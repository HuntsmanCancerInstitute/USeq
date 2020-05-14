package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.data.Graphite;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Tool for contrasting two replicas for common variants.  Uses Graphite for calculating actual counts. */
public class VCFReplicaComparator {

	//user fields
	private File vcfA;
	private File vcfB;
	private File bamA;
	private File bamB;
	private File graphiteExe;
	private File outputDir;
	private File fasta;
	private float minRatio = 0.5f;
	private int minAltCounts = 3;
	
	//internal 
	private ArrayList<VCFMatch> vcfMatches = new ArrayList<VCFMatch>();
	private ArrayList<VCFRecord> vcfNonMatchesA = new ArrayList<VCFRecord>();
	private ArrayList<VCFRecord> vcfNonMatchesB = new ArrayList<VCFRecord>();
	private ArrayList<VCFRecord> vcfNonMatchesANowMatched = null;
	private ArrayList<VCFRecord> vcfNonMatchesBNowMatched = null;

	private VCFParser vcfParserA = null;
	private VCFParser vcfParserB = null;
	private double numVcfA = 0;
	private double numVcfB = 0;
	private Graphite graphite = null;
	private File graphiteTempDir = null;
	
	
	//constructor
	public VCFReplicaComparator(String[] args){
		try {
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		parseVcfs();
		
		compareInitialCalls();
		
		compareNoMatchCalls();
		
		writeVcfs();

		//finish and calc run time
		IO.deleteDirectory(graphiteTempDir);
		
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem comparing vcf replicas.");
		}
	}



	private void writeVcfs() throws IOException {
		IO.pl("\nWriting vcfs...");
		// vcfNonMatchesA
		if (vcfNonMatchesA.size() > 0) {
			File v = new File (outputDir, Misc.removeExtension(vcfA.getName())+"_NoMatch.vcf.gz");
			writeVcf(v,vcfNonMatchesA, vcfParserA.getStringComments());
		}
		// vcfNonMatchesB
		if (vcfNonMatchesB.size() > 0) {
			File v = new File (outputDir, Misc.removeExtension(vcfB.getName())+"_NoMatch.vcf.gz");
			writeVcf(v,vcfNonMatchesB, vcfParserB.getStringComments());
		}
		// matches
		ArrayList<VCFRecord> matches = new ArrayList<VCFRecord>();
		matches.addAll(vcfNonMatchesANowMatched);
		matches.addAll(vcfNonMatchesBNowMatched);
		// pick the one with the highest QUAL
		for (VCFMatch m: vcfMatches) {
			if (m.getKey().getQuality() > m.getTest().getQuality()) matches.add(m.getKey());
			else matches.add(m.getTest());
		}
		if (matches.size()>0) {
			VCFRecord[] r = new VCFRecord[matches.size()];
			matches.toArray(r);
			Arrays.sort(r);
			matches.clear();
			for (VCFRecord v: r) matches.add(v);
			String name = Misc.removeExtension(vcfA.getName())+"_"+ Misc.removeExtension(vcfB.getName()); 
			File f = new File (outputDir, name+"_Match.vcf.gz");
			writeVcf(f,matches, vcfParserA.getStringComments());
		}
		
		double numMatch = matches.size();
		double fracANonMatch = ((double)vcfNonMatchesA.size())/numVcfA;
		double fracBNonMatch = ((double)vcfNonMatchesB.size())/numVcfB;
		double fracAMatch = 1.0 - fracANonMatch;
		double fracBMatch = 1.0 - fracBNonMatch;
		IO.pl("\nFinal Intersection Statistics:");
		IO.pl("\t"+ (int)numMatch+ "\tMatch - A: "+Num.formatPercentOneFraction(fracAMatch)+"  B: "+Num.formatPercentOneFraction(fracBMatch));
		IO.pl("\t"+vcfNonMatchesA.size()+"\tUnique to A "+Num.formatPercentOneFraction(fracANonMatch));
		IO.pl("\t"+vcfNonMatchesB.size()+"\tUnique to B "+Num.formatPercentOneFraction(fracBNonMatch));
	}

	private void writeVcf(File vcfFile, ArrayList<VCFRecord> records, String[] header) throws IOException {
		Gzipper out = new Gzipper(vcfFile);
		for (String s: header) out.println(s);
		for (VCFRecord r: records) out.println(r.getOriginalRecord());
		out.close();		
	}



	private void compareNoMatchCalls() throws IOException {
		IO.pl("\nAdjudicating non matching variants with Graphite...");
		//might be null if no nonMatches
		File noMatchA =  writeOutNoMatchVcf(vcfNonMatchesA, vcfParserA);
		//For A vcfs not found in B
		if (noMatchA != null) {
			//Fetch AF and DP for A vars in A bam
			IO.pl("\tPulling AF DP from Bam A");
			HashMap<String, float[]> keyAfDp_AVarsInABam = graphite.annotate(noMatchA, bamA);
			IO.pl("\tPulling AF DP from Bam B");
			HashMap<String, float[]> keyAfDp_AVarsInBBam = graphite.annotate(noMatchA, bamB);
			ArrayList<VCFRecord>[] matchNoMatchA = adjudicate(keyAfDp_AVarsInABam, keyAfDp_AVarsInBBam, vcfNonMatchesA);
			// assign no matches
			vcfNonMatchesA = matchNoMatchA[1];
			IO.pl("\t"+matchNoMatchA[1].size()+"\tStill Unique to A");
			// assign matches
			vcfNonMatchesANowMatched = matchNoMatchA[0];
			IO.pl("\t"+matchNoMatchA[0].size()+"\tA vars assigned as matching in B");	
		}
		
		File noMatchB =  writeOutNoMatchVcf(vcfNonMatchesB, vcfParserB);
		//For B vcfs not found in A
		if (noMatchB != null) {
			//Fetch AF and DP for B vars in B bam
			IO.pl("\tPulling AF DP from Bam B");
			HashMap<String, float[]> keyAfDp_BVarsInBBam = graphite.annotate(noMatchB, bamB);
			IO.pl("\tPulling AF DP from Bam A");
			HashMap<String, float[]> keyAfDp_BVarsInABam = graphite.annotate(noMatchB, bamA);
			ArrayList<VCFRecord>[] matchNoMatchB = adjudicate(keyAfDp_BVarsInBBam, keyAfDp_BVarsInABam, vcfNonMatchesB);
			// assign no matches
			vcfNonMatchesB = matchNoMatchB[1];
			IO.pl("\t"+matchNoMatchB[1].size()+"\tStill Unique to B");
			// assign matches
			vcfNonMatchesBNowMatched = matchNoMatchB[0];
			IO.pl("\t"+matchNoMatchB[0].size()+"\tB vars assigned as matching in A");	
		}
		
	}

	private void parseVcfs() {
		IO.pl("\nInitial Intersection Statistics:");
		vcfParserA = new VCFParser(vcfA, true, false, false);
		numVcfA = vcfParserA.getVcfRecords().length;
		IO.pl("\t"+vcfParserA.getVcfRecords().length+"\tin A");
		vcfParserB = new VCFParser(vcfB, true, false, false);
		numVcfB = vcfParserB.getVcfRecords().length;
		IO.pl("\t"+vcfParserB.getVcfRecords().length+"\tin B");
	}

	private ArrayList<VCFRecord>[] adjudicate(HashMap<String, float[]> keyAfDpCalled, HashMap<String, float[]> keyAfDpTest, ArrayList<VCFRecord> vcfNonMatches) throws IOException {
		 ArrayList<VCFRecord> stillNoMatch = new  ArrayList<VCFRecord>();
		 ArrayList<VCFRecord> match = new  ArrayList<VCFRecord>();
	
		//for each variant
		for (int i=0; i< vcfNonMatches.size(); i++) {
			VCFRecord v = vcfNonMatches.get(i);
			//#CHROM POS ID REF ALT 0 1 2 3 4
			String key = v.getChrPosRefAlt(false);
			float[] afDpCalled =  keyAfDpCalled.get(key);
			float[] afDpTest =  keyAfDpTest.get(key);
			if (afDpCalled == null || afDpTest == null) throw new IOException("Failed to find the AF and DP from the the called or test graphite results, see "+key);
			
			//pass min alt counts?
			int altCountsInTest = Math.round(afDpTest[0] * afDpTest[1]);
			if (altCountsInTest< minAltCounts) stillNoMatch.add(v);
			else {
				//pass test/called AF ratio 
				if (afDpCalled[0] <= 0.0f) throw new IOException("Error, the AF in the called variant is <= 0 "+key+"\t"+afDpCalled[0]+","+afDpCalled[1]+"\t"+
						afDpTest[0]+","+afDpTest[1]+"\t"+altCountsInTest+"\t"+(afDpTest[0]/afDpCalled[0]));
				float ratio = afDpTest[0]/afDpCalled[0];
				if (ratio < minRatio) stillNoMatch.add(v);
				else match.add(v);
			}	
		}
		return new ArrayList[] {match, stillNoMatch};
	}

	private File writeOutNoMatchVcf(ArrayList<VCFRecord> vcfs, VCFParser vcfParser) throws IOException {
		if (vcfs.size() == 0) return null;
		String name = vcfParser.getVcfFile().getName();
		name = Misc.removeExtension(name);
		File file = new File(outputDir, name+ "_NoMatchTemp.vcf");
		file.deleteOnExit();
		PrintWriter out = new PrintWriter(new FileWriter(file ));
		for (String s: vcfParser.getStringComments()) out.println(s);
		for (VCFRecord r: vcfs) out.println(r.getOriginalRecord());
		out.close();		
		return file;
	}


	public void compareInitialCalls(){	
		//intersect and split test into matching and non matching
		intersectVCF();
		double numMatch = vcfMatches.size();
		
		double fracANonMatch = ((double)vcfNonMatchesA.size())/numVcfA;
		double fracBNonMatch = ((double)vcfNonMatchesB.size())/numVcfB;
		double fracAMatch = 1.0 - fracANonMatch;
		double fracBMatch = 1.0 - fracBNonMatch;
		
		IO.pl("\t"+ (int)numMatch+ "\tMatch - A: "+Num.formatPercentOneFraction(fracAMatch)+"  B: "+Num.formatPercentOneFraction(fracBMatch));
		IO.pl("\t"+vcfNonMatchesA.size()+"\tUnique to A "+Num.formatPercentOneFraction(fracANonMatch));
		IO.pl("\t"+vcfNonMatchesB.size()+"\tUnique to B "+Num.formatPercentOneFraction(fracBNonMatch));	
	}

	public void intersectVCF(){
		//set all records to fail
		vcfParserA.setFilterFieldOnAllRecords(VCFRecord.FAIL);
		vcfParserB.setFilterFieldOnAllRecords(VCFRecord.FAIL);

		//for each A record
		for (String chr: vcfParserA.getChromosomeVCFRecords().keySet()){
			VCFLookUp luA = vcfParserA.getChromosomeVCFRecords().get(chr);
			VCFLookUp luB = vcfParserB.getChromosomeVCFRecords().get(chr);
			//any vcf records in B?
			if (luB != null) countMatches(luB, luA, vcfMatches);
		} 
		
		//record non matches
		VCFRecord[] vcfs = vcfParserA.getVcfRecords();
		for (int i=0; i< vcfs.length; i++) if (vcfs[i].getFilter().equals(VCFRecord.FAIL)) vcfNonMatchesA.add(vcfs[i]);
		vcfs = vcfParserB.getVcfRecords();
		for (int i=0; i< vcfs.length; i++) if (vcfs[i].getFilter().equals(VCFRecord.FAIL)) vcfNonMatchesB.add(vcfs[i]);
	}


	public void countMatches(VCFLookUp luA, VCFLookUp luB, ArrayList<VCFMatch> matches){
		
		int[] posB = luB.getBasePosition();
		VCFRecord[] vcfB = luB.getVcfRecord();
		
		//for each B record 
		for (int i=0; i< vcfB.length; i++){
			
			//fetch the A records that intersect
			int start = posB[i];
			int stop = posB[i]+1;
			VCFRecord[] aRecords = luA.fetchVCFRecords(start, stop);	
			
			//no records found
			if (aRecords == null) continue;
			
			//for each A that matches pos
			else {				
				for (int x=0; x< aRecords.length; x++){	
					//check to see if it matches by position, ref and by one of the alts
					if (vcfB[i].matchesPosRefAlt(aRecords[x])) {
						aRecords[x].setFilter(VCFRecord.PASS);
						vcfB[i].setFilter(VCFRecord.PASS);
						matches.add(new VCFMatch(aRecords[x], vcfB[i]));
						break;
					}
				}
			}
		}
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFReplicaComparator(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		Pattern pat = Pattern.compile("-[a-z0-9]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': vcfA = new File(args[++i]); break;
					case 'b': vcfB = new File(args[++i]); break;
					case '1': bamA = new File(args[++i]); break;
					case '2': bamB = new File(args[++i]); break;
					case 'g': graphiteExe = new File(args[++i]); break;
					case 'o': outputDir = new File(args[++i]); break;
					case 'f': fasta = new File(args[++i]); break;
					case 'c': minAltCounts = Integer.parseInt(args[++i]); break;
					case 'r': minRatio = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check files
		if (vcfA == null || vcfB == null || vcfA.exists()==false || vcfB.exists()==false) Misc.printErrAndExit("\nError -a or -b: please provide paths to the first vcf ("+vcfA+") and second vcf ("+vcfB+") files.");
		if (bamA == null || bamB == null || bamA.exists()==false || bamB.exists()==false) Misc.printErrAndExit("\nError -1 or -2: please provide paths to the first bam ("+bamA+") and second bam ("+bamB+") files.");
		if (graphiteExe == null || graphiteExe.canExecute() == false) Misc.printErrAndExit("\nError -g: please provide path to the graphite executable.");
		if (outputDir == null) Misc.printErrAndExit("\nError -o: please provide a path to a directory to write the adjudicated vcf files.\n");
		if (fasta == null || fasta.exists()==false) Misc.printErrAndExit("\nError -f: please provide an indexed fasta file.\n");
		if (outputDir.exists() == false) outputDir.mkdirs();
		
		//make temp dir for Graphite
		graphiteTempDir = new File (outputDir, "GraphiteTempDir_"+Misc.getRandomString(6));
		graphiteTempDir = graphiteTempDir.getCanonicalFile();
		graphiteTempDir.mkdirs();
		graphite = new Graphite(graphiteExe, fasta, graphiteTempDir);
		
		printOptions();
	}	
	
	public void printOptions() {
		IO.pl("Parameters:");
		IO.pl(" -a VcfA "+vcfA);
		IO.pl(" -b VcfB "+vcfB);
		IO.pl(" -1 BamA "+bamA);
		IO.pl(" -2 BamB "+bamB);
		IO.pl(" -o OutputDir\t"+outputDir);
		IO.pl(" -c MinAltCounts\t"+minAltCounts);
		IO.pl(" -r MinAFRatio Test/RealCall "+minRatio);
		IO.pl(" -f FastaIndex "+fasta);
		IO.pl(" -g Graphite "+graphiteExe);
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           VCF Replica Comparator : April 2020                    **\n" +
				"**************************************************************************************\n" +
				"Compares two vcf files derived from a multi replica experiment. Variants unique to\n"+
				"either are run through Graphite (https://github.com/dillonl/graphite) to accurately\n"+
				"calculate alignment support and potentially add more matches based on the -c and -r\n"+
				"thresholds.\n"+

				"\nRequired Options:\n"+
				"-a Vcf file for replica A (xxx.vcf(.gz/.zip OK))\n"+
				"-b Vcf file for replica B\n"+
				"-1 Bam file for replica A (e.g. tumor) with its index\n"+
				"-2 Bam file for replica B\n"+
				"-o Output directory for adjudicated variant calls\n"+
				"-f Fasta index used in aligning the reads\n"+
				"-g Graphite executable, see https://github.com/dillonl/graphite install e.g.:\n"+
				"     module load gcc/4.9.2\n" + 
				"     git clone https://github.com/dillonl/graphite.git\n" + 
				"     mkdir graphite/bin\n" + 
				"     cd graphite/bin\n" + 
				"     cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++)\n" + 
				"     make\n"+
				
				"\nDefault Options:\n"+
				"-c Min alt allele variant observations to match, defaults to 3\n"+
				"-r Min test AF/ called AF ratio to match, defaults to 0.5\n"+
				
				"\nExample: module load gcc/4.9.2; java -Xmx10G -jar USeq/Apps/VCFReplicaComparator\n" +
				"       -a repA.vcf.gz -b repB.vcf.gz -1 repA.bam -2 repB.bam -o AdjVcfs \n" +
				"       -g ~/Graphite/bin/graphite -f ~/Hg38Ref/hs38DH.fa -c 2 -r 0.3 \n\n"+

		"**************************************************************************************\n");

	}
}