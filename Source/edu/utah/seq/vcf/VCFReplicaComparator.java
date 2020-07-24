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
import edu.utah.seq.parsers.jpileup.BamPileup;
import edu.utah.seq.parsers.jpileup.BamPileupTabixLoader;
import edu.utah.seq.parsers.jpileup.BaseCount;
import edu.utah.seq.parsers.jpileup.BpileupLine;
import edu.utah.seq.query.QueryIndexFileLoader;
import htsjdk.tribble.readers.TabixReader;
import util.apps.MergeRegions;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Tool for contrasting two replicas for common variants.  Uses Graphite or BamPileup for calculating actual counts. */
public class VCFReplicaComparator {

	//user fields
	private File vcfA;
	private File vcfB;
	private File bamA;
	private File bamB;
	private File graphiteExe;
	private File outputDir;
	private File fasta;
	private File tabix;
	private File bgzip;
	private double minRatio = 0.5;
	private int minAltCounts = 3;
	private float minQualFilt = 0;
	private String filterTxtToExclude = null;
	
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
	private int bpPadding = 200;
	private TabixReader tabixReader = null;
	
	
	//constructor
	public VCFReplicaComparator(String[] args) {
		try {
			//start clock
			long startTime = System.currentTimeMillis();

			//process args
			processArgs(args);

			parseVcfs();

			compareInitialCalls();

			compareNoMatchCallsBamPileup();

			//compareNoMatchCallsGraphite();

			writeVcfs();

			//finish and calc run time
			//IO.deleteDirectory(graphiteTempDir);

			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
		} catch (Exception e) {
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
		IO.pl("\t"+ (int)numMatch+ "\tMatch w/ extras - A: "+Num.formatPercentOneFraction(fracAMatch)+"  B: "+Num.formatPercentOneFraction(fracBMatch));
		IO.pl("\t"+vcfNonMatchesA.size()+"\tUnique to A "+Num.formatPercentOneFraction(fracANonMatch));
		IO.pl("\t"+vcfNonMatchesB.size()+"\tUnique to B "+Num.formatPercentOneFraction(fracBNonMatch));
	}

	private void writeVcf(File vcfFile, ArrayList<VCFRecord> records, String[] header) throws IOException {
		Gzipper out = new Gzipper(vcfFile);
		for (String s: header) out.println(s);
		for (VCFRecord r: records) out.println(r.getOriginalRecord());
		out.close();		
	}

	private void compareNoMatchCallsGraphite() throws IOException {
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
		int numAPre = vcfParserA.getVcfRecords().length;
		vcfParserA.getVcfRecords();
		if (minQualFilt !=0.0f) vcfParserA.filterVCFRecordsOnQUAL(minQualFilt);
		if (filterTxtToExclude != null) vcfParserA.filterVCFRecordsOnFILTER(filterTxtToExclude);
		numVcfA = vcfParserA.getVcfRecords().length;
		IO.pl("\t"+(int)numVcfA+"\tin A ("+numAPre+" before filtering)");
		
		vcfParserB = new VCFParser(vcfB, true, false, false);
		int numBPre = vcfParserB.getVcfRecords().length;
		if (minQualFilt !=0.0f) vcfParserB.filterVCFRecordsOnQUAL(minQualFilt);
		if (filterTxtToExclude != null) vcfParserB.filterVCFRecordsOnFILTER(filterTxtToExclude);
		numVcfB = vcfParserB.getVcfRecords().length;
		IO.pl("\t"+(int)numVcfB+"\tin B ("+numBPre+" before filtering)");
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
	
	public void compareNoMatchCallsBamPileup() throws Exception {
		IO.pl("\nAdjudicating non matching variants with BamPileup...");
		
		if (vcfNonMatchesA.size() == 0 && vcfNonMatchesA.size() == 0) return;
			
			//write out bed file of regions around each vcf record
			File mergedBed = writeBedFromVcfs();

			//run BamPileup with sample order A then B, 0 and 1 respectively
			File bpTempResults = createBamPileup(mergedBed);
			
			//create tabix index reader
			tabixReader = new TabixReader(bpTempResults.getCanonicalPath());
			
			//for A vars not in B
			if (vcfNonMatchesA.size() != 0) {
				ArrayList<VCFRecord>[] matchNoMatch = scoreVcfs(vcfNonMatchesA, 0);
				// assign no matches
				vcfNonMatchesA = matchNoMatch[1];
				IO.pl("\t"+matchNoMatch[1].size()+"\tStill Unique to A");
				// assign matches
				vcfNonMatchesANowMatched = matchNoMatch[0];
				IO.pl("\t"+matchNoMatch[0].size()+"\tA vars assigned as matching in B");
			}
			
			//for B vars not in A
			if (vcfNonMatchesB.size() != 0) {
				ArrayList<VCFRecord>[] matchNoMatch = scoreVcfs(vcfNonMatchesB, 1);
				// assign no matches
				vcfNonMatchesB = matchNoMatch[1];
				IO.pl("\t"+matchNoMatch[1].size()+"\tStill Unique to B");
				// assign matches
				vcfNonMatchesBNowMatched = matchNoMatch[0];
				IO.pl("\t"+matchNoMatch[0].size()+"\tB vars assigned as matching in A");	
			}
			
			//cleanup
			bpTempResults.deleteOnExit();
			tabixReader.close();
				
	
	}
	
	private File createBamPileup(File mergedBed) throws IOException {
		File bpTempResults = new File (outputDir, "tempNonMatch.bp.txt.gz");
		bpTempResults.deleteOnExit();
		File tbi = new File (outputDir, "tempNonMatch.bp.txt.gz.tbi");
		tbi.deleteOnExit();
		new BamPileup(mergedBed, new File[] {bamA, bamB}, outputDir, fasta, bpTempResults, bgzip, tabix, false);
		if (bpTempResults.exists() == false) throw new IOException("ERROR: failed to find "+bpTempResults);
		return bpTempResults;
	}

	private File writeBedFromVcfs() throws IOException {
		File bed = new File(outputDir, "tempNonMatchVcf_"+Misc.getRandomString(6)+".bed");
		bed.deleteOnExit();
		PrintWriter out = new PrintWriter( new FileWriter (bed));
		for (VCFRecord v: vcfNonMatchesA) {
			int[] startStop = QueryIndexFileLoader.fetchEffectedBps(Misc.TAB.split(v.getOriginalRecord()), false);
			int start = (startStop[0]-bpPadding);
			if (start < 0) start = 0;
			out.println(v.getChromosome()+"\t"+ start+ "\t"+(startStop[1]+bpPadding));
		}
		for (VCFRecord v: vcfNonMatchesB) {
			int[] startStop = QueryIndexFileLoader.fetchEffectedBps(Misc.TAB.split(v.getOriginalRecord()), false);
			int start = (startStop[0]-bpPadding);
			if (start < 0) start = 0;
			out.println(v.getChromosome()+"\t"+ start+ "\t"+(startStop[1]+bpPadding));
		}
		out.close();
		
		//merge and sort the regions
		File mergedBed = new File(outputDir, "tempNonMatchVcfMerged_"+Misc.getRandomString(6)+".bed");
		mergedBed.deleteOnExit();
		new MergeRegions (new File[] {bed}, mergedBed, false);
		return mergedBed;
	}

	private ArrayList<VCFRecord>[]  scoreVcfs(ArrayList<VCFRecord> vcfs, int index) throws Exception {
		// if index is 0 these are A variants not found in B so look in bp sample 1 (B)
		// if index is 1 these are B variants not found in A so look in bp sample 0 (A)
		ArrayList<VCFRecord> stillNoMatch = new  ArrayList<VCFRecord>();
		ArrayList<VCFRecord> match = new  ArrayList<VCFRecord>();

		for (VCFRecord v: vcfs) {
			//what kind of variant, returns GATCID
			String[] fields = Misc.TAB.split(v.getOriginalRecord());
			char allele = BamPileupTabixLoader.fetchAllele(fields);
			Boolean pass;
			if (allele == 'I' || allele == 'D') pass = scoreIndel(fields, allele, index);
			else pass = scoreSnv(fields, allele, index);

			if (pass == null) {
				IO.pl("\nWARNING: failed to score, assigning to still no match "+v);
				stillNoMatch.add(v);
			}
			else if (pass == true) {
//IO.pl("\t\t\tPass");
				match.add(v);
			}
			else {
//IO.pl("\t\t\tFail");
				stillNoMatch.add(v);
			}
		}
		return new ArrayList[] {match, stillNoMatch};
	}

	private Boolean scoreIndel(String[] fields, char allele, int index) throws Exception {
//IO.pl("Scoring INDEL "+Misc.stringArrayToString(fields, " ")+" "+allele +" "+index);

		String tabixCoor = null;
		//fetch interbase coordinates for del
		if (allele == 'D') {
			int[] startStop = QueryIndexFileLoader.fetchEffectedBps(fields, true);
			if (startStop == null) return null; 
			//pull bpileup records over deleted bps
			tabixCoor = fields[0]+":"+(startStop[0]+2)+"-"+ (startStop[1]);
		}
		else if (allele ==  'I') {
			//single downstream base
			int start = Integer.parseInt(fields[1])+1;
			int stop = start;
			tabixCoor = fields[0]+":"+start+"-"+stop;
		}

		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) return null;

		//for each bpileup record, check minAltCounts and minRatio
		String bpileupLine = null;

		while ((bpileupLine = it.next()) != null){
//IO.pl("\tLine "+bpileupLine);

			//parse bpileup line and pull filtered set of samples
			BpileupLine ml = new BpileupLine(bpileupLine);
			BaseCount[] abSamples = ml.getSamples();

			//check minAltCount
//IO.pl("\t\tAltCounts "+abSamples[0].getIndelCount(allele)+"\t"+abSamples[1].getIndelCount(allele));
			boolean pass = false;
			if (index == 0) {
				int bAltCount = abSamples[1].getIndelCount(allele);
				if (bAltCount >= minAltCounts) pass = true;
			}
			else {
				int aAltCount = abSamples[0].getIndelCount(allele);
				if (aAltCount >= minAltCounts) pass = true;
			}
			if (pass == false) continue;

			pass = false;
			//check AF ratio
			double aAF = abSamples[0].getIndelAlleleFreq(allele);
			double bAF = abSamples[1].getIndelAlleleFreq(allele);

			if (index == 0) {
//IO.pl("\t\tAltAFs xxx "+aAF+"\t"+bAF+" b/a "+bAF/aAF);
				if ((bAF/aAF) >= minRatio) return new Boolean (true);
			}
			else {
//IO.pl("\t\tAltAFs xxx "+aAF+"\t"+bAF+" a/b "+aAF/bAF);
				if ((aAF/bAF) >= minRatio) return new Boolean (true);
			}
		}
		//got through all the lines and all failed
		return new Boolean (false);

	}

	/**Returns null if it cannot score the variant*/
	private Boolean scoreSnv(String[] fields, char allele, int index) throws Exception {
//IO.pl("Scoring Snv "+Misc.stringArrayToString(fields, " ")+" "+allele+" "+index);
		// if index is 0 these are A variants not found in B so look in bp sample 1 (B)
		// if index is 1 these are B variants not found in A so look in bp sample 0 (A)

		//pull bpileup record, if none return null;
		String tabixCoor = fields[0]+":"+fields[1]+"-"+fields[1];
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) return null;

		String bpileupLine = it.next();
		if (bpileupLine == null) return null;

		//parse bpileup line 
		BpileupLine ml = new BpileupLine(bpileupLine);
		BaseCount[] abSamples = ml.getSamples();

		//check minAltCount


//IO.pl("\t\tAltCounts "+abSamples[0].getSnvCount(allele)+"\t"+abSamples[1].getSnvCount(allele));

		if (index == 0) {
			int bAltCount = abSamples[1].getSnvCount(allele);
			if (bAltCount < minAltCounts) return new Boolean (false);
		}
		else {
			int aAltCount = abSamples[0].getSnvCount(allele);
			if (aAltCount < minAltCounts) return new Boolean (false);
		}
		//check AF ratio
		double aAF = abSamples[0].getSnvAlleleFreq(allele);
		double bAF = abSamples[1].getSnvAlleleFreq(allele);

		if (index == 0) {
//IO.pl("\t\tAltAFs "+aAF+"\t"+bAF+" b/a "+bAF/aAF);
			if ((bAF/aAF) < minRatio) return new Boolean (false);
		}
		else {
//IO.pl("\t\tAltAFs "+aAF+"\t"+bAF+" a/b "+aAF/bAF);
			if ((aAF/bAF) < minRatio) return new Boolean (false);
		}
		return new Boolean (true);
	}
	
	private TabixReader.Iterator fetchInteratorOnCoordinates(String coordinates) {
		TabixReader.Iterator it = null;
		//watch out for no retrieved data error from tabix
		try {
			it = tabixReader.query(coordinates);
		} catch (ArrayIndexOutOfBoundsException e){
		}
		return it;
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
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File tabixBinDirectory = null;
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
					//case 'g': graphiteExe = new File(args[++i]); break;
					case 'o': outputDir = new File(args[++i]); break;
					case 'i': fasta = new File(args[++i]); break;
					case 'c': minAltCounts = Integer.parseInt(args[++i]); break;
					case 'r': minRatio = Double.parseDouble(args[++i]); break;
					case 'f': filterTxtToExclude = args[++i]; break;
					case 'q': minQualFilt = Float.parseFloat(args[++i]); break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
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
		//if (graphiteExe == null || graphiteExe.canExecute() == false) Misc.printErrAndExit("\nError -g: please provide path to the graphite executable.");
		if (outputDir == null) Misc.printErrAndExit("\nError -o: please provide a path to a directory to write the adjudicated vcf files.\n");
		if (fasta == null || fasta.exists()==false) Misc.printErrAndExit("\nError -i: please provide an indexed fasta file.\n");
		if (outputDir.exists() == false) outputDir.mkdirs();
		
		//pull tabix and bgzip
		if (tabixBinDirectory == null) Misc.printExit("\nError: please point to your HTSlib directory containing the tabix and bgzip executables (e.g. ~/BioApps/HTSlib/1.10.2/bin/ )\n");
		bgzip = new File (tabixBinDirectory, "bgzip");
		tabix = new File (tabixBinDirectory, "tabix");
		if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+tabixBinDirectory);

		//make temp dir for Graphite
		//int numThreads = Runtime.getRuntime().availableProcessors() - 1;
		//graphiteTempDir = new File (outputDir, "GraphiteTempDir_"+Misc.getRandomString(6));
		//graphiteTempDir = graphiteTempDir.getCanonicalFile();
		//graphiteTempDir.mkdirs();
		//graphite = new Graphite(graphiteExe, fasta, graphiteTempDir, numThreads);
		
		printOptions();
	}	
	
	public void printOptions() {
		IO.pl("Parameters:");
		IO.pl(" -a VcfA "+vcfA);
		IO.pl(" -b VcfB "+vcfB);
		IO.pl(" -1 BamA "+bamA);
		IO.pl(" -2 BamB "+bamB);
		IO.pl(" -o OutputDir\t"+outputDir);
		IO.pl(" -i FastaIndex\t"+fasta);
		IO.pl(" -t TabixDir\t"+tabix.getParent());
		IO.pl(" -c MinAltCounts\t"+minAltCounts);
		IO.pl(" -r MinAFRatio Test/RealCall "+minRatio);
		IO.pl(" -q Min QUAL "+minQualFilt);
		IO.pl(" -f FILTER exclude "+filterTxtToExclude);

		//IO.pl(" -g Graphite app "+graphiteExe);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           VCF Replica Comparator : July 2020                     **\n" +
				"**************************************************************************************\n" +
				"Compares two vcf files derived from a multi replica experiment. Variants unique to\n"+
				"either are run through BamPilup to estimate alignment support and potentially add more\n"+
				"replica matches based on the -c and -r thresholds. Thus filter the input VCFs for min\n"+
				"QUAL and min alt allele coverage and allele freq before running.\n"+

				"\nRequired Options:\n"+
				"-a Vcf file for replica A (xxx.vcf(.gz/.zip OK))\n"+
				"-b Vcf file for replica B\n"+
				"-1 Path to a merged bam file for replica A (e.g. tumor) with its index\n"+
				"-2 Path to a merged bam file for replica B\n"+
				"-o Output directory for adjudicated variant calls\n"+
				"-i Fasta index used in aligning the reads\n"+
				//"-g Graphite executable, see https://github.com/dillonl/graphite install e.g.:\n"+
				//"     module load gcc/4.9.2\n" + 
				//"     git clone https://github.com/dillonl/graphite.git\n" + 
				//"     mkdir graphite/bin\n" + 
				//"     cd graphite/bin\n" + 
				//"     cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++)\n" + 
				//"     make\n"+
				"-t Path to a directory containing the bgzip and tabix executables to compress and index\n"+
				"      the bp file, see htslib.org \n"+
				
				"\nDefault Options:\n"+
				"-c Min alt allele variant observations in replica to match, defaults to 3\n"+
				"-r Min test AF/ called AF ratio to match, defaults to 0.5\n"+
				"-q Min QUAL score for vcf record import, defaults to 0\n"+
				"-f FILTER field value to exclude records from import, defaults to no filtering\n"+
				
				"\nExample: module load gcc/4.9.2; java -Xmx25G -jar USeq/Apps/VCFReplicaComparator\n" +
				"       -a repA.vcf.gz -b repB.vcf.gz -1 repA.bam -2 repB.bam -o AdjVcfs -f BKAF\n" +
				"       -t ~/BioApps/HTSlib/1.10.2/bin/ -i ~/Hg38Ref/hs38DH.fa -c 5 -r 0.75 -q 3\n\n"+
				//"       -g ~/Graphite/bin/graphite -i ~/Hg38Ref/hs38DH.fa -c 5 -r 0.75 -q 3\n\n"+

		"**************************************************************************************\n");

	}
}