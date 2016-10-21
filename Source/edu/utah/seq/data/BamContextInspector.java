package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.Bed;
import util.gen.*;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamLayoutForMutation;

/** Application for identifying flanking bps with a potential variant that could be confounding the call.
 * @author Nix
 * */
public class BamContextInspector {

	//user defined fields
	private File bedIn;
	private File bamFile;
	private int bpToScan = 25;
	private double maxNonRefBaseFreq = 0.03;
	private int minReadCov = 6;
	private int minNonRef = 2;
	private int minBaseQual = 13;
	private int minMapQual = 13;
	private int maxNumNonRefBps = 1;
	private boolean skipIndels = true;
	
	private IndexedFastaSequenceFile fasta; 
	private Bed[] regions;
	private SamReader samReader = null;
	private ArrayList<Bed> passAL = new ArrayList<Bed>();
	private ArrayList<Bed> failAL = new ArrayList<Bed>();
	private File geneiAseFile;

	//constructor
	public BamContextInspector(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			printParams();

			//load the bed file of regions to scan
			if (bedIn != null) regions = Bed.parseFile(bedIn, 0, 0);
			else parseGeneiAseTable();

			//create reader
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			samReader = factory.open(bamFile);

			walkRegions();
			
			//close the readers
			fasta.close();
			samReader.close();
			
			printRegions();
			System.out.println(passAL.size()+"\tPassing");
			System.out.println(failAL.size()+"\tFailing");
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes\n");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void parseGeneiAseTable() {
		try {
			BufferedReader in = IO.fetchBufferedReader(geneiAseFile);
			String line;
			String[] t;
			String[] c;
			ArrayList<Bed> al = new ArrayList<Bed>();
			while ((line = in.readLine()) != null){
				if (line.startsWith("gene")) continue;
				t = Misc.TAB.split(line);
				c = Misc.UNDERSCORE.split(t[1]);
				int start = Integer.parseInt(c[1]);
				Bed b = new Bed(c[0], start, start+1, line, 0f, '.');
				al.add(b);
			}
			regions = new Bed[al.size()];
			al.toArray(regions);
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private void printParams(){
		System.out.println("Params:");
		if (bedIn != null) System.out.println(bedIn.getName()+"\t Regions file");
		else System.out.println(geneiAseFile.getName()+"\t Regions file");
		System.out.println(bamFile.getName()+"\t Bam file");
		System.out.println(bpToScan+"\t BPs of flank to scan");
		System.out.println(maxNumNonRefBps+ "\t Max num non ref BPs in flanks");
		System.out.println(maxNonRefBaseFreq+"\tMax AF");
		System.out.println(minReadCov+"\tMin read cov");
		System.out.println(minNonRef+"\tMin non ref cov");
		System.out.println(minBaseQual+"\tMin base qual");
		System.out.println(minMapQual+"\tMin map qual");
		System.out.println(skipIndels +"\tFail regions with 1+ INDEL(s) in flanks\n");
	}
	
	private int score(BaseCount[] bc) {
		int failing = 0;
		//for each base
		for (int i=0; i< bc.length; i++){
			double depth = bc[i].getReadDepth();
			if (depth < minReadCov) continue;
			double numNonRef = bc[i].getReadDepthNonRef();
			if (numNonRef < minNonRef) continue;
			double af = numNonRef/depth;
			if (af >= maxNonRefBaseFreq) failing++;
			else if (skipIndels) {
				//check indel status
				double numIndel = bc[i].getIndelCount();
				double afIndel = numIndel/depth;
				if (numIndel > minNonRef && afIndel >= maxNonRefBaseFreq) return Integer.MAX_VALUE;
			}
		}
		return failing;
	}


	private void printRegions() {
		String name;
		File parent;
		if (bedIn != null) {
			parent = bedIn.getParentFile();
			name = bedIn.getName();
		}
		else {
			parent = geneiAseFile.getParentFile();
			name = geneiAseFile.getName();
		}
		if (passAL.size() !=0) write (new File (parent, "pass_"+name), passAL);
		if (failAL.size() !=0) write (new File (parent, "fail_"+name), failAL);
	}

	private void write(File f, ArrayList<Bed> al) {
		try {
			Gzipper out = new Gzipper(f);
			for (Bed b: al) {
				if (bedIn != null) out.println(b.toString());
				else out.println(b.getName());
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void walkRegions() throws Exception {
		//for each region

		
		for (Bed region : regions){
			
			//scan upstream
			int upStart = region.getStart()- bpToScan;
			int upStop = region.getStart();
			
			BaseCount[] up = pileup(region.getChromosome(), upStart, upStop);
			int numNonRef = score(up);
			
			//System.out.println("P\tR\tG\tA\tT\tC\tI\tD");
			//for (BaseCount bc : up) System.out.println(bc);
			if (numNonRef > maxNumNonRefBps) failAL.add(region);
			
			//scan downstream if needed
			else {
				int downStart = region.getStop();
				int downStop = downStart+ bpToScan;
				BaseCount[] down = pileup(region.getChromosome(), downStart, downStop);
				numNonRef += score(down);
				if (numNonRef > maxNumNonRefBps) failAL.add(region);
				else passAL.add(region);
			}
		}
		
	}



	private class BaseCount{
		int bpPosition;
		int g = 0;
		int a = 0;
		int t = 0;
		int c = 0;
		int ins = 0;
		int del = 0;
		char ref;
		
		BaseCount(int bpPosition, char ref){
			this.bpPosition = bpPosition;
			this.ref = ref;
		}
		
		public double getReadDepthNonRef() {
			double depth = getReadDepth();
			if (ref == 'G') return depth - g;
			if (ref == 'A') return depth - a;
			if (ref == 'T') return depth - t;
			if (ref == 'C') return depth - c;
			return 0;
		}
		public double getReadDepth() {
			return g+a+t+c+ins+del;
		}
		public double getIndelCount() {
			return ins+del;
		}
		public void increment(char x) throws Exception{
			if (x == 'G') g++;
			else if (x == 'A') a++;
			else if (x == 'T') t++;
			else if (x == 'C') c++;
			else if (x == 'I') ins++;
			else if (x == 'D') del++;
			else throw new Exception("Unrecognized base "+x);
		}
		/**pos ref G A T C INS DEL*/
		public String toString(){
			return bpPosition+"\t"+ref+"\t"+g+"\t"+a+"\t"+t+"\t"+c+"\t"+ins+"\t"+del;
		}
	}
	
	private BaseCount[] pileup(String chr, int start, int stop) throws Exception{
		
		//create container for counts
		ReferenceSequence p = fasta.getSubsequenceAt(chr, start+1, stop+1);
		char[] refSeq = new String(p.getBases()).toCharArray();
		BaseCount[] bc = new BaseCount[stop-start];
		int counter = 0;
		for (int i=start; i< stop; i++) bc[counter] = new BaseCount(i, refSeq[counter++]);

		//fetch alignments
		SAMRecordIterator it = samReader.queryOverlapping(chr, start-1, stop+1);
		if (it.hasNext() == false) {
			it.close();
			return bc;
		}

		//for each record
		while (it.hasNext()){

			SAMRecord sam = it.next();
//if (sam.getCigarString().contains("I")) {
	//System.out.println("\nProc "+sam.getSAMString());
//}
		//}
		
	//return null;}

			//pass map quality?
			if (sam.getMappingQuality() < minMapQual) {
				//System.out.println("\tPoor qual skipping");
				continue;
			}

			//make a layout
			SamAlignment sa = new SamAlignment(sam.getSAMString().trim(), true);
			SamLayoutForMutation layout = new SamLayoutForMutation(sa);
			//layout.print();

			//for each base in the region see if it overlaps
			counter = 0;
			for (int i = start; i< stop; i++){
				int index = layout.findBaseIndex(i);
				//System.out.println(i+ " i \t index "+index);

				//present in the alignment?
				if (index == -1) {
					//System.out.println("\tNotFound");
					counter++;
					continue;
				}

				//good score? watch out for -1, not set
				int qual = layout.getQual()[index];
				if (qual != -1 && qual < minBaseQual) {
					//System.out.println("\tPoorQual");
					counter++;
					continue;
				}

				//a match?
				char call = layout.getCall()[index];
				if (call == 'M') {
					//System.out.println("\tMatch "+layout.getSeq()[index]);
					bc[counter].increment(layout.getSeq()[index]);
				}
				//an indel
				else if (call == 'I' || call == 'D') {
					bc[counter].increment(call);
					//System.out.println("\tINDEL "+call);
				}
				//nope must be a masked base (S,H) or N
				//else System.out.println("\tMask or N "+call);
				counter++;
			}
		}
		it.close();
		return bc;
	}

	public static void main(String[] args) {

		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamContextInspector(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args) throws Exception{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File indexedFastaFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bamFile = new File(args[++i]); break;
					case 'r': bedIn = new File(args[++i]); break;
					case 'f': indexedFastaFile = new File(args[++i]); break;
					case 'p': bpToScan = Integer.parseInt(args[++i]); break;
					case 'c': minReadCov = Integer.parseInt(args[++i]); break;
					case 'n': minNonRef = Integer.parseInt(args[++i]); break;
					case 'q': minBaseQual = Integer.parseInt(args[++i]); break;
					case 'm': minMapQual = Integer.parseInt(args[++i]); break;
					case 'a': maxNonRefBaseFreq = Double.parseDouble(args[++i]); break;
					
					case 'x': maxNumNonRefBps = Integer.parseInt(args[++i]); break;
					case 'i': skipIndels = false; break;
					//hidden
					case 'g': geneiAseFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//Create fasta fetcher
		if (indexedFastaFile == null || indexedFastaFile.exists() == false)  Misc.printErrAndExit("\nError: please provide an indexed reference genome fasta file and it's index.");
		fasta = new IndexedFastaSequenceFile(indexedFastaFile);
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFastaFile);
		
		if (bedIn == null && geneiAseFile == null) Misc.printErrAndExit("\nError: please provide a file of regions to scan.");
		
		if (bamFile == null ) Misc.printErrAndExit("\nError: please provide a sorted and indexed bam alignment file.");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Bam Context Inspector:  Oct 2016                       **\n" +
				"**************************************************************************************\n" +
				"Application for scanning the surrounding context of a set of regions for non reference\n"+
				"bps.  Use to flag variants with adjacent potentially confounding changes.\n"+

				"\nRequired Options:\n"+
				"-b Sorted bam alignment file with associated index.\n"+
				"-r Bed file of regions to split into pass (no non refs) or fail (with a non ref bp).\n"+
				"-f Path to the reference fasta with and xxx.fai index\n"+

				"\nDefault Options:\n"+
				"-p Bp of flanking bases to scan for non ref bases, defaults to 25\n"+
				"-c Minimum alignment coverage of a base before scanning, defaults to 6\n"+
				"-n Minimum non reference base count, defaults to 2\n"+
				"-q Minimum base quality, defaults to 13\n"+
				"-m Minimum alignment mapping quality, defaults to 13\n"+
				"-a Maximum non ref allele frequency, defaults to 0.03\n"+
				"-x Maximum number non ref bps in flanks, defaults to 1\n"+
				"-i Don't fail regions with an indel in the flanks\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/BamContextInspector -b Bam/rPENormal.bam\n"+
				"-r rPENormal_calls.bed -r Ref/human_g1k_v37_decoy.fasta -q 20 -m 20 -a 2 \n\n" +

				"**************************************************************************************\n");

	}
}
