package edu.utah.seq.parsers.jpileup;

import java.io.*;
import java.util.ArrayList;
import java.util.regex.*;
import util.gen.*;

/** 
 * Merged multiple BamPileups into one file.

# MinMapQual	13
# MinBaseQual	10
# IncludeOverlappingBpCounts false
# Bed /uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner/Bed/AvatarV2IDTBedNov2019/hg38IdtProbesPad150bp.bed.gz
# Bam/Cram 0 /scratch/general/pe-nfs1/u0028003/Avatar/AJobs/1101000/Alignment/1101000_NormalDNA/Bam/1101000_NormalDNA_Hg38_final.bam
# Chr	1BasePos	Ref	A,C,G,T,N,Del,Ins,FailBQ
chr1	68964	G	0,0,1,0,0,0,0,0
chr1	68965	C	0,1,0,0,0,0,0,0
chr1	68966	A	1,0,0,0,0,0,0,0

 * @author Nix
 * */
public class BamPileupMerger {

	//user defined fields
	private File[] pileups;
	private File mergedPileup;
	private File tabix;
	private File bgzip;
	private PrintWriter out;
	private BufferedReader[] ins;
	private boolean skipHeaderCheck = false;
	

	//constructor
	public BamPileupMerger(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
		
			buildHeader();
			
			mergeBPs();
			
			
			//close IO
			out.close();
			for (int i=0; i< pileups.length; i++) ins[i].close();
			
			compressAndIndex(bgzip, tabix, mergedPileup, new String[] {"-f", "-s", "1", "-b", "2", "-e", "2"}, true);
			
			IO.pl("\tMerged pileup saved to "+mergedPileup+".gz");
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes");
			
		} catch (Exception e) {
			if (out !=null) out.close();
			if (ins !=null) for (int i=0; i< pileups.length; i++) IO.closeNoException(ins[i]);
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR making bpileup! Correct and restart.");
		}
	}

	public static void compressAndIndex(File bgzip, File tabix, File uncompressedPileupTxtFile, String[] tabixCmds, boolean verbose) throws IOException {
		if (verbose) IO.pl("Bgzipping...");
		//compress with bgzip
		String[] cmd = { bgzip.getCanonicalPath(), "-f", "--threads", new Integer(Runtime.getRuntime().availableProcessors()).toString(), uncompressedPileupTxtFile.getCanonicalPath()};
		String[] output = IO.executeCommandLineReturnAll(cmd);
		File bgzipped = new File(uncompressedPileupTxtFile.toString()+".gz");
		if (output == null || output.length != 0 || bgzipped.exists() == false){
			throw new IOException("\nFailed to bgzip compress "+uncompressedPileupTxtFile+"\nError: "+Misc.stringArrayToString(output, "\n"));
		}

		//tabix
		if (verbose) IO.pl("Tabixing...");
		
		cmd = new String[tabixCmds.length+2];
		cmd[0]= tabix.getCanonicalPath();
		cmd[cmd.length-1] = bgzipped.getCanonicalPath();
		int i = 1;
		for (String c: tabixCmds)cmd[i++] = c;
		//cmd = new String[]{ tabix.getCanonicalPath(), "-f", "-s", "1", "-b", "2", "-e", "2", bgzipped.getCanonicalPath() };
		
		output = IO.executeCommandLineReturnAll(cmd);
		File indexed = new File (bgzipped+".tbi");
		if (output == null || output.length != 0 || indexed.exists() == false){
			throw new IOException("\nFailed to tabix index "+bgzipped+"\nError: "+Misc.stringArrayToString(output, "\n"));
		}
	}

	private void mergeBPs() throws IOException {
		IO.pl("Merging bam pileup files...");
		
		//all files should be of the same length
		String[] f = null;
		String[] s = null;
		String fLine = null;
		String sLine = null;
		StringBuilder sb = null;
		boolean countsFound;
		long linesWithCounts = 0;
		long totalLines = 0;
		
		//read in first bp file lines
		while ((fLine = ins[0].readLine()) != null) {
			totalLines++;
			sb = new StringBuilder();
			countsFound = false;
			//chr1	68964	G	0,0,1,0,0,0,0,0
			//load first
			f = Misc.TAB.split(fLine);
			sb.append(fLine);
			if (f[3].equals("0,0,0,0,0,0,0,0") == false) countsFound = true;
			
			//remainder
			for (int i=1; i< ins.length; i++) {
				sLine = ins[i].readLine();
				s = Misc.TAB.split(sLine);
				if ( (s[0].equals(f[0])==false) || (s[1].equals(f[1])==false)) throw new IOException("BamPileup file lines don't match?!\n\t"+fLine+"\n\t"+sLine);
				sb.append("\t");
				sb.append(s[3]);
				if (s[3].equals("0,0,0,0,0,0,0,0") == false) countsFound = true;
			}
			
			//save it?
			if(countsFound) {
				out.println(sb.toString());
				linesWithCounts++;
			}
		}
		
		//check that remainder are all null
		for (int i=1; i< ins.length; i++) {
			sLine = ins[i].readLine();
			if (sLine != null) throw new IOException("This BamPileup file has additional lines. "+pileups[i]);
		}
		
		IO.pl("\t"+totalLines+" bps parsed, "+linesWithCounts+" with counts.");
		
	}

	private void buildHeader() {
		try {
			if (skipHeaderCheck == false) IO.pl("Checking bam pileup file headers...");
			
			//load headers
			BpHeader[] headers = new BpHeader[ins.length];
			for (int i=0; i<headers.length; i++) headers[i] = new BpHeader(ins[i], pileups[i]);
			
			//check headers
			for (int i=1; i<headers.length; i++) headers[0].compare(headers[i]);
			
			//write out first header
			out.println(headers[0].minMap);
			out.println(headers[0].minBase);
			out.println(headers[0].overlap);
			out.println(headers[0].bed);
			for (int i=0; i<headers.length; i++) out.println("# BamCram "+i+" "+headers[i].bam);
			out.println(headers[0].chrom);
		
		} catch (IOException e) {
			out.close();
			for (int i=0; i< pileups.length; i++) IO.closeNoException(ins[i]);
			Misc.printErrAndExit("\nERROR: headers in bamPileup files differ.");
		}
	}
	
	private class BpHeader {
		String chrom = null;
		String bed = null;
		String minMap = null;
		String minBase = null;
		String overlap = null;
		String bam= null;
		String printAll= null;

		private BpHeader (BufferedReader in, File bp) throws IOException {
			String line = null;
			for (int i=0; i< 1000; i++) {
				line = in.readLine();
				if (line.startsWith("# Chr")) {
					chrom = line;
					break;
				}
				if (line.startsWith("# MinMapQual")) minMap = line;
				else if (line.startsWith("# MinBaseQual")) minBase = line;
				else if (line.startsWith("# IncludeOverlappingBpCounts")) overlap = line;
				else if (line.startsWith("# PrintAll")) {
					printAll = line;
					if (printAll.contains("true") == false) throw new IOException("\nLooks like the PrintAll flag was not true. "
							+ "BamPileupMerged can only merge single bam/cram bamPileup files. Re run the BamPileup app on just one bam or cram file.");
				}
				else if (line.startsWith("# Bed")) bed = line;
				else if (line.startsWith("# Bam")) {
					if (bam != null) throw new IOException("\nFound more than one bam files in "+bp);
					String[] f = Misc.WHITESPACE.split(line);
					bam = (f[f.length-1]);
				}
			}
			if (chrom==null || bed==null || minMap==null || minBase==null || overlap==null || bam==null || printAll==null) throw new IOException ("\nFailed to find a proper bp header in "+bp);
		}

		public void compare(BpHeader other) throws IOException {
			if (chrom.equals(other.chrom) == false) throw new IOException("Header chrom lines differ "+other.chrom);
			if (skipHeaderCheck) return;
			if (bed.equals(other.bed) == false) throw new IOException("Header bed lines differ "+other.bed);
			if (minMap.equals(other.minMap) == false) throw new IOException("Header minMap lines differ "+other.minMap);
			if (minBase.equals(other.minBase) == false) throw new IOException("Header minBase lines differ "+other.minBase);
			if (overlap.equals(other.overlap) == false) throw new IOException("Header overlap lines differ "+other.overlap);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamPileupMerger(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args) throws Exception{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File tabixBinDirectory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': pileups = IO.extractFiles(new File(args[++i]), ".gz"); break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
					case 'p': mergedPileup = new File(args[++i]); break;
					case 's': skipHeaderCheck = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull tabix and bgzip
		if (tabixBinDirectory == null) Misc.printExit("\nError: please point to your HTSlib directory containing the tabix and bgzip executables (e.g. ~/BioApps/HTSlib/1.10.2/bin/ )\n");
		bgzip = new File (tabixBinDirectory, "bgzip");
		tabix = new File (tabixBinDirectory, "tabix");
		//look for bgzip and tabix executables
		if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+bgzip+" "+tabix);

		//check pileups and results
		if (pileups == null || pileups.length <2) Misc.printErrAndExit("\nPlease provide one or more xxx.gz bamPileup files to merge.\n");
		
		//set results file
		if (mergedPileup == null || mergedPileup.getName().endsWith(".gz")== false) Misc.printErrAndExit("\nPlease provide a path to a file ending with xxx.gz to save the merged pileup file, e.g. bp.txt.gz\n");
		String path = mergedPileup.getCanonicalPath();
		File noGz = new File (path.substring(0, path.length()-3));
		mergedPileup = noGz;
		
		//make IO
		out = new PrintWriter(new BufferedWriter(new FileWriter(mergedPileup)));
		ins = new BufferedReader[pileups.length];
		for (int i=0; i< pileups.length; i++) ins[i] = IO.fetchBufferedReader(pileups[i]);
		
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Bam Pileup Meger:  Sept 2021                       **\n" +
				"**************************************************************************************\n" +
				"BPM merges BamPileup files, compresses, and indexes them. \n"+
				
				"\nOptions:\n"+
				"-d Path to a directory of single sample bgzipped BamPileup files with the tbi indexes.\n"+
				"-p Path to save the output pileup file, must end in gz.\n"+
				"-t Path to the directory containing the compiled bgzip and tabix executables. See\n" +
				"     https://www.htslib.org\n"+
				"-s Skip header check, use with caution!\n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/BamPileupMerger -d BPFiles/\n"+
				"     -t ~/BioApps/HTSlib/1.10.2/bin -p abl.bp.txt.gz\n\n" +

				"**************************************************************************************\n");
	}

}
