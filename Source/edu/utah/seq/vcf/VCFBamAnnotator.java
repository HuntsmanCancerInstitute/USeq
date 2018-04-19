package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.parsers.SamAlignmentLoader;
import edu.utah.seq.parsers.SamSorter;
import edu.utah.seq.parsers.mpileup.MpileupTabixLoaderRs;
import edu.utah.seq.query.QueryIndexFileLoader;
import util.apps.MergeRegions;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

/** Beta for checking read 1 and read 2 bias
 * @author Nix*/
public class VCFBamAnnotator {
	
	//user defined fields
	private File vcfFile;
	private File[] bamFiles;
	private String[] sampleNames;
	private File vcfResultFile;
	private File bedFileToMerge = null;
	private File mergedBedFile = null;
	private File indexedFasta = null;
	private File samtools = null;
	private File bgzip = null;
	private File tabix = null;
	private int numberThreads = 0;
	private File tempDir = null;

	//internal
	private File mpileup = null;
	public static final Pattern END_POSITION = Pattern.compile(".*END=(\\d+).*", Pattern.CASE_INSENSITIVE);
	private int bpPadding = 250;
	private BufferedReader vcfIn = null;
	private int numVcfToProc = 0;
	private int minBaseQual = 20;
	private boolean verbose = true;
	private int numberShortRegions = 0;
	
	public VCFBamAnnotator(String[] args){
		try {	
			processArgs(args);
			
			//load vcf records and write out expanded bed for merging
			System.out.println("Parsing vcf file...");
			parseShortRegions();
			
			//merge them
			System.out.println("Merging regions...");
			String resName = Misc.removeExtension(vcfResultFile.getName());
			mergedBedFile = new File (tempDir, resName +".tempMerged.bed");
			new MergeRegions(new File[]{bedFileToMerge}, mergedBedFile);
			
			
			//pull alignments over these regions and split by R1 and R2
			System.out.println("Extracting R1 and R2 alignments over regions...");
			SamAlignmentLoader[] sals = new SamAlignmentLoader[bamFiles.length];
			File[] splitBamFiles = new File[bamFiles.length * 2];
			int counter = 0;
			for (int i=0; i< sals.length; i++) {
				sals[i] = new SamAlignmentLoader(bamFiles[i], mergedBedFile, tempDir, numberThreads);
				File[] r1r2 = sals[i].getParsedBams();
				splitBamFiles[counter++] = r1r2[0];
				splitBamFiles[counter++] = r1r2[1];
			}
			
			//make indexed mpileup file for the alignments that overlap these regions
			System.out.println("Building mpileup...");
			mpileup = new File (tempDir, resName +"_mpileup.gz");
			makeMPileup(mergedBedFile, mpileup, splitBamFiles, indexedFasta, samtools, bgzip, tabix, numberThreads);

			//open reader on vcfFile	
			vcfIn = IO.fetchBufferedReader(vcfFile);
			
			//load and process records
			loadAndParseMpileup();
			
			
			
			vcfIn.close();
			
			//write out merged multi sample vcf with DP, R1, R2, AF, NRAF
			
			//IO.deleteDirectory(tempDir);
			//IO.deleteDirectoryViaCmdLine(tempDir);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void loadAndParseMpileup() throws IOException {
		numVcfToProc = numberShortRegions/numberThreads;
		if (numVcfToProc < 25) numVcfToProc = 25;
		MpileupTabixLoaderRs[] loaders = new MpileupTabixLoaderRs[numberThreads];
		ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
		for (MpileupTabixLoaderRs l: loaders) {
			l = new MpileupTabixLoaderRs(mpileup, this);
			executor.execute(l);
		}
		executor.shutdown();
		while (!executor.isTerminated()) {}  //wait here until complete
		
		//check loaders 
		for (MpileupTabixLoaderRs l: loaders) if (l != null && l.isFailed()) throw new IOException("ERROR: Failed loading mpileup data!");
	}

	/**Loads the AL with vcf lines up to the numVcfToProc, if none could be read, returns false, otherwise true.*/
	public synchronized boolean loadVcfRecords(ArrayList<String> chunk) throws IOException{
		String record;
		int count = 0;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.startsWith("#") || record.length() == 0) continue;
			chunk.add(record);
			if (count++ > numVcfToProc) break;
		}
		
		if (chunk.size() == 0) return false;
		return true;
	}
	
	public synchronized void addProcessedVcfs(ArrayList<String> modVcfRecords, ArrayList<String> notAnnotated) {
		//System.out.println("NumProc "+modVcfRecords.size()+"    NotAnno "+notAnnotated.size());
	}
	
	public static void makeMPileup(File bedRegions, File tabixedMpileup, File[] bamFiles, File indexedFasta, File samtools, File bgzip, File tabix, int numberThreads) throws IOException {
		//tempDir
		String name = Misc.removeExtension(tabixedMpileup.getName());
		File tempDir = tabixedMpileup.getParentFile();
		tempDir.mkdirs();
		
		//mpileup
		File output = new File(tempDir, name);
		ArrayList<String> al = new ArrayList<String>();
		al.add(samtools.getCanonicalPath());
		al.add("mpileup -R -B -q 13 -d 100000 -f");
		al.add(indexedFasta.getCanonicalPath());
		al.add("-l"); al.add(bedRegions.getCanonicalPath());
		for (File b: bamFiles) al.add(b.getCanonicalPath());
		al.add(">");
		al.add(output.getCanonicalPath());
		al.add("\n");
		
		//bgzip
		al.add(bgzip.getCanonicalPath());
		al.add("--threads "+ numberThreads);
		al.add(output.getCanonicalPath());
		al.add("\n");
		
		//tabix
		al.add(tabix.getCanonicalPath());
		al.add("-s 1 -b 2 -e 2");
		al.add(output.getCanonicalPath()+".gz");
		al.add("\n");
		
		//execute the cmd
		String cmd = Misc.stringArrayListToString(al, " ");
		String[] shellOut = IO.executeShellScript(cmd, tempDir);
		
		//mpileup throws one line
		if (shellOut.length > 1) throw new IOException ("Failed to execute "+cmd+"\nerror, "+Misc.stringArrayToString(shellOut, "\n"));
		
		//look for files
		File gz = new File(output.getCanonicalPath()+".gz");
		File tbi = new File(output.getCanonicalPath()+".gz.tbi");
		if (tbi.exists() == false)  throw new IOException ("Failed to successfully execute "+cmd+"\nno tbi index file.");
		if (gz.exists() == false)  throw new IOException ("Failed to successfully execute "+cmd+"\nno bgzip results file.");
	}


	private void parseShortRegions() {
		try {
			bedFileToMerge = new File (tempDir, Misc.removeExtension(vcfResultFile.getName())+".tempToMerge.bed");
			PrintWriter out = new PrintWriter(new FileWriter(bedFileToMerge));
			BufferedReader in = IO.fetchBufferedReader(vcfFile);
			String line;
			String[] tokens;
			while ((line=in.readLine())!= null){
				if (line.startsWith("#") || line.trim().length() == 0) continue;
				tokens = Misc.TAB.split(line);
				//watch out for structural variants
				if (tokens[4].contains("<")) continue;
				int[] startStop = null;
				try {
					startStop = QueryIndexFileLoader.fetchEffectedBps(tokens, false); 
				} catch (Exception e){
					System.err.println("\tProblem parsing effected bps for the following variant, skipping!\n\t\t"+e.getMessage());
				}
				if (startStop != null){
					int start = startStop[0] - bpPadding;
					if (start < 0) start = 0;
					out.println(tokens[0]+"\t"+ start+ "\t"+ (startStop[1]+bpPadding));
					numberShortRegions++;
				}
			}
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error loading vcf files!");
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFBamAnnotator(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcfFile = new File(args[++i]); break;
					case 'b': bamFiles = IO.extractFiles(args[++i]); break;
					case 'n': sampleNames = Misc.COMMA.split(args[++i]); break;
					case 'r': vcfResultFile = new File(args[++i]); break;
					case 'f': indexedFasta = new File(args[++i]); break;
					case 's': samtools = new File(args[++i]); break;
					case 'g': bgzip = new File(args[++i]); break;
					case 't': tabix = new File(args[++i]); break;
					case 'p': numberThreads = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//check files
		if (vcfFile == null || vcfFile.exists() == false) Misc.printExit("\nError: please provide a vcf to annotate with count data.\n");
		if (bamFiles == null || bamFiles.length == 0) Misc.printExit("\nError: please provide a comma delimited list of bam files to use in extracting variant information.\n");
		if (vcfResultFile == null) Misc.printExit("\nError: please provide a file name to write the merged vcf records.\n");
		if (indexedFasta == null || indexedFasta.exists() == false) Misc.printExit("\nError: please provide an indexed fasta genome file.\n");
		if (samtools == null || samtools.exists() == false) Misc.printExit("\nError: please provide a path to your samtools application.\n");
		if (bgzip == null || bgzip.exists() == false) Misc.printExit("\nError: please provide a path to your bgzip application.\n");
		if (tabix == null || tabix.exists() == false) Misc.printExit("\nError: please provide a path to your tabix application.\n");
		
		if (sampleNames == null){
			sampleNames = new String[bamFiles.length];
			for (int i=0; i< sampleNames.length; i++){
				sampleNames[i] = Misc.removeExtension(bamFiles[i].getName());
			}
		}
		
		if (numberThreads == 0) numberThreads = Runtime.getRuntime().availableProcessors() - 1;
		if (numberThreads == 0) numberThreads = 1;
		tempDir = new File (vcfResultFile.getParentFile(), Misc.getRandomString(6)+"_TempDeleteMe");
		tempDir.mkdirs();
				
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Bam Annotator : Aug 2017                         **\n" +
				"**************************************************************************************\n" +
				"This app pulls read coverage information from each bam over the vcf variants and adds\n"+
				"the data as sample information in the vcf. Of particular use is the breakdown of\n"+
				"alignment depth by R1 and R2. \n"+

				"\nOptions:\n"+
				"-v Vcf file (xxx.vcf(.gz/.zip OK)) to annotate.\n"+
				"-b Comma delimited list of indexed bam files to extract alignment info.\n"+
				"-n Comma delimited list of 'Sample Names' corresponding to the bams, defaults to the\n"+
				"     bam file names.\n" +
				"-r Annotated vcf results file.\n"+
				"-f Indexed reference fasta used in aligning the bams.\n"+
				"-s Path to the samtools executable.\n"+
				"-g Path to the HTSlib bgzip executable.\n"+
				"-t Path to the HTSlib tabix executable.\n"+
				"-p Number of threads to use, defaults to all available.\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFConsensus -p illumina.vcf -q Strelka\n"+
				"-s stnd.indel.vcf.gz -t Scalpel -o indelCalls.vcf.gz \n\n"+

		"**************************************************************************************\n");

	}
	
	public static Bed makeBedFromVcfRecord(String vcfLine) throws Exception{
		String[] t = Misc.TAB.split(vcfLine);
		
			//fetch start and stop for effected bps.
			int[] startStop = QueryIndexFileLoader.fetchEffectedBps(t, false);
			if (startStop != null) {
				if (startStop[0] < 0) startStop[0] = 0;
				return new Bed(t[0], startStop[0], startStop[1], vcfLine, 0, '.');
			}
			
		return null;
	}

	public int getMinBaseQuality() {
		return minBaseQual;
	}

	public boolean isVerbose() {
		return verbose;
	}


	

}
