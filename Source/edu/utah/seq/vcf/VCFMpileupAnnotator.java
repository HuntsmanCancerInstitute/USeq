package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.parsers.mpileup.MpileupTabixLoaderAFDP;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Given a vcf file compares each record against a tabix indexed mpileup file and adds/ overwrites AF and DP info.
 * @author Nix*/
public class VCFMpileupAnnotator {
	
	//user defined fields
	private File vcfFile;
	private File mpileup;
	private File vcfOut;
	private int minBaseQuality = 0;
	private boolean verbose = false;
	private int numberThreads = 0;
	private int numDecimals = 4;
	
	//internal
	private static final int numVcfToProc = 100;
	private MpileupTabixLoaderAFDP[] loaders;

	//working
	private BufferedReader vcfIn;
	private int numRecords = 0;
	private ArrayList<String> vcfHeader = new ArrayList<String>();
	private ArrayList<String> vcfRecords = new ArrayList<String>();
	
	//constructor
	public VCFMpileupAnnotator(String[] args){
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);
			
			//make Loaders
			loaders = new MpileupTabixLoaderAFDP[numberThreads];
			for (int i=0; i< loaders.length; i++) loaders[i] = new MpileupTabixLoaderAFDP(mpileup, this);

			System.out.print("\nParsing "+vcfFile.getName()+"\t");

			createReaderSaveHeader();

			//load mpileup data
			ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
			for (MpileupTabixLoaderAFDP l: loaders) executor.execute(l);
			executor.shutdown();
			while (!executor.isTerminated()) {}

			//check loaders 
			for (MpileupTabixLoaderAFDP l: loaders) {
				if (l.isFailed()) throw new IOException("ERROR: File Loader issue! \n");
			}
			vcfIn.close();

			VcfSorter[] finalVcfs = fetchVcfSorterArray();
			Arrays.sort(finalVcfs);
			writeOutModifiedVcf (finalVcfs, vcfOut);

			//print summary stats
			System.out.println(numRecords);


		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			//shut down loaders
			for (MpileupTabixLoaderAFDP m : loaders) m.getTabixReader().close();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
		}
	}

	private void writeOutModifiedVcf(VcfSorter[] finalVcfs, File file) throws IOException {
		Gzipper out = new Gzipper(file);
		//add header
		for (String s: vcfHeader) out.println(s);
		for (VcfSorter v : finalVcfs) out.println(v.fields, "\t");
		out.close();
	}
	
	private VcfSorter[] fetchVcfSorterArray() throws Exception {
		int num= vcfRecords.size();

		VcfSorter[] v = new VcfSorter[num];
		for (int x = 0; x < num; x++){
			String vcf = vcfRecords.get(x);
			v[x] = new VcfSorter(Misc.TAB.split(vcf));
		}
		return v;
	}
	
	private class VcfSorter implements Comparable<VcfSorter>{
		private String chr;
		private int pos;
		//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES
		String[] fields;
		
		public VcfSorter(String[] fields){
			this.fields = fields;
			chr = fields[0];
			pos = Integer.parseInt(fields[1]);
		}
		
		public int compareTo(VcfSorter other) {
			//sort by chromosome
			int x = chr.compareTo(other.chr);
			if (x !=0) return x;
			//sort by position
			if (pos < other.pos) return -1;
			if (pos > other.pos) return 1;
			return 0;
		}
	}

	public synchronized void save(ArrayList<String> vcf) {
		vcfRecords.addAll(vcf);
	}

	/**Loads the AL with vcf lines up to the numVcfToProc, if none could be read, returns false, otherwise true.*/
	public synchronized boolean loadVcfRecords(ArrayList<String> chunk) throws IOException{
		String record;
		int count = 0;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			chunk.add(record);
			if (count++ > numVcfToProc) break;
		}
		
		if (chunk.size() == 0) return false;
		numRecords += chunk.size();
		return true;
	}
	
	private void createReaderSaveHeader() throws Exception {
		vcfIn = IO.fetchBufferedReader(vcfFile);
		vcfHeader.clear();
		String record;
		boolean addInfo = true;
		boolean endFound = false;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			if (record.startsWith("#")){
				//add new info?
				if (addInfo && record.startsWith("##INFO=")) {
					vcfHeader.add("##INFO=<ID=AF,Number=1,Type=String,Description=\"Comma delimited list of allele frequencies for each sample in the mpileup data.\">");
					vcfHeader.add("##INFO=<ID=DP,Number=1,Type=String,Description=\"Comma delimited list of read depths for each sample in the mpileup data.\">");
					addInfo = false;
				}
				//skip prior AF and DP infos
				if (record.startsWith("##INFO=<ID=AF,") || record.startsWith("##INFO=<ID=DP,")){}
				else vcfHeader.add(record);
				
				//at the end? exit
				if (record.startsWith("#CHROM")){
					endFound = true;
					break;
				}
			}
			else {
				endFound = true;
				break;
			}
		}
		if (endFound == false) Misc.printErrAndExit("\nERROR: failed to find the #CHROM line, aborting\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFMpileupAnnotator(args);
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
					case 'o': vcfOut = new File(args[++i]); break;
					case 'm': mpileup = new File(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'e': numDecimals = Integer.parseInt(args[++i]); break;
					case 'd': verbose = true; break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull vcf files
		if (vcfFile == null || vcfFile.canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");
		if (vcfOut == null) Misc.printExit("\nError: please provide a path to a xxx.vcf.gz file to write out the gzipped results!\n");

		//tabix indexed mpileup?
		if (mpileup == null || mpileup.canRead() == false) Misc.printExit("\nError: please provide a path to a bgzipped tabix indexed mpileup file.\n");
		File index = new File (mpileup.toString()+".tbi");
		if (index.exists() == false) Misc.printExit("\nError: cannot find the '"+index.getName()+"' index file corresponding to this indexed mpileup file "+mpileup.getName());	
		
		//threads to use
		double gigaBytesAvailable = ((double)Runtime.getRuntime().maxMemory())/ 1073741824.0;
		int numPossCores = (int)Math.round(gigaBytesAvailable/5.0);
		if (numPossCores < 1) numPossCores = 1;
		int numPossThreads = Runtime.getRuntime().availableProcessors();
		if (numberThreads == 0){
			if (numPossCores <= numPossThreads) numberThreads = numPossCores;
			else numberThreads = numPossThreads;
			System.out.println("Core usage:\n\tGB available to Java:\t"+ Num.formatNumber(gigaBytesAvailable, 1));
			System.out.println("\tAvailable cores:\t"+numPossThreads);
			System.out.println("\tNumber cores to use @ 5GB/core:\t"+numberThreads);
		}
	}	


	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          VCF Mpileup Annotator : April 2019                      **\n" +
				"**************************************************************************************\n" +
				"VMA estimates the AF and DP of a vcf record from a single sample mpileup file.  It \n"+
				"replaces the AF or DP INFO values in the vcf records if present. For INDELs, the\n"+
				"region effected is scanned and the maximum AF and DP assigned. Provide the max memory\n"+
				"available to the app to maximize cpu usage.\n"+

				"\nRequired:\n"+
				"-v Path to a xxx.vcf(.gz/.zip OK) file.\n" +
				"-m Path to a bgzip compressed and tabix indexed single sample mpileup file. e.g.:\n"+
				"      1) Mpileup: 'samtools mpileup -B -R -A -d 1000000 -f $fastaIndex -l\n"+ 
				"         $bedFile $indexedBamFile > bam.mpileup'\n"+
				"      2) Bgzip: 'tabix-0.2.6/bgzip bam.mpileup'\n"+
                "         Tabix: 'tabix-0.2.6/tabix -s 1 -b 2 -e 2 bam.mpileup.gz'\n"+
				"-o Path to a xxx.vcf.gz file to save the modified vcf records.\n"+
						
				"\nOptional:\n" +
				"-q Minimum mpileup sample bp quality, defaults to 0\n"+
				"-e Number of decimals in the AF, defaults to 4\n"+
				"-d Print verbose debugging output.\n" +
				"-t Number of threads to use, defaults to all/ 5GB.\n"+
				

				"\n"+

				"Example: java -Xmx64G -jar pathTo/USeq/Apps/VCFMpileupAnnotator -v spikes.vcf -q 13\n"+
				"-m bam.mpileup.gz -o spikes.mod.vcf.gz\n\n"+

		        "**************************************************************************************\n");
	}

	//getters and setters

	public int getMinBaseQuality() {
		return minBaseQuality;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public int getNumDecimals() {
		return numDecimals;
	}
}
