package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.data.Graphite;
import util.gen.FisherExact;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Tool for rescoring variants using Graphite for calculating actual counts. This catches some complex indels. */
public class VCFGraphiteFilter {

	//user fields
	private File vcf;
	private File tumorBam;
	private File normalBam;
	private File graphiteExe;
	private File outputDir;
	private File fasta;
	private boolean verbose = false;

	//pre filters
	private float minQualFilt = 0;
	private String filterTxtToExclude = null;

	//post filters
	//tumor
	private int minimumTumorReadDepth = 0;
	private int minimumAltReadDepth = 0;
	private double minimumTumorAltFraction = 0;
	private double maximumTumorAltFraction = 1;

	//normal
	private double maximumNormalAltFraction = 1;
	private int minimumNormalReadDepth = 0;
	private double minimumTNFractionDiff = 0;
	private double minimumTNRatio = 0;
	private double maximumFisherPVal = 1;

	//internal 
	private String userArgs = null;
	private int maxCountSum = 10000;
	private FisherExact fisher = new FisherExact(maxCountSum);
	private Gzipper outPass = null;
	private Gzipper outFail = null;
	private VCFParser vcfParser = null;
	private int numVcf;
	private int numPostPreFiltVcf;
	private int numPassVcf;
	private int numFailVcf;
	private Graphite graphite = null;
	private File graphiteTempDir = null;
	private String tAfInfo = "##INFO=<ID=T_AF,Number=1,Type=Float,Description=\"Graphite allele frequency for tumor\">";
	private String tDpInfo = "##INFO=<ID=T_DP,Number=1,Type=Integer,Description=\"Graphite read depth for tumor\">";
	private String nAfInfo = "##INFO=<ID=N_AF,Number=1,Type=Float,Description=\"Graphite allele frequency for normal\">";
	private String nDpInfo = "##INFO=<ID=N_DP,Number=1,Type=Integer,Description=\"Graphite read depth for normal\">";
	private String fpvInfo = "##INFO=<ID=FPV,Number=1,Type=Float,Description=\"Left tailed Fisher's exact 2x2 ref alt tumor normal observation, -10Log10(p-val)\">";
	private Pattern iT_AF = Pattern.compile("T_AF=[\\d/.]+;");
	private Pattern iT_DP = Pattern.compile("T_DP=\\d+;");
	private Pattern iN_AF = Pattern.compile("N_AF=[\\d/.]+;");
	private Pattern iN_DP = Pattern.compile("N_DP=\\d+;");
	private Pattern iFPV = Pattern.compile("FPV=[\\d/.]+;");

	//constructor
	public VCFGraphiteFilter(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			parseVcf();

			scoreCalls();

			IO.deleteDirectory(graphiteTempDir);
			

			IO.pl("\nNumInputRecords\tAfterPreFiltering\tPass\tFail\tFileName:");
			IO.pl(numVcf+"\t"+numPostPreFiltVcf+"\t"+numPassVcf+"\t"+numFailVcf+"\t"+vcf);

			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
			
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem comparing vcf replicas.");
		}
	}


	private void scoreCalls() throws IOException {
		IO.pl("Extracting AF and DP with Graphite...");
		//Fetch AF and DP for A vars in tumorBam
		if (verbose) IO.pl("\tTumor Bam");
		HashMap<String, float[]> keyAfDpTumor = graphite.annotate(vcf, tumorBam);
		if (verbose) IO.pl("\tNormal Bam");
		HashMap<String, float[]> keyAfDpNormal = graphite.annotate(vcf, normalBam);

		IO.pl("Scoring Vcf Records...");
		scoreVcf(keyAfDpTumor, keyAfDpNormal);

		outPass.close();
		outFail.close();
	}

	private void addHeader() throws IOException {
		//strip out any T_AF, T_DP, N_AF, N_DP
		ArrayList<String> filtered = new ArrayList<String>();
		for (String s: vcfParser.getStringComments()) {
			if (s.contains("INFO=<ID=T_AF")) continue;
			if (s.contains("INFO=<ID=T_DP")) continue;
			if (s.contains("INFO=<ID=N_AF")) continue;
			if (s.contains("INFO=<ID=N_DP")) continue;
			if (s.contains("INFO=<ID=FPV")) continue;
			filtered.add(s);
		}
		boolean printed = false;
		String[] infos = new String[] {tAfInfo, tDpInfo, nAfInfo, nDpInfo, fpvInfo};
		for (String s: filtered) {
			if (printed == false && s.startsWith("##INFO=")) {
				outPass.println("##VCFGraphiteFilter= "+userArgs);
				outFail.println("##VCFGraphiteFilter= "+userArgs);
				printed = true;
				for (String i: infos) {
					outPass.println(i);
					outFail.println(i);
				}
			}
			outPass.println(s);
			outFail.println(s);
		}
	}

	private void parseVcf() throws FileNotFoundException, IOException {
		//create writers and add headers
		String name = Misc.removeExtension(vcf.getName());
		File pass = new File(outputDir,name+".pass.vcf.gz");
		File fail = new File(outputDir,name+".fail.vcf.gz");
		outPass = new Gzipper(pass);
		outFail = new Gzipper(fail);
		
		IO.pl("\nParsing and pre filtering vcf records...");
		vcfParser = new VCFParser(vcf, true, false, false);
		numVcf = vcfParser.getVcfRecords().length;
		addHeader();
		ArrayList<VCFRecord> failing = new ArrayList<VCFRecord>();
		if (minQualFilt !=0.0f) failing.addAll(vcfParser.filterVCFRecordsOnQUAL(minQualFilt));
		if (filterTxtToExclude != null) failing.addAll(vcfParser.filterVCFRecordsOnFILTER(filterTxtToExclude));
		numPostPreFiltVcf = vcfParser.getVcfRecords().length;
		
		//any filtered out?
		if (failing.size()>0) for (VCFRecord v: failing) outFail.println(v.getOriginalRecord());
	}

	private void scoreVcf(HashMap<String, float[]> keyAfDpTumor, HashMap<String, float[]> keyAfDpNormal) throws IOException {

		//for each variant
		for (VCFRecord v: vcfParser.getVcfRecords()) {
			if (verbose) IO.pl(v.getOriginalRecord());
			//#CHROM POS ID REF ALT 0 1 2 3 4
			String key = v.getChrPosRefAlt(false);
			//AF, DP, fRef, rRef, fAlt, rAlt
			//0   1    2      3     4     5
			float[] afDpTumor =  keyAfDpTumor.get(key);
			float[] afDpNormal =  keyAfDpNormal.get(key);
			if (afDpTumor == null || afDpNormal == null) throw new IOException("Failed to find the AF and DP from the the tumor or normal graphite results, see "+v.getOriginalRecord());

			//calc pval
			int totalSum = (int)(afDpNormal[1]+ afDpTumor[1]);
			if (totalSum > maxCountSum) {
				fisher = new FisherExact(totalSum);
				maxCountSum = totalSum;
			}
			int tRef = (int)(afDpTumor[2]+ afDpTumor[3]);
			int tAlt = (int)(afDpTumor[4]+ afDpTumor[5]);
			int nRef = (int)(afDpNormal[2]+ afDpNormal[3]);
			int nAlt = (int)(afDpNormal[4]+ afDpNormal[5]);
			double pval = fisher.getLeftTailedP(tRef,tAlt, nRef, nAlt);
			double phred = Num.minus10log10(pval);

			//add in new info replacing old if it exists
			String vcfCorr = addInfo(afDpTumor, afDpNormal, Num.formatNumber(phred, 3), v);
			if (verbose) IO.pl(vcfCorr);

			//maximumFisherPVal
			if(pval > maximumFisherPVal) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing pval "+pval);
			}
			//minimumAltReadDepth
			else if (tAlt < minimumAltReadDepth) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing minimumAltReadDepth "+tAlt);
			}
			//minimumTumorReadDepth
			else if ((int)afDpTumor[1]< minimumTumorReadDepth ) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing minimumTumorReadDepth "+(int)afDpTumor[1]);
			}
			//minimumNormalReadDepth
			else if ((int)afDpNormal[1]< minimumNormalReadDepth ) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing minimumTumorReadDepth "+((int)afDpNormal[1]));
			}
			//minimumTumorAltFraction
			else if (afDpTumor[0] < minimumTumorAltFraction) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing minimumTumorAltFraction "+afDpTumor[0]);
			}
			//maximumTumorAltFraction
			else if (afDpTumor[0] > maximumTumorAltFraction) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing maximumTumorAltFraction "+afDpTumor[0]);
			}
			//maximumNormalAltFraction
			else if (afDpNormal[0]> maximumNormalAltFraction ) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing maximumNormalAltFraction "+afDpNormal[0]);
			}
			//minimumTNFractionDiff
			else if ((afDpTumor[0] - afDpNormal[0]) < minimumTNFractionDiff ) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing minimumTNFractionDiff "+(afDpTumor[0] - afDpNormal[0]));
			}
			//minimumTNRatio
			else if ((afDpTumor[0] / afDpNormal[0]) < minimumTNRatio ) {
				outFail.println(vcfCorr);
				numFailVcf++;
				if (verbose) IO.pl("\tFailing minimumTNRatio "+(afDpTumor[0] / afDpNormal[0]));
			}
			else {
				outPass.println(vcfCorr);
				numPassVcf++;
				if (verbose) IO.pl("\tPASS");
			}
		}
	}

	private String addInfo(float[] afDpTumor, float[] afDpNormal, String pval, VCFRecord v) {
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
		String[] fields = Misc.TAB.split(v.getOriginalRecord());
		String info = fields[7];

		String tAF = "T_AF="+formatAf(afDpTumor[0])+";";
		String tDP = "T_DP="+(int)afDpTumor[1]+";";
		String nAF = "N_AF="+formatAf(afDpNormal[0])+";";
		String nDP = "N_DP="+(int)afDpNormal[1]+";";
		String fpv = "FPV="+pval+";";

		info = swap(iT_AF, info, tAF);
		info = swap(iT_DP, info, tDP);
		info = swap(iN_AF, info, nAF);
		info = swap(iN_DP, info, nDP);
		info = swap(iFPV, info, fpv);

		fields[7] = info;
		return Misc.stringArrayToString(fields, "\t");
	}

	private String formatAf(double af){
		if (af == 0.0) return "0";
		return Num.formatNumberJustMax(af, 4);
	}

	public static String swap(Pattern pat, String info, String replacement) {
		Matcher mat = pat.matcher(info);
		return mat.replaceFirst(replacement);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFGraphiteFilter(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		Pattern pat = Pattern.compile("-[a-z0-9]");
		userArgs = IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ");
		IO.pl("\n"+userArgs +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcf = new File(args[++i]); break;
					case 'b': tumorBam = new File(args[++i]); break;
					case 'm': normalBam = new File(args[++i]); break;
					case 'g': graphiteExe = new File(args[++i]); break;
					case 's': outputDir = new File(args[++i]); break;
					case 'i': fasta = new File(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'x': maximumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'a': minimumAltReadDepth = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
					case 'f': maximumFisherPVal = Double.parseDouble(args[++i]); break;
					case 'e': filterTxtToExclude = args[++i]; break;
					case 'q': minQualFilt = Float.parseFloat(args[++i]); break;
					case 'z': verbose = true; break;
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
		if (vcf == null || vcf.exists()==false) Misc.printErrAndExit("\nError -v: please provide a path to the vcf ("+vcf+") file.");
		if (tumorBam == null || normalBam == null || tumorBam.exists()==false || normalBam.exists()==false) Misc.printErrAndExit("\nError -t or -n: please provide paths to the tumor or normal bam files.");
		if (graphiteExe == null || graphiteExe.canExecute() == false) Misc.printErrAndExit("\nError -g: please provide path to the graphite executable.");
		if (outputDir == null) Misc.printErrAndExit("\nError -o: please provide a path to a directory to write the adjudicated vcf files.\n");
		if (fasta == null || fasta.exists()==false) Misc.printErrAndExit("\nError -i: please provide an indexed fasta file.\n");
		if (outputDir.exists() == false) outputDir.mkdirs();

		//make temp dir for Graphite
		int numThreads = Runtime.getRuntime().availableProcessors() - 1;
		graphiteTempDir = new File (outputDir, "GraphiteTempDir_"+Misc.getRandomString(6));
		graphiteTempDir = graphiteTempDir.getCanonicalFile();
		graphiteTempDir.mkdirs();
		graphite = new Graphite(graphiteExe, fasta, graphiteTempDir, numThreads);

		printArgs();
	}	
	
	private void printArgs() {
		IO.pl("Parameters:");
		IO.pl(" -v Vcf "+vcf);
		IO.pl(" -b TumorBam "+tumorBam);
		IO.pl(" -m NormalBam "+normalBam);
		IO.pl(" -s OutputDir\t"+outputDir);
		IO.pl(" -i FastaIndex "+fasta);
		IO.pl(" -g Graphite app "+graphiteExe);
		IO.pl(" -q "+minQualFilt+ "\tMin QUAL prefilter ");
		IO.pl(" -e "+filterTxtToExclude+ "\tFILTER exclude prefilter ");
		
		IO.pl(" -u "+minimumTumorReadDepth+"\tMin T alignment depth");
		IO.pl(" -a "+minimumAltReadDepth+"\tMin T alt count");
		IO.pl(" -t "+minimumTumorAltFraction+"\tMin T allelic fraction");
		IO.pl(" -x "+maximumTumorAltFraction+"\tMax T allelic fraction");
		IO.pl(" -o "+minimumNormalReadDepth+"\tMin N alignment depth");
		IO.pl(" -n "+maximumNormalAltFraction+"\tMax N allelic fraction");
		IO.pl(" -d "+minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		IO.pl(" -r "+minimumTNRatio+"\tMin T/N allelic fraction ratio");
		IO.pl(" -f "+maximumFisherPVal+"\tMax Fisher's one-sided pvalue");
		IO.pl(" -z "+verbose+"\tVerbose debugging output");
	}
	
	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Graphite Filter : May 2020                       **\n" +
				"**************************************************************************************\n" +
				"Uses Graphite (https://github.com/dillonl/graphite) to accurately calculate tumor and\n"+
				"normal alignment support and apply a set of filters to split a vcf file into passing\n"+
				"and failing. Apply pre filters to reduce the size of the vcf records to limit the slow\n"+
				"re scoring process. Be sure to run the USeq MergePairedAlignments app on the bam\n"+
				"files, otherwise Graphite will double count variant observations from overlapping\n"+
				"pairs. Uses all threads available.\n"+

				"\nRequired Options:\n"+	
				"-v Path to a vcf file xxx.vcf(.gz/.zip OK) to rescore.\n" +
				"-b Path to the merged pair alignment tumor bam (with bai index in the same dir).\n"+
				"-m Path to the merged pair normal bam\n"+
				"-s Path to a director to save the split vcf files.\n"+
				"-i Fasta index used in aligning the reads\n"+
				"-g Graphite executable, see https://github.com/dillonl/graphite install e.g.:\n"+
				"     module load gcc/4.9.2\n" + 
				"     git clone https://github.com/dillonl/graphite.git\n" + 
				"     mkdir graphite/bin\n" + 
				"     cd graphite/bin\n" + 
				"     cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++)\n" + 
				"     make\n"+
				
				"\nPre filters:\n"+
				"-q Min QUAL score for vcf record import, defaults to 0\n"+
				"-e FILTER field value to exclude records from import, defaults to no filtering\n"+
				
				"\nPost filters:\n"+
				"-t Minimum tumor allele frequency (AF)\n"+
				"-x Maximum tumor AF\n"+
				"-n Maximum normal AF\n"+
				"-u Minimum tumor alignment depth (DP)\n"+
				"-a Minimum tumor alt count\n"+
				"-o Minimum normal DP\n"+
				"-d Minimum T-N AF difference\n"+
				"-r Minimum T/N AF ratio\n"+
				"-f Maximum Fisher's TN pvalue, defaults to 1, no threshold\n"+
				"-z Verbose debugging output\n"+

				"\nExample: module load gcc/4.9.2; java -Xmx25G -jar USeq/Apps/VCFGraphiteFilter\n" +
				"-v ssc.vcf.gz -b mergedTumor.bam -m mergedNormal.bam -s VGF/ -i ~/Ref/hg38.fa \n"+
				"-g ~/Apps/graphite -q 4 -e BKAF -f 0.001 -t 0.005 -n 0.1 -u 100 -a 4 -o 100 -d 0.0045 \n"+
				"-r 3\n\n"+

				"**************************************************************************************\n");

	}
}