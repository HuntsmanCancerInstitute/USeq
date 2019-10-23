package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.FisherExact;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Takes output like the following and filters it. 
 * 
bcftools=/uufs/chpc.utah.edu/common/HIPAA/u0028003/BioApps/Bcftools/1.9/bin/bcftools
fasta=/uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa
regions=/scratch/mammoth/serial/u0028003/Underhill/SpikeBuild/VcfToInject/Bed/HSV1_GBM_IDT_Probes_Hg38Pad150bps_91K.bed.gz

$bcftools mpileup \
--count-orphans \
--no-BAQ \
--regions-file $regions \
--max-depth 100000 \
--max-idepth 100000 \
--gap-frac 0.001 \
--per-sample-mF \
--ignore-RG \
--min-MQ 13 \
--fasta-ref $fasta \
--ff UNMAP,SECONDARY,QCFAIL \
--annotate FORMAT/AD \
-Ou \
MergedTumor/15352X1_0.01.bam MergedNormal/15352X17_Hg38_final.bam | \
$bcftools norm \
--fasta-ref $fasta \
--multiallelics -any \
-Oz --threads 5 --output rawCallsNoBAQ.vcf.gz

 */
public class SimpleSomaticCaller {

	private File[] vcfFiles;

	//tumor
	private int minimumTumorReadDepth = 0;
	private int minimumAltReadDepth = 0;
	private double minimumTumorAltFraction = 0;
	private double maximumTumorAltFraction = 1;
	private boolean saveJustOneMultiAlt = false;

	//normal
	private double maximumNormalAltFraction = 1;
	private int minimumNormalReadDepth = 0;

	private double minimumTNFractionDiff = 0;
	private double minimumTNRatio = 0;
	private double maximumFisherPVal = 1;

	private String afInfo = "##INFO=<ID=T_AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=T_DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private String nafInfo = "##INFO=<ID=N_AF,Number=1,Type=Float,Description=\"Allele Frequency for normal\">";
	private String ndpInfo = "##INFO=<ID=N_DP,Number=1,Type=Integer,Description=\"Read depth for normal\">";
	private String fpvInfo = "##INFO=<ID=FPV,Number=1,Type=Float,Description=\"Left tailed Fisher's exact 2x2 ref alt tumor normal observation, -10Log10(p-val)\">";
	private boolean printInfo = true;

	private File saveDirectory = null;
	private Gzipper pass = null;
	private Gzipper fail = null;
	private int numPass = 0;
	private int numFail = 0;
	private boolean saveFail = false;
	private int vcfIndex = 0;

	private HashMap<String, Integer> formatIndexes = null;
	private int maxCountSum = 10000;
	private FisherExact fisher = new FisherExact(maxCountSum);
	private String userArgs = null;
	private BufferedReader in = null;
	private File workingVcf = null;
	private double workingAFDiff = 0;
	private String workingChromPosition = null;
	private String workingRecord = null;

	public SimpleSomaticCaller (String[] args) { 

		try {
			processArgs(args);

			printArgs();

			//Parse files or std in?
			IO.pl("\nName\tPassing\tFailing");
			if (vcfFiles != null) parseVcfFiles();
			else parseVcfStream();
		}
		catch (Exception e) {
			if (vcfFiles != null) IO.el("\nProblem parsing "+workingVcf);
			e.printStackTrace();
		}
		finally {
			try {
				if (pass != null) pass.close();
				if (fail != null) fail.close();
				if (in != null) in.close();
			} catch (IOException e) {}	
		}
		IO.pl("\nComplete!");
	}

	private void parseVcfStream() throws Exception {
		//save dir?
		if (saveDirectory == null) {
			File workingDir = new File(System.getProperty("user.dir")).getCanonicalFile();
			saveDirectory = new File(workingDir, "SSC");
			saveDirectory.mkdirs();
		}
		System.out.print("StandardIn\t");
		File passGz = new File (saveDirectory, "stdIn.ssc.vcf.gz");
		pass = new Gzipper(passGz);
		if (saveFail) {
			File failGz = new File (saveDirectory, "stdIn.ssc.fail.vcf.gz");
			fail = new Gzipper(failGz);
		}
		in = new BufferedReader(new InputStreamReader(System.in));
		parse();
		
		//save last?
		if (workingRecord !=null) {
			numPass++;
			pass.println(workingRecord);
		}

		IO.pl(numPass+"\t"+numFail);
	}

	private void parseVcfFiles() throws Exception {
		for (File vcf: vcfFiles){
			workingVcf = vcf;
			System.out.print(workingVcf.getName()+"\t");
			File passGz = new File (saveDirectory, Misc.removeExtension(vcf.getName())+".ssc.vcf.gz");
			pass = new Gzipper(passGz);
			if (saveFail) {
				File failGz = new File (saveDirectory, Misc.removeExtension(vcf.getName())+".ssc.fail.vcf.gz");
				fail = new Gzipper(failGz);
			}
			in = IO.fetchBufferedReader(vcf);
			parse();
			//save last?
			if (workingRecord !=null) {
				numPass++;
				pass.println(workingRecord);
			}
			IO.pl(numPass+"\t"+numFail);
		}
	}

	private void printArgs() {
		IO.pl("Thresholds for Tumor and Normal:");

		IO.pl("  -u "+minimumTumorReadDepth+"\tMin T alignment depth");
		IO.pl("  -a "+minimumAltReadDepth+"\tMin T alt count");
		IO.pl("  -t "+minimumTumorAltFraction+"\tMin T allelic fraction");
		IO.pl("  -x "+maximumTumorAltFraction+"\tMax T allelic fraction");

		IO.pl("  -o "+minimumNormalReadDepth+"\tMin N alignment depth");
		IO.pl("  -n "+maximumNormalAltFraction+"\tMax N allelic fraction");

		IO.pl("  -d "+minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		IO.pl("  -r "+minimumTNRatio+"\tMin T/N allelic fraction ratio");
		IO.pl("  -f "+maximumFisherPVal+"\tMax Fisher's one-sided pvalue");
		IO.pl("  -y "+ saveFail+ "\tSave failing records to file");
		IO.pl("  -m "+ saveJustOneMultiAlt+ "\tSave just one multi alt allele with largest tAF-nAF");
	}

	public void parse() throws Exception{

		String line = null;
		numPass = 0;
		numFail = 0;
		vcfIndex = 0;
		try {
			while ((line = in.readLine())!= null) {
				if (line.trim().length() == 0) continue;
				if (line.startsWith("##INFO") && printInfo) {
					addHeader(pass);
					pass.println(line);
					if (fail!=null) fail.println(line);
					printInfo = false;
				}
				else if (line.startsWith("#")) {
					pass.println(line);
					if (fail!=null) fail.println(line);
				}
				else parseRecord(line);
			}
		} catch (IOException e) {
			e.printStackTrace();
			throw new Exception ("\nProblem parsing line "+line);
		}
	}



	private void addHeader(Gzipper out) throws IOException {
		out.println("##SimpleSomaticCaller "+userArgs);
		out.println(afInfo);
		out.println(dpInfo);
		out.println(nafInfo);
		out.println(ndpInfo);
		out.println(fpvInfo);
	}

	private void parseRecord(String line) throws IOException {
		//#CHROM	POS	   ID	REF	ALT	QUAL FILTER	INFO	      FORMAT	MergedTumor/15352X1_0.01.bam	MergedNormal/15352X17_Hg38_final.bam
		//chr1	11199396	.	AG	A	0	  .	    INDEL;IDV=75  PL:AD	     0,255,239:5977,59	            0,255,237:13580,12
		//  0      1        2   3   4   5     6       7             8            9                          10
		String[] fields = Misc.TAB.split(line);
		if (formatIndexes == null) loadFormatIndexes(Misc.COLON.split(fields[8]));

		//parse and calculate info
		String[] tSample = Misc.COLON.split(fields[9]);
		double[] tRefAlt = parseAD(tSample);
		double tDP = tRefAlt[1]+tRefAlt[0];
		double tAF = tRefAlt[1]/tDP;
		String[] nSample = Misc.COLON.split(fields[10]);
		double[] nRefAlt = parseAD(nSample);
		double nDP = nRefAlt[1]+nRefAlt[0];
		double nAF = nRefAlt[1]/nDP;
		double afDiff = tAF-nAF;
		double afRto = tAF/nAF;
		int totalSum = (int)(nDP+ tDP);
		if (totalSum > this.maxCountSum) fisher = new FisherExact(totalSum);
		double pval = fisher.getLeftTailedP((int)tRefAlt[0],(int)tRefAlt[1], (int)nRefAlt[0], (int)nRefAlt[1]);

		ArrayList<String> flags = new ArrayList<String>();
		//pass min alt count?
		if (tRefAlt[1] < minimumAltReadDepth) flags.add("minAltReads="+tRefAlt[1]);

		//pass min max tumor AF
		if (tAF < minimumTumorAltFraction || tAF > maximumTumorAltFraction) flags.add("tAF="+tAF);

		//pass min tumor DP
		if (tDP < minimumTumorReadDepth) flags.add("tDP="+tDP);

		//pass max normal AF
		if (nAF > maximumNormalAltFraction) flags.add("nAF="+nAF);

		//pass min normal DP
		if (nDP < minimumNormalReadDepth) flags.add("nDP="+nDP);

		//pass af diff

		if (afDiff < minimumTNFractionDiff) flags.add("TNDiff="+afDiff);
			
		//pass af ratio
		if (afRto < minimumTNRatio) flags.add("TNRatio="+afRto);

		//calc fisher's left side pval
		if (pval > maximumFisherPVal) flags.add("pVal="+pval);
			
		//fail?
		if (flags.size() !=0) {
			fields[6] = Misc.stringArrayListToString(flags, ";");
			printFail(buildNewRecord(tAF,tDP,nAF,nDP,pval,fields, vcfIndex++));
		}
		
		else printPass(buildNewRecord(tAF,tDP,nAF,nDP,pval,fields, vcfIndex++), afDiff,fields);
		
	}
	



	private void printPass(String currentRecord, double diffAF, String[] fields) throws IOException {
		if (saveJustOneMultiAlt) {
			String key = fields[0]+"_"+fields[1];
			if (workingChromPosition == null) {
				workingChromPosition = key;
				workingAFDiff = diffAF;
				workingRecord = currentRecord;
			}
			//same position
			else if (workingChromPosition.equals(key)){
				//check scores
				if (diffAF> workingAFDiff) {
					//fail working and reset
					printFail(workingRecord);
					workingAFDiff = diffAF;
					workingRecord = currentRecord;
				}
				else {
					//fail current
					printFail(currentRecord);
				}
			}
			//diff position
			else {
				//print old
				numPass++;
				pass.println(workingRecord);
				//set current
				workingChromPosition = key;
				workingAFDiff = diffAF;
				workingRecord = currentRecord;
			}
		}
		//save everything
		else {
			numPass++;
			pass.println(currentRecord);
		}
	}

	private void printFail(String record) throws IOException {
		//only save if they want it
		if (fail!=null) fail.println(record);
		numFail++;
	}

	private String buildNewRecord(double tAF, double tDP, double nAF, double nDP, double pval, String[] fields, int index) {
		//#CHROM	POS	   ID	REF	ALT	QUAL FILTER	INFO	      FORMAT	MergedTumor/15352X1_0.01.bam	MergedNormal/15352X17_Hg38_final.bam
		//chr1	11199396	.	AG	A	0	  .	    INDEL;IDV=75  PL:AD	     0,255,239:5977,59	            0,255,237:13580,12
		//  0      1        2   3   4   5     6       7             8            9                          10
		StringBuilder sb = new StringBuilder();
		//chrom and pos
		sb.append(fields[0]); sb.append("\t");
		sb.append(fields[1]); sb.append("\t");
		//id
		sb.append("SSC_"); sb.append(index); sb.append("\t");
		//ref alt qual filter
		sb.append(fields[3]); sb.append("\t");
		sb.append(fields[4]); sb.append("\t");
		sb.append(fields[5]); sb.append("\t");
		sb.append(fields[6]); sb.append("\t");
		//info
		sb.append("T_AF=");
		sb.append(Num.formatNumber(tAF, 5)); sb.append(";");
		sb.append("T_DP=");
		sb.append((int)tDP); sb.append(";");
		sb.append("N_AF=");
		sb.append(Num.formatNumber(nAF, 5)); sb.append(";");
		sb.append("N_DP=");
		sb.append((int)nDP); sb.append(";");
		sb.append("FPV=");
		//sb.append(Num.formatNumber(pval, 5));sb.append(";");
		double phred = Num.minus10log10(pval);
		sb.append(Num.formatNumber(phred, 3));sb.append(";");
		sb.append(fields[7]);
		//samples
		for (int i=8; i< fields.length; i++) {
			sb.append("\t");
			sb.append(fields[i]);
		}
		return sb.toString();
	}

	private double[] parseAD(String[] sample) {
		String[] refAlt = Misc.COMMA.split(sample[formatIndexes.get("AD")]);
		return new double[] {Double.parseDouble(refAlt[0]), Double.parseDouble(refAlt[1])};
	}

	private void loadFormatIndexes(String[] f) {
		formatIndexes = new HashMap<String, Integer>();
		for (int i=0; i< f.length; i++)formatIndexes.put(f[i], i);
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SimpleSomaticCaller(args);
	}	



	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'x': maximumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'a': minimumAltReadDepth = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
					case 'f': maximumFisherPVal = Double.parseDouble(args[++i]); break;
					case 'y': saveFail = true; break;
					case 'm': saveJustOneMultiAlt = true; break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'z': printBash(); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		userArgs = IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ");
		IO.pl("\n"+userArgs +"\n");

		//pull vcf files
		if (forExtraction != null ) {
			File[][] tot = new File[3][];
			tot[0] = IO.extractFiles(forExtraction, ".vcf");
			tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
			tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
			vcfFiles = IO.collapseFileArray(tot);
			if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");
			if (saveDirectory == null) saveDirectory = vcfFiles[0].getParentFile();
		}
		if (saveDirectory != null) saveDirectory.mkdirs();

	}	
	
	public static void printBash() {
		IO.pl("\n# Be sure to use bcftools 1.9 and extensively test any changes to the param values\n"+
				"\n"+
	            "tumorBam=0.005.bam\n"+
				"normalBam=normal.bam\n"+
				"vcf=0.005.raw.vcf.gz\n"+
				"\n"+
				"bcftools=/BioApps/Bcftools1.9/bin/bcftools\n"+
				"fasta=/Indexes/H38IndexForBwa-0.7.17/hs38DH.fa\n"+
				"regions=/Bed/IDT_Probes_Hg38Pad150bps_91K.bed\n"+
				"\n"+
				"$bcftools mpileup --count-orphans --no-BAQ --regions-file $regions --max-depth 100000 --max-idepth 100000 \n"+
				"--gap-frac 0.001 --per-sample-mF --ignore-RG --min-MQ 13 --fasta-ref $fasta --ff UNMAP,SECONDARY,QCFAIL \n"+
				"--annotate FORMAT/AD -Ou $tumorBam $normalBam | \n"+
				"$bcftools norm --fasta-ref $fasta --multiallelics -any --threads 5 | grep -w -v '<\\*>' | gzip > $vcf\n");
		System.exit(0);
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                          Simple Somatic Caller: Oct 2019                         **\n" +
				"**************************************************************************************\n" +
				"Takes vcf output from Bcftools mpileup and norm applications (http://www.htslib.org)\n"+
				"run on a paired tumor and normal bam file set and filters the variants for somatic\n"+
				"calls. Select -z to print an example bash script. A one-tailed Fisher's Exact test is\n"+
				"performed on each alt allele count to select somatic variants.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-s Directory in which to save the parsed files, defaults to parent dir of the vcfs.\n"+
				"-t Minimum tumor allele frequency (AF)\n"+
				"-x Maximum tumor AF\n"+
				"-n Maximum normal AF\n"+
				"-u Minimum tumor alignment depth (DP)\n"+
				"-a Minimum tumor alt count\n"+
				"-o Minimum normal DP\n"+
				"-d Minimum T-N AF difference\n"+
				"-r Minimum T/N AF ratio\n"+
				"-f Maximum Fisher's TN pvalue, defaults to 1, no threshold\n"+
				"-y Save failing records\n"+
				"-m For multi alt alleles, save just the biggest T-N AF difference.\n"+
				"-z Print example bcftools bash script.\n"+

				"\nExample: java -jar pathToUSeq/Apps/SimpleSomaticCaller -v VCFFiles/ -f 0.001 \n"+
				"    -t 0.005 -x 0.2 -n 0.1 -u 100 -a 3 -o 100 -d 0.0025 -r 3 -m \n\n"+


				"**************************************************************************************\n");

	}

}
