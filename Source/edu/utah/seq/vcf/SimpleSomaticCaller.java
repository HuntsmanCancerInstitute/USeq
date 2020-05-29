package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.jpileup.BamPileupTabixLoaderSingle;
import edu.utah.seq.parsers.jpileup.BaseCount;
import edu.utah.seq.parsers.jpileup.BpileupLine;
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

	private File vcfFile;
	private String fileName;
	private File normalBamPileup = null;
	private int normalBamPileupIndex = 0;
	private int bpPadding = 3;

	//tumor
	private int minimumTumorReadDepth = 0;
	private int minimumAltReadDepth = 0;
	private double minimumTumorAltFraction = 0;
	private double maximumTumorAltFraction = 1;
	private boolean saveJustOneMultiAlt = false;

	//normal
	private double maximumNormalAltFraction = 1;
	private int minimumNormalReadDepth = 0;
	private double maximumNormalIndelAF = 0.01;
	private double maximumNormalIndelAFRatio = 0.5;

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

	private HashMap<String, Integer> formatIndexes = null;
	private int maxCountSum = 10000;
	private FisherExact fisher = new FisherExact(maxCountSum);
	private String userArgs = null;
	private BufferedReader in = null;
	private String workingChromPosition = null;
	private double[] workingScores = null;
	private String[] workingFields = null;
	private int workingIndex = 0;
	private BamPileupTabixLoaderSingle bamPileupLoader = null;

	public SimpleSomaticCaller (String[] args) { 

		try {
			processArgs(args);

			printArgs();

			//Parse files or std in?
			IO.pl("\nName\tPassing\tFailing");
			if (vcfFile != null) parseVcfFile();
			else parseVcfStream();
			
			IO.pl(numPass+"\t"+numFail);
			
		}
		catch (Exception e) {
			if (vcfFile != null) IO.el("\nProblem parsing "+vcfFile);
			e.printStackTrace();
		}
		finally {
			try {
				if (pass != null) pass.close();
				if (fail != null) fail.close();
				if (in != null) in.close();
				if (normalBamPileup != null) bamPileupLoader.getTabixReader().close();
			} catch (IOException e) {}	
		}
		IO.pl("\nComplete!");
	}

	private void parseVcfStream() throws Exception {
		fileName = "stdIn";
		//save dir?
		if (saveDirectory == null) {
			File workingDir = new File(System.getProperty("user.dir")).getCanonicalFile();
			saveDirectory = new File(workingDir, "SSC_stdIn");
			saveDirectory.mkdirs();
		}
		IO.p("StandardIn\t");
		File passGz = new File (saveDirectory, "stdIn.ssc.pass.vcf.gz");
		pass = new Gzipper(passGz);
		File failGz = new File (saveDirectory, "stdIn.ssc.fail.vcf.gz");
		fail = new Gzipper(failGz);
		in = new BufferedReader(new InputStreamReader(System.in));
		parse();
	}

	private void parseVcfFile() throws Exception {
		fileName = Misc.removeExtension(vcfFile.getName());
		IO.p(vcfFile.getName()+"\t");
		File passGz = new File (saveDirectory, fileName+ ".ssc.pass.vcf.gz");
		pass = new Gzipper(passGz);
		File failGz = new File (saveDirectory, fileName+ ".ssc.fail.vcf.gz");
		fail = new Gzipper(failGz);
		in = IO.fetchBufferedReader(vcfFile);
		parse();
	}

	private void printArgs() {
		IO.pl("Thresholds for Tumor and Normal:");

		IO.pl("  -u "+minimumTumorReadDepth+"\tMin T alignment depth");
		IO.pl("  -a "+minimumAltReadDepth+"\tMin T alt count");
		IO.pl("  -t "+minimumTumorAltFraction+"\tMin T allelic fraction");
		IO.pl("  -x "+maximumTumorAltFraction+"\tMax T allelic fraction");

		IO.pl("\n  -o "+minimumNormalReadDepth+"\tMin N alignment depth");
		IO.pl("  -n "+maximumNormalAltFraction+"\tMax N allelic fraction");

		IO.pl("\n  -d "+minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		IO.pl("  -r "+minimumTNRatio+"\tMin T/N allelic fraction ratio");
		IO.pl("  -f "+maximumFisherPVal+"\tMax Fisher's one-sided pvalue");
		IO.pl("  -m "+ saveJustOneMultiAlt+ "\tSave just one multi alt allele with largest tAF-nAF");
		
		if (normalBamPileup != null) {
			IO.pl("\n  -p "+ bpPadding+ "\tBP padding for scanning normal alignments for INDEL background artifacts");
			IO.pl("  -i "+ maximumNormalIndelAF+ "\tMaximum normal INDEL AF");
			IO.pl("  -j "+ maximumNormalIndelAFRatio+ "\tMaximum normal INDEL AF/ tumor AF");
			IO.pl("  -c "+ normalBamPileupIndex+"\tNormal BamPileup sample index");
		}
	}

	public void parse() throws Exception{

		String line = null;
		numPass = 0;
		numFail = 0;
		int index = 0;
		try {
			while ((line = in.readLine())!= null) {
				if (line.trim().length() == 0) continue;
				if (line.startsWith("##INFO") && printInfo) {
					addHeader(pass);
					addHeader(fail);
					pass.println(line);
					fail.println(line);
					printInfo = false;
				}
				else if (line.startsWith("#")) {
					pass.println(line);
					fail.println(line);
				}
				else {
					index++;
					parseRecord(line, index);
				}
			}
			//save last one?
			if (saveJustOneMultiAlt) {
				numPass++;
				pass.println(buildNewRecord(workingFields, workingScores, workingIndex));
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			throw new Exception ("\nProblem parsing line "+line);
		}
	}

	private void addHeader(Gzipper out) throws IOException {
		out.println("##SimpleSomaticCaller= "+userArgs);
		out.println(afInfo);
		out.println(dpInfo);
		out.println(nafInfo);
		out.println(ndpInfo);
		out.println(fpvInfo);
	}

	private void parseRecord(String line, int index) throws Exception {
		//#CHROM	POS	   ID	REF	ALT	QUAL FILTER	INFO	      FORMAT	MergedTumor/15352X1_0.01.bam	MergedNormal/15352X17_Hg38_final.bam
		//chr1	11199396	.	AG	A	0	  .	    INDEL;IDV=75  PL:AD	     0,255,239:5977,59	            0,255,237:13580,12
		//  0      1        2   3   4   5     6       7             8            9                          10
		String[] fields = Misc.TAB.split(line);
		if (formatIndexes == null) loadFormatIndexes(Misc.COLON.split(fields[8]));

		//parse and calculate info
		String[] tSample = Misc.COLON.split(fields[9]);
		double[] tRefAlt = parseAD(tSample);
		String[] nSample = Misc.COLON.split(fields[10]);
		double[] nRefAlt = parseAD(nSample);	
		
		//tDP, tAF, nDP, nAF, afDiff, afRto, pval
		//0     1    2    3      4      5     6
		ThresholdCheck tc = passThresholds(tRefAlt, nRefAlt);
		
		//did it fail a threshold?
		if (tc.filter != null) printFail(fields, tc.scores, index, tc.filter);
		
		else {
			double[] scores = tc.scores;
			// do they want to do an indel scan?
			if (normalBamPileup != null) {
				
				//calculate the max normal indel AF and check thresholds
				double maxIndelAF = fetchMaxIndelAF(fields);
				double testAF = maximumNormalIndelAFRatio * scores[1];
				
				if (maxIndelAF >= maximumNormalIndelAF || maxIndelAF > testAF) printFail(fields, scores, index, "indelBkg_"+nF(maxIndelAF));
				else printPass(fields, scores, index);	
				
			}
			// nope, just pass it
			else printPass(fields, scores, index);
		}
	}
	
	private static String nF(double d) {
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(4);
		return f.format(d);
	}
	
	private class ThresholdCheck {
		String filter = null;
		double[] scores = null;
		private ThresholdCheck (String filter, double[] scores) {
			this.filter = filter;
			this.scores = scores;
		}
	}
	
	private ThresholdCheck passThresholds(double[] tRefAlt, double[] nRefAlt) {
		ArrayList<String> flags = new ArrayList<String>();
		
		//pass min alt count?
		if (tRefAlt[1] < minimumAltReadDepth) flags.add("minAltDp_"+ (int)tRefAlt[1]);

		//pass min max tumor AF
		double tDP = tRefAlt[1]+tRefAlt[0];
		double tAF = tRefAlt[1]/tDP;
		if (tAF < minimumTumorAltFraction || tAF > maximumTumorAltFraction) flags.add("tAF_"+ nF(tAF));
		
		//pass min tumor DP
		if (tDP < minimumTumorReadDepth) flags.add("tDP_"+ (int)tDP);

		//pass max normal AF
		double nDP = nRefAlt[1]+nRefAlt[0];
		double nAF = nRefAlt[1]/nDP;
		if (nAF > maximumNormalAltFraction)  flags.add("nAF_"+ nF(nAF));

		//pass min normal DP
		if (nDP < minimumNormalReadDepth)  flags.add("nDP_"+ (int)nDP);

		//pass af diff
		double afDiff = tAF-nAF;
		if (afDiff < minimumTNFractionDiff)  flags.add("afDiff_"+ nF(afDiff));
			
		//pass af ratio
		double afRto = tAF/nAF;
		if (afRto < minimumTNRatio)  flags.add("afRto_"+ nF(afRto));

		//calc fisher's left side pval
		int totalSum = (int)(nDP+ tDP);
		if (totalSum > maxCountSum) {
			fisher = new FisherExact(totalSum);
			maxCountSum = totalSum;
		}
		double pval = fisher.getLeftTailedP((int)tRefAlt[0],(int)tRefAlt[1], (int)nRefAlt[0], (int)nRefAlt[1]);

		if (pval > maximumFisherPVal)  flags.add("pVal_"+ nF(pval));
		
		String cFlag = null;
		if (flags.size() != 0) cFlag = Misc.stringArrayListToString(flags, ";");
		return new ThresholdCheck(cFlag, new double[] {tDP, tAF, nDP, nAF, afDiff, afRto, pval});
	}

	private double fetchMaxIndelAF(String[] fields) throws Exception {
		ArrayList<BpileupLine> bpl = bamPileupLoader.fetchBpileupRecords(fields);
		double maxIndelAf = -1;
		for (BpileupLine bp: bpl) {
			BaseCount bc = bp.getSamples()[normalBamPileupIndex];
			double iAf = bc.getAlleleFreqIndel();
			if (iAf > maxIndelAf) {
				maxIndelAf = iAf;
			}
		}
		return maxIndelAf;
	}

	private void printPass(String[] fields, double[] scores, int index) throws IOException {
		//tDP, tAF, nDP, nAF, afDiff, afRto, pval
		//0     1    2    3      4      5     6
		
		//save one multi alt
		if (saveJustOneMultiAlt) {
			//				chrom		pos
			String key = fields[0]+"_"+fields[1];
			
			//first record?
			if (workingChromPosition == null) {
				workingChromPosition = key;
				workingScores = scores;
				workingFields = fields;
				workingIndex = index;
			}
			
			//same position
			else if (workingChromPosition.equals(key)){
				//check scores
				if (workingScores[4] < scores[4]) {
					//fail working and reset
					printFail(workingFields, workingScores, workingIndex, "multiAlt_"+workingScores[4]);
					workingFields = fields;
					workingScores = scores;
					workingIndex = index;
				}
				else {
					//fail current
					printFail(fields, scores, index, "multiAlt_"+scores[4]);
				}
			}
			//diff position
			else {
				//print old
				numPass++;
				pass.println(buildNewRecord(workingFields, workingScores, workingIndex));
				//set current
				workingChromPosition = key;
				workingFields = fields;
				workingScores = scores;
				workingIndex = index;
			}
		}
		
		//save everything
		else {
			numPass++;
			pass.println(buildNewRecord(fields, scores, index));
		}
	}

	private void printFail(String[] fields, double[] scores, int index, String filter) throws IOException {
		//modify fields[6] the filter field
		if (fields[6].length()==0 || fields[6].equals(".")) fields[6] = filter;
		else fields[6] = filter+";"+fields[6];
		fail.println( buildNewRecord(fields, scores, index) );
		numFail++;
	}

	private String buildNewRecord(String[] fields, double[] scores, int index) {
		//tDP, tAF, nDP, nAF, afDiff, afRto, pval
		//0     1    2    3      4      5     6
		double tAF =scores[1];
		double tDP =scores[0];
		double nAF =scores[3];
		double nDP =scores[2];
		double pval=scores[6];
		
		//#CHROM	POS	   ID	REF	ALT	QUAL FILTER	INFO	      FORMAT	MergedTumor/15352X1_0.01.bam	MergedNormal/15352X17_Hg38_final.bam
		//chr1	11199396	.	AG	A	0	  .	    INDEL;IDV=75  PL:AD	     0,255,239:5977,59	            0,255,237:13580,12
		//  0      1        2   3   4   5     6       7             8            9                          10
		StringBuilder sb = new StringBuilder();
		//chrom and pos
		sb.append(fields[0]); sb.append("\t");
		sb.append(fields[1]); sb.append("\t");
		//id
		sb.append("SSC_"); sb.append(fileName); sb.append("_"); sb.append(index); sb.append("\t");
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



	/**This method will process each argument and assign new variables
	 * @throws IOException 
	 * */
	public void processArgs(String[] args) throws IOException{
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcfFile = new File(args[++i]); break;
					case 'b': normalBamPileup = new File(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'x': maximumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'a': minimumAltReadDepth = Integer.parseInt(args[++i]); break;
					case 'p': bpPadding = Integer.parseInt(args[++i]); break;
					case 'c': normalBamPileupIndex = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
					case 'f': maximumFisherPVal = Double.parseDouble(args[++i]); break;
					case 'i': maximumNormalIndelAF = Double.parseDouble(args[++i]); break;
					case 'j': maximumNormalIndelAFRatio = Double.parseDouble(args[++i]); break;
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

		//any vcf file or are they streaming?
		if (vcfFile != null) {
			if (vcfFile.canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file!\n");
			if (saveDirectory == null) {
				String name = Misc.removeExtension(vcfFile.getName());
				saveDirectory = new File (vcfFile.getParentFile(), "SSC_"+name);
			}
			saveDirectory.mkdirs();
		}
		if (normalBamPileup != null) bamPileupLoader = new BamPileupTabixLoaderSingle(normalBamPileup, bpPadding);

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
				"**                          Simple Somatic Caller: May 2020                         **\n" +
				"**************************************************************************************\n" +
				"Takes vcf output from Bcftools mpileup and norm applications (http://www.htslib.org)\n"+
				"run on a paired tumor and normal bam file set and filters the variants for somatic\n"+
				"calls. Select -z to print an example bash script. A one-tailed Fisher's Exact test is\n"+
				"performed on each alt allele count to select somatic variants. Vt normalization isn't\n"+
				"needed for the output vcfs. An indel scan is performed to down weight snvs and indels\n"+
				" that occur in noisy indel regions.\n"+

				"\nRequired Params:\n"+
				"-v Path to the bcftools xxx.vcf(.gz/.zip OK) file.\n" +
				"-b Normal BamPileup file, compressed and tabix indexed, see the USeq BamPileup app.\n"+
				
				"\nDefault Params:\n"+
				"-s Directory in which to save the parsed file, defaults to parent dir of the vcf.\n"+
				"-c Normal BamPileup index, 0\n"+
				"-t Minimum tumor allele frequency (AF), 0\n"+
				"-x Maximum tumor AF, 1\n"+
				"-n Maximum normal AF, 1\n"+
				"-u Minimum tumor alignment depth (DP), 0\n"+
				"-a Minimum tumor alt count, 0\n"+
				"-o Minimum normal DP, 0\n"+
				"-i Maximum normal INDEL AF, 0\n"+
				"-j Maximum normal INDEL AF/ tumor AF ratio, 0\n"+
				"-d Minimum T-N AF difference, 0\n"+
				"-r Minimum T/N AF ratio, 0\n"+
				"-f Maximum Fisher's TN pvalue, defaults to 1\n"+
				"-m For multi alt alleles, save just the biggest T-N AF difference\n"+
				"-p BP padding for indel scan, default 3\n"+
				"-z Print example bcftools bash script.\n"+

				"\nExample: java -jar pathToUSeq/Apps/SimpleSomaticCaller -v bcfTools.vcf.gz -f 0.001 \n"+
				" -t 0.005 -x 0.2 -n 0.1 -u 100 -a 3 -o 100 -d 0.0025 -r 3 -m -s SSC/ -b norm.bp.txt.gz\n"+
				"-i 0.01 -j 0.5\n\n"+


				"**************************************************************************************\n");

	}

}
