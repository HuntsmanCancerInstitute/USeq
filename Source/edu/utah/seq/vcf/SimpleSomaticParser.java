package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Takes output like the following and filters it. Be sure to use bcftools 1.9, 1.10 won't work, not sure about subsequent releases, 1.10 contains bugs.
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
public class SimpleSomaticParser {

	//user defined
	private File[] vcfFiles;
	private File vcfOut;
	private int minimumReadDepth = 12;
	private int minimumAltObservations = 2;
	private double minimumAltFraction = 0.01;

	//internal
	private Gzipper pass = null;
	private int numPass = 0;
	private int numFail = 0;
	private int numSkipped = 0;
	private HashMap<String, Integer> formatIndexes = null;
	private int fieldCount = 0;
	private BufferedReader in = null;
	private boolean printHeader = true;

	public SimpleSomaticParser (String[] args) { 

		try {
			processArgs(args);
			pass = new Gzipper(vcfOut);

			for (File vcf: vcfFiles) {
				IO.pl("Parsing "+vcf.getName());
				parse(vcf);
				IO.pl("\t#Pass: "+numPass+"  #Fail: "+numFail);
				printHeader = false;
			}


		}
		catch (Exception e) {
			vcfOut.delete();
			e.printStackTrace();
			System.exit(1);
		}
		finally {
			try {
				if (pass != null) pass.close();
			} catch (IOException e) {}	
		}
		IO.pl("\nComplete!");
	}


	public void parse(File vcf) throws Exception{

		in = IO.fetchBufferedReader(vcf);
		numPass = 0;
		numFail = 0;
		String line = null;
		try {
			while ((line = in.readLine())!= null) {
				if (line.trim().length() == 0) continue;
				else if (line.startsWith("#")) {
					if (printHeader) pass.println(line);
					else continue;
				}
				else parseRecord(line);
			}
		} catch (IOException e) {
			e.printStackTrace();
			throw new Exception ("\nProblem parsing line "+line);
		}
		in.close();
	}


	private void parseRecord(String line) throws Exception {
		//#CHROM	POS	   ID	REF	ALT	QUAL FILTER	INFO	      FORMAT	21310X10_Hg38_final.bam  21310X1_Hg38_final.bam      ...
		//chr1	11199396	.	AG	A	0	  .	    INDEL;IDV=75  PL:AD	       0,255,239:5977,59	       0,255,237:13580,12
		//  0      1        2   3   4   5     6       7             8                9                          10                11...
		// AD  reference counts: alt counts
		// need to potentially match multiple tabs since bcltools is inserting blank cells
		String[] fields = Misc.TAB.split(line);
		if (formatIndexes == null) loadFormatIndexes(Misc.COLON.split(fields[8]));

		//check field count
		if (fieldCount == 0) fieldCount = fields.length;
		else if (fieldCount != fields.length) {
			IO.el("\tMissing fields, skipping:\t"+line);
			numSkipped++;
			return;
		}

		//for each sample
		for (int i=9; i< fields.length; i++) {
			String[] sample = Misc.COLON.split(fields[i]);
			double[] refAlt = parseAD(sample);
			if (passThresholds(refAlt)) {
				numPass++;
				pass.println(line);
				return;
			}
		}
		numFail++;
	}

	private boolean passThresholds(double[] refAlt) {

		//pass min alt count?
		if (refAlt[1] < minimumAltObservations) return false;

		//pass DP
		double dp = refAlt[1]+refAlt[0];
		if (dp < minimumReadDepth) return false;

		//pass AF
		double af = refAlt[1]/dp;
		if (af < minimumAltFraction ) return false;

		return true;
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
		new SimpleSomaticParser(args);
	}	



	/**This method will process each argument and assign new variables
	 * @throws IOException 
	 * */
	public void processArgs(String[] args) throws IOException{
		String useqVersion = IO.fetchUSeqVersion();
		String source = useqVersion+" Args: USeq/SimpleSomaticParser "+ Misc.stringArrayToString(args, " ");
		IO.pl("\n"+ source +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': forExtraction = new File(args[++i]); break;
					case 'o': vcfOut = new File(args[++i]); break;
					case 'a': minimumAltFraction = Double.parseDouble(args[++i]); break;
					case 'd': minimumReadDepth = Integer.parseInt(args[++i]); break;
					case 't': minimumAltObservations = Integer.parseInt(args[++i]); break;
					case 'b': printBash(); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//check IO
		if (vcfOut == null) Misc.printExit("\nError: please provide a file path xxx.vcf.gz to write out the parsed variants.\n");

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");
	}	

	public static void printBash() {
		IO.pl("\n# Be sure to use bcftools 1.9 and extensively test any changes to the param values\n"+
				"\n"+
				"bcftools=/BioApps/Bcftools1.9/bin/bcftools\n"+
				"fasta=/Indexes/H38IndexForBwa-0.7.17/hs38DH.fa\n"+
				"regions=/Bed/IDT_Probes_Hg38Pad150bps_91K.bed\n"+
				"\n"+
				"$bcftools mpileup --count-orphans --no-BAQ --regions-file $regions --max-depth 100000 --max-idepth 100000 \n"+
				"--gap-frac 0.001 --per-sample-mF --ignore-RG --min-MQ 13 --fasta-ref $fasta --ff UNMAP,SECONDARY,QCFAIL \n"+
				"--annotate FORMAT/AD -Ou *.bam | \n"+
				"$bcftools norm --fasta-ref $fasta --multiallelics -any --threads 5 | grep -w -v '<\\*>' | gzip > raw.vcf.gz\n");
		System.exit(0);
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                          Simple Somatic Parser: Oct 2023                         **\n" +
				"**************************************************************************************\n" +
				"Takes vcf output from Bcftools mpileup and norm applications (http://www.htslib.org)\n"+
				"run on one or more bam files and filters the variants for those meeting several hard\n"+
				"thresholds regarding allele fraction, read depth, and number of alt observations.\n"+
				"Only one of the samples must pass these criteria to save the entire vcf record.\n"+

				"\nParams:\n"+
				"-i Path to the bcftools xxx.vcf(.gz/.zip OK) file.\n" +
				"-o Path to write out the passing variants, xxx.vcf.gz\n"+
				"-a Minimum allele fraction (AF), defaults to 0.01\n"+
				"-d Minimum alignment depth (DP), defaults to 12\n"+
				"-t Minimum alt count, defaults to 2\n"+
				"-b Print example bcftools bash script.\n"+

				"\nExample: java -jar pathToUSeq/Apps/SimpleSomaticParser -i bcfTools.vcf.gz -a 0.025 \n"+
				" -u 100 -t 3 -o passing.vcf.gz\n\n"+


				"**************************************************************************************\n");

	}

}
