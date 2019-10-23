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
import util.gen.Num;

/**Splits and filters a joint genotyped vcf file.*/
public class JointGenotypeVCFParser {

	private File vcfFile;

	private File saveDirectory;
	private double alleleFreq = 0.2;
	private double readDepth = 10;
	private double genotypeQuality = 20;
	private double qual = 20;
	private boolean debug = false;
	private HashMap<String, String> formatValues = new HashMap<String,String>();
	private String workingHeader = null;
	private String[] sampleNames = null;
	private int[] numPass = null;
	private int[] numFail = null;
	private int numFailingQual = 0;
	private int numRecords;
	
	public static final Pattern AF = Pattern.compile("AF=[\\d\\.]+;");
	public static final Pattern DP = Pattern.compile("DP=[\\d\\.]+;");	
	
	public JointGenotypeVCFParser (String[] args) {

		processArgs(args);
		
		System.out.println("Thresholds:");
		System.out.println(alleleFreq + "\tMinimum allele frequency based on AD");
		System.out.println((int)readDepth + "\tMinimum read depth based on AD");
		System.out.println((int)genotypeQuality + "\tMinimum genotype quality, GQ");
		System.out.println(qual + "\tMinimum QUAL value");

		parse();
		
	}
	
	public BufferedReader loadHeader() throws Exception{
		BufferedReader in = IO.fetchBufferedReader(vcfFile);
		String line = null;
		StringBuilder sb = new StringBuilder();
		//for each line in the file
		while ((line = in.readLine()) != null){
			line = line.trim();
			//header? just print out
			if (line.startsWith("#")) {
				if (line.startsWith("#CHROM")){
					sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
					workingHeader = sb.toString();
					
					//save sample names
					String[] columns = Misc.TAB.split(line);
					int numSamples = columns.length- 9;
					sampleNames = new String[numSamples];
					int index = 0;
					for (int i= 9; i< columns.length; i++) sampleNames[index++] = columns[i];
					return in;
				}
				sb.append(line);
				sb.append("\n");
			}
		}
		throw new Exception("\nERROR: failed to find the #CHROM line?");
	}

	public void parse() {
		String line = null;
		boolean failed = false;
		BufferedReader in = null;
		Gzipper[] out = null;
		try {
			//load the header and get the sample count
			in = loadHeader();
			int numSamples = sampleNames.length;

			out = new Gzipper[numSamples];
			//sampleName 1218982_NormalExome_Hg38 plus _JointGenotyped.vcf.gz
			for (int i=0; i< numSamples; i++){
				out[i] = new Gzipper( new File(saveDirectory, sampleNames[i]+"_JointGenotyped.vcf.gz"));
				out[i].println(workingHeader+sampleNames[i]);
			}
			numPass = new int[numSamples];
			numFail = new int[numSamples];

			//for each line in the file
			while ((line = in.readLine()) != null){
				if (debug) System.out.println("\n"+line);
				numRecords++;
				
				//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1, Sample2.....
				//   0    1   2  3   4   5     6      7     8      9       10
				String[] fields = Misc.TAB.split(line);

				//check whole line QUAL
				double q = Double.parseDouble(fields[5]);
				if (q < qual) {
					numFailingQual++;
					continue;
				}

				String[] format = Misc.COLON.split(fields[8]);
				//build partial line
				StringBuilder sb = new StringBuilder();
				for (int i=0; i< 7; i++){
					sb.append(fields[i]);
					sb.append("\t");
				}
				String partialVcf = sb.toString();

				//for each sample
				int index = 0;
				for (int i=9; i< fields.length; i++){
					double[] dpAf = checkSample(format, fields[i]);
					if (dpAf == null) numFail[index]++;
					else {
						numPass[index]++;
						String modInfo = replaceAfDb(fields[7], dpAf);
						StringBuilder sub = new StringBuilder(partialVcf);
						sub.append(modInfo);
						sub.append("\t");
						sub.append(fields[8]);
						sub.append("\t");
						sub.append(fields[i]);
						out[index].println(sub);
					}
					index++;
				}
			}

			//print stats
			printStats();
			
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("ERROR: parsing line "+line);
			failed = true;
		} finally {
			//close io
			try {
				for (Gzipper g: out) {
					g.closeNoException();
					if (failed) g.getGzipFile().delete();
				}
				if (in != null) in.close();
			} catch (IOException e) {}
			if (failed) System.exit(1);
		}
	}

	private void printStats() {
		System.out.println("\nParsing Statistics:");
		System.out.println(numRecords+ "\tNumber records");
		System.out.println(numFailingQual+ "\tNumber failing QUAL");
		System.out.println("\nPassing\tFailing\tSampleName");
		for (int i=0; i< sampleNames.length; i++){
			System.out.println(numPass[i] + "\t"+ numFail[i] + "\t"+ sampleNames[i] );
		}
		
	}

	private String replaceAfDb(String info, double[] dpAf) throws Exception {
		String afString = "AF="+ Num.formatNumberMinOne(dpAf[1], 4)+";";
		String dpString = "DP="+ (int)dpAf[0] +";";
		
		//check it contains AF= and DP=
		if (info.contains("AF=") == false || info.contains("DP=") == false) throw new Exception("Failed to find AF= and/or DP= in "+info);
		
		String modAF = AF.matcher(info).replaceFirst(afString);
		String finalMod = DP.matcher(modAF).replaceFirst(dpString);
		return finalMod;
	}

	/**Returns null if fails, otherwise the DP and AF.*/
	private double[] checkSample(String[] ids, String sample) throws Exception {
		if (debug) System.out.println("\t"+sample);
		//no values?
		if (sample.startsWith("./.")) {
			if (debug) System.out.println ("\t\tEmpty");
			return null;
		}
		
		formatValues.clear();
		String[] values = Misc.COLON.split(sample);
		if (ids.length != values.length) throw new Exception ("\t\tThe FORMAT length and sample values do not match");
		
		for (int i=0; i< ids.length; i++) formatValues.put(ids[i], values[i]);
		
		//check GQ?
		if (genotypeQuality != 0){
			String gqString = formatValues.get("GQ");
			if (gqString == null || gqString.equals(".")) {
				if (debug) System.out.println ("\t\tNo GQ");
				return null;
			}
			else {
				double gq = Double.parseDouble(gqString);
				if (gq < genotypeQuality) {
					if (debug) System.out.println ("\t\tFailing GQ");
					return null;
				}
			}
		}
		
		//check depth and allele freq, should be only one alt!
		String ad = formatValues.get("AD");
		if (ad == null) throw new Exception ("\t\tNo AD?");

		String[] refAltCounts = Misc.COMMA.split(ad);
		if (refAltCounts.length != 2) throw new Exception ("\nThe AD for this sample does not contain 2 values, did you vt decompose? "+ad+" "+sample);
		double refCounts = Double.parseDouble(refAltCounts[0]);
		double altCounts = Double.parseDouble(refAltCounts[1]);
		double dp = refCounts + altCounts;
		double af = altCounts/dp;
		if (dp >= readDepth && af >= alleleFreq) return new double[]{dp, af};
		
		if (debug) System.out.println ("\t\tFailing dp ("+dp+") or af ("+af+")");
		return null;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new JointGenotypeVCFParser(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcfFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'd': readDepth = Double.parseDouble(args[++i]); break;
					case 'a': alleleFreq = Double.parseDouble(args[++i]); break;
					case 'q': qual = Double.parseDouble(args[++i]); break;
					case 'g': genotypeQuality = Double.parseDouble(args[++i]); break;
					case 'f': debug = true; break;
					
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

		if (vcfFile == null || vcfFile.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file to split and filter.\n");	
		if (saveDirectory == null) Misc.printErrAndExit("\nCannot find your save directory?! "+saveDirectory);
		saveDirectory.mkdirs();
		
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Joint Genotype VCF Parser: Oct 2018                     **\n" +
				"**************************************************************************************\n" +
				"Splits and filters GATK joint genotyped multi sample vcf files. Use vt to decompose \n"+
				"the multi alts. See https://genome.sph.umich.edu/wiki/Vt#Decompose . Replaces the AF\n"+
				"and DP INFO fields with the sample level values.\n"+

				"\nRequired Params:\n"+
				"-v Path to vt decomposed GATK joint genotyped multi sample vcf file, gz/zip OK.\n"+
				"     ~/BioApps/vt decompose -s jointGenotyped.vcf.gz -o jointGenotyped.decomp.vcf.gz\n" +
				"-s Path to a directory to save the split files.\n"+
				
				"\nOptional Params:\n"+
				"-q Minimum QUAL value, defaults to 20\n"+
				"-d Minimum read depth based on the AD sample values, defaults to 10\n"+
				"-a Minimum AF allele freq, defaults to 0.2\n"+
				"-g Minimum GT genotype quality, defaults to 20\n"+
				"-f Print debugging output to screen\n"+

				"\nExample: java -jar -Xmx2G pathToUSeq/Apps/HaplotypeVCFParser -d 20 -a 0.25 -g 30 -f\n"+
				"      -v jointGenotyped.decomp.vcf.gz -s SplitFilteredVcfs/ -q 30 \n\n" +


				"**************************************************************************************\n");

	}

}
