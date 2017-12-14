package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Lofreq formatter and parser. */
public class HaplotypeVCFParser {

	private File[] vcfFiles;

	private File saveDirectory;
	private double alleleFreq = 0;
	private double readDepth = 0;
	private double genotypeQuality = 0;
	private boolean debug = false;
	HashMap<String, String> formatValues = new HashMap<String,String>();
	
	public HaplotypeVCFParser (String[] args) {

		processArgs(args);
		
		System.out.println("Thresholds:");
		System.out.println(alleleFreq + "\tMinimum allele frequency based on AD");
		System.out.println((int)readDepth + "\tMinimum read depth based on AD");
		System.out.println((int)genotypeQuality + "\tMinimum genotype quality, GQ");
		
		System.out.println("\nParsing:");
		for (File vcf: vcfFiles){
			System.out.println(vcf.getName());
			parse(vcf);
			System.out.println();
		}
	}

	public void parse(File vcf) {
		Gzipper outPass = null;
		Gzipper outFail = null;
		String line = null;
		try {
			//counters
			int numFail = 0;
			int numPass = 0;

			//IO
			String name = Misc.removeExtension(vcf.getName());
			outPass = new Gzipper( new File (saveDirectory, name + "_Pass.vcf.gz"));
			outFail = new Gzipper( new File (saveDirectory, name + "_Fail.vcf.gz"));
			BufferedReader in = IO.fetchBufferedReader(vcf);
			
			//for each line in the file
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#")) {
					outPass.println(line);
					outFail.println(line);
				}
				
				//data line
				else {
					//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1, Sample2.....
					//   0    1   2  3   4   5     6      7     8      9       10
					String[] tokens = Misc.TAB.split(line);
					String[] ids = Misc.COLON.split(tokens[8]);
					if (debug) System.out.println("\n"+line);
					
					//for each sample, only one must pass
					boolean pass = false;
					for (int i=9; i< tokens.length; i++){
						pass = checkSample(ids, tokens[i]);
						if (pass) break;
					}
					if (debug) System.out.println("\t"+ pass);
					if (pass == false) {
						numFail++;
						outFail.println(line);
					}
					else {
						numPass++;
						outPass.println(line);
					}
				}
			}
			in.close();
			outPass.close();
			outFail.close();
			System.out.println("\tNumPass: "+numPass+"\tNumFail: "+numFail);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: parsing haplotype caller vcf file "+vcf+"\nline "+line);
			if (outPass != null){
				outPass.closeNoException();
				outPass.getGzipFile().delete();
			}
			if (outFail != null){
				outFail.closeNoException();
				outFail.getGzipFile().delete();
			}
		} 
	}

	private boolean checkSample(String[] ids, String sample) throws Exception {
		if (debug) System.out.println("\t"+sample);
		formatValues.clear();
		String[] values = Misc.COLON.split(sample);
		if (ids.length != values.length) {
			if (debug) System.out.println ("\t\tThe FORMAT length and sample values do not match");
			return false;
		}
		for (int i=0; i< ids.length; i++) formatValues.put(ids[i], values[i]);
		
		//check GQ?
		if (genotypeQuality != 0){
			String gqString = formatValues.get("GQ");
			if (gqString == null || gqString.equals(".")) {
				if (debug) System.out.println ("\t\tNo GQ");
				return false;
			}
			else {
				double gq = Double.parseDouble(gqString);
				if (gq < genotypeQuality) {
					if (debug) System.out.println ("\t\tFailing GQ");
					return false;
				}
			}
		}
		
		//check depth and allele freq for each alternate
		String ad = formatValues.get("AD");
		if (ad == null) {
			if (debug) System.out.println ("\t\tNo AD, skipping");
			return false;
		}
		String[] refAltCounts = Misc.COMMA.split(ad);
		if (refAltCounts.length < 2) throw new Exception ("\nThe AD for this sample does not contain 2 or more values");
		double refCounts = Double.parseDouble(refAltCounts[0]);
		
		//for each alt count
		for (int i=1; i < refAltCounts.length; i++){
			double altCounts = Double.parseDouble(refAltCounts[i]);
			double dp = refCounts + altCounts;
			double af = altCounts/dp;
			if (dp >= readDepth && af >= alleleFreq) return true;
		}
		
		if (debug) System.out.println ("\t\tFailing dp or af");
		return false;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new HaplotypeVCFParser(args);
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
					case 's': saveDirectory = new File(args[++i]); break;
					case 'd': readDepth = Double.parseDouble(args[++i]); break;
					case 'a': alleleFreq = Double.parseDouble(args[++i]); break;
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

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");
		
		if (saveDirectory != null){
			saveDirectory.mkdirs();
			if (saveDirectory.isDirectory() == false || saveDirectory.exists() == false) Misc.printErrAndExit("\nCannot find or make your save directory?! "+saveDirectory);
		}
		else saveDirectory = vcfFiles[0].getParentFile();
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Haplotype VCF Parser: Nov 2017                        **\n" +
				"**************************************************************************************\n" +
				"Parses GATK Haplotype called single and multi sample vcf files. Run GATKs standard gVCF\n"+
				"to joint genotyping and split the multi sample vcf by sample with their SelectVariants.\n"+
				"The FORMAT field must contain AD and GQ . For multi alts, only one must pass. For multi\n"+
				"sample vcf, only one sample must pass.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				
				"\nOptional Params:\n"+
				"-s File path to a directory for saving the modified vcfs\n"+
				"-d Minimum read depth based on the AD values, defaults to 0\n"+
				"-a Minimum AF allele freq, defaults to 0\n"+
				"-g Minimum GT genotype quality, defaults to 0\n"+
				"-f Print debugging output to screen\n"+

				"\nExample: java -jar -Xmx2G pathToUSeq/Apps/HaplotypeVCFParser -v SplitVCF -d 20 -a 0.2\n"+
				"      -g 20 -s FilteredVcfs/ \n\n" +


				"**************************************************************************************\n");

	}

}
