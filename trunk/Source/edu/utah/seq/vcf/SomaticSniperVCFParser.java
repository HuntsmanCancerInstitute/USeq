package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

/**Simple SomaticSniper vcf formatter.
 * takes the tumor ssc score and replaces the qual. Can also filter for minimum read depth count in Tum and Normal, as well as minimum Alt Allele fraction change. */
public class SomaticSniperVCFParser {

	private File[] vcfFiles;
	private float minimumSSC = 0;
	private int minimumCount = 0;
	private double minimumAbsFractionChange = 0;
	private double maximumNormalAltFraction = 1;

	public SomaticSniperVCFParser (String[] args) {

		processArgs(args);
		
		System.out.println("Thresholds:");
		System.out.println(minimumSSC+"\tMin SSC score");
		System.out.println(minimumCount+"\tMin alignment depth");
		System.out.println(minimumAbsFractionChange+"\tMin Allelic fraction change");
		System.out.println(maximumNormalAltFraction+"\tMax Normal alt allelic fraction");
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);
			if (parser.getSampleNames()[1].equals("TUMOR") == false) Misc.printErrAndExit("TUMOR doesn't appear to be the second sample in the VCF file?! "+vcf.getName());
			for (VCFRecord r: parser.getVcfRecords()){
				VCFSample[] normTum = r.getSample();
				//check depth
				int normDepth = normTum[0].getReadDepthDP();
				int tumDepth = normTum[1].getReadDepthDP();
				if (normDepth < minimumCount || tumDepth < minimumCount) {
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				//check allelic ratio shift, SomSni has DP4
				if (minimumAbsFractionChange != 0.0){
					double normRto = normTum[0].getAltRatio();
					double tumRto = normTum[1].getAltRatio();
					double change = Math.abs(normRto-tumRto);
					if (change < minimumAbsFractionChange){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check normal alt fraction?
				if (maximumNormalAltFraction !=1){
					double normRto = normTum[0].getAltRatio();
					if (normRto > maximumNormalAltFraction){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				
				//set SSC score
				VCFSample tumor = normTum[1];
				String[] format = tumor.getFormat();
				if (format[format.length-1].equals("SSC") == false) Misc.printErrAndExit("SSC isn't the last FORMAT in record "+r.getOriginalRecord());
				r.setQuality(Float.parseFloat(tumor.getData()[format.length-1]));
			}
			File out = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"_SSC.vcf");
			printRecords(parser, out);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void printRecords(VCFParser parser, File f) throws Exception{
		PrintWriter out = new PrintWriter (new FileWriter (f));
		int numFail = 0;
		int numPass = 0;
		//write out header
		for (String h : parser.getStringComments()) out.println(h);
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(VCFRecord.FAIL) || vcf.getQuality() < minimumSSC) numFail++;
			else {
				numPass++;
				String orig = vcf.toString();
				String[] fields = VCFParser.TAB.split(orig);
				//reset ID
				fields[2] = "SomSni"+numPass;
				//reset score
				fields[5] = Integer.toString((int)vcf.getQuality());
				out.println(Misc.stringArrayToString(fields, "\t"));
			}
		}
		out.close();
		System.out.println(numPass+"\t"+numFail);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SomaticSniperVCFParser(args);
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
					case 'm': minimumSSC = Float.parseFloat(args[++i]); break;
					case 'a': minimumCount = Integer.parseInt(args[++i]); break;
					case 'r': minimumAbsFractionChange = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
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

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Somatic Sniper VCF Parser: Jan 2015                  **\n" +
				"**************************************************************************************\n" +
				"Parses Somatic Sniper VCF files, replacing the QUAl score with the SSC score. Also \n"+
				"filters for minimum tumor normal read depth, difference in alt allelic ratios, and on\n"+
				"normal alt allelic rations.\n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-m Minimum SSC score, defaults to 0.\n"+
				"-a Minimum alignment depth for both tumor and normal samples, defaults to 0.\n"+
				"-r Minimum absolute difference in alt allelic ratios, defaults to 0.\n"+
				"-n Maximum normal alt allelic fraction, defaults to 1.\n"+

				"\nExample: java -jar pathToUSeq/Apps/SomaticSniperVCFParser -v /VCFFiles/ -m 32 -a 15\n"+
				"      -r 0.25 -n 0.02\n\n" +


				"**************************************************************************************\n");

	}

}
