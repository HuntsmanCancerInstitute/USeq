package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

/**Simple MuTech vcf parser. Filters for minimum read depth count in Tum and Normal, as well as minimum Alt Allele fraction change and normal Alt fraction */
public class MuTectVCFParser {

	private File[] vcfFiles;
	private int minimumCount = 0;
	private double minimumAbsFractionChange = 0;
	private double maximumNormalAltFraction = 1;
	private double minimumTumorAltFraction = 0;
	
	public MuTectVCFParser (String[] args) {

		processArgs(args);
		
		System.out.println("Thresholds:");
		System.out.println(minimumCount+"\tMin alignment depth");
		System.out.println(minimumAbsFractionChange+"\tMin T-N allelic fraction change");
		System.out.println(maximumNormalAltFraction+"\tMax Normal alt allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin Tumor alt allelic fraction");
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);
			for (VCFRecord r: parser.getVcfRecords()){			
				VCFSample[] tumNorm = r.getSample();
				//check depth
				int normDepth = tumNorm[1].getReadDepthDP();
				int tumDepth = tumNorm[0].getReadDepthDP();
				if (normDepth < minimumCount || tumDepth < minimumCount) {
					r.setFilter(VCFRecord.FAIL);					
					continue;
				}
				double normRto = tumNorm[1].getAltRatio();
				double tumRto = tumNorm[0].getAltRatio();
				
				//check allelic ratio shift, SomSni has DP4
				if (minimumAbsFractionChange != 0){
					double change = Math.abs(normRto-tumRto);
					if (change < minimumAbsFractionChange){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check normal alt fraction?
				if (maximumNormalAltFraction !=1){
					if (normRto > maximumNormalAltFraction){	
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check tumor alt fraction?
				if (minimumTumorAltFraction !=0){
					if (tumRto < minimumTumorAltFraction){	
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
			}
			printRecords(parser, vcf);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void printRecords(VCFParser parser, File f) throws Exception{
		
		String name = Misc.removeExtension(f.getName());
		File passFile = new File(f.getParentFile(), name+ "_pass.vcf");
		PrintWriter outPass = new PrintWriter (new FileWriter (passFile));
		File failFile = new File(f.getParentFile(), name+ "_fail.vcf");
		PrintWriter outFail = new PrintWriter (new FileWriter (failFile));
		
		int numFail = 0;
		int numPass = 0;

		for (String h : parser.getStringComments()) {
			outPass.println(h);
			outFail.println(h);
		}
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(VCFRecord.FAIL)) {
				numFail++;
				outFail.println(vcf.getOriginalRecord());
			}
			else {
				numPass++;
				outPass.println(vcf.getOriginalRecord());
			}
		}
		outPass.close();
		outFail.close();
		if (numPass == 0) passFile.deleteOnExit();
		if (numFail == 0) failFile.deleteOnExit();
		
		System.out.println(numPass+"\t"+numFail);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MuTectVCFParser(args);
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
					case 'a': minimumCount = Integer.parseInt(args[++i]); break;
					case 'r': minimumAbsFractionChange = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'm': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
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
				"**                               MuTect VCF Parser: Aug 2016                        **\n" +
				"**************************************************************************************\n" +
				"Parses MuTech2 VCF files, filtering for minimum read depth, difference in allelic\n"+
				"ratios, and allelic fractions. Splits files into pass and fail with no modifications.\n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-a Minimum alignment depth for tumor and normal samples, defaults to 0.\n"+
				"-r Minimum absolute difference in alt allelic ratios, defaults to 0.\n"+
				"-n Maximum normal alt allelic fraction, defaults to 1.\n"+
				"-m Minimum tumor alt allelic fraction, defaults to 0.\n"+

				"\nExample: java -jar pathToUSeq/Apps/MuTechVCFParser -v /VCFFiles/ -a 15 -r 0.10\n\n"+


				"**************************************************************************************\n");

	}

}
