package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Scalpel filtering app. */
public class ScalpelVCFParser {

	private File[] vcfFiles;
	private float minimumQual = 0;
	private int minimumCount = 100;
	private double minimumAltTNRatio = 1.6;
	private double maximumNormalAltFraction = 0.1;
	private double minimumTumorAltFraction = 0.05;
	
	//private String[] keys = {"13:28608243", "13:28608218", "13:28608262", "13:28608259", "13:28608257", "13:28608226", "13:28608239", "13:28608223", "13:28608232", "13:28608266", "13:28608273", "13:28608262"};
	//private HashSet<String> key = null;
	
	public ScalpelVCFParser (String[] args) { 

		processArgs(args);
		
		System.out.println("Thresholds:");
		System.out.println(minimumQual+"\tMin QUAL score");
		System.out.println(minimumCount+"\tMin alignment depth");
		System.out.println(minimumAltTNRatio+"\tMin Allelic fraction change tAltRto/nAltRto");
		System.out.println(maximumNormalAltFraction+"\tMax Normal alt allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin Tumor alt allelic fraction");
		
		//key = fetchKey();
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}
	
	/*public HashSet<String> fetchKey(){
		HashSet<String> key = new HashSet<String>();
		for (String x : keys) key.add(x);
		return key;
	}*/
	

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);
			
			//set all to pass
			parser.setFilterFieldOnAllRecords(VCFRecord.PASS);
			
			for (VCFRecord r: parser.getVcfRecords()){	
				VCFSample[] normTum = r.getSample();
				
				//check depth
				int normDepth = normTum[0].getReadDepthDP();
				int tumDepth = normTum[1].getReadDepthDP();
				if (normDepth < minimumCount || tumDepth < minimumCount) {
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				
				//check QUAL
				if (r.getQuality() < minimumQual){
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				
				double normRto = normTum[0].getAltRatio();
				double tumRto = normTum[1].getAltRatio();
				
				//check alt allelic ratio 
				if (minimumAltTNRatio != 0.0){
					if (normRto == 0) normRto = 0.001;
					if (tumRto == 0) tumRto = 0.001;
					double rto = tumRto/normRto;
					if (rto < minimumAltTNRatio){
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
			File out = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"_Filtered.vcf.gz");
			printRecords(parser, out);
			
			//scoreRecords(parser);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/*public void scoreRecords(VCFParser parser) throws Exception{
		int numFail = 0;
		int numPass = 0;
		int numOnTarget = 0;
		
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(VCFRecord.FAIL)) numFail++;
			else {
				numPass++;
				String test = vcf.getChromosome()+":"+(vcf.getPosition()+1);
				if (key.contains(test)) numOnTarget++;
			}
		}
		System.out.println(numPass+"\t"+numFail+"\t"+numOnTarget);
	}*/
	
	public void printRecords(VCFParser parser, File f) throws Exception{
		Gzipper out = new Gzipper (f);
		
		//write out header inserting AF
		writeHeaderWithExtraInfo(out, parser);
		
		int numFail = 0;
		int numPass = 0;
		
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(VCFRecord.FAIL)) numFail++;
			else {
				numPass++;
				String orig = vcf.toString();
				out.println(orig);
			}
		}
		out.close();
		System.out.println(numPass+"\t"+numFail);
	}
	
	private void writeHeaderWithExtraInfo(Gzipper out, VCFParser parser) throws IOException {
		for (String h : parser.getStringComments()) {
			out.println(h);
		}
	}


	

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ScalpelVCFParser(args);
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
					case 'm': minimumQual = Float.parseFloat(args[++i]); break;
					case 'a': minimumCount = Integer.parseInt(args[++i]); break;
					case 'r': minimumAltTNRatio = Double.parseDouble(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
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
				"**                            Scalpel VCF Parser: May 2016                          **\n" +
				"**************************************************************************************\n" +
				"Filters Scalpel VCF INDEL files.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				"-m Minimum QUAL score, defaults to 0\n"+
				"-a Minimum alignment depth for both tumor and normal samples, defaults to 100\n"+
				"-t Minimum tumor alt allelic fraction, defaults to 0.05\n"+
				"-n Maximum normal alt allelic fraction, defaults to 0.1\n"+
				"-r Minimum alt Tum/Norm allelic ratio, defaults to 1.6.\n"+

				"\nExample: java -jar pathToUSeq/Apps/Scalpel/VCFParser -v /VCFFiles/ -m 150 -a 200\n"+
				"      -r 2 -n 0.2 -t 0.03  \n\n" +


				"**************************************************************************************\n");

	}
}
