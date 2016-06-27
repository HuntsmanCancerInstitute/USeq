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

/**App to score Pindel results with filtering. */
public class PindelVCFParser {

	private File[] vcfFiles;
	private int minimumReadDepth = 100;
	private double minimumAllelicRatio = 0.05;
	private String[] keys = {"13:28608243", "13:28608218", "13:28608262", "13:28608259", "13:28608257", "13:28608226", "13:28608239", "13:28608223", "13:28608232", "13:28608266", "13:28608273", "13:28608262"};
	private HashSet<String> key = null;
	
	public PindelVCFParser (String[] args) { 

		processArgs(args);
		
		System.out.println("Thresholds:");
		System.out.println(minimumAllelicRatio+"\tMin allelic freq AF");
		System.out.println(minimumReadDepth+"\tMin depth DP");
		
		key = fetchKey();
		
		System.out.println("\nName\tPassing\tFailing\tInKey");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}
	

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);

			//set all to pass
			parser.setFilterFieldOnAllRecords(VCFRecord.PASS);
			HashSet<String> vars = new HashSet<String>();
			
			for (VCFRecord r: parser.getVcfRecords()){	
				VCFSample[] tumor = r.getSample();
				
				
				//check if dup
				String test = r.getChromosome()+":"+(r.getPosition()+1);
//System.out.println(test);
				
boolean x = test.equals("xxx");
if (x) System.out.println("Found it!" +r.getOriginalRecord());
				
				if (vars.contains(test)){
if (x) System.out.println("Already seen!");
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				vars.add(test);
				
				//check depth
				if (minimumReadDepth !=0 && tumor[0].getReadDepthDP() < minimumReadDepth) {
					r.setFilter(VCFRecord.FAIL);
if (x) System.out.println("Failed min depth");
					continue;
				}
				
				//check alt allelic ratio 
				if (minimumAllelicRatio != 0.0 && tumor[0].getAltRatio() < minimumAllelicRatio){
					r.setFilter(VCFRecord.FAIL);
if (x) System.out.println("Failed min AF");
					continue;
				}
			}
			
			scoreRecords(parser);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public HashSet<String> fetchKey(){
		HashSet<String> key = new HashSet<String>();
		for (String x : keys) key.add(x);
		return key;
	}

	public void scoreRecords(VCFParser parser) throws Exception{
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
	}
	
	

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PindelVCFParser(args);
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
					case 'd': minimumReadDepth = Integer.parseInt(args[++i]); break;
					case 'a': minimumAllelicRatio = Double.parseDouble(args[++i]); break;
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
				"**                             Pindel VCF Parser: April 2016                        **\n" +
				"**************************************************************************************\n" +
				"Beta. \n\n" +


				"**************************************************************************************\n");

	}

}
