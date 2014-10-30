package edu.utah.seq.vcf;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

/**Simple VarScan SOMATIC vcf formatter.*/
public class VarScanVCFParser {
	
	private File[] vcfFiles;
	
	public VarScanVCFParser (String[] args) {
		
		processArgs(args);
		
		System.out.println("VCFFile\t#Somatic\t#NonSomatic");
		for (File vcf: vcfFiles){
			parse(vcf);
		}
	}
	
	public void parse(File vcf){
		VCFParser parser = new VCFParser (vcf);
		int numSomatic = 0;
		int numNon = 0;
		for (VCFRecord record: parser.getVcfRecords()){
			int ssc = -1;
			try {
				if (record.getInfoObject().doesInfoEntryExist("SOMATIC")){
					ssc = record.getInfoObject().getInfoInt("SSC");
				}
			} catch (Exception e) {}
			if (ssc != -1){
				record.setFilter(VCFRecord.PASS);
				record.setQuality((float)ssc);
				numSomatic++;
			}
			else {
				record.setFilter(VCFRecord.FAIL);
				numNon++;
			}
		}
		File out = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"_Som.vcf");
		parser.printFilteredRecords(out, VCFRecord.PASS);
		out = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"_NonSom.vcf");
		parser.printFilteredRecords(out, VCFRecord.FAIL);
		System.out.println(vcf.getName()+"\t"+numSomatic+"\t"+numNon);
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VarScanVCFParser(args);
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
				"**                             VarScan VCFParser: Oct 2014                          **\n" +
				"**************************************************************************************\n" +
				"Parses and filters VarScan VCF files for those called SOMATIC.  Replaces the QUAl\n" +
				"score with the ssc score.\n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s). Recursive!\n" +
				
				"\nExample: java -jar pathToUSeq/Apps/VarScanVCFParser -v /VarScan2/VCFFiles/\n\n" +
				

		"**************************************************************************************\n");

	}

}
