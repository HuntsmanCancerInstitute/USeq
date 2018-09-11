package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Simple Strelka INDEL vcf formatter and parser.
 * Takes the tumor QSI or QSS score and replaces the qual. 
 */
public class VCFAlleleFreqBinner {

	private File[] vcfFiles;
	private File saveDirectory;
	private double minAF = 0;
	private double maxAF = 1;
	
	
	public VCFAlleleFreqBinner (String[] args) { 

		processArgs(args);
		System.out.println("Thresholds:\n\tMinAF\t"+minAF+"\n\tMaxAF\t"+maxAF);
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}
	

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);
			String name=Misc.removeExtension(vcf.getName());
			Gzipper out = new Gzipper(new File(saveDirectory, name+"_"+minAF+"_"+maxAF+".vcf.gz"));
			//print header
			for (String s: parser.getStringComments()) out.println(s);
			//set counters
			int numPass = 0;
			int numFail = 0;
			//for each record
			for (VCFRecord r: parser.getVcfRecords()){	
				String afString = r.getInfoObject().getInfo("AF");
				double af = Double.parseDouble(afString);
				if (af >= minAF && af < maxAF) {
					out.println(r.getOriginalRecord());
					numPass++;
				}
				else numFail++;
			}
			out.close();
			System.out.println(numPass+"\t"+numFail);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFAlleleFreqBinner(args);
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
					case 'm': minAF = Double.parseDouble(args[++i]); break;
					case 'x': maxAF = Double.parseDouble(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
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
		if (saveDirectory == null) saveDirectory = vcfFiles[0].getParentFile();
		else saveDirectory.mkdirs();
		
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Strelka VCF Parser: June 2018                         **\n" +
				"**************************************************************************************\n" +
				"Parses Strelka VCF INDEL and SNV files, replacing the QUAl score with the QSI or QSS\n"+
				"score. Also filters for read depth, T/N alt allelic ratio and diff, ref/alt with '.',\n"+
				"and tumor and normal alt allelic ratios. Lastly, it inserts the tumor DP and AF info.\n"+
				"For somatic exome datasets sequenced at >100X unique observation read depth, try the \n"+
				"tuned tier filtering to select for lists with 1-3%, 4-10%, and 10-20% FDR. Follow the\n"+
				"example.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-m Minimum QSI or QSS score, defaults to 0.\n"+
				"-e Apply a tuned 100X exome QSI/S stringency tier: 1 (10-20%FDR), 2 (4-10%FDR),\n"+
				"         3 (1-3%FDR), defaults to 0, no tiered filtiering (25-65%FDR).\n"+
				"-t Minimum tumor allele frequency (AF), defaults to 0.\n"+
				"-n Maximum normal AF, defaults to 1.\n"+
				"-u Minimum tumor alignment depth, defaults to 0.\n"+
				"-a Minimum tumor alt count, defaults to 0.\n"+
				"-o Minimum normal alignment depth, defaults to 0.\n"+
				"-d Minimum T-N AF difference, defaults to 0.\n"+
				"-r Minimum T/N AF ratio, defaults to 0.\n"+
				"-p Remove non PASS filter field records.\n"+
				"-s Print spreadsheet variant summary.\n"+
				"-f Directory in which to save the parsed files, defaults to the parent dir of the vcfs.\n"+

				"\nExample: java -jar pathToUSeq/Apps/StrelkaVCFParser -v /VCFFiles/ -t 0.03 -n 0.6 \n"+
				"-u 25 -o 10 -a 3 -d 0.03 -r 2 -e 2 \n\n"+


				"**************************************************************************************\n");

	}

}
