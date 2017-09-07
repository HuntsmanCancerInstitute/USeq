package edu.utah.seq.vcf;

import java.io.File;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Converts a vcf file into a simple bed w/ w/o padding.*/
public class VCF2Bed {

	private File[] vcfFiles;
	private File saveDirectory;
	private int padding = 0;
	private boolean onlyEnds = false;
	
	public VCF2Bed (String[] args) {

		processArgs(args);
		for (File vcf: vcfFiles) {
			System.out.print(vcf.getName());
			parse(vcf);
		}
		
		
	}

	public void parse(File vcf){
	
			File bed = new File (saveDirectory, Misc.removeExtension(vcf.getName())+"Pad"+ padding +"bp.bed.gz");
			try {
				Gzipper out = new Gzipper(bed);
				HashMap<String,RegionScoreText[]> chrReg = null;
				if (onlyEnds) chrReg = Bed.parseVcfFileForENDVars(vcf);
				else chrReg = Bed.parseVcfFile(vcf, padding, false);
				int numPrinted = 0;
				for (String chr: chrReg.keySet()){
					RegionScoreText[] regions = chrReg.get(chr);
					numPrinted += regions.length;
					for (RegionScoreText r: regions) out.println(r.getBedLine(chr));
				}
				out.close();
				if (numPrinted == 0) {
					bed.deleteOnExit();
					System.out.println("\nNo records to print?");
				}
				else System.out.println("\t"+numPrinted);
			} catch (Exception e) {
				e.printStackTrace();
				Misc.printErrAndExit("\nProblem parsing and saving bed file from "+vcf);
			} 
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCF2Bed(args);
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
					case 'p': padding = Integer.parseInt(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'e': onlyEnds = true; break;
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

		//save dir?
		if (saveDirectory == null) saveDirectory = vcfFiles[0].getParentFile();
		else saveDirectory.mkdirs();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  VCF 2 Bed: June 2017                            **\n" +
				"**************************************************************************************\n" +
				"Converts a vcf file to bed format.\n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-p Padding to expand each variants size, defaults to 0\n"+
				"-s Directory to save the bed files, defaults to the parent of the vcf\n"+
				"-e Print out only the END=xxx containing vcf records as the bed based on the end value.\n"+

				"\nExample: java -jar pathToUSeq/Apps/VCF2Bed -v /VCFFiles/ -p 25 \n\n" +
				"**************************************************************************************\n");

	}

}
