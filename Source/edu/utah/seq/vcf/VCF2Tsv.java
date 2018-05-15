package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Converts a vcf file into a simple chr pos ref alt %AF for Illumina CE/CA uploads.*/
public class VCF2Tsv {

	private File[] vcfFiles;
	private File saveDirectory;
	private double minimumQual = 0;
	private boolean ignoreNoBKZThresholding = false;
	private double maxBKAF = 1;
	private Pattern keepId = null;
	private String keepIdString = "";

	public VCF2Tsv (String[] args) {

		processArgs(args);

		printSettings();

		System.out.println("\nName\tPass\tFail");
		for (File vcf: vcfFiles) {
			System.out.print(vcf.getName());
			parse(vcf);
		}
		System.out.println("\nDone!");

	}

	public void parse(File vcfFile){

		File bed = new File (saveDirectory, Misc.removeExtension(vcfFile.getName())+".tsv.gz");
		String line = null;
		int numPrinted = 0;
		int numNotPrinted = 0;
		VCFParser vp = new VCFParser();

		try {
			Gzipper out = new Gzipper(bed);
			BufferedReader in = IO.fetchBufferedReader(vcfFile);

			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				VCFRecord vcf = new VCFRecord(line, vp, false, true);
				boolean parse = false;
				double af = 0;

				//skip structural vars even in Id match calls
				if (vcf.getAlternate()[0].startsWith("<") == false) {
					
					
					//check quality
					float qual = vcf.getQuality();
					if (qual ==0 && line.contains("BKZ=") == false && ignoreNoBKZThresholding) parse = true;
					if (parse == false){
						if (qual >= minimumQual) parse = true;
					}
					
					//check for multi alts
					if (vcf.getAlternate().length != 1){
						parse = false;
					}

					//check frac BKAF?
					af = vcf.getInfoObject().getInfoFloat("AF");
					if (parse == true && maxBKAF !=1){
						String bkafs = vcf.getInfoObject().getInfo("BKAF");
						if (bkafs != null){
							double[] bks = Num.parseDoubles(bkafs, Misc.COMMA);
							double exceeds = 0;
							for (double d: bks) if (d>=af) exceeds++;
							exceeds = exceeds/(double)bks.length;
							//IO.p("\nBKAFs check "+exceeds+" "+Num.doubleArrayToString(bks, ":"));
							if (exceeds >= maxBKAF) parse = false;
						}
					}


					//check id field, this overrides prior score thresholds
					if (keepId != null && parse == false){
						if (keepId.matcher(vcf.getRsNumber()).matches()) parse = true;
					}
					
					//check for . in ref or alt, these are junk vars from the old Strelka
					if (vcf.getReference().equals(".") || vcf.getAlternate()[0].equals(".")) parse = false;

				}

				if (parse){
					//check for bad AFs
					if (af > 1.0) af = 1.0;
					String afPercent = Num.formatNumber(af*100, 1);
					
					if (vcf.getAlternate().length !=1) throw new Exception("More than one alt, use vt decompose "+line);
					out.println(vcf.getChromosome()+"\t"+(vcf.getPosition()+1)+ "\t"+ vcf.getReference()+"\t"+vcf.getAlternate()[0]+"\t"+afPercent);
					numPrinted++;
					//System.out.println("Pass\t"+line);
				}
				else {
					numNotPrinted++;
					//System.out.println("Fail\t"+line);
				}
			}
			out.close();
			System.out.println("\t"+numPrinted+"\t"+numNotPrinted);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing this record "+line);
		} 
	}

	public void printSettings(){
		IO.p("Thresholds:");
		IO.p("Save Dir:\t"+saveDirectory);
		IO.p("Min Qual:\t"+minimumQual);
		IO.p("Ignore NoBKZ:\t"+ignoreNoBKZThresholding);
		IO.p("Max BKAF Frac:\t"+maxBKAF);
		IO.p("Keep ID String:\t"+keepIdString);
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCF2Tsv(args);
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
					case 'q': minimumQual = Double.parseDouble(args[++i]); break;
					case 'b': maxBKAF = Double.parseDouble(args[++i]); break;
					case 'i': keepIdString = args[++i]; break;
					case 'z': ignoreNoBKZThresholding = true; break;
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

		if (keepIdString.equals("") == false) keepId = Pattern.compile(".*"+keepIdString+".*", Pattern.CASE_INSENSITIVE);
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 VCF 2 Tsv: March 2018                            **\n" +
				"**************************************************************************************\n" +
				"Converts vcf files' SNVs and INDELs to tsv Illumina CE/CA format. \n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-s Directory to save the tvs files, defaults to the parent of the vcf\n"+
				"-q Minimum QUAL, defaults to 0, no threshold.\n"+
				"-z Ignore QUAL thresholding when BKZ= is absent from INFO field.\n"+
				"-b Max frac BKAFs >= AF, defaults to 1, no threshold.\n"+
				"-i ID column string forcing export. Case insensitive.\n"+

				"\nExample: java -jar pathToUSeq/Apps/VCF2Tsv -v VCFss/ -q 3 -z -b 0.2 -i Foundation \n\n" +
				"**************************************************************************************\n");

	}

}
