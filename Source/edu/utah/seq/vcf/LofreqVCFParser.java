package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Lofreq formatter and parser. */
public class LofreqVCFParser {

	private File[] vcfFiles;
	private float minimumScore = 0;
	private boolean appendFNT = false;
	private boolean removeIndels = false;
	private boolean clearFilter = false;
	private File saveDirectory;
	private float alleleFreq = 0;
	private float readDepth = 0;
	private boolean markFail = false;
	
	public LofreqVCFParser (String[] args) {

		processArgs(args);
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
			System.out.println();
		}
	}

	public void parse(File vcf) {
		try {
			//counters
			int numFail = 0;
			int numPass = 0;
			int counter = 0;

			//IO
			String name = Misc.removeExtension(vcf.getName());
			Gzipper modVcf = new Gzipper( new File (saveDirectory, name + "_Filtered.vcf.gz"));
			BufferedReader in = IO.fetchBufferedReader(vcf);
			
			//for each line in the file
			String line;
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#")) {
					if (appendFNT && line.startsWith("#CHROM")) {
						//add GT line, then concluding comment line
						modVcf.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
						modVcf.println(line+"\tFORMAT\tNORMAL\tTUMOR");
					}
					else modVcf.println(line);
				}
				//data line
				else {
					//#CHROM POS ID REF ALT QUAL FILTER INFO ......
					//   0    1   2  3   4   5     6      7
					String[] tokens = Misc.TAB.split(line);
					boolean pass = true;
					//skip indels?
					if (removeIndels){
						if (tokens[3].length() !=1 || tokens[4].length() !=1) {
							numFail++;
							pass= false;
						}
					}
					//pass score?
					if (pass && minimumScore !=0){
						int score = Integer.parseInt(tokens[5]);
						if (score < minimumScore) {
							numFail++;
							pass= false;
						}
					}
					//read depth or allele freq?
					if (pass && readDepth !=0 || alleleFreq !=0){
						float[] dpAf = parseDpAf(tokens[7]);
						if (dpAf[0] < readDepth || dpAf[1] < alleleFreq){
							numFail++;
							pass= false;
						}
					}
					if (pass) numPass++;
					//so it failed, do they want to still print it?
					else if (markFail == false) continue;
					
					
					//reset FILTER?
					if (clearFilter) tokens[6] = ".";
					if (pass == false && markFail) tokens[6] = "FAIL";
					
					//modify id
					tokens[2] = "Lofreq_"+counter;
					counter++;
					line = Misc.stringArrayToString(tokens, "\t");
					
					//print
					if (appendFNT) modVcf.println(line+"\tGT\t./.\t./.");
					else modVcf.println(line);
				}
			}
			in.close();
			modVcf.close();
			System.out.println(numPass+"\t"+numFail);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: parsing lofreq vcf "+vcf);
		} 
	}



	private float[] parseDpAf(String info) {
		//DP=400;AF=0.222500;SB=1;DP4=137,173,37,52
		String[] t = Misc.SEMI_COLON.split(info);
		float dp = -1f;
		float af = -1f;
		for (int i=0; i< t.length; i++){
			if (t[i].startsWith("DP=")) dp = Float.parseFloat(t[i].substring(3));
			else if (t[i].startsWith("AF=")) af = Float.parseFloat(t[i].substring(3));
		}
		if (dp == -1 || af == -1) Misc.printErrAndExit("\nError: failed to parse a DP or AF from "+info);
		return new float[]{dp, af};
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new LofreqVCFParser(args);
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
					case 'i': removeIndels = true; break;
					case 'f': clearFilter = true; break;
					case 'a': appendFNT = true; break;
					case 'n': markFail = true; break;
					case 'm': minimumScore = Float.parseFloat(args[++i]); break;
					case 'd': readDepth = Float.parseFloat(args[++i]); break;
					case 't': alleleFreq = Float.parseFloat(args[++i]); break;
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
				"**                             Lofreq VCF Parser: July 2016                         **\n" +
				"**************************************************************************************\n" +
				"Parses Lofreq vcf files with options for filtering for minimum QUAL, modifying the\n"+
				"FILTER field, removing non SNVs, and appending FORMAT info for downstream merging.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				
				"\nOptional Params:\n"+
				"-s File path to a directory for saving the modified vcfs\n"+
				"-m Minimum QUAL score, defaults to 0\n"+
				"-d Minimum DP read depth, defaults to 0\n"+
				"-t Minimum AF allele freq, defaults to 0\n"+
				"-i Remove non SNV records\n"+
				"-f Replace the FILTER field with '.'\n"+
				"-a Append FORMAT NORMAL TUMOR to #CHROM line and add empty columns to records\n"+
				"-n Mark variants failing thresholds FAIL instead of not printing\n"+

				"\nExample: java -jar pathToUSeq/Apps/LofreqVCFParser -v VCFFiles/ -m 32 -i -f -a\n"+
				"      -s FilteredLofreqVcfs/ \n\n" +


				"**************************************************************************************\n");

	}

}
