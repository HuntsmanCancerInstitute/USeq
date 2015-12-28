package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;

/**Lofreq formatter and parser. */
public class GenericMungeVCFParser {

	private File[] vcfFiles;
	private float alleleFreq = 0.04f;
	private float readDepth = 100;
	private Histogram af = new Histogram(alleleFreq, 1.05, 25);
	private Gzipper out;

	public GenericMungeVCFParser (String[] args) throws FileNotFoundException, IOException {

		processArgs(args);
		out = new Gzipper (new File (vcfFiles[0].getParentFile(), "passing.vcf"));
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
		
		System.out.println();
		af.printScaledHistogram();
		out.close();
	}

	public void parse(File vcf) {
		String line = null;
		try {
			//counters
			int numFail = 0;
			int numPass = 0;

			//IO
			BufferedReader in = IO.fetchBufferedReader(vcf);

			//for each line in the file
			
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#")) continue;
				//data line

				//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT	NORMAL	TUMOR
				//   0    1   2  3   4   5     6      7     8     9       10
				String[] tokens = Misc.TAB.split(line);


				float[] dpAf;

				//indel or snp
				//if (tokens[3].length() != tokens[4].length()) dpAf = parseDpAfStrelka(tokens[10]);
				if (tokens[2].startsWith("Strelka")) dpAf = parseDpAfStrelka(tokens[10]);
				else dpAf = parseDpAf(tokens[7]);

				//pass readDepth and alleleFreq
				if (dpAf[0] < readDepth || dpAf[1] < alleleFreq){
					numFail++;
				}
				else {
					numPass++;
					out.println(line);
					af.count(dpAf[1]);
				}

			}
			in.close();
			System.out.println(numPass+"\t"+numFail);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: parsing lofreq vcf "+vcf+" \n"+line);
		} 
	}



	private float[] parseDpAfStrelka(String string) {
		// DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:AF
		String[] t = Misc.COLON.split(string);
		float dp = Float.parseFloat(t[0]);
		float af = Float.parseFloat(t[t.length-1]);
		return new float[]{dp,af};
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
		try {
			new GenericMungeVCFParser(args);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
				"**                             VCF Parser: Dec 2015                          **\n" +
				"**************************************************************************************\n" +


				"**************************************************************************************\n");

	}

}
