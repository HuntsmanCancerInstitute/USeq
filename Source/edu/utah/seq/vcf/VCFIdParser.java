package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Quickie ID app parser and counter. */
public class VCFIdParser {

	private File[] vcfFiles;
	
	public VCFIdParser (String[] args) {

		processArgs(args);

		for (File vcf: vcfFiles){
			parse(vcf);
			System.out.println();
		}
	}

	public void parse(File vcf) {
		try {
			System.out.println(vcf.getName());
			BufferedReader in = IO.fetchBufferedReader(vcf);
			HashMap<String, Integer> typeCount = new HashMap<String, Integer>();
			Pattern toDrop = Pattern.compile("[\\d_;]+");
			
			//for each line in the file
			String line;
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#")) continue;
				//data line
				else {
					//#CHROM POS ID REF ALT QUAL FILTER INFO ......
					//   0    1   2  3   4   5     6      7
					String[] tokens = Misc.TAB.split(line);
					String key = toDrop.matcher(tokens[2]).replaceAll("");
					Integer count = typeCount.get(key);
					if (count == null) count = new Integer(1);
					else count = new Integer(1+count.intValue());
					typeCount.put(key, count);
				}
			}
			in.close();
			//System.out.println(typeCount+"");
			Misc.printHashMap(typeCount, "\t");
			
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: parsing lofreq vcf "+vcf);
		} 
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFIdParser(args);
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
				"**                             Lofreq VCF Parser: March 2017                        **\n" +
				"**************************************************************************************\n" +
				"NOT for distribution.  Parses vcf files ID column and counts what app called each var.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				


				"**************************************************************************************\n");

	}

}
