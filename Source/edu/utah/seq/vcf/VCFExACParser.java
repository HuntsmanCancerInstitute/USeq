package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Parses the AF field from a VCF files Info field, outputs just this and the other vcf fields.*/
public class VCFExACParser {

	private File vcfFile;
	private static final Pattern AF = Pattern.compile(".+AF=([\\d\\.e-]+);.+");

	public VCFExACParser (String[] args) {

		processArgs(args);

		String line = null;
		long numVars = 0;
		try {
			File parsedVcf = new File ( Misc.removeExtension(vcfFile.getCanonicalPath())+".parsed.vcf.gz");
			Gzipper out = new Gzipper(parsedVcf);
			BufferedReader in = IO.fetchBufferedReader(vcfFile);
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0) continue;
				else if (line.startsWith("#")) {
					if (line.startsWith("##INFO=")) {
						if (line.startsWith("##INFO=<ID=AF,")) out.println(line);
					}
					else out.println(line);
				}
				else {
					numVars++;
					//#CHROM0  POS1  ID2  REF3  ALT4     QUAL5    FILTER6  INFO7
					String[] fields = Misc.TAB.split(line);
					Matcher mat = AF.matcher(fields[7]);
					if (mat.matches()) fields[7] = "AF="+mat.group(1);
					else Misc.printErrAndExit("Did you Vt normalize the VCF? Failed to find the AF= in "+line);
					out.println(Misc.stringArrayToString(fields, "\t"));
				}
			}
			out.close();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing this record "+line);
		} 
		System.out.println("Done! "+numVars+" records parsed.");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFExACParser(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcfFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF ExAC Parser: April 2021                         **\n" +
				"**************************************************************************************\n" +
				"Extracts the AF= value from the INFO field and uses it for the INFO.\n"+

				"\nRequired Options:\n"+
				"-v Path to the xxx.vcf(.gz/.zip OK) file.\n" +

				"\nExample: java -jar pathToUSeq/Apps/VCFExACParser -v ExAC.r1.sites.vep.hg38.vcf\n\n" +
				"**************************************************************************************\n");

	}

}
