package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Recurses through finding vcf files and calculates a median AF on the T_DP= INFO column value.*/
public class VCFMedianAF {

	private File[] vcfFiles;
	//;T_DP=2411;
	private Pattern taf = Pattern.compile(".+;T_DP=(\\d+);.+");
	
	public VCFMedianAF (String[] args) {

		processArgs(args);
		for (File vcf: vcfFiles) {
			IO.p(vcf+"\t");
			parse(vcf);
		}
		
		
	}

	public void parse(File vcf){
		BufferedReader in = null;
			try {
				String line;
				in = IO.fetchBufferedReader(vcf);
				ArrayList<Double> afs = new ArrayList<Double>();
				while ((line = in.readLine())!= null) {
					if (line.startsWith("#") == false) {
						Matcher mat = taf.matcher(line);
						if (mat.matches()) {
							Double af = new Double(mat.group(1));
							afs.add(af);
						}
						else throw new Exception("\nFailed to parse AF from "+line);
					}
				}
				//calc median AF
				double[] dAFs = Num.arrayListOfDoubleToArray(afs);
				Arrays.sort(dAFs);
				if (dAFs.length > 2)IO.pl(Num.median(dAFs));
				else if (dAFs.length >0) IO.pl(Num.mean(dAFs));
				else IO.pl("NA");
				
			} catch (Exception e) {
				if (in != null) IO.closeNoException(in);
				e.printStackTrace();
				Misc.printErrAndExit("\nProblem parsing and saving bed file from "+vcf);
			} finally {
				if (in != null) IO.closeNoException(in);
			}
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFMedianAF(args);
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
					case 'v': vcfFiles = IO.fetchFilesRecursively(new File(args[++i]), ".vcf.gz"); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

		
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf.gz file(s)!\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 VCF Median AF: April 2020                        **\n" +
				"**************************************************************************************\n" +
				"Recurses through directories finding vcf files and calculates a median AF on the\n"+
				"T_DP= INFO column value.\n"+

				"\nRequired Options:\n"+
				"-v Directory containing xxx.vcf.gz file(s), recursive.\n" +
				

				"\nExample: java -jar pathToUSeq/Apps/VCFMedianAF -v VCFFiles/  \n\n" +
				"**************************************************************************************\n");

	}

}
