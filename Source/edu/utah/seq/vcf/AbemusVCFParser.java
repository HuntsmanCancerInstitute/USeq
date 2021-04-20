package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class AbemusVCFParser {

	private File[] tsvFiles;
	private String afInfo = "##INFO=<ID=T_AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=T_DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private String nafInfo = "##INFO=<ID=N_AF,Number=1,Type=Float,Description=\"Allele Frequency for normal\">";
	private String ndpInfo = "##INFO=<ID=N_DP,Number=1,Type=Integer,Description=\"Read depth for normal\">";
	
	public AbemusVCFParser (String[] args) throws Exception {

		processArgs(args);

		for (File tsv: tsvFiles){
			System.out.println(tsv.getName()+"\t");
			BufferedReader in = IO.fetchBufferedReader(tsv);
			String name = Misc.removeExtension(tsv.getName())+".vcf.gz";
			Gzipper out = new Gzipper (new File(tsv.getParentFile(), name));
			String[] f = null;
			String line = null;
			int recNum=1;
			while ((line = in.readLine())!=null) {
				line = line.trim();
				if (line.length()==0) continue;
				if (line.startsWith("group")) {
					//add header
					out.println("##fileformat=VCFv4.1\n##reference=GRCh38");
					out.println(afInfo);
					out.println(dpInfo);
					out.println(nafInfo);
					out.println(ndpInfo);
					out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
				}
				else {
					f = Misc.TAB.split(line);
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
					//   1        2       3   5                   T_AF=10; T_DP=11;N_AF=26; N_DP=27
					out.print(f[1]); out.print("\t");
					out.print(f[2]); out.print("\t");
					out.print(recNum++); out.print("\t");
					out.print(f[3]); out.print("\t");
					out.print(f[5]); 
					out.print("\t.\t.\tT_AF="); out.print(f[10]); 
					out.print(";T_DP="); out.print(f[11]); 
					out.print(";N_AF="); out.print(f[26]); 
					out.print(";N_DP="); out.println(f[27]);
				}
			}
			in.close();
			out.close();
		}
		
		IO.pl("\nDone!");
	}

	

	public static void main(String[] args) throws Exception {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AbemusVCFParser(args);
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
					case 't': tsvFiles = IO.extractFiles(new File(args[++i]), ".tsv"); break;
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
				"**                             Abemus VCF Parser: April 2021                        **\n" +
				"**************************************************************************************\n" +
				"Parses ABEMUS tsv files to a minimal vcf format.\n"+

				"\nRequired Options:\n"+
				"-t Path to an ABEMUS tsv file or directory containing such, xxx.tsv(.gz/.zip OK)\n" +

				"\nExample: java -jar pathToUSeq/Apps/AbemusVcfParser -t TSVFiles/ \n\n"+


				"**************************************************************************************\n");

	}

}
