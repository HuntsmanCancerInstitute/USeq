
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.gen.*;


/**
 * @author david.nix@hci.utah.edu 
 **/
public class CompareParsedAlignments {
	
	//fields
	private File alleleFile1;
	private File alleleFile2;
	private AlleleWithAlignments[] alleles1;
	private AlleleWithAlignments[] alleles2;
	private File tempDirectory;
	private File rApp = new File ("/usr/bin/R");
	
	//constructors
	public CompareParsedAlignments(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		//check lengths
		if (alleles1.length != alleles2.length) Misc.printErrAndExit("Allele array lengths differ, did you run the ParseIntersectingAlignments app on the same allele table? Aborting!\n");
		
		//make maxtricies
		int[][] baseCounts1 = new int[alleles1.length][4];
		int[][] baseCounts2 = new int[alleles2.length][4];
		
		for (int i=0; i< alleles1.length; i++){
			baseCounts1[i] = alleles1[i].getGATCCounts();
			baseCounts2[i] = alleles2[i].getGATCCounts();
		}
		
		//calc fisher exact p-value approximations
		double[] pvalues = Num.fisherTest(baseCounts1, baseCounts2, tempDirectory, rApp);
		
		//print
		for (int i=0; i< alleles1.length; i++){
			System.out.println(pvalues[i]+"\t"+alleles1[i]+"\t"+alleles2[i]);
		}
		

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
		public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CompareParsedAlignments(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': alleleFile1 = new File(args[++i]); break;
					case 'b': alleleFile2 = new File(args[++i]); break;
					case 'r': rApp = new File(args[++i]); break;
					case 'd': tempDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (alleleFile1 == null || alleleFile1 == null) Misc.printErrAndExit("\nCouldn't find one or both of your xxx.alleles files?\n");
		alleles1 = (AlleleWithAlignments[])IO.fetchObject(alleleFile1);
		alleles2 = (AlleleWithAlignments[])IO.fetchObject(alleleFile2);
		tempDirectory.mkdir();
		if (rApp == null || rApp.exists() == false) Misc.printErrAndExit("Cannot find R?!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Compare Parsed Alignments: Nov 2009                      **\n" +
				"**************************************************************************************\n" +
				"Compares two parsed alignments for a common distribution of snps using R's Fisher's\n" +
				"Exact. Run the ParseIntersectingAlignments with the same snp table first.\n" +

				"\nOptions:\n"+
				"-a Full path file name for the first xxx.alleles file.\n" +
				"-b Full path file name for the first xxx.alleles file.\n" +
				"-d Full path directory name for writing temporary files.\n"+
				"-r Full path file name for R, defaults to '/usr/bin/R'\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/CompareParsedAlignments. \n" +
				"     -a /SeqData/lymphSNPs.alleles -b /SeqData/normalSNPs.alleles -b /temp/\n\n" +

		"**************************************************************************************\n");

	}	

}
