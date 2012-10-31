
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import edu.utah.seq.data.*;
import util.gen.*;
import java.util.*;
import util.bio.annotation.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class PrbFilter {
	
	//fields
	private File[] directories;
	private Pattern bad = Pattern.compile(".+-.+-.+-.+-.+");
	private Pattern tab = Pattern.compile("\t");
	private double maximumFractionBadBases = 0.33333;
	

	//constructors
	public PrbFilter(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		System.out.println("Parsing and filtering...");

		//for each directory, parse prb files and 
		for (int i=0; i< directories.length; i++){
			double numBadReads = 0;
			double numTotalReads = 0;
			System.out.println("\t"+directories[i]);
			//fetch the files
			File[] prbFiles = IO.extractFiles(directories[i], "_prb.txt.gz");
			if (prbFiles.length == 0) System.err.println("\t\tNo xxx_prb.txt.gz files found! Skipping.");
			//make the results file
			File parsedPrbFile = new File(directories[i], "parsed_prb.txt");
			try {
				PrintWriter out = new PrintWriter (new FileWriter (parsedPrbFile));
				//for each file
				for (int j=0; j< prbFiles.length; j++){
					BufferedReader in = IO.fetchBufferedReader(prbFiles[j]);
					String line;
					String[] bases;
					while ((line = in.readLine()) !=null){
						bases = tab.split(line);
						double numBadBases = 0;
						for (int k=0; k<bases.length; k++){
							Matcher mat = bad.matcher(bases[k]);
							if (mat.matches()) numBadBases++;
						}
						double fractionBad = numBadBases/(double)bases.length;
						if (fractionBad< maximumFractionBadBases) out.println(line);
						else {
							numBadReads++;
							System.out.println("BadRead "+line);
						}
						numTotalReads++;
					}
					in.close();
				}
				out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
			double fractionBad = numBadReads/numTotalReads;
			String fractionBadString = Num.formatNumber(fractionBad, 2);
			System.out.println("\t\t"+fractionBadString+" Fraction bad reads ("+(int)numBadReads+"/"+(int)numTotalReads+")");
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
		new PrbFilter(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': directories = IO.extractFiles(args[++i]); break;
					case 'm': maximumFractionBadBases = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (directories == null || directories.length ==0) Misc.printErrAndExit("\nCouldn't find your directory(s) of xxx_prb.txt.gz files?\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Prb Filter: March 2009                          **\n" +
				"**************************************************************************************\n" +
				"Parses and filters xxx_prb.txt.gz files for reads with many all negative base scores.\n" +
				"Writes the good reads to a single file for each directory.\n" +

				"\nOptions:\n"+
				"-d Comma delimited list of full path directory names containing xxx_prb.txt.gz files to\n" +
				"      parse.\n" +
				"-m Maximum fraction of bad bases in a read, defaults to 0.33\n" +

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/PrbFilter -m 0.25 \n" +
				"     -d /GAIIData/Lane1/,/GAIIData/Lane2/,/OldGAIData/Lane2/ \n\n" +

		"**************************************************************************************\n");

	}	

}
