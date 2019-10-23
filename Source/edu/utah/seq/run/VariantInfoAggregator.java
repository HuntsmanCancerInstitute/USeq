package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class VariantInfoAggregator {
	
	private File jobDir = null;
	private File resultSpreadsheet = null;

	public VariantInfoAggregator (String[] args){
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);
			
			parseAndPrint();
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
	}
	

	private void parseAndPrint() throws IOException {
		File[] dirs = IO.extractOnlyDirectories(jobDir);
		if (dirs.length == 0) Misc.printErrAndExit("\nFailed to find any sample dirs in "+jobDir);
		
		PrintWriter out = new PrintWriter(resultSpreadsheet);
		out.println(VariantInfoSample.header);
		//for each sample dir
		IO.pl("Processing...");
		for (File dir: dirs) {
			IO.pl("\t"+dir);
			VariantInfoSample s = new VariantInfoSample(dir);
			out.println(s.toString());
		}
		
		out.close();
		
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VariantInfoAggregator(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': resultSpreadsheet = new File(args[++i]); break;
					case 'j': jobDir = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (jobDir == null || resultSpreadsheet == null) Misc.printErrAndExit("\nPlease provide both a directory containing TNRunner results directories and the name of a spreadsheet results file.\n"); 
		jobDir.mkdirs();
		
	}


	public static void printDocs(){
		
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                            VariantInfoAggregator : Sept 2019                     **\n" +
				"**************************************************************************************\n" +
				"Parses and aggregates somatic and germline variant calls from many samples precessed\n"+
				"by the USeq TNRunner application.\n"+

				"\nOptions:\n"+
				"-j Job directory containing directories for each patient results TNRunner. \n" +
				"-s Spreadsheet file name in which to write the results.\n"+
				
				
                "\nExample: java -jar -Xmx2G ~/USeqApps/VariantInfoAggregator -j AJobs/ -s \n"+
                "    aggregatedVariantInfo.txt\n\n"+


				"**************************************************************************************\n");
	}
}
