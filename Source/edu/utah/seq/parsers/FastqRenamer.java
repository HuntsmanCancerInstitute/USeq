package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.gen.*;

/**Takes paired gzipped fastq files and outputs fastq with names prepended with a line row.
 * @author david.nix@aruplab.com
 **/
public class FastqRenamer{

	//user defined fields
	private File firstReadFastq;
	private File secondReadFastq;
	private File saveDirectory;

	//internal fields
	private BufferedReader firstFastqIn;
	private BufferedReader secondFastqIn;
	private Gzipper firstOut;
	private Gzipper secondOut;
	private String[] first = new String[4];
	private String[] second = new String[4];
	private long lineNumber = 0;
	private long fastqCount = 0;
	private String plus = "+";

	//constructors
	public FastqRenamer(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		try {
			makeReadersAndWriters();

			parse();

			closeIO();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			String res = "Done! "+Math.round(diffTime)+" min for "+fastqCount+" fastq records.";
			System.out.println(res);

		} catch (IOException e) {
			System.err.println("\n\nProblems found parsing fastq data!");
			Misc.printErrAndExit(e.getMessage());
		} 




	}

	public boolean loadARead(){
		//attempt to load 4 lines
		int i=0;
		try {
			for (; i<4; i++){
				first[i] = firstFastqIn.readLine();
				second[i] = secondFastqIn.readLine();

				//check all are good
				if (first[i] == null || second[i] == null) throw new Exception();
				lineNumber++;
			}
			return true;
		} catch (Exception e){	
			//on index 0
			if (i!=0) Misc.printErrAndExit("\nError: looks like one of your fastq files is truncated! Only loaded part of a four line fastq record, line number "+lineNumber);
			//all null?
			if (first[i]!=null || second[i]!=null){
				Misc.printErrAndExit("\nError: looks like one of your fastq files is shorter than another? line number "+lineNumber);
			}
			return false;
		}
	}

	public void parse() throws IOException {
		while (loadARead()){
			modifyReads();
			printReads();
		}
	}

	private void modifyReads() {
		fastqCount++;
		
		//replace with fastqCount
		first[0] = "@"+fastqCount;
		second[0] = first[0];
		
		//wipe any name in 3rd line
		first[2] = plus;
		second[2] = plus;
	}

	private void printReads() throws IOException {
		for (int i=0; i< 4; i++) firstOut.println(first[i]);
		for (int i=0; i< 4; i++) secondOut.println(second[i]);	
	}

	public void makeReadersAndWriters() throws IOException{
		//readers
		firstFastqIn = IO.fetchBufferedReader(firstReadFastq);
		secondFastqIn = IO.fetchBufferedReader(secondReadFastq);
		firstOut = new Gzipper(new File(saveDirectory, firstReadFastq.getName()));
		secondOut = new Gzipper(new File(saveDirectory, secondReadFastq.getName()));
	}

	public void closeIO() throws IOException{
		//close readers
		firstFastqIn.close();
		secondFastqIn.close();
		firstOut.close();
		secondOut.close();
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FastqRenamer(args);
	}		



	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': firstReadFastq = new File(args[++i]); break;
					case 's': secondReadFastq = new File(args[++i]); break;
					case 'd': saveDirectory = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}

			}
		}

		String params = "\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n";
		System.out.println(params);

		//check fields
		if (firstReadFastq == null || firstReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your first read fastq file -> "+firstReadFastq);
		if (secondReadFastq == null || secondReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your second read fastq file -> "+secondReadFastq);
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save directory?");
		saveDirectory.mkdirs();
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Fastq Renamer: April 2018                             **\n" +
				"**************************************************************************************\n" +
				"Takes paired fastq files and replaces the header with the record count. \n"+

				"\nOptions:\n"+
				"-f First fastq file, .gz/.zip OK.\n" +
				"-s Second fastq file, .gz/.zip OK.\n" +
				"-d Path to a directory for saving the modified fastq files.\n"+

				"\nExample: java -Xmx1G -jar pathToUSeq/Apps/FastqRenamer -f lob_1.fastq.gz\n" +
				"     -s lob_2.fastq.gz -d UniquifiedFastq/ \n\n" +

				"**************************************************************************************\n");

	}
}
