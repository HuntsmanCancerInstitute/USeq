package edu.utah.seq.barcodes;

import java.io.*;
import java.util.regex.*;
import util.gen.*;

/**Takes paired gzipped fastq files and outputs interlaced fastq which checking for fragment name integrity.
 * @author david.nix@aruplab.com
 **/
public class FastqInterlacer{

	//user defined fields
	private File firstReadFastq;
	private File secondReadFastq;

	//internal fields
	private BufferedReader firstFastqIn;
	private BufferedReader secondFastqIn;
	private String[] first = new String[4];
	private String[] second = new String[4];
	private long lineNumber = 0;

	//constructors
	public FastqInterlacer(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		try {
			makeReadersAndWriters();

			parse();

			closeIO();

		} catch (IOException e) {
			System.err.println("\n\nProblems found parsing fastq data!");
			Misc.printErrAndExit(e.getMessage());
		} 

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		String res = "Done! "+Math.round(diffTime)+" min for "+(int)((double)lineNumber/4.0)+" fastq records.";
		System.err.println(res);


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
			checkAndAssignFastqNames();
			printReads();
		}
	}

	private void printReads() throws IOException {
		for (int i=0; i< 4; i++) System.out.println(first[i]);
		for (int i=0; i< 4; i++) System.out.println(second[i]);	
	}

	private void checkAndAssignFastqNames() {
		//split on tab, sometimes the sample barcode is appended with whitespace divider, might be others; the frag name is always first
		String[] f = Misc.WHITESPACE.split(first[0]);
		String[] s = Misc.WHITESPACE.split(second[0]);

		//check that name is the same from first, second, and barcode
		if (f[0].equals(s[0]) == false) Misc.printErrAndExit("\nError, looks like your first and second fastq names differ? \nFirst\t"+f[0]+"\nSecond\t"+s[0]);

		//make new frag name with appended barcode
		String fragmentName = f[0];

		//make new header line for first
		StringBuilder sb = new StringBuilder(fragmentName);
		sb.append("/1");
		for (int i=1; i< f.length; i++){
			sb.append(" ");
			sb.append(f[i]);
		}
		first[0] = sb.toString();

		//make new header for second
		sb = new StringBuilder(fragmentName);
		sb.append("/2");
		for (int i=1; i< s.length; i++){
			sb.append(" ");
			sb.append(s[i]);
		}
		second[0] = sb.toString();

		//watch out for non + names for other 1/2 of fastq record
		if (first[2].equals("+") == false) first[2] = first[0];
		if (second[2].equals("+") == false) second[2] = second[0];
	}

	public void makeReadersAndWriters() throws IOException{
		//readers
		firstFastqIn = IO.fetchBufferedReader(firstReadFastq);
		secondFastqIn = IO.fetchBufferedReader(secondReadFastq);
	}

	public void closeIO() throws IOException{
		//close readers
		firstFastqIn.close();
		secondFastqIn.close();
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FastqInterlacer(args);
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
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}

			}
		}

		String params = "\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n";
		System.err.println(params);

		//check fields
		if (firstReadFastq == null || firstReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your first read fastq file -> "+firstReadFastq);
		if (secondReadFastq == null || secondReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your second read fastq file -> "+secondReadFastq);
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Fastq Interlacer: April 2016                          **\n" +
				"**************************************************************************************\n" +
				"Takes paired fastq files and writes interlaced/ interleaved fastq to stndOut. \n"+

				"\nOptions:\n"+
				"-f First fastq file, .gz/.zip OK.\n" +
				"-s Second fastq file, .gz/.zip OK.\n" +

				"\nExample: java -Xmx1G -jar pathToUSeq/Apps/FastqInterlacer -f lob_1.fastq.gz\n" +
				"     -s lob_2.fastq.gz | cutadapt | bwa | samblaster .... \n\n" +

				"**************************************************************************************\n");

	}
}
