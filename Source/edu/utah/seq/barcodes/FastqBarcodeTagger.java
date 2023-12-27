
package edu.utah.seq.barcodes;

import java.io.*;
import java.util.regex.*;
import util.gen.*;

/**Takes three fastq files (paired end reads and a third containing barcodes), appends the barcode and quality to the fastq header.
 * Can also parse inline barcodes and trim the fastq seq and qual appropriately
 * @author david.nix@hci.utah.edu
 **/
public class FastqBarcodeTagger{
	
	//user defined fields
	private File firstReadFastq;
	private File secondReadFastq;
	private File barcodeFastq;
	private File resultsDirectory;
	private boolean interlace = false;
	private boolean appendLineNumber = false;
	private int barcodeLength = -1;
	private int bpToTrim = -1;
	
	//internal fields
	private BufferedReader firstFastqIn;
	private BufferedReader secondFastqIn;
	private BufferedReader barcodeFastqIn;
	private Gzipper firstFastqOut;
	private Gzipper secondFastqOut;
	private String[] first = new String[4];
	private String[] second = new String[4];
	private String[] barcode = new String[4];
	private long lineNumber = 0;
	private boolean extractInline = false;
	private Pattern badEnding = Pattern.compile(".+/\\d");

	//constructors
	public FastqBarcodeTagger(String[] args){
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
		if (interlace == false) System.out.println(res);
		else System.err.println(res);
		
		
	}
	
	public boolean loadARead(){
		//attempt to load 4 lines
		int i=0;
		try {
			for (; i<4; i++){
				first[i] = firstFastqIn.readLine();
				second[i] = secondFastqIn.readLine();
				barcode[i] = barcodeFastqIn.readLine();
				//check all are good
				if (first[i] == null || second[i] == null || barcode[i] == null) throw new Exception();
				lineNumber++;
			}
			return true;
		} catch (Exception e){	
			//on index 0
			if (i!=0) Misc.printErrAndExit("\nError: looks like one of your fastq files is truncated! Only loaded part of a four line fastq record, line number "+lineNumber);
			//all null?
			if (first[i]!=null || second[i]!=null || barcode[i]!=null){
					Misc.printErrAndExit("\nError: looks like one of your fastq files is shorter than another? line number "+lineNumber);
			}
			return false;
		}
	}
	
	public boolean loadAndExtractARead(){
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
			//extract the umi
			extract3MerUMI();
			
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
	

	
	/**For parsing IDTs 3N,skip2,insert from read 1 and 2 and creating the 6mer UMI*/
	private void extract3MerUMI() {
		//line 0, name
		barcode[0] = first[0];
		
		//line 1, seq
		//String readOneSeq = first[1].substring(0, 3);
		//String readTwoSeq = second[1].substring(0, 3);
		String readOneSeq = first[1].substring(0, barcodeLength);
		String readTwoSeq = second[1].substring(0, barcodeLength);
		barcode[1] = readOneSeq + readTwoSeq;
		
		//change first and second seqs
		//first[1] = first[1].substring(5);
		//second[1] = second[1].substring(5);
		first[1] = first[1].substring(bpToTrim);
		second[1] = second[1].substring(bpToTrim);
		
		//line 2, name
		barcode[2] = first[2];
		
		//line 3, qual
		//String readOneQual = first[3].substring(0, 3);
		//String readTwoQual = second[3].substring(0, 3);
		String readOneQual = first[3].substring(0, barcodeLength);
		String readTwoQual = second[3].substring(0, barcodeLength);
		barcode[3] = readOneQual + readTwoQual;

		//change first and second quals
		//first[3] = first[3].substring(5);
		//second[3] = second[3].substring(5);	
		first[3] = first[3].substring(bpToTrim);
		second[3] = second[3].substring(bpToTrim);
	}

	public void parse() throws IOException {		
		if (extractInline == false){
			while (loadARead()){
				checkAndAssignFastqNames();
				printReads();
			}
		}
		else {
			while (loadAndExtractARead()){
				checkAndAssignFastqNames();
				printReads();
			}
		}
	}
	
	private void printReads() throws IOException {
		if (interlace){
			for (int i=0; i< 4; i++) System.out.println(first[i]);
			for (int i=0; i< 4; i++) System.out.println(second[i]);
			
		}
		else {
			for (int i=0; i< 4; i++){
				firstFastqOut.println(first[i]);
				secondFastqOut.println(second[i]);
			}
		}
	}

	private void checkAndAssignFastqNames() {
		//split on tab, sometimes the sample barcode is appended with whitespace divider, might be others; the frag name is always first
		String[] f = Misc.WHITESPACE.split(first[0]);
		String[] s = Misc.WHITESPACE.split(second[0]);
		String[] b = Misc.WHITESPACE.split(barcode[0]);
		
		//check that name is the same from first, second, and barcode
		if (f[0].equals(s[0]) == false) Misc.printErrAndExit("\nError, looks like your first and second fastq names differ? \nFirst\t"+f[0]+"\nSecond\t"+s[0]);
		if (f[0].equals(b[0]) == false) Misc.printErrAndExit("\nError, looks like your first and barcode fastq names differ? \nFirst\t"+f[0]+"\nBarcode\t"+b[0]);
		
		//make new frag name with appended barcode
		StringBuilder fn = new StringBuilder(f[0]);
		if (appendLineNumber) {
			fn.append(":");
			fn.append(lineNumber);
		}
		fn.append(":BMF:");
		fn.append(barcode[1]);
		
		fn.append(checkForBadEnds(barcode[3]));
		String fragmentName = fn.toString();
		
		//make new header line for first
		StringBuilder sb = new StringBuilder(fragmentName);
		if (interlace) sb.append("/1");
		for (int i=1; i< f.length; i++){
			sb.append(" ");
			sb.append(f[i]);
		}
		first[0] = sb.toString();

		//make new header for second
		sb = new StringBuilder(fragmentName);
		if (interlace) sb.append("/2");
		for (int i=1; i< s.length; i++){
			sb.append(" ");
			sb.append(s[i]);
		}
		second[0] = sb.toString();

		//watch out for non + names for other 1/2 of fastq record
		if (first[2].equals("+") == false) first[2] = first[0];
		if (second[2].equals("+") == false) second[2] = second[0];
	}

	/* Must swap /0-9 with 00-9 otherwise bwa strips off the /0 /1 /2 etc
	 */
	private String checkForBadEnds(String seq) {
		Matcher mat = badEnding.matcher(seq);
		if (mat.matches()) {
			int len = seq.length();
			return seq.substring(0, len-2)+"0"+seq.substring(len-1);
		}
		return seq;
	}

	public void makeReadersAndWriters() throws IOException{
		//readers
		firstFastqIn = IO.fetchBufferedReader(firstReadFastq);
		secondFastqIn = IO.fetchBufferedReader(secondReadFastq);
		if (extractInline == false) barcodeFastqIn = IO.fetchBufferedReader(barcodeFastq);

		//writers
		if (interlace == false){
			firstFastqOut = new Gzipper(new File (resultsDirectory, "barcoded_"+firstReadFastq.getName()));
			secondFastqOut = new Gzipper(new File (resultsDirectory, "barcoded_"+secondReadFastq.getName()));	
		}
	}

	public void closeIO() throws IOException{
		//close readers
		firstFastqIn.close();
		secondFastqIn.close();
		if (extractInline == false) barcodeFastqIn.close();
		//close writers
		if (interlace == false){
			firstFastqOut.close();
			secondFastqOut.close();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FastqBarcodeTagger(args);
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
					case 'b': barcodeFastq = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'l': barcodeLength = Integer.parseInt(args[++i]); break;
					case 't': bpToTrim = Integer.parseInt(args[++i]); break;
					case 'i': interlace = true; break;
					case 'a': appendLineNumber = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		String params = "\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n";
		if (interlace == false) System.out.println(params);
		else System.err.println(params);
		
		//check fields
		if (firstReadFastq == null || firstReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your first read fastq file -> "+firstReadFastq);
		if (secondReadFastq == null || secondReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your second read fastq file -> "+secondReadFastq);
		
		//barcode file?
		if (barcodeFastq != null) {
			extractInline = false;
			if (barcodeLength != -1 || bpToTrim != -1) Misc.printErrAndExit("\nError: cannot provide a barcode file AND set inline barcode length info");
			if (barcodeFastq.canRead()==false) Misc.printErrAndExit("\nError: cannot find or read your barcode fastq file -> "+barcodeFastq);
		}
		//must be inline
		else {
			extractInline = true;
			if (barcodeLength == -1) Misc.printErrAndExit("\nError: if a barcode file is not provided, then please set a barcode length value for inline UMI extraction.");
			if (bpToTrim == -1) bpToTrim = barcodeLength;
		}
		
		if (interlace == false) {
			if (resultsDirectory != null) resultsDirectory.mkdirs();
			else resultsDirectory = firstReadFastq.getParentFile();
		}
		else {
			if (resultsDirectory != null) Misc.printErrAndExit("\nError: it doesn't make sense to designate a results directory when writing output to stdout!");
		}
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Fastq Barcode Tagger: Oct 2023                          **\n" +
				"**************************************************************************************\n" +
				"Takes 2 or 3 fastq files (paired end reads and possibly a third containing unique \n"+
				"molecular indexes/ barcodes), appends the barcode and quality to the fastq header, and\n"+
				"writes out the modified records. For inline UMI datasets, provide a length to extract\n"+
				"the UMI from beginning of each fastq.\n"+

				"\nOptions:\n"+
				"-f First fastq file, .gz/.zip OK.\n" +
				"-s Second fastq file, .gz/.zip OK.\n" +
				"-b Barcode fastq file, .gz/.zip OK, or set -e\n" + 
				"-l For inline, length of barcode to parse from the beginning of each fastq file.\n"+
				"-t For inline, length to trim the fastqs, defaults to -l\n"+
				"-i Write interlaced fastq to stdout for direct piping to other apps\n"+
				"-r Directory to save the modified fastqs, defaults to the parent of -f\n"+
				"-a Append the line number to the read name to uniquify.\n"+

				"\nExample: java -Xmx1G -jar pathToUSeq/Apps/FastqBarcodeTagger -f lob_1.fastq.gz\n" +
				"     -s lob_2.fastq.gz -b lob_barcode.fastq.gz -i | bwa mem -p /ref/hg19.fa \n\n" +

		"**************************************************************************************\n");

	}
}
