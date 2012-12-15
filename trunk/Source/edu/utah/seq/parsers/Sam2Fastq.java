
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.data.*;
import edu.utah.seq.data.sam.MalformedSamAlignmentException;
import edu.utah.seq.data.sam.SamAlignment;
import util.gen.*;

import java.util.*;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecord.SAMTagAndValue;


/**Takes an unsorted sam file and the original fastq files and pulls the original fastq data lines.  Assumes the order of the reads is preserved in all files.
 * @author david.nix@hci.utah.edu 
 **/
public class Sam2Fastq{
	//fields
	private File samFile;
	private File firstReadFastq;
	private File secondReadFastq;
	private BufferedReader samIn;
	private BufferedReader firstFastqIn;
	private BufferedReader secondFastqIn;
	private Gzipper firstFastqOut;
	private Gzipper secondFastqOut;
	private Gzipper samOut;
	private int readNameIndex = 0;
	private boolean paired = false;
	private static final Pattern TAB = Pattern.compile("\t");
	private int firstFastqLineNumber = 0;
	private int secondFastqLineNumber = 0;
	private boolean replaceSamSeqWithOriginal = false;

	//constructors

	public Sam2Fastq(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		try {
			
			makeReadersAndWriters();
			
			System.out.print("Parsing");
			if (replaceSamSeqWithOriginal) replaceEm();
			else parseEm();
			System.out.println();

			closeIO();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
			
		} catch (IOException e) {
			System.err.println("\n\nProblems found extracting fastq data!");
			System.err.println(e.getMessage());
		} 
		
		
		
	}
	
	public void replaceEm() throws IOException{
		
		String line;
		String oldReadNameFirst = "";
		String oldReadNameSecond = "";
		String oldFirstSeq = null;
		String oldSecondSeq = null;
		//for each alignment line
		int dotCounter = 0;
		int numBadLines = 0;
		while ((line = samIn.readLine()) !=null){
			if (++dotCounter > 1000000){
				System.out.print(".");
				dotCounter = 0;
			}
			line = line.trim();

			//skip blank lines
			if (line.length() == 0) continue;

			//print header lines
			if (line.startsWith("@")) {
				samOut.println(line);
				continue;
			}
			
			SamAlignment sa;
			try {
				sa = new SamAlignment(line, false);
			} catch (MalformedSamAlignmentException e) {
				System.out.println("\nSkipping malformed sam alignment -> "+e.getMessage());
				if (numBadLines++ > 100) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
				continue;
			}
			
			//fetch read name
			String readName = sa.getName();
System.out.println("\nNew alignment for "+readName + "  fp? "+sa.isFirstPair()+"  rs? "+sa.isReverseStrand());
System.out.println("OriSam "+line);
			
			//is it first?
			if (sa.isFirstPair()){
				//new?
				if (readName.equals(oldReadNameFirst) == false) {
					oldReadNameFirst = readName;
					oldFirstSeq = fetchFirstSeq(oldReadNameFirst);
					oldFirstSeq = sa.processOriginalSequence(oldFirstSeq);
				}
				else System.out.println("\tAlready parsed");
				sa.setSequence(oldFirstSeq);
			}
			
			//second read
			else {
				//new?
				if (readName.equals(oldReadNameSecond) == false) {
					oldReadNameSecond = readName;
					oldSecondSeq = fetchSecondSeq(oldReadNameSecond);
					oldSecondSeq = sa.processOriginalSequence(oldSecondSeq);
				}
				else System.out.println("\tAlready parsed");
				sa.setSequence(oldSecondSeq);
			}
			
			//print it
System.out.println("FinSam "+line);			
			samOut.println(sa.toStringNoMDField());
			
		}
		
	}

	public void parseEm() throws IOException{
		
		String line;
		String oldReadName = "";
		//for each alignment line
		int dotCounter = 0;
		int numLines = 0;
		while ((line = samIn.readLine()) !=null){
			numLines++;
			//print dot?
			if (++dotCounter > 50000){
				System.out.print(".");
				dotCounter = 0;
			}
			
			//skip comment lines
			if (line.startsWith("@")) continue;
			
			//fetch read name
			String readName = "@"+(TAB.split(line)[readNameIndex]);
			
			//check it's not the name
			if (readName.equals(oldReadName)) continue;
			oldReadName = readName;
			
			//find and print fastq data
			findPrintFirstFastq(readName);
			if (paired){
				//find and print fastq data
				findPrintSecondFastq(readName);
			}
		}
		
	}
	
	public String fetchFirstSeq(String readName) throws IOException{
		String line;
		String headerLine = "@"+readName;
		//walk through file
		while ((line = firstFastqIn.readLine()) != null){
			if (line.startsWith(headerLine)){
				//read in seq
				line = firstFastqIn.readLine();
				firstFastqIn.readLine();
				firstFastqIn.readLine();
				return line;
			}
		}
		//should never end
		throw new IOException("\nError: reached the end of the first fastq file without finding the read name -> "+readName);
	}
	public String fetchSecondSeq(String readName) throws IOException{
		String line;
		String headerLine = "@"+readName;
		//walk through file
		while ((line = secondFastqIn.readLine()) != null){
			if (line.startsWith(headerLine)){
				//read in seq
				line = secondFastqIn.readLine();
				secondFastqIn.readLine();
				secondFastqIn.readLine();
				return line;
			}
		}
		//should never end
		throw new IOException("\nError: reached the end of the second fastq file without finding the read name -> "+readName);
	}
	
	public void findPrintFirstFastq(String readName) throws IOException{
		String line;
		//walk through file
		while ((line = firstFastqIn.readLine()) != null){
			if (line.startsWith(readName)){
				firstFastqOut.println(line);
				firstFastqOut.println(firstFastqIn.readLine());
				firstFastqOut.println(firstFastqIn.readLine());
				firstFastqOut.println(firstFastqIn.readLine());
				firstFastqLineNumber+=4;
				return;
			}
			firstFastqLineNumber++;
		}
		//should never end
		throw new IOException("\nError: reached the end of the first fastq file without finding the read name -> "+readName);
	}
	
	public void findPrintSecondFastq(String readName) throws IOException{
		String line;
		//walk through file
		while ((line = secondFastqIn.readLine()) != null){
			if (line.startsWith(readName)){
				secondFastqOut.println(line);
				secondFastqOut.println(secondFastqIn.readLine());
				secondFastqOut.println(secondFastqIn.readLine());
				secondFastqOut.println(secondFastqIn.readLine());
				secondFastqLineNumber+=4;
				//check line numbers
				if (secondFastqLineNumber != firstFastqLineNumber) throw new IOException("\nError: line numbers between the first and the second fastq files differ for read name -> "+readName);
				return;
			}
			secondFastqLineNumber++;
		}
		//should never end
		throw new IOException("\nError: reached the end of the second fastq file without finding the read name -> "+readName);
	}
	
	
	public void makeReadersAndWriters() throws IOException{
			//readers
			samIn = IO.fetchBufferedReader(samFile);
			firstFastqIn = IO.fetchBufferedReader(firstReadFastq);
			if (paired) secondFastqIn = IO.fetchBufferedReader(secondReadFastq);
			
			//writers
			if (replaceSamSeqWithOriginal){
				File parsedSam = new File ("oriSeq_"+samFile);
				samOut = new Gzipper(parsedSam);
			}
			else {
				File firstParsed = new File (firstReadFastq.getParentFile(), "parsed_"+firstReadFastq.getName());
				firstFastqOut = new Gzipper(firstParsed);
				if (paired){
					File secondParsed = new File (secondReadFastq.getParentFile(), "parsed_"+secondReadFastq.getName());
					secondFastqOut = new Gzipper(secondParsed);
				}
			}
			
	}
	
	public void closeIO() throws IOException{
			samIn.close();
			firstFastqIn.close();
			
			if (replaceSamSeqWithOriginal) {
				samOut.close();
				if (paired) secondFastqIn.close();
			}
			else {
				firstFastqOut.close();
				if (paired){
					secondFastqIn.close();
					secondFastqOut.close();
				}
			}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Sam2Fastq(args);
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
					case 'a': samFile = new File(args[++i]); break;
					case 'f': firstReadFastq = new File(args[++i]); break;
					case 's': secondReadFastq = new File(args[++i]); paired = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (samFile == null || samFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your alignment file -> "+samFile);
		if (firstReadFastq == null || firstReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your first read fastq file -> "+firstReadFastq);
		if (paired){
			if (secondReadFastq == null || secondReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your second read fastq file -> "+secondReadFastq);
		}

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Sam 2 Fastq: March 2012                            **\n" +
				"**************************************************************************************\n" +
				"Extracts the original Illumina fastq data from single or paired end sam alignments.\n" +
				"Assumes alignments and reads are in the same order. In novoalign, set -oSync .\n"+

				"\nOptions:\n"+
				"-a Sam alignment txt file, full path, .gz/.zip OK.\n"+
				"-f First read fastq file, ditto.\n" +
				"-s (Optional) Second read fastq file, from paired read sequencing, ditto.\n" +


				"\nExample: java -Xmx1G -jar pathToUSeq/Apps/Sam2Fastq -a /SAM/unaligned.sam.gz -f \n" +
				"     /Fastq/X1_110825_SN141_0377_AD06YNACXX_1_1.txt.gz -s \n" +
				"     /Fastq/X1_110825_SN141_0377_AD06YNACXX_1_2.txt.gz\n\n" +

		"**************************************************************************************\n");

	}
}
