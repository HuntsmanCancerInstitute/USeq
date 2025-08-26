package edu.utah.seq.barcodes;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;


public class UmiWorker extends Thread{

	private FastqUmiTagger fut = null;
	private ArrayList<String[]> reads = null;
	private ArrayList<String> processedReads = null;
	private boolean failed = false;
	private boolean loading3Reads;
	private static final Pattern badEnding = Pattern.compile(".+/\\d");
	private String[] first = new String[4];
	private String[] second = new String[4];
	private String[] umi = new String[4];
	private int umiLength = -1;
	private int bpToTrim = -1;
	private int maxNsInUMI = -1;
	private long numUmiNFailures = 0;

	//constructor
	public UmiWorker(FastqUmiTagger fut) throws IOException{
		this.fut = fut;
		loading3Reads = fut.getUmiFastq() != null;
		if (loading3Reads) {
			reads = new ArrayList<String[]>(fut.getChunkSize()*3);
		}
		else {
			reads = new ArrayList<String[]>(fut.getChunkSize()*2);
			umiLength = fut.getUmiLength();
			bpToTrim = fut.getBpToTrim();
		}
		processedReads = new ArrayList<String>(fut.getChunkSize()*4*2);
		maxNsInUMI = fut.getMaxNsInUMI();
	}

	public void run(){

		try {
			if (loading3Reads) {
				while (fut.load3Reads(reads)) {
					//process reads
					int num = reads.size();
					for (int i=0; i< num; i++) {
						first = reads.get(i++);
						second = reads.get(i++);
						umi = reads.get(i);
						
						if (checkUMI(umi[1])) {
							checkAndAssignFastqNames();
							saveRecords();
						}
					}
					
					fut.printReads(processedReads);
					
					//clear the worker data
					reads.clear();
					processedReads.clear();
				}
			}
			else {
				while (fut.load2Reads(reads)) {
					//process reads
					int num = reads.size();
					for (int i=0; i< num; i++) {
						first = reads.get(i++);
						second = reads.get(i);
						//for in line
						if (extractInlineUMI()) {
							checkAndAssignFastqNames();
							saveRecords();
						}
					}
					fut.printReads(processedReads);
					//clear the worker data
					reads.clear();
					processedReads.clear();
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			failed = true;
		}
	}
	
	private void saveRecords() {
		for (int i=0; i< 4; i++) processedReads.add(first[i]);
		for (int i=0; i< 4; i++) processedReads.add(second[i]);
	}

	/**For parsing IDTs 3N,skip2, insert from read 1 and 2 and creating the 6mer UMI*/
	private boolean extractInlineUMI() {
		//line 0, name
		umi[0] = first[0];
		
		//line 1, seq
		String readOneSeq = first[1].substring(0, umiLength);
		String readTwoSeq = second[1].substring(0, umiLength);
		umi[1] = readOneSeq + readTwoSeq;
		
		if (checkUMI(umi[1])==false) return false;
		
		//change first and second seqs
		first[1] = first[1].substring(bpToTrim);
		second[1] = second[1].substring(bpToTrim);
		
		//line 2, name
		umi[2] = first[2];
		
		//line 3, qual
		String readOneQual = first[3].substring(0, umiLength);
		String readTwoQual = second[3].substring(0, umiLength);
		umi[3] = readOneQual + readTwoQual;

		//change first and second quals	
		first[3] = first[3].substring(bpToTrim);
		second[3] = second[3].substring(bpToTrim);
		
		return true;
	}
	
	private boolean checkUMI(String umi) {
		//look for Ns?
		if (maxNsInUMI != -1) {
			
			int numNs = 0;
			for (int x=0; x< umi.length(); x++) {
				if (umi.charAt(x) == 'N') numNs++;
			}
			//Misc.printErrAndExit("Checking UMI "+numNs+"\t"+umi);
			if (numNs > maxNsInUMI) {
				numUmiNFailures++;
				return false;
			}
		}
		return true;
	}
	
	private void checkAndAssignFastqNames() {
		//split on tab, sometimes the sample umi is appended with whitespace divider, might be others; the frag name is always first
		String[] f = Misc.WHITESPACE.split(first[0]);
		String[] s = Misc.WHITESPACE.split(second[0]);
		String[] b = Misc.WHITESPACE.split(umi[0]);
		
		//check that name is the same from first, second, and umi
		if (f[0].equals(s[0]) == false) {
			IO.el("\nError, looks like your first and second fastq names differ? \nFirst\t"+f[0]+"\nSecond\t"+s[0]);
			failed = true;
		}
		if (f[0].equals(b[0]) == false) {
			IO.el("\nError, looks like your first and umi fastq names differ? \nFirst\t"+f[0]+"\nUmi\t"+b[0]);
			failed = true;
		}
		
		//make new frag name with appended umi
		StringBuilder fn = new StringBuilder(f[0]);
		fn.append(":BMF:");
		fn.append(umi[1]);
		
		fn.append(checkForBadEnds(umi[3]));
		String fragmentName = fn.toString();
		
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

	public boolean isFailed() {
		return failed;
	}

	public long getNumUmiNFailures() {
		return numUmiNFailures;
	}

	
}
