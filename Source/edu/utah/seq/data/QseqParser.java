package edu.utah.seq.data;
import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.zip.*;
import util.gen.IO;
import java.util.regex.*;

public class QseqParser extends Thread{
	
	//fields
	private boolean complete = false;
	private boolean passed = false;
	private File saveDirectory;
	private ArrayList<File> qseqFiles;
	private boolean pairsPresent = false;
	private boolean printFullHeaders;
	private boolean keepAllReads;
	private static final Pattern name = Pattern.compile("(.+)_qseq.*");
	private static final Pattern tab = Pattern.compile("\\t");
	private static final Pattern period = Pattern.compile("\\.");
	private int numberPassingReads = 0;
	private int numberFailedReads = 0;
	private int numberMalformedQseqLines = 0;

	//constructor
	public QseqParser (File saveDirectory, ArrayList<File> qseqFiles, boolean printFullHeaders, boolean pairsPresent, boolean keepAllReads){
		this.saveDirectory = saveDirectory;
		this.qseqFiles = qseqFiles;
		this.printFullHeaders = printFullHeaders;
		this.pairsPresent = pairsPresent;
		this.keepAllReads = keepAllReads;
	}
	
	//methods
	public int getTotalNumberReads(){
		return numberPassingReads+ numberFailedReads;
	}
	
	public void run(){
		if (pairsPresent) parsePairedReads();
		else parseSingleReads();
	}
	
	public void parseSingleReads(){
		File file = null;
		BufferedReader in = null;
		Writer out = null;
		String line = null;;
		String[] parts;
		StringBuilder sb;	
		try{
			//for each file
			int num = qseqFiles.size();
			for (int i=0; i< num; i++){
				
				//make input stream
				file = qseqFiles.get(i);
				in = IO.fetchBufferedReader(file);
				
				//make output writer
				Matcher mat = name.matcher(file.getName());
				mat.matches();
				String fastqName = mat.group(1) +".fastq.gz";
				File fastqFile = new File (saveDirectory, fastqName);
				out = new OutputStreamWriter ( new GZIPOutputStream ( new FileOutputStream (fastqFile)));
				
				//parse it!
				while ((line = in.readLine()) != null){
					//split on tab and check
					parts = tab.split(line);
					if (parts.length != 11){
						numberMalformedQseqLines++;
						continue;
					}
					//check if it passed filters
					if (keepAllReads == false && parts[10].equals("0")){
						numberFailedReads++;
						continue;
					}
					numberPassingReads++;
					
					//make fastq block
					sb = new StringBuilder();
					
					//build seq header
					sb.append("@");
					String h;
					if (printFullHeaders == false){
						h = Integer.toString(numberPassingReads);
					}
					else {
						h= parts[0] + ":" + parts[2] + ":" + parts[3] + ":" + parts[4] + ":" + parts[5] + "#" + parts[6] + "/" + parts[7];
					}
					sb.append(h);
					sb.append("\n");
					
					//add seq replacing any . with N
					sb.append(period.matcher(parts[8]).replaceAll("N"));
					sb.append("\n");
					
					//build quality header
					sb.append("+");
					if (printFullHeaders) sb.append(h);
					sb.append("\n");
					
					//add quality string
					sb.append(parts[9]);
					sb.append("\n");
					
					out.write(sb.toString());
					
				}
				in.close();
				out.close();
			}
			passed = true;
		} catch (Exception e){
			System.err.println("Error: Problem parsing qseq file -> "+file.getName()+" line -> "+line);
			e.printStackTrace();
			passed = false;
		} finally {
			complete = true;
			try {
				if (in != null) in.close();
				if (out != null) out.close();
			} catch (IOException e) {}
		}
	}

	
	public void parsePairedReads(){
		File file1 = null;
		BufferedReader in1 = null;
		Writer out1 = null;
		String line1 = null;;
		String[] parts1;
		StringBuilder sb1;
		
		File file2 = null;
		BufferedReader in2 = null;
		Writer out2 = null;
		String line2 = null;;
		String[] parts2;
		StringBuilder sb2;		
		
		try{
			//for each pair of files file
			int num = qseqFiles.size();
			for (int i=0; i< num; i+=2){
				
				//make input streams
				file1 = qseqFiles.get(i);
				in1 = IO.fetchBufferedReader(file1);
				file2 = qseqFiles.get(i+1);
				in2 = IO.fetchBufferedReader(file2);
				
				//make output stream
				Matcher mat = name.matcher(file1.getName());
				mat.matches();
				String fastqName1 = mat.group(1) +".fastq.gz";
				File fastqFile1 = new File (saveDirectory, fastqName1);
				out1 = new OutputStreamWriter ( new GZIPOutputStream ( new FileOutputStream (fastqFile1)));
				
				mat = name.matcher(file2.getName());
				mat.matches();
				String fastqName2 = mat.group(1) +".fastq.gz";
				File fastqFile2 = new File (saveDirectory, fastqName2);
				out2 = new OutputStreamWriter ( new GZIPOutputStream ( new FileOutputStream (fastqFile2)));
				
				//parse it!
				while ((line1 = in1.readLine()) != null){
					line2 = in2.readLine();
					if (line2 == null) throw new Exception("Error: File two has less lines than file one?!");
					//split on tab and check
					parts1 = tab.split(line1);
					parts2 = tab.split(line2);
					if (parts1.length != 11 || parts2.length != 11){
						numberMalformedQseqLines++;
						continue;
					}
					//check if both failed
					if (keepAllReads == false && parts1[10].equals("0") && parts2[10].equals("0")){
						numberFailedReads++;
						continue;
					}
					numberPassingReads++;
					
					//make fastq block
					sb1 = new StringBuilder();
					sb2 = new StringBuilder();
					
					//build seq header
					sb1.append("@");
					sb2.append("@");
					String h1;
					String h2;
					if (printFullHeaders == false){
						h1 = Integer.toString(numberPassingReads);
						h2 = h1;
					}
					else {
						h1= parts1[0] + ":" + parts1[2] + ":" + parts1[3] + ":" + parts1[4] + ":" + parts1[5] + "#" + parts1[6] + "/" + parts1[7];
						h2= parts2[0] + ":" + parts2[2] + ":" + parts2[3] + ":" + parts2[4] + ":" + parts2[5] + "#" + parts2[6] + "/" + parts2[7];
					}
					sb1.append(h1);
					sb2.append(h2);
					sb1.append("\n");
					sb2.append("\n");
					
					//add seq replacing any . with N
					sb1.append(period.matcher(parts1[8]).replaceAll("N"));
					sb1.append("\n");
					sb2.append(period.matcher(parts2[8]).replaceAll("N"));
					sb2.append("\n");
					
					//build quality header
					sb1.append("+");
					sb2.append("+");
					if (printFullHeaders) {
						sb1.append(h1);
						sb2.append(h2);
					}
					sb1.append("\n");
					sb2.append("\n");
					
					//add quality string
					sb1.append(parts1[9]);
					sb2.append(parts2[9]);
					sb1.append("\n");
					sb2.append("\n");
					
					out1.write(sb1.toString());
					out2.write(sb2.toString());
					
				}
				//check that file two is out of lines
				line2 = in2.readLine();
				if (line2 != null) throw new Exception("Error: File one has less lines than file two?!");
				
				//close readers
				in1.close();
				in2.close();
				out1.close();
				out2.close();
			}
			passed = true;
		} catch (Exception e){
			System.err.println("Error: Problem parsing qseq file -> "+file1.getName()+" line -> "+line1);
			System.err.println("\t Or problem parsing qseq file -> "+file2.getName()+" line -> "+line2);
			e.printStackTrace();
			passed = false;
		} finally {
			complete = true;
			try {
				//close the streams
				if (out1 != null) out1.close();
				if (out2 != null) out2.close();
				if (in1 != null) in1.close();
				if (in2 != null) in2.close();
			} catch (IOException e) {}
		}
	}

	
	public boolean isComplete() {
		return complete;
	}

	public boolean isPassed() {
		return passed;
	}

	public int getNumberPassingReads() {
		return numberPassingReads;
	}

	public int getNumberFailedReads() {
		return numberFailedReads;
	}

	public int getNumberMalformedQseqLines() {
		return numberMalformedQseqLines;
	}

	public File getSaveDirectory() {
		return saveDirectory;
	}
}
