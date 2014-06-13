package edu.utah.ames.bioinfo;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import util.gen.Misc;
/**
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import util.gen.IO;
import util.gen.Num;
 **/
public class MeasureContaminatingSeqs {

	//fields
	private static String inputFastq;
	private static String zooDB = "/tomato/dev/app/contSeqs/zooBWAdb/zooBwaIdx";
	private static String bwaPath = "/tomato/dev/app/bwa/0.7.5/bwa";
	private static int splitLines;
	private static int totalReads;
	private String splitFastq;
	private static String outputSam;
	private static String outName;
	private String MATCHMOUSE = "\\tchrMmouse\\t";
	private String MATCHHUMAN = "\\tchrMhuman\\t";
	private String MATCHECOLI = "\\tchrAllEColiK12\\t";
	private String MATCHSALMON = "\\tchrMsalmon\\t";
	private String MATCHZFISH = "\\tchrMzebrafish\\t";
	private String MATCHWORM = "\\tchrMworm\\t";
	private String MATCHFLY = "\\tchrMfly\\t";
	private String MATCHALLM = "\\t(chr[MA]\\w+)\\t";
	
	private int totalZooReads = 0; //only reads mapped to any of the chrM seqs included in the index
	private int totalReadsSplit = 0; //both mapped and unmapped reads
	private int humanReads = 0;
	private int mouseReads = 0;
	private int salmonReads = 0;
	private int zfishReads = 0;
	private int flyReads = 0;
	private int wormReads = 0;
	private int ecoliReads = 0;
	private int percentMappedReadsAlignedToX = 0;
	private int percentZooReads = 0;
	private int totalReadsAll = 0;

	public static void main(String[] args) throws Exception {
		//check for args
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		MeasureContaminatingSeqs mcs = new MeasureContaminatingSeqs(args);
		//mcs.getFileLines();
		mcs.makeAligns();
	}

	//constructor
	public MeasureContaminatingSeqs(String[] args) {
		//process args
		processArgs(args);
	}

	/**
	 * This method calls bwa to make the alignments and convert to sams
	 * @throws Exception
	 */
	public void makeAligns() throws Exception {

		//make reader for input fastq
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream
				(new FileInputStream(inputFastq))));

		//get #reads in original fastq file
		String line;
		long z = 0; //total reads in input file
		while ((line = br.readLine()) != null) {
			z++;
			this.setTotalReads(z/4);
		}
		br.close();
		
		//set bwa output file name
		//String dest = Misc.removeExtension(inputFastq) + "." + Integer.toString(numInputReads) + "Reads.MCS.bwa";
		String dest = Misc.removeExtension(inputFastq) + ".MCS.bwa";
		
		//build process for execution
		//ProcessBuilder pb1 = new ProcessBuilder(bwaPath, "aln", zooDB, this.getSplitFastq(), dest); 
		ProcessBuilder pb1 = new ProcessBuilder(bwaPath, "aln", "-t 10", zooDB, inputFastq, dest);

		//System.out.println("Aligning against zoo contamination db...");
		Process p1 = pb1.start();

		//make buffered input and output streams
		BufferedInputStream bis1 = new BufferedInputStream(p1.getInputStream());
		BufferedOutputStream bos1 = new BufferedOutputStream(new FileOutputStream(dest));

		//create buffer to read in chunks
		byte[] buffer = new byte[1024*1024*10];
		int n = -1;

		//write output
		while ((n = bis1.read(buffer)) != -1) {
			bos1.write(buffer,0,n);
		}

		//wait for process to finish
		int val = p1.waitFor();

		//close buffered streams
		bos1.close();
		bis1.close();

		//catch error stream if something weird happens
		if (val != 0) {
			System.out.println("Something happened with the alignments!");
			BufferedReader br1 = new BufferedReader(new InputStreamReader(p1.getErrorStream()));
			String line1 = null;
			while ((line1 = br1.readLine()) != null) {
				System.out.println(line1);
			}
			System.exit(1);
		}

		//set bwa sam file output
		//outputSam = Misc.removeExtension(inputFastq) + "." + Integer.toString(numInputReads) + "Reads.MCS.sam";
		outputSam = Misc.removeExtension(inputFastq) + ".MCS.sam";


		//ProcessBuilder pb2 = new ProcessBuilder(bwaPath, "samse", zooDB, dest, this.getSplitFastq(),
		//	outputSam);
		ProcessBuilder pb2 = new ProcessBuilder(bwaPath, "samse", zooDB, dest, inputFastq, outputSam);
		System.out.println("Converting alignment to sam format...");
		Process p2 = pb2.start();

		BufferedInputStream bis2 = new BufferedInputStream(p2.getInputStream());
		BufferedOutputStream bos2 = new BufferedOutputStream(new FileOutputStream(outputSam));
		
		byte[] buf = new byte[1024*1024*10];
		int m = -1;

		while ((m = bis2.read(buf)) != -1) {
			bos2.write(buf,0,m);
		}

		int v = p2.waitFor();
		bos2.close();
		bis2.close();

		if (v != 0) {
			System.out.println("Something happened converting alignments to sam format!");
			BufferedReader br2 = new BufferedReader(new InputStreamReader(p2.getErrorStream()));
			String line2 = null;
			while ((line2 = br2.readLine()) != null) {
				System.out.println(line2);
			}
			System.exit(1);
		}
		//parse output
		parseData();
	}

	/**
	 * This method generates the contamination stats file
	 * @throws Exception
	 */
	public void parseData() throws Exception {
		
		//System.out.println("Parsing the sam file...");
		//outName = (Misc.removeExtension(this.getSplitFastq()) +  ".xls");
		outName = Misc.removeExtension(inputFastq) + ".MCS.xls";

		//create buffered reader to read input sam file
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(outputSam)));

		String line;

		//output stats writer
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outName)));

		while ((line = br.readLine()) != null) {

			//skip header lines in sam file
			if (line.contains("@")) {
				continue;
			}

			//look for contaminating reads
			Pattern p0 = Pattern.compile(MATCHMOUSE);
			Matcher m0 = p0.matcher(line);
			if (m0.find()) {
				mouseReads++;
				this.setMouseReads(mouseReads);
			}
			Pattern p1 = Pattern.compile(MATCHHUMAN);
			Matcher m1 = p1.matcher(line);
			if (m1.find()) {
				humanReads++;
				this.setHumanReads(humanReads);
			}
			Pattern p2 = Pattern.compile(MATCHSALMON);
			Matcher m2 = p2.matcher(line);
			if (m2.find()) {
				salmonReads++;
				this.setSalmonReads(salmonReads);
			}
			Pattern p3 = Pattern.compile(MATCHECOLI);
			Matcher m3 = p3.matcher(line);
			if (m3.find()) {
				ecoliReads++;
				this.setEcoliReads(ecoliReads);
			}
			Pattern p4 = Pattern.compile(MATCHZFISH);
			Matcher m4 = p4.matcher(line);
			if (m4.find()) {
				zfishReads++;
				this.setZfishReads(zfishReads);
			}
			Pattern p5 = Pattern.compile(MATCHWORM);
			Matcher m5 = p5.matcher(line);
			if (m5.find()) {
				wormReads++;
				this.setWormReads(wormReads);
			}
			Pattern p6 = Pattern.compile(MATCHFLY);
			Matcher m6 = p6.matcher(line);
			if (m6.find()) {
				flyReads++;
				this.setFlyReads(flyReads);
			}
			Pattern p7 = Pattern.compile(MATCHALLM);
			Matcher m7 = p7.matcher(line);
			if (m7.find()) {
				totalZooReads++;
				this.setTotalZooReads(totalZooReads);
			}
		}
		/**
		//write output
		bw.write(String.format("%s\nMouse reads", mouseZooReads));
		bw.write(String.format("%s\nHuman reads", humanZooReads));
		bw.write(String.format("%s\nE. coli reads", ecoliZooReads));
		bw.write(String.format("%s\nFly reads", flyZooReads));
		bw.write(String.format("%s\nWorm reads", wormZooReads));
		bw.write(String.format("%s\nSalmon reads", salmonZooReads));
		bw.write(String.format("%s\nZebrafish reads", zfishZooReads));
		bw.write(String.format("%s\nTotal reads", totalReadsSplit));
		bw.write(String.format("%s\nTotal contaminating reads", totalZooReads));
		**/
		//print a bunch of stuff
		System.out.println(this.inputFastq + "\n\n");
		System.out.println("\nGenomes\t" + "ReadCounts\t" + "% Total");

		System.out.println("Total input reads:\t" + (float)this.getTotalReads());
		
		System.out.println("\nMouse:\t" + this.getMouseReads() + "\t" +
				(float)this.getMouseReads()/(this.getTotalReads()));

		System.out.println("Human:\t" + this.getHumanReads() + "\t" + 
				(float)this.getHumanReads()/(this.getTotalReads()));

		System.out.println("E. coli:\t" + this.getEcoliReads() + "\t" +
				(float)this.getEcoliReads()/(this.getTotalReads()));

		System.out.println("Fly:\t" + this.getFlyReads() + "\t" +
				(float)this.getFlyReads()/(this.getTotalReads()));

		System.out.println("Worm:\t" + this.getWormReads() + "\t" +
				(float)this.getWormReads()/(this.getTotalReads()));

		System.out.println("Salmon:\t" + this.getSalmonReads() + "\t" +
				(float)this.getSalmonReads()/(this.getTotalReads()));

		System.out.println("Zebrafish:\t" + this.getZfishReads() + "\t" +
				(float)this.getZfishReads()/(this.getTotalReads()));

		System.out.println("Total potentially contaminating reads:\t" + this.getTotalZooReads() + "\t" +
				(float)this.getTotalZooReads()/(this.getTotalReads()));

		//close reader/writer
		br.close();
		bw.close();
	}
	
	/**
	 * This method processes each argument and assigns new variables
	 * @param args
	 */
	public void processArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		String programArgs = Misc.stringArrayToString(args, ",");
		boolean verbose = false;
		if (verbose) System.out.println("\nArguments: " + programArgs + "\n");
		for (int i = 0; i < args.length; i++) {
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()) {
				char test = args[i].charAt(1);
				try {
					switch (test) {
					case 'f': inputFastq = new String(args[++i]); break;
					case 'z': zooDB = new String(args[++i]); break;
					case 'b': bwaPath = new String(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem--unknown option used!" + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
	}

	public static void printDocs() {
		System.out.println("\n" +
				"*********************************************************************\n" +
				"**                MeasureContaminatingSeqs: May 2014               **\n" +
				"*********************************************************************\n");
	}

	public int getTotalZooReads() {
		return totalZooReads;
	}

	public int getTotalZooReadsAll() {
		return totalReadsSplit;
	}

	public int getMouseReads() {
		return mouseReads;
	}

	public int getSalmonReads() {
		return salmonReads;
	}

	public int getZfishReads() {
		return zfishReads;
	}

	public int getWormReads() {
		return wormReads;
	}

	public int getPercentMappedReadsAlignedToX() {
		return percentMappedReadsAlignedToX;
	}

	public int getPercentZooReads() {
		return percentZooReads;
	}

	public int getHumanReads() {
		return humanReads;
	}

	public int getFlyReads() {
		return flyReads;
	}

	public int getEcoliReads() {
		return ecoliReads;
	}

	public void setTotalZooReads(int totalZooReads) {
		this.totalZooReads = totalZooReads;
	}

	public void setTotalZooReadsAll(int totalZooReadsAll) {
		this.totalReadsSplit = totalZooReadsAll;
	}

	public void setMouseReads(int mouseReads) {
		this.mouseReads = mouseReads;
	}

	public void setSalmonReads(int salmonReads) {
		this.salmonReads = salmonReads;
	}

	public void setZfishReads(int zfishReads) {
		this.zfishReads = zfishReads;
	}

	public void setWormReads(int wormReads) {
		this.wormReads = wormReads;
	}

	public void setPercentMappedReadsAlignedToX(int percentMappedReadsAlignedToX) {
		this.percentMappedReadsAlignedToX = percentMappedReadsAlignedToX;
	}

	public void setPercentZooReads(int percentZooReads) {
		this.percentZooReads = percentZooReads;
	}

	public void setHumanReads(int humanReads) {
		this.humanReads = humanReads;
	}

	public void setFlyReads(int flyReads) {
		this.flyReads = flyReads;
	}

	public void setEcoliReads(int ecoliReads) {
		this.ecoliReads = ecoliReads;
	}

	public static int getTotalReads() {
		return totalReads;
	}

	public static void setTotalReads(long z) {
		MeasureContaminatingSeqs.totalReads = (int) z;
	}

	public String getSplitFastq() {
		return splitFastq;
	}

	public void setSplitFastq(String splitFastq) {
		this.splitFastq = splitFastq;
	}
}