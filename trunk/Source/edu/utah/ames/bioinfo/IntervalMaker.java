package edu.utah.ames.bioinfo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Misc;

/**
 * Generates chromosome-specific intervals of user-defined length, with output as
 * .bed file (chr, start, stop). Zero-based.
 * @author darren.ames@hci.utah.edu
 *
 */
public class IntervalMaker {
	
	//fields
	private static String inputFile;
	private static int binSize = 1000000;
	
	//constructor
	public IntervalMaker(String[] args) throws Exception {
		processArgs(args);
	}
	
	public static void main(String[] args) throws Exception {
		
		//check for args
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		IntervalMaker im = new IntervalMaker(args);
		im.makeIntervals();
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	public void makeIntervals() throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		String line;
		
		File f = new File(inputFile);
		File fNew = new File((f.getAbsolutePath().replace(".bed", "_intervals.bed")));
		
		//if file doesn't exist, create it
		if (!fNew.exists()) {
			fNew.createNewFile();
		}
		
		FileWriter fw = new FileWriter(fNew.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		//loop the read block until all lines in the file are read
		while ((line = br.readLine()) != null) {

			//split contents of tab-delimited file 
			String[] dataValue = line.split("\t");

			//check for correct number of columns in input file
			if (dataValue.length < 3) {
				continue;
			}
			else {
				Region r = new Region(dataValue);
				//System.out.println(r.getStop());
				double d = Double.parseDouble(r.getStop());
				double c = d/binSize;
				//calculate the required number of intervals to cover length of each chromosome
				int numInts = roundUp(c);
				//make the regions
				makeRegions(r, numInts, bw);
			}
		}
		//close the reader and writer
		br.close();
		bw.close();
	}
	
	/**
	 * creates the regions and prints to .bed file
	 * @param r
	 * @param numInts
	 */
	public void makeRegions(Region r, int numInts, BufferedWriter bw) {

		try {

			for (int n = Integer.parseInt(r.getStart()); n < numInts; n++) {
				if (n == numInts -1) {
					bw.write(String.format("chr%s\t%s\t%s\n", r.getChromNum(), String.valueOf(n*binSize), String.valueOf(r.getStop())));
				}
				else {
					bw.write(String.format("chr%s\t%s\t%s\n", r.getChromNum(), String.valueOf(n*binSize), String.valueOf((n+1)*binSize)));
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * rounds up a double by 1 to return an int
	 * @param d
	 * @return
	 */
	private int roundUp(double d) {
		return (d > (int) d) ? (int) d + 1 : (int) d;
	}
	
	/**
	 * This method will process each argument and assign new variables.
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
					case 'f': inputFile = new String(args[++i]); break;
					case 'n': binSize = new Integer(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem--unknown option used!" + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
		//print error and exit if no input file specified
		if (inputFile == null) Misc.printErrAndExit("\nPlease provide full path to input file to parse.\n");
	}
	
	public static void printDocs() {
		System.out.println("\n" +
				"**********************************************************************************\n" +
				"**                           IntervalMaker: July 2013                            **\n" +
				"**********************************************************************************\n" + 
				"This class contains methods to feed in a .bed file of chromosome lengths (chr, start," +
				" stop) and generate an output .bed chromosomes broken into regions of defined intervals." +
				"\nParameters: \n\n" +
				"-f filename for autoalign report to process\n\n" +
				"-n interval size (default 1Mb)\n" +

				"Usage:\n\n" +
				"java -jar pathTo/IntervalMaker.jar -f pathTo/input.bed -n 1000000\n\n" +
				"Questions or comments? Contact: darren.ames@hci.utah.edu\n" +
				"**********************************************************************************\n");
	}
}
