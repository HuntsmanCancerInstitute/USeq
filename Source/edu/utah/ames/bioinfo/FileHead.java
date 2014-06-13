package edu.utah.ames.bioinfo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import util.gen.Misc;

public class FileHead {

	private static String inFile;
	private static String outFile;
	private static int numLines;

	//constructor
	public FileHead(String[] args) {
		processArgs(args);
	}

	public static void main(String[] args) {
		//check for args
		if (args.length == 0) {
			System.exit(0);
		}
		FileHead fh = new FileHead(args);
		try {
			fh.getFileLines();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void getFileLines() throws Exception {
		BufferedReader br = null;
		BufferedWriter bw = null;
		outFile = Misc.removeExtension(inFile) + "." + numLines + ".txt.gz";

		
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream
					(new FileInputStream(inFile))));
			bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream
					(new FileOutputStream(outFile))));
			//bw = new BufferedWriter(new FileWriter(new File(outFile)));
			
			String line;
			int x = 0; //num lines read
			//get specified lines
			while ((line = br.readLine()) != null && x < numLines) {
				System.out.println("blar");
				x++;
				bw.write(line);
				bw.newLine();	
				System.out.println("stuff");
			}
		System.out.println("Successfully copied " + numLines + " lines to " + outFile.toString());
		bw.close();
		br.close();
	}

	public static void processArgs(String[] args) {
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
					case 'f': inFile = new String(args[++i]); break;
					case 'n': numLines = new Integer(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem--unknown option used!" + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
		//print error and exit if no input file specified
		if (inFile == null) Misc.printErrAndExit("\nPlease provide full path to input file to parse.\n");
	}
}
