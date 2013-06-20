package edu.utah.diagnostics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.sf.samtools.*;

import util.gen.Misc;

public class CountChromosomes {
	private File inputFile = null;
	private File path = null;
	private SAMFileHeader header = null;
	private File outputFile = null;
	private String reference = null;
	private ArrayList<String> supportedReferences = new ArrayList<String>(){{
		add("hg19");
	}};
	private HashMap<String,ArrayList<ArrayList<String>>> referenceChroms = new HashMap<String, ArrayList<ArrayList<String>>>();

	public CountChromosomes(String[] args) {
		this.parseArgs(args);
		
		SAMFileReader sr = new SAMFileReader(inputFile);
		header = sr.getFileHeader();
		
		referenceChroms.put("hg19", new ArrayList<ArrayList<String>>());
		ArrayList<String> stdChroms = new ArrayList<String>(){{
			add("chr1");
			add("chr2");
			add("chr3");
			add("chr4");
			add("chr5");
			add("chr6");
			add("chr7");
			add("chr8");
			add("chr9");
			add("chr10");
			add("chr11");
			add("chr12");
			add("chr13");
			add("chr14");
			add("chr15");
			add("chr16");
			add("chr17");
			add("chr18");
			add("chr19");
			add("chr20");
			add("chr21");
			add("chr22");
			add("chrX");
			add("chrY");
			add("chrM");
		}};
		
		ArrayList<String> phiX = new ArrayList<String>() {{
			add("chrPhiX_Illumina");
		}};
		
		ArrayList<String> adapter = new ArrayList<String>() {{
			add("chrAdapter");
		}};
		
		referenceChroms.get("hg19").add(stdChroms);
		referenceChroms.get("hg19").add(phiX);
		referenceChroms.get("hg19").add(adapter);
		
		
	}
	
	public static void main(String[] args) {
		CountChromosomes cc = new CountChromosomes(args);
		cc.run();
	}
	
	private void run() {
		String allReads = null;
		String alignReads = null;
		String hg19Reads = null;
		String phixReads = null;
		String adaptReads = null;
		
		ProcessBuilder pb = null;
		
		//Grab the total number of reads
		pb = new ProcessBuilder(path.toString(),"view","-c",inputFile.toString());
		allReads = this.runSystemCommand(pb, "Count all Reads");
		
		//Create standard command
		ArrayList<String> stdCommand = new ArrayList<String>();
		stdCommand.add(path.toString());
		stdCommand.add("view");
		stdCommand.add("-F");
		stdCommand.add("0x4");
		stdCommand.add("-c");
		stdCommand.add(inputFile.toString());
		
		//Grab the number of aligned reads 
		pb = new ProcessBuilder(stdCommand);
		alignReads = this.runSystemCommand(pb, "Count Aligned Reads");
		
		//Grab the number of standard hg19 reads
		if (!checkChromosomeExistence(referenceChroms.get(reference).get(0))) {
			ArrayList<String> hg19Command = new ArrayList<String>();
			hg19Command.addAll(stdCommand);
			hg19Command.addAll(this.referenceChroms.get(reference).get(0));
			
			pb = new ProcessBuilder(hg19Command);
			hg19Reads = this.runSystemCommand(pb, "Count hg19 reads");
		} else {
			System.out.println("Some or all of the hg19 chromsomes are missing from the BAM file, reporting zero reads");
			hg19Reads = "0";
		}
		
		
		//Grab the number of phix Reads
		if (!checkChromosomeExistence(referenceChroms.get(reference).get(1))) {
			ArrayList<String> phixCommand = new ArrayList<String>();
			phixCommand.addAll(stdCommand);
			phixCommand.addAll(this.referenceChroms.get(reference).get(1));
			
			pb = new ProcessBuilder(phixCommand);
			phixReads = this.runSystemCommand(pb, "Count phiX reads");
		} else {
			System.out.println("Some or all of the phiX chromsomes are missing from the BAM file, reporting zero reads");
			phixReads = "0";
		}
		
		if (!checkChromosomeExistence(referenceChroms.get(reference).get(1))) {
			ArrayList<String> adapCommand = new ArrayList<String>();
			adapCommand.addAll(stdCommand);
			adapCommand.addAll(this.referenceChroms.get(reference).get(2));
			
			pb = new ProcessBuilder(adapCommand);
			adaptReads = this.runSystemCommand(pb, "Count adapter reads");
		} else {
			System.out.println("Some or all of the adapter chromsomes are missing from the BAM file, reporting zero reads");
			adaptReads = "0";
		}
		
		
		BufferedWriter br = null;
		
		try {
			br = new BufferedWriter(new FileWriter(outputFile));
			br.write(allReads + "\n" + alignReads + "\n" + hg19Reads + "\n" + phixReads + "\n" + adaptReads + "\n");
			br.close();
		} catch (IOException ioex) {
			System.out.println("Could not write to file, exiting");
			System.exit(1);
		} finally {
			try {
				br.close();
			} catch (Exception ex) {};
		}
	}
	
	private boolean checkChromosomeExistence(ArrayList<String> chroms) {
		boolean missing = false;
		for (String c: chroms) {
			if (header.getSequenceDictionary().getSequence(c) == null) {
				missing = true;
				break;
			}
		}
		
		return missing;
	}
	
	private String runSystemCommand(ProcessBuilder pb, String commandName) {
		StringBuilder errorOut = new StringBuilder();
		StringBuilder standardOut = new StringBuilder();
		
		
		System.out.println("Running: " + commandName);
		StringBuilder cmd = new StringBuilder(); 
		for (String c: pb.command()) {
			cmd.append(c + " ");
		}
		System.out.println(cmd.toString().trim());
		
		
		
		try {
			Process pa = pb.start();
			
			BufferedReader bro = new BufferedReader(new InputStreamReader(pa.getInputStream()));
			BufferedReader bre = new BufferedReader(new InputStreamReader(pa.getErrorStream()));
			String temp = null;
			while ((temp = bro.readLine()) != null) {
				standardOut.append(temp + "\n");
			}
			
			while((temp = bre.readLine()) != null) {
				errorOut.append(temp + "\n");
			}
			
			int retVal = pa.waitFor();
			
			if (retVal != 0) {	
				System.out.println("Error running " + commandName + ", reporting zero: \n" + errorOut.toString());
				standardOut.append("0");
			}
			
			
			
		} catch (IOException ioex) {
			System.out.println("Could not run command, exiting: " + ioex.getMessage());
			ioex.printStackTrace();
			System.exit(1);
		} catch (InterruptedException iex) {
			System.out.println(commandName + " was interrupted, exiting: " + iex.getMessage());
			iex.printStackTrace();
			System.exit(1);
		}
		
		return standardOut.toString().trim();
	}
	
	private void parseArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': inputFile = new File(args[++i]); break;
					case 'o': outputFile = new File(args[++i]); break;
					case 'r': reference = args[++i]; break;
					case 'p': path = new File(args[++i]); break;
					case 'h': usage(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		
		
		if (outputFile == null) {
			exitMessage("Output File not specified, exiting");
		}
		
		if (reference == null) {
			exitMessage("Reference File not specified, exiting");
		}
		
		if (inputFile == null) {
			exitMessage("Input File not specified, exiting");
		}
		
		if (!inputFile.exists()) {
			exitMessage("Input file can't be found, exiting: " + inputFile.getName());
		}
		
		if (!this.supportedReferences.contains(reference)) {
			exitMessage("Reference sequence not currently supported. Contact developers for addition");
		}
		
	}
	
	private void usage() {
		System.out.println("***************************************************************");
		System.out.println("*                      CountChromosomes                       *");
		System.out.println("*                                                             *");
		System.out.println("* This script drives samtools view command.  It will create   *");
		System.out.println("* a report that lists counds to standard chroms, extra        *");
		System.out.println("* chroms, phiX and adatpter.  This data will be used in the   *");
		System.out.println("* ParseMetrics App.                                           *");
		System.out.println("*                                                             *");
		System.out.println("* -i Input file (bam format)                                  *");
		System.out.println("* -o Output file (.txt format)                                *");
		System.out.println("* -r Reference (hg19, hg18, mm10, mm9 etc.                    *");
		System.out.println("* -p path to samtools                                         *");   
		System.out.println("***************************************************************\n");
	}
	
	private void exitMessage(String message) {
		usage();
		System.out.println(message);
		System.exit(0);
	}

}
