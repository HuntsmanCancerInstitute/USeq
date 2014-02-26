package edu.utah.seq.analysis;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;


import net.sf.samtools.*;


public class PullMatchingAlignments {
	private File readNameFile = null;
	private File alignmentFile = null;
	private File outputFile = null;
	

	public PullMatchingAlignments(String[] args) {
		parseArgs(args);
		run();
		
	}
	
	public void run() {
		//Read in read names
		System.out.println("Reading in read names");
		HashSet<String> readNames = this.readInReadNames();
		
		//Find matches
		System.out.println("Scanning alignments for matches");
		findMatches(readNames);
		
		System.out.println("Finished");
	}
	
	private void findMatches(HashSet<String> readNames) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(this.outputFile));
			
			SAMFileReader sfr = null;
			if (this.alignmentFile.toString().endsWith(".sam") || this.alignmentFile.toString().endsWith(".bam")) {
				sfr = new SAMFileReader(this.alignmentFile);
			} else if (this.alignmentFile.toString().endsWith(".sam.gz")) {
				sfr = new SAMFileReader(new GZIPInputStream(new FileInputStream(this.alignmentFile)));
			}  else {
				System.out.println("File does not end with a known suffix (bam, sam, sam.gz), exiting:  " + this.alignmentFile.toString() );
				System.exit(1);
			}
			
			 
			sfr.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			
			SAMRecordIterator sri = sfr.iterator();
			SAMRecord sr = null;
			
			
			while(sri.hasNext()) {
				sr = sri.next();
				
				if (readNames.contains(sr.getReadName())) {
					String aligned = null;
					if (sr.getReadUnmappedFlag()) {
						aligned = "false";
					} else {
						aligned = "true";
					}
					bw.write(String.format("%s\t%s\n",aligned,sr.getSAMString()));
				}	
			}
			
			sfr.close();
			bw.close();
		} catch (IOException ioex) {
			System.out.println("I/O Error: " + ioex.getMessage());
			System.exit(1);
		}
	}
	
	private HashSet<String> readInReadNames() {
		
		HashSet<String> readNames = new HashSet<String>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.readNameFile));
			String temp = null;
			
			while ((temp = br.readLine()) != null) {
				readNames.add(temp);
			}
			
			br.close();
		} catch (FileNotFoundException fnfe) {
			System.out.println("Read name file not found: " + this.readNameFile.toString());
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Error reading file, exiting: " + ioex.getMessage());
			System.exit(1);
		} 
		
		return readNames;
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
					case 'r': readNameFile = new File(args[++i]); break;
					case 'a': alignmentFile = new File(args[++i]); break;
					case 'o': outputFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}

		if (readNameFile ==null){
			System.out.println("\nMissing read name file (-r).\n");
			System.exit(1);
		}

		if (alignmentFile==null){
			System.out.println("\nMissing alignment file (-a).\n");
			System.exit(1);
		}
		
		if (outputFile==null) {
			System.out.println("\nMissing output file (-o)\n");
			System.exit(1);
		}
	}
	
	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(1);
		}
		
		new PullMatchingAlignments(args);
		
		
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                  Pull Matching Alignments: Feb 2014                              **\n" +
				"**************************************************************************************\n" +
				"Read in a list of read names and then write out alignments that match.  The output is \n" + 
				"true/false depending on alignment status and then the standard SAM entry.\n\n" +
				
				"-r The full path to a file listing read names.  One read name per line\n" +
				"-a The full path to an alignment file. The app can handle bam/sam/sam.gz.\n"+
				"-o Path to the output file.  The file will be any alignment that matches\n"+
				"     a read name in (-r).\n" +
				"\n" +
				"Example: java -Xmx500M -jar pathTo/Apps/PullMatchingAlignments -r /path/to/readnames.txt\n" +
				"      -a /path/to/alignment.bam -o /path/to/matching.txt \n\n" +

		"**************************************************************************************\n");
	}

}
