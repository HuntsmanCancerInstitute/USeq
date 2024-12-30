package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;

import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import edu.utah.seq.analysis.ScanSeqs;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamAlignmentFlags;
import edu.utah.seq.useq.data.IntersectingRegions;
import edu.utah.seq.useq.data.Region;

/**
 * @author david.nix@hci.utah.edu
 * 
 * FT-SA159358_R1.fastq.gz
 * SL350323_1.fastq.gz
 * 24426X1_xxxxxx_L001_R1_001.fastq.gz
 **/
public class MergePairedFastqs{
	
	//user defined fields
	private File[] allGzippedFiles;
	private String nameSeparator = "-";
	private boolean deleteFastqs = false;
	
	//internal
	private Pattern internalR1 = Pattern.compile(".+(_R1_).+");
	private Pattern internalR2 = Pattern.compile(".+(_R2_).+");
	
	private String[] end1s = {"_1.fastq.gz", ".1.fastq.gz","_R1.fastq.gz", ".R1.fastq.gz", "_1.fq.gz", ".1.fq.gz"};
	private String[] end2s = {"_2.fastq.gz", ".2.fastq.gz","_R2.fastq.gz", ".R2.fastq.gz", "_2.fq.gz", ".2.fq.gz"};

	private ArrayList<File> toMergeFile1s = new ArrayList<File>();
	private ArrayList<File> toMergeFile2s = new ArrayList<File>();
	private ArrayList<String> file1Names = new ArrayList<String>();
	private ArrayList<String> file2Names = new ArrayList<String>();
	
	private ArrayList<String> filesNotMatched = new ArrayList<String>();
	
	private String mergedName = null;
	
	//constructor
	public MergePairedFastqs(String[] args){
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		
		processFiles();
		
		checkFiles();
		
		catFiles();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("Done! "+Math.round(diffTime)+" min\n");
	}

	private void catFiles() {
		ArrayList<String> al1 = new ArrayList<String>();
		for (File f: toMergeFile1s) al1.add(f.toString());
		ArrayList<String> al2 = new ArrayList<String>();
		for (File f: toMergeFile2s) al2.add(f.toString());
		String parent = toMergeFile1s.get(0).getParent();
		
		IO.pl("\nCat merging:\n\t"+ 
				Misc.stringArrayListToString(al1, ", ")+"\n\t"+
				Misc.stringArrayListToString(al2, ", ")+"\n");
		
		String files1Space = Misc.stringArrayListToString(al1, " ");
		String files2Space = Misc.stringArrayListToString(al2, " ");
		
		String cmd = "cat "+files1Space+" > "+parent+"/"+mergedName+"_1.fastq.gz && "+
				     "cat "+files2Space+" > "+parent+"/"+mergedName+"_2.fastq.gz && ";
		if (deleteFastqs) cmd = cmd +"rm -f "+files1Space+" "+files2Space+" && ";
		cmd = cmd + "echo COMPLETE\n";
		
		File tmp = new File(System.getProperty("java.io.tmpdir"));
		String[] results = IO.executeShellScript(cmd, tmp);
		String resultsMerged = Misc.stringArrayToString(results, "\n");
		
		if (resultsMerged.contains("COMPLETE") == false) {
			Misc.printErrAndExit("\nFailed to cat files! Aboring!\n"+resultsMerged+"\n"+cmd);
		}
				
		
		
	}

	private void checkFiles() {
		//any files not matched
		if (filesNotMatched.size()!=0) Misc.printErrAndExit("\nSome files not matched to read 1 or 2!\n"+ Misc.stringArrayListToString(filesNotMatched, ", ")+"\nAborting\n");
		
		//diff number
		if (toMergeFile1s.size() != toMergeFile2s.size()) {
			Misc.printErrAndExit("\nDifferent number of fastq read 1 and 2!\n\t"+ 
					Misc.stringArrayListToString(toMergeFile1s, ", ")+"\n\t"+
					Misc.stringArrayListToString(toMergeFile2s, ", ")+"\n");
		}
		
		//check names
		mergedName = Misc.stringArrayListToString(file1Names, nameSeparator);
		String name2 = Misc.stringArrayListToString(file2Names, nameSeparator);
		if (mergedName.equals(name2)==false) Misc.printErrAndExit("\nMerged file names differ!\n\t"+mergedName+"\n\t"+name2+"\n\tAborting.");
	}

	private void processFiles() {
		for (File f: allGzippedFiles) {
			String fileName = f.getName();
			boolean notFound = true;
			// end in 1?
			for (String e1: end1s) {
				if (fileName.endsWith(e1)) {
					toMergeFile1s.add(f);
					file1Names.add(fileName.replace(e1, ""));
					notFound = false;
					break;
				}
			}
			// end in 2?
			if (notFound) {
				for (String e2: end2s) {
					if (fileName.endsWith(e2)) {
						toMergeFile2s.add(f);
						file2Names.add(fileName.replace(e2, ""));
						notFound = false;
						break;
					}
				}
			}
			if (notFound) {
				Matcher mat = internalR1.matcher(fileName);
				if (mat.matches()) {
					toMergeFile1s.add(f);
					int start = mat.start(1);
					file1Names.add(fileName.substring(0, start));
					notFound = false;
				}
			}
			if (notFound) {
				Matcher mat = internalR2.matcher(fileName);
				if (mat.matches()) {
					toMergeFile2s.add(f);
					int start = mat.start(1);
					file2Names.add(fileName.substring(0, start));
					notFound = false;
				}
			}
			if (notFound) filesNotMatched.add(fileName);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergePairedFastqs(args);
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[++i]); break;
					case 'd': deleteFastqs = true; break;
					case 's': nameSeparator = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					 
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		allGzippedFiles = IO.extractFiles(forExtraction,"q.gz");
		if (allGzippedFiles.length == 0) Misc.printErrAndExit("\nFailed to find any files ending in q.gz in "+forExtraction+", aborting.\n");
		if (allGzippedFiles.length == 2) Misc.printExit("\nJust two files ending in q.gz in "+forExtraction+", nothing to do!\n");
		
		
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Merge Paired Fastqs: Sept 2024                       **\n" +
				"**************************************************************************************\n" +
				"Uses cat to combine multiple sets of paired fastq.gz files.\n"+

				"\nOption:\n"+
				"-f Directory containing the fastqs to cat\n" +
				"-s Merged name separator, defaults to '-'\n"+
				"-d Delete the starting fastqs after successfully merging. Will delete softlinks but\n"+
				"      not the sources.\n"+
				
				"\nExample: java -jar pathTo/USeq_xxx/Apps/MergePairedFastqs -f .\n\n" +


				"**************************************************************************************\n");

	}

}
