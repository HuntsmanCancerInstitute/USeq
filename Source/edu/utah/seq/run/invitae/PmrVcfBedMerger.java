package edu.utah.seq.run.invitae;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.hci.misc.Gzipper;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

/**Merges vcf and bed files with the same staring PMR_XXXX ID, from Invitae etc. Renames them PMRID_merged.vcf.gz or PMRID_merged.bed.gz */
public class PmrVcfBedMerger {

	private File[] vcfFiles;
	private File vcfOutputDir;
	private File[] bedFiles;
	private File bedOutputDir;

	public PmrVcfBedMerger (String[] args) {
		try {

			long startTime = System.currentTimeMillis();

			processArgs(args);

			if (vcfFiles != null) {
				IO.pl(vcfFiles.length+"\tVcf files to merge...");
				mergeVcfs(splitFilesByPmr(vcfFiles));
			}

			if (bedFiles !=null) {
				IO.pl(bedFiles.length+"\tBed files to merge...");
				mergeBedFiles(splitFilesByPmr(bedFiles));
			}

			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Min\n");
		} catch (Exception e) {
			IO.el("Problem merging the vcf and bed files...");
			e.printStackTrace();
			System.exit(1);
		}
	}

	public static TreeMap<String, ArrayList<File>> splitFilesByPmr(File[] files) {
		TreeMap<String, ArrayList<File>> pmrFiles = new TreeMap<String, ArrayList<File>>();
		for (File f: files) {
			String name = f.getName();
			String pmr = name.substring(0, name.indexOf("_"));
			ArrayList<File> al = pmrFiles.get(pmr);
			if (al == null) {
				al = new ArrayList<File>();
				pmrFiles.put(pmr, al);
			}
			al.add(f);
		}
		return pmrFiles;
	}

	private void mergeBedFiles(TreeMap<String, ArrayList<File>> pmrFiles) throws IOException {
		//for each pmr
		for (String pmr : pmrFiles.keySet()) {
			File mergedBed = new File(bedOutputDir, pmr+"_merged.bed.gz");
			ArrayList<File> beds = pmrFiles.get(pmr);
			if (beds.size()==1) {
				IO.copyViaFileChannel(beds.get(0), mergedBed);
			}
			else {
				mergeSortBeds(beds, mergedBed);
			}
		}
	}

	private void mergeVcfs(TreeMap<String, ArrayList<File>> pmrFiles) throws IOException {
		//for each pmr
		for (String pmr : pmrFiles.keySet()) {
			File mergedVcf = new File(vcfOutputDir, pmr+"_merged.vcf.gz");
			ArrayList<File> vcfs = pmrFiles.get(pmr);
			if (vcfs.size()==1) {
				IO.copyViaFileChannel(vcfs.get(0), mergedVcf);
			}
			else {
				mergeSortVcfs(vcfs, mergedVcf);
			}
		}
	}

	private void mergeSortBeds(ArrayList<File> beds, File mergedBed) throws IOException {
		//sort the files 
		File[] bedFiles = new File[beds.size()];
		beds.toArray(bedFiles);
		Arrays.sort(bedFiles);

		//merge the headers, if any
		LinkedHashSet<String> header = new LinkedHashSet<String>();
		for (File f: bedFiles) {
			String h = IO.fetchHeader(f);
			String[] lines = Misc.RETURN.split(h);
			for (String l: lines) {
				if (header.contains(l) == false) header.add(l);
			}
		}

		ArrayList<String> fileNames = new ArrayList<String>();
		TreeSet<Bed> bedRecords = new TreeSet<Bed>();
		//for each bed file
		//merge the records, the TreeSet should sort these
		for (int i=0; i< bedFiles.length; i++){
			for (Bed b: Bed.parseFile(bedFiles[i], 0, 0)) bedRecords.add(b);
			fileNames.add(bedFiles[i].getName());
		}

		//print header and records
		Gzipper out = new Gzipper(mergedBed);
		out.println("# mergedBedFiles="+Misc.stringArrayListToString(fileNames, ","));
		for (String l : header) out.println(l);
		for (Bed b : bedRecords) out.println(b.toString());
		out.close();
		
	}

	private void mergeSortVcfs(ArrayList<File> vcfs, File mergedVcf) throws IOException {
		//sort the files 
		File[] vcfFiles = new File[vcfs.size()];
		vcfs.toArray(vcfFiles);
		Arrays.sort(vcfFiles);

		VCFParser[] vcfParsers = new VCFParser[vcfFiles.length];
		ArrayList<String> fileNames = new ArrayList<String>();
		//walk them in reverse so the most recent result is processed
		for (int i=0; i< vcfParsers.length; i++){
			vcfParsers[i] = new VCFParser(vcfFiles[i], true, false, false);
			fileNames.add(vcfFiles[i].getName());
		}

		String[] mergedHeader = VCFParser.mergeHeaders(vcfParsers, false, false);
		if (mergedHeader == null) Misc.printErrAndExit("\nError: hmm something is wrong when merging headers, are the #CHROM lines different?\n");

		//merge the records taking first occurance
		LinkedHashMap<String, VCFRecord> mergedRecords = new LinkedHashMap<String, VCFRecord>();
		for (int i=0; i< vcfParsers.length; i++){
			for (VCFRecord r : vcfParsers[i].getVcfRecords()) {
				String key = r.getChrPosRefAlt(false);
				if (mergedRecords.containsKey(key) == false) mergedRecords.put(key, r);
			}
		}

		VCFRecord[] mergedVcfs = new VCFRecord[mergedRecords.size()];
		mergedRecords.values().toArray(mergedVcfs);
		Arrays.sort(mergedVcfs);				

		//print header and records
		Gzipper out = new Gzipper(mergedVcf);
		for (String l : mergedHeader) {
			//add header line for merged files
			if (l.startsWith("#CHROM")) out.println("##MergedVcfFiles="+Misc.stringArrayListToString(fileNames, ","));
			out.println(l);
		}
		for (int i=0; i< mergedVcfs.length; i++) {
			String[] fields = Misc.TAB.split(mergedVcfs[i].getOriginalRecord());
			//replace IDs so these are unique to the file
			fields[2] = "Invitae_"+i;
			out.println(Misc.stringArrayToString(fields, "\t"));
		}
		out.close();
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PmrVcfBedMerger(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcfFiles = IO.extractFiles(new File(args[++i]),".vcf.gz"); break;
					case 'b': bedFiles = IO.extractFiles(new File(args[++i]),".bed.gz"); break;
					case 'c': vcfOutputDir = new File(args[++i]); break;
					case 'e': bedOutputDir = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//pull vcf files
		if (vcfFiles != null) vcfOutputDir.mkdirs();
		//pull bed files
		if (bedFiles != null) bedOutputDir.mkdirs();
		//missing both?
		if (vcfFiles == null && bedFiles == null) Misc.printErrAndExit("\nError: please provide a directory of vcfs and or directory of bed files to merge.\n");

	}	

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                             Pmr Vcf Bed Merger: Oct 2023                         **\n" +
				"**************************************************************************************\n" +
				"Parses directories of Vcf and Bed files named such that the first underscore separated\n"+
				"field represents a subject's ID.  These will be used to group the vcf and bed files\n"+
				"and merge them by removing duplicate lines.\n"+

				"\nOptions:\n"+
				"-v Directory containing xxx.vcf.gz files\n" +
				"-b Directory containing xxx.bed.gz files\n" +
				"-c Directory to save the merged vcf files\n"+
				"-e Directory to save the merged bed files\n"+

				"\nExample: java -jar pathToUSeq/Apps/PmrVcfBedMerger -v VcfFiles/ -c MergedVcfs/\n"+
				"        -b BedFiles/ -e MergedBeds/ \n\n"+


				"**************************************************************************************\n");

	}

}
