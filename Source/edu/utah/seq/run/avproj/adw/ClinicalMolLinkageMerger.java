package edu.utah.seq.run.avproj.adw;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

/**Tool for merging clin mol linkage files. M2Gen is deleting tumor and normal samples that fail qc from the most recent reports.*/
public class ClinicalMolLinkageMerger {
	
	//fields
	private String header = null;
	private String[] headerKeys = null;
	private HashSet<String> uniqueLines = new HashSet<String>();
	private HashMap<String, Integer> headerKeyIndex = new HashMap<String, Integer>();
	private HashMap<String, ArrayList<String[]>> avatarIdDataLines = new HashMap<String, ArrayList<String[]>>();
	
	public static final Pattern commaQuote = Pattern.compile("\",\"");
	public static final Pattern quote = Pattern.compile("\"");
	public static final Pattern dash = Pattern.compile("-");
	
	private Integer avatarIdIndex = null;
	private Integer tumNormIndex = null;
	private Integer specimenIdIndex = null;
	private Integer wesIndex = null;
	private Integer rnaSeqIndex = null;
	
	private PrintWriter out = null;

	public ClinicalMolLinkageMerger( File[] csvFiles, File output) throws IOException {
		
		//parse all of the individual files
		for (int i= csvFiles.length-1; i>=0; i--) {
			IO.pl("Parsing "+csvFiles[i].getName());
			parseHeader(csvFiles[i]);
			parseIt(csvFiles[i]);
		}
		out = new PrintWriter( new FileWriter(output));
		save(headerKeys);
		collapseRecords(output);
		out.close();
	}

	private void collapseRecords(File output) throws IOException {
		
		
		//for each avatarId_tumNorm_specimenId 
		//         A007039_Germline_16-0061475.1b
		for (String key: avatarIdDataLines.keySet()) {
			
			IO.pl("\n"+ key);
			ArrayList<String[]> al = avatarIdDataLines.get(key);
			
			//just one? no prob!
			if (al.size()==1) {
				IO.pl("Ju1\t"+ Misc.stringArrayToString(al.get(0), "\t"));
				save(al.get(0));
			}
			
			//scan em to collapse
			else {
				IO.pl("ToCollapse:");
				for (String[] l: al) IO.pl("\t"+ Misc.stringArrayToString(l, "\t"));
				//Germline or Tumor
				String[] best = null;
				if (key.contains("Germline")) best = collapseGermline(al);
				else best = collapseTumor(al);
				IO.pl("Col\t"+ Misc.stringArrayToString(best, "\t"));
				save(best);
			}
			
		}
	}
	
	private void save(String[] f) {
		String line = Misc.stringArrayToString(f, "\",\"");
		out.println("\""+ line + "\"");
	}

	private String[] collapseTumor(ArrayList<String[]> al) throws IOException {
		String[] mostRecent = al.get(0);
		HashSet<String> wesIds = new HashSet<String>();
		HashSet<String> rnaIds = new HashSet<String>();
		
		// add in wes and rna SL ids
		for (String[] line: al) {
			if (line[wesIndex].length()!=0) wesIds.add(line[wesIndex]);
			if (line[rnaSeqIndex].length()!=0) rnaIds.add(line[rnaSeqIndex]);
		}
		// more than one?
		if (wesIds.size()>1) throw new IOException ("\nERROR, found multiple WES IDs?! Manually fix "+wesIds);
		if (rnaIds.size()>1) throw new IOException ("\nERROR, found multiple RNASeq IDs?! Manually fix "+rnaIds);
		
		//add in the ids to the most recent
		if (wesIds.size()==1) mostRecent[wesIndex] = wesIds.iterator().next();
		if (rnaIds.size()==1) mostRecent[rnaSeqIndex] = rnaIds.iterator().next();
		
		return mostRecent;
	}

	private String[] collapseGermline(ArrayList<String[]> al) {
		// take first
		return al.get(0);
		
	}

	private void parseHeader(File cvsFromM2Gen) {
		String line = null;
		
		//clear any prior info
		header = null;
		headerKeyIndex.clear();
		avatarIdIndex= null;
		tumNormIndex= null;
		specimenIdIndex= null;
		wesIndex= null;
		rnaSeqIndex = null;
				
		try {
			BufferedReader in = IO.fetchBufferedReader(cvsFromM2Gen);
			//parse header line
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() == 0) continue;
				//swap out "," for tabs
				line = commaQuote.matcher(line).replaceAll("\t");
				line = quote.matcher(line).replaceAll("");
				if (line.contains("ORIENAvatarKey")) {
					header = line;
					headerKeys = Misc.TAB.split(header);
					for (int i=0; i<headerKeys.length; i++) headerKeyIndex.put(headerKeys[i], i);
					break;
				}
			}
			in.close();
			if (header == null) throw new IOException("Failed to find a line with 'ORIENAvatarKey'");
			
			//pull constant keys      ORIENAvatarKey, Tumor/Germline, ORIENSpecimenID
			avatarIdIndex = headerKeyIndex.get("ORIENAvatarKey");
			tumNormIndex = headerKeyIndex.get("Tumor/Germline");
			specimenIdIndex = headerKeyIndex.get("ORIENSpecimenID");
			wesIndex = headerKeyIndex.get("WES");
			rnaSeqIndex = headerKeyIndex.get("RNASeq");
			
			if (avatarIdIndex == null || tumNormIndex == null || specimenIdIndex == null || wesIndex == null || rnaSeqIndex == null) {
				throw new IOException("Failed to find one of the header keys, ORIENAvatarKey, Tumor/Germline, ORIENSpecimenID, WES, or RNASeq");
			}
			
		} catch (IOException e) {
			IO.el("Failed to parse the header, "+ line + " in "+cvsFromM2Gen);
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void parseIt(File cvsFromM2Gen) {
		String line = null;
		try {
			BufferedReader in = IO.fetchBufferedReader(cvsFromM2Gen);
			//parse header line
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() == 0 || line.contains("ORIENAvatarKey")) continue;
				//swap out "," for tabs
				line = commaQuote.matcher(line).replaceAll("\t");
				line = quote.matcher(line).replaceAll("");
				addToHash(line);
			}
			in.close();
		} catch (IOException e) {
			IO.el("Failed to parse "+ line + " in "+cvsFromM2Gen);
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void addToHash(String line) {
		if (uniqueLines.contains(line)) return;
		uniqueLines.add(line);
		String[] f = Misc.TAB.split(line);
		String avatarId = f[avatarIdIndex];
		String tumNorm = f[tumNormIndex];
		String specimenId = f[specimenIdIndex];
		String idKey = avatarId+"_"+tumNorm+"_"+specimenId;
		ArrayList<String[]> al = avatarIdDataLines.get(idKey);
		if (al == null) {
			al = new ArrayList<String[]>();
			avatarIdDataLines.put(idKey, al);
		}
		al.add(f);
	}
	
	public static void main (String[] args) throws IOException {
		if (args.length !=2) Misc.printExit("Provide a path to directory containing xxx.csv ClinicalMolLinkage files "
				+ "named so that they sort by age with the most recent last.\nAlso provide a path to a xxx.csv file for saving the results.");
		File[] toMerge = IO.extractFiles(new File(args[0]), ".csv");
		File results = new File(args[1]);
		new ClinicalMolLinkageMerger(toMerge, results);
	}
}
