package edu.utah.seq.cnv.wgs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class CopyAnalysisBedMerger {

	//User fields
	private File[] bedFilesToParse = null;
	private File geneFileToParse = null;
	private File saveFile = null;
	private float minAbsLg2Tum = 0;

	//Internal fields
	private HashMap<String, String> aliasesToOncoKBGene = new HashMap<String, String>();
	private TreeMap<String, Patient> patients = new TreeMap<String, Patient>();
	private String[] genes = null;

	//https://www.oncokb.org/cancer-genes download "Cancer Gene List" 1231 genes, cancerGeneList.tsv
	//https://www.oncokb.org/gene/AR/somatic for url link 

	public CopyAnalysisBedMerger(String[] args) {

		try {

			processArgs(args);

			parseBedFiles();

			parseGeneFile();

			//printJustPanelCopyAlteredGenes();
			//printGeneStats();

			printGenesSpreadsheet();

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: running the CopyAnalysisBedMerger\n");
		}
		IO.pl("Done");

	}

	private void parseGeneFile() throws IOException{
		BufferedReader in = IO.fetchBufferedReader(geneFileToParse);
		String line = null;
		String[] fields = null;
		ArrayList<String> oncoKBGenesAL = new ArrayList<String>();
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.length()==0) continue;
			fields = Misc.TAB.split(line);
			if (line.startsWith("Hugo")) {
				//check Gene Aliases
				if (fields[15].startsWith("Gene") == false) throw new IOException("\nError: failed to find 'Gene Aliases' at index 15 in the OncoKB gene table?\n");
			}
			else {
				String oncoKBGene = fields[0];
				oncoKBGenesAL.add(oncoKBGene);
				aliasesToOncoKBGene.put(oncoKBGene, oncoKBGene);
				//any aliases?
				if (fields.length > 15) {
					String aliases = fields[15].trim();
					if (aliases.length()!=0) {
						for (String a: Misc.COMMA_WHITESPACE.split(aliases)) {
							//don't overwrite if it exists
							if (aliasesToOncoKBGene.containsKey(a)==false) aliasesToOncoKBGene.put(a, oncoKBGene);
						}
					}
				}
			}
		}
		genes = new String[oncoKBGenesAL.size()];
		oncoKBGenesAL.toArray(genes);
		in.close();
	}

	private void parseBedFiles() {
		//parse bed files
		for (File bedFile: bedFilesToParse) {
			//126006-01-001.C4D15_25KB_Hg38.called.seg.pass.bed
			//   patientId   type  window
			String[] underscore = Misc.UNDERSCORE.split(bedFile.getName());
			String[] period = Misc.PERIOD.split(underscore[0]);
			String patientId = period[0];
			String type = "NA";
			if (period.length>1) type = period[1];

			Patient p = patients.get(patientId);
			if (p == null) {
				p = new Patient(patientId);
				patients.put(patientId, p);
			}
			p.addCallSet(bedFile, type);
		}
	}

	private void printGeneStats() {
		TreeMap<String, Integer> panelGeneCount = new TreeMap<String, Integer>();
		//for each Patient
		for (String patientId: patients.keySet()) {
			//for each call set
			TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
			for (MergedCallSet mcs: mergedCalls.values()) {
				String geneString = mcs.getPanelGeneCalls();
				for (String geneStrand: Misc.WHITESPACE.split(geneString)) {
					Integer i = panelGeneCount.get(geneStrand);
					if (i==null) panelGeneCount.put(geneStrand, 1);
					else panelGeneCount.put(geneStrand, 1+i);
				}
			}
		}
		for (String geneStrand: panelGeneCount.keySet()) {
			IO.pl(geneStrand+"\t"+panelGeneCount.get(geneStrand));
		}
	}

	private void printJustPanelCopyAlteredGenes() {
		//Print results
		IO.pl("Patient\tCondition\t#CopyCalls\t#AlteredGenes\tAlteredGeneNames: - deletion, + amplification\n");
		//for each Patient
		for (String patientId: patients.keySet()) {
			IO.pl(patientId);
			//for each call set
			TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
			for (String type: mergedCalls.keySet()) {
				int totalNumCalls = mergedCalls.get(type).geneStrandCalls.size();
				int numGenes = mergedCalls.get(type).getGenePanelCalls().size();
				IO.pl("\t"+type+"\t"+totalNumCalls+"\t"+ numGenes+"\t"+ mergedCalls.get(type).getPanelGeneCalls());
			}
		}
	}

	private void printGenesSpreadsheet() throws IOException {
		PrintWriter out = new PrintWriter (new FileWriter(saveFile));
		//create header
		String[] stats = new String[] {"Lg2Tum","Lg2Norm","NumOb"};
		out.print("OncoKB Gene Link\tDataset Hits");
		//for each patient
		for (String patientId: patients.keySet()) {
			//for each type
			TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
			for (String type: mergedCalls.keySet()) {
				for (String s: stats) {
					//add header
					out.print("\t");
					out.print(patientId);
					if (type.equals("NA") == false) {
						out.print(" ");
						out.print(type);
					}
					out.print(" ");
					out.print(s);
				}
			}
		}
		out.println();

		//for each gene
		for (String g: genes) {
			String gP = g+"+";
			String gM = g+"-";
			String[] hyperlinks = fetchHyperLinks(g,gP, gM);
			ArrayList<String> amp = new ArrayList<String>();
			amp.add(hyperlinks[0]);
			ArrayList<String> del = new ArrayList<String>();
			del.add(hyperlinks[1]);
			int numAmp = 0;
			int numDel = 0;
			//for each patient
			for (String patientId: patients.keySet()) {
				//for each type
				TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
				for (String type: mergedCalls.keySet()) {
					MergedCallSet mcs = mergedCalls.get(type);
					CopyStats mergedStatsP = mcs.geneStrandCalls.get(gP);
					if (mergedStatsP == null) amp.add(".\t.\t.");
					else {
						amp.add(mergedStatsP.toStringTabs());
						numAmp++;
					}
					CopyStats mergedStatsG = mcs.geneStrandCalls.get(gM);
					if (mergedStatsG == null) del.add(".\t.\t.");
					else {
						del.add(mergedStatsG.toStringTabs());
						numDel++;
					}
				}
			}

			//insert count
			amp.add(1, new Integer(numAmp).toString());
			del.add(1, new Integer(numDel).toString());

			//print
			out.println(Misc.stringArrayListToString(amp, "\t"));
			out.println(Misc.stringArrayListToString(del, "\t"));
		}
		out.close();
	}

	private String[] fetchHyperLinks(String g, String gP, String gM) {
		StringBuilder sb = new StringBuilder();
		// https://www.oncokb.org/gene/AR/somatic
		sb.append("=HYPERLINK(\"https://www.oncokb.org/gene/");
		sb.append(g);		
		sb.append("/somatic\",\"");
		String lead = sb.toString();
		String a = lead+gP+"\")";
		String b = lead+gM+"\")";
		return new String[] {a,b};
	}

	private class CopyStats {
		// always present
		int numOb;
		float lg2Tum;
		// set to zero for no norm analysis
		float lg2Norm;

		//constructor
		public CopyStats(String stats) {
			// numOb=4041;lg2Tum=0.2932;lg2Norm=0;
			String[] s = Misc.SEMI_COLON.split(stats);
			numOb = Integer.parseInt(s[0].substring(6));
			lg2Tum = Float.parseFloat(s[1].substring(7));
			lg2Norm = Float.parseFloat(s[2].substring(8));
		}

		public String toStringTabs() {
			return lg2Tum+"\t"+lg2Norm + "\t"+numOb;
		}
	}

	private void printAllCopyGenes() {
		//Print results
		//for each Patient
		for (String patientId: patients.keySet()) {
			IO.pl(patientId);
			//for each call set
			TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
			for (String type: mergedCalls.keySet()) {
				IO.pl("\t"+type+"\t"+mergedCalls.get(type).geneStrandCalls);
			}
			IO.pl();
		}

	}

	private class Patient {
		String patientId = null;
		TreeMap<String, MergedCallSet> mergedCalls = new TreeMap<String, MergedCallSet>();

		private Patient (String patientId) {
			this.patientId = patientId;
		}

		private void addCallSet(File bedFile, String type) {
			MergedCallSet csCallSet = mergedCalls.get(type);
			if (csCallSet == null) {
				csCallSet = new MergedCallSet(type);
				mergedCalls.put(type, csCallSet);
			}
			csCallSet.addCalls(bedFile);
		}
	}

	private class MergedCallSet {
		String type = null;
		TreeMap<String, CopyStats> geneStrandCalls = new TreeMap<String, CopyStats>();
		ArrayList<String> panelGenePanelCalls = null;

		private MergedCallSet(String type) {
			this.type = type;
		}

		private void addCalls(File bedFile) {
			Bed[] bedRegions = Bed.parseFile(bedFile, 0, 0);
			for (Bed b: bedRegions) {
				//numOb=64;lg2Tum=0.2725;lg2Norm=-0.0285;genes=TTTY17C,TTTY17B,TTTY17A
				String[] split = b.getName().split("genes=");
				String genes = split[1];
				CopyStats cs = new CopyStats(split[0]);
				for (String g: Misc.COMMA.split(genes)) {
					//exclude any antisense stuff or empty call sets '.'
					if (g.contains("-AS") == false && g.contains(".") == false) {
						//try to get oncoKB gene symbol
						String oncoKBGeneSym = aliasesToOncoKBGene.get(g);
						if (oncoKBGeneSym == null) oncoKBGeneSym = g;
						//strand is used to represent + amp, - del
						String key = oncoKBGeneSym+b.getStrand();
						CopyStats oldCS = geneStrandCalls.get(key);

						if (oldCS == null) geneStrandCalls.put(key, cs);
						else {
							float oldCSTum = Math.abs(oldCS.lg2Tum);
							float newCSTum = Math.abs(cs.lg2Tum);
							//only put new cs if it has a bigger abs(lg2Tum)
							if (oldCSTum < newCSTum) geneStrandCalls.put(key, cs);
						}

					}
				}
			}
		}

		private String getPanelGeneCalls() {
			return Misc.stringArrayListToString(getGenePanelCalls(), " ");
		}

		private ArrayList<String> getGenePanelCalls(){
			if (panelGenePanelCalls == null) {
				panelGenePanelCalls = new ArrayList<String>();
				for (String g: genes) {
					String gP = g+"+";
					String gM = g+"-";
					if (geneStrandCalls.get(gP) !=null) panelGenePanelCalls.add(gP);
					if (geneStrandCalls.get(gM) !=null) panelGenePanelCalls.add(gM);
				}
			}
			return panelGenePanelCalls;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CopyAnalysisBedMerger(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File bedDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedDir = new File (args[++i]); break;
					case 'g': geneFileToParse = new File (args[++i]); break;
					case 'r': saveFile = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull bed files
		if (bedDir == null || bedDir.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a bed file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(bedDir, ".bed");
		tot[1] = IO.extractFiles(bedDir,".bed.gz");
		tot[2] = IO.extractFiles(bedDir,".vcf.zip");
		bedFilesToParse = IO.collapseFileArray(tot);
		if (bedFilesToParse == null || bedFilesToParse.length ==0 || bedFilesToParse[0].canRead() == false) {
			Misc.printExit("\nError: cannot find your xxx.bed(.zip/.gz OK) file(s)!\n");
		}

		//gene file
		if (geneFileToParse == null || geneFileToParse.exists()==false) {
			Misc.printExit("\nError: cannot find your tab delimited Hugo gene symbol (.zip/.gz OK) file(s)!\n");
		}

		if (saveFile == null || saveFile.getName().endsWith(".xls")==false) Misc.printExit("\nError: please provide a file path ending in xls for saving the results.\n");
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           CopyAnalysisBedMerger: March 2026                      **\n" +
				"**************************************************************************************\n" +
				"Merges bed files from the GATK-USeq workflow called with different window sizes. See\n"+
				"https://github.com/HuntsmanCancerInstitute/Workflows/tree/master/Hg38RunnerWorkflows2\n"+
				"/CopyAnalysis Window calling stats with the largest abs(Lg2Tum) value are selected to\n"+
				"represent the alteration. Only genes and aliases in OncoKB are analyzed. See -g\n"+

				"\nOptions:\n"+
				"-b Directory containing the xxx.called.seg.pass.bed.gz files. The file name will be\n"+
				"     split on '_' then '.' to define the patient ID, and if present, the condition\n"+
				"     e.g.  126006-01-001.Cisplatin_25KB_Hg38.called.seg.pass.bed\n"+
				"               patientID.condition_xxx datasets with the same ID.condition are merged\n"+
				"-g OncoKB gene list file. See https://www.oncokb.org/cancer-genes and the 'Cancer Gene\n"+
				"     List' download link. Only 'Hugo Symbol' and 'Gene Aliases' columns are parsed.\n"+
				"-r Path to an xls spreadsheet file for saving the results.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/CopyAnalysisBedMerger -b PassingBeds/ \n" +
				"     -g cancerGeneList.tsv -r mergedResults.xls \n\n"+ 


				"************************************************************************************\n");

	}
}
