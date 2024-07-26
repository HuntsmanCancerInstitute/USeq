package edu.utah.seq.data.cbio;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.gen.IO;
import util.gen.Misc;

public class CBio2Vcf {

	//fields
	private File mutationTable = null;
	private File indexedFasta = null;
	private IndexedFastaSequenceFile fasta; 
	private String workingChromosomeName = "";
	private String workingSequence = null;
	private HashMap<String, Integer> headerNameIndex = null;
	private HashMap<String, ArrayList<String[]>> chromosomeLines = new HashMap<String, ArrayList<String[]>>(); 
	private String[] toFind = new String[]{"Gene", "Protein Change", "Variant Type", "Chromosome", "Start Pos", "End Pos", "Ref", "Var", "HGVSc", "ClinVar"};
	//										 0           1                 2                3            4           5       6      7       8           9
	private int geneIndex = -1;
	private int proteinChangeIndex = -1;
	private int variantTypeIndex = -1;
	private int chromosomeIndex = -1;
	private int startIndex = -1;
	private int endIndex = -1;
	private int refIndex = -1;
	private int altIndex = -1;
	private int hgvscIndex = -1;
	private int clinvarIndex = -1;
	
	private ArrayList<CBioVar> variantsAL = new ArrayList<CBioVar>();
	private CBioVar[] variants = null;
	
	public CBio2Vcf(String[] args){
		processArgs(args);

		parseMutationTable();
		
		//make an array and sort
		variants = new CBioVar[variantsAL.size()];
		variantsAL.toArray(variants);
		Arrays.sort(variants);
		
		fixRefAlts();
		
		appendCounts();
		
		writeVcf();
		
		IO.pl("\nCOMPLETE!");
		
	}
	
	private void writeVcf() {
		IO.pl("Writing vcf...");
		String name = mutationTable.getName();
		name = Misc.removeExtension(name);
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter(new File(mutationTable.getParentFile(), name+".vcf")));
			out.println(CBioVar.getVcfHeader(indexedFasta));
			for (int i=0; i< variants.length; i++) out.println(variants[i].toVcf(i+""));
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: writing vcf file.");
		} finally {
			IO.closeNoException(out);
		}
	}

	private void appendCounts() {
		IO.pl("Appending counts...");
		int numTotal = 0;
		LinkedHashMap<String, ArrayList<CBioVar>> keyVars = new LinkedHashMap<String, ArrayList<CBioVar>>();
		for (CBioVar v: variants) {
			String key = v.fetchKey();
			ArrayList<CBioVar> al = keyVars.get(key);
			if (al == null) {
				al = new ArrayList<CBioVar>();
				keyVars.put(key, al);
			}
			al.add(v);
		}
		variants = new CBioVar[keyVars.size()];
		int index = 0;
		for (ArrayList<CBioVar> al: keyVars.values()) {
			int num = al.size();
			numTotal+= num;
			variants[index] = al.get(0);
			variants[index].setCount(num);
			index++;
		}
		IO.pl("\tMerging "+numTotal+" to "+variants.length);
		
	}

	private void fixRefAlts() {
		IO.pl("Converting ref alt and pos to vcf format...");
		//all have chrom, startMinOne, ref, alt
		for (CBioVar var: variants) {
			loadChromosomeData(var.getChrom());
			var.checkFixRefAlt(workingSequence);
		}
	}

	private void loadChromosomeData(String chrom) {
		try {
			if (workingChromosomeName.equals(chrom)) return;
			//set new name
			workingChromosomeName = chrom;

			//find and load sequence
			ReferenceSequence p = fasta.getSequence(workingChromosomeName);
			if (p == null ) throw new IOException ("\n\nFailed to find or load a fasta sequence for '"+workingChromosomeName+"', aborting.\n");
			workingSequence = new String(p.getBases());
			workingSequence = workingSequence.toUpperCase();

		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nDo your chromosome names and the fasta match? Watch for 'chr'. "+workingChromosomeName);
		}
	}

	private void parseMutationTable() {
		IO.pl("Parsing mutation table...");
		BufferedReader in = null;
		int fail = 0;
		try {
			in = IO.fetchBufferedReader(mutationTable);
			String line = null;
			while ((line = in.readLine())!=null) {
				String[] f = Misc.TAB.split(line);
				//header line? parse it
				if (headerNameIndex == null && line.contains("Chromosome")) {
					headerNameIndex = new HashMap<String, Integer>();
					for (int i=0; i< f.length; i++) headerNameIndex.put(f[i], i);
					//check it
					for (String s: toFind) if (headerNameIndex.get(s) == null) throw new IOException("Failed to find '"+s+"' in the header! Correct and restart.");
					//assign indexes
					geneIndex = headerNameIndex.get(toFind[0]);
					proteinChangeIndex = headerNameIndex.get(toFind[1]);
					variantTypeIndex = headerNameIndex.get(toFind[2]);
					chromosomeIndex = headerNameIndex.get(toFind[3]);
					startIndex = headerNameIndex.get(toFind[4]);
					endIndex = headerNameIndex.get(toFind[5]);
					refIndex = headerNameIndex.get(toFind[6]);
					altIndex = headerNameIndex.get(toFind[7]);
					hgvscIndex = headerNameIndex.get(toFind[8]);
					clinvarIndex = headerNameIndex.get(toFind[9]);
				}
				//save it? must have a Chromosome, Start Pos, Ref, Var, 
				else {
					//public CBioVar (String gene, String chrom, int startMinOne, int stop, String ref, String alt, String proteinChange, String variantType, String hgvsC, String sampleType) {
					CBioVar var = new CBioVar(f[geneIndex], f[chromosomeIndex], f[startIndex], f[endIndex], f[refIndex], f[altIndex], f[proteinChangeIndex], f[variantTypeIndex], f[hgvscIndex], f[clinvarIndex]);
					if (var.isGoodForVcf()) variantsAL.add(var);
					else fail++;
				}
			}
			IO.pl("\tPassing: "+variantsAL.size()+" Failing: "+fail);
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			IO.closeNoException(in);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CBio2Vcf(args);
	}

	/**This method will process each argument and assign new variables
	 * @throws FileNotFoundException */
	public void processArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': indexedFasta = new File(args[++i]); break;
					case 'm': mutationTable = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//mutation file
		if (mutationTable == null || mutationTable.exists() == false) Misc.printErrAndExit("\nError: cannot find your mutation table file.");

		//Create fasta fetcher
		if (indexedFasta == null || indexedFasta.exists() == false)  Misc.printErrAndExit("\nError: cannot find indexed fasta file? -> "+ indexedFasta);
		try {
			fasta = new IndexedFastaSequenceFile(indexedFasta);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: loading fasta lookup.");
		}
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                CBio 2 Vcf : July 2024                            **\n" +
				"**************************************************************************************\n" +
				"Converts a cBioPortal mutation table download to a vcf file for snvs and indels.\n"+
				"Does not parse fusions or copy alterations. Add a column called Gene and populate it.\n"+

				"\nOptions:\n"+
				"-m Downloaded mutation xxx.tsv table for parsing.\n"+
				"-f Indexed multi chromosome fasta file.\n"+

				"\n"+

				"Example: java -Xmx1G -jar pathTo/USeq/Apps/CBio2Vcf -m table.tsv -f ~/Indexes/hg19.fa\n\n"+

				"**************************************************************************************\n");

	}		


}
