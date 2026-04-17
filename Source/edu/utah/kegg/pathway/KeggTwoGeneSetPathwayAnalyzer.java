package edu.utah.kegg.pathway;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.*;

import edu.utah.hci.bioinfo.smm.Util;
import util.gen.*;

public class KeggTwoGeneSetPathwayAnalyzer {

	//user fields
	private File keggNetworkDirectory;
	private File keggIdsFile;
	private File resultsDirectory = null;
	private File interrogatedGeneList;
	private File selectGeneListOne;
	private File selectGeneListTwo;
	private File fullPathToR;
	private int minimumNumberGenes = 4;
	private HashSet<String> networkTypesToExclude = null;
	private String typesToExclude = null;
	private File tempDirectory = null;
	private double maximumFdr = 0.15;
	private boolean replaceNetworksWithPathways = false;
	
	private KeggGenePathwayAnalyzer genesA = null;
	private KeggGenePathwayAnalyzer genesB = null;
	private HashMap<String, ArrayList<String>> gs2ki = null;
	
	//constructor
	public KeggTwoGeneSetPathwayAnalyzer(String[] args){
		try {
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		checkFiles();
		
		IO.pl("\nLaunching pathway analysis...");
		runAnalysis();
		gs2ki = genesA.getGs2ki();
		
		IO.pl("\nSaving combine pathway results...");
		buildAndSaveGenePathways();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" seconds\n");
		
		} catch (Exception e) {
			e.printStackTrace();
			IO.el("ERROR running the KeggGeneAndVariantPathwayAnalyzer!");
			System.exit(1);
		}
	}
	
	private void buildAndSaveGenePathways() throws IOException {
		CombinePathwayRoot cpr = new CombinePathwayRoot(genesA.getAnalyzedNetworks(), genesB.getAnalyzedNetworks(),"AB");	
		cpr.makeCombinePathways(maximumFdr);
		cpr.saveTwoGeneSetPathways(maximumFdr, gs2ki, resultsDirectory, minimumNumberGenes);
	}

	private void runAnalysis() throws IOException {
		//create the objects
		boolean verbose = false;
		genesA = new KeggGenePathwayAnalyzer(keggIdsFile, keggNetworkDirectory, fullPathToR, resultsDirectory, tempDirectory, interrogatedGeneList, 
				selectGeneListOne, maximumFdr, minimumNumberGenes, typesToExclude, networkTypesToExclude, verbose, replaceNetworksWithPathways, "A");
		IO.pl("\tGene Set A: "+selectGeneListOne.getName());
		genesB = new KeggGenePathwayAnalyzer(keggIdsFile, keggNetworkDirectory, fullPathToR, resultsDirectory, tempDirectory, interrogatedGeneList, 
				selectGeneListTwo, maximumFdr, minimumNumberGenes, typesToExclude, networkTypesToExclude, verbose, replaceNetworksWithPathways,"B");
		IO.pl("\tGene Set B: "+selectGeneListTwo.getName());
		
		//run the analysis in separate threads
		ExecutorService executor = Executors.newFixedThreadPool(2);
		executor.execute(genesA);
		executor.execute(genesB);
		executor.shutdown();

		//spins here until the executer is terminated, e.g. all threads complete
		while (!executor.isTerminated()) {}
		
		//write out logs
		File geneALog = new File (resultsDirectory, "geneSetAAnalysisLog.txt");
		IO.writeArrayList(genesA.getLog(), geneALog);
		File geneBLog = new File (resultsDirectory, "geneSetBAnalysisLog.txt");
		IO.writeArrayList(genesB.getLog(), geneBLog);
		
		//check both 
		if (genesA.isFailed() || genesB.isFailed()) throw new IOException();
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KeggTwoGeneSetPathwayAnalyzer(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'k': keggIdsFile = new File(args[++i]); break;
					case 'n': keggNetworkDirectory = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'm': minimumNumberGenes = Integer.parseInt(args[++i]); break;
					case 't': typesToExclude = args[++i]; break;
					case 'x': maximumFdr = Double.parseDouble(args[++i]); break;
					case 'i': interrogatedGeneList = new File(args[++i]); break;
					case 'a': selectGeneListOne = new File(args[++i]); break;
					case 'b': selectGeneListTwo = new File(args[++i]); break;
					case 'e': fullPathToR = new File(args[++i]); break;
					case 'p': replaceNetworksWithPathways = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		tempDirectory = new File(System.getProperty("java.io.tmpdir"));
		if (resultsDirectory != null) resultsDirectory.mkdirs();

		if (typesToExclude != null) {
			String[] split = Misc.COMMA.split(typesToExclude);
			networkTypesToExclude = new HashSet<String>();
			for (String s: split) networkTypesToExclude.add(s.toLowerCase());
		}
	}
	
	private void checkFiles() throws IOException {
		IO.pl("Checking required files and directories...");
		File[] toCheck = {keggIdsFile, keggNetworkDirectory, interrogatedGeneList, selectGeneListOne, selectGeneListTwo, fullPathToR, resultsDirectory};
		String[] names = {"-k keggIdsLookupFile", "-n keggNetworkDir", "-i interrogatedGeneList", "-a selectGeneListOne", "-b selectGeneListTwo", "-e fullPathToR", "-r resultsDir"};
		boolean notFound = false;
		for (int i=0; i< toCheck.length; i++) {
			if (toCheck[i] == null || toCheck[i].exists()== false) {
				IO.el("\tMissing "+names[i]+" : "+toCheck[i]);
				notFound = true;
			}
		}
		if (notFound) throw new IOException("ERROR: Correct issues with required files and restart.\n");
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                  Kegg Two Gene Set Pathway Analyzer : Mar 2026                   **\n" +
				"**************************************************************************************\n" +
				"Runs two KeggGenePathwayAnalyzer analysis and combines the results at the KEGG Pathway\n"+
				"level coloring each gene appropriately for interacive exploration in the KEGG Viewer\n"+
				"(https://www.kegg.jp). Calculates a combine Pathway p-value (Fisher's method) and FDR\n"+
				"(Benjamini-Hochberg) for Pathways with multiple significant Networks. Best to use this\n"+
				"tool if you wish to compare variant data (copy or snv/indel) with differential gene\n"+
				"expression for a single condition.\n"+
				
				"\nRelated Applications:\n\n"+
				"KeggGeneSymbolIdExtractor - Uses the Kegg API to look up and parse Kegg Gene Ids that\n"+
				"   match each of the provided HUGO Gene Symbols. Use this tool to generate the\n"+
				"   required '-k KEGG gene Id HUGO gene symbol lookup file'.\n"+
				"KeggResourceExtractor - Pulls the latest list of KEGG Networks, their associated\n"+
				"   genes, and pathways. Use this to generate the required '-n KEGG network dir'.\n"+
				"KeggGenePathwayAnalyzer - Differential gene expression KEGG Network and Pathway\n"+
				"   analysis.\n"+
				"KeggVariantPathwayAnalyzer - Differential gene mutation cohort KEGG Network and\n"+
				"   Pathway analysis.\n"+
				"MergeKeggNetworkResults and MergeKeggPathwayResults - Merges network or pathway xls\n"+
				"   results from multiple USeq KEGG analysis.\n"+
				"AnnotatedVcfParser - App to select high impact, gain/loss of function gene mutations.\n"+
				"DESeq2, edgeR - R packages for selecting differentially expressed gene sets.\n"+
				
				"\nKEGG Gene Color Key:\n\n"+
				"   Negative Genes Set A - Light Blue - #87CEEB\n"+
				"   Positive Genes Set A - Light Red, Pink - #FFB6C1\n"+
				"   Negative Genes Set B - Light Violet - #D6B4FC\n"+
				"   Positive Genes Set B - Light Orange - #FFD580\n"+
		        "   Kegg colors genes in red text that are disease pathway associated\n"+
				
				"\nRequired Parameters:\n\n"+
				"-i File containing all of the interrogated genes in your gene expression study, one\n"+
				"     gene per line. Typically > 15K genes\n"+
				"-a File containing the selected genes of interest A with their 'log2Rto', one per\n"+
				"     line, tab delimited, from your study, e.g. differentially expressed, copy\n"+
				"     mutated (loss or gain of function), etc. typically <2K.\n"+
				"     altered, 'Log2Rtos' should be + or -, not 0. Only the sign is used here, not\n"+
				"     the value.\n"+
				"-b File containing the selected genes of interest B, ditto.\n"+
				"-e File path to the R executable\n"+
				"-r Directory to save the spreadsheet results\n"+
				"-k KEGG gene Id HUGO gene symbol lookup file, see above\n"+
				"-n KEGG network directory, see above\n"+
				
				"\nOptional Parameters:\n\n"+
				"-t Network TYPEs to exclude from testing, comma delimited no spaces.\n"+
				"-m Minimum number of interrogated genes in a network for analysis,\n"+
				"     defaults to 4\n"+
				"-x Maximum FDR for including networks into the combine Kegg Pathway\n"+
				"     spreadsheet, defaults to 0.15\n"+
				"-p Replace networks with composite pathways.\n"+
	
				"\nExample:\n\n"+
				"java -Xmx1G -jar pathTo/USeq/Apps/KeggTwoGeneSetPathwayAnalyzer\n"+
				"   -i allTestedGenes.txt -a diffExpGenes.txt -b copyAlteredGenes.txt \n"+
				"   -k geneSymbol2KeggGeneInfo.txt.gz -n Networks/ -o -r CombineResults -t\n"+
				"   'Pathogen,Ev factor' -e /usr/local/bin/R\n"+

		"\n**************************************************************************************\n");
	}
}
