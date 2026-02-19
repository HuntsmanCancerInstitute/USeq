package edu.utah.kegg.pathway;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.kegg.api.KeggApiGene;
import edu.utah.kegg.api.KeggApiNetwork;
import edu.utah.kegg.api.KeggApiPathway;
import edu.utah.kegg.api.KeggGeneSymbolIdExtractor;
import edu.utah.kegg.api.KeggResourceExtractor;
import util.gen.*;

public class KeggGenePathwayAnalyzer implements Runnable{

	//user fields
	private File interrogatedGeneList;
	private File selectGeneList;
	private File keggNetworkDirectory;
	private File keggIdsFile;
	private File fullPathToR;
	private File resultsDirectory = null;
	private double maximumFdr = 0.15;
	private int minimumNumberGenes = 4;
	private HashSet<String> networkTypesToExclude = null;
	private String typesToExclude = null;

	//internal
	private File tempDirectory = null;
	private KeggApiNetwork[] allNetworks = null;
	private HashMap<String, AnalyzedNetwork> networkIdAnalyzedNetwork = null;
	private AnalyzedNetwork[] analyzedNetworks = null;
	private HashMap<String, ArrayList<String>> gs2ki = null;
	private HashSet<String> uniqueInterrogatedGenes = null;
	private HashSet<SelectGene> uniqueSelectGenes = null;
	private boolean failed = false;
	private boolean verbose = true;
	private ArrayList<String> log = new ArrayList<String>();
	private long startTime = -1;
	
	//constructor for cmd line
	public KeggGenePathwayAnalyzer(String[] args) {
		try {

			processArgs(args);
			
			checkFiles();

			run();
			
			if (failed) throw new Exception();
			
		} catch (Exception e) {
			e.printStackTrace();
			IO.el("ERROR running the KeggGenePathwayAnalyzer!");
			failed = true;
			System.exit(1);
		}
	}
	

	//constructor for the joint analysis
	public KeggGenePathwayAnalyzer(File keggIdsFile, File keggNetworkDirectory, File fullPathToR, File resultsDirectory, File tempDirectory, File interrogatedGeneList, 
			File selectGeneList, double maximumFdr , int minimumNumberGenes, String typesToExclude, HashSet<String> networkTypesToExclude, boolean verbose) {
			this.interrogatedGeneList = interrogatedGeneList;
			this.keggIdsFile = keggIdsFile;
			this.keggNetworkDirectory = keggNetworkDirectory;
			this.tempDirectory = tempDirectory;
			this.fullPathToR = fullPathToR;
			this.resultsDirectory = resultsDirectory;
			this.maximumFdr = maximumFdr;
			this.minimumNumberGenes = minimumNumberGenes;
			this.typesToExclude = typesToExclude;
			this.networkTypesToExclude = networkTypesToExclude;
			this.selectGeneList = selectGeneList;
			this.verbose = verbose;
	}

	public void run() {
		startTime = System.currentTimeMillis();
		
		try {
			lg("\nLoading KEGG gene symbol <-> id lookup tables...");
			loadKeggIdLookupHashes();

			lg("\nLoading your interrogated gene symbols...");
			loadInterrogatedGenes();

			lg("\nLoading your select genes...");
			loadSelectedGenes();

			lg("\nLoading and filtering KEGG Medicus networks for genes...");
			loadFilterKeggApiNetworks();

			lg("\nIntersect select genes with networks...");
			intersectGeneNetworks();

			lg("\nComparing gene networks...");
			compareGeneNetworks();

			lg("\nSaving network results...");
			saveGeneNetworks();

			lg("\nSaving significant pathways...\n");
			buildAndSaveGenePathways();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			verbose = true;
			lg("\nDone - Gene Pathway Analysis! "+Math.round(diffTime)+" seconds");

		} catch (Exception e) {
			e.printStackTrace();
			IO.el("ERROR running the KeggGenePathwayAnalyzer!");
			failed = true;
		}

	}
	
	private void lg(String message) {
		if (verbose) System.out.println(message);
		log.add(message);
	}

	private void buildAndSaveGenePathways() throws IOException {
		CombinePathwayRoot cpr = new CombinePathwayRoot(analyzedNetworks, null);
		cpr.makeCombinePathways(maximumFdr);
		cpr.saveGenePathways(maximumFdr, gs2ki, resultsDirectory, minimumNumberGenes);
	}
	
	private void loadKeggIdLookupHashes() throws IOException {
		HashMap<String, ArrayList<String>>[] hashes = KeggGeneSymbolIdExtractor.loadGeneLookupHashes(keggIdsFile);
		gs2ki = hashes[0];
	}

	private void loadFilterKeggApiNetworks() throws IOException {
		allNetworks = KeggResourceExtractor.loadNetworks(keggNetworkDirectory);
		networkIdAnalyzedNetwork = new HashMap<String, AnalyzedNetwork>(allNetworks.length);
		for (int i=0; i< allNetworks.length; i++) {
			//check type?
			if (networkTypesToExclude != null) {
				String type = allNetworks[i].getNetworkType().toLowerCase();
				if (networkTypesToExclude.contains(type)) continue; 
			}

			//check minimum # genes
			KeggApiGene[]  genes = allNetworks[i].getGenes();
			if (genes!=null && genes.length >= minimumNumberGenes) {
				networkIdAnalyzedNetwork.put(allNetworks[i].getNetworkId(), new AnalyzedNetwork(allNetworks[i], uniqueSelectGenes.size()));
			}
		}
		lg("\t"+networkIdAnalyzedNetwork.size()+"\tNetworks loaded that pass minimum # genes ("+minimumNumberGenes+") and excluded types: "+typesToExclude+"\n");
		
		//for each network gene, look for them in the user's interrogated genes
		for (AnalyzedNetwork an: networkIdAnalyzedNetwork.values()) {
			TreeSet<String> found = an.getGeneNetworkGeneNamesInInterrogatedGenes();
			for (KeggApiGene kag: an.getKeggApiNetwork().getGenes()) {
				if (uniqueInterrogatedGenes.contains(kag.getName())) found.add(kag.getName());
			}
			int numNetGenes = an.getKeggApiNetwork().getGenes().length;
			int numFound = found.size();
			if (numFound < minimumNumberGenes) {
				lg("\tSkipping '"+an.getKeggApiNetwork().getNetworkIdName(" : ")+"' too few of its genes ("+numFound+") intersected your interrogated gene list.");
				an.setGeneAnalyzeIt(false);
			}
			else if (numFound != numNetGenes) {
				//lg("\t"+numFound+" of "+numNetGenes+" were found in your interrogated gene list for '"+ an.getKeggApiNetwork().getNetworkName()+"'");
			}
		}

		//merge networks since the # genes has changed, don't want to test these twice
		lg("\nMerging networks with the same gene sets...");
		HashMap<String, ArrayList<AnalyzedNetwork>> geneNameStringNetwork = new HashMap<String, ArrayList<AnalyzedNetwork>>();
		for (AnalyzedNetwork an: networkIdAnalyzedNetwork.values()) {
			if (an.isGeneAnalyzeIt() == false) continue;
			TreeSet<String> found = an.getGeneNetworkGeneNamesInInterrogatedGenes();
			String key = Misc.treeSetToString(found, ",");
			ArrayList<AnalyzedNetwork> al = geneNameStringNetwork.get(key);
			if (al == null) {
				al = new ArrayList<AnalyzedNetwork>();
				geneNameStringNetwork.put(key, al);
			}
			al.add(an);	
		}

		//make final set to intersect
		analyzedNetworks = new AnalyzedNetwork[geneNameStringNetwork.size()];
		int index = 0;
		for (String genes: geneNameStringNetwork.keySet()) {
			ArrayList<AnalyzedNetwork> al = geneNameStringNetwork.get(genes);
			//lg("\t"+ al.size()+"\tNetworks with "+genes);
			analyzedNetworks[index] = al.get(0);
			//add in networks
			ArrayList<KeggApiNetwork> n = analyzedNetworks[index].getAnalizedKeggApiNetworks(); 
			for (AnalyzedNetwork an: al) n.add(an.getKeggApiNetwork());
			index++;
		}
		lg("\t"+analyzedNetworks.length+"\tMerged networks.");
		
	}

	private void loadSelectedGenes() throws IOException {
		uniqueSelectGenes = new HashSet<SelectGene>();

		String[] selectGenes = IO.loadFileIntoStringArray(selectGeneList);
		ArrayList<String> failedToFind = new ArrayList<String>();
		for (String s: selectGenes) {
			//GeneName Log2Rto
			String[] nameRto = Misc.TAB.split(s);
			if (nameRto.length!=2) throw new IOException ("Failed to parse the gene name and it's log2Rto for the select gene "+s);
			SelectGene sg = new SelectGene(nameRto);
			if (uniqueInterrogatedGenes.contains(sg.getGeneSymbol())) uniqueSelectGenes.add(sg);
			else failedToFind.add(sg.getGeneSymbol());
		}
		if (failedToFind.size()!=0) lg("\tWARNING: failed to find "+failedToFind.size()+" in the interrogated gene list, skipping: "+Misc.stringArrayListToString(failedToFind, ","));
		lg("\t"+selectGenes.length+" -> "+uniqueSelectGenes.size()+ " unique selected genes");

	}

	private void loadInterrogatedGenes() {
		String[] allinterrogatedGenes = IO.loadFileIntoStringArray(interrogatedGeneList);
		uniqueInterrogatedGenes = new HashSet<String>();
		for (String s: allinterrogatedGenes) {
			if (gs2ki.containsKey(s)) uniqueInterrogatedGenes.add(s);
		}
		lg("\t"+allinterrogatedGenes.length+" -> "+uniqueInterrogatedGenes.size()+ " unique HUGO gene symbols that are recognized by KEGG.");
	}

	private void saveGeneNetworks() throws IOException {
		PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "geneNetworksMinGen"+minimumNumberGenes+".xls")));
		//name descLink pval, adjPval, #UniqueNetworkGenes, #FoundUniqueNetworkGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/PathGenes, IntersectingGenes
		out.println("# Network Name(s)\tNetwork Desc Link\tPval\tAdj Pval\tAll Network Genes\tFound Network Genes\tSelect Genes\tIntersect\tInt/Net\tInt Gene Symbols\tInt Genes LgRtos\tPathway Map Links...");
		for (AnalyzedNetwork an: analyzedNetworks) out.print(an.toStringGene(gs2ki));
		out.println("\n"+CombinePathwayRoot.fetchColorKey(true, false));
		out.close();
	}


	private void intersectGeneNetworks() {
		//for each AnalyzedNetwork
		ArrayList<AnalyzedNetwork> withHits = new ArrayList<AnalyzedNetwork>();
		for (AnalyzedNetwork p: analyzedNetworks) {
			//pull the genes from the network after filtered against the users interrogated gene list
			TreeSet<String> networkGeneNames = p.getGeneNetworkGeneNamesInInterrogatedGenes();
			ArrayList<SelectGene> sharedGenes = p.getGeneSharedGenes();
			//for each select gene from the user, if it intersects the network genes save it.
			Iterator<SelectGene> it = uniqueSelectGenes.iterator();
			while (it.hasNext()) {
				SelectGene sg = it.next();
				if (networkGeneNames.contains(sg.getGeneSymbol())) sharedGenes.add(sg);
			}
			if (sharedGenes.size()!=0) withHits.add(p);
		}
		//replace AnalyzedNetworks
		AnalyzedNetwork[] finalAn = new AnalyzedNetwork[withHits.size()];
		withHits.toArray(finalAn);
		lg("\t"+finalAn.length+"\tNetworks found with genes that intersect your select list.");
		analyzedNetworks = finalAn;	
	}
		
	private void compareGeneNetworks() throws IOException {
		int[][] nabt = new int[analyzedNetworks.length][];
		//for each AnalyzedNetwork calculate a hyperG pval for the intersection
		//n = total number of genes
		int n = uniqueInterrogatedGenes.size();
		//a = number of genes in set A
		int a = uniqueSelectGenes.size();
		for (int i=0; i< analyzedNetworks.length; i++) {
			//b = number of genes in set B
			//int b = p.genes.size();
			int b = analyzedNetworks[i].getGeneNetworkGeneNamesInInterrogatedGenes().size();
			//t = number of genes in common to A and B
			//int t = p.sharedGenes.size();
			int t = analyzedNetworks[i].getGeneSharedGenes().size();
			nabt[i] = new int[]{n,a,b,t};
			//if (i<30) lg(analyzedNetworks[i].getKeggApiNetwork().getNetworkName()+"\t"+Num.intArrayToString(nabt[i], "-"));
		}
		//Calculate pvalues using hg dist, fast way to do a fisher's exact on a 2x2 table
		IntersectListsHypergeometric ih = new IntersectListsHypergeometric(tempDirectory, fullPathToR);
		double[] pvals = ih.calculateOverRepresentationPValues(nabt);
		//for (int i=0; i< pathways.length; i++) pathways[i].pval = pvals[i];
		for (int i=0; i< analyzedNetworks.length; i++) {
			analyzedNetworks[i].setPValue(pvals[i]);
			analyzedNetworks[i].setSortValue(pvals[i]);
		}
		//convert the pvals to fdrs
		Arrays.sort(analyzedNetworks); //smallest to largest
		double[] pvalsLargeToSmall = new double[analyzedNetworks.length];
		int counter = 0;
		for (int i=analyzedNetworks.length-1; i>=0; i--) pvalsLargeToSmall[counter++] = analyzedNetworks[i].getPValue();
		Num.benjaminiHochbergCorrect(pvalsLargeToSmall);
		counter=0;
		for (int i=analyzedNetworks.length-1; i>=0; i--) analyzedNetworks[i].setFdr(pvalsLargeToSmall[counter++]); 
	}

	public static void main(String[] args) throws IOException {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KeggGenePathwayAnalyzer(args);
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
					case 'i': interrogatedGeneList = new File(args[++i]); break;
					case 'k': keggIdsFile = new File(args[++i]); break;
					case 'g': selectGeneList = new File(args[++i]); break;
					case 'e': fullPathToR = new File(args[++i]); break;
					case 'n': keggNetworkDirectory = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'x': maximumFdr = Double.parseDouble(args[++i]); break;
					case 'm': minimumNumberGenes = Integer.parseInt(args[++i]); break;
					case 't': typesToExclude = args[++i]; break;
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
		lg("Checking required files and directories...");
		File[] toCheck = {keggIdsFile, interrogatedGeneList, selectGeneList, fullPathToR, keggNetworkDirectory, resultsDirectory};
		String[] names = {"-k keggIdsFile", "-i interrogatedGeneList", "-g selectGeneList", "-e fullPathToR", "-n keggNetworkDir", "-r resultsDir"};
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
				"**                       Kegg Gene Pathway Analyzer : Feb 2026                      **\n" +
				"**************************************************************************************\n" +
				"KGPA uses your interrogated gene list to 1) filter your select gene list and the KEGG\n"+
				"Network gene sets then 2) intersects these lists, 3) calculates hypergeometric \n"+
				"p-values, 4) controls for multiple testing (Benjamini-Hochberg FDR method), and\n"+
				"5) outputs two Excel spreadsheets focused on each tested KEGG Network (sub pathways)\n"+
				"and KEGG Pathway. Use the URL links in the spreadsheets to interactively explort the\n"+
				"results in the KEGG Pathway Viewer (https://www.kegg.jp).\n"+
				
				"\nRelated Apps:\n\n"+
				
				"KeggGeneSymbolIdExtractor - Uses the Kegg API to look up and parse Kegg Gene Ids that\n"+
				"   match each of the provided HUGO Gene Symbols. Use this tool to generate the\n"+
				"   required '-k KEGG gene Id HUGO gene symbol lookup file'.\n"+
				"KeggResourceExtractor - Pulls the latest list of KEGG Networks, their associated\n"+
				"   genes, and pathways. Use this to generate the required '-n KEGG network dir'.\n"+
				"KeggVariantPathwayAnalyzer - Differential gene mutation cohort KEGG Network and\n"+
				"   Pathway analysis.\n"+
				"KeggGeneAndVariantPathway Analyzer - Runs a joint gene expression and gene mutation\n"+
				"   KEGG Network and Pathway analysis. Recommended if both available.\n"+
				"MergeKeggNetworkResults - Merges network xls results from multiple pathway analysis.\n"+
				"DESeq2, edgeR - R packages for selecting differentially expressed gene sets.\n"+
				
				"\nKEGG Viewer Gene Color Key:\n\n"+
				
				"   Negative Diff Exp LgRto - Light Blue - #87CEEB\n"+
				"   Positive Diff Exp LgRto - Light Red, Pink - #FFB6C1\n"+
		        "   Genes in red text are disease associated\n"+

				"\nApp Parameters:\n\n"+
				
				"-i File containing all of the interrogated genes in your study, one gene per line.\n"+
				"     Typically > 15K gene HUGO gene symbols.\n"+
				"-g File containing the selected genes of interest with their log2Rto, one per line,\n"+
				"     tab delimited, from your study, e.g. differentially expressed, typically <2K\n"+
				"-e File path to the R executable\n\n"+
				
				"-r Directory to save the spreadsheet results\n"+
				"-k KEGG gene Id HUGO gene symbol lookup file, see above\n"+
				"-n KEGG network directory, see above\n"+
				"-t (Optional) Network TYPEs to exclude from testing, comma delimited no spaces.\n"+
				"-m (Optional) Minimum number of interrogated genes in a network for analysis,\n"+
				"     defaults to 4\n"+
				"-x (Optional) Maximum FDR for including networks into the combine Kegg Pathway\n"+
				"     spreadsheet, defaults to 0.15\n"+

				"\nExample:\n\n"+ 
				"java -Xmx1G -jar pathTo/USeq/Apps/KeggGenePathwayAnalyzer -i \n"+
				"   allTestedGenes.txt -g diffExpGenes.txt -k geneSymbol2KeggGeneInfo.txt.gz -n Networks/\n"+
				"   -r GeneAnalyzerResults -t 'Pathogen,Ev factor' -m 5 -e /usr/bin/R\n"+

				"\n**************************************************************************************\n");

	}


	public boolean isFailed() {
		return failed;
	}


	public ArrayList<String> getLog() {
		return log;
	}


	public KeggApiNetwork[] getAllNetworks() {
		return allNetworks;
	}


	public AnalyzedNetwork[] getAnalyzedNetworks() {
		return analyzedNetworks;
	}


	public HashMap<String, ArrayList<String>> getGs2ki() {
		return gs2ki;
	}

}
