package edu.utah.kegg.pathway;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.*;
import edu.utah.kegg.api.KeggApiNetwork;
import util.gen.*;

public class KeggGeneAndVariantPathwayAnalyzer {

	//user fields
	private File keggNetworkDirectory;
	private File keggIdsFile;
	private File resultsDirectory = null;
	private File interrogatedGeneList;
	private File selectGeneList;
	private File fullPathToR;
	private File groupAGeneHits;
	private File groupBGeneHits;
	private int minimumNumberGenes = 4;
	private HashSet<String> networkTypesToExclude = null;
	private String typesToExclude = null;
	private boolean addOne = false;
	private File tempDirectory = null;
	private double maximumFdr = 0.15;
	private boolean replaceNetworksWithPathways = false;
	
	private KeggGenePathwayAnalyzer genes = null;
	private KeggVariantPathwayAnalyzer variants = null;
	private HashMap<String, ArrayList<String>> gs2ki = null;
	
	//constructor
	public KeggGeneAndVariantPathwayAnalyzer(String[] args){
		try {
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		checkFiles();
		
		IO.pl("\nLaunching gene and variant pathway analysis...");
		runAnalysis();
		gs2ki = genes.getGs2ki();

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
		CombinePathwayRoot cpr = new CombinePathwayRoot(genes.getAnalyzedNetworks(), variants.getAnalyzedNetworks(), "");
		cpr.makeCombinePathways(maximumFdr);
		cpr.saveGeneAndVariantPathways(maximumFdr, gs2ki, resultsDirectory, minimumNumberGenes);
	}

	private void runAnalysis() throws IOException {
		//create the objects
		boolean verbose = false;
		genes = new KeggGenePathwayAnalyzer(keggIdsFile, keggNetworkDirectory, fullPathToR, resultsDirectory, tempDirectory, interrogatedGeneList, 
				selectGeneList, maximumFdr, minimumNumberGenes, typesToExclude, networkTypesToExclude, verbose, replaceNetworksWithPathways, "");

		variants = new KeggVariantPathwayAnalyzer(keggIdsFile, keggNetworkDirectory, resultsDirectory, groupAGeneHits, groupBGeneHits, minimumNumberGenes, 
				maximumFdr, addOne, typesToExclude, networkTypesToExclude, verbose, replaceNetworksWithPathways);
		
		//run the analysis in separate threads
		ExecutorService executor = Executors.newFixedThreadPool(2);
		executor.execute(variants);
		executor.execute(genes);
		executor.shutdown();

		//spins here until the executer is terminated, e.g. all threads complete
		while (!executor.isTerminated()) {}
		
		//write out logs
		File geneLog = new File (resultsDirectory, "geneAnalysisLog.txt");
		IO.writeArrayList(genes.getLog(), geneLog);
		File variantLog = new File (resultsDirectory, "variantAnalysisLog.txt");
		IO.writeArrayList(variants.getLog(), variantLog);
		
		//check both 
		if (genes.isFailed() || variants.isFailed()) throw new IOException();
	}




	private void saveNetworks() {
		
		//make a hash of analyzed networks
		AnalyzedNetwork[] geneAnNet = genes.getAnalyzedNetworks();
		AnalyzedNetwork[] variantAnNet = variants.getAnalyzedNetworks();
		TreeMap<String, AnalyzedNetwork[]> netIdGenVarAnNet = new TreeMap<String, AnalyzedNetwork[]>();
		
		for (AnalyzedNetwork an: geneAnNet) {
			ArrayList<KeggApiNetwork> networks = an.getAnalizedKeggApiNetworks();
			for (KeggApiNetwork n: networks) {
				AnalyzedNetwork[] geneVar = netIdGenVarAnNet.get(n.getNetworkId());
				if (geneVar == null) {
					geneVar = new AnalyzedNetwork[2];
					netIdGenVarAnNet.put(n.getNetworkId(), geneVar);
				}
				geneVar[0]=an;
			}
		}
		
		for (AnalyzedNetwork an: variantAnNet) {
			ArrayList<KeggApiNetwork> networks = an.getAnalizedKeggApiNetworks();
			for (KeggApiNetwork n: networks) {
				AnalyzedNetwork[] geneVar = netIdGenVarAnNet.get(n.getNetworkId());
				if (geneVar == null) {
					geneVar = new AnalyzedNetwork[2];
					netIdGenVarAnNet.put(n.getNetworkId(), geneVar);
				}
				geneVar[1]=an;
			}
		}
		
		//combine the pvalues
		CombinePValues cp = new CombinePValues();
		double[] toCombine = new double[2];
		ArrayList<Float> toCorrect = new ArrayList<Float>();

		for (String netId: netIdGenVarAnNet.keySet()) {
			AnalyzedNetwork[] geneVar = netIdGenVarAnNet.get(netId);
			if (geneVar[0]!=null && geneVar[1]!=null) {
				toCombine[0] = geneVar[0].getPValue();
				toCombine[1] = geneVar[1].getPValue(); 
				//watch out with 1's can cause non convergence
				double combinePVal = 1.0;
				if (toCombine[0]!=1.0 && toCombine[1]!=1.0) combinePVal = cp.calculateCombinePValues(toCombine);
				geneVar[0].setCombinePValue(combinePVal);
				geneVar[1].setCombinePValue(combinePVal);
				toCorrect.add(Num.minus10log10Float(combinePVal));
			}
		}
		
		//correct the pvals and add back
		float[] f = Num.benjaminiHochbergCorrectUnsorted(Num.arrayListOfFloatToArray(toCorrect));
		double[] fdrs = Num.antiNeg10log10(f);
		
		int counter = 0;
		for (String netId: netIdGenVarAnNet.keySet()) {
			AnalyzedNetwork[] geneVar = netIdGenVarAnNet.get(netId);
			if (geneVar[0]!=null && geneVar[1]!=null) {
				geneVar[0].setCombineFdr(fdrs[counter]);
				geneVar[1].setCombineFdr(fdrs[counter]);
				counter++;
			}
		}
		
		//print out a combine analysis for those hitting both
		try {
			PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "combineNetworksMinGen"+minimumNumberGenes+".xls")));
			//genes
			//name descLink pval, adjPval, #UniqueNetworkGenes, #FoundUniqueNetworkGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/PathGenes, IntersectingGenes
			//out.println("# Gene Network Name(s)\tNetwork Desc Link\tPval\tAdj Pval\tAll Network Genes\tFound Network Genes\tSelect Genes\tIntersect\tInt/Net\tInt Gene Symbols\tInt Genes LgRtos\tPathway Map Links...");
			//for (AnalyzedNetwork an: analyzedNetworks) out.print(an.toStringGene(gs2ki));
			
			//variants
			//out.println("# Var Network Name(s)\tNetwork Desc Link\tPval\tAdjPval\tAHits\tANoHits\tFracAHits\tAGeneHits\tBHits\tBNoHits\tFracBHits\tBGeneHits\tLog2(fracA/fracB)\tAllGeneHits\tPathwayMapLinksWithTopMapDescription...");
			//for (AnalyzedNetwork an: analyzedNetworks) out.print(an.toStringVariant(addOne, gs2ki));
			
			out.println("# Network Name\tNetwork Desc Link\tCombine Pval\tCombine Adj Pval\tGene Pval\tGene Adj Pval\tGene #IntGenes/ #NetGenes\tInt Gene Symbols\tVariant Pval\tVariant AdjPval\tVariant Log2(fracA/fracB)\tVariant AllGeneHits\tPathway Map Links...");
			
			
			for (String netId: netIdGenVarAnNet.keySet()) {
				AnalyzedNetwork[] geneVar = netIdGenVarAnNet.get(netId);
				if (geneVar[0]!=null && geneVar[1]!=null) {
					out.println(netId+"\t"+geneVar[0].getCombinePValue()+"\t"+geneVar[0].getCombineFdr());
					out.print(geneVar[0].toStringGeneCombo(netId, gs2ki));
					out.print(geneVar[1].toStringVariantCombo(netId, geneVar[0], addOne, gs2ki));
					out.println();
				}
			}
			
			out.close();
			IO.pl("Combo print complete");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	

	



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KeggGeneAndVariantPathwayAnalyzer(args);
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
					case 'a': groupAGeneHits = new File(args[++i]); break;
					case 'b': groupBGeneHits = new File(args[++i]); break;
					case 'k': keggIdsFile = new File(args[++i]); break;
					case 'n': keggNetworkDirectory = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'm': minimumNumberGenes = Integer.parseInt(args[++i]); break;
					case 't': typesToExclude = args[++i]; break;
					case 'x': maximumFdr = Double.parseDouble(args[++i]); break;
					case 'o': addOne = true; break;
					case 'i': interrogatedGeneList = new File(args[++i]); break;
					case 'g': selectGeneList = new File(args[++i]); break;
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
		File[] toCheck = {keggIdsFile, keggNetworkDirectory, interrogatedGeneList, selectGeneList, fullPathToR, resultsDirectory, groupAGeneHits, groupBGeneHits};
		String[] names = {"-k keggIdsLookupFile", "-n keggNetworkDir", "-i interrogatedGeneList", "-g selectGeneList", "-e fullPathToR", "-r resultsDir", "-a groupAGeneHitsFile", "-b groupBGeneHitsFile"};
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
				"**                 Kegg Gene and Variant Pathway Analyzer : Feb 2026                **\n" +
				"**************************************************************************************\n" +
				"Runs both the KeggGenePathwayAnalyzer and KeggVariantPathwayAnalyzer applications.\n"+
				"Combines the results at the KEGG Pathway level, coloring each gene for interacive\n"+
				"exploration in the KEGG Pathway Viewer (https://www.kegg.jp). Calculates\n"+
				"a combine Pathway p-value (Fisher's method) and FDR (Benjamini-Hochberg) for Pathways\n"+
				"with multiple significant Networks. Best to use this tool if you are comparing two\n"+
				"large cohorts with both differential gene expression and somatic mutation datasets.\n"+
				
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
				"   Negative Diff Exp LgRto - Light Blue - #87CEEB\n"+
				"   Positive Diff Exp LgRto - Light Red, Pink - #FFB6C1\n"+
				"   Negative Mutation Freq LgRto - Light Violet - #D6B4FC\n"+
				"   Positive Mutation Freq LgRto - Light Orange - #FFD580\n"+
				"   Abs(Mut Freq) < 1.25x - Light Yellow - #FFFFAC\n"+
		        "   Genes in red text are disease associated\n"+
				
				"\nApp Parameters:\n\n"+
				"-i File containing all of the interrogated genes in your gene expression study, one\n"+
				"     gene per line. Typically > 15K genes\n"+
				"-g File containing the selected genes of interest with their log2Rto, one per line,\n"+
				"     tab delimited, from your study, e.g. differentially expressed, typically <2K\n"+
				"-e File path to the R executable\n\n"+
				
				"-a File containing cohort A gene sets, each line represents a subject's genes of\n"+
				"     interest (e.g. those with HIGH impact mutations), tab delimited, the first cell\n"+
				"     is the subject ID, subsequent cells are the gene names.\n"+
				"-b File containing cohort B gene sets, ditto.\n"+
				"-o (Optional) Add one to zero count A or B fractions when calculating the\n"+
				"     log2Rto(fracA/fracB)\n\n"+
				
				"-r Directory to save the spreadsheet results\n"+
				"-k KEGG gene Id HUGO gene symbol lookup file, see above\n"+
				"-n KEGG network directory, see above\n"+
				"-t (Optional) Network TYPEs to exclude from testing, comma delimited no spaces.\n"+
				"-m (Optional) Minimum number of interrogated genes in a network for analysis,\n"+
				"     defaults to 4\n"+
				"-x (Optional) Maximum FDR for including networks into the combine Kegg Pathway\n"+
				"     spreadsheet, defaults to 0.15\n"+
				"-p Replace networks with composite pathways, thus no network analysis.\n"+
	
				"\nExample:\n\n"+
				"java -Xmx1G -jar pathTo/USeq/Apps/KeggGeneAndVariantPathwayAnalyzer\n"+
				"   -i allTestedGenes.txt -g diffExpGenes.txt -a earlyCRC.txt -b lateCRC.txt \n"+
				"   -k geneSymbol2KeggGeneInfo.txt.gz -n Networks/ -o -r CombineResults -t\n"+
				"   'Pathogen,Ev factor'\n"+

		"\n**************************************************************************************\n");
	}
}
