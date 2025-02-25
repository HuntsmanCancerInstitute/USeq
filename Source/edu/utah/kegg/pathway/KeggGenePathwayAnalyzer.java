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

public class KeggGenePathwayAnalyzer {

	//user fields
	private File interrogatedGeneList;
	private File selectGeneList;
	private File keggNetworkDirectory;
	private File keggIdsFile;
	private File fullPathToR;
	private File resultsDirectory = null;
	private double maximumFdr = 0.15;
	private boolean downloadSigPathways = false;  //not working so keep off
	private int minimumNumberGenes = 4;
	private HashSet<String> networkTypesToExclude = null;
	private String typesToExclude = null;

	//internal
	private File tempDirectory = null;
	private HashMap<String, AnalyzedNetwork> networkIdAnalyzedNetwork = null;
	private AnalyzedNetwork[] analyzedNetworks = null;
	private HashMap<String, ArrayList<String>> gs2ki = null;
	private HashMap<String, ArrayList<String>> ki2gs = null;
	private HashSet<String> uniqueInterrogatedGenes = null;
	private HashSet<SelectGene> uniqueSelectGenes = null;
	
	//constructor
	public KeggGenePathwayAnalyzer(String[] args) {
		try {
			//start clock
			long startTime = System.currentTimeMillis();

			//process args
			processArgs(args);
			checkFiles();

			IO.pl("\nLoading KEGG gene symbol <-> id lookup tables...");
			loadKeggIdLookupHashes();

			IO.pl("\nLoading your interrogated gene symbols...");
			loadInterrogatedGenes();

			IO.pl("\nLoading your select genes...");
			loadSelectedGenes();

			//parse network info
			IO.pl("\nLoading and filtering KEGG Medicus networks...");
			loadKeggApiNetworks();

			IO.pl("\nIntersect select genes with networks...");
			intersectNetworks();

			IO.pl("\nComparing networks...");
			compareNetworks();

			IO.pl("\nSaving all network results...");
			saveNetworks();

			IO.pl("\nSaving significant pathways...\n");
			savePathways();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
		} catch (Exception e) {
			e.printStackTrace();
			IO.el("ERROR running the KeggGenePathwayAnalyzer!");
		}
	}

	private void savePathways() throws IOException {
		//the idea here is to find all of the significant networks with a pathway
		//	then for each pathway decorate it with all of the network genes colored, and highlight the networks
		//	really only useful when multiple diff networks hit the same pathway
		//		output pathway view xls file

		//find significant networks 
		//pathway id and the ANs it came from
		TreeMap<String, ArrayList<AnalyzedNetwork>> sigPath = new TreeMap<String, ArrayList<AnalyzedNetwork>>();
		HashMap<String, KeggApiPathway> pathwayIdKeggApiPathway = new HashMap<String, KeggApiPathway>();
		for (AnalyzedNetwork an: analyzedNetworks) {
			if (an.getGeneAdjPVal()<= maximumFdr) {
				//for each network in the an
				for (KeggApiNetwork kan: an.getGeneAnalizedNetworks()){
					//any pathways?
					if (kan.getPathways() != null) {
						for (KeggApiPathway kap: kan.getPathways()) {
							pathwayIdKeggApiPathway.put(kap.getId(), kap);
							ArrayList<AnalyzedNetwork> al = sigPath.get(kap.getId());
							if (al == null) {
								al = new ArrayList<AnalyzedNetwork>();
								sigPath.put(kap.getId(), al);
							}
							al.add(an);
						}
					}
				}
			}
		}

		//set up pathmap dir? seems to be blocked by kegg so not working!
		File pathwayDir = null;
		ArrayList<String> downloadCmds = new ArrayList<String>();
		ArrayList<File> downloadFiles = new ArrayList<File>();
		if (sigPath.size() !=0 && downloadSigPathways) {
			pathwayDir = new File(resultsDirectory, "PathwayImageMaps");
			pathwayDir.mkdirs();
		}

		//for each pathway, create a results obj to sort by pvalue
		StringValueSort[] results = new StringValueSort[sigPath.size()];
		int counter = 0;
		for (String pathwayId: sigPath.keySet()) {
			StringBuilder xls = new StringBuilder();
			xls.append("PATHWAY:\t");
			xls.append(pathwayId); 
			xls.append("\t"); 
			xls.append(pathwayIdKeggApiPathway.get(pathwayId).getName()); 
			xls.append("\t");  
			xls.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);		
			xls.append("\",\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);	
			xls.append("\")\n");

			//fetch all of the genes and networks
			TreeMap<String, SelectGene> allGenes = new TreeMap<String,SelectGene>();
			HashSet<String> networkIds = new HashSet<String>();
			xls.append("Networks:\n");
			double minPVal = Double.MAX_VALUE;
			for (AnalyzedNetwork an : sigPath.get(pathwayId)) {
				if (an.getGenePVal()< minPVal) minPVal = an.getGenePVal();
				for (SelectGene sg: an.getGeneSharedGenes()) allGenes.put(sg.getGeneSymbol(), sg);
				for (KeggApiNetwork net: an.getGeneAnalizedNetworks()) networkIds.add(net.getNetworkId());

				xls.append("\t"); xls.append(an.getGeneAnalizedNetworkIdNames()); xls.append("\n");
				xls.append("\t\tAdjPval: "); xls.append(Num.formatNumber(an.getGeneAdjPVal(),3)); xls.append("\n");
				xls.append("\t\tGenes: "); xls.append(an.getSharedGeneSymbols()); xls.append("\n");

			}
			SelectGene[] sg = new SelectGene[allGenes.size()];
			int i = 0;
			for (SelectGene g: allGenes.values()) sg[i++] = g;
			xls.append("AllGenes:\t"); xls.append(Misc.stringSetToString(allGenes.keySet(), ", ")); xls.append("\n");
			
			//networks
			String[] networkIdsStringArray = Misc.hashSetToStringArray(networkIds);
			//"skyblue" for negColor, "pink" for posColor and zeroColor
			String url = AnalyzedNetwork.fetchKeggPathwayMapLink(pathwayId, networkIdsStringArray, sg, gs2ki,AnalyzedNetwork.geneNegColor,AnalyzedNetwork.genePosColor,AnalyzedNetwork.genePosColor);
			xls.append("KeggLink:\t"); xls.append(url); xls.append("\n");
			if (url.length()<250) {
				xls.append("ExcelLink:\t=HYPERLINK(\"");
				xls.append(url);		
				xls.append("\",\""+url);	
				xls.append("\")\n\n");
			}
			else xls.append("ExcelLink:\tToo big\n\n");
			results[counter++] = new StringValueSort(xls, minPVal);
			//download the network image? not working, seems blocked by KEGG?
			/*
			if (downloadSigPathways) {
				File pathwayImage = new File(pathwayDir, pathwayId+"_Fdr"+ minimumFdr+".png");
				String cmd = "curl -s -o "+pathwayImage +" '"+url+ "#downloadImage2x' &";
				downloadFiles.add(pathwayImage);
				downloadCmds.add(cmd);
			}*/
		}

		//executeCmds(downloadCmds, downloadFiles);
		
		Arrays.sort(results);
		PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "genePathwaysMinGen"+minimumNumberGenes+"MaxFdr"+maximumFdr+".xls")));
		//name descLink pval, adjPval, #UniqueNetworkGenes, #FoundUniqueNetworkGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/PathGenes, IntersectingGenes
		out.println("#Composite view of significant networks in KEGG pathways\n");
		for (StringValueSort s: results) out.print(s.getCargo().toString());
		out.close();
	}

	/*
	private void executeCmds(ArrayList<String> downloadCmds, ArrayList<File> downloadFiles) throws IOException {
		if (downloadCmds.size()==0) return;
		File complete = new File(downloadFiles.get(0).getParent(), "COMPLETE");
		downloadFiles.add(complete);
		downloadCmds.add("sleep 5");
		downloadCmds.add("touch "+complete);
		String shellScript = "set -e \n" + Misc.stringArrayListToString(downloadCmds, "\n");
		IO.pl("Executing:\n"+shellScript);
		IO.executeShellScript(shellScript, tempDirectory);
		//look for files
		for (File f: downloadFiles) {
			if (f.exists()==false) throw new IOException("ERROR: problem executing\n"+shellScript);
		}
	}*/

	
	private void loadKeggIdLookupHashes() throws IOException {
		HashMap<String, ArrayList<String>>[] hashes = KeggGeneSymbolIdExtractor.loadGeneLookupHashes(keggIdsFile);
		gs2ki = hashes[0];
		ki2gs = hashes[1];
	}

	private void loadKeggApiNetworks() throws IOException {
		KeggApiNetwork[] networks = KeggResourceExtractor.loadNetworks(keggNetworkDirectory);
		networkIdAnalyzedNetwork = new HashMap<String, AnalyzedNetwork>(networks.length);
		for (int i=0; i< networks.length; i++) {
			//check type?
			if (networkTypesToExclude != null) {
				String type = networks[i].getNetworkType().toLowerCase();
				if (networkTypesToExclude.contains(type)) continue; 
			}

			//check minimum # genes
			KeggApiGene[]  genes = networks[i].getGenes();
			if (genes!=null && genes.length >= minimumNumberGenes) {
				networkIdAnalyzedNetwork.put(networks[i].getNetworkId(), new AnalyzedNetwork(networks[i], uniqueSelectGenes.size()));
			}
		}
		IO.pl("\t"+networkIdAnalyzedNetwork.size()+"\tNetworks loaded that pass minimum # genes ("+minimumNumberGenes+") and excluded types: "+typesToExclude+"\n");

		//for each network gene, look for them in the user's interrogated genes
		for (AnalyzedNetwork an: networkIdAnalyzedNetwork.values()) {
			TreeSet<String> found = an.getGeneNetworkGeneNamesInInterrogatedGenes();
			for (KeggApiGene kag: an.getKeggApiNetwork().getGenes()) {
				if (uniqueInterrogatedGenes.contains(kag.getName())) found.add(kag.getName());
			}
			int numNetGenes = an.getKeggApiNetwork().getGenes().length;
			int numFound = found.size();
			if (numFound < minimumNumberGenes) {
				IO.pl("\tSkipping '"+an.getKeggApiNetwork().getNetworkIdName()+"' too few of its genes ("+numFound+") intersected your interrogated gene list.");
				an.setGeneAnalyzeIt(false);
			}
			else if (numFound != numNetGenes) {
				//IO.pl("\t"+numFound+" of "+numNetGenes+" were found in your interrogated gene list for '"+ an.getKeggApiNetwork().getNetworkName()+"'");
			}
		}

		//merge networks since the # genes has changed, don't want to test these twice
		IO.pl("\nMerging networks with the same gene sets...");
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

		//make final set to analyze
		analyzedNetworks = new AnalyzedNetwork[geneNameStringNetwork.size()];
		int index = 0;
		for (String genes: geneNameStringNetwork.keySet()) {
			ArrayList<AnalyzedNetwork> al = geneNameStringNetwork.get(genes);
			//IO.pl("\t"+ al.size()+"\tNetworks with "+genes);
			analyzedNetworks[index] = al.get(0);
			//add in networks
			ArrayList<KeggApiNetwork> n = analyzedNetworks[index].getGeneAnalizedNetworks(); 
			for (AnalyzedNetwork an: al) n.add(an.getKeggApiNetwork());
			index++;
		}
		IO.pl("\t"+analyzedNetworks.length+"\tGene sets to test");
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
		if (failedToFind.size()!=0) IO.pl("\tWARNING: failed to find "+failedToFind.size()+" in the interrogated gene list, skipping: "+Misc.stringArrayListToString(failedToFind, ","));
		IO.pl("\t"+selectGenes.length+" -> "+uniqueSelectGenes.size()+ " unique selected genes");

	}

	private void loadInterrogatedGenes() {
		String[] allinterrogatedGenes = IO.loadFileIntoStringArray(interrogatedGeneList);
		uniqueInterrogatedGenes = new HashSet<String>();
		for (String s: allinterrogatedGenes) {
			if (gs2ki.containsKey(s)) uniqueInterrogatedGenes.add(s);
		}
		IO.pl("\t"+allinterrogatedGenes.length+" -> "+uniqueInterrogatedGenes.size()+ " unique HUGO gene symbols that are recognized by KEGG.");
	}

	private void saveNetworks() throws IOException {
		PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "geneNetworksMinGen"+minimumNumberGenes+".xls")));
		//name descLink pval, adjPval, #UniqueNetworkGenes, #FoundUniqueNetworkGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/PathGenes, IntersectingGenes
		out.println("# Network Name(s)\tNetwork Desc Link\tPval\tAdj Pval\tAll Network Genes\tFound Network Genes\tSelect Genes\tIntersect\tInt/Net\tInt Gene Symbols\tInt Genes LgRtos\tPathway Map Links...");
		for (AnalyzedNetwork an: analyzedNetworks) out.print(an.toStringGene(gs2ki));
		out.close();
	}


	private void intersectNetworks() {
		//for each AnalyzedNetwork
		for (AnalyzedNetwork p: networkIdAnalyzedNetwork.values()) {

			//pull the genes from the network after filtered against the users interrogated gene list
			TreeSet<String> networkGeneNames = p.getGeneNetworkGeneNamesInInterrogatedGenes();
			ArrayList<SelectGene> sharedGenes = p.getGeneSharedGenes();

			//for each select gene from the user, if it intersects the network genes save it.
			Iterator<SelectGene> it = uniqueSelectGenes.iterator();
			while (it.hasNext()) {
				SelectGene sg = it.next();
				if (networkGeneNames.contains(sg.getGeneSymbol())) sharedGenes.add(sg);
			}
		}
	}

	private void compareNetworks() throws IOException {
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

			//if (i<30) IO.pl(analyzedNetworks[i].getKeggApiNetwork().getNetworkName()+"\t"+Num.intArrayToString(nabt[i], "-"));
		}

		//Calculate pvalues using hg dist, fast way to do a fisher's exact on a 2x2 table
		IntersectListsHypergeometric ih = new IntersectListsHypergeometric(tempDirectory, fullPathToR);
		double[] pvals = ih.calculatePValues(nabt);
		//for (int i=0; i< pathways.length; i++) pathways[i].pval = pvals[i];
		for (int i=0; i< analyzedNetworks.length; i++) {
			analyzedNetworks[i].setGenePVal(pvals[i]);
			analyzedNetworks[i].setSortValue(pvals[i]);
		}

		//convert the pvals to fdrs
		Arrays.sort(analyzedNetworks); //smallest to largest
		double[] pvalsLargeToSmall = new double[analyzedNetworks.length];
		int counter = 0;
		for (int i=analyzedNetworks.length-1; i>=0; i--) pvalsLargeToSmall[counter++] = analyzedNetworks[i].getGenePVal();
		Num.benjaminiHochbergCorrect(pvalsLargeToSmall);
		counter=0;
		for (int i=analyzedNetworks.length-1; i>=0; i--) analyzedNetworks[i].setGeneAdjPVal(pvalsLargeToSmall[counter++]); 
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
		IO.pl("Checking required files and directories...");
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
				"**                      Kegg Gene Set Pathway Analyzer : Feb 2025                   **\n" +
				"**************************************************************************************\n" +
				"KGSPA uses your interrogated gene list to 1) filter your select gene list and the KEGG\n"+
				"Network gene sets then 2) intersects these lists, 3) calculates hypergeometric \n"+
				"p-values, 4) controls for multiple testing using the Benjamini-Hochberg FDR method, and\n"+
				"5) output two Excel spreadsheets focused on each tested KEGG Network (sub pathways)\n"+
				"and KEGG Pathway. Use the URL links in the spreadsheets to color code genes and gene\n"+
				"pathways in the interacive KEGG Pathway Viewer (https://www.kegg.jp). Run the\n"+
				"USeq KeggJointPathwayAnalyzer to combine these analysis with those from the USeq\n"+
				"KeggVariantPathwayAnalyzer.\n"+

				"\nRequired Parameters:\n"+
				"-i File containing all of the interrogated genes in your study, one gene per line.\n"+
				"-g File containing the selected genes of interest with their log2Rto, one per line,\n"+
				"     tab delimited, from your study.\n"+
				"-k KEGG gene Id HUGO gene symbol lookup file generated by running the USeq\n"+
				"     KeggGeneSymbolIdExtractor.\n"+
				"-n KEGG network directory created by running the USeq KeggResourceExtractor.\n"+
				"-r Directory to save the spreadsheet results.\n"+
				"-t Network TYPEs to exclude from testing, comma delimited no spaces.\n"+
				"-m Minimum number of interrogated genes in a network for analysis, defaults to 5\n"+
				"-x Maximum FDR for including networks into the Kegg Pathway spreadsheet, defaults\n"+
				"     to 0.15\n"+
				"-e File path to the R executable.\n"+

				"\nExample: java -Xmx1G -jar pathTo/USeq/Apps/GeneSetPathwayComparator -i \n"+
				"   allTestedGenes.txt -g diffExpGenes.txt -k keggIdLookup.txt -n KeggNetworks/ \n"+
				"   -r GenePathwayAnalyzerResults -t 'Pathogen,Ev factor' -e /usr/bin/R\n"+

				"\n**************************************************************************************\n");

	}
}
