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

public class KeggVariantPathwayAnalyzer {

	//user fields
	private File groupAGeneHits;
	private File groupBGeneHits;
	private boolean addOne = false;
	
	//new stuff
	private double maximumFdr = 0.15;
	private HashMap<String, ArrayList<String>> gs2ki = null;
	private HashMap<String, ArrayList<String>> ki2gs = null;
	private File keggNetworkDirectory;
	private File keggIdsFile;
	private File resultsDirectory = null;
	private int minimumNumberGenes = 4;
	private HashSet<String> networkTypesToExclude = null;
	private String typesToExclude = null;
	private HashMap<String, AnalyzedNetwork> networkIdAnalyzedNetwork = null;
	private AnalyzedNetwork[] analyzedNetworks = null;
	private TreeMap<String, Integer> missingGeneSymbols = new TreeMap<String, Integer>();
	private TreeMap<String, Integer> foundGeneSymbols = new TreeMap<String, Integer>();
	
	
	//constructor
	public KeggVariantPathwayAnalyzer(String[] args){
		try {
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		checkFiles();
		
		IO.pl("\nLoading KEGG HUGO gene symbol <-> id lookup tables...");
		loadKeggIdLookupHashes();
		
		//parse network info
		IO.pl("\nLoading and filtering KEGG Medicus networks...");
		loadKeggApiNetworks();
		
		IO.pl("\nLoading group A gene hits...");
		loadGroup(groupAGeneHits, true);
		
		IO.pl("\nLoading group B gene hits...");
		loadGroup(groupBGeneHits, false);
		
		IO.pl("\nFilter network hits...");
		filterNetworks();
		
		IO.pl("\nComparing pathway hits...");
		comparePathways();
		
		IO.pl("\nSaving individual network results...");
		saveNetworks();
		
		IO.pl("\nSaving combine pathway results...");
		savePathways();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
		
		} catch (Exception e) {
			e.printStackTrace();
			IO.el("ERROR running the KeggVariantPathwayAnalyzer!");
		}
	}
	
	private void filterNetworks() {
		HashMap<String, ArrayList<AnalyzedNetwork>> geneNameKeyNetwork = new HashMap<String, ArrayList<AnalyzedNetwork>>();


		//first toss any with no hits
		for (AnalyzedNetwork an: networkIdAnalyzedNetwork.values()) {
			//any hits?
			if (an.getGroupAHits().size() > 0 || an.getGroupBHits().size() > 0) {
				TreeSet<String> found = an.getVarinatGeneNameHits();
				String key = Misc.treeSetToString(found, ",");
				ArrayList<AnalyzedNetwork> al = geneNameKeyNetwork.get(key);
				if (al == null) {
					al = new ArrayList<AnalyzedNetwork>();
					geneNameKeyNetwork.put(key, al);
				}
				al.add(an);
			}
		}

		//second merge networks with the same gene name hits, 
		//several networks sets are variants of the same genes so don't want to score them separately if the gene hits are the same
		analyzedNetworks = new AnalyzedNetwork[geneNameKeyNetwork.size()];
		int index = 0;
		for (String genes: geneNameKeyNetwork.keySet()) {
			ArrayList<AnalyzedNetwork> al = geneNameKeyNetwork.get(genes);
			//IO.pl("\t"+ al.size()+"\tNetworks with "+genes);
			analyzedNetworks[index] = al.get(0);
			//add in networks
			ArrayList<KeggApiNetwork> n = analyzedNetworks[index].getVariantAnalizedNetworks(); 
			for (AnalyzedNetwork an: al) n.add(an.getKeggApiNetwork());
			index++;
		}

		IO.pl("\t"+analyzedNetworks.length+" networks with dataset hits from "+networkIdAnalyzedNetwork.size()+" total.");
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
				//networkIdAnalyzedNetwork.put(networks[i].getNetworkId(), new AnalyzedNetwork(networks[i], uniqueSelectGenes.size()));
				networkIdAnalyzedNetwork.put(networks[i].getNetworkId(), new AnalyzedNetwork(networks[i], 0));
			}
		}
		IO.pl("\t"+networkIdAnalyzedNetwork.size()+"\tNetworks loaded that pass minimum # genes ("+minimumNumberGenes+") and excluded types: "+typesToExclude);
	}
	
	private void loadKeggIdLookupHashes() throws IOException {
		HashMap<String, ArrayList<String>>[] hashes = KeggGeneSymbolIdExtractor.loadGeneLookupHashes(keggIdsFile);
		gs2ki = hashes[0];
		ki2gs = hashes[1];
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
			if (an.getVariantAdjPVal()<= maximumFdr) {
				//for each network in the an
				for (KeggApiNetwork kan: an.getVariantAnalizedNetworks()){
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
				if (an.getVariantPVal()< minPVal) minPVal = an.getVariantPVal();
				for (SelectGene sg: an.fetchVariantGenes()) allGenes.put(sg.getGeneSymbol(), sg);
				for (KeggApiNetwork net: an.getVariantAnalizedNetworks()) networkIds.add(net.getNetworkId());

				xls.append("\t"); xls.append(an.getVariantAnalizedNetworkIdNames()); xls.append("\n");
				xls.append("\t\tAdjPval: "); xls.append(Num.formatNumber(an.getVariantAdjPVal(),3)); xls.append("\n");
				xls.append("\t\tLog2Rto: "); xls.append(Num.formatNumber(an.getVariantLog2Rto(),3)); xls.append("\n");
				xls.append("\t\tGenes: "); xls.append(Misc.stringSetToString(an.getVarinatGeneNameHits(), ", ")); xls.append("\n");
			}
			SelectGene[] sg = new SelectGene[allGenes.size()];
			int i = 0;
			for (SelectGene g: allGenes.values()) sg[i++] = g;
			xls.append("AllGenes:\t"); xls.append(Misc.stringSetToString(allGenes.keySet(), ", ")); xls.append("\n");

			//networks
			String[] networkIdsStringArray = Misc.hashSetToStringArray(networkIds);
			String url = AnalyzedNetwork.fetchKeggPathwayMapLink(pathwayId, networkIdsStringArray, sg, gs2ki,AnalyzedNetwork.variantNegColor,AnalyzedNetwork.variantPosColor,AnalyzedNetwork.variantZeroColor);
			xls.append("KeggLink:\t"); xls.append(url); xls.append("\n");
			if (url.length()<250) {
				xls.append("ExcelLink:\t=HYPERLINK(\"");
				xls.append(url);		
				xls.append("\",\""+url);	
				xls.append("\")\n\n");
			}
			else xls.append("ExcelLink:\tToo big\n\n");
			results[counter++] = new StringValueSort(xls, minPVal);
		}
		Arrays.sort(results);
		PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "variantPathwaysMinGen"+minimumNumberGenes+"MaxFdr"+maximumFdr+".xls")));
		out.println("#Composite view of significant networks in KEGG pathways\n");
		for (StringValueSort s: results) out.print(s.getCargo().toString());
		out.close();
	}


	private void saveNetworks() {
		Arrays.sort(analyzedNetworks);
		try {
			PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "variantNetworksMinGen"+minimumNumberGenes+".xls")));
			out.println("# Network Name(s)\tNetwork Desc Link\tPval\tAdjPval\tAHits\tANoHits\tFracAHits\tAGeneHits\tBHits\tBNoHits\tFracBHits\tBGeneHits\tLog2(fracA/fracB)\tAllGeneHits\tPathwayMapLinksWithTopMapDescription...");
			for (AnalyzedNetwork an: analyzedNetworks) out.print(an.toStringVariant(addOne, gs2ki));
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void comparePathways() {
		FisherExact fe = null;
		float[] minLog10Pvals = new float[analyzedNetworks.length];
		for (int i=0; i< analyzedNetworks.length; i++) { 
			AnalyzedNetwork an = analyzedNetworks[i];
			int numGroupAHits = an.getGroupAHits().size();
			
			int numGroupANoHits = an.getGroupANoHits().size();
			
			int numGroupBHits = an.getGroupBHits().size();
			
			int numGroupBNoHits = an.getGroupBNoHits().size();
			
			if (fe == null) fe = new FisherExact(numGroupAHits+ numGroupANoHits+ numGroupBHits+ numGroupBNoHits);
			double pVal = fe.getTwoTailedP(numGroupAHits, numGroupANoHits, numGroupBHits, numGroupBNoHits);
			an.setVariantPVal(pVal);
			an.setSortValue(pVal);
			minLog10Pvals[i] = (float)Num.minus10log10(pVal);
		}
		
		//adjust the pvalues
		float[] adjPVals = Num.benjaminiHochbergCorrectUnsorted(minLog10Pvals);
		for (int i=0; i< analyzedNetworks.length; i++) {
			AnalyzedNetwork an = analyzedNetworks[i];
			double adjPval = Num.antiNeg10log10(adjPVals[i]);
			if (adjPval>1.0) adjPval = 1.0;
			an.setVariantAdjPVal(adjPval);
		}
	}
	
	//TODO: multithread, slow step
	private void loadGroup(File groupGeneHits, boolean isGroupA) {
		BufferedReader in;
		try {
			in = IO.fetchBufferedReader(groupGeneHits);
			String line;
			String[] cells;
			while ((line=in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("#") || line.length() ==0) continue;
				// 10DTM2THWO_IDTv1_A58476_A58475_SL430712	ABCC6	ADGRG1	APC	CREB1	DTD2	EFCAB1	EML6	FBXW7	KRAS	SLAIN1	TP53	TRPM3	ZFP36L2
				//                  0                         1        2     3    4       5        6      7    ....
				cells = Misc.TAB.split(line);
		
				//for each pathway network, count the number of gene hits
				for (AnalyzedNetwork p: networkIdAnalyzedNetwork.values()) {

					boolean hit = false;
					HashMap<String, KeggApiGene> networkGenes = p.getKeggApiNetwork().getGeneNameKeggApiGene();
					
					//for each gene in the sample
					for (int i=1; i< cells.length; i++) {
						//is it a known to kegg?
						if (lookupKeggGeneSymbol(cells[i]) == false) continue;
						
						if (networkGenes.containsKey(cells[i])) {
							if(isGroupA) {
								p.getGroupAHits().add(cells[0]);
								Integer geneCount = p.getaGenes().get(cells[i]);
								if (geneCount == null) p.getaGenes().put(cells[i], new Integer(1));
								else p.getaGenes().put(cells[i], new Integer(geneCount+1));
							}
							else {
								p.getGroupBHits().add(cells[0]);
								Integer geneCount = p.getbGenes().get(cells[i]);
								if (geneCount == null) p.getbGenes().put(cells[i], new Integer(1));
								else p.getbGenes().put(cells[i], new Integer(geneCount+1));
							}
							hit = true;
							break;
						}
					}
					if (hit == false) {
						if(isGroupA) p.getGroupANoHits().add(cells[0]);
						else p.getGroupBNoHits().add(cells[0]);
					}
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the gene hit file "+groupGeneHits);
		}
	}

	private boolean lookupKeggGeneSymbol(String geneSymbol) {
		if (gs2ki.containsKey(geneSymbol)) {
			Integer counter = foundGeneSymbols.get(geneSymbol);
			if (counter == null) foundGeneSymbols.put(geneSymbol, 1);
			else foundGeneSymbols.put(geneSymbol, counter++);
			return true;
		}
		else {
			Integer counter = missingGeneSymbols.get(geneSymbol);
			if (counter == null) missingGeneSymbols.put(geneSymbol, 1);
			else missingGeneSymbols.put(geneSymbol, counter++);
			return false;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KeggVariantPathwayAnalyzer(args);
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
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (resultsDirectory != null) resultsDirectory.mkdirs();

		if (typesToExclude != null) {
			String[] split = Misc.COMMA.split(typesToExclude);
			networkTypesToExclude = new HashSet<String>();
			for (String s: split) networkTypesToExclude.add(s.toLowerCase());
		}
	}
	
	private void checkFiles() throws IOException {
		IO.pl("Checking required files and directories...");
		File[] toCheck = {keggIdsFile, keggNetworkDirectory, resultsDirectory, groupAGeneHits, groupBGeneHits};
		String[] names = {"-k keggIdsLookupFile", "-n keggNetworkDir", "-r resultsDir", "-a groupAGeneHitsFile", "-b groupBGeneHitsFile"};
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
				"**                 Variant Pathway Comparator Kegg Medicus : June 2024              **\n" +
				"**************************************************************************************\n" +
				"For each KEGG MEDICUS pathway, VKPC creates a 2x2 contingency table and calculates a\n"+
				"Fisher's exact p-value that is subsequently multiple test corrected using Benjamini-\n"+
				"Hochberg's method. The contingency table is the number of subjects from cohort A with\n"+
				"one or more matching gene names, the number from A without any gene matches, and\n"+
				"likewise for the subjects in cohort B. A variety of statistics, including the\n"+
				"matching gene frequency, degree and directionality of change, and html links to each\n"+
				"network and pathway are saved in two spreadsheets. For TNRunner processed somatic\n"+
				"variant files, use the USeq AnnotatedVcfParser to select high impact, loss of\n"+
				"function/ CLINVAR patho/likely-pathogenic variants.\n"+
				
				"\nRequired Parameters:\n"+
				"-a File containing cohort A gene sets, each line represents a subject's genes of\n"+
				"     interest (e.g. those with HIGH impact mutations), tab delimited, the first cell\n"+
				"     is the subject ID, subsequent cells are the gene names.\n"+
				"-b File containing cohort B gene sets, ditto.\n"+
				"-k KEGG gene Id HUGO gene symbol lookup file generated by running the USeq\n"+
				"     KeggGeneSymbolIdExtractor.\n"+
				"-n KEGG network directory created by running the USeq KeggResourceExtractor.\n"+
				"-r Directory to save the spreadsheet results.\n"+
				"-t Network TYPEs to exclude from testing, comma delimited no spaces.\n"+
				"-m Minimum number of interrogated genes in a network for analysis, defaults to 5\n"+
				"-x Maximum FDR for including networks into the Kegg Pathway spreadsheet, defaults\n"+
				"     to 0.15\n"+
				"-o Add one to zero count A or B fractions when calculating the log2Rto(fracA/fracB)\n"+
				
				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VariantPathwayComparatorKegg -a \n"+
				"   earlyCRC.txt -b lateCRC.txt -k keggIdLookup.txt -n KeggNetworks/ -o\n"+
				"   -r VariantPathwayAnalyzerResults -t 'Pathogen,Ev factor' -x 0.1 \n"+

		"\n**************************************************************************************\n");

	}
}
