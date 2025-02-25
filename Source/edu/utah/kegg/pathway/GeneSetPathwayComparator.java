package edu.utah.kegg.pathway;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;


public class GeneSetPathwayComparator {

	//user fields
	private File interrogatedGeneList;
	private File selectGeneList;
	private File pathwayGenes;
	private File pathwayLinks;
	private File tempDirectory;
	private File fullPathToR;
	private File resultsSpreadsheet = null;
	private double minimumFdr = 0.1;

	//internal
	private HashSet<String> uniqueInterrogatedGenes = null;
	private HashSet<SelectGene> uniqueSelectGenes = null;
	private Pathway[] pathways = null;
	private HashMap<String, Pathway> networkIdPathway = new HashMap<String, Pathway>();
	
	
	//constructor
	public GeneSetPathwayComparator(String[] args) throws IOException{
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		checkFiles();
		
		IO.pl("\nLoading interrogated genes...");
		loadInterrogatedGenes();
		
		//parse pathways
		IO.pl("Loading KEGG Medicus pathways and links...");
		loadPathways();
		loadLinks();
		
		IO.pl("\nLoading select genes...");
		loadSelectedGenes();
		
		IO.pl("\nIntersect select genes with pathways...");
		intersectPathways();
		
		IO.pl("\nComparing pathway hits...");
		comparePathways();
		
		savePathways();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void loadSelectedGenes() throws IOException {
		String[] selectGenes = IO.loadFileIntoStringArray(selectGeneList);
		uniqueSelectGenes = new HashSet<SelectGene>();
		ArrayList<String> failedToFind = new ArrayList<String>();
		for (String s: selectGenes) {
			//GeneName Log2Rto
			String[] nameRto = Misc.TAB.split(s);
			if (nameRto.length!=2) throw new IOException ("Failed to parse the gene name and it's log2Rto for the select gene "+s);
			SelectGene sg = new SelectGene(nameRto);
			if (uniqueInterrogatedGenes.contains(sg.geneSymbol)) uniqueSelectGenes.add(sg);
			else failedToFind.add(sg.geneSymbol);
		}
		if (failedToFind.size()!=0) IO.pl("\tWARNING: failed to find "+failedToFind.size()+" in the interrogated gene list, skipping: "+Misc.stringArrayListToString(failedToFind, ","));
		IO.pl("\t"+selectGenes.length+" -> "+uniqueSelectGenes.size()+ " unique selected genes");
	}
	
	private class SelectGene implements Comparable<SelectGene>{
		String geneSymbol;
		String log2Rto;
		
		SelectGene(String[] geneRto){
			geneSymbol = geneRto[0];
			log2Rto = geneRto[1];
		}
		public int compareTo(SelectGene o) {
			return geneSymbol.compareTo(o.geneSymbol);
		}
	}
	
	private void loadInterrogatedGenes() {
		String[] allinterrogatedGenes = IO.loadFileIntoStringArray(interrogatedGeneList);
		uniqueInterrogatedGenes = new HashSet<String>();
		for (String s: allinterrogatedGenes) uniqueInterrogatedGenes.add(s);
		if (allinterrogatedGenes.length!= uniqueInterrogatedGenes.size()) IO.pl("\tWARNING: your interrogated gene list contains duplicates, these have been removed.");
		IO.pl("\t"+allinterrogatedGenes.length+" -> "+uniqueInterrogatedGenes.size()+ " unique genes\n");
	}

	private void savePathways() {
		try {
			PrintWriter out = new PrintWriter( new FileWriter(resultsSpreadsheet));
			//name pval, adjPval, #UniquePathwayGenes, #FoundUniquePathwayGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/PathGenes, IntersectingGenes
			out.println("#PathwayLink\tPval\tAdjPval\tPathGenes\tFoundPathGenes\tSelectGenes\tIntersect\tInt/Path\tIntGeneSymbols\tIntGenesLgRtos\tPathwayMapLinksWithTopMapDescription...");
			for (Pathway p: pathways) out.println(p.toString());
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void walkPathwayMaps() {
		
		//for each pathway(hsa05220) in the network(N00002), build a URL
		//	color the intersected genes pink for upregulated, skyblue for downregulated
		//  add selector for all of the significant networks (adjPVal <= 0.1)
		//e.g. https://www.kegg.jp/kegg-bin/show_pathway?map=hsa05220&multi_query=673%20skyblue%0A3265%20pink&network=N00002
		//  %20 between gene# and bkg color
		//  %0A between genes
		// pink:skyblue purple:orange
		/*Present the data by pathway map!
			For each map with a significant network
				Select the networks
				Aggregate and color the genes
				
		Rewrite using common objects from the Extractor?
		
		Here*/
		
			//find significant networks with a map
			HashSet<String> sigPathwayIds = new HashSet<String>();
			HashSet<String> sigNetworkIds = new HashSet<String>();
			for (Pathway p: pathways) {
				if (p.adjPval <= minimumFdr) {
					sigPathwayIds.add(p.name);
					sigNetworkIds.add(p.networkId);
				}
			}
			
			//for each sig network
		
	}

	private void comparePathways() throws IOException {
		int[][] nabt = new int[pathways.length][];
		
		//for each pathway calculate a hyperG pval for the intersection
		
		//n = total number of genes
		int n = uniqueInterrogatedGenes.size();
		//a = number of genes in set A
		int a = uniqueSelectGenes.size();
		
		for (int i=0; i< pathways.length; i++) { 
			Pathway p = pathways[i];
			
			//b = number of genes in set B
			int b = p.genes.size();
			
			//t = number of genes in common to A and B
			int t = p.sharedGenes.size();
			
			nabt[i] = new int[]{n,a,b,t};
			
			//if (i<30) IO.pl(p.name+"\t"+p.genes.size()+"\t"+Num.intArrayToString(nabt[i], "-")+"\t"+p.sharedGenes);
		}
		
		//Calculate pvalues using hg dist, fast way to do a fisher's exact on a 2x2 table
		IntersectListsHypergeometric ih = new IntersectListsHypergeometric(tempDirectory, fullPathToR);
		double[] pvals = ih.calculatePValues(nabt);
		for (int i=0; i< pathways.length; i++) pathways[i].pval = pvals[i];
		
		//convert the pvals to fdrs
		Arrays.sort(pathways); //smallest to largest
		double[] pvalsLargeToSmall = new double[pathways.length];
		int counter = 0;
		for (int i=pathways.length-1; i>=0; i--) pvalsLargeToSmall[counter++] = pathways[i].pval;
		Num.benjaminiHochbergCorrect(pvalsLargeToSmall);
		counter=0;
		for (int i=pathways.length-1; i>=0; i--) pathways[i].adjPval = pvalsLargeToSmall[counter++];
	}

	private void intersectPathways() {
		//for each pathway
		for (Pathway p: pathways) {
			//for each select gene
			Iterator<SelectGene> it = uniqueSelectGenes.iterator();
			while (it.hasNext()) {
				SelectGene sg = it.next();
				
				if (p.genes.contains(sg.geneSymbol)) p.sharedGenes.add(sg);
			}
		}
	}
	
	private void loadPathways() {
		BufferedReader in;
		ArrayList<Pathway> pathwaysAL = new ArrayList<Pathway>();
		try {
			in = IO.fetchBufferedReader(pathwayGenes);
			String line;
			String[] cells;
			while ((line=in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("#") || line.length() == 0) continue;
				cells = Misc.TAB.split(line);
				if (cells.length > 1) {
					Pathway p = new Pathway(cells);
					pathwaysAL.add(p);
					networkIdPathway.put(p.networkId, p);
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the pathway file "+pathwayGenes);
		}
		
		pathways = new Pathway[pathwaysAL.size()];
		pathwaysAL.toArray(pathways);
	}
	
	private class Pathway implements Comparable<Pathway>{
		String name;
		String networkId;
		int numInputGenes = -1;
		HashSet<String> genes;
		ArrayList<SelectGene> sharedGenes = new ArrayList<SelectGene>();
		double pval = -1;
		double adjPval = -1;
		ArrayList<String> failedToFind = new ArrayList<String>();
		ArrayList<String[]> urlLinkDescriptions = new ArrayList<String[]>();
		
		Pathway(String[] cells){
			name = cells[0];
			networkId = cells[1];
			IO.pl("\t"+name +"\t"+networkId);
			genes = new HashSet<String>(cells.length-2);
			for (int i=2; i< cells.length; i++) {
				if (uniqueInterrogatedGenes.contains(cells[i]))genes.add(cells[i]);
				else failedToFind.add(cells[i]);
				
			}
			
			if (failedToFind.size()!=0) IO.pl("\t\tWARNING: failed to find "+failedToFind.size()+" in the interrogated gene list, skipping: "+Misc.stringArrayListToString(failedToFind, ","));
			numInputGenes = genes.size()+failedToFind.size();
			IO.pl("\t\t"+ numInputGenes +" -> "+genes.size());
		}
		
		public String toString() {
			//name pval, adjPval, #UniquePathwayGenes, #FoundUniquePathwayGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/FoundUniPathGenes, IntersectingGenes, IntersectingGeneLogRtos, PathwayMapLinks...
			StringBuilder sb = new StringBuilder();
			
			sb.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			sb.append(networkId);
			sb.append("\",\"");
			sb.append(name);
			sb.append("\")");
			String networkCell = sb.toString();
			sb.append("\t");
			
			sb.append(pval); sb.append("\t");
			sb.append(Num.formatNumber(adjPval, 5)); sb.append("\t");
			sb.append(numInputGenes); sb.append("\t");
			sb.append(genes.size()); sb.append("\t");
			sb.append(uniqueSelectGenes.size()); sb.append("\t");
			sb.append(sharedGenes.size()); sb.append("\t");
			double frac = (double)sharedGenes.size()/ (double)genes.size();
			sb.append(Num.formatNumber(frac, 5)); sb.append("\t");
			
			//sort the sharedGenes by gene symbol so easier to group
			SelectGene[] sortedSGs = new SelectGene[sharedGenes.size()];
			sharedGenes.toArray(sortedSGs);
			Arrays.sort(sortedSGs);
			//just gene symbols
			for (SelectGene sg: sortedSGs) {
				sb.append(sg.geneSymbol);
				sb.append(" ");
			}
			sb.append("\t"); 
			
			//just log2Rtos have to use , to keep excel from converting it to a formula
			if (sortedSGs.length!=0) {
				sb.append(sortedSGs[0].log2Rto);
				for (int i=1; i< sortedSGs.length; i++) {
					sb.append(",");
					sb.append(sortedSGs[i].log2Rto);
				}
			}
			
			if (urlLinkDescriptions.size()!=0) {
				for (int i=0; i< urlLinkDescriptions.size(); i++) {
					sb.append("\t");
					String[] toLinkDesc = urlLinkDescriptions.get(i);
					sb.append("=HYPERLINK(\"https://www.kegg.jp/pathway/");
					sb.append(toLinkDesc[0]);
					sb.append("\",\"");
					sb.append(toLinkDesc[1]);
					sb.append("\")");
				}
			}
			else {
				sb.append("\t");
				sb.append(networkCell);
			}

			return sb.toString();
		}
		
		public String toStringNew() {
			//name pval, adjPval, #UniquePathwayGenes, #FoundUniquePathwayGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/FoundUniPathGenes, IntersectingGenes, IntersectingGeneLogRtos, PathwayMapLinks...
			StringBuilder sb = new StringBuilder();
			
			sb.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			sb.append(networkId);
			sb.append("\",\"");
			sb.append(name);
			sb.append("\")");
			String networkCell = sb.toString();
			sb.append("\t");
			
			sb.append(pval); sb.append("\t");
			sb.append(Num.formatNumber(adjPval, 5)); sb.append("\t");
			sb.append(numInputGenes); sb.append("\t");
			sb.append(genes.size()); sb.append("\t");
			sb.append(uniqueSelectGenes.size()); sb.append("\t");
			sb.append(sharedGenes.size()); sb.append("\t");
			double frac = (double)sharedGenes.size()/ (double)genes.size();
			sb.append(Num.formatNumber(frac, 5)); sb.append("\t");
			
			//sort the sharedGenes by gene symbol so easier to group
			SelectGene[] sortedSGs = new SelectGene[sharedGenes.size()];
			sharedGenes.toArray(sortedSGs);
			Arrays.sort(sortedSGs);
			//just gene symbols
			for (SelectGene sg: sortedSGs) {
				sb.append(sg.geneSymbol);
				sb.append(" ");
			}
			sb.append("\t"); 
			
			//just log2Rtos have to use , to keep excel from converting it to a formula
			if (sortedSGs.length!=0) {
				sb.append(sortedSGs[0].log2Rto);
				for (int i=1; i< sortedSGs.length; i++) {
					sb.append(",");
					sb.append(sortedSGs[i].log2Rto);
				}
			}
			
			if (urlLinkDescriptions.size()!=0) {
				for (int i=0; i< urlLinkDescriptions.size(); i++) {
					sb.append("\t");
					String[] toLinkDesc = urlLinkDescriptions.get(i);
					sb.append("=HYPERLINK(\"https://www.kegg.jp/pathway/");
					sb.append(toLinkDesc[0]);
					sb.append("\",\"");
					sb.append(toLinkDesc[1]);
					sb.append("\")");
				}
			}
			else {
				sb.append("\t");
				sb.append(networkCell);
			}

			return sb.toString();
		}


		public int compareTo(Pathway o) {
			if (this.pval< o.pval) return -1;
			if (this.pval> o.pval) return 1;
			return 0;
		}
		
	}
	
	private void loadLinks() {
		BufferedReader in;
		try {
			in = IO.fetchBufferedReader(pathwayLinks);
			String line;
			String[] cells;
			while ((line=in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("#") || line.length() == 0 || line.contains("NADA")) continue;
				//N01486	hsa04110+N01486	 Cell cycle
				cells = Misc.TAB.split(line);
				if (cells.length > 1) {
					Pathway kp = networkIdPathway.get(cells[0]);
					if (kp == null) throw new IOException("Missing Kegg Network ID in the gene pathways file "+cells[0]);
					String[] pd = new String[] {cells[1], cells[2]};
					kp.urlLinkDescriptions.add(pd);
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the pathway link file "+pathwayLinks);
		}
		
	}


	public static void main(String[] args) throws IOException {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new GeneSetPathwayComparator(args);
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
					case 'g': selectGeneList = new File(args[++i]); break;
					case 't': tempDirectory = new File(args[++i]); break;
					case 'e': fullPathToR = new File(args[++i]); break;
					
					case 'p': pathwayGenes = new File(args[++i]); break;
					case 'l': pathwayLinks = new File(args[++i]); break;
					case 'r': resultsSpreadsheet = new File(args[++i]); break;
					
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (tempDirectory != null) tempDirectory.mkdirs();

	}	
	
	private void checkFiles() {
		
		IO.pl("Checking required files...");
		File[] toCheck = {interrogatedGeneList, selectGeneList, tempDirectory, fullPathToR, pathwayGenes, pathwayLinks};
		String[] names = {"-i interrogatedGeneList", "-g selectGeneList", "-t tempDirectory", "-e fullPathToR", "-p pathwayGenes", "-l pathwayLinks"};
		boolean notFound = false;
		for (int i=0; i< toCheck.length; i++) {
			if (toCheck[i] == null || toCheck[i].exists()== false) {
				IO.el("\tMissing "+names[i]+" : "+toCheck[i]);
				notFound = true;
			}
		}
		if (resultsSpreadsheet == null) {
			IO.el("\tMissing -r results spreadsheet output file?");
			notFound = true;
		}
		if (notFound) Misc.printErrAndExit("Correct issues and restart.\n");
		
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Gene Set Pathway Comparator : June 2023                    **\n" +
				"**************************************************************************************\n" +
				"GSPC uses the interrogatedGeneList to filter the selectGeneList and pathwayGenes sets,\n"+
				"intersects the filtered selectGeneList against each filtered pathwayGene sets and\n"+
				"calculates a p-value for the degree of intersection using a hypergeometric\n"+
				"distribution. These are subsequently multiple test corrected using Benjamini-\n"+
				"Hochberg FDR method. Use the same KEGG MEDICUS reference files as used by the USeq\n"+
				"VariantPathwayComparatorKeggMedicus app then run the USeq JointPathwayComparator to\n"+
				"calculate combine p-values for each.\n"+
				
				"\nRequired Parameters:\n"+
				"-i File containing all of the interrogated genes in the study, one gene per line.\n"+
				"-g File containing the selected genes of interest with their log2Rto, one per line,\n"+
				"   tab delimited.\n"+
				"-p File containing KEGG MEDICUS pathways to compare, each line represents one pathway,\n"+
				"     tab delimited: unique pathway name, network ID, and associated gene symbols (e.g.\n"+
		        "     REFERENCE_ATR_SIGNALING  N01451 ATR ATRIP CHEK1 HUS1 ...)\n"+
		        "-l  File containing KEGG network IDs, pathway links, and their descriptions (e.g.\n"+
		        "    N01486 hsa04110+N01486 CellCycle), one network per row, tab delimited.\n"+
				"-r File to save the txt results, should end with .xls\n"+
				"-e File path to the R executable.\n"+
				"-t Directory in which to save temporary files.\n"+
				
				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/GeneSetPathwayComparator -i allGenes.txt\n"+
				"   -g diffExpGenes.txt -l keggLinks.txt -p keggPathways.txt -r earlyVsLateRnaSeq.xls\n"+
				"   -t Temp -e /usr/bin/R \n"+

		"\n**************************************************************************************\n");

	}
}
