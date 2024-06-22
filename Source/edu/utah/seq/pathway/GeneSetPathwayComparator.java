package edu.utah.seq.pathway;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;


public class GeneSetPathwayComparator {

	//user fields
	private File interrogatedGeneList;
	private File selectGeneList;
	private File pathwayGenes;
	private File tempDirectory;
	private File fullPathToR;
	private File resultsSpreadsheet = null;

	//internal
	private HashSet<String> uniqueInterrogatedGenes = null;
	private HashSet<String> uniqueSelectGenes = null;
	private Pathway[] pathways = null;
	
	
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
		IO.pl("Loading pathways...");
		loadPathways();
		
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
	
	private void loadSelectedGenes() {
		String[] selectGenes = IO.loadFileIntoStringArray(selectGeneList);
		uniqueSelectGenes = new HashSet<String>();
		ArrayList<String> failedToFind = new ArrayList<String>();
		for (String s: selectGenes) {
			if (uniqueInterrogatedGenes.contains(s)) uniqueSelectGenes.add(s);
			else failedToFind.add(s);
		}
		if (failedToFind.size()!=0) IO.pl("\tWARNING: failed to find "+failedToFind.size()+" in the interrogated gene list, skipping: "+Misc.stringArrayListToString(failedToFind, ","));
		IO.pl("\t"+selectGenes.length+" -> "+uniqueSelectGenes.size()+ " unique selected genes");
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
			out.println("#Name\tPval\tAdjPval\tPathGenes\tFoundPathGenes\tSelectGenes\tIntersect\tInt/Path\tIntGenes");
			for (Pathway p: pathways) out.println(p.toString());
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
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
			Iterator<String> it = uniqueSelectGenes.iterator();
			while (it.hasNext()) {
				String g = it.next();
				if (p.genes.contains(g)) p.sharedGenes.add(g);
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
				if (cells.length > 1) pathwaysAL.add(new Pathway(cells));
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
		int numInputGenes = -1;
		HashSet<String> genes;
		ArrayList<String> sharedGenes = new ArrayList<String>();
		double pval = -1;
		double adjPval = -1;
		ArrayList<String> failedToFind = new ArrayList<String>();
		
		Pathway(String[] cells){
			name = cells[0];
			IO.pl("\t"+name);
			genes = new HashSet<String>(cells.length-1);
			for (int i=1; i< cells.length; i++) {
				if (uniqueInterrogatedGenes.contains(cells[i]))genes.add(cells[i]);
				else failedToFind.add(cells[i]);
				
			}
			if (failedToFind.size()!=0) IO.pl("\t\tWARNING: failed to find "+failedToFind.size()+" in the interrogated gene list, skipping: "+Misc.stringArrayListToString(failedToFind, ","));
			numInputGenes = genes.size()+failedToFind.size();
			IO.pl("\t\t"+ numInputGenes +" -> "+genes.size());
		}
		
		private final Pattern hsaPath = Pattern.compile("(^hsa\\d+)\\s.+");
		
		public String toString() {
			//name pval, adjPval, #UniquePathwayGenes, #FoundUniquePathwayGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/FoundUniPathGenes, IntersectingGenes
			StringBuilder sb = new StringBuilder();
			//attempt to parse the hsa0000 pathway number
			Matcher mat = hsaPath.matcher(name);
			if (mat.matches()) {
				sb.append("=HYPERLINK(\"https://www.kegg.jp/pathway/");
				sb.append(mat.group(1));
				sb.append("\",\"");
				sb.append(name);
				sb.append("\")");
			}
			else sb.append(name); sb.append("\t");
			sb.append(pval); sb.append("\t");
			sb.append(Num.formatNumber(adjPval, 5)); sb.append("\t");
			sb.append(numInputGenes); sb.append("\t");
			sb.append(genes.size()); sb.append("\t");
			sb.append(uniqueSelectGenes.size()); sb.append("\t");
			sb.append(sharedGenes.size()); sb.append("\t");
			double frac = (double)sharedGenes.size()/ (double)genes.size();
			sb.append(Num.formatNumber(frac, 5)); sb.append("\t");
			if (sharedGenes.size()!=0) sb.append(Misc.stringArrayListToString(sharedGenes, " ")); 

			return sb.toString();
		}

		public int compareTo(Pathway o) {
			if (this.pval< o.pval) return -1;
			if (this.pval> o.pval) return 1;
			return 0;
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
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'p': pathwayGenes = new File(args[++i]); break;
					case 's': resultsSpreadsheet = new File(args[++i]); break;
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
		File[] toCheck = {interrogatedGeneList, selectGeneList, tempDirectory, fullPathToR, pathwayGenes, resultsSpreadsheet};
		String[] names = {"-i interrogatedGeneList", "-g selectGeneList", "-t tempDirectory", "-r fullPathToR", "-p pathwayGenes", "-s resultsSpreadsheet"};
		boolean notFound = false;
		for (int i=0; i< toCheck.length; i++) {
			if (toCheck[i] == null || toCheck[i].exists()== false) {
				IO.el("\tMissing "+names[i]+" : "+toCheck[i]);
				notFound = true;
			}
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
				"Hochberg FDR method.\n"+
				
				"\nRequired Parameters:\n"+
				"-i File containing all of the interrogated genes in the study, one gene symbol per\n"+
				"     line.\n"+
				"-g File containing the selected genes of interest, ditto.\n"+
				"-p File containing pathways to compare, each line represents one pathway, tab\n"+
				"     delimited, the first cell is the pathway ID and description (e.g.\n"+
		        "     'hsa05210  Colorectal cancer'), subsequent cells, the associated genes.\n"+
				"-s File to save the results spreadsheet, should end with .txt\n"+
				"-r File path to the R application.\n"+
				"-t Directory in which to save temporary files.\n"+
				
				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/GeneSetPathwayComparator -i \n"+
				"   allGenes.txt -g diffExpGenes.txt -p keggPathways.txt -s earlyVsLate.xls -t Temp\n"+
				"   -r /usr/bin/R \n"+

		"\n**************************************************************************************\n");

	}
}
