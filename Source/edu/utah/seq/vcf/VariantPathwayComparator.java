package edu.utah.seq.vcf;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;


public class VariantPathwayComparator {

	//user fields
	private File groupAGeneHits;
	private File groupBGeneHits;
	private File pathwayGenes;
	private File geneKeggIdName;
	
	//internal
	private ArrayList<Pathway> pathways = new ArrayList<Pathway>();
	private File resultsSpreadsheet = null;
	private HashMap<String,String> keggGeneNameID = null;
	private double minKeggFreq = 0.025;
	
	
	//constructor
	public VariantPathwayComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//parse pathways
		IO.pl("Loading pathways...");
		loadPathways();
		
		IO.pl("\nLoading group A gene hits...");
		loadGroup(groupAGeneHits, true);
		
		IO.pl("\nLoading group B gene hits...");
		loadGroup(groupBGeneHits, false);
		
		IO.pl("\nComparing pathway hits...");
		comparePathways();
		
		sortSavePathways();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void sortSavePathways() {
		Pathway[] loadedPathways = new Pathway[pathways.size()];
		pathways.toArray(loadedPathways);
		Arrays.sort(loadedPathways);
		
		try {
			PrintWriter out = new PrintWriter( new FileWriter(this.resultsSpreadsheet));
			if (geneKeggIdName == null) out.println("NameHyperLink\tPval\tAdjPval\tAHits\tANoHits\tFracAHits\tAGeneHits\tBHits\tBNoHits\tFracBHits\tBGeneHits\tLog2(fracA/fracB)");
			else out.println("NameHyperLink\tPval\tAdjPval\tAHits\tANoHits\tFracAHits\tAGeneHits\tAGeneLinks\tBHits\tBNoHits\tFracBHits\tBGeneHits\tBGeneLinks\tLog2(fracA/fracB)");
			for (Pathway p: loadedPathways) out.println(p.toString(keggGeneNameID));
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	private void comparePathways() {
		FisherExact fe = null;
		int numPathways = pathways.size();
		float[] minLog10Pvals = new float[numPathways];
		for (int i=0; i< numPathways; i++) { 
			Pathway p = pathways.get(i);
//IO.pl(p.name);
//IO.pl("\t"+p.genes.size()+"\t"+Misc.hashSetToString(p.genes, ","));
			int numGroupAHits = p.groupAHits.size();
//IO.pl("\t"+numGroupAHits +"\tGroupAHits\t"+Misc.stringArrayListToString(p.groupAHits, ","));
			
			int numGroupANoHits = p.groupANoHits.size();
//IO.pl("\t"+numGroupANoHits +"\tGroupANoHits");
			
			int numGroupBHits = p.groupBHits.size();
//IO.pl("\t"+numGroupBHits +"\tGroupBHits\t"+Misc.stringArrayListToString(p.groupBHits, ","));
			
			int numGroupBNoHits = p.groupBNoHits.size();
//IO.pl("\t"+numGroupBNoHits +"\tGroupBNoHits");
			
			if (fe == null) fe = new FisherExact(numGroupAHits+ numGroupANoHits+ numGroupBHits+ numGroupBNoHits);
			double pVal = fe.getTwoTailedP(numGroupAHits, numGroupANoHits, numGroupBHits, numGroupBNoHits);
//IO.pl("\t"+pVal+"\t"+p.name);
			p.pval = pVal;
			minLog10Pvals[i] = (float)Num.minus10log10(pVal);
		}
		
		//adjust the pvalues
		float[] adjPVals = Num.benjaminiHochbergCorrectUnsorted(minLog10Pvals);
		for (int i=0; i< numPathways; i++) {
			Pathway p = pathways.get(i);
			p.adjPval = Num.antiNeg10log10(adjPVals[i]);
			if (p.adjPval>1) p.adjPval = 1.0;
		}
	}

	private void loadGroup(File groupGeneHits, boolean isGroupA) {
		BufferedReader in;
		try {
			in = IO.fetchBufferedReader(groupGeneHits);
			String line;
			String[] cells;
			while ((line=in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("#") || line.length() ==0) continue;
				cells = Misc.TAB.split(line);
//IO.pl("\nSample "+line);
				//for each pathway
				for (Pathway p: pathways) {
//IO.pl("Pathway "+p.name+" -> "+p.genes);
					boolean hit = false;
					//for each gene in the sample
					for (int i=1; i< cells.length; i++) {
//IO.pl("\tSampleGene\t"+cells[i]);
						if (p.genes.contains(cells[i])) {
							if(isGroupA) {
								p.groupAHits.add(cells[0]);
								Integer geneCount = p.aGenes.get(cells[i]);
								if (geneCount == null) p.aGenes.put(cells[i], new Integer(1));
								else p.aGenes.put(cells[i], new Integer(geneCount+1));
							}
							else {
								p.groupBHits.add(cells[0]);
								Integer geneCount = p.bGenes.get(cells[i]);
								if (geneCount == null) p.bGenes.put(cells[i], new Integer(1));
								else p.bGenes.put(cells[i], new Integer(geneCount+1));
							}
							hit = true;
//IO.pl("\t\t\tHIT");
							break;
						}
					}
					if (hit == false) {
//IO.pl("\t\tNoHit");
						if(isGroupA) p.groupANoHits.add(cells[0]);
						else p.groupBNoHits.add(cells[0]);
					}
					else {
//IO.pl("\t\tHit");
//IO.pl("\nSample "+line);
//System.exit(0);
					}
				}
				
				
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the gene hit file "+groupGeneHits);
		}
		
	}

	private void loadPathways() {
		BufferedReader in;
		try {
			in = IO.fetchBufferedReader(pathwayGenes);
			String line;
			String[] cells;
			while ((line=in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("#") || line.length() == 0) continue;
				cells = Misc.TAB.split(line);
				if (cells.length > 1) pathways.add(new Pathway(cells));
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the pathway file "+pathwayGenes);
		}
		
		if (geneKeggIdName!=null) keggGeneNameID = IO.loadFileIntoHash(geneKeggIdName, 1, 0);
	}
	
	private class Pathway implements Comparable<Pathway>{
		String name;
		HashSet<String> genes;
		ArrayList<String> groupAHits = new ArrayList<String>();
		ArrayList<String> groupANoHits = new ArrayList<String>();
		ArrayList<String> groupBHits = new ArrayList<String>();
		ArrayList<String> groupBNoHits = new ArrayList<String>();
		double pval = -1;
		double adjPval = -1;
		HashMap<String, Integer> aGenes = new HashMap<String, Integer>();
		HashMap<String, Integer> bGenes = new HashMap<String, Integer>();
		
		Pathway(String[] cells){
			name = cells[0];
			genes = new HashSet<String>(cells.length-1);
			for (int i=1; i< cells.length; i++) genes.add(cells[i]);
		}
		
		private final Pattern hsaPath = Pattern.compile("(^hsa\\d+)\\s.+");
		
		
		
		public String toString(HashMap<String,String>keggGeneNameId) {
			//name pval, adjPval, AHits, ANoHits, fracAHits, AGeneHits, (ALinks,) BHits, BNoHits, fracBHits, BGeneHits, (BLinks), log2(fracA/fracB)
			StringBuilder sb = new StringBuilder();
			//attempt to parse the hsa0000 pathway number
			Matcher mat = hsaPath.matcher(name);
			String partialHyperLink = null;
			if (mat.matches()) {
				sb.append("=HYPERLINK(\"https://www.kegg.jp/pathway/");
				sb.append(mat.group(1));
				partialHyperLink = sb.toString();
				sb.append("\",\"");
				sb.append(name);
				sb.append("\")");
			}
			else sb.append(name); sb.append("\t");
			sb.append(Num.formatNumber(pval, 5)); sb.append("\t");
			sb.append(Num.formatNumber(adjPval, 5)); sb.append("\t");
			sb.append(groupAHits.size()); sb.append("\t");
			sb.append(groupANoHits.size()); sb.append("\t");
			double aFrac = (double)groupAHits.size()/ (double)(groupAHits.size()+groupANoHits.size());
			sb.append(Num.formatNumber(aFrac, 5)); sb.append("\t");
			sb.append(convertToSortedFrequencies(aGenes, partialHyperLink)); sb.append("\t");
			
			sb.append(groupBHits.size()); sb.append("\t");
			sb.append(groupBNoHits.size()); sb.append("\t");
			double bFrac = (double)groupBHits.size()/ (double)(groupBHits.size()+groupBNoHits.size());
			sb.append(Num.formatNumber(bFrac, 5)); sb.append("\t");
			sb.append(convertToSortedFrequencies(bGenes, partialHyperLink)); sb.append("\t");

			double rto = aFrac/bFrac;
			sb.append(Num.formatNumber(Num.log2(rto), 5));
			return sb.toString();
		}

		public int compareTo(Pathway o) {
			if (this.pval< o.pval) return -1;
			if (this.pval> o.pval) return 1;
			return 0;
		}
		
	}
	
	public String convertToSortedFrequencies(HashMap<String,Integer> keyCounts, String pathwayName) {
		double totalCounts = 0;
		for (String key: keyCounts.keySet()) totalCounts+= keyCounts.get(key);
		GeneCount[] gc = new GeneCount[keyCounts.size()];
		int count = 0;
		for (String key: keyCounts.keySet()) {
			double geneCount = keyCounts.get(key);
			if (totalCounts ==0) geneCount = 0;
			else geneCount = geneCount/totalCounts;
			gc[count++] = new GeneCount(key, geneCount);
		}
		Arrays.sort(gc);
		
		StringBuilder namesToView = new StringBuilder();
		ArrayList<String> namesToLink = new ArrayList<String>();
		int last = gc.length-1;
		for (int i=0; i<gc.length; i++) {
			namesToView.append(gc[i].name);
			namesToView.append("=");
			namesToView.append(Num.formatNumber(gc[i].fraction, 3));
			if (i!=last)namesToView.append(",");
			//any for linking?
			if (gc[i].fraction >= minKeggFreq) namesToLink.add(gc[i].name);
		}
		String geneNamesToView = namesToView.toString();
		if (keggGeneNameID != null) return geneNamesToView+ "\t"+ createGeneHyperLink(pathwayName, namesToLink, geneNamesToView);
		return geneNamesToView;
	}
	private String createGeneHyperLink(String pathwayLink, ArrayList<String> namesToLink, String geneNamesToView) {
		//https://www.kegg.jp/pathway/hsa03440
		//must watch size, max 250
		StringBuilder sb = new StringBuilder(pathwayLink);
		for (String geneName: namesToLink) {
			String id = keggGeneNameID.get(geneName);
			int size = sb.length() + id.length()+1;
			if (size > 250) break;
			if (id != null) {
				sb.append("+");
				sb.append(id);
			
			}
		}
		sb.append("\",\"");
		if (geneNamesToView.length()> 250) sb.append(geneNamesToView.subSequence(0, 247)+"...");
		else sb.append(geneNamesToView);
		sb.append("\")");
		return sb.toString();
	}

	public static class GeneCount implements Comparable<GeneCount>{
		double fraction = 0;
		String name = null;
		
		public GeneCount(String name, double fraction) {
			this.name = name;
			this.fraction = fraction;
		}
		public int compareTo(GeneCount o) {
			if (this.fraction> o.fraction) return -1;
			if (this.fraction< o.fraction) return 1;
			return this.name.compareTo(o.name);
		}
		
	}
	


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VariantPathwayComparator(args);
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
					case 'p': pathwayGenes = new File(args[++i]); break;
					case 'k': geneKeggIdName = new File(args[++i]); break;
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
		
		if (groupAGeneHits == null || groupBGeneHits == null || pathwayGenes == null || resultsSpreadsheet == null) {
			Misc.printErrAndExit("\nERROR: failed to find one or more of the required input files. See the help menu.\n");
		}
	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Variant Pathway Comparator : Feb 2022                  **\n" +
				"**************************************************************************************\n" +
				"For each pathway, VPC creates a 2x2 contingency table and calculates a Fisher's exact\n"+
				"p-value that is subsequently multiple test corrected using Benjamini-Hochberg's\n"+
				"method. The contingency table is the number of subjects from cohort A with one or more\n"+
				"matching gene names, the number from A without any gene matches, and likewise for the\n"+
				"subjects in cohort B. A variety of statistics, including the matching gene frequency,\n"+
				"degree and directionality of change, and html links to each pathway are saved.\n"+
				
				"\nRequired Parameters:\n"+
				"-a File containing cohort A gene sets, each line represents a subject's genes of\n"+
				"     interest (e.g. those with HIGH impact mutations), tab delimited, the first cell\n"+
				"     is the subject ID, subsequent cells are the gene names.\n"+
				"-b File containing cohort B gene sets, ditto.\n"+
				"-p File containing pathways to compare, each line represents one pathway, tab\n"+
				"     delimited, the first cell is the pathway ID and description (e.g.\n"+
		        "     'hsa05210  Colorectal cancer'), subsequent cells, the associated genes.\n"+
				"-r File to save results, should end with .xls\n"+
		        "-k (Optional) Tab delimited file containing KEGG gene IDs and gene Names, one row\n"+
				"      per gene, for highlighting impacted genes in each KEGG pathway hyperlink.\n"+
				
				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VariantPathwayComparator -a earlyCRC.txt\n"+
				"   -b lateCRC.txt -p keggHumanPathways.txt -r earlyVsLateVPC.xls -k keggGeneIDName.txt\n"+

		"\n**************************************************************************************\n");

	}
}
