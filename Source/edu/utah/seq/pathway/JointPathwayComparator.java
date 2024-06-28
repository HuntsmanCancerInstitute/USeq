package edu.utah.seq.pathway;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;

public class JointPathwayComparator {

	//user fields
	private File variantResults;
	private File geneResults;
	private File resultsSpreadsheet = null;

	//internal
	private HashMap<String, Pathway> pathways = null;
	
	//constructor
	public JointPathwayComparator(String[] args) throws IOException{
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		IO.pl("\nLoading gene set pathway comparator results...");
		loadGeneResults();
		
		IO.pl("Loading variant pathway comparator results...");
		loadVariantResults();
		
		IO.pl("Combining PValues...");
		combinePValues();
		
		savePathways();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	

	


	private void savePathways() {
		try {
			PrintWriter out = new PrintWriter( new FileWriter(resultsSpreadsheet));
			String header = "#Pathway Page Link\tDiff Exp PVal\tDiff Exp FDR\tDiff Exp Gene Symb\tDiff Exp Gene LgRtos\t"+
              "Som Var PVal\tSom Var FDR\tSom Var Genes\tSom Var LgRto\t"+
              "Combine PVal\tCombine FDR\tCombine Genes\tPathway Map Links...";			
			out.println(header);
			for (Pathway p: pathways.values()) out.println(p.toString());
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void combinePValues() throws IOException {
		CombinePValues cp = new CombinePValues();
		double[] toCombine = new double[2];
		Pathway[] pathwayArray = new Pathway[pathways.size()];
		int counter = 0;
		for (Pathway p: pathways.values()) {
			toCombine[0] = p.diffExpPVal;
			toCombine[1] = p.variantPVal; 
			//watch out with 1's can cause non convergence
			if (p.diffExpPVal ==1.0 || p.variantPVal==1.0) p.combinePVal = 1.0;
			else p.combinePVal = cp.calculateCombinePValues(toCombine);
			pathwayArray[counter++] = p;
		}
		//convert the pvals to fdrs
		Arrays.sort(pathwayArray); //largest to smallest
		
		double[] pvalsLargeToSmall = new double[pathwayArray.length];
		for (int i=0; i< pathwayArray.length; i++) pvalsLargeToSmall[i] = pathwayArray[i].combinePVal;
		
		Num.benjaminiHochbergCorrect(pvalsLargeToSmall);
		for (int i=0; i< pathwayArray.length; i++) pathwayArray[i].combineAdjPval = pvalsLargeToSmall[i];
	}
	
	private void loadVariantResults() {
		BufferedReader in;
		ArrayList<String> missingPathways = new ArrayList<String>();
		try {
			in = IO.fetchBufferedReader(variantResults);
			String line;
			String[] cells;
			while ((line=in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("#") || line.length() == 0) continue;
				//#NameHyperLink\tPval
				cells = Misc.TAB.split(line);
				Pathway p = pathways.get(cells[0]);
				if (p==null) missingPathways.add(cells[0]);
				else p.addSomVarInfo(cells);;
			}
			in.close();
			
			//checkem
			if (missingPathways.size()!=0) Misc.printErrAndExit("\tMissing the following pathways from the gene results: "+missingPathways);
			missingPathways.clear();
			for (Pathway p : pathways.values()) {
				if (p.variantPVal == -1.0) missingPathways.add(p.pathwayPageLink);
			}
			if (missingPathways.size()!=0) Misc.printErrAndExit("\tMissing the following pathways from the variant results: "+missingPathways);
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the variant set pathway results file "+geneResults);
		}
	}
	
	private void loadGeneResults() {
		BufferedReader in;
		pathways = new HashMap<String, Pathway>();
		try {
			in = IO.fetchBufferedReader(geneResults);
			String line;
			String[] cells;
			while ((line=in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("#") || line.startsWith("Name") || line.length() == 0) continue;
				cells = Misc.TAB.split(line);
				pathways.put(cells[0], new Pathway(cells));
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the gene set pathway results file "+geneResults);
		}
	}
	
	private class Pathway implements Comparable<Pathway>{
		String pathwayPageLink;
		ArrayList<String> subPathwayLinks = new ArrayList<String>();
		
		//diffExp
		double diffExpPVal = -1.0;
		double diffExpAdjPVal = -1.0;
		String diffExpGenes = null;
		String diffExpGeneLgRtos = null;
		
		//somVar
		double variantPVal = -1.0;
		double variantAdjPVal = -1.0;
		String variantLog2Rto = null;
		String variantGenes = null;
		
		
		//combine
		double combinePVal = -1.0;
		double combineAdjPval = -1.0;
		String combineGenes = null;
		
		Pathway(String[] cells){
			//diffExp
			//	      0        1        2           8               9                        10...
			//#PathwayLink	Pval	AdjPval	IntGeneSymbols	IntGenesLgRtos	PathwayMapLinksWithTopMapDescription...
			pathwayPageLink = cells[0];
			diffExpPVal = Double.parseDouble(cells[1]);
			diffExpAdjPVal = Double.parseDouble(cells[2]);
			diffExpGenes = cells[8];
			diffExpGeneLgRtos = cells[9];
			for (int i=10; i< cells.length; i++) subPathwayLinks.add(cells[i]);
			
		}
		
		
		
		public void addSomVarInfo (String[] cells) {
			//	      0      1     2           11             12                     13
			//#PathwayLink Pval AdjPval Log2(fracA/fracB) AllGeneHits PathwayMapLinksWithTopMapDescription...
			variantPVal = Double.parseDouble(cells[1]);
			variantAdjPVal = Double.parseDouble(cells[2]);
			variantLog2Rto = cells[11];
			variantGenes = cells[12];
			
			//make merged gene symbols
			TreeSet<String> allGenes = new TreeSet<String>();
			if (diffExpGenes.trim().length()!=0) {
				for (String s: Misc.WHITESPACE.split(diffExpGenes)) allGenes.add(s);
			}
			if (variantGenes.trim().length()!=0) {
				for (String s: Misc.WHITESPACE.split(variantGenes)) allGenes.add(s);
			}
			combineGenes = Misc.treeSetToString(allGenes, " ");

		}
		
		public int compareTo(Pathway o) {
			if (this.combinePVal > o.combinePVal) return -1;
			if (this.combinePVal < o.combinePVal) return 1;
			return 0;
		}
		
		public String toString() {
            
			StringBuilder sb = new StringBuilder();
			sb.append(pathwayPageLink); sb.append("\t");
			
			sb.append(diffExpPVal); sb.append("\t");
			sb.append(diffExpAdjPVal); sb.append("\t");
			sb.append(diffExpGenes); sb.append("\t");
			sb.append(diffExpGeneLgRtos); sb.append("\t");;
			
			sb.append(variantPVal); sb.append("\t");
			sb.append(variantAdjPVal); sb.append("\t");
			sb.append(variantGenes); sb.append("\t");
			sb.append(variantLog2Rto); sb.append("\t");
			
			
			sb.append(combinePVal); sb.append("\t");
			sb.append(combineAdjPval);sb.append("\t");
			sb.append(combineGenes);
			
			for (String ml: this.subPathwayLinks) {
				sb.append("\t");
				sb.append(ml);
			}
			
			return sb.toString();
		}
	}
	
	public static void main(String[] args) throws IOException {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new JointPathwayComparator(args);
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
					case 'v': variantResults = new File(args[++i]); break;
					case 'g': geneResults = new File(args[++i]); break;
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
		checkFiles();
	}
	
	private void checkFiles() {
		IO.pl("Checking required files...");
		File[] toCheck = {variantResults, geneResults, resultsSpreadsheet};
		boolean[] exists = {true, true, false};
		String[] names = {"-v variantResults", "-g geneResults", "-s spreadsheetOut"};
		boolean notFound = false;
		for (int i=0; i< toCheck.length; i++) {
			if (toCheck[i] == null || (exists[i] && toCheck[i].exists()== false)) {
				IO.el("\tMissing "+names[i]+" : "+toCheck[i]);
				notFound = true;
			}
		}
		if (notFound) Misc.printErrAndExit("Correct issues and restart.\n");
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Joint Pathway Comparator : June 2024                      **\n" +
				"**************************************************************************************\n" +
				"JPC parses the output of the GeneSet and Variant Pathway Comparators when run on the\n"+
				"same pathway set, combines the pvalues using Fisher's method, and multiple test\n"+
				"corrects these using the Benjamini-Hochberg FDR method.\n"+
				
				"\nRequired Parameters:\n"+
				"-g File containing the GeneSetPathwayComparator txt results.\n"+
				"-v File containing the VariantPathwayComparator txt results.\n"+
				"-s File to save the spreadsheet results, should end with .txt or .xls\n"+

				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/JointPathwayComparator -s \n"+
				"   joint.kegg.txt -g gene.kegg.txt -v var.kegg.txt \n"+

		"\n**************************************************************************************\n");

	}
}
