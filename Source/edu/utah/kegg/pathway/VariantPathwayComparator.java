package edu.utah.kegg.pathway;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;


public class VariantPathwayComparator {

	//user fields
	private File groupAGeneHits;
	private File groupBGeneHits;
	private File pathwayGenes;
	private File pathwayLinks;
	private boolean addOne = false;
	
	//internal
	private ArrayList<KeggPathway> pathways = new ArrayList<KeggPathway>();
	private HashMap<String, KeggPathway> networkIdPathway = new HashMap<String, KeggPathway>();
	private File resultsSpreadsheet = null;
	private double minKeggFreq = 0.0005;
	
	
	//constructor
	public VariantPathwayComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//parse pathways
		IO.pl("Loading KEGG Medicus pathways and links...");
		loadPathways();
		loadLinks();
		
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
		KeggPathway[] loadedPathways = new KeggPathway[pathways.size()];
		pathways.toArray(loadedPathways);
		Arrays.sort(loadedPathways);
		
		try {
			PrintWriter out = new PrintWriter( new FileWriter(this.resultsSpreadsheet));
			out.println("#PathwayLink\tPval\tAdjPval\tAHits\tANoHits\tFracAHits\tAGeneHits\tBHits\tBNoHits\tFracBHits\tBGeneHits\tLog2(fracA/fracB)\tAllGeneHits\tPathwayMapLinksWithTopMapDescription...");

			for (KeggPathway p: loadedPathways) out.println(p.toString(addOne));
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
			KeggPathway p = pathways.get(i);
//IO.pl(p.name);
//IO.pl("\t"+p.genes.size()+"\t"+Misc.hashSetToString(p.genes, ","));
			int numGroupAHits = p.getGroupAHits().size();
//IO.pl("\t"+numGroupAHits +"\tGroupAHits\t"+Misc.stringArrayListToString(p.groupAHits, ","));
			
			int numGroupANoHits = p.getGroupANoHits().size();
//IO.pl("\t"+numGroupANoHits +"\tGroupANoHits");
			
			int numGroupBHits = p.getGroupBHits().size();
//IO.pl("\t"+numGroupBHits +"\tGroupBHits\t"+Misc.stringArrayListToString(p.groupBHits, ","));
			
			int numGroupBNoHits = p.getGroupBNoHits().size();
//IO.pl("\t"+numGroupBNoHits +"\tGroupBNoHits");
			
			if (fe == null) fe = new FisherExact(numGroupAHits+ numGroupANoHits+ numGroupBHits+ numGroupBNoHits);
			double pVal = fe.getTwoTailedP(numGroupAHits, numGroupANoHits, numGroupBHits, numGroupBNoHits);
//IO.pl("\t"+pVal+"\t"+p.name);
			p.setPval(pVal);
			minLog10Pvals[i] = (float)Num.minus10log10(pVal);
		}
		
		//adjust the pvalues
		float[] adjPVals = Num.benjaminiHochbergCorrectUnsorted(minLog10Pvals);
		for (int i=0; i< numPathways; i++) {
			KeggPathway p = pathways.get(i);
			double adjPval = Num.antiNeg10log10(adjPVals[i]);
			if (adjPval>1.0) adjPval = 1.0;
			p.setAdjPval(adjPval);
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
				for (KeggPathway p: pathways) {
//IO.pl("Pathway "+p.name+" -> "+p.genes);
					boolean hit = false;
					//for each gene in the sample
					for (int i=1; i< cells.length; i++) {
//IO.pl("\tSampleGene\t"+cells[i]);
						if (p.getGenes().contains(cells[i])) {
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
//IO.pl("\t\t\tHIT");
							break;
						}
					}
					if (hit == false) {
//IO.pl("\t\tNoHit");
						if(isGroupA) p.getGroupANoHits().add(cells[0]);
						else p.getGroupBNoHits().add(cells[0]);
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
				if (cells.length > 1) {
					KeggPathway kp = new KeggPathway(cells, minKeggFreq);
					pathways.add(kp);
					if (networkIdPathway.containsKey(kp.getNetworkId())) throw new IOException("Duplicate Kegg Network ID "+kp.getNetworkId());
					networkIdPathway.put(kp.getNetworkId(), kp);
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the pathway file "+pathwayGenes);
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
					KeggPathway kp = networkIdPathway.get(cells[0]);
					if (kp == null) throw new IOException("Missing Kegg Network ID in the gene pathways file "+cells[0]);
					String[] pd = new String[] {cells[1], cells[2]};
					kp.addPathwayMapDescription(pd);
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse the pathway link file "+pathwayLinks);
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
					case 'l': pathwayLinks = new File(args[++i]); break;
					case 'r': resultsSpreadsheet = new File(args[++i]); break;
					case 'm': minKeggFreq = Double.parseDouble(args[++i]); break;
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
		
		if (groupAGeneHits == null || groupBGeneHits == null || pathwayGenes == null || pathwayLinks == null || resultsSpreadsheet == null) {
			Misc.printErrAndExit("\nERROR: failed to find one or more of the required input files. See the help menu.\n");
		}
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
				"pathway are saved. For TNRunner processed somatic variant files, use the USeq \n"+
				"AnnotatedVcfParser to select high impact, loss of function/ CLINVAR patho/likely-\n"+
				"pathogenic variants.\n"+
				
				"\nRequired Parameters:\n"+
				"-a File containing cohort A gene sets, each line represents a subject's genes of\n"+
				"     interest (e.g. those with HIGH impact mutations), tab delimited, the first cell\n"+
				"     is the subject ID, subsequent cells are the gene names.\n"+
				"-b File containing cohort B gene sets, ditto.\n"+
				"-p File containing KEGG MEDICUS pathways to compare, each line represents one pathway,\n"+
				"     tab delimited: unique pathway name, network ID, and associated gene symbols (e.g.\n"+
		        "     REFERENCE_ATR_SIGNALING  N01451 ATR ATRIP CHEK1 HUS1 ...)\n"+
		        "-l  File containing KEGG network IDs, pathway links, and their descriptions (e.g.\n"+
		        "    N01486 hsa04110+N01486 CellCycle), one network per row, tab delimited.\n"+
				"-r File to save the txt results, should end with .xls\n"+
		        "-m Minimum gene hit frequency for inclusion in output, defaults to 0.0005\n"+
				"-o Add one to zero count A or B fractions when calculating the log2Rto(fracA/fracB)\n"+
				
				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VariantPathwayComparatorKegg -a \n"+
				"   earlyCRC.txt -l keggLinks.txt -p keggPathways.txt -a eCRC.txt -b lCRC.txt -o\n"+
				"   -r earlyVsLateVPCKM.xls \n"+

		"\n**************************************************************************************\n");

	}
}
