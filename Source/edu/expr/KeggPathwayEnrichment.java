package edu.expr;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.parsers.*;
import util.bio.seq.Seq;
import util.gen.*;


/**Takes a list of select genes from Ensembl (ENSG00...) compares them to the total list and identifies Kegg pathways statistically overrepresented in the select list.
 * Needs three files from Kegg http://www.genome.jp/kegg/pathway.html
 * 1) KeggGeneIDs : EnsemblGeneIDs - ftp://ftp.genome.jp/pub/kegg/genes/organisms/hsa/hsa_ensembl-hsa.list 
 * 2) KeggPathwayIDs : TextDescription - ftp://ftp.genome.jp/pub/kegg/pathway/map_title.tab 
 * 3) KeggGeneIDs : KeggPathwayIDs - ftp://ftp.genome.jp/pub/kegg/pathway/organisms/hsa/hsa_gene_map.tab
 * Need two user files:
 * 1) All Ensemble genes interrogated in experiment, each text must start with ENSG...
 * 2) Select list of genes found interesting in experiment (enriched, reduced, mutated...)*/
public class KeggPathwayEnrichment {

	private File keggEnsemblMapFile;
	private HashMap<String,String> ensembleGeneID2keggGeneID;
	private HashMap<String,String> keggGeneID2ensembleGeneID;
	private File keggPathwayTitleFile;
	private File keggGeneIDsPathwayIDsFile;
	private HashMap<String,String[]> keggGeneID2Pathways;
	private File allGeneListFile;
	private HashSet<String> allGenes;
	private File selectGeneListFile;
	private HashSet<String> selectGenes;
	private int numberOfIterations = 10000;
	private int numberOfIterationsDouble;
	private Pattern tab = Pattern.compile("\\t");
	private Pattern tabColon = Pattern.compile("\\t|:");
	private Pattern whiteSpace = Pattern.compile("\\s");
	private HashMap<String,Pathway> pathways;
	private String organismCode = null;



	public KeggPathwayEnrichment(String[] args){

		processArgs(args);

		loadData();

		scoreRealData();

		scoreRandomData();

		printResults();

		System.out.println("\nDone!\n");
	}

	public void printResults(){
		//order Pathways
		ArrayList<Pathway> pathwaysAL = new ArrayList<Pathway>();
		Iterator<String> pIt = pathways.keySet().iterator();
		while (pIt.hasNext()){
			Pathway p = pathways.get(pIt.next());
			if (p.numberRealHits != 0) {
				p.score();
				pathwaysAL.add(p);
			}
		}
		Pathway[] pathwaysArray = new Pathway[pathwaysAL.size()];
		pathwaysAL.toArray(pathwaysArray);
		Arrays.sort(pathwaysArray);
		//write
		try {
			File results = new File (Misc.removeExtension(selectGeneListFile.getCanonicalPath())+"_KeggEnrich.xls");
			PrintWriter out = new PrintWriter (new FileWriter(results));
			out.println("PathwayID_Title\t#RealHits\tMeanRandomHits\tLog2((1+#RealHits)/(1+MeanRandom))\tpValue\tGeneHits2Pathway");
			for (int i=0; i< pathwaysArray.length; i++){
				out.println(pathwaysArray[i]);
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}

	}

	public void scoreRandomData(){
		System.out.println("\nScoring random data...");
		String[] allEnsGeneNames = Misc.hashSetToStringArray(allGenes);
		int numGenesToDraw = selectGenes.size();
		for (int i=0; i< numberOfIterations; i++){
			//draw set of random gene names
			Misc.randomize(allEnsGeneNames, System.currentTimeMillis());
			String[] randomGenes = new String[numGenesToDraw];
			System.arraycopy(allEnsGeneNames, 0, randomGenes, 0, numGenesToDraw);
			//convert to keggIDs
			String[] keggGeneIDs = convertGeneIDs(randomGenes);
			String[] pathwayIDs = findMemberPathways(keggGeneIDs, false);
			Arrays.sort(pathwayIDs);
			HashMap<String, Integer> pathwayIDsHitCount = hitCount(pathwayIDs);
			//add hits to pathways
			Iterator<String> it = pathwayIDsHitCount.keySet().iterator();
			while (it.hasNext()){
				String pathwayID= it.next();
				int numHits = pathwayIDsHitCount.get(pathwayID).intValue();
				pathways.get(pathwayID).numberRandomHits[i] = numHits;
			}
		}
	}

	public void scoreRealData(){
		System.out.println("\nScoring real data...");
		String[] realSelectGenes = Misc.hashSetToStringArray(selectGenes);
		String[] keggGeneIDs = convertGeneIDs(realSelectGenes);
		String[] pathwayIDs = findMemberPathways(keggGeneIDs, true);
		Arrays.sort(pathwayIDs);
		HashMap<String, Integer> pathwayIDsHitCount = hitCount(pathwayIDs);
		//add hits to pathways
		Iterator<String> it = pathwayIDsHitCount.keySet().iterator();
		while (it.hasNext()){
			String pathwayID= it.next();
			int numHits = pathwayIDsHitCount.get(pathwayID).intValue();
			pathways.get(pathwayID).numberRealHits = numHits;
		}
	}

	public HashMap<String, Integer> hitCount(String[] sorted){
		HashMap<String, Integer> hitCount = new HashMap<String, Integer>();
		String current = sorted[0];
		int numCounts = 1;
		for (int i=1; i< sorted.length; i++){
			String test = sorted[i];
			if (test.equals(current)) numCounts++;
			else {
				//add count
				hitCount.put(current, numCounts);
				//reset
				numCounts = 1;
				current = test;
			}
		}
		//add last
		hitCount.put(current, numCounts);
		return hitCount;
	}

	public String[] findMemberPathways(String[] keggGeneIDs, boolean setHitsInPathways){
		ArrayList<String> pathwaysAL = new ArrayList<String>(keggGeneIDs.length * 3);
		if (setHitsInPathways){
			for (int i=0; i< keggGeneIDs.length; i++){
				String[] ps = keggGeneID2Pathways.get(keggGeneIDs[i]);
				if (ps == null) continue;
				for (int j=0; j< ps.length; j++) {
					pathwaysAL.add(ps[j]);
					pathways.get(ps[j]).intersectingKeggGeneIDs.add(keggGeneIDs[i]);
				}
			}
		}
		else {
			for (int i=0; i< keggGeneIDs.length; i++){
				String[] ps = keggGeneID2Pathways.get(keggGeneIDs[i]);
				if (ps == null) continue;
				for (int j=0; j< ps.length; j++) {
					pathwaysAL.add(ps[j]);
				}
			}
		}
		return Misc.stringArrayListToStringArray(pathwaysAL);
	}

	public String[] convertGeneIDs(String[] ensGeneIDs){
		ArrayList<String> ids = new ArrayList<String>(ensGeneIDs.length);
		for (int i=0; i< ensGeneIDs.length; i++){
			String keggGeneID = ensembleGeneID2keggGeneID.get(ensGeneIDs[i]);
			if (keggGeneID != null) ids.add(keggGeneID);
		}
		return Misc.stringArrayListToStringArray(ids);
	}

	private class Pathway implements Comparable{
		String keggPathwayID;	
		String keggPathwayTitle;
		int numberRealHits;
		HashSet<String> intersectingKeggGeneIDs = new HashSet<String>();
		int[] numberRandomHits;
		double meanNumberHits;
		double pValue;
		double enrichment;

		public Pathway(String keggPathwayID, String keggPathwayTitle){
			this.keggPathwayID = keggPathwayID;
			this.keggPathwayTitle = keggPathwayTitle;
			numberRandomHits = new int[numberOfIterations];
		}

		public void score(){
			//calculate number of exceeding hits
			double numBigHits = 0;
			double totalNumberHits = 0;
			for (int i=0; i< numberRandomHits.length; i++){
				if (numberRandomHits[i] >= numberRealHits) numBigHits++;
				totalNumberHits += numberRandomHits[i];
			}
			meanNumberHits = totalNumberHits/numberOfIterationsDouble;
			pValue = numBigHits/numberOfIterationsDouble;
			enrichment = Num.log2((1.0 + (double)numberRealHits)/(1+meanNumberHits));
		}

		/**ID_Title NumRealHits MeanRandomHits EnrichmentLog2((1+NumReal)/(1+MeanRandom)) pVal GeneHits
		 * Only to be called after calling score() and with real hits.*/
		public String toString(){
			StringBuilder sb = new StringBuilder();
			//ID_Title =HYPERLINK("http://www.genome.jp/dbget-bin/www_bget?pathway+map00290","00290")
			sb.append("=HYPERLINK(\"http://www.genome.jp/dbget-bin/www_bget?pathway+map");
			sb.append(keggPathwayID);
			sb.append("\",\"");
			sb.append(keggPathwayID);
			sb.append("_");
			sb.append(keggPathwayTitle);
			sb.append("\")\t");
			//num real
			sb.append(numberRealHits);
			sb.append("\t");
			//num mean random
			sb.append(Num.formatNumber(meanNumberHits, 2));
			sb.append("\t");
			//enrichment
			sb.append(Num.formatNumber(enrichment, 3));
			sb.append("\t");
			//pValue
			sb.append(Num.formatNumber(pValue, 6));
			sb.append("\t");
			//genes  =HYPERLINK("http://www.genome.jp/dbget-bin/www_bfind_sub?max_hit=1000&dbkey=genes&keywords=hsa%3A1001&mode=bget","text")
			Iterator<String> it = intersectingKeggGeneIDs.iterator();
			while (it.hasNext()){
				String keggGeneID = it.next();
				String ensemblGeneID = keggGeneID2ensembleGeneID.get(keggGeneID);
				sb.append("=HYPERLINK(\"http://www.genome.jp/dbget-bin/www_bfind_sub?max_hit=1000&dbkey=genes&keywords=");
				sb.append(organismCode);
				sb.append("%3A");
				sb.append(keggGeneID);
				sb.append("&mode=bget\",\"");
				sb.append(ensemblGeneID);
				sb.append("\")\t");
			}

			return sb.toString();
		}

		public int compareTo(Object o){
			Pathway other = (Pathway)o;
			//sort by pval
			if (pValue < other.pValue) return -1;
			if (pValue > other.pValue) return 1;
			if (pValue == 0){
				if (enrichment > other.enrichment) return -1;
				if (enrichment < other.enrichment) return 1;
			}
			if (pValue == 1){
				if (enrichment > other.enrichment) return 1;
				if (enrichment < other.enrichment) return 1;
			}
			return 0;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KeggPathwayEnrichment(args);
	}		

	public void loadData(){
		System.out.println("Loading data...");
		loadEnsembleID2KeggGeneID();
		//System.out.println("\nensembleGeneID2keggGeneID "+ensembleGeneID2keggGeneID);

		loadKeggPathway2Title();
		//System.out.println("\nkeggPathwayID2Title "+keggPathwayID2Title);

		loadKeggGeneID2PathwayIDs();
		//System.out.println("\nkeggGeneID2Pathways "+keggGeneID2Pathways);

		int numberAllGenesInKeggPathways = loadAllGenes();
		//System.out.println("\nallGenes "+allGenes);

		int numberSelectGenesInKeggPathways = loadSelectGenes();
		//System.out.println("\nselectGenes "+selectGenes);

		//print stats
		System.out.println(selectGenes.size()+"\tUnique select genes");
		System.out.println(numberSelectGenesInKeggPathways+"\tUnique select genes in Kegg pathways");
		System.out.println(allGenes.size()+"\tUnique all genes");
		System.out.println(numberAllGenesInKeggPathways+"\tUnique all genes in Kegg pathways");

		if (numberSelectGenesInKeggPathways == 0) Misc.printExit("\nNo intersecting pathways in your select gene list!");

	}

	public int loadSelectGenes(){
		//load it
		String line = null;
		int numberMappedToPathway = 0;
		try{
			selectGenes = new HashSet<String>();
			BufferedReader in = IO.fetchBufferedReader(selectGeneListFile);
			// ENSG00000121410
			while ((line=in.readLine()) != null){
				if (line.startsWith("#")) continue;
				line = line.trim();
				if (allGenes.contains(line)) {
					selectGenes.add(line);
					if (ensembleGeneID2keggGeneID.containsKey(line)) {
						String keggGeneID = ensembleGeneID2keggGeneID.get(line);
						if (keggGeneID2Pathways.containsKey(keggGeneID)) numberMappedToPathway++;
					}
				}
				else Misc.printErrAndExit("\nThe following select gene is not found in your all gene list?! "+line);
			}
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem loading your select Ensembl gene file, see option -a "+selectGeneListFile);
		}
		return numberMappedToPathway;
	}

	public int loadAllGenes(){
		//load it
		String line = null;
		int numberMappedToPathway = 0;
		try{
			allGenes = new HashSet<String>();
			BufferedReader in = IO.fetchBufferedReader(allGeneListFile);
			// ENSG00000121410
			while ((line=in.readLine()) != null){
				if (line.startsWith("#")) continue;
				line = line.trim();
				if (ensembleGeneID2keggGeneID.containsKey(line)) {
					String keggGeneID = ensembleGeneID2keggGeneID.get(line);
					if (keggGeneID2Pathways.containsKey(keggGeneID)) numberMappedToPathway++;
				}
				allGenes.add(line.trim());
			}
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem loading your all interrogated Ensembl gene file, see option -a "+allGeneListFile);
		}
		return numberMappedToPathway;
	}

	public void loadKeggGeneID2PathwayIDs(){
		String line = null;
		try{
			keggGeneID2Pathways = new HashMap<String,String[]>();
			BufferedReader in = IO.fetchBufferedReader(keggGeneIDsPathwayIDsFile);
			String[] tokens;
			// 100130247	00190 01100 04260 05010 05012 05016
			while ((line=in.readLine()) != null){
				if (line.startsWith("#")) continue;
				tokens = tab.split(line);
				if (tokens.length != 2) throw new Exception();
				keggGeneID2Pathways.put(tokens[0], whiteSpace.split(tokens[1]));
			}
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem loading KeggGeneIDs : KeggPathwayIDs file, see option -g "+keggGeneIDsPathwayIDsFile+"\n\tLine "+line);
		}
	}

	public void loadKeggPathway2Title(){
		String line = null;
		try{
			pathways = new HashMap<String,Pathway>();
			BufferedReader in = IO.fetchBufferedReader(keggPathwayTitleFile);
			String[] tokens;
			// 07013	Cephalosporins - oral agents
			while ((line=in.readLine()) != null){
				if (line.startsWith("#")) continue;
				tokens = tab.split(line);
				if (tokens.length != 2) throw new Exception();
				pathways.put(tokens[0], new Pathway(tokens[0], tokens[1]));
			}
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem loading KeggPathwayIDs : TextTitle file, see option -p "+keggPathwayTitleFile+"\n\tLine "+line);
		}
	}

	public void loadEnsembleID2KeggGeneID(){
		String line = null;
		try{
			ensembleGeneID2keggGeneID = new HashMap<String,String>();
			keggGeneID2ensembleGeneID = new HashMap<String,String>();
			BufferedReader in = IO.fetchBufferedReader(keggEnsemblMapFile);
			String[] tokens = null;
			// hsa:100008586	ensembl-hsa:ENSG00000215274
			while ((line=in.readLine()) != null){
				if (line.startsWith("#")) continue;
				tokens = tabColon.split(line);
				if (tokens.length != 4) throw new Exception();
				ensembleGeneID2keggGeneID.put(tokens[3], tokens[1]);
				keggGeneID2ensembleGeneID.put(tokens[1], tokens[3]);
			}
			//set organismCode
			organismCode = tokens[0];
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem loading KeggGeneIDs : EnsemblGeneIDs file, see option -e "+keggEnsemblMapFile+"\n\tLine "+line);
		}
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
					case 'e': keggEnsemblMapFile = new File(args[++i]); break;
					case 'p': keggPathwayTitleFile = new File(args[++i]); break;
					case 'g': keggGeneIDsPathwayIDsFile = new File(args[++i]); break;
					case 'a': allGeneListFile = new File(args[++i]); break;
					case 's': selectGeneListFile = new File(args[++i]); break;
					case 'n': numberOfIterations = Integer.parseInt(args[++i]);
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check files
		if (keggEnsemblMapFile == null || keggEnsemblMapFile.canRead()==false) Misc.printErrAndExit("\nCan't find or read your KeggGeneIDs : EnsemblGeneIDs file, see option -e "+keggEnsemblMapFile);
		if (keggPathwayTitleFile == null || keggPathwayTitleFile.canRead()==false) Misc.printErrAndExit("\nCan't find or read your KeggPathwayIDs : TextTitle file, see option -p "+keggPathwayTitleFile);
		if (keggGeneIDsPathwayIDsFile == null || keggGeneIDsPathwayIDsFile.canRead()==false) Misc.printErrAndExit("\nCan't find or read your KeggGeneIDs : KeggPathwayIDs file, see option -g "+keggGeneIDsPathwayIDsFile);
		if (allGeneListFile == null || allGeneListFile.canRead()==false) Misc.printErrAndExit("\nCan't find or read your all interrogated Ensembl gene file, see option -a "+allGeneListFile);
		if (selectGeneListFile == null || selectGeneListFile.canRead()==false) Misc.printErrAndExit("\nCan't find or read your select Ensembl gene file file, see option -s "+selectGeneListFile);

		numberOfIterationsDouble = numberOfIterations;

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Kegg Pathway Enrichment:  Aug 2009                     **\n" +
				"**************************************************************************************\n" +
				"KPE looks for overrepresentation of genes from a user's list in Kegg pathways using a\n" +
				"random permutation test. Several files are needed from http://www.genome.jp/kegg \n" +
				"Gene names must be in Ensembl Gene notation and begin with ENSG.\n\n" +

				"Options:\n"+
				"-e Full path file text for a KeggGeneIDs : EnsemblGeneIDs file (e.g. Human \n" +
				"      ftp://ftp.genome.jp/pub/kegg/genes/organisms/hsa/hsa_ensembl-hsa.list)\n"+
				"-p Full path file text for a KeggPathwayIDs : TextDescription file (e.g. Human \n" +
				"      ftp://ftp.genome.jp/pub/kegg/pathway/map_title.tab)\n"+
				"-g Full path file text for a KeggGeneIDs : KeggPathwayIDs file (e.g. Human \n" +
				"      ftp://ftp.genome.jp/pub/kegg/pathway/organisms/hsa/hsa_gene_map.tab)\n"+
				"-a Full path file text for your all interrogated Ensembl gene list (e.g. ENSG00...)\n"+
				"      One gene per line.\n"+
				"-s Full path file text for your select gene list.\n" +
				"-n Number of random iterations, defaults to 10000\n"+
				
				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/KeggPathwayEnrichment -e \n" +
				"      /Kegg/hsa_ensembl-hsa.list -p /Kegg/map_title.tab -g /Kegg/hsa_gene_map.tab\n" +
				"      -a /HCV/ensemblGenesWith20OrMoreReads.txt -s /HCV/upRegInHCV_Norm.txt\n\n" +

		"**************************************************************************************\n");

	}
}
