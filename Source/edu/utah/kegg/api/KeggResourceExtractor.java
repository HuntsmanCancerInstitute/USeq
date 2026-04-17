package edu.utah.kegg.api;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.hci.misc.Util;
import util.gen.IO;
import util.gen.Misc;

public class KeggResourceExtractor {

	//fields
	private int minimumNumberGenes = 4;
	private boolean requirePathways = true;
	private File resultsDirectory = null;
	private String networkListUrl = "https://rest.kegg.jp/list/network";
	private String networkGetUrl = "https://rest.kegg.jp/get/ne:";
	private Random random = new Random();
	private File tmpDir = new File(System.getProperty("java.io.tmpdir"));
	
	//results
	private ArrayList<String> networkIds = new ArrayList<String>();
	private ArrayList<KeggApiNetwork> networks = new ArrayList<KeggApiNetwork>();
	private ArrayList<KeggApiNetwork> passingNetworks = new ArrayList<KeggApiNetwork>();
	
	public KeggResourceExtractor(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);
		
		try {
			fetchNetworkIds();
			fetchNetworkInfo();
			filterNetworks();
			saveResourceFilesLegacy();
			
		} catch (Exception e) {
			IO.el("\nERROR running KeggResourceExtractor\n");
			e.printStackTrace();
		}

		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		IO.pl("\nDone! "+Math.round(diffTime)+" minutes\n");
	}
	

	private void saveResourceFilesLegacy() throws FileNotFoundException {
		IO.pl("\nSaving Legacy Pathway Links...");
		String reqPath = "T_";
		if (this.requirePathways == false) reqPath = "F_";
		
		String date = "Min"+minimumNumberGenes+ reqPath +Misc.getDateNoSpaces()+".txt";
		
		//write out pathway links
		File links = new File(resultsDirectory, "keggMedicusPathwayLinks_"+date);
		PrintWriter out = new PrintWriter(links);
		out.println("# Kegg Network TxtInfo: https://www.kegg.jp/entry/hsa05221+N00003 VisPath: https://www.kegg.jp/pathway/hsa05221+N00003");
		for (KeggApiNetwork ni: passingNetworks) {
			//N00003	hsa05221+N00003	Acute myeloid leukemia, one or more
			out.print(ni.getPathwayLinks());
		}
		out.close();
		IO.pl("\t"+links);
		
		//write out pathway genes
		File genes = new File(resultsDirectory, "keggMedicusPathwayGenes_"+date);
		out = new PrintWriter(genes);
		out.println("# NameOfPathwayTested https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/KEGG_MEDICUS_	KeggNetwork https://www.kegg.jp/entry/	Genes…");
		for (KeggApiNetwork ni: passingNetworks) {
			//REFERENCE_WNT_SIGNALING_PATHWAY	N00056	APC	AXIN1	AXIN2	CCND1	CTNNB1	DVL1...
			out.println(ni.getPathwayGenes());
		}
		out.close();
		IO.pl("\t"+genes);
	}


	private void filterNetworks() {
		IO.pl("\nFiltering pathways for those with a minimum # genes and whether they have associated pathways...");
		TreeSet<String> uniquePathways = new TreeSet<String>();
		TreeSet<String> uniqueGenes = new TreeSet<String>();
		
		for (KeggApiNetwork ni: networks) {
			if (ni.getGenes()== null || ni.getGenes().length < minimumNumberGenes) continue;
			if (requirePathways == true && ni.getPathways() == null ) continue;
			passingNetworks.add(ni);
			//IO.pl(ni+"\n");
			if (ni.getPathways()!=null) for (KeggApiPathway kp: ni.getPathways()) uniquePathways.add(kp.getId());
			for (KeggApiGene kp: ni.getGenes()) uniqueGenes.add(kp.getId());
		}
		IO.pl("\t"+passingNetworks.size()+" of "+networks.size()+" networks passed (>="+minimumNumberGenes+" #Genes, "+ requirePathways+" reqPath)");
		IO.pl("\t"+uniquePathways.size()+" unique pathways");
		IO.pl("\t"+uniqueGenes.size()+" unique genes");
	}

	private void fetchNetworkInfo() throws IOException {
		IO.pl("\nFetching and parsing network info...");
		File networkInfoDir = new File(resultsDirectory, "Networks");
		networkInfoDir.mkdirs();
		
		int numCached = 0;
		int numDownloaded = 0;
		
		//for each network id
		for (String nId: networkIds) {
			File tmp = new File(networkInfoDir, nId+ ".txt");
			//does it already exist?
			if (tmp.exists()) {
				numCached++;
				//IO.pl("\t"+nId+" exists using cached file. Delete to reload.");
			}
			//fetch it
			else {
				String[] cmd = new String[] {"curl","-o",tmp.getCanonicalPath(), networkGetUrl+nId};
				int exit = IO.executeViaProcessBuilderReturnExit(cmd);
				if (exit !=0) throw new IOException("ERROR: failed to fetch the network info:\n"+Misc.stringArrayToString(cmd, " ")+"\n");
				numDownloaded++;
			}
			KeggApiNetwork ni = new KeggApiNetwork(tmp, nId);
			networks.add(ni);
		}
		IO.pl("\t"+numDownloaded+" Downloaded, "+numCached+" Cached (delete "+networkInfoDir+" to force a fresh download and restart.)");
	}
	
	public static KeggApiNetwork[] loadNetworks(File networkDir) throws IOException {
		File[] networkFiles = IO.extractFilesStartingWith(networkDir, "N");
		KeggApiNetwork[] networks = new KeggApiNetwork[networkFiles.length];
		for (int i=0; i< networks.length; i++) networks[i]= new KeggApiNetwork(networkFiles[i], networkFiles[i].getName().replaceAll(".txt", ""));
		return networks;
	}
	
	public static KeggApiNetwork[] buildCompositePathways (KeggApiNetwork[] networks) {
		
		HashMap<String, ArrayList<KeggApiNetwork>> pathwayIdNetworks = new HashMap<String, ArrayList<KeggApiNetwork>>();
		HashMap<String, KeggApiPathway> pathwayIdPathway = new HashMap<String, KeggApiPathway>();
		
		//for each network
		for (int i=0; i< networks.length; i++) {
			//for each pathway
			KeggApiPathway[] pathways = networks[i].getPathways();
			if (pathways !=null) {
				for (int j=0; j < pathways.length; j++) {
					ArrayList<KeggApiNetwork> al = pathwayIdNetworks.get(pathways[j].getId());
					if (al == null) {
						al = new ArrayList<KeggApiNetwork>();
						pathwayIdNetworks.put(pathways[j].getId(), al);
					}
					al.add(networks[i]);
					pathwayIdPathway.put(pathways[j].getId(), pathways[j]);
				}
			}
		}
		
		//for each pathway, make a composite 'network'
		KeggApiNetwork[] compNets = new KeggApiNetwork[pathwayIdNetworks.size()];
		
		int counter = 0;
		for (String pathwayId: pathwayIdNetworks.keySet()) {
			ArrayList<KeggApiNetwork> nets = pathwayIdNetworks.get(pathwayId);
			
			String id = pathwayId;
			String name = nets.get(0).getPathways()[0].getName();
			HashSet<String> types = new HashSet<String>();
			KeggApiPathway[] pathways = new KeggApiPathway[] {pathwayIdPathway.get(pathwayId)};
			HashMap<String, KeggApiGene> geneNameKeggApiGene = new HashMap<String, KeggApiGene>();
			
			//for each network, add in genes and types
			for (KeggApiNetwork n: nets) {
				if (n.getGeneNameKeggApiGene()!=null) geneNameKeggApiGene.putAll(n.getGeneNameKeggApiGene());
				if (n.getNetworkType()!= null) types.add(n.getNetworkType());
			}
			
			KeggApiGene[] genes = new KeggApiGene[geneNameKeggApiGene.size()];
			int i = 0;
			for (KeggApiGene g: geneNameKeggApiGene.values()) genes[i++] = g;
			
			KeggApiNetwork n = new KeggApiNetwork(id, name, Misc.hashSetToString(types, ","), pathways, geneNameKeggApiGene, genes);
			compNets[counter++] = n;
//IO.pl(n+"\n");
		}
		return compNets;
	}


	/*
N00001	EGF-EGFR-RAS-ERK signaling pathway
N00002	BCR-ABL fusion kinase to RAS-ERK signaling pathway
N00003	Mutation-activated KIT to RAS-ERK signaling pathway
...
nt06031	Citrate cycle and pyruvate metabolism
nt06017	Glycogen metabolism
nt06023	Galactose degradation
	 */
	private void fetchNetworkIds() throws IOException {
		IO.pl("Fetching and parsing network IDs...");
		File tmpList = new File(tmpDir,"tmpNetworkList."+ random.nextInt(10000)+ ".txt");
		tmpList.deleteOnExit();
		
		String[] cmd = new String[] {"curl","-o",tmpList.getCanonicalPath(), networkListUrl};
		int exit = IO.executeViaProcessBuilderReturnExit(cmd);
		if (exit !=0) throw new IOException("ERROR: failed to fetch the network list:\n"+Misc.stringArrayToString(cmd, " ")+"\n");
		
		String[] list = IO.loadFileIntoStringArray(tmpList);
		
		for (int i=0; i< list.length; i++) {
			String id = Misc.TAB.split(list[i])[0];
			if (id.startsWith("N")) networkIds.add(id);
		}
		IO.pl("\t"+networkIds.size()+"\tNetworks found");
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KeggResourceExtractor(args);
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
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'm': minimumNumberGenes = Integer.parseInt(args[++i]); break;
					case 'p': requirePathways = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		
		//check results dir
		if (resultsDirectory == null) Misc.printErrAndExit("\nError: provide a directory to save results.\n");
	}	


	public static void printDocs(){
		
		
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Kegg Resource Extractor: Jan 2025                        **\n" +
				"**************************************************************************************\n" +
				"Pulls the latest list of KEGG Networks, their associated genes, and pathways.\n"+
				"Writes out reference files for the USeq KEGG pathway applications.\n"+

				"\nOptions:\n"+
				"-r Results directory to save the Kegg resource files.\n"+
				"-m (Optional) Minimum # genes in network, defaults to 4.\n"+
				"-p (Optional) Save networks lacking pathways, defaults to excluding.\n"+
				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/KeggResourceExtractor -r KeggResources\n"+
				" \n\n" +

				"**************************************************************************************\n");
	}
}
