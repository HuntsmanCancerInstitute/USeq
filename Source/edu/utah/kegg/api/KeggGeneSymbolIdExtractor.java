package edu.utah.kegg.api;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class KeggGeneSymbolIdExtractor {

	//fields
	private File resultsDirectory = null;
	private File interrogatedGenes = null;
	
	//internal fields
	private String findUrl = "https://rest.kegg.jp/find/genes/hsa%20";
	private ArrayList<KeggApiGeneMatch> geneMatches = new ArrayList<KeggApiGeneMatch>();
	private ArrayList<KeggApiGeneMatch> geneNoMatches = new ArrayList<KeggApiGeneMatch>();
	private File lookupFile = null;
	
	public KeggGeneSymbolIdExtractor(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);
		
		try {
			findGenes();
			saveResourceFiles();
			loadCheckHashes();
			
		} catch (Exception e) {
			IO.el("\nERROR running KeggGeneSymbolIDExtractor\n");
			e.printStackTrace();
		}

		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		IO.pl("\nDone! "+Math.round(diffTime)+" minutes\n");
	}
	
	private void loadCheckHashes() throws IOException {
		IO.pl("\nLoading and checking the lookup hashes...");
		HashMap<String, ArrayList<String>>[] gsTKi = KeggGeneSymbolIdExtractor.loadGeneLookupHashes(lookupFile);

		IO.pl("Dups in GeneSymbol 2 KeggID:");
		for (String geneSymbol: gsTKi[0].keySet()) {
			if (gsTKi[0].get(geneSymbol).size()>1) IO.pl("\t"+geneSymbol +" -> "+gsTKi[0].get(geneSymbol));
		}
		IO.pl("Dups in KeggID 2 GeneSymbol:");
		for (String ki: gsTKi[1].keySet()) {
			if (gsTKi[1].get(ki).size()>1) IO.pl("\t"+ki +" -> "+gsTKi[1].get(ki));
		}
		
	}

	/**Run on the geneSymbol2KeggGeneInfo27Jan2025.txt file, returns a GeneSymbol2KeggGeneId and a KeggGeneId2GeneSymbol HashMaps
	 * Won't load any null info.  
	 * @throws IOException */
	public static HashMap<String, ArrayList<String>>[] loadGeneLookupHashes( File geneSym2KegGenInfo) throws IOException{
		BufferedReader in = new BufferedReader( new FileReader( geneSym2KegGenInfo));
		String[] f;
		String l;
		HashMap<String, ArrayList<String>> gs2ki = new HashMap<String,ArrayList<String>>();
		HashMap<String, ArrayList<String>> ki2gs = new HashMap<String,ArrayList<String>>();
		while ((l=in.readLine()) != null) {
			//AAR2	25980	hsa:25980	AAR2, C20orf4, CGI-23; AAR2 splicing factor
			// 0       1       2		  3
			if (l.startsWith("#") == false) {
				f = Misc.TAB.split(l);
				if (f[1].equals("null")) continue;
				ArrayList<String> al = gs2ki.get(f[0]);
				if (al==null) {
					al = new ArrayList<String>();
					gs2ki.put(f[0], al);
				}
				al.add(f[1]);
				
				ArrayList<String> al2 = ki2gs.get(f[1]);
				if (al2==null) {
					al2 = new ArrayList<String>();
					ki2gs.put(f[1], al2);
				}
				al2.add(f[0]);
			}
		}
		in.close();
		@SuppressWarnings("unchecked")
		HashMap<String, ArrayList<String>>[] toReturn = (HashMap<String, ArrayList<String>>[]) new HashMap[2];
		toReturn[0] = gs2ki;
		toReturn[1] = ki2gs;
		return toReturn;
		
	}
	

	private void saveResourceFiles() throws FileNotFoundException {
		String fileName = "geneSymbol2KeggGeneInfo"+Misc.getDateNoSpaces()+".txt";
		
		//write out no match, then matches
		lookupFile = new File(resultsDirectory, fileName);
		IO.pl("\nSaving Gene Symbol 2 Kegg Gene Info: "+lookupFile);
		PrintWriter out = new PrintWriter(lookupFile);
		out.println("#GeneSymbol\tKeggGeneId\tMatchingKeggQueryResponse");
		for (KeggApiGeneMatch gm: geneNoMatches) out.println(gm.toStringSummary());
		for (KeggApiGeneMatch gm: geneMatches) out.println(gm.toStringSummary());
		out.close();
	}

	private void findGenes() throws IOException {
		IO.pl("\nFetching kegg genes...");
		String[] geneSymbolsToLookUp = IO.loadFileIntoStringArray(interrogatedGenes);
		
		File genesDir = new File(resultsDirectory, "Genes");
		genesDir.mkdirs();
		
		int numCached = 0;
		int numDownloaded = 0;
		
		//for each network id
		for (String nId: geneSymbolsToLookUp) {
			File tmp = new File(genesDir, nId+ ".txt");
			KeggApiGeneMatch ni = null;
			//does it already exist?
			if (tmp.exists()) {
				numCached++;
				IO.pl("\t"+nId+" exists using cached file. Delete to reload.");
				ni = new KeggApiGeneMatch(tmp, nId);
			}
			//fetch it
			else {
				String[] cmd = new String[] {"curl","-o",tmp.getCanonicalPath(), findUrl+nId};
				int exit = IO.executeViaProcessBuilderReturnExit(cmd);
				if (exit !=0) throw new IOException("ERROR: failed to fetch the gene info:\n"+Misc.stringArrayToString(cmd, " ")+"\n");
				numDownloaded++;
				ni = new KeggApiGeneMatch(tmp, nId);
				IO.pl("\n"+ ni.toString());
			}
			
			if (ni.isOk()) geneMatches.add(ni);
			else geneNoMatches.add(ni);
		}
		IO.pl("\t"+numDownloaded+" Downloaded, "+numCached+" Cached (delete "+genesDir+" to force a fresh download and restart.)");
		IO.pl("\t"+geneMatches.size()+" Gene symbols matched to Kegg gene ids.");
		IO.p("\t"+geneNoMatches.size()+" Gene symbols not matched:\n\t\t");
		for (KeggApiGeneMatch g: geneNoMatches) IO.p(g.getGeneSymbol()+" ");
		IO.pl();
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KeggGeneSymbolIdExtractor(args);
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
					case 'g': interrogatedGenes = new File(args[++i]); break;
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
		if (interrogatedGenes == null || interrogatedGenes.exists() == false) Misc.printErrAndExit("\nError: provide a file containing gene symbols, one per line, to match to Kegg Gene IDs.\n");
	}	


	public static void printDocs(){
		
		
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Kegg Gene Symbol Id Extractor: Jan 2025                  **\n" +
				"**************************************************************************************\n" +
				"Uses the Kegg API to look up and then parse the Kegg Gene Id that most closely\n"+
				"matches each of the provided HUGO Gene Symbols. Delete the xxx/Genes directory or \n"+
				"files within to load the most recent gene info. See https://rest.kegg.jp/find/genes/\n"+

				"\nOptions:\n"+
				"-g File containing HUGO Gene Symbols to look up, one per line.\n"+
				"-r Results directory to save the Kegg resource files.\n"+
				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/KeggGeneSymbolIdExtractor -r KeggResources\n"+
				"         -g interrogatedGeneSymbolNames.txt\n\n" +

				"**************************************************************************************\n");
	}
}
