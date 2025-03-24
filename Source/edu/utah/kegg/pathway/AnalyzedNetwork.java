package edu.utah.kegg.pathway;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

import edu.utah.kegg.api.KeggApiNetwork;
import edu.utah.kegg.api.KeggApiPathway;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class AnalyzedNetwork implements Comparable<AnalyzedNetwork>{
	
	//Either Gene or Variant Info
	private KeggApiNetwork keggApiNetwork = null;
	private ArrayList<KeggApiNetwork> analizedKeggApiNetworks = new ArrayList<KeggApiNetwork>();
	private double pValue = -1;
	private double fdr = -1;
	public double sortValue = -1;
	private double combinePValue = -1;
	private double combineFdr = -1;
	
	//Gene Specific Info
	private int totalNumDiffExpGenes = 0;
	// users interrogated genes that intersect this network's genes
	private TreeSet<String> geneNetworkGeneNamesInInterrogatedGenes = new TreeSet<String>();
	// users select genes that intersect this network's genes
	private ArrayList<SelectGene> geneSharedGenes = new ArrayList<SelectGene>();
	private boolean geneAnalyzed = true;

	
	//Variant Specific Info
	private ArrayList<String> groupAHits = new ArrayList<String>();
	private ArrayList<String> groupANoHits = new ArrayList<String>();
	private ArrayList<String> groupBHits = new ArrayList<String>();
	private ArrayList<String> groupBNoHits = new ArrayList<String>();
	private double variantLog2Rto = -1;
	private HashMap<String, Integer> aGenes = new HashMap<String, Integer>();
	private HashMap<String, Integer> bGenes = new HashMap<String, Integer>();
	private double minKeggFreq = 0;
	private HashMap<String, Double> geneFractionsGroupA = null;
	private HashMap<String, Double> geneFractionsGroupB = null;
	
	//Kegg info, hex color coding, replace # with %23	
	public static String geneNegColor = "%2387CEEB";	// light blue
	public static String genePosColor = "%23FFB6C1";	// light red, pink
	public static String geneZeroColor = "%23FFB6C1";	// zero value genes don't exist, so not ever used
	public static String variantNegColor = "%23D6B4FC";	// light violet
	public static String variantPosColor = "%23FFD580";	// light orange
	public static String variantZeroColor = "%23FFFFAC";	// light yellow, sometimes set when the pos/neg rto is very small

	public static final String showPathwayUrl = "https://www.kegg.jp/kegg-bin/show_pathway?map=";
	
	
	public AnalyzedNetwork (KeggApiNetwork keggApiNetwork, int totalNumDiffExpGenes) {
		this.keggApiNetwork = keggApiNetwork;
		this.totalNumDiffExpGenes = totalNumDiffExpGenes;
	}
	
	public String toStringVariant(boolean addOne, HashMap<String, ArrayList<String>> geneSymbolKeggId) {
		//name descLink \tPval\tAdjPval\tAHits\tANoHits\tFracAHits\tAGeneHits\tBHits\tBNoHits\tFracBHits\tBGeneHits\tLog2(fracA/fracB)\tAllGeneHits\tPathwayMapLinks.."
		
		String[] networkIds = fetchNetworkIds();
		String mergedNetIds = Misc.stringArrayToString(networkIds, ",");
		
		StringBuilder sb = new StringBuilder();
		
		//for each network
		for (KeggApiNetwork n: analizedKeggApiNetworks) {

			//name
			sb.append(mergedNetIds);
			sb.append("\t");

			sb.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			sb.append(n.getNetworkId());		
			sb.append("\",\"");
			sb.append(n.getNetworkId());
			sb.append(": ");
			sb.append(n.getNetworkName());	
			sb.append("\")");
			sb.append("\t");

		sb.append(Num.formatNumber(pValue, 5)); sb.append("\t");

		sb.append(Num.formatNumber(fdr, 5)); sb.append("\t");

		sb.append(groupAHits.size()); sb.append("\t");

		sb.append(groupANoHits.size()); sb.append("\t");

		double aFrac = (double)groupAHits.size()/ (double)(groupAHits.size()+groupANoHits.size());

		sb.append(Num.formatNumber(aFrac, 5)); sb.append("\t");
		sb.append(convertToSortedFrequencies(aGenes, true)); sb.append("\t");

		sb.append(groupBHits.size()); sb.append("\t");
		sb.append(groupBNoHits.size()); sb.append("\t");
		double bFrac = (double)groupBHits.size()/ (double)(groupBHits.size()+groupBNoHits.size());
		sb.append(Num.formatNumber(bFrac, 5)); sb.append("\t");
		sb.append(convertToSortedFrequencies(bGenes, false)); sb.append("\t");

		//correct for zero log2Rto?
		if (addOne) {
			if (aFrac == 0) {
				double a = 1;
				aFrac = a/ (a+(double)groupANoHits.size());
			}
			if (bFrac == 0) {
				double b = 1;
				bFrac = b/ (b+(double)groupBNoHits.size());
			}
		}
		variantLog2Rto = aFrac/bFrac;
		sb.append(Num.formatNumber(Num.log2(variantLog2Rto), 5));
		sb.append("\t");

		TreeSet<String> allGenes = new TreeSet<String>();
		allGenes.addAll(aGenes.keySet());
		allGenes.addAll(bGenes.keySet());
		for (String geneName: allGenes) {
			sb.append(geneName);
			sb.append(" ");
		}

		
		
		//PathwayMapLinks...
		if (n.getPathways()!= null) {
			KeggApiPathway[] pathways = n.getPathways(); 
			SelectGene[] selectGenes = fetchVariantGenes();

			for (int i=0; i< pathways.length; i++) {
				sb.append("\t");
				String link = fetchKeggPathwayMapLink(pathways[i].getId(), networkIds, selectGenes, geneSymbolKeggId,variantNegColor ,variantPosColor ,variantZeroColor);
				if (link.length()<250) {
					sb.append("=HYPERLINK(\"");
					sb.append(link);
					sb.append("\",\"");
					sb.append(pathways[i].getName());
					sb.append("\")");
				}
				else {
					sb.append(pathways[i].getName());
					sb.append(" : ");
					sb.append(link);
				}
			}
		}
		else sb.append("\t");
		sb.append("\n");
		}
		return sb.toString();
	}

	public String toStringVariantCombo(String selectNetId, AnalyzedNetwork geneAnalyzedNetowrk, boolean addOne, HashMap<String, ArrayList<String>> geneSymbolKeggId) {
		//Pval AdjPval Log2(fracA/fracB) AllGeneHits PathwayMapLinks.."
		
		StringBuilder sb = new StringBuilder();
		
		//find select network
		for (KeggApiNetwork n: analizedKeggApiNetworks) {
			if (selectNetId.equals(n.getNetworkId())== false) continue;

		sb.append(Num.formatNumber(pValue, 5)); sb.append("\t");
		sb.append(Num.formatNumber(fdr, 5)); sb.append("\t");

		double aFrac = (double)groupAHits.size()/ (double)(groupAHits.size()+groupANoHits.size());
		double bFrac = (double)groupBHits.size()/ (double)(groupBHits.size()+groupBNoHits.size());
		
		//correct for zero log2Rto?
		if (addOne) {
			if (aFrac == 0) {
				double a = 1;
				aFrac = a/ (a+(double)groupANoHits.size());
			}
			if (bFrac == 0) {
				double b = 1;
				bFrac = b/ (b+(double)groupBNoHits.size());
			}
		}
		variantLog2Rto = aFrac/bFrac;
		sb.append(Num.formatNumber(Num.log2(variantLog2Rto), 5));
		sb.append("\t");

		TreeSet<String> allGenes = new TreeSet<String>();
		allGenes.addAll(aGenes.keySet());
		allGenes.addAll(bGenes.keySet());
		for (String geneName: allGenes) {
			sb.append(geneName);
			sb.append(" ");
		}

		//PathwayMapLinks...
		if (n.getPathways()!= null) {
			KeggApiPathway[] pathways = n.getPathways(); 
			
			SelectGene[] varSelectGenes = fetchVariantGenes();
			
			ArrayList<SelectGene> gsgAl = geneAnalyzedNetowrk.getGeneSharedGenes();
			SelectGene[] geneSelectGenes = new SelectGene[gsgAl.size()];
			gsgAl.toArray(geneSelectGenes);

			for (int i=0; i< pathways.length; i++) {
				sb.append("\t");
				
				String link = fetchCombineKeggPathwayMapLink(pathways[i].getId(), new String[] {selectNetId}, geneSelectGenes, varSelectGenes, geneSymbolKeggId);
				if (link.length()<250) {
					sb.append("=HYPERLINK(\"");
					sb.append(link);
					sb.append("\",\"");
					sb.append(pathways[i].getName());
					sb.append("\")");
				}
				else {
					sb.append(pathways[i].getName());
					sb.append(" : ");
					sb.append(link);
				}
			}
		}
		else sb.append("\t");
		sb.append("\n");
		}
		return sb.toString();
	}

	public SelectGene[] fetchVariantGenes() {
		//get all Variant genes
		TreeSet<String> allGeneHits = getVarinatGeneNameHits();
		SelectGene[] selectGenes = new SelectGene[allGeneHits.size()];
		int i = 0;
		for (String geneSymbol: allGeneHits) {
			Double aFrac = geneFractionsGroupA.get(geneSymbol);
			Double bFrac = geneFractionsGroupB.get(geneSymbol);
			if (aFrac == null) aFrac = 0d;
			if (bFrac == null) bFrac = 0d;
			
			//are either zero? thus no observations
			if (aFrac==0 || bFrac ==0) selectGenes[i++] = new SelectGene(geneSymbol,"0");
			else {
				double rto = Num.log2(aFrac/bFrac);
				//too small a diff? 1.25x
				if (Math.abs(rto) < 0.322 ) selectGenes[i++] = new SelectGene(geneSymbol,"0");
				else if(rto > 0) selectGenes[i++] = new SelectGene(geneSymbol,"1");
				else selectGenes[i++] = new SelectGene(geneSymbol,"-1");
			}
		}
		return selectGenes;
	}

	public String convertToSortedFrequencies(HashMap<String,Integer> keyCounts, boolean isGroupA) {
		double totalCounts = 0;
		for (String key: keyCounts.keySet()) totalCounts+= keyCounts.get(key);
		GeneCount[] gc = new GeneCount[keyCounts.size()];
		HashMap<String, Double> geneFraction = new HashMap<String, Double>(); 
		int count = 0;
		for (String key: keyCounts.keySet()) {
			double geneCount = keyCounts.get(key);
			if (totalCounts ==0) geneCount = 0;
			else geneCount = geneCount/totalCounts;
			gc[count++] = new GeneCount(key, geneCount);
			geneFraction.put(key, geneCount);
		}
		Arrays.sort(gc);
		if (isGroupA) geneFractionsGroupA = geneFraction;
		else geneFractionsGroupB = geneFraction;

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
		return namesToView.toString();
	}

	public String toStringGene(HashMap<String, ArrayList<String>> geneSymbolKeggId) {
		//name descLink pval, adjPval, #UniquePathwayGenes, #FoundUniquePathwayGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/FoundUniPathGenes, 
		//IntersectingGenes, IntersectingGeneLogRtos, PathwayMapLinks...
		
		String[] networkIds = fetchNetworkIds();
		String mergedNetIds = Misc.stringArrayToString(networkIds, ",");
		
		StringBuilder sb = new StringBuilder();
		
		//for each network
		for (KeggApiNetwork n: analizedKeggApiNetworks) {

			//name
			sb.append(mergedNetIds);
			sb.append("\t");

			sb.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			sb.append(n.getNetworkId());		
			sb.append("\",\"");
			sb.append(n.getNetworkId());
			sb.append(": ");
			sb.append(n.getNetworkName());	
			sb.append("\")");
			sb.append("\t");

			//pval, adjPval
			sb.append(pValue); sb.append("\t");
			sb.append(Num.formatNumber(fdr, 5)); sb.append("\t");

			// #UniquePathwayGenes
			sb.append(n.getGenes().length); sb.append("\t");

			//#FoundUniquePathwayGenes
			sb.append(geneNetworkGeneNamesInInterrogatedGenes.size()); sb.append("\t");

			//#DiffExpGenes
			sb.append(totalNumDiffExpGenes); sb.append("\t");

			//#GenesIntersect
			sb.append(geneSharedGenes.size()); sb.append("\t");

			//GenesIntersect/FoundUniPathGenes
			double frac = (double)geneSharedGenes.size()/ (double)geneNetworkGeneNamesInInterrogatedGenes.size();
			sb.append(Num.formatNumber(frac, 5)); sb.append("\t");

			//IntersectingGenes
			//sort the sharedGenes by gene symbol so easier to group
			SelectGene[] sortedSGs = new SelectGene[geneSharedGenes.size()];
			geneSharedGenes.toArray(sortedSGs);
			Arrays.sort(sortedSGs);
			//just gene symbols
			for (SelectGene sg: sortedSGs) {
				sb.append(sg.getGeneSymbol());
				sb.append(" ");
			}
			sb.append("\t"); 

			//IntersectingGeneLogRtos
			//just log2Rtos have to use , to keep excel from converting it to a formula
			if (sortedSGs.length!=0) {
				sb.append(sortedSGs[0].getLog2Rto());
				for (int i=1; i< sortedSGs.length; i++) {
					sb.append(",");
					sb.append(sortedSGs[i].getLog2Rto());
				}
			}
			//PathwayMapLinks...
			if (n.getPathways()!= null) {
				KeggApiPathway[] pathways = n.getPathways(); 
				SelectGene[] intersectingGenes = new SelectGene[geneSharedGenes.size()];
				geneSharedGenes.toArray(intersectingGenes);

				for (int i=0; i< pathways.length; i++) {
					sb.append("\t");
					String link = fetchKeggPathwayMapLink(pathways[i].getId(), networkIds, intersectingGenes, geneSymbolKeggId,geneNegColor ,genePosColor ,geneZeroColor);
					if (link.length()<250) {
						sb.append("=HYPERLINK(\"");
						sb.append(link);
						sb.append("\",\"");
						sb.append(pathways[i].getName());
						sb.append("\")");
					}
					else {
						sb.append(pathways[i].getName());
						sb.append(" : ");
						sb.append(link);
					}
				}
			}
			else sb.append("\t");
			sb.append("\n");
		}
		return sb.toString();
	}
	
	public String toStringGeneCombo(String selectNetworkId, HashMap<String, ArrayList<String>> geneSymbolKeggId) {
		//name, descLink, combinePval, combineFDR, pval, adjPval, GenesIntersect/FoundUniPathGenes, IntersectingGenes
		
		StringBuilder sb = new StringBuilder(selectNetworkId);
		sb.append("\t");
		
		//find network
		for (KeggApiNetwork n: analizedKeggApiNetworks) {
			if (selectNetworkId.equals(n.getNetworkId()) == false) continue;

			sb.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			sb.append(n.getNetworkId());		
			sb.append("\",\"");
			sb.append(n.getNetworkId());
			sb.append(": ");
			sb.append(n.getNetworkName());	
			sb.append("\")");
			sb.append("\t");
			
			//combine pval, adjPval
			sb.append(Num.formatNumber(combinePValue, 5)); sb.append("\t");
			sb.append(Num.formatNumber(combineFdr, 5)); sb.append("\t");

			//gene pval, adjPval
			sb.append(pValue); sb.append("\t");
			sb.append(Num.formatNumber(pValue, 5)); sb.append("\t");

			//GenesIntersect/FoundUniPathGenes
			double frac = (double)geneSharedGenes.size()/ (double)geneNetworkGeneNamesInInterrogatedGenes.size();
			sb.append(Num.formatNumber(frac, 5)); sb.append("\t");

			//IntersectingGenes
			//sort the sharedGenes by gene symbol so easier to group
			SelectGene[] sortedSGs = new SelectGene[geneSharedGenes.size()];
			geneSharedGenes.toArray(sortedSGs);
			Arrays.sort(sortedSGs);
			//just gene symbols
			for (SelectGene sg: sortedSGs) {
				sb.append(sg.getGeneSymbol());
				sb.append(" ");
			}
			sb.append("\t"); 
			break;
		}
			
		return sb.toString();
	}


	public TreeSet<String> getVarinatGeneNameHits() {
		TreeSet<String> allSymbols = new TreeSet<String>();
		allSymbols.addAll(aGenes.keySet());
		allSymbols.addAll(bGenes.keySet());
		return allSymbols;
	}
	
	public String[] fetchNetworkIds() {
		String[] ids = new String[analizedKeggApiNetworks.size()];
		for (int i=0; i< ids.length; i++) ids[i] = analizedKeggApiNetworks.get(i).getNetworkId();
		return ids;
	}
	
	public static String fetchKeggPathwayMapLink(String pathwayId, String[] networkIds, SelectGene[] genes, 
			HashMap<String, ArrayList<String>> geneSymbolKeggId, String negColor, String posColor, String zeroColor) {
		//for each pathway(hsa05220) in the network(N00002), build a URL
		//	color the intersected genes pink for upregulated, skyblue for downregulated
		//  add selector for all of the significant networks (adjPVal <= 0.1)
		//e.g. https://www.kegg.jp/kegg-bin/show_pathway?map=hsa05220&multi_query=673%20skyblue%0A3265%20pink&network=N00002+N00115
		//  %20 between gene and bkg color
		//  %0A between genes
		//  %23 in place of #
		//see https://www.kegg.jp/kegg/webapp/color_url.html and use hex colors
		/*
		
		# Diff exp
		Decrease Neg light blue  - Increase Pos
		blue 						red
		#87CEEB 					#FFB6C1

		# Variant
		Decrease Neg 	Increase Pos 	No change
		light violet	light orange	lightyellow
		#D6B4FC 		#FFD580 		#FFFFAC
		
		*/
		
		StringBuilder sb = new StringBuilder("https://www.kegg.jp/kegg-bin/show_pathway?map=");
		sb.append(pathwayId); 
		sb.append("&multi_query=");
		//for each gene
		
		int lastIndex = genes.length-1;
		for (int i=0; i< genes.length; i++) {
			//just take first id, very rarely there is more than one
			String geneId = geneSymbolKeggId.get(genes[i].getGeneSymbol()).get(0);
			sb.append(geneId);
			sb.append("%20");
			if (genes[i].getLog2Rto().startsWith("-")) sb.append(negColor);
			else if (genes[i].getLog2Rto().equals("0")) sb.append(zeroColor);
			else sb.append(posColor);
			if (i < lastIndex) sb.append("%0A");
		}
		
		//for each network, it's OK to include networks that are not in the pathway
		sb.append("&network=");
		sb.append(networkIds[0]);
		for (int i=1; i< networkIds.length; i++) {
			sb.append("+");
			sb.append(networkIds[i]);
		}
		return sb.toString();
	}

	public static String fetchCombineKeggPathwayMapLink(String pathwayId, String[] networkIds, SelectGene[] geneSelectGenes, 
			SelectGene[] variantSelectGenes, HashMap<String, ArrayList<String>> geneSymbolKeggId) {
		
		StringBuilder sb = new StringBuilder("https://www.kegg.jp/kegg-bin/show_pathway?map=");
		sb.append(pathwayId); 
		sb.append("&multi_query=");
		
		//find those specific to genes, variants, and in common
		HashMap<String,SelectGene> genesHash = new HashMap<String,SelectGene>();
		HashMap<String,SelectGene> variantsHash = new HashMap<String,SelectGene>();
		for (SelectGene s: geneSelectGenes) genesHash.put(s.getGeneSymbol(), s);
		for (SelectGene s: variantSelectGenes) variantsHash.put(s.getGeneSymbol(), s);
		TreeSet<String> allGeneSymbols = new TreeSet<String>();
		allGeneSymbols.addAll(genesHash.keySet());
		allGeneSymbols.addAll(variantsHash.keySet());
		
		//for each gene
		int lastIndex = allGeneSymbols.size()-1;
		int j = 0;
		for (String gs: allGeneSymbols) {
			String geneId = geneSymbolKeggId.get(gs).get(0);
			sb.append(geneId);
			sb.append("%20");
			
			//gene result present?
			SelectGene gene = genesHash.get(gs);
			String geneColor = null;
			if (gene!=null) {
				if (gene.getLog2Rto().startsWith("-")) geneColor = geneNegColor;
				else if (gene.getLog2Rto().equals("0")) geneColor = geneZeroColor; //should never happen
				else geneColor = genePosColor;
			}
			
			//variant result present?
			SelectGene var = variantsHash.get(gs);
			String varColor = null;
			if (var!=null) {
				if (var.getLog2Rto().startsWith("-")) varColor = variantNegColor;
				else if (var.getLog2Rto().equals("0")) varColor = variantZeroColor;
				else varColor = variantPosColor;
			}
			
			//both present
			if (geneColor!=null && varColor!=null) {
				sb.append(geneColor);
				sb.append("%20");
				sb.append(varColor);
			}
			else if (geneColor!=null) sb.append(geneColor);
			else sb.append(varColor);
			
			//if not the last
			if (j < lastIndex) sb.append("%0A");
			j++;
		}

		//for each network, it's OK to include networks that are not in the pathway
		sb.append("&network=");
		sb.append(networkIds[0]);
		for (int i=1; i< networkIds.length; i++) {
			sb.append("+");
			sb.append(networkIds[i]);
		}
		return sb.toString();
	}
	
	public KeggApiNetwork getKeggApiNetwork() {
		return keggApiNetwork;
	}
	public void setKeggApiNetwork(KeggApiNetwork keggApiNetwork) {
		this.keggApiNetwork = keggApiNetwork;
	}
	public ArrayList<SelectGene> getGeneSharedGenes() {
		return geneSharedGenes;
	}
	public TreeSet<String> getGeneNetworkGeneNamesInInterrogatedGenes() {
		return geneNetworkGeneNamesInInterrogatedGenes;
	}
	public boolean isGeneAnalyzeIt() {
		return geneAnalyzed;
	}
	public void setGeneAnalyzeIt(boolean geneAnalyzeIt) {
		this.geneAnalyzed = geneAnalyzeIt;
	}

	/*Sorts on sortValue*/
	public int compareTo(AnalyzedNetwork o) {
		if (this.sortValue< o.sortValue) return -1;
		if (this.sortValue> o.sortValue) return 1;
		return 0;
	}

	public String getSharedGeneSymbols() {
		StringBuilder sb = new StringBuilder(geneSharedGenes.get(0).getGeneSymbol());
		for (int i=1; i< geneSharedGenes.size(); i++) {
			sb.append(" ");
			sb.append(geneSharedGenes.get(i).getGeneSymbol());
		}
		return sb.toString();
	}
	
	public String getAnalizedNetworkIdNames() {
		StringBuilder sb = new StringBuilder(analizedKeggApiNetworks.get(0).getNetworkIdName(" : "));
		for (int i=1; i< analizedKeggApiNetworks.size(); i++) {
			sb.append(" ");
			sb.append(analizedKeggApiNetworks.get(i).getNetworkIdName(" : "));
		}
		return sb.toString();
		
	}
	
	public String getAnalizedNetworkIds() {
		StringBuilder sb = new StringBuilder(analizedKeggApiNetworks.get(0).getNetworkId());
		for (int i=1; i< analizedKeggApiNetworks.size(); i++) {
			sb.append(",");
			sb.append(analizedKeggApiNetworks.get(i).getNetworkId());
		}
		return sb.toString();
		
	}
	
	public boolean isGeneAnalyzedNetwork() {
		return totalNumDiffExpGenes > 0;
	}
	public double getPValue() {
		return pValue;
	}
	public void setPValue(double pValue) {
		this.pValue = pValue;
	}
	public double getFdr() {
		return fdr;
	}
	public void setFdr(double fdr) {
		this.fdr = fdr;
	}
	public ArrayList<String> getGroupAHits() {
		return groupAHits;
	}
	public ArrayList<String> getGroupANoHits() {
		return groupANoHits;
	}
	public ArrayList<String> getGroupBHits() {
		return groupBHits;
	}
	public ArrayList<String> getGroupBNoHits() {
		return groupBNoHits;
	}
	public HashMap<String, Integer> getaGenes() {
		return aGenes;
	}
	public HashMap<String, Integer> getbGenes() {
		return bGenes;
	}
	public void setMinKeggFreq(double minKeggFreq) {
		this.minKeggFreq = minKeggFreq;
	}
	public ArrayList<KeggApiNetwork> getAnalizedKeggApiNetworks() {
		return analizedKeggApiNetworks;
	}
	public double getSortValue() {
		return sortValue;
	}
	public void setSortValue(double sortValue) {
		this.sortValue = sortValue;
	}
	public double getVariantLog2Rto() {
		return variantLog2Rto;
	}
	public void setCombinePValue(double combinePValue) {
		this.combinePValue = combinePValue;
		
	}
	public void setCombineFdr(double combineFdr) {
		this.combineFdr = combineFdr;
	}
	public double getCombinePValue() {
		return combinePValue;
	}
	public double getCombineFdr() {
		return combineFdr;
	}
}
