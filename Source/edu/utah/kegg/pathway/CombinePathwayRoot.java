package edu.utah.kegg.pathway;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.utah.kegg.api.KeggApiNetwork;
import edu.utah.kegg.api.KeggApiPathway;
import util.gen.CombinePValues;
import util.gen.Misc;
import util.gen.Num;

public class CombinePathwayRoot {

	//fields
	private AnalyzedNetwork[] geneAnalyzedNetworks = null;
	private AnalyzedNetwork[] variantAnalyzedNetworks = null;
	private HashMap<String, CombinePathway> pathwayIdCombinePathway = null;
	
	/**Constructor, either may be null*/
	public CombinePathwayRoot(AnalyzedNetwork[] geneAnalyzedNetworks, AnalyzedNetwork[] variantAnalyzedNetworks) {
		this.geneAnalyzedNetworks = geneAnalyzedNetworks;
		this.variantAnalyzedNetworks = variantAnalyzedNetworks;
	};
	
	
	//methods
	public void makeCombinePathways(double maximumFdr) {
		pathwayIdCombinePathway = new HashMap<String, CombinePathway>();
		
		//merge all the AnalyzedNetworks
		ArrayList<AnalyzedNetwork> all = new ArrayList<AnalyzedNetwork>();
		if (geneAnalyzedNetworks!=null) for (AnalyzedNetwork an: geneAnalyzedNetworks) all.add(an);
		if (variantAnalyzedNetworks!=null)  for (AnalyzedNetwork an: variantAnalyzedNetworks) all.add(an);

		//for each Analyzed Network
		for (AnalyzedNetwork an: all) {
			
			// does it pass the FDR?
			if (an.getFdr()<= maximumFdr) {
				
				//for each kegg api network in the an
				for (KeggApiNetwork kan: an.getAnalizedKeggApiNetworks()){
					
					//any associated pathways?
					if (kan.getPathways() != null) {
						//for each kegg pathway, add the Analyzed Network
						for (KeggApiPathway kap: kan.getPathways()) {
							CombinePathway cp = pathwayIdCombinePathway.get(kap.getId());
							if (cp==null) {
								cp = new CombinePathway(kap);
								pathwayIdCombinePathway.put(kap.getId(), cp);
							}
							//must watch for duplicates
							String networkIds = an.getAnalizedNetworkIds();
							cp.getNetworkIdsAnalyzedNetworks().put(networkIds, an);
						}
					}
				}
			}
		}
		calculateCombinePValues();
	}
	
	public void calculateCombinePValues() {
		
		//calc combine pvalues for pathways with more than one analysis
		ArrayList<Float> combinePValues = new ArrayList<Float>();
		CombinePValues combPVal = new CombinePValues();
		for (CombinePathway cp: pathwayIdCombinePathway.values()) {
			int numAns = cp.getNetworkIdsAnalyzedNetworks().size();
			if (numAns > 1) {
				double[] toCombine = new double[numAns];
				int i = 0;
				for (AnalyzedNetwork an: cp.getNetworkIdsAnalyzedNetworks().values()) toCombine[i++] = an.getPValue();
				double cPVal = combPVal.calculateCombinePValues(toCombine);
				cp.setCombinePValue(cPVal);
				combinePValues.add(Num.minus10log10Float(cPVal));
			}
		}
		
		//adjust pvalues
		float[] f = Num.benjaminiHochbergCorrectUnsorted(Num.arrayListOfFloatToArray(combinePValues));
		double[] fdrs = Num.antiNeg10log10(f);
		int i=0;
		for (CombinePathway cp: pathwayIdCombinePathway.values()) {
			if (cp.getNetworkIdsAnalyzedNetworks().size() > 1) cp.setCombineFdr(fdrs[i++]);
		}

	}
	
	public void saveGenePathways(double maximumFdr, HashMap<String, ArrayList<String>> gs2ki, File resultsDirectory, int minimumNumberGenes) throws IOException {
		//the idea here is to find all of the significant networks with the same pathway
		//	then for each pathway decorate it with all of the network genes colored, and highlight the networks
		//	really only useful when multiple diff networks hit the same pathway
		//		output pathway view xls file or json for web apps
		
		//for each pathway, create a results obj to sort by pvalue
		StringValueSort[] results = new StringValueSort[pathwayIdCombinePathway.size()];
		int counter = 0;
		for (String pathwayId: pathwayIdCombinePathway.keySet()) {
			CombinePathway cp = pathwayIdCombinePathway.get(pathwayId);
			
			StringBuilder xls = new StringBuilder();
			xls.append("PATHWAY\t");
			xls.append(pathwayId); 
			xls.append("\t"); 
			xls.append(cp.getKeggPathway().getName()); 
			xls.append("\t");  
			xls.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);		
			xls.append("\",\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);	
			xls.append("\")\n");

			//fetch all of the genes and networks
			TreeMap<String, SelectGene> allGenes = new TreeMap<String,SelectGene>();
			HashSet<String> networkIds = new HashSet<String>();
			xls.append("Networks\n");
			double minPVal = Double.MAX_VALUE;
			double fdr = Double.MAX_VALUE;
			for (AnalyzedNetwork an : cp.getNetworkIdsAnalyzedNetworks().values()) {
				if (an.isGeneAnalyzedNetwork() == false) continue;
				if (an.getPValue()< minPVal) {
					minPVal = an.getPValue();
					fdr = an.getFdr();
				}
				for (SelectGene sg: an.getGeneSharedGenes()) allGenes.put(sg.getGeneSymbol(), sg);
				for (KeggApiNetwork net: an.getAnalizedKeggApiNetworks()) {
					networkIds.add(net.getNetworkId());
					xls.append("\t"); xls.append(net.getNetworkIdName("\t")); xls.append("\n");
				}
				xls.append("\t\tAdjPval\t"); xls.append(Num.formatNumber(an.getFdr(),3)); xls.append("\n");
				xls.append("\t\tGenes\t"); xls.append(an.getSharedGeneSymbols()); xls.append("\n");
			}
			
			SelectGene[] sg = new SelectGene[allGenes.size()];
			int i = 0;
			for (SelectGene g: allGenes.values()) sg[i++] = g;
			xls.append("AllGenes\t"); xls.append(Misc.stringSetToString(allGenes.keySet(), ", ")); xls.append("\n");
			
			//networks
			String[] networkIdsStringArray = Misc.hashSetToStringArray(networkIds);
			//"skyblue" for negColor, "pink" for posColor and zeroColor
			String url = AnalyzedNetwork.fetchKeggPathwayMapLink(pathwayId, networkIdsStringArray, sg, gs2ki,AnalyzedNetwork.geneNegColor,AnalyzedNetwork.genePosColor,AnalyzedNetwork.genePosColor);
			xls.append("KeggLink\t"); xls.append(url); xls.append("\n");
			if (url.length()<250) {
				xls.append("ExcelLink\t=HYPERLINK(\"");
				xls.append(url);		
				xls.append("\",\""+url);	
				xls.append("\")\n");
			}
			else xls.append("ExcelLink\tToo big\n");
			double pathwayPValue = minPVal;
			double pathwayFdr = fdr;
			if (cp.getCombinePValue()!=-1) {
				pathwayPValue = cp.getCombinePValue();
				pathwayFdr = cp.getCombineFdr();
			}
			xls.append("PathwayPValue\t"); xls.append(pathwayPValue);
			xls.append("\nPathwayFDR\t"); xls.append(pathwayFdr);
			xls.append("\n\n");
			
			results[counter++] = new StringValueSort(xls, minPVal);
		}
		
		Arrays.sort(results);
		PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "genePathwaysMinGen"+minimumNumberGenes+"MaxFdr"+maximumFdr+".xls")));
		//name descLink pval, adjPval, #UniqueNetworkGenes, #FoundUniqueNetworkGenes, #DiffExpGenes, #GenesIntersect, GenesIntersect/PathGenes, IntersectingGenes
		out.println("#Composite view of significant networks in KEGG pathways from the gene set analysis\n");
		out.println(fetchColorKey(true, false)+"\n");
		for (StringValueSort s: results) out.print(s.getCargo().toString());
		out.close();
		
	}
	
	public static String fetchColorKey(boolean includeGenes, boolean includeVariants) {
		StringBuilder sb = new StringBuilder("Kegg Link Gene Color Key:\n");
		if (includeGenes) {
			sb.append("Negative Diff Exp LgRto\tLight Blue\t#87CEEB\n");
			sb.append("Positive Diff Exp LgRto\tLight Red, Pink\t#FFB6C1\n");
		}
		if (includeVariants) {
			sb.append("Negative Mutation Freq LgRto\tLight Violet\t#D6B4FC\n");
			sb.append("Positive Mutation Freq LgRto\tLight Orange\t#FFD580\n");
			sb.append("Abs(Mut Freq) < 1.25x\tLight Yellow\t#FFFFAC\n");
		}
		sb.append("Genes in red text are disease associated\t");
		return sb.toString();
	}


	public void saveVariantPathways(double maximumFdr, HashMap<String, ArrayList<String>> gs2ki, File resultsDirectory, int minimumNumberGenes) throws IOException {
		//for each pathway, create a results obj to sort by pvalue
		StringValueSort[] results = new StringValueSort[pathwayIdCombinePathway.size()];
		int counter = 0;
		
		for (String pathwayId: pathwayIdCombinePathway.keySet()) {
			CombinePathway cp = pathwayIdCombinePathway.get(pathwayId);
			
			StringBuilder xls = new StringBuilder();
			xls.append("PATHWAY\t");
			xls.append(pathwayId); 
			xls.append("\t"); 
			xls.append(cp.getKeggPathway().getName()); 
			xls.append("\t");  
			xls.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);		
			xls.append("\",\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);	
			xls.append("\")\n");

			//fetch all of the genes and networks
			TreeMap<String, SelectGene> allGenes = new TreeMap<String,SelectGene>();
			HashSet<String> networkIds = new HashSet<String>();
			xls.append("Networks\n");
			double minPVal = Double.MAX_VALUE;
			double fdr = Double.MAX_VALUE;
			
			// for each AN
			for (AnalyzedNetwork an : cp.getNetworkIdsAnalyzedNetworks().values()) {
				// check if this is a variant AN 
				if (an.isGeneAnalyzedNetwork() == true) continue;
				if (an.getPValue()< minPVal) {
					minPVal = an.getPValue();
					fdr = an.getFdr();
				}
				for (SelectGene sg: an.fetchVariantGenes()) allGenes.put(sg.getGeneSymbol(), sg);

				for (KeggApiNetwork net: an.getAnalizedKeggApiNetworks()) {
					networkIds.add(net.getNetworkId());
					xls.append("\t"); xls.append(net.getNetworkIdName("\t")); xls.append("\n");
				}
				xls.append("\t\tAdjPval\t"); xls.append(Num.formatNumber(an.getFdr(),3)); xls.append("\n");
				xls.append("\t\tLog2Rto\t"); xls.append(Num.formatNumber(an.getVariantLog2Rto(),3)); xls.append("\n");
				xls.append("\t\tGenes\t"); xls.append(Misc.stringSetToString(an.getVarinatGeneNameHits(), ", ")); xls.append("\n");
			}
			
			// pathway info
			SelectGene[] sg = new SelectGene[allGenes.size()];
			int i = 0;
			for (SelectGene g: allGenes.values()) sg[i++] = g;
			xls.append("AllGenes\t"); xls.append(Misc.stringSetToString(allGenes.keySet(), ", ")); xls.append("\n");

			//pathway link
			String[] networkIdsStringArray = Misc.hashSetToStringArray(networkIds);
			String url = AnalyzedNetwork.fetchKeggPathwayMapLink(pathwayId, networkIdsStringArray, sg, gs2ki, AnalyzedNetwork.variantNegColor, AnalyzedNetwork.variantPosColor,AnalyzedNetwork.variantZeroColor);
			xls.append("KeggLink\t"); xls.append(url); xls.append("\n");
			if (url.length()<250) {
				xls.append("ExcelLink\t=HYPERLINK(\"");
				xls.append(url);		
				xls.append("\",\""+url);	
				xls.append("\")\n");
			}
			else xls.append("ExcelLink\tToo big\n");
			
			double pathwayPValue = minPVal;
			double pathwayFdr = fdr;
			if (cp.getCombinePValue()!=-1) {
				pathwayPValue = cp.getCombinePValue();
				pathwayFdr = cp.getCombineFdr();
			}
			xls.append("PathwayPValue\t"); xls.append(pathwayPValue);
			xls.append("\nPathwayFDR\t"); xls.append(pathwayFdr);
			xls.append("\n\n");
			
			results[counter++] = new StringValueSort(xls, minPVal);
		}
		
		Arrays.sort(results);
		PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "variantPathwaysMinGen"+minimumNumberGenes+"MaxFdr"+maximumFdr+".xls")));
		out.println("#Composite view of significant networks in KEGG pathways from the variant set analysis\n");
		out.println(fetchColorKey(false, true)+"\n");
		for (StringValueSort s: results) out.print(s.getCargo().toString());
		out.close();
	}

	public void saveGeneAndVariantPathways(double maximumFdr, HashMap<String, ArrayList<String>> gs2ki, File resultsDirectory, int minimumNumberGenes) throws IOException {
		//for each pathway, create a results obj to sort by pvalue
		StringValueSort[] results = new StringValueSort[pathwayIdCombinePathway.size()];
		int counter = 0;
		
		for (String pathwayId: pathwayIdCombinePathway.keySet()) {
			CombinePathway cp = pathwayIdCombinePathway.get(pathwayId);
			
			StringBuilder xls = new StringBuilder();
			xls.append("PATHWAY\t");
			xls.append(pathwayId); 
			xls.append("\t"); 
			xls.append(cp.getKeggPathway().getName()); 
			xls.append("\t");  
			xls.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);		
			xls.append("\",\"https://www.kegg.jp/entry/");
			xls.append(pathwayId);	
			xls.append("\")\n"); 
			
			//fetch all of the genes and networks
			TreeMap<String, SelectGene> allGenesGenes = new TreeMap<String,SelectGene>();
			TreeMap<String, SelectGene> allGenesVariants = new TreeMap<String,SelectGene>();
			HashSet<String> networkIds = new HashSet<String>();
			xls.append("Networks\n");
			double minPVal = Double.MAX_VALUE;
			double fdr = Double.MAX_VALUE;
			
			// for each AN
			for (AnalyzedNetwork an : cp.getNetworkIdsAnalyzedNetworks().values()) {
				
				String type = "Diff Exp";
				if (an.isGeneAnalyzedNetwork()==false) type = "Diff Mut";
				
				if (an.getPValue()< minPVal) {
					minPVal = an.getPValue();
					fdr = an.getFdr();
				}
				
				for (KeggApiNetwork net: an.getAnalizedKeggApiNetworks()) {
					networkIds.add(net.getNetworkId());
					xls.append("\t"); 
					xls.append(type);
					xls.append("\t");
					xls.append(net.getNetworkIdName("\t")); xls.append("\n");
				}
				
				// check if this is a variant AN 
				if (an.isGeneAnalyzedNetwork() == true) {
					//gene
					for (SelectGene sg: an.getGeneSharedGenes()) allGenesGenes.put(sg.getGeneSymbol(), sg);
					xls.append("\t\tAdjPval\t"); xls.append(Num.formatNumber(an.getFdr(),3)); xls.append("\n");
					xls.append("\t\tGenes\t"); xls.append(an.getSharedGeneSymbols()); xls.append("\n");
					
				}
				else {
					//variant
					for (SelectGene sg: an.fetchVariantGenes()) allGenesVariants.put(sg.getGeneSymbol(), sg);
					xls.append("\t\tAdjPval\t"); xls.append(Num.formatNumber(an.getFdr(),3)); xls.append("\n");
					xls.append("\t\tLog2Rto\t"); xls.append(Num.formatNumber(an.getVariantLog2Rto(),3)); xls.append("\n");
					xls.append("\t\tGenes\t"); xls.append(Misc.stringSetToString(an.getVarinatGeneNameHits(), ", ")); xls.append("\n");
				}
			}
			
			// pathway info
			TreeSet<String> allGenes = new TreeSet<String>();
			SelectGene[] sgG = new SelectGene[allGenesGenes.size()];
			int i = 0;
			for (SelectGene g: allGenesGenes.values()) sgG[i++] = g;
			allGenes.addAll(allGenesGenes.keySet());
			SelectGene[] sgV = new SelectGene[allGenesVariants.size()];
			i = 0;
			for (SelectGene g: allGenesVariants.values()) sgV[i++] = g;
			allGenes.addAll(allGenesVariants.keySet());
			xls.append("AllGenes\t"); xls.append(Misc.stringSetToString(allGenes, ", ")); xls.append("\n");

			//pathway link
			String[] networkIdsStringArray = Misc.hashSetToStringArray(networkIds);
			String url = AnalyzedNetwork.fetchCombineKeggPathwayMapLink(pathwayId, networkIdsStringArray, sgG, sgV, gs2ki);
			xls.append("KeggLink\t"); xls.append(url); xls.append("\n");
			if (url.length()<250) {
				xls.append("ExcelLink:\t=HYPERLINK(\"");
				xls.append(url);		
				xls.append("\",\""+url);	
				xls.append("\")\n");
			}
			else xls.append("ExcelLink\tToo big\n");
			
			double pathwayPValue = minPVal;
			double pathwayFdr = fdr;
			if (cp.getCombinePValue()!=-1) {
				pathwayPValue = cp.getCombinePValue();
				pathwayFdr = cp.getCombineFdr();
			}
			xls.append("PathwayPValue\t"); xls.append(pathwayPValue);
			xls.append("\nPathwayFDR\t"); xls.append(pathwayFdr);
			xls.append("\n\n");
			
			results[counter++] = new StringValueSort(xls, minPVal);
		}
		
		Arrays.sort(results);
		PrintWriter out = new PrintWriter( new FileWriter(new File(resultsDirectory, "combineGeneVariantPathwaysMinGen"+minimumNumberGenes+"MaxFdr"+maximumFdr+".xls")));
		out.println("#Composite view of significant networks in KEGG pathways from the variant and gene set analysis\n");
		out.println(fetchColorKey(true, true)+"\n");
		for (StringValueSort s: results) out.print(s.getCargo().toString());
		out.close();
	}


}
