package edu.utah.kegg.pathway;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import util.gen.Num;

public class KeggPathway implements Comparable<KeggPathway>{

	private String name;
	private String networkId;
	private HashSet<String> genes;
	private ArrayList<String> groupAHits = new ArrayList<String>();
	private ArrayList<String> groupANoHits = new ArrayList<String>();
	private ArrayList<String> groupBHits = new ArrayList<String>();
	private ArrayList<String> groupBNoHits = new ArrayList<String>();
	private double pval = -1;
	private double adjPval = -1;
	private HashMap<String, Integer> aGenes = new HashMap<String, Integer>();
	private HashMap<String, Integer> bGenes = new HashMap<String, Integer>();
	private ArrayList<String[]> urlLinkDescriptions = new ArrayList<String[]>();
	private double minKeggFreq = 0;



	public KeggPathway(String[] cells, double minKeggFreq){
		name = cells[0];
		networkId = cells[1];
		genes = new HashSet<String>(cells.length-2);
		for (int i=2; i< cells.length; i++) genes.add(cells[i]);
		this.minKeggFreq = minKeggFreq;
	}

	public String toString(boolean addOne) {
		//"#NameHyperLink\tPval\tAdjPval\tAHits\tANoHits\tFracAHits\tAGeneHits\tBHits\tBNoHits\tFracBHits\tBGeneHits\tLog2(fracA/fracB)\tAllGeneHits\tPathwayMapLinks.."
		StringBuilder sb = new StringBuilder();

		sb.append("=HYPERLINK(\"https://www.kegg.jp/entry/");
		sb.append(networkId);
		sb.append("\",\"");
		sb.append(name);
		sb.append("\")");
		String networkCell = sb.toString();
		sb.append("\t");

		sb.append(Num.formatNumber(pval, 5)); sb.append("\t");

		sb.append(Num.formatNumber(adjPval, 5)); sb.append("\t");

		sb.append(groupAHits.size()); sb.append("\t");

		sb.append(groupANoHits.size()); sb.append("\t");

		double aFrac = (double)groupAHits.size()/ (double)(groupAHits.size()+groupANoHits.size());

		sb.append(Num.formatNumber(aFrac, 5)); sb.append("\t");
		sb.append(convertToSortedFrequencies(aGenes)); sb.append("\t");

		sb.append(groupBHits.size()); sb.append("\t");
		sb.append(groupBNoHits.size()); sb.append("\t");
		double bFrac = (double)groupBHits.size()/ (double)(groupBHits.size()+groupBNoHits.size());
		sb.append(Num.formatNumber(bFrac, 5)); sb.append("\t");
		sb.append(convertToSortedFrequencies(bGenes)); sb.append("\t");

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
		double rto = aFrac/bFrac;
		sb.append(Num.formatNumber(Num.log2(rto), 5));
		sb.append("\t");

		TreeSet<String> allGenes = new TreeSet<String>();
		allGenes.addAll(aGenes.keySet());
		allGenes.addAll(bGenes.keySet());
		for (String geneName: allGenes) {
			sb.append(geneName);
			sb.append(" ");
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

	public String convertToSortedFrequencies(HashMap<String,Integer> keyCounts) {
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
		return namesToView.toString();
	}

	public int compareTo(KeggPathway o) {
		if (this.pval< o.pval) return -1;
		if (this.pval> o.pval) return 1;
		return 0;
	}

	public String getNetworkId() {
		return networkId;
	}

	public void addPathwayMapDescription(String[] pd) {
		urlLinkDescriptions.add(pd);
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

	public double getPval() {
		return pval;
	}

	public void setPval(double pval) {
		this.pval = pval;
	}

	public void setAdjPval(double adjPval) {
		this.adjPval = adjPval;
	}

	public HashSet<String> getGenes() {
		return genes;
	}

	public HashMap<String, Integer> getaGenes() {
		return aGenes;
	}

	public HashMap<String, Integer> getbGenes() {
		return bGenes;
	}

}
