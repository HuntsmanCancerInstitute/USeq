package edu.expr;
import util.bio.parsers.*;
import util.gen.Num;

import java.util.*;
import java.util.regex.Pattern;

/**Container for holding alleles affecting a particular gene model.*/
public class AffectedGene {
	
	//fields
	private UCSCGeneLine geneModel;
	private Allele[] alleles;
	private String[] effects;
	public static final Pattern whiteSpace = Pattern.compile("\\s+");
	public static final Pattern numberWhiteSpace = Pattern.compile("\\d+|\\s+");
	public HashSet<String> printedTranscripts = null;
	
	//constructor
	public AffectedGene (UCSCGeneLine geneModel, Allele[] alleles, String chromSeq, HashSet<String> printedTranscripts){
		this.geneModel = geneModel;
		this.alleles = alleles;
		this.printedTranscripts = printedTranscripts;
		mutate(chromSeq);
	}

	//methods
	/**Characterizes effects of each allele on the gene model.*/
	public void mutate(String chromSeq){
		effects = new String[alleles.length];
		for (int i=0; i< alleles.length; i++){
			effects[i] = alleles[i].affect(geneModel, chromSeq);
		}
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(geneModel.simpleToString());
		sb.append("\n");
		//warn if no cds
		if (geneModel.getCdsEnd() == geneModel.getCdsStart()) {
			sb.append("\tWARNING: No coding sequence indicated!\n");
		}
		for (int i=0; i< alleles.length; i++){
			sb.append("\t");
			sb.append(alleles[i]);
			sb.append("\n");
			sb.append(effects[i]);
		}
		return sb.toString();
	}
	
	public String toStringBedLine(boolean collapse){
		//get chrom
		String chrom = geneModel.getChrom()+"\t";
		StringBuilder sb = new StringBuilder();
		for (int i=0; i< alleles.length; i++){
			//check and see if different?
			if (collapse){
				String test = geneModel.getChrom()+"_"+ alleles[i].getStart()+ "_"+ alleles[i].getStop() +new String (numberWhiteSpace.matcher(effects[i]).replaceAll(""));
				if (printedTranscripts.contains(test) == false) {
					printedTranscripts.add(test);
				}
				else continue;
			}
			//add chrom
			sb.append(chrom);
			//watch out for INDELS
			int start = alleles[i].getStart();
			int stop = alleles[i].getStop();
			if (start == stop){
				start--;
				stop++;
			}
			sb.append(start);
			sb.append("\t");
			sb.append(stop);
			sb.append("\t");
			String simple = geneModel.getNames("_") +"_|_"+ alleles[i].getNotes()+"_|_"+effects[i];
			simple = simple.trim();
			simple = whiteSpace.matcher(simple).replaceAll("_");
			sb.append(simple);
			sb.append ("\t1000\t");
			sb.append(geneModel.getStrand());
			sb.append("\n");
		}
		return sb.toString();
	}
	
	/**Just return potentially deleterious alleles ('Non-synonymous or splice') in bed file format.
	 * Returns null if not deleterious.*/
	public String toStringDeleteriousBedLine(boolean collapse){
		//scan for deleterious mutations
		ArrayList<Integer> bad = new ArrayList<Integer>();
		for (int i=0; i< effects.length; i++){
			if (effects[i].contains("Non-synonymous") || effects[i].contains("splice")){
				if (collapse){
					String test = geneModel.getChrom()+"_"+ alleles[i].getStart()+ "_"+ alleles[i].getStop() +new String (numberWhiteSpace.matcher(effects[i]).replaceAll(""));
					if (printedTranscripts.contains(test) == false) {
						bad.add(new Integer(i));
						printedTranscripts.add(test);					
					}
				}
				else bad.add(new Integer(i));
			}
		}
		if (bad.size()==0) return null;
		//get chrom
		String chrom = geneModel.getChrom()+"\t";
		StringBuilder sb = new StringBuilder();

		int[] badAlleles = Num.arrayListOfIntegerToInts(bad);
		for (int i=0; i< badAlleles.length; i++){
			
			//add chrom
			sb.append(chrom);
			//watch out for INDELS
			int start = alleles[badAlleles[i]].getStart();
			int stop = alleles[badAlleles[i]].getStop();
			if (start == stop){
				start--;
				stop++;
			}
			sb.append(start);
			sb.append("\t");
			sb.append(stop);
			sb.append("\t");
			String simple = geneModel.getNames("_") +"_|_"+ alleles[badAlleles[i]].getNotes()+"_|_"+effects[badAlleles[i]];
			simple = simple.trim();
			simple = whiteSpace.matcher(simple).replaceAll("_");
			sb.append(simple);			
			sb.append("\t1000\t");
			sb.append(geneModel.getStrand());
			sb.append("\n");
		}
		return sb.toString();
	}
	
		/**Just return potentially deleterious alleles ('Non-synonymous or splice').*/
	public String toStringDeleterious(boolean collapse){
		//scan for deleterious mutations
		ArrayList<Integer> bad = new ArrayList<Integer>();
		for (int i=0; i< effects.length; i++){
			if (effects[i].contains("Non-synonymous") || effects[i].contains("splice")){
				if (collapse){
					String test = geneModel.getChrom()+"_"+ alleles[i].getStart()+ "_"+ alleles[i].getStop() +new String (numberWhiteSpace.matcher(effects[i]).replaceAll(""));
					if (printedTranscripts.contains(test) == false) {
						bad.add(new Integer(i));
						printedTranscripts.add(test);
					}
				}
				else bad.add(new Integer(i));
			}
		}
		if (bad.size()==0) return null;
		StringBuilder sb = new StringBuilder();
		sb.append(geneModel.simpleToString());
		sb.append("\n");
		//warn if no cds
		if (geneModel.getCdsEnd() == geneModel.getCdsStart()) {
			sb.append("\n\tWARNING: No coding sequence indicated!\n\n");
		}
		int[] badAlleles = Num.arrayListOfIntegerToInts(bad);
		for (int i=0; i< badAlleles.length; i++){
			sb.append("\t");
			sb.append(alleles[badAlleles[i]]);
			sb.append("\n");
			sb.append(effects[badAlleles[i]]);
		}
		return sb.toString();
	}
	
	//getters and setters
	public UCSCGeneLine getGeneModel() {
		return geneModel;
	}
	public Allele[] getAlleles() {
		return alleles;
	}
	
	
}
