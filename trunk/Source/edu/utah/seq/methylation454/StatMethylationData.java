package edu.utah.seq.methylation454;
import java.io.*;

import util.bio.calc.Alignment;
import util.gen.*;
import java.util.*;
import flanagan.analysis.Stat;
import trans.main.*;

public class StatMethylationData {
	//fields
	private int minimumNumberReads = 50;
	private String samplesName;
	private String ampliconName;
	private String cpGName;
	
	
	public StatMethylationData(String[] args){
		File file = new File ("/Users/nix/HCI/PIs/Yu/MethylSeqData/20071005.SV10855.Utah.Yu.amp/Amplicons/Results/Align/parsedCpGDataI85I87.txt");
		//stat samples
		Sample[] samples = new MethylationParser(file).getSamples();
		System.out.println("#Interbase coordinates, UCSC build hg18");
		System.out.println("#Minimum number reads per CpG"+minimumNumberReads);
		System.out.println("#Samples\tAmplicon\tCpG\t#ReadsOne\t#ReadsTwo\t%COne\t%CTwo\tLog2(Ratio)\t-10Log10(Chi-SqrPVal)\tLog2(AmpRatio)\t-10Log10(AmpWilcoxPVal)");
		statCpGs(samples[0],samples[1]);
	}
	
	public static void main (String[] args){
		new StatMethylationData (args);
	}
	
	/**Generates two tailed Chi-square p-values for differences in %C (a proxy for %methylation) for each CpG in each Amplicon between two samples
	 * and a Wilcoxon Signed Rank p-value for all of the CpG %Cs in an amplicon.*/
	public void statCpGs (Sample one, Sample two){
		samplesName = one.getId()+"Vs"+two.getId();
		//find common set of amplicons, paired one, two
		ArrayList commonAmplicons = findCommonAmplicons(one, two);
		//for each amplicon pair find common set of CpGs
		for (int i=0; i< commonAmplicons.size(); i++){
			Amplicon oneAmp = (Amplicon)commonAmplicons.get(i);
			i++;
			Amplicon twoAmp = (Amplicon)commonAmplicons.get(i);
			ampliconName = oneAmp.getId();
			ArrayList commonCpGs = findCommonCpGs(oneAmp, twoAmp);
			ArrayList oneCs = new ArrayList();
			ArrayList twoCs = new ArrayList();
			//perform chi-square test on CpGs, one's %C is the expected, two is the observed
			for (int j=0; j< commonCpGs.size(); j++){
				
				CpG oneCpG = (CpG)commonCpGs.get(j++);
				CpG twoCpG = (CpG)commonCpGs.get(j);
				
				//collect terms
				cpGName = oneCpG.getId();
				double numReadsOne = oneCpG.getNumberReads();
				double numReadsTwo = twoCpG.getNumberReads();
				
				//check number of reads
				if (numReadsOne < minimumNumberReads || numReadsTwo < minimumNumberReads) continue;
				
				//calculate fractions
				double fractionCOne = oneCpG.fractionC();
				double fractionCTwo = twoCpG.fractionC();
				
				//save for wilcoxon
				oneCs.add(new Float(fractionCOne));
				twoCs.add(new Float(fractionCTwo));
				
				//format for printing
				String perCOne = Num.formatNumber(fractionCOne * 100.0, 3);
				String perCTwo = Num.formatNumber(fractionCTwo * 100.0, 3);
				String ratio; 

				//watch out for zero fractions, adjust for chi-square
				if (fractionCOne !=0 && fractionCTwo !=0) ratio = Num.formatNumber(Num.log2(fractionCOne/fractionCTwo), 3);
				else {
					if (fractionCOne == 0.0){
						double cs = oneCpG.getNumberCs() + 1.0;
						fractionCOne = cs/numReadsOne;
					}
					if (fractionCTwo == 0.0){
						double cs = twoCpG.getNumberCs() + 1.0;
						fractionCTwo = cs/numReadsTwo;
					}
					ratio = Num.formatNumber(Num.log2(fractionCOne/fractionCTwo), 3)+"*";
				}

				//chi-square stuff
				//assume that the fraction C One is the Expected fraction for Two
				//what's the number of expected reads to be C in two?
				double expectedNumberCsInTwo = fractionCOne * numReadsTwo;
				double observedNumberCsInTwo = twoCpG.getNumberCs();
				double expectedNumberTsInTwo = numReadsTwo - expectedNumberCsInTwo;
				double observedNumberTsInTwo = numReadsTwo - observedNumberCsInTwo;
				double x = Stat.chiSquareFreq(new double[]{observedNumberCsInTwo, observedNumberTsInTwo}, new double[]{expectedNumberCsInTwo, expectedNumberTsInTwo});
				//two tailed!
				double pval = Alignment.transform10Log10(1- Stat.chiSquareCDF(x, 1));
				if (Double.isInfinite(pval) || pval > 150) pval = 150;
				String pvalString = Num.formatNumber(pval, 1);
				System.out.println(samplesName+"\t"+ampliconName+"\t"+cpGName+"\t"+numReadsOne+"\t"+numReadsTwo+"\t"+perCOne+"\t"+perCTwo+"\t"+ratio+"\t"+pvalString);
			}
			//calculate and print wilcoxon signed rank test for whole amplicon
			if (oneCs.size() > 4){
				float[] oneCsFloat = Num.arrayListOfFloatToArray(oneCs);
				float[] twoCsFloat = Num.arrayListOfFloatToArray(twoCs);
				String ratioAmp = Num.formatNumber(Num.log2(Num.mean(oneCsFloat)/Num.mean(twoCsFloat)), 3);
				double wilcoxonPVal = calculateWilcoxon(oneCsFloat, twoCsFloat);
				if (Double.isInfinite(wilcoxonPVal) || wilcoxonPVal > 150) wilcoxonPVal = 150;
				String pvalString = Num.formatNumber(wilcoxonPVal, 1);
				System.out.print(samplesName+"\t"+ampliconName+"\t\t\t\t\t\t\t\t");
				if (wilcoxonPVal == -1) System.out.print(ratioAmp+"\tNA");
				else System.out.println(ratioAmp+"\t"+pvalString);
			}
		}
		
	}
	
	public double calculateWilcoxon(float[] one, float[] two){
		if (one.length > 4){
			WilcoxonSignedRankTest wst = new WilcoxonSignedRankTest (one, two);
			return wst.getTransformedPValue();
		}
		return -1;
		
	}
	
	public static ArrayList findCommonCpGs (Amplicon one, Amplicon two){
		ArrayList al = new ArrayList();
		CpG[] cpGOne = one.getCpGs();
		CpG[] cpGTwo = two.getCpGs();
		for (int i=0; i< cpGOne.length; i++){
			String oneName = cpGOne[i].getId();
			for (int j=0; j< cpGTwo.length; j++){
				if (oneName.equals(cpGTwo[j].getId())){
					al.add(cpGOne[i]);
					al.add(cpGTwo[j]);
					break;
				}
			}
		}
		return al;
	}
	
	public static ArrayList findCommonAmplicons (Sample one, Sample two){
		ArrayList al = new ArrayList();
		Amplicon[] ampOne = one.getAmplicons();
		Amplicon[] ampTwo = two.getAmplicons();
		for (int i=0; i< ampOne.length; i++){
			String ampOneName = ampOne[i].getId();
			for (int j=0; j< ampTwo.length; j++){
				if (ampOneName.equals(ampTwo[j].getId())){
					al.add(ampOne[i]);
					al.add(ampTwo[j]);
					break;
				}
			}
		}
		return al;
	}
}
