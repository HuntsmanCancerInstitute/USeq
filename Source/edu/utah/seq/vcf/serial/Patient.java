package edu.utah.seq.vcf.serial;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;

import trans.main.WilcoxonSignedRankTest;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class Patient {

	//fields
	private String id = null;
	private HashMap<Integer, Sample> samples = new HashMap<Integer, Sample>();
	private double minimumVaf = 0.0;
	private boolean justVariantsInA;
	private File saveDirectory = null;
	private double minimumFracChange = 0;
	private File rScript = null;
	
	public Patient(String id, double minimumVaf, boolean justVariantsInA, File saveDirectory, double minimumFracChange, File rScript) {
		this.id = id;
		this.minimumVaf = minimumVaf;
		this.justVariantsInA = justVariantsInA;
		this.saveDirectory = saveDirectory;
		this.minimumFracChange = minimumFracChange;
		this.rScript = rScript;
	}

	public String getId() {
		return id;
	}

	public HashMap<Integer, Sample> getSamples() {
		return samples;
	}

	public void compareSamples(boolean verbose) {
		Integer[] sampleNumbers = new Integer[samples.size()];
		int index = 0;
		for (Integer i: samples.keySet()) sampleNumbers[index++] = i;
		Arrays.sort(sampleNumbers);
		
		for (int a=0; a< sampleNumbers.length-1; a++) {
			Sample sampleA = samples.get(sampleNumbers[a]);
			for (int b=a+1; b< sampleNumbers.length; b++) {
				Sample sampleB = samples.get(sampleNumbers[b]);
				compareSamplePair(sampleA, sampleB, verbose);
				//compareSamplePairSaveForR (sampleA, sampleB);
			}
		}
		
	}
	
	//uses first sample as baseline then compares it to following
	public void compareSamples2(boolean verbose) {
		Integer[] sampleNumbers = new Integer[samples.size()];
		int index = 0;
		for (Integer i: samples.keySet()) sampleNumbers[index++] = i;
		Arrays.sort(sampleNumbers);


		Sample sampleA = samples.get(sampleNumbers[0]);
		for (int b=1; b< sampleNumbers.length; b++) {
			Sample sampleB = samples.get(sampleNumbers[b]);
			//compareSamplePair(sampleA, sampleB, verbose);
			//compareSamplePair2(sampleA, sampleB, verbose);
			//compareSamplePairSaveForR (sampleA, sampleB);
			printSamplePairMeanVaf(sampleA, sampleB, verbose);
		}


	}
	
	private void printSamplePairMeanVaf(Sample sampleA, Sample sampleB, boolean verbose) {
		if (verbose) IO.pl("Calculating mean VAFs "+sampleA.getNumber()+" ("+sampleA.getName()+")"+" vs "+sampleB.getNumber()+" ("+sampleB.getName()+")");
		
		//find variants to compare
		ArrayList<Double> foundVariantsA = new ArrayList<Double>();
		for (Variant v: sampleA.getVariants().values()) if (v.getVaf()>=minimumVaf) foundVariantsA.add(v.getVaf());
		if (verbose) IO.pl("\t"+foundVariantsA.size()+" A variants found and >= min VAF, set size");
		
		ArrayList<Double> foundVariantsB = new ArrayList<Double>();
		for (Variant v: sampleB.getVariants().values()) if (v.getVaf()>=minimumVaf) foundVariantsB.add(v.getVaf());
		if (verbose) IO.pl("\t"+foundVariantsB.size()+" B variants found and >= min VAF, set size");

		//meanVAF in A and meanVAF in B
		double meanVafA = Num.meanDouble(foundVariantsA);
		double meanVafB = Num.meanDouble(foundVariantsB);

		StringBuilder sb = new StringBuilder();
		sb.append(id); sb.append("\t");
		sb.append(sampleA.getName()); sb.append("\t");
		sb.append(sampleB.getName()); sb.append("\t");
		sb.append(Num.formatNumber(meanVafA, 4)); sb.append("\t");
		sb.append(Num.formatNumber(meanVafB, 4)); sb.append("\t");
		IO.pl(sb);
	}

	
	private void compareSamplePair2(Sample sampleA, Sample sampleB, boolean verbose) {
		if (verbose) IO.pl("Comparing "+sampleA.getNumber()+" ("+sampleA.getName()+")"+" vs "+sampleB.getNumber()+" ("+sampleB.getName()+")");
		
		//find variants to compare
		TreeSet<String> foundVariants = new TreeSet<String>();
		for (Variant v: sampleA.getVariants().values()) if (v.isFoundInSample() && v.getVaf()>=minimumVaf) foundVariants.add(v.getKey());
		if (justVariantsInA==false) {
			for (Variant v: sampleB.getVariants().values()) if (v.isFoundInSample() && v.getVaf()>=minimumVaf) foundVariants.add(v.getKey());
		}
		if (verbose) IO.pl("\t"+foundVariants.size()+" variants found and >= min VAF, set size");
		
		//collect pairs
		float[] aVaf = new float[foundVariants.size()];
		float[] bVaf = new float[foundVariants.size()];
		int index = 0;
		for (String key: foundVariants) {
			aVaf[index] = (float) sampleA.getVariants().get(key).getVaf();
			bVaf[index] = (float) sampleB.getVariants().get(key).getVaf();
			if (verbose) IO.pl("\t\t"+key+"\t"+aVaf[index]+"\t"+bVaf[index]);
			index++;
		}
		
		WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(aVaf, bVaf);
		int diff = w.getNumberWMinusRanks() - w.getNumberWPlusRanks();
		int sum = w.getNumberWPlusRanks() + w.getNumberWMinusRanks();

		//meanVAF in A and meanVAF in B
		double meanVafA = Num.mean(aVaf);
		double meanVafB = Num.mean(bVaf);
		double meanFracB = meanVafB/meanVafA;

		if (verbose) IO.pl("\tPVal\t"+w.getPValue()+"\tNum+\t"+w.getNumberWPlusRanks()+"\tNum-\t"+w.getNumberWMinusRanks()+"\tTrendB-A\t"+diff+"/"+sum);
		StringBuilder sb = new StringBuilder();
		sb.append(id); sb.append("\t");
		sb.append(sampleA.getName()); sb.append("\t");
		sb.append(sampleB.getName()); sb.append("\t");
		sb.append(w.getNumberWPlusRanks()); sb.append("\t");
		sb.append(w.getNumberWMinusRanks()); sb.append("\t");
		sb.append(w.getPValue()); sb.append("\t");
		sb.append(Num.formatNumber(meanVafA, 4)); sb.append("\t");
		sb.append(Num.formatNumber(meanVafB, 4)); sb.append("\t");
		sb.append(Num.formatNumber(meanFracB, 4)); sb.append("\t");
		sb.append(Num.formatNumber((meanFracB - 1.0), 4));
		IO.pl(sb);
	}


	private void compareSamplePair(Sample sampleA, Sample sampleB, boolean verbose) {
		if (verbose) IO.pl("Comparing "+sampleA.getNumber()+" ("+sampleA.getName()+")"+" vs "+sampleB.getNumber()+" ("+sampleB.getName()+")");
		
		//find variants to compare
		TreeSet<String> foundVariants = new TreeSet<String>();
		for (Variant v: sampleA.getVariants().values()) if (v.isFoundInSample() && v.getVaf()>=minimumVaf) foundVariants.add(v.getKey());
		if (justVariantsInA==false) {
			for (Variant v: sampleB.getVariants().values()) if (v.isFoundInSample() && v.getVaf()>=minimumVaf) foundVariants.add(v.getKey());
		}
		if (verbose) IO.pl("\t"+foundVariants.size()+" variants found and >= min VAF, set size");
		
		//collect pairs
		float[] aVaf = new float[foundVariants.size()];
		float[] bVaf = new float[foundVariants.size()];
		int index = 0;
		for (String key: foundVariants) {
			aVaf[index] = (float) sampleA.getVariants().get(key).getVaf();
			bVaf[index] = (float) sampleB.getVariants().get(key).getVaf();
			if (verbose) IO.pl("\t\t"+key+"\t"+aVaf[index]+"\t"+bVaf[index]);
			index++;
		}
		
		WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(aVaf, bVaf);
		int diff = w.getNumberWMinusRanks() - w.getNumberWPlusRanks();
		int sum = w.getNumberWPlusRanks() + w.getNumberWMinusRanks();
		double frac = (double)diff/(double)sum;
		if (verbose) IO.pl("\tPVal\t"+w.getPValue()+"\tNum+\t"+w.getNumberWPlusRanks()+"\tNum-\t"+w.getNumberWMinusRanks()+"\tTrendB-A\t"+diff+"/"+sum);
		StringBuilder sb = new StringBuilder();
		sb.append(id); sb.append("\t");
		sb.append(sampleA.getName()); sb.append("\t");
		sb.append(sampleB.getName()); sb.append("\t");
		sb.append(w.getNumberWPlusRanks()); sb.append("\t");
		sb.append(w.getNumberWMinusRanks()); sb.append("\t");
		sb.append(Num.formatPercentOneFraction(frac)); sb.append("\t");
		sb.append(w.getPValue());
		IO.pl(sb);
	}
	
	private void compareSamplePairSaveForR(Sample sampleA, Sample sampleB) {
		IO.pl("Comparing "+ sampleA.getName()+ " vs "+sampleB.getName());
		ArrayList<String> forR = new ArrayList<String>();
		forR.add("VariantName\tVafA\tVafB");
		
		ArrayList<String> oneLine = new ArrayList<String>();
		oneLine.add("\tOneLine");
		oneLine.add(id);
		oneLine.add(sampleA.getName());
		oneLine.add(sampleB.getName());
		
		//find variants to compare
		TreeSet<String> foundVariants = new TreeSet<String>();
		for (Variant v: sampleA.getVariants().values()) if (v.isFoundInSample() && v.getVaf()>=minimumVaf) foundVariants.add(v.getKey());
		if (justVariantsInA==false) {
			for (Variant v: sampleB.getVariants().values()) if (v.isFoundInSample() && v.getVaf()>=minimumVaf) foundVariants.add(v.getKey());
		}
		IO.pl("\t"+foundVariants.size()+" variants found and >= min VAF ("+minimumVaf+")");
		oneLine.add(""+foundVariants.size());
		
		//collect pairs
		double[] aVafAll = new double[foundVariants.size()];
		double[] bVafAll = new double[foundVariants.size()];
		ArrayList<Double> aVafPassingMin = new ArrayList<Double>();
		ArrayList<Double> bVafPassingMin = new ArrayList<Double>();
		ArrayList<Double> diffsPassingMin = new ArrayList<Double>();
		ArrayList<Double> fracChangePassingMin = new ArrayList<Double>();
		
		int index = 0;
		for (String key: foundVariants) {
			aVafAll[index] = sampleA.getVariants().get(key).getVaf();
			bVafAll[index] = sampleB.getVariants().get(key).getVaf();
			forR.add(key+"\t"+aVafAll[index]+"\t"+bVafAll[index]);
			double diff = aVafAll[index] - bVafAll[index];
			double fracChange = Num.calculateNormalizedFractionChange(aVafAll[index],bVafAll[index]);
			if (Math.abs(fracChange)>= minimumFracChange) {
				aVafPassingMin.add(aVafAll[index]);
				bVafPassingMin.add(bVafAll[index]);
				diffsPassingMin.add(diff);
				fracChangePassingMin.add(fracChange);
			}
			index++;
		}
		
		//save all the output for R
		File f = new File(saveDirectory, id+"-"+sampleA.getName()+ "-"+sampleB.getName()+".txt");
		IO.writeArrayList(forR, f);
		
		//output cmd
		File h = new File(saveDirectory, id+"-"+sampleA.getName()+ "-"+sampleB.getName()+".html");
		String cmd = "Rscript "+rScript+" "+f+" "+h+" "+sampleA.getName()+" "+sampleB.getName()+" 0.025";
		IO.pl("\t"+cmd);
		
		//calculate stats on the min diff set
		IO.pl("\t"+diffsPassingMin.size()+" variants >= min VAF frac change ("+minimumFracChange+") for statistics");
		oneLine.add(""+diffsPassingMin.size());
		
		WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(Num.arrayListOfDoubleToFloatArray(aVafPassingMin), Num.arrayListOfDoubleToFloatArray(bVafPassingMin));
		IO.pl("\t\tPVal\t"+w.getPValue()+"\tNum+\t"+w.getNumberWPlusRanks()+"\tNum-\t"+w.getNumberWMinusRanks());
		oneLine.add(""+w.getPValue());
		oneLine.add(""+w.getNumberWPlusRanks());
		oneLine.add(""+w.getNumberWMinusRanks());
		
		double meanDiff = Num.meanDouble(diffsPassingMin);
		IO.pl("\t\tMean VAF diff\t"+meanDiff);
		oneLine.add(""+meanDiff);
		
		double meanFracChange = Num.meanDouble(fracChangePassingMin);
		IO.pl("\t\tMean VAF frac change\t"+meanFracChange);
		oneLine.add(""+meanFracChange);
		
		oneLine.add(cmd);
		
		IO.pl(Misc.stringArrayListToString(oneLine, "\t"));
		
		
	}


}
