package edu.utah.seq.cnv;

import java.util.ArrayList;

import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class GatkSegment implements Comparable<GatkSegment>{
	
	//fields
	private String chr;
	private int start;
	private int end;
	private int numPoints;
	private float geoMeanCopyRto;
	private String call;
	private ArrayList<Float> tumorCopyRatios = new ArrayList<Float>();
	private ArrayList<Float> normalCopyRatios = new ArrayList<Float>();
	private ArrayList<Float> tumorAlleleFrequencies = new ArrayList<Float>();
	private ArrayList<Float> normalAlleleFrequencies = new ArrayList<Float>();
	private double logMeanTumorCopyRatios = 0;
	private double logMeanNormalCopyRatios = 0;
	private double meanTumorAlleleFrequencies = 0;
	private double meanNormalAlleleFrequencies = 0;
	private double lgMeanTNRatios = 0;
	private double tnRatioAFMeans = 0;
	private Bed bedRegion;
	
	public GatkSegment (String line){
		//CONTIG	START	END	NUM_POINTS_COPY_RATIO	MEAN_LOG2_COPY_RATIO	CALL
		//chr1	14016	16502273	1747	-0.009973	0
		//chr4	58936	49197716	2290	-0.243664	-
		//chr4	49488068	49560163	6	0.923619	+
		String[] tokens = Misc.TAB.split(line);
		chr = tokens[0];
		start = Integer.parseInt(tokens[1]);
		end = Integer.parseInt(tokens[2]);
		numPoints = Integer.parseInt(tokens[3]);
		geoMeanCopyRto = Float.parseFloat(tokens[4]);
		call = tokens[5];
	}
	
	public void calculateMeans(){
		float[] tumorCRs = Num.antiLog(Num.arrayListOfFloatToArray(tumorCopyRatios), 2);
		double tumorCRMean = Num.mean(tumorCRs);
		if (tumorCRMean != 0) logMeanTumorCopyRatios = Num.log2(tumorCRMean);
		
		float[] normalCRs = Num.antiLog(Num.arrayListOfFloatToArray(normalCopyRatios), 2);
		double normalCRMean = Num.mean(normalCRs);
		if (normalCRMean != 0) logMeanNormalCopyRatios = Num.log2(normalCRMean);
		
		double[] tnRatios = Num.ratio(tumorCRs, normalCRs);
		lgMeanTNRatios = Num.log2(Num.mean(tnRatios));
		//could do a t-test here if >2 obs...
		
		//probably should match the AFs and calc individual ratios then the mean of the ratios.  Prob is these aren't always present in each set.
		if (tumorAlleleFrequencies.size()!=0){
			float[] floats = Num.arrayListOfFloatToArray(tumorAlleleFrequencies);
			meanTumorAlleleFrequencies = Num.mean(floats);
		}
		if (normalAlleleFrequencies.size()!=0){
			float[] floats = Num.arrayListOfFloatToArray(normalAlleleFrequencies);
			meanNormalAlleleFrequencies = Num.mean(floats);
		}
		if (meanNormalAlleleFrequencies !=0) tnRatioAFMeans = meanTumorAlleleFrequencies/meanNormalAlleleFrequencies;
		
		/*if (start == 46887794){
			this.printData();
			Misc.printArray(tumorCopyRatios);
			IO.pl();
			Misc.printArray(normalCopyRatios);
			IO.pl();
			Misc.printArray(tumorAlleleFrequencies);
			IO.pl();
			Misc.printArray(normalAlleleFrequencies);
			System.exit(0);
		}*/
	}
	
	public void printData(){
		IO.pl(chr+":"+start+"-"+end);
		IO.pl("\t"+numPoints +" "+geoMeanCopyRto+" "+call);
		IO.pl("\tTCR: "+tumorCopyRatios.size()+" "+logMeanTumorCopyRatios);
		IO.pl("\tNCR: "+normalCopyRatios.size()+" "+logMeanNormalCopyRatios);
		IO.pl("\tLgMeanTNRatioCR: "+lgMeanTNRatios);
		IO.pl("\tTAF: "+tumorAlleleFrequencies.size()+" "+meanTumorAlleleFrequencies);
		IO.pl("\tNAF: "+normalAlleleFrequencies.size()+" "+meanNormalAlleleFrequencies);
		IO.pl("\tRatioAF: "+tnRatioAFMeans);
	}

	public String toSpreadSheet(boolean pass){
		StringBuilder sb = new StringBuilder();
		//create IGV link  =HYPERLINK('http://localhost:60151/goto?locus=17:7379035-7395388','1')
		sb.append("=HYPERLINK(\"http://localhost:60151/goto?locus=");
		sb.append(chr); 
		sb.append(":");
		int s = start-10000;
		if (s < 0) s = 0;
		sb.append(s); sb.append("-");
		sb.append(end+10000); sb.append("\",\"");
		sb.append(chr); sb.append(":");
		sb.append(start); sb.append("-");
		sb.append(end); sb.append("\")\t");
		sb.append(pass); sb.append("\t");
		sb.append(call); sb.append("\t");
		sb.append(numPoints); sb.append("\t");
		sb.append(f(geoMeanCopyRto)); sb.append("\t");
		sb.append(f(logMeanTumorCopyRatios)); sb.append("\t");
		sb.append(f(logMeanNormalCopyRatios)); sb.append("\t");
		sb.append(f(lgMeanTNRatios)); sb.append("\t");
		sb.append(tumorAlleleFrequencies.size()); sb.append("\t");
		sb.append(f(meanTumorAlleleFrequencies));sb.append("\t");
		sb.append(normalAlleleFrequencies.size()); sb.append("\t");
		sb.append(f(meanNormalAlleleFrequencies));sb.append("\t");
		sb.append(f(tnRatioAFMeans)); 
		sb.append(bedRegion.getName());
		return sb.toString();
	}
	
	public String f(double num){
		return Num.formatNumber(num, 4);
	}
	
	
	public String toSeg(String sampleName){
		//Seg info: Sample Chr Start End
		StringBuilder sb = new StringBuilder(sampleName); sb.append("\t");
		sb.append(chr); sb.append("\t");
		sb.append(start); sb.append("\t");
		sb.append(end); sb.append("\t");
		
		//CR info: CR_NumPoints CR_TumorMean CR_NormalMean CR_TNRatio
		sb.append(numPoints); sb.append("\t");
		sb.append(f( logMeanTumorCopyRatios )); sb.append("\t");
		sb.append(f( logMeanNormalCopyRatios )); sb.append("\t");
		sb.append(f( lgMeanTNRatios)); sb.append("\t");
		
		//AF info: AF_TumorSize AF_TumorMean AF_NormalSize AF_NormalMean AF_TNRatio
		sb.append(tumorAlleleFrequencies.size()); sb.append("\t");
		sb.append(f( meanTumorAlleleFrequencies));sb.append("\t");
		sb.append(normalAlleleFrequencies.size()); sb.append("\t");
		sb.append(f( meanNormalAlleleFrequencies)); sb.append("\t");
		sb.append(f( tnRatioAFMeans)); sb.append("\t");
		
		//GATK info: GATK_Call GATK_TumorGeoMean
		sb.append(call); sb.append("\t");
		sb.append(f( geoMeanCopyRto));
		return sb.toString();
	}
	
	public boolean copyRatiosOK(){
		if (tumorCopyRatios.size() != normalCopyRatios.size() || tumorCopyRatios.size() != numPoints) return false;
		return true;
	}

	public int compareTo(GatkSegment o) {
		if (this.chr.equals(o.chr) == false) return chr.compareTo(o.chr);
		if (this.start > o.start) return 1;
		if (this.start < o.start) return -1;
		return 0;
	}

	public String getChr() {
		return chr;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public int getNumPoints() {
		return numPoints;
	}
	public float getMeanCopyRto() {
		return geoMeanCopyRto;
	}
	public String getCall() {
		return call;
	}
	public ArrayList<Float> getTumorCopyRatios() {
		return tumorCopyRatios;
	}
	public ArrayList<Float> getNormalCopyRatios() {
		return normalCopyRatios;
	}
	public ArrayList<Float> getTumorAlleleFrequencies() {
		return tumorAlleleFrequencies;
	}
	public ArrayList<Float> getNormalAlleleFrequencies() {
		return normalAlleleFrequencies;
	}

	public float getGeoMeanCopyRto() {
		return geoMeanCopyRto;
	}

	public double getLogMeanTumorCopyRatios() {
		return logMeanTumorCopyRatios;
	}

	public double getLogMeanNormalCopyRatios() {
		return logMeanNormalCopyRatios;
	}

	public double getMeanTumorAlleleFrequencies() {
		return meanTumorAlleleFrequencies;
	}

	public double getMeanNormalAlleleFrequencies() {
		return meanNormalAlleleFrequencies;
	}

	public double getLgMeanTNRatios() {
		return lgMeanTNRatios;
	}

	public double getTnRatioAFMeans() {
		return tnRatioAFMeans;
	}

	public Bed getBedRegion() {
		return bedRegion;
	}

	public void setBedRegion(Bed bedRegion) {
		this.bedRegion = bedRegion;
	}
}
