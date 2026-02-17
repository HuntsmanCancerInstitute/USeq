package edu.utah.seq.maf;

import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeSet;

import util.gen.IntersectListsHypergeometric;
import util.gen.Misc;

public class CBioProjStatRes {
	
	//fields
	private String oncoTreeTumorCode = null; // SKCM
	private String cBioPortalStudyIdentifier = null; // skcm_tcga_pan_can_atlas_2018
	
	private int numberSamples = -1;
	
	//MAP2K1 and MAP2K2
	private ArrayList<String[]> mafMAP2KLines = null;
	private int numberSamplesMAP2K = -1;
	private TreeSet<String> observedMAP2KpDots = null;
	
	//RAS
	private ArrayList<RasStats> rasStats = new ArrayList<RasStats>();
	
	//constructor
	public CBioProjStatRes(String fileName) {
		String[] n = fileName.split("--");
		oncoTreeTumorCode = n[0].toUpperCase();
		cBioPortalStudyIdentifier = n[1];
	}
	

	public void calculateIntersectionPValues(IntersectListsHypergeometric ih) throws IOException {
		/*make array of nabt for each ras
		n = total number of samples
		a = number of samples in set A
		b = number of samples in set B
		t = number of samples in common to A and B
		*/
		int[][] nabt = new int[rasStats.size()*2][4];
		int counter = 0;
		for (int i=0; i< rasStats.size(); i++) {
			RasStats r = rasStats.get(i);
			//canonical ras
			int[] cNABT = {numberSamples, numberSamplesMAP2K, r.getNumberSamplesRasCanonical(), r.getNumberSamplesMAP2K_RasCanonical()}; 
			int[] ncNABT = {numberSamples, numberSamplesMAP2K, r.getNumberSamplesRasNonCanonical(), r.getNumberSamplesMAP2K_RasNonCanonical()};
			nabt[counter++] =cNABT;
			nabt[counter++] =ncNABT;
		}
		double[] overUnderPVals = ih.calculateOverUnderRepresentationPValues(nabt);
		
		//set values into each ras stat
		counter = 0;
		for (int i=0; i< rasStats.size(); i++) {
			RasStats r = rasStats.get(i);
			double canOvPval = overUnderPVals[counter++];
			double canUnPval = overUnderPVals[counter++];
			double nonCanOvPval = overUnderPVals[counter++];
			double nonCanUnPval = overUnderPVals[counter++];
			r.setPValues(canOvPval, canUnPval, nonCanOvPval, nonCanUnPval);
		}
	}
	
	public String getHeader() {
		StringBuilder sb = new StringBuilder();
		sb.append("OT Code"); sb.append("\t");
		sb.append("cBio Study ID"); sb.append("\t");
		sb.append("Num Samples"); sb.append("\t");
		sb.append("Num Samples MAP2K Muts"); sb.append("\t");
		sb.append("Obs MAP2K pDots"); 
		for (RasStats rs: rasStats) {
			sb.append("\t");
			sb.append(rs.getHeader());
		}
		return sb.toString();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(oncoTreeTumorCode); sb.append("\t");
		sb.append(cBioPortalStudyIdentifier); sb.append("\t");
		sb.append(numberSamples); sb.append("\t");
		sb.append(numberSamplesMAP2K); sb.append("\t");
		sb.append(Misc.treeSetToString(observedMAP2KpDots, ", "));
		
		for (RasStats rs: rasStats) {
			 sb.append("\t");
			 rs.appendInfo(sb, numberSamples);
		}
		return sb.toString();
	}

	public int getNumberSamples() {
		return numberSamples;
	}

	public void setNumberSamples(int numberSamples) {
		this.numberSamples = numberSamples;
	}

	public ArrayList<String[]> getMafMAP2KLines() {
		return mafMAP2KLines;
	}

	public void setMafMAP2KLines(ArrayList<String[]> mafMAP2KLines) {
		this.mafMAP2KLines = mafMAP2KLines;
	}

	public int getNumberSamplesMAP2K() {
		return numberSamplesMAP2K;
	}

	public void setNumberSamplesMAP2K(int numberSamplesMAP2K) {
		this.numberSamplesMAP2K = numberSamplesMAP2K;
	}

	public TreeSet<String> getObservedMAP2KpDots() {
		return observedMAP2KpDots;
	}

	public void setObservedMAP2KpDots(TreeSet<String> observedMAP2KpDots) {
		this.observedMAP2KpDots = observedMAP2KpDots;
	}

	public String getOncoTreeTumorCode() {
		return oncoTreeTumorCode;
	}

	public String getcBioPortalStudyIdentifier() {
		return cBioPortalStudyIdentifier;
	}

	public ArrayList<RasStats> getRasStats() {
		return rasStats;
	}
	
}
