package edu.utah.seq.cnv.wgs;

import java.io.File;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class CopyAnalysisBedMerger {

	public TreeMap<String, Patient> patients = new TreeMap<String, Patient>();
	
	public String[] panelGenes = {"AKT1", "APC", "AR", "ARID1A", "ATM", "BARD1", "BRAF", "BRCA1", "BRCA2", "BRIP1", "CCND1", "CDK12", "CDK4", "CDK6", "CDKN1B", "CDKN2A", "CHD1", "CHEK1", "CHEK2", "COL22A1", "CTNNB1", "FANCL", "FOXA1", "FOXP1", "HOXB13", "KMT2C", "KMT2D", "MDM2", "MLH1", "MSH2", "MSH6", "MYC", "NBN", "NCOA2", "NCOR1", "NKX3-1", "NOTCH1", "PALB2", "PIK3CA", "PIK3CB", "PIK3R1", "PMS2", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "RB1", "RNF43", "SPOP", "TP53", "ZBTB16", "ZFHX3"};
	
	public CopyAnalysisBedMerger(String bedDirString) {
		
		//parse bed files
		File bedDir = new File(bedDirString);
		for (File bedFile: IO.extractFiles(bedDir, ".bed")) {
			//126006-01-001.C4D15_25KB_Hg38.called.seg.pass.bed
			String[] underscore = Misc.UNDERSCORE.split(bedFile.getName());
			String window = underscore[1];
			String [] period = Misc.PERIOD.split(underscore[0]);
			String patientId = period[0];
			String type = period[1];
			
			Patient p = patients.get(patientId);
			if (p == null) {
				p = new Patient(patientId);
				patients.put(patientId, p);
			}
			p.addCallSet(bedFile, type);
		}
		
		//printAllCopyGenes();
		
		//printJustPanelCopyAlteredGenes();
		
		printPanelGeneStats();
		
		IO.pl("Done");
		
	}

	private void printPanelGeneStats() {
		TreeMap<String, Integer> panelGeneCount = new TreeMap<String, Integer>();
		
		//for each Patient
		for (String patientId: patients.keySet()) {
			//for each call set
			TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
			for (MergedCallSet mcs: mergedCalls.values()) {
				String geneString = mcs.getPanelGeneCalls();
				for (String geneStrand: Misc.WHITESPACE.split(geneString)) {
					Integer i = panelGeneCount.get(geneStrand);
					if (i==null) panelGeneCount.put(geneStrand, 1);
					else panelGeneCount.put(geneStrand, 1+i);
				}
			}
		}
		
		for (String geneStrand: panelGeneCount.keySet()) {
			IO.pl(geneStrand+"\t"+panelGeneCount.get(geneStrand));
		}

		
	}

	private void printJustPanelCopyAlteredGenes() {
		//Print results
		//for each Patient
		for (String patientId: patients.keySet()) {
			IO.pl(patientId);
			//for each call set
			TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
			for (String type: mergedCalls.keySet()) {
				int totalNumCalls = mergedCalls.get(type).geneStrandCalls.size();
				IO.pl("\t"+type+"\t"+totalNumCalls+"\t"+ mergedCalls.get(type).getPanelGeneCalls());
			}
		}
		
	}

	private void printAllCopyGenes() {
		//Print results
		//for each Patient
		for (String patientId: patients.keySet()) {
			IO.pl(patientId);
			//for each call set
			TreeMap<String, MergedCallSet> mergedCalls = patients.get(patientId).mergedCalls;
			for (String type: mergedCalls.keySet()) {
				IO.pl("\t"+type+"\t"+mergedCalls.get(type).geneStrandCalls);
			}
			IO.pl();
		}
		
	}

	public static void main(String[] args) {
		new CopyAnalysisBedMerger("/Users/u0028003/HCI/Labs/Maughan_Benjamin/ProstateCancerProject2025/26250R/WGS/ByPatient/PassingBeds/ReallyFixed");

	}
	
	
	private class Patient {
		String patientId = null;
		TreeMap<String, MergedCallSet> mergedCalls = new TreeMap<String, MergedCallSet>();
		
		private Patient (String patientId) {
			this.patientId = patientId;
		}
		
		private void addCallSet(File bedFile, String type) {
			MergedCallSet csCallSet = mergedCalls.get(type);
			if (csCallSet == null) {
				csCallSet = new MergedCallSet(type);
				mergedCalls.put(type, csCallSet);
			}
			csCallSet.addCalls(bedFile);
		}
	}
	
	private class MergedCallSet {
		String type = null;
		TreeSet<String> geneStrandCalls = new TreeSet<String>();
		
		private MergedCallSet(String type) {
			this.type = type;
		}
		
		private void addCalls(File bedFile) {
			Bed[] bedRegions = Bed.parseFile(bedFile, 0, 0);
			for (Bed b: bedRegions) {
				//numOb=64;lg2Tum=0.2725;lg2Norm=-0.0285;genes=TTTY17C,TTTY17B,TTTY17A
				String genes = b.getName().split("genes=")[1];
				for (String g: Misc.COMMA.split(genes)) {
					//exclude any antisense stuff or empty call sets '.'
					if (g.contains("-AS") == false && g.contains(".") == false) geneStrandCalls.add(g+b.getStrand());
				}
				
			}
		}
		
		private String getPanelGeneCalls() {
			ArrayList<String> panelGenePanelCalls = new ArrayList();
			for (String g: panelGenes) {
				String gP = g+"+";
				String gM = g+"-";
				if (geneStrandCalls.contains(gP)) panelGenePanelCalls.add(gP);
				if (geneStrandCalls.contains(gM)) panelGenePanelCalls.add(gM);
			}
			return Misc.stringArrayListToString(panelGenePanelCalls, " ");
		}
	}

}
