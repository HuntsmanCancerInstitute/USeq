package edu.utah.kohli;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;

import org.freehep.graphicsio.swf.SWFAction.Call;

import edu.utah.kohli.KohliGene.CopyTestResult;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class CompareGATKDatasets {
	
	private TreeMap<String, KohliPatient> patients = new TreeMap<String,KohliPatient>();
	private TreeMap<String, KohliSample> samples = new TreeMap<String,KohliSample>();
	private LinkedHashMap<String, KohliGene> genes = null;
	private File resultsDirectory = null;
	private String[] testedGenes = new String[]{"AR-Enh","ARID1A","PTEN","CCND1","ATM","ZBTB16","CDKN1B","KMT2D","CDK4","MDM2","BRCA2","RB1","FOXA1","AKT1","ZFHX3","TP53","NCOR1","CDK12","BRCA1","SPOP","RNF43","MSH2","ERG","TMPRSS2","CHEK2","CTNNB1","FOXP1","PIK3CB","PIK3CA","PIK3R1","CHD1","APC","CDK6","BRAF","KMT2C","NKX3-1","NCOA2","MYC","COL22A1","CDKN2A","NOTCH1","AR","OPHN1"};
	//private String[] testedGenes = new String[]{"AR","BRCA2","CHEK2","MYC","NKX3-1","OPHN1","PIK3CA","PIK3CB","TP53","ZBTB16"};
	private HashMap<String, String> geneNameInfo = null;
	private Bed arEnhancer = new Bed("chrX", 66899317, 66907698);
	private HashMap<String, String> vcfKeyRecord = null;
	private File fullPathToR = new File("/usr/local/bin/R");
	private File saveDir = null;

	public static void main(String[] args) throws IOException {
		File patientSampleInfo = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt");
		File geneInfo = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNADesign/CBioReviewProstateCancer23July2024/targetGeneSummary.txt");
		File gatkResultsDir = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/WgsAnalysis/GATKForWGS/SubSamplingComparision/SubSampled");
		File gatkKeyResultsDir = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/WgsAnalysis/GATKForWGS/SubSamplingComparision/S0_100_PassBed");
		new CompareGATKDatasets(patientSampleInfo, geneInfo, gatkResultsDir, gatkKeyResultsDir);
	}
	
	public CompareGATKDatasets(File patientSampleInfo, File geneInfo, File gatkResultsDir, File gatkKeyResultsDir) throws IOException {

		//load patients and samples
		loadPatients(patientSampleInfo);
		for (KohliPatient kp: patients.values()) IO.pl(kp);

		//load gene info
		loadGeneInfo(geneInfo);

		//load the key
		loadGatkCopyRatioCalls(gatkKeyResultsDir, "Key");

		//load dirs of cnv results
		File[] splitDirs = IO.extractOnlyDirectories(gatkResultsDir);
		for (File sd: splitDirs) loadGatkCopyRatioCalls(sd, sd.getName());

		compareResultsWithMergedGatkWgsConfusionMatrix();

		IO.pl("\nCOMPLETE!");

	}
	
	




	private void loadGeneInfo(File geneInfo) {
		String[] lines = IO.loadFileIntoStringArray(geneInfo);
		geneNameInfo = new HashMap<String,String>();
		for (String l: lines) {
			if (l.startsWith("#") || l.length() ==0) continue;
			String[] f = Misc.TAB.split(l);
			IO.pl(l);
			geneNameInfo.put(f[0], f[1]+" -> "+f[2]);
		}
		
	}

	private void compareResultsWithMergedGatkWgsConfusionMatrix() {
		IO.pl("\nComparing copy number results...");
		IO.pl("Header:\t"+ConfusionMatrix.toStringHeader());
		//for each patient
		for (KohliPatient kp: patients.values()) {

			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.p(kp.getHciPatientId()+"\t"+ks.getSampleId());
				//find the key
				CnvCallSet key = null;
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("Key")) {
						key = cc;
						break;
					}
				}
				if (key == null) {
					IO.pl("\tKey not found....");
					continue;
				}

				String keyGenes = key.getCopyAlteredGeneString(testedGenes);
				TreeSet<String> keyGeneCalls = mergeGeneStrings(keyGenes);
				IO.p("\tKey: "+keyGenes);
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("Key")) continue;
					IO.p("\t"+cc.getMethodName()+": "+cc.getCopyAlteredGeneString(testedGenes));
					ArrayList<String> genesDirection = cc.getCopyAlteredGenes(testedGenes);
					
					//create confusion matrix
					HashSet<String> testGeneCalls = new HashSet<String>();
					for (String s: genesDirection) testGeneCalls.add(s);
					ConfusionMatrix cm = new ConfusionMatrix(testedGenes, keyGeneCalls, testGeneCalls);
					
					IO.p("\t"+ cm.toString());
					
				}
				IO.pl();
			}	
		}
	}


	public String formatFraction(int a, int b) {
		double f = (double)a / (double)b;
		return a+"/"+b+"("+Num.formatNumber(f, 4)+")";
	}



	public TreeSet<String> mergeGeneStrings(String one) {
		TreeSet<String> merge = new TreeSet<String>();
		String[] oneSplit = Misc.COMMA.split(one);
		for (String x: oneSplit) merge.add(x);
		if (merge.size()>1 && merge.contains("None")) merge.remove("None");
		return merge;
	}

	private void loadGatkCopyRatioCalls(File gatkResultsDir, String name) {
		IO.pl("Loading GATK "+name+" results...");
		IO.pl("\tDataset\tNumBedLines\tNumInterrogatedGenesCalled");
		File[] beds = IO.extractFiles(gatkResultsDir, "called.seg.pass.bed.gz");
		Pattern genePattern = Pattern.compile(";genes=");
		//for each results bed file results set
		for (File b: beds) {
			String[] patient_cfDNA_germline = Misc.UNDERSCORE.split(b.getName().replaceAll("G", "X"));
			
			//pull KohliSample
			KohliSample ks = samples.get(patient_cfDNA_germline[1]);
			if (ks == null) Misc.printErrAndExit("\nFailed to find the sample associated with "+b);
			
			//load all of the bed results
			Bed[] bedResults = Bed.parseFile(b, 0, 0);
			int numCalls = bedResults.length;
			TreeMap<String, GeneCallResult> affectedGeneStats = new TreeMap<String, GeneCallResult>();
			
			//for each bed line
			for (Bed bed: bedResults) {
				//name    numOb=28;lg2Tum=-0.191;lg2Norm=0;genes=RALGAPA1P1,SNORD121B,LOC101928775
				String[] f = genePattern.split(bed.getName());
				if (f.length!=2 || f[0].startsWith("numOb=")==false) Misc.printErrAndExit("Failed to split on genes the bed result: "+bed.toString());
				String[] genesAffected = Misc.COMMA.split(f[1]);
				boolean containsMinus = f[0].contains("lg2Tum=-");
				HashSet<String> gas = Misc.loadHashSet(genesAffected);
				//for each interrogated gene
				for (String intGene: testedGenes) {
					if (gas.contains(intGene)) {
						affectedGeneStats.put(intGene, new GeneCallResult(true, containsMinus==false, f[0]));
					}
				}
				//does it intersect the AR Enhancer?
				if (bed.intersects(arEnhancer)) {
					affectedGeneStats.put("AR-Enh", new GeneCallResult(true, containsMinus==false, f[0]));
				}
			}
			IO.pl("\t"+b.getName()+"\t"+numCalls+"\t"+affectedGeneStats.size()+"\t"+affectedGeneStats.keySet());
			ks.getCnvCallSets().add(new CnvCallSet(name, ks.getSampleId(), numCalls, affectedGeneStats));
		}
		
	}
	
	private void loadPatients(File anno) {
		IO.pl("Loading patient and sample meta data...");
		String[] lines = IO.loadFileIntoStringArray(anno);
		if (lines[0].equals("#HCIPersonID\tSampleID\tType\tDateDrawn") == false) Misc.printErrAndExit("First line in the anno file isn't '#HCIPersonID SampleID Type	DateDrawn'? -> "+lines[0]);
		for (int i=1; i< lines.length; i++) {
			String[] tokens = Misc.TAB.split(lines[i]);
			KohliPatient kp = patients.get(tokens[0]);
			if (kp == null) {
				kp = new KohliPatient(tokens);
				patients.put(tokens[0], kp);
			}
			else {
				boolean isGermline = false;
				if (tokens[2].toLowerCase().contains("germline")) isGermline = true;
				kp.getSamples().add(new KohliSample(kp, tokens[1], isGermline, tokens[3]));
			}
		}
		//load the sample treemap
		for (KohliPatient kp: patients.values()) {
			for (KohliSample ks: kp.getSamples()) {
				if (samples.containsKey(ks.getSampleId())) Misc.printErrAndExit("Duplicate sample found! "+ks.getSampleId());
				samples.put(ks.getSampleId(), ks);
			}
		}
	}

}
