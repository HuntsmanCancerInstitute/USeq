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

public class CompareKohliDatasets {
	
	private TreeMap<String, KohliPatient> patients = new TreeMap<String,KohliPatient>();
	private TreeMap<String, KohliSample> samples = new TreeMap<String,KohliSample>();
	private LinkedHashMap<String, KohliGene> genes = null;
	private File resultsDirectory = null;
	//private String[] testedGenes = new String[]{"AR-Enh","ARID1A","PTEN","CCND1","ATM","ZBTB16","CDKN1B","KMT2D","CDK4","MDM2","BRCA2","RB1","FOXA1","AKT1","ZFHX3","TP53","NCOR1","CDK12","BRCA1","SPOP","RNF43","MSH2","ERG","TMPRSS2","CHEK2","CTNNB1","FOXP1","PIK3CB","PIK3CA","PIK3R1","CHD1","APC","CDK6","BRAF","KMT2C","NKX3-1","NCOA2","MYC","COL22A1","CDKN2A","NOTCH1","AR","OPHN1"};
	private String[] testedGenes = new String[]{"AR","BRCA2","CHEK2","MYC","NKX3-1","OPHN1","PIK3CA","PIK3CB","TP53","ZBTB16"};
	private HashMap<String, String> geneNameInfo = null;
	private Bed arEnhancer = new Bed("chrX", 66899317, 66907698);
	private HashMap<String, String> vcfKeyRecord = null;
	private File fullPathToR = new File("/usr/local/bin/R");
	private File saveDir = null;

	public static void main(String[] args) throws IOException {
		File patientSampleInfo = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt");
		File geneInfo = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNADesign/CBioReviewProstateCancer23July2024/targetGeneSummary.txt");
		File gatkPanelResultsDir = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/GatkCopyRatioPassingBeds");
		File gatkWGSResultsDir = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/WgsAnalysis/GATKForWGS/SubSamplingComparision/S0_100_PassBed");
		File smallPanelResults = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/KohliCopyRatio/CopyAlterationCalls/fourNormMethodsFinal.txt");
		File cBioVcf = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNADesign/CBioReviewProstateCancer23July2024/Mutations_Hg19/cBio23July2024MutationsHg38.vcf.gz");
		File somVarCallDir = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/SomVarCalling/Anno");
		File snaqSeqResults = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/AccuGenomicsCopyCalls/snaqSeqStats.txt");
		new CompareKohliDatasets(patientSampleInfo, geneInfo, gatkPanelResultsDir, gatkWGSResultsDir, smallPanelResults, cBioVcf, somVarCallDir, snaqSeqResults);
	}
	
	public CompareKohliDatasets(File patientSampleInfo, File geneInfo, File gatkPanelResultsDir, File gatkWGSResultsDir, File smallPanelResults, File cBioVcf, File somVarCallDir, File snaqSeqResults) throws IOException {

		//load patients and samples
		loadPatients(patientSampleInfo);
		for (KohliPatient kp: patients.values()) IO.pl(kp);
		
		//load gene info
		loadGeneInfo(geneInfo);

		if (true) {
			//load AccuGenomics results
			loadSnaqSeqCopyCalls(snaqSeqResults);
			
			//load dir of cnv results
			loadGatkCopyRatioCalls(gatkPanelResultsDir, "GatkPanel");
			
			//load dir of cnv results
			loadGatkCopyRatioCalls(gatkWGSResultsDir, "GatkWGS");

			//load the small panel method results
			loadSmallPanelCopyRatioCalls(smallPanelResults);
			
			compareResultsWithMergedGatkWgsConfusionMatrix();
		}
		
		if (false) {
			//load cBio vcf
			loadCBioVcf(cBioVcf);
			
			//load som var variant calls
			loadSomVarCalls(somVarCallDir);
			
			//print results
			printVariantResults();
			
		}
		
		IO.pl("\nCOMPLETE!");

	}
	
	
	private void printVariantResults() {
		IO.pl("Printing variant results...");
		//for each patient
		for (KohliPatient kp: patients.values()) {
			IO.pl("\n"+kp.getHciPatientId());
			TreeSet<String> modifiedGenes = new TreeSet<String>();
			
			//for each cfDNA sample
			for (KohliSample ks: kp.getNonGermlineSamples()) {
				if (ks.getIntersectingVariants().size()!=0) IO.pl(ks.toStringVariants());
			}
			/*
			IO.pl("GeneInfo:");
			for (String g: modifiedGenes) IO.pl(g+"\t"+geneNameInfo.get(g));
			IO.pl();
			*/
		}
		
	}


	private void loadSomVarCalls(File somVarCallDir) {
		BufferedReader in = null;
		try {
			//for each 1274089_22712X12_23802X4.merged.anno.vcf.gz
			for (File v: IO.extractFiles(somVarCallDir, ".vcf.gz")) {
				String[] c = Misc.UNDERSCORE.split(v.getName());
				
				//pull sample
				KohliSample ks = samples.get(c[1]);
				if (ks == null) Misc.printErrAndExit("Failed to find a KohliSample for "+ v.getName());
				
				//load the vcfs
				in = IO.fetchBufferedReader(v);
				String[] f = null;
				String line = null;
				boolean foundVars = false;
				while ((line = in.readLine())!=null) {
					line = line.trim();
					if (line.length()==0 || line.startsWith("#")) continue;
					f = Misc.TAB.split(line);
					String key = f[0]+"_"+f[1]+"_"+f[3]+"_"+f[4];
					if (vcfKeyRecord.containsKey(key)){
						ks.getIntersectingVariants().add("\tcBio\t"+vcfKeyRecord.get(key)+"\n\tsVar\t"+line);
						foundVars = true;
					}
				}
				if (foundVars == false) ks.getIntersectingVariants().add("\tNo intersecting variants");
			}

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		} finally {
			IO.closeNoException(in);
		}

	}

	private void loadCBioVcf(File cBioVcf) {
		IO.pl("Loading cBioPortal snvs and indels...");
		//#CHROM	POS	   ID	REF	ALT	QUAL FILTER	INFO	      FORMAT	21310X10_Hg38_final.bam  21310X1_Hg38_final.bam      ...
		//chr1	11199396	.	AG	A	0	  .	    INDEL;IDV=75  PL:AD	       0,255,239:5977,59	       0,255,237:13580,12
		//  0      1        2   3   4   5     6       7             8                9                          10                11...
		vcfKeyRecord = new HashMap<String,String>();
		BufferedReader in = null;
		try {
			in = IO.fetchBufferedReader(cBioVcf);
			String[] f = null;
			String line = null;
			while ((line = in.readLine())!=null) {
				line = line.trim();
				if (line.length()==0 || line.startsWith("#")) continue;
				f = Misc.TAB.split(line);
				String key = f[0]+"_"+f[1]+"_"+f[3]+"_"+f[4];
				vcfKeyRecord.put(key, line);
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		} finally {
			IO.closeNoException(in);
		}

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

	private void compareResultsAll() {
		IO.pl("Comparing copy number results...");
		//for each patient
		for (KohliPatient kp: patients.values()) {
			TreeSet<String> modifiedGenes = new TreeSet<String>();
			
			
			//get samples, for the germline sample, might be null
			KohliSample germline = kp.getGermlineSample();
			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			
			//header
			IO.pl("Patient:\t"+kp.getHciPatientId());
			

			if (germline != null) {
				IO.pl("Germline:\t"+germline.getSampleId());
				for (CnvCallSet cc: germline.getCnvCallSets()) {
					IO.pl("\t"+cc.getMethodName()+":\t"+cc.getCopyAlteredGeneString(testedGenes));
					for (String gd: cc.getCopyAlteredGenes(testedGenes)) modifiedGenes.add(gd.substring(0, gd.length()-1));
				}
			}
			
			
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.pl("cfDNA:\t"+ks.getSampleId());
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					IO.pl("\t"+cc.getMethodName()+":\t"+cc.getCopyAlteredGeneString(testedGenes));
					for (String gd: cc.getCopyAlteredGenes(testedGenes)) modifiedGenes.add(gd.substring(0, gd.length()-1));
				}
			}
			
			//gene info

			if (modifiedGenes.size()!=0) {
				IO.pl("GeneInfo:");
				for (String g: modifiedGenes) IO.pl(g+"\t"+geneNameInfo.get(g));
			}
			
			IO.pl();
		}
		
	}
	
	private void compareResults() {
		IO.pl("Comparing copy number results...");
		//for each patient
		for (KohliPatient kp: patients.values()) {
			TreeSet<String> modifiedGenes = new TreeSet<String>();
			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.p(kp.getHciPatientId()+"\t"+ks.getSampleId());
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					IO.p("\t"+cc.getMethodName()+": "+cc.getCopyAlteredGeneString(testedGenes));
					for (String gd: cc.getCopyAlteredGenes(testedGenes)) modifiedGenes.add(gd.substring(0, gd.length()-1));
				}
				IO.pl();
			}	
		}
	}
	
	private void compareResultsWithIChor() {
		IO.pl("Comparing copy number results to ichor...");
		//for each patient
		for (KohliPatient kp: patients.values()) {
			TreeSet<String> modifiedGenes = new TreeSet<String>();
			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.p(kp.getHciPatientId()+"\t"+ks.getSampleId());
				//find the ichor
				CnvCallSet ichor = null;
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("ichor")) {
						ichor = cc;
						break;
					}
				}
				if (ichor == null) {
					IO.pl("\tichor not found....");
					continue;
				}
				
				
				//add ichor results
				IO.p("\t"+ichor.getMethodName()+": "+ichor.getCopyAlteredGeneString(testedGenes));
				
				//make has of ichor results
				//AR+2, BRCA2+1, CHEK2+2
				ArrayList<String> ichorGenesDirection = ichor.getCopyAlteredGenes(testedGenes);
				HashSet<String> ichorResults = new HashSet<String>();
				for (String gd: ichorGenesDirection) ichorResults.add(gd.substring(0, gd.length()-1));
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("ichor")) continue;
					IO.p("\t"+cc.getMethodName()+": "+cc.getCopyAlteredGeneString(testedGenes));
					ArrayList<String> genesDirection = cc.getCopyAlteredGenes(testedGenes);
					int numMatchIchor = 0;
					int numNoMatchIchor = 0;
					for (String gd: genesDirection) {
						modifiedGenes.add(gd.substring(0, gd.length()-1));
						if (ichorResults.contains(gd)) numMatchIchor++;
						else numNoMatchIchor++;
					}
					IO.p("\t"+numMatchIchor+"/"+ichorResults.size()+"\t"+numNoMatchIchor+"/"+ichorResults.size()+"\t"+numNoMatchIchor+"/"+(numMatchIchor+numNoMatchIchor));
				}
				IO.pl();
			}	
		}
	}
	
	private void compareResultsWithGatkWGS() {
		IO.pl("Comparing copy number results to GatkWGS...");
		//for each patient
		for (KohliPatient kp: patients.values()) {
			TreeSet<String> modifiedGenes = new TreeSet<String>();
			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.p(kp.getHciPatientId()+"\t"+ks.getSampleId());
				//find the GatkWGS
				CnvCallSet gatkWGS = null;
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("GatkWGS")) {
						gatkWGS = cc;
						break;
					}
				}
				if (gatkWGS == null) {
					IO.pl("\tGatkWGS not found....");
					continue;
				}
				
				
				//add gatkWGS results
				IO.p("\t"+gatkWGS.getMethodName()+": "+gatkWGS.getCopyAlteredGeneString(testedGenes));
				
				//make hash of gatkWGS results
				//AR+2, BRCA2+1, CHEK2+2
				ArrayList<String> gatkWGSGenesDirection = gatkWGS.getCopyAlteredGenes(testedGenes);
				HashSet<String> gatkWGSResults = new HashSet<String>();
				for (String gd: gatkWGSGenesDirection) gatkWGSResults.add(gd);
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("GatkWGS")) continue;
					IO.p("\t"+cc.getMethodName()+": "+cc.getCopyAlteredGeneString(testedGenes));
					ArrayList<String> genesDirection = cc.getCopyAlteredGenes(testedGenes);
					int numMatchGatkWGS = 0;
					int numNoMatchGatkWGS = 0;
					for (String gd: genesDirection) {
						modifiedGenes.add(gd.substring(0, gd.length()-1));
						if (gatkWGSResults.contains(gd)) numMatchGatkWGS++;
						else numNoMatchGatkWGS++;
					}
					IO.p("\t"+numMatchGatkWGS+"/"+gatkWGSResults.size()+"\t"+numNoMatchGatkWGS+"/"+gatkWGSResults.size()+"\t"+numNoMatchGatkWGS+"/"+(numMatchGatkWGS+numNoMatchGatkWGS));
				}
				IO.pl();
			}	
		}
	}
	
	private void compareResultsWithGatkWGSMerged() {
		IO.pl("Comparing copy number results to GatkWGS 2.5KB and 10KB...");
		//for each patient
		for (KohliPatient kp: patients.values()) {
			TreeSet<String> modifiedGenes = new TreeSet<String>();
			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.p(kp.getHciPatientId()+"\t"+ks.getSampleId());
				//find the GatkWGS2.5 and GatkWGS10
				CnvCallSet gatkWGS2 = null;
				CnvCallSet gatkWGS10 = null;
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("GatkWGS2.5")) gatkWGS2 = cc;
					else if (cc.getMethodName().equals("GatkWGS10")) gatkWGS10 = cc;
				}
				if (gatkWGS2 == null || gatkWGS10 == null) {
					IO.pl("\tGatkWGSs not found....");
					continue;
				}
				

				//add gatkWGS results
				String wgs2Genes = gatkWGS2.getCopyAlteredGeneString(testedGenes);
				String wgs10Genes = gatkWGS10.getCopyAlteredGeneString(testedGenes);
				TreeSet<String> mergedGenes = mergeGeneStrings(wgs2Genes, wgs10Genes);
				String merged = Misc.treeSetToString(mergedGenes, ",");
				IO.p("\tGATK WGS: "+merged);
				
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("GatkWGS2.5") || cc.getMethodName().equals("GatkWGS10")) continue;
					IO.p("\t"+cc.getMethodName()+": "+cc.getCopyAlteredGeneString(testedGenes));
					ArrayList<String> genesDirection = cc.getCopyAlteredGenes(testedGenes);
					int numMatchGatkWGS = 0;
					int numNoMatchGatkWGS = 0;
					for (String gd: genesDirection) {
						modifiedGenes.add(gd.substring(0, gd.length()-1));
						if (mergedGenes.contains(gd)) numMatchGatkWGS++;
						else numNoMatchGatkWGS++;
					}
					int numMergedGenes = 0;
					if (mergedGenes.contains("None") == false) numMergedGenes = mergedGenes.size();
					//IO.p("\t"+numMatchGatkWGS+"/"+numMergedGenes+"\t"+numNoMatchGatkWGS+"/"+numMergedGenes+"\t"+numNoMatchGatkWGS+"/"+(numMatchGatkWGS+numNoMatchGatkWGS));

					
					//FDR=TestNonMatches/TotalTest
					String fdr = formatFraction(numNoMatchGatkWGS, genesDirection.size());
					
					//Precision PPV=TestMatches/TotalTest
					String ppv = formatFraction(numMatchGatkWGS, genesDirection.size());
					double precision = (double)numMatchGatkWGS/ (double)genesDirection.size();
					
					//Recall TPR=KeyMatches/TotalKey
					String tpr = formatFraction(numMatchGatkWGS, numMergedGenes);
					double recall = (double)numMatchGatkWGS / (double)numMergedGenes;
					
					//FPR=KeyNonMatchs/TotalKey, wrong!!!!
					String fpr = formatFraction(numMergedGenes-numMatchGatkWGS, numMergedGenes);

					//F-score=HarmonicMean(Precision,Recall)
					String fscore = Num.formatNumber(Num.harmonicMean(new double[] {precision, recall}), 4);
					
					IO.p("\t"+numMatchGatkWGS+"\t"+numNoMatchGatkWGS+ "\t"+fdr+"\t"+ppv+"\t"+tpr+"\t"+fpr+"\t"+fscore);
					

				
				}
				IO.pl();
			}	
		}
	}
	
	private void compareResultsWithMergedGatkWgsConfusionMatrix() {
		IO.pl("\nComparing copy number results to GatkWGS 2.5KB and 10KB confusion matrix...");
		IO.pl("Header:\t"+ConfusionMatrix.toStringHeader());
		//for each patient
		for (KohliPatient kp: patients.values()) {

			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.p(kp.getHciPatientId()+"\t"+ks.getSampleId());
				
				//find the GatkWGS2.5 and GatkWGS10
				CnvCallSet gatkWGS2 = null;
				CnvCallSet gatkWGS10 = null;
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("GatkWGS2.5")) gatkWGS2 = cc;
					else if (cc.getMethodName().equals("GatkWGS10")) gatkWGS10 = cc;
				}
				if (gatkWGS2 == null || gatkWGS10 == null) {
					IO.pl("\tGatkWGSs not found....");
					continue;
				}

				//add gatkWGS results
				String wgs2Genes = gatkWGS2.getCopyAlteredGeneString(testedGenes);
				String wgs10Genes = gatkWGS10.getCopyAlteredGeneString(testedGenes);
				TreeSet<String> keyGeneCalls = mergeGeneStrings(wgs2Genes, wgs10Genes);
				String merged = Misc.treeSetToString(keyGeneCalls, ",");
				IO.p("\tGATK WGS: "+merged);
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("GatkWGS2.5") || cc.getMethodName().equals("GatkWGS10")) continue;
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



	public TreeSet<String> mergeGeneStrings(String one, String two) {
		TreeSet<String> merge = new TreeSet<String>();
		String[] oneSplit = Misc.COMMA.split(one);
		String[] twoSplit = Misc.COMMA.split(two);
		for (String x: oneSplit) merge.add(x);
		for (String x: twoSplit) merge.add(x);
		if (merge.size()>1 && merge.contains("None")) merge.remove("None");
		return merge;
	}

	private void loadSmallPanelCopyRatioCalls(File smallPanelResults) {
		IO.pl("Loading Small Panel CR results...");
		IO.pl("\tDataset\tNumCalls\tNumWithCopyAlteration\t");

			
			//load all of the results
			String[] lines = IO.loadFileIntoStringArray(smallPanelResults);
			//parse the header
			String[] header = Misc.TAB.split(lines[0]);
			HashMap<String, Integer> nameIndex = new HashMap<String,Integer>();
			for (int i=0; i< header.length; i++) nameIndex.put(header[i], i);
			if (nameIndex.containsKey("Patient") == false) Misc.printErrAndExit("Failed to find the header 'Patient' in "+lines[0]+" from "+smallPanelResults);
			int geneIndex = nameIndex.get("Gene");
			int sampleIndex = nameIndex.get("Test Sample");
			int callIndex = nameIndex.get("CopyAltered");
			int log2Ratio = nameIndex.get("Log2 (Test/ PoN)");
			int normMethIndex = nameIndex.get("Norm Method");
			
			HashMap<String, CnvCallSet> appPckBstScCnvCallSets = new HashMap<String, CnvCallSet>();
			HashMap<String, CnvCallSet> appPckGnNmCnvCallSets = new HashMap<String, CnvCallSet>();
			HashMap<String, CnvCallSet> manCnvCallSets = new HashMap<String, CnvCallSet>();
			HashMap<String, CnvCallSet> accuPckCnvCallSets = new HashMap<String, CnvCallSet>();
			
			//for each line
			for (String line: lines) {
				//Patient, Test Sample, Germline, Used To Norm, Gene, # Capture Regions, # CRs Test Outside PoN, Mean CR Z-Score, 
				//    0         1           2           3         4          5                     6                   7
				//PseudoMedian CR Z-Score, Mean CR Scaled Test, Mean CR Scaled PoN, Log2(Test/PoN), TTestPVals, 
				//             8                  9                     10                11             12     
				//TTestCombineAdjPVal, TTestFracPassAdjPVals, CopyAltered, NormMethod		
				//          13                  14                 15          16
				if (line.startsWith("Patient")) continue;

				String[] f = Misc.TAB.split(line);
				
				//pull KohliSample
				KohliSample ks = samples.get(f[sampleIndex]);
				if (ks == null) Misc.printErrAndExit("\nFailed to find the sample associated with "+line);
				
				//pull the CnvCallSet
				HashMap<String, CnvCallSet> hm = null;
				String normMeth = null;
				
				if (f[normMethIndex].equals("AppPckBstSc")) {
					hm = appPckBstScCnvCallSets;
					normMeth="AppPckBstSc";
				}
				else if (f[normMethIndex].equals("AppPckGnNm")) {
					hm = appPckGnNmCnvCallSets;
					normMeth="AppPckGnNm";
				}
				else if (f[normMethIndex].equals("ManPck")) {
					hm = manCnvCallSets;
					normMeth="ManPck";
				}
				else if (f[normMethIndex].equals("AccuPck")) {
					hm = accuPckCnvCallSets;
					normMeth="AccuPck";
				}
				else Misc.printErrAndExit("\nFailed to find the norm method with line "+line);
				
				CnvCallSet ccs = hm.get(ks.getSampleId());
				if (ccs == null) {
					ccs = new CnvCallSet(normMeth, ks.getSampleId(), 0, new TreeMap<String, GeneCallResult>());
					hm.put(ks.getSampleId(), ccs);
				}
				
				//pull the gene call results
				TreeMap<String, GeneCallResult> gcr = ccs.getGeneCallResults();
				
				//add a new results, boolean isCopyAltered, boolean isAmplified, String statistics
				String call = f[callIndex]; //Pass or Fail
				if (call.equals("Pass") == false && call.equals("Fail") == false) Misc.printErrAndExit("\nFailed to parse 'Pass' or 'Fail' copy call from "+line);
				double ratio = Double.parseDouble(f[log2Ratio]);
				GeneCallResult res = new GeneCallResult(call.equals("Pass"), ratio>0, line);
				gcr.put(f[geneIndex], res);
				
			}

			//fix the count, set it in the sample
			HashMap<String, CnvCallSet>[] hms = new HashMap[4];
			hms[0]=appPckBstScCnvCallSets;
			hms[1]=appPckGnNmCnvCallSets;
			hms[2]=manCnvCallSets;
			hms[3]=accuPckCnvCallSets;
			for (HashMap<String, CnvCallSet> hm: hms) {
				for (CnvCallSet ccs: hm.values()) {
					int numPass = 0;
					for (GeneCallResult res: ccs.getGeneCallResults().values()) if (res.isCopyAltered()) numPass++;
					ccs.setNumberCalls(numPass);
					KohliSample ks = samples.get(ccs.getSampleId());
					ks.getCnvCallSets().add(ccs);
					IO.pl("\t"+ks.getSampleId()+"\t"+ccs.getGeneCallResults().size()+"\t"+numPass);
				}
			}
	}
	
	private void loadSnaqSeqCopyCalls(File snaqSeqStatsFile) {
		IO.pl("\nLoading SnaqSeq CR results...");
		IO.pl("\tDataset\tNumCalls\tNumWithCopyAlteration\t");

		//load file
		String[] lines = IO.loadFileIntoStringArray(snaqSeqStatsFile);
		//parse the header
		String[] header = Misc.TAB.split(lines[0]);
		HashMap<String, Integer> nameIndex = new HashMap<String,Integer>();
		for (int i=0; i< header.length; i++) nameIndex.put(header[i], i);
		if (nameIndex.containsKey("gene") == false) Misc.printErrAndExit("Failed to find the header 'gene' in in "+lines[0]+" from "+snaqSeqStatsFile);
		int geneIndex = nameIndex.get("gene");
		int sampleIndex = nameIndex.get("short_name");
		int callIndex = nameIndex.get("comparitor-snaq");
		
		HashMap<String, CnvCallSet> sampleCnvCallSets = new HashMap<String, CnvCallSet>();
		
		for (int i=1; i< lines.length; i++) {
			String[] f = Misc.TAB.split(lines[i]);

			//pull KohliSample
			KohliSample ks = samples.get(f[sampleIndex]);
			if (ks == null) Misc.printErrAndExit("\nFailed to find the sample associated with "+lines[i]);
			
			//pull the CnvCallSet
			CnvCallSet ccs = sampleCnvCallSets.get(ks.getSampleId());
			if (ccs == null) {
				ccs = new CnvCallSet("SnaqSeq", ks.getSampleId(), 0, new TreeMap<String, GeneCallResult>());
				sampleCnvCallSets.put(ks.getSampleId(), ccs);
			}
			
			//pull the gene call results
			TreeMap<String, GeneCallResult> gcr = ccs.getGeneCallResults();
			
			//add a new results, boolean isCopyAltered, boolean isAmplified, String statistics
			int call = Integer.parseInt(f[callIndex]);
			GeneCallResult res = new GeneCallResult(call!=0, call>0, lines[i]);
			gcr.put(f[geneIndex], res);
		}
		
		//fix the count, set it in the sample
		for (CnvCallSet ccs: sampleCnvCallSets.values()) {
			int numPass = 0;
			for (GeneCallResult res: ccs.getGeneCallResults().values()) if (res.isCopyAltered()) numPass++;
			ccs.setNumberCalls(numPass);
			KohliSample ks = samples.get(ccs.getSampleId());
			ks.getCnvCallSets().add(ccs);
			IO.pl("\t"+ks.getSampleId()+"\t"+ccs.getGeneCallResults().size()+"\t"+numPass);
		}
	}


	private void loadGatkCopyRatioCalls(File gatkResultsDir, String name) {
		IO.pl("Loading GATK "+name+" results...");
		IO.pl("\tDataset\tNumBedLines\tNumInterrogatedGenesCalled");
		File[] beds = IO.extractFiles(gatkResultsDir, "called.seg.pass.bed");
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
	
	private void loadIChorCopyCalls(File ichorResultsDir) {
		IO.pl("Loading ichorCna results...");
		IO.pl("\tDataset\tNumDataLines\tNumInterrogatedGenesCalled");
		File[] annoSegs = IO.extractFiles(ichorResultsDir, ".anno.seg.txt.gz");
		HashMap<String, Integer> callDirection = fetchIChorCalls();
		HashMap<String, String> geneNameDirection = fetchPrioritizedCnvDirection();
				
		//for each results file
		for (File b: annoSegs) {
			//1110039_22712G28_23802G10.anno.seg.txt.gz and 1110039_22712G28_NA.anno.seg.txt.gz
			String[] patient_cfDNA_germline = Misc.UNDERSCORE.split(b.getName());
			
			//pull KohliSample
			KohliSample ks = samples.get(patient_cfDNA_germline[1].replaceFirst("G", "X"));
			if (ks == null) Misc.printErrAndExit("\nFailed to find the sample associated with "+b);
			//KohliSample ks = new KohliSample(null, null, false, null);
			
			//load all of the lines
			String[] lines = IO.loadFile(b);
			int numCalls = lines.length;
			TreeMap<String, GeneCallResult> affectedGeneStats = new TreeMap<String, GeneCallResult>();
			
			//for each bed line
			for (String line: lines) {
				//ID chrom start end num.mark seg.median.logR copy.number call subclone.status logR_Copy_Number Corrected_Copy_Number Corrected_Call Gene_Names 
				//0     1     2   3      4          5              6        7          8              9                    10                 11         12
				String[] f = Misc.TAB.split(line);
				String[] genesAffected = Misc.COMMA.split(f[12]);
				Integer call = callDirection.get(f[7]);
				boolean cnvCalled = call !=0;
				boolean amplified = call > 0;
				HashSet<String> gas = Misc.loadHashSet(genesAffected);
				String stats =  f[4]+","+ f[5]+","+ f[7]+","+ f[9]; 
				//for each interrogated gene
				for (String intGene: testedGenes) {
					if (gas.contains(intGene)) {
						GeneCallResult gcr = new GeneCallResult(cnvCalled, amplified, stats);
						//is an existing call already present?
						if (affectedGeneStats.containsKey(intGene)) {
							GeneCallResult prior = affectedGeneStats.get(intGene);
							String plusOrMinus = geneNameDirection.get(intGene);
							int priorCallNumber = Integer.parseInt(prior.getIchorCode()); 
							if (priorCallNumber < 0 && plusOrMinus.equals("+")) affectedGeneStats.put(intGene, gcr);
							else if (priorCallNumber > 0 && plusOrMinus.equals("-")) affectedGeneStats.put(intGene, gcr);
						}
						else affectedGeneStats.put(intGene, gcr);
					}
				}
				//does it intersect the AR Enhancer?
				Bed bed = new Bed (f[1], Integer.parseInt(f[2]), Integer.parseInt(f[3]));
				if (bed.intersects(arEnhancer)) {
					String intGene = "AR-Enh";
					GeneCallResult gcr = new GeneCallResult(cnvCalled, amplified, stats);
					if (affectedGeneStats.containsKey(intGene)) {
						GeneCallResult prior = affectedGeneStats.get(intGene);
						String plusOrMinus = geneNameDirection.get(intGene);
						int priorCallNumber = Integer.parseInt(prior.getIchorCode()); 
						if (priorCallNumber < 0 && plusOrMinus.equals("+")) affectedGeneStats.put(intGene, gcr);
						else if (priorCallNumber > 0 && plusOrMinus.equals("-")) affectedGeneStats.put(intGene, gcr);
					}
					else affectedGeneStats.put(intGene, gcr);
				}
			}
			IO.pl("\t"+b.getName()+"\t"+numCalls+"\t"+affectedGeneStats.size()+"\t"+affectedGeneStats.keySet());
			ks.getCnvCallSets().add(new CnvCallSet("ichor", ks.getSampleId(), numCalls, affectedGeneStats));
			
		}	
	}	
	
	private static HashMap<String, Integer> fetchIChorCalls(){
		HashMap<String, Integer> callDir = new HashMap<String, Integer>();
		callDir.put("HETD", -1);
		callDir.put("HOMD", -2);
		callDir.put("NEUT", 0);
		callDir.put("GAIN", 1);
		callDir.put("AMP", 2);
		callDir.put("HLAMP", 3);
		return callDir;
	}

	public HashMap<String, String> fetchPrioritizedCnvDirection() {
		String prioritizedCnvDirection= "CDKN1B -,CDKN2A -,CHD1 -,RNF43 -,PIK3R1 -,ARID1A -,BRCA1 -,BRCA2 -,CDK12 -,CHEK2 -,FOXP1 -,NCOR1 -,PTEN -,TP53 -,APC -,ATM -,KMT2C -,KMT2D -,MSH2 -,NKX3-1 -,RB1 -,ZFHX3 -,ZBTB16 -,AKT1 +,BRAF +,CCND1 +,CDK4 +,COL22A1 +,CTNNB1 +,FOXA1 +,MDM2 +,NOTCH1 +,PIK3CA +,PIK3CB +,SPOP +,ERG +,TMPRSS2 +,MYC +,AR +,AR-Enh +,OPHN1 +,CDK6 +,NCOA2 +";
		HashMap<String, String> geneNameDirection = new HashMap<String, String>();
		String[] genes = Misc.COMMA.split(prioritizedCnvDirection);
		for (String s: genes) {
			String[] geneDir = Misc.WHITESPACE.split(s);
			geneNameDirection.put(geneDir[0], geneDir[1]);
		}
		return geneNameDirection;
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
