package edu.utah.seq.vcf.anno; 

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Parser of vcf files annotated with SnpEff (dbNSFP,ClinVar), VCFBackgroundChecker, and VCFSpliceScanner.  */
public class AnnotatedVcfParser {

	//user input and default settings, don't modify defaults here
	private File[] vcfFiles;
	private File saveDirectory;
	private File transcriptList = null;
	private String clinvarDate = null;
	private int minimumDP = 0;
	private double minimumAF = 0;
	private double minimumAlt = 0;
	private double maximumAF = 1;
	private double maximumCF = 1;
	private double maximumPopAF = 0;
	private double maxFracBKAFs = 0; //maximum number of background samples with >= AF, measure as fraction of total samples, 50 for foundation so 5 max
	private double minimumBKZ = 0;
	private double minFractionPathogenic = 0.51;
	private String[] passingIDKeys = null;
	private TreeSet<String> passingAnnImpact = null;
	private TreeSet<String> passingClinSig = null;
	private TreeSet<String> excludeClinSig = null;
	private TreeSet<String> passingGermlineGenes = null;
	private TreeSet<String> drugResClinSigGenes = null;
	private String[] passingVCFSS = null;
	private double minimumVCFSSDiff = 0; //difference between ref and alt for splice junction score
	private boolean verbose = false;
	private boolean somaticProcessing = false;
	private boolean orAnnos = false;
	private boolean skipWarningTrans = true;
	private boolean onlyProteinCoding = true;
	private boolean justFrameShiftStartStop = false;
	private String annInfo = "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' \">";
	private Gzipper sumarySpreadSheet = null;
	private static final Pattern COMMA_SLASH = Pattern.compile(",|/");
	private String appSettings = null;
	private HashSet<String> transcriptFilter = null;
	private Gzipper impactedGenes = null;
	
	//trackers
	private Histogram afs = new Histogram(0, 1.01, 101);
	private TreeMap<String, Integer> clinsig = new TreeMap<String, Integer>();
	private int[] dps = new int[10000];
	private TreeMap<String, Integer> impacts = new TreeMap<String, Integer>();
	private TreeMap<String, Integer> effects = new TreeMap<String, Integer>();
	private String acmgGenes = "ACTA2,ACTC1,APC,APOB,ATM,ATP7B,BARD1,BMPR1A,BRCA1,BRCA2,BRIP1,CACNA1S,CDH1,CDK4,CDKN2A,CHEK2,COL3A1,DSC2,DSG2,DSP,"
			+ "EPCAM,FBN1,GLA,GREM1,KCNH2,KCNQ1,LDLR,LMNA,MEN1,MLH1,MSH2,MSH6,MUTYH,MYBPC3,MYH7,MYH11,MYL2,MYL3,NBN,NF2,OTC,PALB2,PCSK9,"
			+ "PKP2,PMS2,POLD1,POLE,PRKAG2,PTEN,RAD51C,RAD51D,RB1,RET,RYR1,RYR2,SCN5A,SDHAF2,SDHB,SDHC,SDHD,SMAD3,SMAD4,STK11,TGFBR1,"
			+ "TGFBR2,TMEM43,TNNI3,TNNT2,TP53,TPM1,TSC1,TSC2,VHL,WT1";
	private TreeMap<String, Integer> observedGermlineGenes = new TreeMap<String, Integer>();
	
	//counters
	private int numRecords = 0;
	private int numPassingDP = 0;
	private int numPassingAlt = 0;
	private int numPassingMaxAF = 0;
	private int numPassingMinAF = 0;
	private int numWithPopAF = 0;
	private int numPassingPopAF = 0;
	private int numWithAnn = 0;
	private int numPassingAnnImpact = 0;
	private int numWithBKAFs = 0;
	private int numPassingBKAFs = 0;
	private int numWithClin = 0;
	private int numPassingClinSing = 0;
	private int numWithCF = 0;
	private int numPassingCF = 0;
	private int numPassingExcludeClinSing = 0;
	private int numWithSplice = 0;
	private int numPassingSplice = 0;
	private int numPassIDs = 0;
	private int numPass = 0;
	
	//working fields
	private File workingVcf = null;
	private String vcfLine = null;
	private double obsAF;
	private AnnotatedVcfParserDataLine dataLine = null;
	private String trimmedFileName = null;
	private HashMap<String,String> infoKeyValue = new HashMap<String,String>();
		
	public AnnotatedVcfParser (String[] args) {
		long startTime = System.currentTimeMillis();
		try {
			
		processArgs(args);
		
		if (somaticProcessing) modifySettingsForFoundation();
		
		openSpreadSheet();
		impactedGenes = new Gzipper( new File (saveDirectory, "impactedGenes.txt.gz"));
		
		IO.pl("Parsing (#fail #pass):");
		for (File vcf: vcfFiles){
			workingVcf = vcf;
			System.out.print("\t"+workingVcf.getName());
			if (verbose) IO.pl();
			parse();
		}
		
		printSettings();
		printStats();
		printDistributions();
		closeSpreadsheet();
		impactedGenes.close();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" sec\n");
		
		} catch (Exception e) {
			System.err.println("\nERROR: parsing vcf file: "+workingVcf.getName()+"\nLine: "+vcfLine+"\n");
			e.printStackTrace();
			System.exit(1);
		} 
	}
	
	
	private void closeSpreadsheet() throws IOException {
		sumarySpreadSheet.println(appSettings);
		sumarySpreadSheet.println(AnnotatedVcfParserDataLine.legend);
		
		sumarySpreadSheet.close();
		
	}


	private void openSpreadSheet() throws FileNotFoundException, IOException {
		String date = Misc.getDateNoSpaces();
		File sss = new File(saveDirectory, "annotatedVcfParser."+date+".xls.gz");
		sumarySpreadSheet = new Gzipper(sss);
		sumarySpreadSheet.print(AnnotatedVcfParserDataLine.headerSpreadSheet);
		if (transcriptFilter != null) sumarySpreadSheet.println(AnnotatedGene.headerSpreadSheetMatch);
		else sumarySpreadSheet.println(AnnotatedGene.headerSpreadSheet);
	}


	public void parse() throws Exception{
		Gzipper passVcf = null;
		Gzipper failVcf = null;

			//IO
			trimmedFileName = Misc.removeExtension(workingVcf.getName());			
			passVcf = new Gzipper( new File (saveDirectory, trimmedFileName + "_Pass.vcf.gz"));
			failVcf = new Gzipper( new File (saveDirectory, trimmedFileName + "_Fail.vcf.gz"));
			
			
			BufferedReader in = IO.fetchBufferedReader(workingVcf);
			int numPassingVcf = 0;
			int numFailingVcf = 0;
			TreeSet<String> passingGeneNames = new TreeSet<String>();
			//for each line in the file
			while ((vcfLine = in.readLine()) != null){
				vcfLine = vcfLine.trim();
				if(vcfLine.length() == 0) continue;
				//header? just print out
				if (vcfLine.startsWith("#")) {
					passVcf.println(vcfLine);
					failVcf.println(vcfLine);
					//ANN?
					if (vcfLine.startsWith("##INFO=<ID=ANN,") && vcfLine.equals(annInfo) == false) throw new Exception("Your ##INFO=<ID=ANN line  doesn't match\n"+vcfLine+"\n"+annInfo);
					continue;
				}
				numRecords++;
				if (verbose) IO.pl("\n"+vcfLine);

				//#CHROM POS ID REF ALT QUAL FILTER INFO ......
				//   0    1   2  3   4   5     6      7
				String[] cells = Misc.TAB.split(vcfLine);
				loadInfoHash(cells[7]);
				dataLine = new AnnotatedVcfParserDataLine(trimmedFileName, cells, clinvarDate);
				
				boolean passDP=true, passAF=true, passImpact=true, passSplice=true, passID = true, passGermlineGenes = false, passAlt=true;
				boolean[] passCF= null;
				boolean[] passBKAF= null; 
				boolean[] passPop= null;
				boolean[] passClinSig= null;
				
				
				if (passingGermlineGenes != null) passGermlineGenes = checkGermlineGenes();
				if (minimumDP != 0) 			passDP = checkDP();
				if (maximumCF != 1) 			passCF = checkCF();
				if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) passAF = checkAF(passGermlineGenes);
				if (minimumAlt !=0)             passAlt = checkAlt(passGermlineGenes);
				if (maximumPopAF != 0) 			passPop = checkPopFreq();
				if (passingAnnImpact != null) 	passImpact = checkImpact();
				if (passingClinSig != null || excludeClinSig != null) 	passClinSig = checkClinSig();
				if (maxFracBKAFs != 0 || minimumBKZ !=0) passBKAF = checkBKAFs();
				if (passingVCFSS != null) 		passSplice = checkSplice();
				if (passingIDKeys != null) 		passID = checkIDs(cells[2]);
				
				boolean pass;
				if (somaticProcessing) pass =  passFoundationFiltering(passID, passDP, passAF, passAlt, passImpact,passSplice, passPop, passBKAF,passClinSig, passCF); 
				else if (orAnnos) pass = passWithOrAnnos(passID, passDP, passAF, passAlt, passImpact, passSplice, passPop, passBKAF,passClinSig, passCF);
				else pass = allPass(passID, passDP, passAF, passAlt, passImpact, passSplice, passPop, passBKAF, passClinSig, passCF);
				
				if (pass){
					numPassingVcf++;
					passingGeneNames.addAll(dataLine.fetchPassingGeneNames());
					passVcf.println(vcfLine);
					dataLine.println(sumarySpreadSheet, transcriptFilter);
				}
				else {
					failVcf.println(vcfLine);
					numFailingVcf++;
				}
				if (verbose) IO.pl("\tPass all filters\t"+pass);	
			}
			numPass+=numPassingVcf;
			
			if (verbose == false) IO.pl("\t"+numFailingVcf+"\t"+numPassingVcf);
			impactedGenes.println(trimmedFileName+"\t"+Misc.treeSetToString(passingGeneNames, "\t")); 
			
			//close io
			in.close();
			passVcf.close();
			failVcf.close();	
	}
	
	private void loadInfoHash(String info) {
		infoKeyValue.clear();
		String[] kvs = Misc.SEMI_COLON.split(info);
		for (String kv: kvs) {
			int first = kv.indexOf('=');
			if (first !=-1) infoKeyValue.put(kv.substring(0, first), kv.substring(first+1));
		}
	}


	/**Requires all are true. Not if these aren't set off defaults then they default to true. Note if passID then returns true regardless of others.*/
	private boolean allPass(boolean passID, boolean passDP, boolean passAF, boolean passAlt, boolean passImpact, boolean passSplice, boolean[] passPop, boolean[] passBKAF, boolean[] passClinSig, boolean[] passCF) {
		if (passingIDKeys != null && passID) return true;
		if (passDP == false) return false;
		if (passAF == false) return false;
		if (passAlt == false) return false;
		if (passCF!= null && passCF[0] == true && passCF[1] == false) return false;
		boolean passBKAF2 = true;
		if (passBKAF != null && passBKAF[0] == true) passBKAF2 = passBKAF[1];
		if (passBKAF2 == false) return false;
		boolean passPop2 = true;
		if (passPop != null && passPop[0] == true) passPop2 = passPop[1];
		if (passPop2 == false) return false;
		//check splice?
		if (passingVCFSS != null && passSplice == false) return false;
		//check impact?
		if (passingAnnImpact != null && passImpact == false) return false;
		//check clinsig?
		if (passingClinSig != null && passClinSig == null ) return false;
		if (passClinSig != null){
			if (excludeClinSig !=null && passClinSig[1] == true) return false;
			if (passingClinSig != null && passClinSig[0] == false) return false;
		}
		return true;
	}
	
	/**Requires DP, AF, Pop, CF, BKAF are true if present. Only one of Clin, Splice, or Impact need be true, Note if passID then returns true regardless of others.*/
	private boolean passWithOrAnnos(boolean passID, boolean passDP, boolean passAF, boolean passAlt, boolean passImpact, boolean passSplice, boolean[] passPop, boolean[] passBKAF, boolean[] passClinSig, boolean[] passCF) {
		if (passingIDKeys != null && passID) return true;
		if (passDP == false) return false;
		if (passAF == false) return false;
		if (passAlt == false) return false;
		if (passCF!= null && passCF[0] == true && passCF[1] == false) return false;
		boolean passBKAF2 = true;
		if (passBKAF != null && passBKAF[0] == true) passBKAF2 = passBKAF[1];
		if (passBKAF2 == false) return false;
		boolean passPop2 = true;
		if (passPop != null && passPop[0] == true) passPop2 = passPop[1];
		if (passPop2 == false) return false;
		
		if (passClinSig != null){
			//are they excluding vars based on clinSig terms? if one was found fail the record
			if (excludeClinSig !=null && passClinSig[1] == true) return false;
			//are they looking for particular clinsig terms
			if (passingClinSig != null && passClinSig[0] == true) return true;
		}
		
		if (passingVCFSS != null && passSplice == true) return true;
		if (passingAnnImpact != null && passImpact == true) return true;
		return false;
	}
	
	private void createGermlineGenes(){
		passingGermlineGenes = new TreeSet<String>();
		for (String s: Misc.COMMA.split(acmgGenes)) passingGermlineGenes.add(s.trim());
	}

	private boolean passFoundationFiltering(boolean passID, boolean passDP, boolean passAF,  boolean passAlt, boolean passImpact, boolean passSplice, boolean[] passPop, boolean[] passBKAF, boolean[] passClinSig, boolean[] passCF){
		//pass it if ID correct
		if (passID) return true;
		//fail it if DP, AF, alt fails
		if (passAF == false || passDP == false || passAlt == false) return false;
		//if present and fails then fail
		if (passPop[0] == true && passPop[1] == false) return false;
		if (passBKAF[0] == true && passBKAF[1] == false) return false;
		if (passCF!= null && passCF[0] == true && passCF[1] == false) return false;
		//check that one is good, these take priority
		if (passClinSig != null){
			//are they looking for particular clinsig terms
			if (passingClinSig != null && passClinSig[0] == true) return true;
			//are they excluding vars based on clinSig terms? if one was found fail the record
			if (excludeClinSig !=null && passClinSig[1] == true) return false;
		}
		if (passingVCFSS != null && passSplice == true) return true;
		if (passingAnnImpact != null && passImpact == true) return true;
		return false;
	}
	
	public double[] calcFracPathoFromConfInterp() {
		double numPathogenic = 0;
		double numBenign = 0;
		String confString = infoKeyValue.get("CLNSIGCONF");
		if (confString != null) {
			String[] type = Misc.COMMA.split(confString);
			dataLine.clinSigConf = confString;
			for (String t: type) {
				String trimmed = t.substring(0, t.length()-1);
				String[] nameNumber = Misc.FORWARD_PARENTHESIS.split(trimmed);
				//IO.pl(nameNumber[0]);
				String nameLC = nameNumber[0].toLowerCase();
				if (nameLC.contains("pathogenic")) numPathogenic+= Double.parseDouble(nameNumber[1]);
				else if (nameLC.contains("benign")) numBenign+= Double.parseDouble(nameNumber[1]);
			}
		}

		//IO.pl("NumPath "+numPathogenic+"\nNumOther "+numBenign);
		return new double[] {numPathogenic, numBenign};
	}
	
	private boolean checkIDs(String id) {
		
		//for each user provided key (all lower cased), check if in id and pass if so
		String lcId = id.toLowerCase();
		for (String key: passingIDKeys){
			if (lcId.contains(key)){
				numPassIDs++;
				if (verbose) IO.pl("\tID Check\ttrue"+key+"\t"+lcId);
				return true;
			}
		}
		if (verbose) IO.pl("\tID Check\tNA\t"+lcId);
		return false;
	}

	private boolean checkSplice() throws Exception {
		
		String vcfSSLine = infoKeyValue.get("VCFSS");
		//any splice info
		if (vcfSSLine == null){
			if (verbose) IO.pl("\tSplice Check\tfalse\tNo VCFSS");
			return false;
		}
		numWithSplice++;
		
		//split by splice effects, may be more than one
		//VCFSS=ENSG00000163554:G3S,158584103,AAAATTCATGTTTTTTTTTTTTCTTTCAGGGACATCAAAGGTGTG,AAAATTCATGTTTTTTTTTTTTTCTTTCAGGGACATCAAAGGTGTG,-7.52,3.04,3.04
		//&
		//ENSG00000163554:G3E,158584103,AAAATTCATGTTTTTTTTTTTTCTTTCAGGGACATCAAAGGTGTG,AAAATTCATGTTTTTTTTTTTTTCTTTCAGGGACATCAAAGGTGTG,-7.52,3.04,3.04
		String[] se = Misc.AMP.split(vcfSSLine);
		
		//for each one look for the G3E, D3S, and parse the delta
		double maxDelta = -1;
		for (String effect: se){
			String[] tokens = Misc.COLON.split(effect);
			//for each of the user defined effects
			for (String code: passingVCFSS){
				if (tokens[1].startsWith(code)){
					String[] items = Misc.COMMA.split(tokens[1]);
					double delta = Double.parseDouble(items[6]);
					if (delta > maxDelta) {
						maxDelta = delta;
						dataLine.spliceGeneName = se[0];
						dataLine.spliceScoreDiff = maxDelta;
					}
				}
			}
		}
		
		boolean passSplice = maxDelta >= minimumVCFSSDiff;
		if (verbose) IO.pl("\tSplice Check\t"+passSplice+"\t"+maxDelta+"\t"+vcfSSLine);	
		if (passSplice) {
			numPassingSplice++;
			dataLine.passesSplice = true;
			return true;
		}
		else return false;
	}

	
	private boolean[] checkBKAFs() throws Exception {
		boolean foundBK = false;
		boolean passBKs = false;
		
		String bkafs = infoKeyValue.get("BKAF");
		String bkzString = infoKeyValue.get("BKZ");
		
		//watch out for ? in the strings
		if (bkafs!= null && bkafs.contains("?")) bkafs = null;
		if (bkzString!= null && bkzString.contains("?")) bkzString = null;
		
		dataLine.bkz = bkzString;
		dataLine.bkzAF = bkafs;
		
		//both of the following must pass if both are being scored
				
		//check max fraction
		if (maxFracBKAFs !=0){
			if (bkafs != null){
				foundBK = true;
				double[] afs = Num.stringArrayToDouble(bkafs, ",");
				//score
				double numFailing = 0;
				for (double af: afs) if (af>= obsAF) numFailing++;
				double fractFailing = numFailing/(double)afs.length;
				boolean passingBkaf = fractFailing <= maxFracBKAFs; 
				if (verbose) IO.pl("\tBKAF Check\t"+passingBkaf+"\t"+fractFailing+"\t"+bkafs);
				if (passingBkaf) passBKs = true;
			}
			else {
				if (verbose) IO.pl("\tBKAF Check\tNA\tNo BKAFs");
			}
		}
		
		//check min BKZ score
		if (minimumBKZ !=0){
			if (bkzString != null){
				double bkz = Double.parseDouble(bkzString);
				boolean pass = bkz >= minimumBKZ;
				if (verbose) IO.pl("\tBKZ Check\t"+pass+"\t"+bkz);
				//prior not scored?
				if (foundBK == false && pass == true) passBKs = true;
				//prior scored but this failed
				if (foundBK == true && pass == false) passBKs = false;
				foundBK = true;
			}
			else {
				if (verbose) IO.pl("\tBKZ Check\tNA\tNo BKZ");
			}
		}
		if (verbose) IO.pl("\tBKZ Checks\t"+foundBK+"\t"+passBKs);
		if (foundBK) numWithBKAFs++;
		if (passBKs) numPassingBKAFs++;
		return new boolean[]{foundBK, passBKs};
	}

	/*Returns null if no clinvar annotation, returns boolean[]{foundRequestedForClinvar, foundRequestedAgainstClinvar}*/
	private boolean[] checkClinSig() throws Exception {
		String csAll = infoKeyValue.get("CLNSIG");
		if (csAll != null){
			//record findings
			numWithClin++;
			dataLine.clinSig = csAll;
			dataLine.clinSigCLNHGVS= infoKeyValue.get("CLNHGVS");
			dataLine.clinSigConf = infoKeyValue.get("CLNSIGCONF");
			dataLine.clinAlleleId = infoKeyValue.get("ALLELEID");
			if (clinsig.containsKey(csAll)) {
				int count = clinsig.get(csAll);
				clinsig.put(csAll, new Integer( ++count ));
			}
			else clinsig.put(csAll, new Integer(1));
			String[] cTerms = COMMA_SLASH.split(csAll);
			boolean passClinSig = false;
			double[] numPathoBenign = null;
			//is there a check for particular clinsig terms?
			if (passingClinSig != null){
				for (String cs: cTerms){
					String cslc = cs.toLowerCase();
					passClinSig = passingClinSig.contains(cslc);
					//Conflicting_interpretations_of_pathogenicity? Check fraction patho
					if (passClinSig && minFractionPathogenic != 0 && cslc.equals("conflicting_interpretations_of_pathogenicity")) {
						numPathoBenign = calcFracPathoFromConfInterp();
						double total = (numPathoBenign[0]+numPathoBenign[1]);
						double fracPatho = 0;
						if (total != 0) fracPatho = numPathoBenign[0]/total;
						passClinSig = (fracPatho >= minFractionPathogenic);
					}
					//drug_response gene name limits?
					if (passClinSig == true && cslc.equals("drug_response") && drugResClinSigGenes != null) {
						boolean found = false;
						for (AnnotatedGene gene: dataLine.annoGenes) {
							if (drugResClinSigGenes.contains(gene.geneName)) {
								found = true;
								break;
							}
						}
						passClinSig = found;
					}
					if (verbose) IO.pl("\tCLINSIG For Check\t"+passClinSig+"\t"+cs);
					if (passClinSig) {
						dataLine.passesCLINVAR = true;
						dataLine.rejectedCLINVAR = false;
						numPassingClinSing++;
						break;
					}
				}
			}
			
			boolean passClinSigAgainst = false;
			//is there a check against particular clinsig terms?
			if (excludeClinSig != null){
				for (String cs: cTerms){
					passClinSigAgainst = excludeClinSig.contains(cs.toLowerCase());
					if (verbose) IO.pl("\tCLINSIG Against Check\t"+passClinSigAgainst+"\t"+cs);
					if (passClinSigAgainst) {
						numPassingExcludeClinSing++;
						dataLine.rejectedCLINVAR = true;
						dataLine.passesCLINVAR = false;
						break;
					}
				}
				//check if just benign in Conflicting_interpretations_of_pathogenicity
				if (passClinSigAgainst == false && numPathoBenign !=null) {
					if (numPathoBenign[0]==0 && numPathoBenign[1]>0) {
						passClinSigAgainst = true;
						numPassingExcludeClinSing++;
						dataLine.rejectedCLINVAR = true;
						dataLine.passesCLINVAR = false;
					}
				}
			}
			return new boolean[]{passClinSig, passClinSigAgainst};
		}
		else {
			if (verbose) IO.pl("\tCLINSIG Check\tfalse\tNo CLINSIG");
			return null;
		}
	}
	
	/**
	 * ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
	 * 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
	 *    0          1                2               3          4            5
	 * ANN=GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20),
	 * GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20)
	 * 
	 * */
	private boolean checkImpact() throws Exception {
		//find ANN=
		String annValue = infoKeyValue.get("ANN");
		if (annValue == null) {
			if (verbose) IO.pl("\tImpact Check\tfalse\tNo ANN");
			return false;
		}
		numWithAnn++;

		//split ANN by comma
		String[] splitAnns = Misc.COMMA.split(annValue);
		
		//for each annotation
		for (String ann: splitAnns){
			//split ANN and save impact
			String[] splitAnn = Misc.PIPE.split(ann);
			String impact = splitAnn[2];
			if (impacts.containsKey(impact)) impacts.put(impact, new Integer( impacts.get(impact)+1 ));
			else impacts.put(impact, new Integer(1));

		}

		//check if impact is one of the ones they want
		if (passingAnnImpact!=null){
			//for each annotation
			
			for (String ann: splitAnns){
				//skip problematic transcripts?
				if (skipWarningTrans == true && ann.contains("WARNING_TRANSCRIPT")== true) continue;
				//protein transcript?
				if (onlyProteinCoding == true && ann.contains("protein_coding") == false) continue;
				
				String[] splitAnn = Misc.PIPE.split(ann);
				String impact = splitAnn[2];
				boolean passImpact = passingAnnImpact.contains(impact.toLowerCase());
				
				//add check for frameshift check?
				if (passImpact == true && justFrameShiftStartStop) {
					String annLC = splitAnn[1].toLowerCase();
					if (annLC.contains("frameshift_variant")==false && 
							annLC.contains("stop_gained")==false &&
							annLC.contains("stop_lost")==false &&
							annLC.contains("start_lost")==false ) passImpact = false;
				}
				
				if (verbose) IO.pl("\tImpact Check\t"+passImpact+"\t"+impact);
				dataLine.annoGenes.add(new AnnotatedGene(passImpact, splitAnn));
				
				if (passImpact) {
					numPassingAnnImpact++;
					dataLine.passesANN = true;
					//Annotation effect
					String effect = splitAnn[1];
					if (effects.containsKey(effect)) effects.put(effect, new Integer( effects.get(effect)+1 ));
					else effects.put(effect, new Integer(1));
					return true;
				}
			}
			//none found so return false
			return false;
		}
		else return false;
	}
	

	
	/**
	 * ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
	 * 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
	 *    0          1                2               3          4            5
	 * ANN=GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20),
	 * GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20)
	 * 
	 * */
	private boolean checkGermlineGenes() throws Exception {
		//find ANN=
		String annValue = infoKeyValue.get("ANN");
		if (annValue == null) {
			if (verbose) IO.pl("\tGermline Check\tNA\tNo ANN");
			return false;
		}

		//split ANN by comma
		String[] splitAnns = Misc.COMMA.split(annValue);

		//Save observed genes, but don't double count!
		HashSet<String> genes = new HashSet<String>();
		for (String ann: splitAnns){
			//split ANN and save gene name
			String[] splitAnn = Misc.PIPE.split(ann);
			genes.add(splitAnn[3]);
		}
		Iterator<String> it = genes.iterator();
		String foundGermlineGene = null;
		while (it.hasNext()){
			String geneName = it.next();
			if (passingGermlineGenes.contains(geneName)) {
				foundGermlineGene = geneName;
				if (observedGermlineGenes.containsKey(geneName)) observedGermlineGenes.put(geneName, new Integer( observedGermlineGenes.get(geneName)+1 ));
				else observedGermlineGenes.put(geneName, new Integer(1));
			}
		}
		
		//for each annotation
		if (foundGermlineGene!= null){
			if (verbose) IO.pl("\tGermGene Check\ttrue\t"+foundGermlineGene);
			return true;
		}
		//none found so return false
		if (verbose) IO.pl("\tGermGene Check\tNA\t"+genes);
		return false;
	}

	
	private boolean[] checkPopFreq() throws Exception {
		boolean foundPop = false;
		boolean passDbNSFPExAC = true;
		boolean passExAC_AF = true;
		boolean pass1K = true;
		double maxPopFreq = -1;
		
		//ExAC from dbNSFP
		String dbExacString = infoKeyValue.get("dbNSFP_ExAC_AF");
		if (dbExacString!= null){
			if (dbExacString.equals(".") == false){
				foundPop = true;
				double popAF = Double.parseDouble(dbExacString);
				if (popAF > maximumPopAF) passDbNSFPExAC = false;
				if (verbose) IO.pl("\tExAC Check\t"+passDbNSFPExAC+"\t"+popAF);
				if (popAF> maxPopFreq) maxPopFreq = popAF;
			}
		}
		//ExAC from just ExAC
		String exacAFString = infoKeyValue.get("ExAC_AF");
		if (exacAFString != null){
			if (exacAFString.equals(".") == false){
				foundPop = true;
				double popAF = Double.parseDouble(exacAFString);
				if (popAF > maximumPopAF) passExAC_AF = false;
				if (verbose) IO.pl("\tExAC AF Check\t"+passExAC_AF+"\t"+popAF);
				if (popAF> maxPopFreq) maxPopFreq = popAF;
			}
		}
		//1K genomes
		String KString = infoKeyValue.get("dbNSFP_1000Gp3_AF");
		if (KString != null){
			if (KString.equals(".") == false){
				foundPop = true;
				double popAF = Double.parseDouble(KString);
				if (popAF > maximumPopAF) pass1K = false;
				if (verbose) IO.pl("\t1KG Check\t"+pass1K+"\t"+popAF);
				if (popAF> maxPopFreq) maxPopFreq = popAF;
			}
		}
		if (foundPop) {
			dataLine.varPopAlleleFreq = maxPopFreq;
			numWithPopAF++;
			if (passDbNSFPExAC == true && pass1K == true && passExAC_AF) {
				numPassingPopAF++;
				return new boolean[]{true, true};
				
			}
			else return new boolean[]{true, false};
		}
		else {
			if (verbose) IO.pl("\tPopFreq Check\tNA\tNo dbNSFP_1000Gp3_AF, dbNSFP_ExAC_AF, or ExAC_AF");
			return new boolean[]{false, false};
		}
	}


	private boolean checkDP() throws Exception {
		String dpString = infoKeyValue.get("T_DP");
		if (dpString == null) dpString = infoKeyValue.get("DP");
		if (dpString != null){
			int obsDP = Integer.parseInt(dpString);
			if (obsDP < dps.length) dps[obsDP]++;
			dataLine.totalUniObDepth = obsDP;
			boolean passDP = obsDP >= minimumDP;
			if (verbose) IO.pl("\tDP Check\t"+passDP+"\t"+obsDP);
			if (passDP == false) return false;
			else {
				numPassingDP++;
				return true;
			}
		}
		else {
			if (verbose) IO.pl("\tDP Check\tfalse\tNo DP");
			return false;
		}
	}
	
	private boolean checkAF(boolean passGermlineGene) {
		try {
			String afString = infoKeyValue.get("T_AF");
			if (afString == null) afString = infoKeyValue.get("AF");
			if (afString != null){
				obsAF = Double.parseDouble(afString);
				afs.count(obsAF);
				dataLine.varAlleleFreq = obsAF;
				boolean passMinAF = obsAF >= minimumAF;
				boolean passMaxAF = obsAF <= maximumAF;
				if (verbose) {
					IO.pl("\tAF Min Check\t"+passMinAF+"\t"+obsAF);
					if( passGermlineGene == false) IO.pl("\tAF Max Check\t"+passMaxAF+"\t"+obsAF);
				}
				if (passMinAF) numPassingMinAF++;
				if (passMaxAF) numPassingMaxAF++;
				if (passGermlineGene) passMaxAF = true;
				if (passMinAF && passMaxAF) return true;
				else return false;
			}
			else {
				if (verbose) IO.pl("\tAF Check\tfalse\tNo AF");
				return false;
			}
		} catch (NumberFormatException e){
			if (verbose) IO.pl("\tAF Check\tfalse\tCan't parse AF");
			return false;
		}
	}
	
	private boolean checkAlt(boolean passGermlineGene) throws Exception {
		//look for db in data line
		if (dataLine.totalUniObDepth == -1) {
			checkDP();
			if (dataLine.totalUniObDepth == -1) {
				if (verbose) IO.pl("\tAlt Check\tfalse\tNo DP");
				return false;
			}
		}
		//look for af
		if (dataLine.varAlleleFreq == -1) {
			checkAF(passGermlineGene);
			if (dataLine.varAlleleFreq == -1) {
				if (verbose) IO.pl("\tAlt Check\tfalse\tNo AF");
				return false;
			}
		}
		
		//calculate it
		dataLine.varUniOb = (int)Math.round(dataLine.varAlleleFreq * (double)dataLine.totalUniObDepth);
		boolean passAlt = dataLine.varUniOb >= minimumAlt;
		if (verbose) IO.pl("\tAlt Check\t"+passAlt+"\t"+dataLine.varUniOb);
		if (passAlt == false) return false;
		else {
			numPassingAlt++;
			return true;
		}
	}
	
	private boolean[] checkCF() throws Exception {
		boolean foundCF = false;
		boolean passCF = false;
		
		String cfString = infoKeyValue.get("CF");
		if (cfString != null) {
			foundCF = true;
			String[] fields = Misc.COMMA.split(cfString);
			dataLine.priorCallFreq = cfString;
			double cf = Double.parseDouble(fields[0]);
			if (cf <= maximumCF) passCF = true;
			else passCF = false;
			
		}
		if (verbose) IO.pl("\tCF Check\t"+foundCF+"\t"+passCF);
		
		if (foundCF) numWithCF++;
		if (passCF) numPassingCF++;
		return new boolean[]{foundCF, passCF};
	}

	public Histogram makeReadDepthHistogram(){
		//find max DP
		int maxDP = 0;
		for (int i=dps.length-1; i>=0; i--){
			if (dps[i] !=0) {
				maxDP = i;
				break;
			}
		}
		//find min DP
		int minDP = 0;
		for (int i=0; i< dps.length; i++){
			if (dps[i] !=0){
				minDP = i;
				break;
			}
		}
		Histogram hist = new Histogram(minDP, maxDP, 100);
		for (int i=minDP; i<=maxDP; i++ ){
			for (int x=0; x < dps[i]; x++) hist.count(i);
		}
		return hist;
	}
	
	private void modifySettingsForFoundation() {
		if (minimumDP == 0) minimumDP = 50;
		if (minimumAlt == 0) minimumAlt = 3;
		if (minimumAF == 0) minimumAF = 0.01;
		if (maximumCF == 1.0) maximumCF = 0.1;
		if (maxFracBKAFs == 0) maxFracBKAFs = 0.1;
		if (maximumPopAF == 0) maximumPopAF = 0.01;
		if (passingAnnImpact == null) {
			passingAnnImpact = new TreeSet<String>();
			passingAnnImpact.add("high");
			passingAnnImpact.add("moderate");
		}
		if (passingClinSig == null){
			passingClinSig = new TreeSet<String>();
			passingClinSig.add("pathogenic");
			passingClinSig.add("conflicting_interpretations_of_pathogenicity");
			passingClinSig.add("likely_pathogenic");
			passingClinSig.add("drug_response");
			passingClinSig.add("risk_factor");
		}
		if (excludeClinSig == null){
			excludeClinSig = new TreeSet<String>();
			excludeClinSig.add("benign");
			excludeClinSig.add("likely_benign");
			
		}
		if (drugResClinSigGenes == null){
			drugResClinSigGenes = new TreeSet<String>();
			drugResClinSigGenes.add("RYR1");
		}
		if (passingIDKeys == null) passingIDKeys = new String[]{"foundation","tempus","caris"};
		if (passingVCFSS == null) passingVCFSS = new String[]{"D5S", "D3S", "G5S", "G3S"};	
		if (minimumVCFSSDiff == 0) minimumVCFSSDiff = 4;
		createGermlineGenes();
	}

	public void printSettings(){
		ArrayList<String> al = new ArrayList<String>();
		
		al.add("Settings:");
		if (somaticProcessing) {
			al.add("\tSomatic processing, unset fields are set to Somatic defaults and the following logic is used to pass or fail each record.\n"
					+ "\t\tIf ID matches pass irregardless of other settings.\n"
					+ "\t\tFail it if DP or AF are false.\n"
					+ "\t\tFail it if ClinSig, Splice, and Impact are all false.\n"
					+ "\t\tFail it if Pop, BKAF, CF are present and they are false.\n");
		}
		if (minimumDP != 0) al.add("\t"+ minimumDP+"\t: Minimum DP read depth");
		if (minimumAlt != 0) al.add("\t"+ minimumAlt+"\t: Minimum ALT observations");
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) {
			al.add("\t"+ minimumAF+"\t: Minimum AF allele frequency");
			al.add("\t"+ maximumAF+"\t: Maximum AF allele frequency");
		}
		if (maximumCF != 1) al.add("\t"+ maximumCF+"\t: Maximum CF prior call frequency");
		if (maximumPopAF != 0) al.add("\t"+ maximumPopAF+"\t: Maximum population AF from dbNSFP_ExAC_AF or dbNSFP_1000Gp3_AF");
		if (maxFracBKAFs != 0) al.add("\t"+ maxFracBKAFs+"\t: Maximum fraction of background samples with an AF >= observed AF");
		if (minimumBKZ != 0) al.add("\t"+ minimumBKZ+"\t: Minimum BKZ quality score");

		al.add("\t"+ skipWarningTrans+"\t: Ignore transcripts labeled WARNING_TRANSCRIPT_XXX");
		al.add("\t"+ onlyProteinCoding+"\t: Only consider protein_coding transcripts");
		if (transcriptFilter != null) al.add("\t"+ transcriptList+"\t: Select annotations that match transcript IDs in this file.");
		
		if (passingAnnImpact != null) {
			al.add("\n\t"+Misc.treeSetToString(passingAnnImpact, ",")+"\t: ANN impact keys");
			al.add("\t"+ justFrameShiftStartStop+"\t: Further restrict ANN impacts to one of these effects: frameshift_variant, stop_gained, stop_lost, or start_lost.");
		}
		if (passingClinSig != null) {
			al.add("\n\t"+Misc.treeSetToString(passingClinSig, ",")+"\t: CLINSIG keep keys");
			if (minFractionPathogenic!=0 && passingClinSig.contains("conflicting_interpretations_of_pathogenicity")) {
				al.add("\t"+minFractionPathogenic+"\t: Minimum fraction pathogenic when CLINSIG is Conflicting_interpretations_of_pathogenicity");
			}
		}
		if (excludeClinSig != null) al.add("\t"+Misc.treeSetToString(excludeClinSig, ",")+"\t: CLINSIG exclude keys");
		if (drugResClinSigGenes != null) al.add("\t"+Misc.treeSetToString(drugResClinSigGenes, ",")+"\t: CLINSIG drug_response restricted genes");
		if (clinvarDate != null) al.add("\t"+ clinvarDate+"\t: CLINVAR file date for spreadsheet.");
		
		if (passingVCFSS != null) {
			al.add("\n\t"+Misc.stringArrayToString(passingVCFSS, ",")+"\t: Splice junction types to scan");
			al.add("\t"+ minimumVCFSSDiff+"\t: Minimum difference in MaxEnt scan scores for a splice junction effect");
		}
		if (somaticProcessing ==false) al.add("\t"+ orAnnos+"\t: Require that only one need pass: ANN Impact or Clinvar or Splice Effect");
		if (passingIDKeys != null) al.add("\n\t"+Misc.stringArrayToString(passingIDKeys, ",")+"\t: VCF ID keys that if present pass it regardless of any filter settings");
		if (passingGermlineGenes != null) al.add("\n\t"+acmgGenes+"\t: ACMG Genes that cause the maxAF filter to be skipped for intersecting variants");
		if (somaticProcessing ==false) 
		al.add("\n\t"+verbose+"\t: Verbose output");
		al.add("\t"+ saveDirectory+"\t: Save directory");
		
		appSettings = Misc.stringArrayListToString(al, "\n");
		IO.pl("\n"+appSettings);
		
	}
	
	public String format(double numerator, double denominator){
		String frac = null;
		if (denominator == 0.0) frac = "0";
		else frac =  Num.formatNumber(100*numerator/denominator, 3);
		return frac +"%\t" +(int)numerator+"/"+(int)denominator+"\t";
	}
	
	public void printStats(){
		IO.pl("\nProcessing statistics for "+vcfFiles.length+" files:");
		IO.pl("\t"+ numRecords+"\tNumber of VCF Records");
		if (minimumDP != 0) IO.pl("\t"+ format(numPassingDP, numRecords)+"Passing DP");
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) {
			IO.pl("\t"+ format(numPassingMinAF, numRecords)+"Passing MinAF");
			IO.pl("\t"+ format(numPassingMaxAF, numRecords)+"Passing MaxAF");
		}
		if (minimumAlt != 0) IO.pl("\t"+ format(numPassingAlt, numRecords)+"Passing ALT Obs");
		if (maximumCF != 1){
			IO.pl("\t"+ format(numWithCF, numRecords)+"With CF");
			IO.pl("\t"+ format(numPassingCF, numWithCF)+"Passing CF");
		}
		if (maximumPopAF != 0){
			IO.pl("\t"+ format(numWithPopAF, numRecords)+"With Pop AF");
			IO.pl("\t"+ format(numPassingPopAF, numWithPopAF)+"Passing Pop AF");
		}
		if (passingAnnImpact != null) {
			IO.pl("\t"+ format(numWithAnn, numRecords)+"With Ann");
			IO.pl("\t"+ format(numPassingAnnImpact, numWithAnn)+"Passing Ann Impact");
		}
		if (maxFracBKAFs != 0) {
			IO.pl("\t"+ format(numWithBKAFs, numRecords)+"With BKAFs BKZs");
			IO.pl("\t"+ format(numPassingBKAFs, numWithBKAFs)+"Passing BKAFs BKZs");
		}
		if (passingClinSig != null){
			IO.pl("\t"+ format(numWithClin, numRecords)+"With Clin");
			IO.pl("\t"+ format(numPassingClinSing, numWithClin)+"Passing Include ClinSig");
			IO.pl("\t"+ format(numPassingExcludeClinSing, numWithClin)+"Passing Exclude ClinSig");
		}
		if (passingVCFSS != null){
			IO.pl("\t"+ format(numWithSplice, numRecords)+"With Splice");
			IO.pl("\t"+ format(numPassingSplice, numWithSplice)+"Passing Splice");
		}
		if (passingIDKeys != null) IO.pl("\t"+ format(numPassIDs, numRecords)+"Pass IDs");
		IO.pl("\t"+ format(numPass, numRecords)+"Records passing all filters");
	}
	
	public void printDistributions(){
		if (passingClinSig != null){
			IO.pl("\nObserved CLINSIG annotations:" );
			for (String key: clinsig.keySet()) IO.pl("\t"+key+"\t"+clinsig.get(key));
		}
		if (passingAnnImpact != null){
			IO.pl("\nObserved ANN impacts:");
			for (String key: impacts.keySet()) IO.pl("\t"+key+"\t"+impacts.get(key));
		}
		if (true){
			IO.pl("\nObserved ANN effects for passing Impacts:");
			for (String key: effects.keySet()) IO.pl("\t"+key+"\t"+effects.get(key));
		}
		if (passingGermlineGenes != null){
			IO.pl("\nObserved ACMG germline genes:");
			for (String key: observedGermlineGenes.keySet()) IO.pl("\t"+key+"\t"+observedGermlineGenes.get(key));
		}
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0){
			IO.pl("\nAF Allele Frequency distribution:");
			afs.printScaledHistogram();
		}
		if (minimumDP != 0) {
			IO.pl("\nDP Read Depth distribution:");
			makeReadDepthHistogram().printScaledHistogram();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AnnotatedVcfParser(args);
	}
	
	public static String[] appendConfigArgs(String[] args, String configArg){
		for (int i=0; i< args.length; i++){
			if (args[i].equals(configArg)){
				File config = new File(args[++i]);
				String[] lines = IO.loadFile(config);
				if (lines == null || lines.length == 0) Misc.printErrAndExit("\nProblem with your config file, this should contain one or more lines of parameteres, none should contain whitespace "+config);
				ArrayList<String> otherArgs = new ArrayList<String>();
				for (String line: lines){
					String[] t = Misc.WHITESPACE.split(line);
					for (String a: t) otherArgs.add(a);
				}
				return Misc.copyAndMerge(args, Misc.stringArrayListToStringArray(otherArgs));
			}
		}
		return args;
	}


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		//look for a config file and append
		args = appendConfigArgs(args,"-y");
		Pattern pat = Pattern.compile("-[a-zA-T]");
		File forExtraction = null;
		String impactString = null;
		String clinString = null;
		String clinStringExclude = null;
		String drugResString = null;
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i];
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'T': transcriptList = new File(args[++i]); break;
					case 'd': minimumDP = Integer.parseInt(args[++i]); break;
					case 'w': minimumAlt = Integer.parseInt(args[++i]); break;
					case 'm': minimumAF = Double.parseDouble(args[++i]); break;
					case 'x': maximumAF = Double.parseDouble(args[++i]); break;
					case 'p': maximumPopAF = Double.parseDouble(args[++i]); break;
					case 'q': maximumCF = Double.parseDouble(args[++i]); break;
					case 'b': maxFracBKAFs = Double.parseDouble(args[++i]); break;
					case 'g': passingVCFSS = Misc.COMMA.split(args[++i]); break;
					case 'n': minimumVCFSSDiff = Double.parseDouble(args[++i]); break;
					case 'z': minimumBKZ = Double.parseDouble(args[++i]); break;
					case 't': minFractionPathogenic = Double.parseDouble(args[++i]); break;
					case 'u': drugResString = args[++i]; break;
					case 'a': impactString = args[++i]; break;
					case 'c': clinString = args[++i]; break;
					case 'C': clinvarDate = args[++i]; break;
					case 'e': clinStringExclude = args[++i]; break;
					case 'j': createGermlineGenes(); break;
					case 'i': passingIDKeys = Misc.COMMA.split(args[++i]); break;
					case 'o': orAnnos = true; break;
					case 'r': verbose = true; break;
					case 'k': onlyProteinCoding = false; break;
					case 'f': somaticProcessing = true; break;
					case 'l': justFrameShiftStartStop = true; break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'y': break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");
		
		if (saveDirectory != null){
			saveDirectory.mkdirs();
			if (saveDirectory.isDirectory() == false || saveDirectory.exists() == false) Misc.printErrAndExit("\nCannot find your save directory?! "+saveDirectory);
		}
		else Misc.printErrAndExit("\nPlease provide a save directory.");
		
		if (impactString != null) {
			impactString = impactString.toLowerCase();
			passingAnnImpact = new TreeSet<String>();
			for (String s: Misc.COMMA.split(impactString)) passingAnnImpact.add(s);
		}
		
		if (clinString != null) {
			clinString = clinString.toLowerCase();
			passingClinSig = new TreeSet<String>();
			for (String s: Misc.COMMA.split(clinString)) passingClinSig.add(s);
		}
		if (clinStringExclude != null) {
			clinStringExclude = clinStringExclude.toLowerCase();
			excludeClinSig = new TreeSet<String>();
			for (String s: Misc.COMMA.split(clinStringExclude)) excludeClinSig.add(s);
		}
		if (drugResString != null) {
			drugResClinSigGenes = new TreeSet<String>();
			for (String s: Misc.COMMA.split(drugResString)) drugResClinSigGenes.add(s);
		}
		if (minimumVCFSSDiff !=0 && passingVCFSS == null) Misc.printErrAndExit("\nError: please provide a comma delimited list of splice junction types (e.g. D5S,D3S,G5S,G3S) to examine.\n");
		
		if (transcriptList != null)  transcriptFilter = IO.loadFileIntoHashSet(transcriptList);
		
	}
	

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                            Annotated Vcf Parser  Jan 2022                        **\n" +
				"**************************************************************************************\n" +
				"Splits VCF files that have been annotated with SnpEff, ExAC, and clinvar, plus the \n"+
				"VCFBkz, VCFCallFrequency, and VCFSpliceScanner USeq apps into passing and failing\n"+
				"records. Use the -r option to inspect the effect of the various filters on each\n"+
				"record. Use the VCFRegionFilter app to restrict variants to particular regions.\n"+
				"A summary spreadsheet is exported with select information and excel hyperlinks for\n"+
				"rapid inspection. Must have AF (or T_AF) and DP in the INFO field for Alt filtering.\n"+

				"\nOptions:\n"+
				"-v File path or directory containing xxx.vcf(.gz/.zip OK) file(s) to filter.\n" +
				"-s Directory for saving the results.\n"+
				"-f Perform a candidate somatic variant processing. Setting the following overrides\n"+
				"        the defaults. \n"+
				"-d Minimum DP alignment depth\n"+
				"-m Minimum AF allele frequency\n"+
				"-x Maximum AF allele frequency\n"+
				"-w Minimum Alt observations\n"+
				"-q Maximum CF prior call frequency\n"+
				"-k Include non protein_coding transcripts in ANN, defaults to excluding\n"+
				"-j Ignore the max AF filter for ACMG incidental germline gene variants.\n"+
                "-p Maximum population allele frequency, only applies if present.\n"+
				"-z Minimum BKZ when present.\n"+
                "-b Maximum fraction of BKAF samples with allele frequency >= VCF AF, only applies\n"+
                "       if present.\n"+
                "-g Splice junction types to scan.\n"+
                "-n Minimum difference in splice junction scores, only applies if present.\n"+
                "-a Comma delimited list of SnpEff ANN impact categories to select for.\n"+
                "-l Further restrict ANN impacts to one of these effects: frameshift_variant,\n"+
                "       stop_gained, stop_lost, or start_lost.\n"+
                "-c Comma delimited list of CLINSIG terms to select for.\n"+
                "-t For CLINSIG=Conflicting_interpretations_of_pathogenicity, minimum fraction\n"+
                "       pathogenicity, defaults to 0\n"+
                "-u For CLINSIG=drug_response, only pass it if associated with these genes. Comma\n"+
                "       delimited, no spaces. Defaults to all.\n"+
                "-e Comma delimited list of CLINSIG terms to select against.\n"+
                "-C CLINVAR file date for spreadsheet output.\n"+

                "-i Comma delimited list of VCF ID keys to select for. If the VCF ID contains one or\n"+
                "       more, the record is passed regardless of other filters. The match is not exact.\n"+
                "-o Only require, if set or present, SnpEff ANN or CLINSIG or Splice to be true to pass.\n"+
                "       Defaults to require that all set pass.\n"+
                "-T Txt file containing a list of transcripts to keep, one per line, to use in filtering\n"+
                "       the summary spreadsheet output.\n "+
				"-r Verbose per record output.\n"+ 
                "-y Path to a config txt file for setting the above.\n"+

                
				
				"\nExample: java -jar pathToUSeq/Apps/AnnotatedVcfParser -v VCFFiles/ -s Parsed/\n"+
				"        -d 74 -m 0.05 -x 0.75 -j -p 0.02 -b 0.1 -z 3 -g D5S,D3S,G5S,G3S -n 4.5 -a\n"+
				"        HIGH,MODERATE -e Benign,Likely_benign -c Pathogenic,Likely_pathogenic,\n"+
				"        Conflicting_interpretations_of_pathogenicity -t 0.51 -u RYR1 -T \n"+
				"        ~/Ref/ACMGTranscripts.txt -C 2021-04-04 -w 3\n\n"+


				"**************************************************************************************\n");
	}
}
