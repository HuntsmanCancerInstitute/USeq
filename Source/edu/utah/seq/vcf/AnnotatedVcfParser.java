package edu.utah.seq.vcf; 

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
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
	private int minimumDP = 0;
	private double minimumAF = 0;
	private double maximumAF = 1;
	private double maximumPopAF = 0;
	private double maxFracBKAFs = 0; //maximum number of background samples with >= AF, measure as fraction of total samples, 50 for foundation so 5 max
	private double minimumBKZ = 0;
	private String[] passingIDKeys = null;
	private TreeSet<String> passingAnnImpact = null;
	private TreeSet<String> passingAnnEffect = new TreeSet<String>();
	private TreeSet<String> passingClinSig = null;
	private TreeSet<String> excludeClinSig = null;
	private TreeSet<String> passingGermlineGenes = null;
	private String[] passingVCFSS = null;
	private double minimumVCFSSDiff = 0; //difference between ref and alt for splice junction score
	private boolean verbose = false;
	private boolean somaticProcessing = false;
	private boolean orAnnos = false;
	
	//Patterns
	private static final Pattern AF = Pattern.compile("AF=([\\d\\.]+);");
	private static final Pattern DP = Pattern.compile("DP=(\\d+);");
	private static final Pattern EXAC = Pattern.compile("dbNSFP_ExAC_AF=([\\d\\.e-]+);");
	private static final Pattern G1000 = Pattern.compile("dbNSFP_1000Gp3_AF=([\\d\\.e-]+);");
	private static final Pattern CLINSIG = Pattern.compile("CLNSIG=([\\w/,]+);");
	private static final Pattern BKAF = Pattern.compile(";*BKAF=([\\d\\.,]+);");
	private static final Pattern BKZ = Pattern.compile("BKZ=([\\d\\.]+);");
	private static final Pattern COMMA_SLASH = Pattern.compile(",|/");
	
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
	private int numPassingExcludeClinSing = 0;
	private int numWithSplice = 0;
	private int numPassingSplice = 0;
	private int numPassIDs = 0;
	private int numPass = 0;
	
	//working fields
	private String vcfLine = null;
	private double obsAF;
		
	public AnnotatedVcfParser (String[] args) {
		long startTime = System.currentTimeMillis();
		processArgs(args);
		if (somaticProcessing) modifySettingsForFoundation();
		IO.pl("Parsing (#pass #fail):");
		for (File vcf: vcfFiles){
			System.out.print("\t"+vcf.getName());
			if (verbose) IO.pl();
			parse(vcf);
		}
		printSettings();
		printStats();
		printDistributions();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" sec\n");
	}
	
	
	public void parse(File vcf) {
		Gzipper passVcf = null;
		Gzipper failVcf = null;
		try {
			//IO
			String name = Misc.removeExtension(vcf.getName());			
			passVcf = new Gzipper( new File (saveDirectory, name + "_Pass.vcf.gz"));
			failVcf = new Gzipper( new File (saveDirectory, name + "_Fail.vcf.gz"));
			BufferedReader in = IO.fetchBufferedReader(vcf);
			int numPassingVcf = 0;
			int numFailingVcf = 0;
			//for each line in the file
			while ((vcfLine = in.readLine()) != null){
				vcfLine = vcfLine.trim();
				if(vcfLine.length() == 0) continue;
				//header? just print out
				if (vcfLine.startsWith("#")) {
					passVcf.println(vcfLine);
					failVcf.println(vcfLine);
					continue;
				}
				numRecords++;
				if (verbose) IO.pl("\n"+vcfLine);

				//#CHROM POS ID REF ALT QUAL FILTER INFO ......
				//   0    1   2  3   4   5     6      7
				String[] cells = Misc.TAB.split(vcfLine);
				String[] info = Misc.SEMI_COLON.split(cells[7]);
				
				boolean passDP=true, passAF=true, passImpact=true, passSplice=true, passID = true, passGermlineGenes = false;
				boolean[] passBKAF= null; 
				boolean[] passPop= null;
				boolean[] passClinSig= null;
				
				if (passingGermlineGenes != null) passGermlineGenes = checkGermlineGenes(info);
				if (minimumDP != 0) 			passDP = checkDP(cells[7]);
				if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) passAF = checkAF(cells[7], passGermlineGenes);
				if (maximumPopAF != 0) 			passPop = checkPopFreq(cells[7]);
				if (passingAnnImpact != null) 	passImpact = checkImpact(info);
				if (passingClinSig != null || excludeClinSig != null) 	passClinSig = checkClinSig(cells[7]);
				if (maxFracBKAFs != 0 || minimumBKZ !=0) passBKAF = checkBKAFs(cells[7]);
				if (passingVCFSS != null) 		passSplice = checkSplice(cells[7], info);
				if (passingIDKeys != null) 		passID = checkIDs(cells[2]);
				
				boolean pass;
				if (somaticProcessing) pass =  passFoundationFiltering(passID, passDP, passAF, passImpact,passSplice, passPop, passBKAF,passClinSig); 
				else if (orAnnos) pass = passWithOrAnnos(passID, passDP, passAF, passImpact, passSplice, passPop, passBKAF,passClinSig);
				else pass = allPass(passID, passDP, passAF, passImpact, passSplice, passPop, passBKAF, passClinSig);
				
				if (pass){
					numPassingVcf++;
					passVcf.println(vcfLine);
				}
				else {
					failVcf.println(vcfLine);
					numFailingVcf++;
				}
				if (verbose) IO.pl("\tPass all filters\t"+pass);	
			}
			numPass+=numPassingVcf;
			
			if (verbose == false) IO.pl("\t"+numPassingVcf+"\t"+numFailingVcf);
			
			//close io
			in.close();
			passVcf.close();
			failVcf.close();
			
		} catch (Exception e) {
			System.err.println("\nERROR: parsing vcf file: "+vcf.getName()+"\nLine: "+vcfLine+"\n");
			e.printStackTrace();
			if (passVcf!= null) {
				passVcf.getGzipFile().delete();
				failVcf.getGzipFile().delete();
			}
			System.exit(1);
		} 
	}
	
	/**Requires all are true. Not if these aren't set off defaults then they default to true. Note if passID then returns true regardless of others.*/
	private boolean allPass(boolean passID, boolean passDP, boolean passAF, boolean passImpact, boolean passSplice, boolean[] passPop, boolean[] passBKAF, boolean[] passClinSig) {
		if (passingIDKeys != null && passID) return true;
		if (passDP == false) return false;
		if (passAF == false) return false;
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
	
	/**Requires DP, AF, Pop, BKAF are true if present. Only one of Clin, Splice, or Impact need be true, Note if passID then returns true regardless of others.*/
	private boolean passWithOrAnnos(boolean passID, boolean passDP, boolean passAF, boolean passImpact, boolean passSplice, boolean[] passPop, boolean[] passBKAF, boolean[] passClinSig) {
		if (passingIDKeys != null && passID) return true;
		if (passDP == false) return false;
		if (passAF == false) return false;
		boolean passBKAF2 = true;
		if (passBKAF != null && passBKAF[0] == true) passBKAF2 = passBKAF[1];
		if (passBKAF2 == false) return false;
		boolean passPop2 = true;
		if (passPop != null && passPop[0] == true) passPop2 = passPop[1];
		if (passPop2 == false) return false;
		//check that one is good
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

	private boolean passFoundationFiltering(boolean passID, boolean passDP, boolean passAF,  boolean passImpact, boolean passSplice, boolean[] passPop, boolean[] passBKAF, boolean[] passClinSig){
		//pass it if ID correct
		if (passID) return true;
		//fail it if DP or AF fails
		if (passAF == false || passDP == false) return false;
		//if present and fails then fail
		if (passPop[0] == true && passPop[1] == false) return false;
		if (passBKAF[0] == true && passBKAF[1] == false) return false;
		//check that one is good
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

	private boolean checkSplice(String info, String[] splitInfo) throws Exception {
		//any splice info
		if (info.contains("VCFSS=") == false){
			if (verbose) IO.pl("\tSplice Check\tfalse\tNo VCFSS");
			return false;
		}
		numWithSplice++;
		
		//find VCFSS
		String vcfSSLine = null;
		for (String s: splitInfo){
			if (s.startsWith("VCFSS=")){
				vcfSSLine = s;
				break;
			}
		}
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
					if (delta > maxDelta) maxDelta = delta;
				}
			}
		}
		
		boolean passSplice = maxDelta >= minimumVCFSSDiff;
		if (verbose) IO.pl("\tSplice Check\t"+passSplice+"\t"+maxDelta+"\t"+vcfSSLine);	
		if (passSplice) {
			numPassingSplice++;
			return true;
		}
		else return false;
	}

	
	private boolean[] checkBKAFs(String info) throws Exception {
		boolean foundBK = false;
		boolean passBKs = false;
		
		//both of the following must pass if both are being scored
				
		//check max fraction
		if (maxFracBKAFs !=0){
			Matcher mat = BKAF.matcher(info);
			if (mat.find()){
				foundBK = true;
				//record findings
				String bkafs = mat.group(1);
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
			Matcher mat = BKZ.matcher(info);
			if (mat.find()){
				double bkz = Double.parseDouble(mat.group(1));
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
	private boolean[] checkClinSig(String info) throws Exception {
		Matcher mat = CLINSIG.matcher(info);
		if (mat.find()){
			//record findings
			numWithClin++;
			String csAll = mat.group(1);
			if (clinsig.containsKey(csAll)) {
				int count = clinsig.get(csAll);
				clinsig.put(csAll, new Integer( ++count ));
			}
			else clinsig.put(csAll, new Integer(1));
			
			String[] cTerms = COMMA_SLASH.split(csAll);
			boolean passClinSig = false;
			//is there a check for particular clinsig terms?
			if (passingClinSig != null){
				for (String cs: cTerms){
					passClinSig = passingClinSig.contains(cs.toLowerCase());
					if (verbose) IO.pl("\tCLINSIG For Check\t"+passClinSig+"\t"+cs);
					if (passClinSig) {
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
						break;
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
	private boolean checkImpact(String[] splitInfo) throws Exception {
		//find ANN=
		String annValue = null;
		for (String i: splitInfo){
			if (i.startsWith("ANN=")) {
				annValue = i;
				break;
			}
		}
		if (annValue == null) {
			if (verbose) IO.pl("\tImpact Check\tfalse\tNo ANN");
			return false;
		}
		numWithAnn++;

		//drop ANN=
		annValue = annValue.substring(4);

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
				String[] splitAnn = Misc.PIPE.split(ann);
				String impact = splitAnn[2];
				boolean passImpact = passingAnnImpact.contains(impact.toLowerCase());
				if (verbose) IO.pl("\tImpact Check\t"+passImpact+"\t"+impact);
				if (passImpact) {
					numPassingAnnImpact++;
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
	private boolean checkGermlineGenes(String[] splitInfo) throws Exception {
		//find ANN=
		String annValue = null;
		for (String i: splitInfo){
			if (i.startsWith("ANN=")) {
				annValue = i;
				break;
			}
		}
		if (annValue == null) {
			if (verbose) IO.pl("\tGermline Check\tNA\tNo ANN");
			return false;
		}

		//drop ANN=
		annValue = annValue.substring(4);

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

	
	private boolean[] checkPopFreq(String info) throws Exception {
		boolean foundPop = false;
		boolean passExAC = true;
		boolean pass1K = true;
		
		//ExAC
		Matcher mat = EXAC.matcher(info);
		if (mat.find()){
			if (mat.group(1).equals(".") == false){
				foundPop = true;
				double popAF = Double.parseDouble(mat.group(1));
				if (popAF > maximumPopAF) passExAC = false;
				if (verbose) IO.pl("\tExAC Check\t"+passExAC+"\t"+popAF);
			}
		}

		//1K genomes
		mat = G1000.matcher(info);
		if (mat.find()){
			if (mat.group(1).equals(".") == false){
				foundPop = true;
				double popAF = Double.parseDouble(mat.group(1));
				if (popAF > maximumPopAF) pass1K = false;
				if (verbose) IO.pl("\t1KG Check\t"+pass1K+"\t"+popAF);
			}
		}
		
		
		if (foundPop) {
			numWithPopAF++;
			if (passExAC == true && pass1K == true) {
				numPassingPopAF++;
				return new boolean[]{true, true};
				
			}
			else return new boolean[]{true, false};
		}
		else {
			if (verbose) IO.pl("\tPopFreq Check\tNA\tNo dbNSFP_1000Gp3_AF or dbNSFP_ExAC_AF");
			return new boolean[]{false, false};
		}
		
	}


	private boolean checkDP(String info) throws Exception {
		Matcher mat = DP.matcher(info);
		if (mat.find()){
			int obsDP = Integer.parseInt(mat.group(1));
			if (obsDP < dps.length) dps[obsDP]++;
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
	
	private boolean checkAF(String info, boolean passGermlineGene) {
		try {
			Matcher mat = AF.matcher(info);
			if (mat.find()){
				obsAF = Double.parseDouble(mat.group(1));
				afs.count(obsAF);
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
		if (minimumAF == 0) minimumAF = 0.01;
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
		if (passingIDKeys == null) passingIDKeys = new String[]{"foundation"};
		if (passingVCFSS == null) passingVCFSS = new String[]{"D5S", "D3S", "G5S", "G3S"};	
		if (minimumVCFSSDiff == 0) minimumVCFSSDiff = 4;
		createGermlineGenes();
		
	}

	public void printSettings(){
		IO.pl("\nSettings:");
		if (somaticProcessing) {
			IO.pl("\tSomatic processing, unset fields are set to Somatic defaults and the following logic is used to pass or fail each record.\n"
					+ "\t\tIf ID matches pass irregardless of other settings.\n"
					+ "\t\tFail it if DP or AF are false.\n"
					+ "\t\tFail it if ClinSig, Splice, and Impact are all false.\n"
					+ "\t\tFail it if Pop or BKAF are present and they are false.\n");
		}
		if (minimumDP != 0) IO.pl("\t"+ minimumDP+"\t: Minimum DP read depth");
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) {
			IO.pl("\t"+ minimumAF+"\t: Minimum AF allele frequency");
			IO.pl("\t"+ maximumAF+"\t: Maximum AF allele frequency");
		}
		if (maximumPopAF != 0) IO.pl("\t"+ maximumPopAF+"\t: Maximum population AF from dbNSFP_ExAC_AF or dbNSFP_1000Gp3_AF");
		if (maxFracBKAFs != 0) IO.pl("\t"+ maxFracBKAFs+"\t: Maximum fraction of background samples with an AF >= observed AF");
		if (minimumBKZ != 0) IO.pl("\t"+ minimumBKZ+"\t: Minimum BKZ quality score");
		if (passingAnnImpact != null) IO.pl("\n\t"+Misc.treeSetToString(passingAnnImpact, ",")+"\t: ANN impact keys");
		if (passingClinSig != null) IO.pl("\n\t"+Misc.treeSetToString(passingClinSig, ",")+"\t: CLINSIG keep keys");
		if (excludeClinSig != null) IO.pl("\t"+Misc.treeSetToString(excludeClinSig, ",")+"\t: CLINSIG exclude keys");
		if (passingVCFSS != null) {
			IO.pl("\n\t"+Misc.stringArrayToString(passingVCFSS, ",")+"\t: Splice junction types to scan");
			IO.pl("\t"+ minimumVCFSSDiff+"\t: Minimum difference in MaxEnt scan scores for a splice junction effect");
		}
		if (somaticProcessing ==false) IO.pl("\t"+ orAnnos+"\t: Require that only one need pass: ANN Impact or Clinvar or Splice Effect");
		if (passingIDKeys != null) IO.pl("\n\t"+Misc.stringArrayToString(passingIDKeys, ",")+"\t: VCF ID keys that if present pass it regardless of any filter settings");
		if (passingGermlineGenes != null) IO.pl("\n\t"+acmgGenes+"\t: ACMG Genes that cause the maxAF filter to be skipped for intersecting variants");
		if (somaticProcessing ==false) 
		IO.pl("\n\t"+verbose+"\t: Verbose output");
		IO.pl("\t"+ saveDirectory+"\t: Save directory");
	}
	
	public String format(double numerator, double denominator){
		String frac =  Num.formatNumber(100*numerator/denominator, 3);
		return frac +"%\t" +(int)numerator+"/"+(int)denominator+"\t";
	}
	
	public void printStats(){
		
		IO.pl("\nProcessing statistics for all files:");
		IO.pl("\t"+ numRecords+"\tNumber of VCF Records");
		if (minimumDP != 0) IO.pl("\t"+ format(numPassingDP, numRecords)+"Passing DP");
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) {
			IO.pl("\t"+ format(numPassingMinAF, numRecords)+"Passing MinAF");
			IO.pl("\t"+ format(numPassingMaxAF, numRecords)+"Passing MaxAF");
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
			IO.pl("\t"+ format(numPassingClinSing, numWithClin)+"Passing ClinSig");
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
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		String impactString = null;
		String clinString = null;
		String clinStringExclude = null;
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'd': minimumDP = Integer.parseInt(args[++i]); break;
					case 'm': minimumAF = Double.parseDouble(args[++i]); break;
					case 'x': maximumAF = Double.parseDouble(args[++i]); break;
					case 'p': maximumPopAF = Double.parseDouble(args[++i]); break;
					case 'b': maxFracBKAFs = Double.parseDouble(args[++i]); break;
					case 'g': passingVCFSS = Misc.COMMA.split(args[++i]); break;
					case 'n': minimumVCFSSDiff = Double.parseDouble(args[++i]); break;
					case 'z': minimumBKZ = Double.parseDouble(args[++i]); break;
					case 'a': impactString = args[++i]; break;
					case 'c': clinString = args[++i]; break;
					case 'e': clinStringExclude = args[++i]; break;
					case 'j': createGermlineGenes(); break;
					case 'i': passingIDKeys = Misc.COMMA.split(args[++i]); break;
					case 'o': orAnnos = true; break;
					case 'r': verbose = true; break;
					case 'f': somaticProcessing = true; break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'y': new File(args[++i]); break;
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
		if (minimumVCFSSDiff !=0 && passingVCFSS == null) Misc.printErrAndExit("\nError: please provide a comma delimited list of splice junction types (e.g. D5S,D3S,G5S,G3S) to examine.\n");
		
	}
	

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                             Annotated Vcf Parser  Sept 2018                      **\n" +
				"**************************************************************************************\n" +
				"Splits VCF files that have been annotated with SnpEff w/ dbNSFP and clinvar, plus the\n"+
				"VCFBackgroundChecker and VCFSpliceScanner USeq apps into passing and failing records.\n"+
				"Use the -e option to inspect the effect of the various filters on each record. Use\n"+
				"the VCFRegionFilter app to restrict variants to particular gene regions.\n"+

				"\nOptions:\n"+
				"-v File path or directory containing xxx.vcf(.gz/.zip OK) file(s) to filter.\n" +
				"-s Directory for saving the results.\n"+
				"-f Perform a candidate somatic variant processing. Setting the following overrides\n"+
				"        the defaults. \n"+
				"-d Minimum DP alignment depth\n"+
				"-m Minimum AF allele frequency\n"+
				"-x Maximum AF allele frequency\n"+
				"-j Ignore the max AF filter for ACMG incidental germline gene variants.\n"+
                "-p Maximum population allele frequency, only applies if present.\n"+
                "-b Maximum fraction of BKAF samples with allele frequency >= VCF AF, only applies\n"+
                "       if present.\n"+
                "-g Splice junction types to scan.\n"+
                "-n Minimum difference in splice junction scores, only applies if present.\n"+
                "-a Comma delimited list of SnpEff ANN impact categories to select for.\n"+
                "-c Comma delimited list of CLINSIG terms to select for.\n"+
                "-e Comma delimited list of CLINSIG terms to select against.\n"+
                "-i Comma delimited list of VCF ID keys to select for. If the VCF ID contains one or\n"+
                "       more, the record is passed regardless of other filters. The match is not exact.\n"+
                "-o Only require, if set or present, SnpEff ANN or CLINSIG or Splice to be true to pass.\n"+
                "       Defaults to require that all set pass.\n"+
				"-r Verbose per record output.\n"+ 
                "-y Path to a config txt file for setting the above.\n"+

				"\nExample: java -jar pathToUSeq/Apps/AnnotatedVcfParser -v VCFFiles/ -s Parsed/\n"+
				"        -d 75 -m 0.05 -x 0.75 -j -p 0.02 -b 0.1 -g D5S,D3S,G5S,G3S -n 3.5 -a\n"+
				"        HIGH,MODERATE -c Pathogenic,Likely_pathogenic -i Foundation,Tempus -v\n\n"+


				"**************************************************************************************\n");
	}
}
