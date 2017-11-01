package edu.utah.seq.vcf; 

import java.io.BufferedReader;
import java.io.File;
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
	private String[] passingIDKeys = null;
	private TreeSet<String> passingAnnImpact = null;
	private TreeSet<String> passingClinSig = null;
	private String[] passingVCFSS = null;
	private double minimumVCFSSDiff = 0; //difference between ref and alt for splice junction score
	private boolean verbose = false;
	private boolean foundationProcessing = false;
	
	//Patterns
	private static final Pattern AF = Pattern.compile(";*AF=([\\d\\.]+);");
	private static final Pattern DP = Pattern.compile(";*DP=(\\d+);");
	private static final Pattern EXAC = Pattern.compile(";*dbNSFP_ExAC_AF=([\\d\\.]+);");
	private static final Pattern G1000 = Pattern.compile(";*dbNSFP_1000Gp3_AF=([\\d\\.]+);");
	private static final Pattern CLINSIG = Pattern.compile(";*CLNSIG=(\\w+);");
	private static final Pattern BKAF = Pattern.compile(";*BKAF=([\\d\\.,]+);");
	
	//trackers
	private Histogram afs = new Histogram(0, 1.01, 101);
	private TreeMap<String, Integer> clinsig = new TreeMap<String, Integer>();
	private int[] dps = new int[10000];
	private TreeMap<String, Integer> impacts = new TreeMap<String, Integer>();
	
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
		if (foundationProcessing) modifySettingsForFoundation();
		IO.p("Parsing:");
		for (File vcf: vcfFiles){
			IO.p("\t"+vcf.getName());
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
				if (verbose) IO.p("\n"+vcfLine);

				//#CHROM POS ID REF ALT QUAL FILTER INFO ......
				//   0    1   2  3   4   5     6      7
				String[] cells = Misc.TAB.split(vcfLine);
				String[] info = Misc.SEMI_COLON.split(cells[7]);
				
				boolean passDP=true, passAF=true, passImpact=true, passClinSig=true, passSplice=true, passID = true;
				boolean[] passBKAF= null; 
				boolean[] passPop= null;
				
				if (minimumDP != 0) 			passDP = checkDP(cells[7]);
				if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) passAF = checkAF(cells[7]);
				if (maximumPopAF != 0) 			passPop = checkPopFreq(cells[7]);
				if (passingAnnImpact != null) 	passImpact = checkImpact(info);
				if (passingClinSig != null) 	passClinSig = checkClinSig(cells[7]);
				if (maxFracBKAFs != 0) 			passBKAF = checkBKAFs(cells[7]);
				if (passingVCFSS != null) 		passSplice = checkSplice(cells[7], info);
				if (passingIDKeys != null) 		passID = checkIDs(cells[2]);
				
				boolean pass;
				if (foundationProcessing) pass =  passFoundationFiltering(passID, passDP, passAF, passImpact, passClinSig, passSplice, passPop, passBKAF); 
				else pass = allPass(passID, passDP, passAF, passImpact, passClinSig, passSplice, passPop, passBKAF);
				
				if (pass){
					numPass++;
					passVcf.println(vcfLine);
				}
				else failVcf.println(vcfLine);
				if (verbose) IO.p("\tPass all filters\t"+pass);
					
			}
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
	private boolean allPass(boolean passID, boolean passDP, boolean passAF, boolean passImpact, boolean passClinSig, boolean passSplice, boolean[] passPop, boolean[] passBKAF) {
		if (passID) return true;
		boolean passBKAF2 = true;
		if (passBKAF != null && passBKAF[0] == true) passBKAF2 = passBKAF[1];
		boolean passPop2 = true;
		if (passPop != null && passPop[0] == true) passPop2 = passPop[1];
		return (passDP && passAF && passImpact && passClinSig && passSplice && passID && passBKAF2 && passPop2);
	}


	private boolean passFoundationFiltering(boolean passID, boolean passDP, boolean passAF,  boolean passImpact, boolean passClinSig, boolean passSplice, boolean[] passPop, boolean[] passBKAF){
		//pass it if ID correct
		if (passID) return true;
		//fail it if DP or AF fails
		if (passAF == false || passDP == false) return false;
		//fail it if all are false
		if (passClinSig == false && passSplice == false && passImpact == false) return false;
		//if present and fails then fail
		if (passPop[0] == true && passPop[1] == false) return false;
		if (passBKAF[0] == true && passBKAF[1] == false) return false;
		return true;
	}
	
	private boolean checkIDs(String id) {
		
		//for each user provided key (all lower cased), check if in id and pass if so
		String lcId = id.toLowerCase();
		for (String key: passingIDKeys){
			if (lcId.contains(key)){
				numPassIDs++;
				if (verbose) IO.p("\tID Check\ttrue"+key+"\t"+lcId);
				return true;
			}
		}
		if (verbose) IO.p("\tID Check\tfalse\t"+lcId);
		return false;
	}

	private boolean checkSplice(String info, String[] splitInfo) throws Exception {
		//any splice info
		if (info.contains("VCFSS=") == false){
			if (verbose) IO.p("\tSplice Check\tfalse\tNo VCFSS");
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
		if (verbose) IO.p("\tSplice Check\t"+passSplice+"\t"+maxDelta+"\t"+vcfSSLine);	
		if (passSplice) {
			numPassingSplice++;
			return true;
		}
		else return false;
	}

	
	private boolean[] checkBKAFs(String info) throws Exception {
		Matcher mat = BKAF.matcher(info);
		if (mat.find()){
			
			//record findings
			numWithBKAFs++;
			String bkafs = mat.group(1);
			double[] afs = Num.stringArrayToDouble(bkafs, ",");
			
			//score
			double numFailing = 0;
			for (double af: afs) if (af>= obsAF) numFailing++;
			double fractFailing = numFailing/(double)afs.length;
			boolean passingBkaf = fractFailing <= maxFracBKAFs; 
			
			if (verbose) IO.p("\tBKAF Check\t"+passingBkaf+"\t"+fractFailing+"\t"+bkafs);
			if (passingBkaf) {
				numPassingBKAFs++;
				return new boolean[]{true, true};
			}
			else return new boolean[]{true, false};
		}
		else {
			if (verbose) IO.p("\tBKAF Check\tfalse\tNo BKAFs");
			return new boolean[]{false, false};
		}
	}

	
	private boolean checkClinSig(String info) throws Exception {
		Matcher mat = CLINSIG.matcher(info);
		if (mat.find()){
			//record findings
			numWithClin++;
			String cs = mat.group(1);
			if (clinsig.containsKey(cs)) {
				int count = clinsig.get(cs);
				clinsig.put(cs, new Integer( ++count ));
			}
			else clinsig.put(cs, new Integer(1));

			boolean passClinSig = passingClinSig.contains(cs.toLowerCase());
			if (verbose) IO.p("\tCLINSIG Check\t"+passClinSig+"\t"+cs);
			if (passClinSig) {
				numPassingClinSing++;
				return true;
			}
			else return false;

		}
		else {
			if (verbose) IO.p("\tCLINSIG Check\tfalse\tNo CLINSIG");
			return false;
		}
	}
	
	/**
	 * ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
	 * 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
	 *    0          1                2               3          4            5
	 * ANN=GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20)

	 * */
	private boolean checkImpact(String[] splitInfo) throws Exception {
		//find ANN=
		String ann = null;
		for (String i: splitInfo){
			if (i.startsWith("ANN=")) {
				ann = i;
				break;
			}
		}
		if (ann == null) {
			if (verbose) IO.p("\tImpact Check\tfalse\tNo ANN");
			return false;
		}
		numWithAnn++;
		//split ANN and save impact
		String[] splitAnn = Misc.PIPE.split(ann);
		String impact = splitAnn[2];
		if (impacts.containsKey(impact)) impacts.put(impact, new Integer( impacts.get(impact)+1 ));
		else impacts.put(impact, new Integer(1));
		
		//check if impact is one of the ones they want
		if (passingAnnImpact!=null){
			boolean passImpact = passingAnnImpact.contains(impact.toLowerCase());
			if (verbose) IO.p("\tImpact Check\t"+passImpact+"\t"+impact);
			if (passImpact) {
				numPassingAnnImpact++;
				return true;
			}
			else return false;
		}
		else return false;
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
				if (verbose) IO.p("\tExAC Check\t"+passExAC+"\t"+popAF);
			}
		}

		//1K genomes
		mat = G1000.matcher(info);
		if (mat.find()){
			if (mat.group(1).equals(".") == false){
				foundPop = true;
				double popAF = Double.parseDouble(mat.group(1));
				if (popAF > maximumPopAF) pass1K = false;
				if (verbose) IO.p("\t1KG Check\t"+pass1K+"\t"+popAF);
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
			if (verbose) IO.p("\tPopFreq Check\tfalse\tNo dbNSFP_1000Gp3_AF or dbNSFP_ExAC_AF");
			return new boolean[]{false, false};
		}
		
	}


	private boolean checkDP(String info) throws Exception {
		Matcher mat = DP.matcher(info);
		if (mat.find()){
			int obsDP = Integer.parseInt(mat.group(1));
			if (obsDP < dps.length) dps[obsDP]++;
			boolean passDP = obsDP >= minimumDP;
			if (verbose) IO.p("\tDP Check\t"+passDP+"\t"+obsDP);
			if (passDP == false) return false;
			else {
				numPassingDP++;
				return true;
			}
		}
		else {
			if (verbose) IO.p("\tDP Check\tfalse\tNo DP");
			return false;
		}
	}
	
	private boolean checkAF(String info) throws Exception {
		Matcher mat = AF.matcher(info);
		if (mat.find()){
			obsAF = Double.parseDouble(mat.group(1));
			afs.count(obsAF);
			boolean passMinAF = obsAF >= minimumAF;
			boolean passMaxAF = obsAF <= maximumAF;
			if (verbose) {
				IO.p("\tAF Min Check\t"+passMinAF+"\t"+obsAF);
				IO.p("\tAF Max Check\t"+passMaxAF+"\t"+obsAF);
			}
			if (passMinAF) numPassingMinAF++;
			if (passMaxAF) numPassingMaxAF++;
			if (passMinAF && passMaxAF) return true;
			else return false;
		}
		else {
			if (verbose) IO.p("\tAF Check\tfalse\tNo AF");
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
		if (maximumAF == 1) maximumAF = 0.5;
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
		
	}

	public void printSettings(){
		IO.p("\nSettings:");
		if (foundationProcessing) {
			IO.p("\tFoundation processing, unset fields are set to Foundation defaults and the following logic is used to pass or fail each record.\n"
					+ "\t\tIf ID matches pass irregardless of other settings.\n"
					+ "\t\tFail it if DP or AF are false.\n"
					+ "\t\tFail it if ClinSig, Splice, and Impact are all false.\n"
					+ "\t\tFail it if Pop or BKAF are present and they are false.\n");
		}
		if (minimumDP != 0) IO.p("\t"+ minimumDP+"\t: Minimum DP read depth");
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) {
			IO.p("\t"+ minimumAF+"\t: Minimum AF allele frequency");
			IO.p("\t"+ maximumAF+"\t: Maximum AF allele frequency");
		}
		if (maximumPopAF != 0) IO.p("\t"+ maximumPopAF+"\t: Maximum population AF from dbNSFP_ExAC_AF or dbNSFP_1000Gp3_AF");
		if (maxFracBKAFs != 0) IO.p("\t"+ maxFracBKAFs+"\t: Maximum fraction of background samples with an AF >= observed AF");
		if (passingAnnImpact != null) IO.p("\n\t"+Misc.treeSetToString(passingAnnImpact, ",")+"\t: ANN impact keys");
		if (passingClinSig != null) IO.p("\n\t"+Misc.treeSetToString(passingClinSig, ",")+"\t: CLINSIG keys");
		if (passingVCFSS != null) {
			IO.p("\n\t"+Misc.stringArrayToString(passingVCFSS, ",")+"\t: Splice junction types to scan");
			IO.p("\t"+ minimumVCFSSDiff+"\t: Minimum difference in MaxEnt scan scores for a splice junction effect");
		}
		if (passingIDKeys != null) IO.p("\n\t"+Misc.stringArrayToString(passingIDKeys, ",")+"\t: VCF ID keys that if present pass it regardless of any filter settings");
		IO.p("\n\t"+verbose+"\t: Verbose output");
		IO.p("\t"+ saveDirectory+"\t: Save directory");
	}
	
	public String format(double numerator, double denominator){
		String frac =  Num.formatNumber(100*numerator/denominator, 3);
		return frac +"%\t" +(int)numerator+"/"+(int)denominator+"\t";
	}
	
	public void printStats(){
		
		IO.p("\nProcessing statistics for all files:");
		IO.p("\t"+ numRecords+"\tNumber of VCF Records");
		if (minimumDP != 0) IO.p("\t"+ format(numPassingDP, numRecords)+"Passing DP");
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0) {
			IO.p("\t"+ format(numPassingMinAF, numRecords)+"Passing MinAF");
			IO.p("\t"+ format(numPassingMaxAF, numRecords)+"Passing MaxAF");
		}
		if (maximumPopAF != 0){
			IO.p("\t"+ format(numWithPopAF, numRecords)+"With Pop AF");
			IO.p("\t"+ format(numPassingPopAF, numWithPopAF)+"Passing Pop AF");
		}
		if (passingAnnImpact != null) {
			IO.p("\t"+ format(numWithAnn, numRecords)+"With Ann");
			IO.p("\t"+ format(numPassingAnnImpact, numWithAnn)+"Passing Ann Impact");
		}
		if (maxFracBKAFs != 0) {
			IO.p("\t"+ format(numWithBKAFs, numRecords)+"With BKAFs");
			IO.p("\t"+ format(numPassingBKAFs, numWithBKAFs)+"Passing BKAFs");
		}
		if (passingClinSig != null){
			IO.p("\t"+ format(numWithClin, numRecords)+"With Clin");
			IO.p("\t"+ format(numPassingClinSing, numWithClin)+"Passing ClinSig");
		}
		if (passingVCFSS != null){
			IO.p("\t"+ format(numWithSplice, numRecords)+"With Splice");
			IO.p("\t"+ format(numPassingSplice, numWithSplice)+"Passing Splice");
		}
		if (passingIDKeys != null) IO.p("\t"+ format(numPassIDs, numRecords)+"Pass IDs");
		IO.p("\t"+ format(numPass, numRecords)+"Records passing all filters");
	}
	
	public void printDistributions(){
		if (passingClinSig != null){
			IO.p("\nObserved CLINSIG annotations and their frequency:" );
			for (String key: clinsig.keySet()) IO.p("\t"+key+"\t"+clinsig.get(key));
		}
		if (passingAnnImpact != null){
			IO.p("\nObserved ANN impacts and their frequency:");
			for (String key: impacts.keySet()) IO.p("\t"+key+"\t"+impacts.get(key));
		}
		if (minimumAF != 0 || maximumAF != 1 || maxFracBKAFs != 0){
			IO.p("\nAF Allele Frequency distribution:");
			afs.printScaledHistogram();
		}
		if (minimumDP != 0) {
			IO.p("\nDP Read Depth distribution:");
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


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		String impactString = null;
		String clinString = null;
		
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
					case 'a': impactString = args[++i]; break;
					case 'c': clinString = args[++i]; break;
					case 'i': passingIDKeys = Misc.COMMA.split(args[++i]); break;
					case 'e': verbose = true; break;
					case 'f': foundationProcessing = true; break;
					case 's': saveDirectory = new File(args[++i]); break;
					
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		IO.p("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

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
		
		if (minimumVCFSSDiff !=0 && passingVCFSS == null) Misc.printErrAndExit("\nError: please provide a comma delimited list of splice junction types (e.g. D5S,D3S,G5S,G3S) to examine.\n");
		
	}
	

	public static void printDocs(){
		IO.p("\n" +
				"**************************************************************************************\n" +
				"**                              Annotated Vcf Parser  Oct 2017                      **\n" +
				"**************************************************************************************\n" +
				"Splits VCF files that have been annotated with SnpEff, dbNSFP, clinvar, plus the USeq\n"+
				"VCFBackgroundChecker and VCFSpliceScanner apps into passing and failing records.\n"+
				"Unless the -f flag is selected, records must pass all of the filters set below.\n"+
				"Use the -e option to inspect the effect of the various filters on each record.\n"+

				"\nOptions:\n"+
				"-v File path or directory containing xxx.vcf(.gz/.zip OK) file(s) to filter.\n" +
				"-s Directory for saving the results.\n"+
				"-f Perform a Foundation dataset processing. Setting the following overrides the\n"+
				"        defaults. \n"+
				"-d Minimum AD alignment depth\n"+
				"-m Minimum AF allele frequency\n"+
				"-x Maximum AF allele frequency\n"+
                "-p Maximum population allele frequency, only applies if present.\n"+
                "-b Maximum fraction of BKAF samples with allele frequency >= VCF AF, only applies\n"+
                "       if present.\n"+
                "-g Splice junction types to scan.\n"+
                "-n Minimum difference in splice junction scores, only applies if present.\n"+
                "-a Comma delimited list of SnpEff ANN impact categories to select for.\n"+
                "-c Comma delimited list of CLINSIG terms to select for.\n"+
                "-i Comma delimited list of VCF ID keys to select for. If the VCF ID contains one or\n"+
                "       more, the record is passed regardless of other filters. The match is not exact.\n "+
				"-v Verbose per record output.\n"+ 

				"\nExample: java -jar pathToUSeq/Apps/AnnotatedVcfParser -v VCFFiles/ -s Parsed/\n"+
				"        -d 75 -m 0.05 -x 0.75 -p 0.02 -b 0.1 -g D5S,D3S,G5S,G3S -n 3.5 -a\n"+
				"        HIGH,MODERATE -c Pathogenic,Likely_pathogenic -i Foundation,Tempus -v\n\n"+


				"**************************************************************************************\n");
	}
	/*case 'v': forExtraction = new File(args[++i]); break;
					case 'a': impactString = args[++i]; break;
					case 'c': clinString = args[++i]; break;
					case 'i': passingIDKeys = Misc.COMMA.split(args[++i]); break;
					case 'e': verbose = true; break;
					*/

}
