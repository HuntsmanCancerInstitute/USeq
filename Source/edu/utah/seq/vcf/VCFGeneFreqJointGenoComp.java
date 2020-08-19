package edu.utah.seq.vcf;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;

/**For SnpEff annotated joint genotyped/ multiple sample vcf file, counts the number of affected genes and compares their counts to those in another group using a fisher's exact test.
 * @author Nix
 * */
public class VCFGeneFreqJointGenoComp {

	//user fields
	private File vcfFile;
	private double minReadDepth = 20;
	private double genotypeQuality = 20;
	private double qual = 20;
	private Pattern groupAPattern = null;
	private Pattern groupBPattern = null;
	
	//VCFSpliceScanner parsing
	private String[] passingVCFSS = new String[]{"D5S", "D3S"};	
	private double minimumVCFSSDiff = 5;
	
	
	private TreeMap<String, ArrayList<JGSample>> geneSamples = new  TreeMap<String, ArrayList<JGSample>>();
	private String[] sampleNames = null;
	private Pattern geneInfo = Pattern.compile(".+;GENEINFO=([a-zA-Z_0-9\\-\\.]+):.+");
	private int maxGeneCount = 0;
	private boolean debug = false;
	private int numRecords = 0;
	private int numFailingQual = 0;
	private HashMap<String, String> formatValues = new HashMap<String,String>();
	private int[] numPass = null;
	private int[] numFail = null;
	private TreeSet<String> groupASampleNames = new TreeSet<String>();
	private TreeSet<String> groupBSampleNames = new TreeSet<String>();
	
	//constructor
	public VCFGeneFreqJointGenoComp(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//parse vcf lines associated with a gene name
		IO.pl("\nParsing Vcf Records...");
		parseVcf();
		
		IO.pl("\nComparing genes...");
		compareGroups();
		
		IO.pl("\nGroup A Sample Names: "+groupASampleNames);
		IO.pl("Group B Sample Names: "+groupBSampleNames);
		if (groupASampleNames.size() == 0 || groupBSampleNames.size() == 0) Misc.printErrAndExit("\nERROR! One or both of the groups had zero associated samples?!");

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public void parseVcf() {
		String line = null;
		boolean failed = false;
		BufferedReader in = null;
		try {
			//load the header and get the sample count
			in = loadHeader();
			
			numPass = new int[sampleNames.length];
			numFail = new int[sampleNames.length];

			//for each line in the file
			while ((line = in.readLine()) != null){
				if (debug) System.out.println("\n"+line);
				numRecords++;
				
				//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1, Sample2.....
				//   0    1   2  3   4   5     6      7     8      9       10
				String[] fields = Misc.TAB.split(line);
				
				if (fields[4].equals("*")) {
					numFailingQual++;
					continue;
				}

				//check whole line QUAL
				double q = Double.parseDouble(fields[5]);
				if (q < qual) {
					numFailingQual++;
					continue;
				}
				
				//pull the clinvar gene for patho
				String gene = fetchClinvarPathoGene(fields[7]);
//if (gene != null) IO.pl("\tCLINVAR "+gene);
				
				//pull the VCFSS gene name
				if (gene == null) gene = fetchVcfSpliceJunctionEffect(fields[7]);
//if (gene != null) IO.pl("\tSPLICE "+gene);
				
				//pull first Anno 
				if (gene == null) {
					//not found so find first ANN
					String[] splitInfo = Misc.SEMI_COLON.split(fields[7]);
					gene = fetchFirstAffectedGene (splitInfo);
//if (gene != null) IO.pl("\tANNO "+gene);
				}
				
				if (gene == null) return;
				
				ArrayList<JGSample> samples = geneSamples.get(gene);
				if (samples == null) {
					samples = new ArrayList<JGSample>();
					geneSamples.put(gene, samples);
				}

				String[] format = Misc.COLON.split(fields[8]);

				//for each sample
				int index = 0;
				for (int i=9; i< fields.length; i++){
					double[] dpAfGt = checkSample(format, fields[i]);
					if (dpAfGt == null) numFail[index]++;
					else {
						numPass[index]++;
						samples.add(new JGSample(sampleNames[index], dpAfGt));
					}
					index++;
				}
			}

			//print stats
			printStats();
			
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("ERROR: parsing line "+line);
			failed = true;
		} finally {
			//close io
			try {
				if (in != null) in.close();
			} catch (IOException e) {}
			if (failed) System.exit(1);
		}
	}
	
	private class JGSample {
		String sampleName;
		double[] dpAfGt;
		
		public JGSample(String sampleName, double[] dpAfGt) {
			this.sampleName = sampleName;
			this.dpAfGt = dpAfGt;
		}

	}
	
	/**Returns null if fails, otherwise the DP, AF, and Gt (0, 0.5, 1).*/
	private double[] checkSample(String[] ids, String sample) throws Exception {
		if (debug) System.out.println("\t"+sample);
		//no values?
		if (sample.startsWith("./.")) {
			if (debug) System.out.println ("\t\tEmpty");
			return null;
		}
		
		formatValues.clear();
		String[] values = Misc.COLON.split(sample);
		if (ids.length != values.length) throw new Exception ("\t\tThe FORMAT length and sample values do not match");
		
		for (int i=0; i< ids.length; i++) formatValues.put(ids[i], values[i]);
		
		//check GQ?
		if (genotypeQuality != 0){
			String gqString = formatValues.get("GQ");
			if (gqString == null || gqString.equals(".")) {
				if (debug) System.out.println ("\t\tNo GQ");
				return null;
			}
			else {
				double gq = Double.parseDouble(gqString);
				if (gq < genotypeQuality) {
					if (debug) System.out.println ("\t\tFailing GQ");
					return null;
				}
			}
		}
		
		//parse genotype
		String gtString = formatValues.get("GT");
		double gt = -1;
		if (gtString == null ) {
			if (debug) System.out.println ("\t\tNo GT");
			return null;
		}
		else {
			if (gtString.equals("1/1")) gt = 1;
			else if (gtString.equals("0/0")) gt = 0;
			else if (gtString.equals("0/1") || gtString.equals("1/0")) gt = 0.5;
			else return null;
		}
		
		
		//check depth, should be only one alt!
		String ad = formatValues.get("AD");
		if (ad == null) throw new Exception ("\t\tNo AD?");

		String[] refAltCounts = Misc.COMMA.split(ad);
		if (refAltCounts.length != 2) throw new Exception ("\nThe AD for this sample does not contain 2 values, did you vt decompose? "+ad+" "+sample);
		double refCounts = Double.parseDouble(refAltCounts[0]);
		double altCounts = Double.parseDouble(refAltCounts[1]);
		double dp = refCounts + altCounts;
		double af = altCounts/dp;
		if (dp >= minReadDepth) return new double[]{dp, af, gt};
		
		if (debug) System.out.println ("\t\tFailing dp ("+dp+")");
		return null;
	}
	
	private void printStats() {
		System.out.println("\nParsing Statistics:");
		System.out.println(numRecords+ "\tNumber records");
		System.out.println(numFailingQual+ "\tNumber failing QUAL");
		System.out.println("\nPassing\tFailing\tSampleName");
		for (int i=0; i< sampleNames.length; i++){
			System.out.println(numPass[i] + "\t"+ numFail[i] + "\t"+ sampleNames[i] );
		}
		
	}
	
	public BufferedReader loadHeader() throws Exception{
		BufferedReader in = IO.fetchBufferedReader(vcfFile);
		String line = null;
		//for each line in the file
		while ((line = in.readLine()) != null){
			line = line.trim();
			//header? just print out
			if (line.startsWith("#")) {
				if (line.startsWith("#CHROM")){
					//save sample names
					String[] columns = Misc.TAB.split(line);
					int numSamples = columns.length- 9;
					sampleNames = new String[numSamples];
					int index = 0;
					for (int i= 9; i< columns.length; i++) sampleNames[index++] = columns[i];
					return in;
				}
			}
		}
		throw new Exception("\nERROR: failed to find the #CHROM line?");
	}

	
	private void compareGroups() {
		try {
			int maxCellCount = 1000;
			FisherExact fe = new FisherExact(maxCellCount);
			double numTests = geneSamples.size();

			//for each gene
			IO.pl("\nGene\t#AMut\t#AWt\t#BMut\t#BWt\tAMutFreq\tBMutFreq\tUnCorrPVal\tCorrPVal");
			for (String gene: geneSamples.keySet()) {
				ArrayList<JGSample> samples = geneSamples.get(gene);
				TreeMap<String, JGSample> sampleName = splitSamplesByName(samples);
				int[] numMutWtAB = splitByPattern(sampleName);
				
				//any mut var samples come through?
				if (numMutWtAB[0] == 0 && numMutWtAB[2] == 0 ) continue;
				
				int total = Num.sumIntArray(numMutWtAB);
				if (total > maxCellCount) {
					maxCellCount = total;
					fe = new FisherExact(maxCellCount);
				}

				double freqMutA = (double)numMutWtAB[0]/(double)(numMutWtAB[0]+numMutWtAB[1]);
				double freqMutB = (double)numMutWtAB[2]/(double)(numMutWtAB[2]+numMutWtAB[3]);
				double pval = fe.getTwoTailedP(numMutWtAB[0], numMutWtAB[1], numMutWtAB[2], numMutWtAB[3]);
				double corrPval = pval * numTests;
				if (pval > 0.999 ) {
					pval = 1;
					corrPval =1;
				}
				else if (corrPval > 1) corrPval = 1;

				IO.pl(gene+"\t"+ numMutWtAB[0]+"\t"+ numMutWtAB[1]+"\t"+ numMutWtAB[2]+"\t"+ numMutWtAB[3]+"\t"+ Num.formatNumber(freqMutA, 4)+"\t"+ Num.formatNumber(freqMutB, 4)+"\t"+ pval+"\t"+ corrPval);
			}
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem comparing groups!");
		}
	}

	private int[] splitByPattern(TreeMap<String, JGSample> samples) throws IOException {
		int numAWt = 0;
		int numAMut = 0;
		int numBWt = 0;
		int numBMut = 0;
		
		for (String sampleName: samples.keySet()) {
			JGSample jgs = samples.get(sampleName);
			
			Matcher matA = groupAPattern.matcher(sampleName);
			if (matA.matches()) {
				groupASampleNames.add(sampleName);
				if (jgs.dpAfGt[2]<0.5) numAWt++;
				else numAMut++;
			}
			else {
				Matcher matB = groupBPattern.matcher(sampleName);
				if (matB.matches()) {
					groupBSampleNames.add(sampleName);
					if (jgs.dpAfGt[2]<0.5) numBWt++;
					else numBMut++;
				}
				else throw new IOException ("Failed to match either of the group patterns for sample name "+sampleName);
			}
		}
		return new int[] {numAMut, numAWt, numBMut, numBWt};
	}

	private TreeMap<String, JGSample> splitSamplesByName(ArrayList<JGSample> samples) {
		TreeMap<String, JGSample> js = new TreeMap<String, JGSample>();
		for (JGSample s: samples) {
			JGSample prior = js.get(s.sampleName);
			//add if null or is a het or hom
			if (prior == null ) {
				prior = s;
				js.put(s.sampleName, prior);
			}
			else if (prior.dpAfGt[2] < s.dpAfGt[2]) js.put(s.sampleName, s);
		}
		return js;
	}

	private int countVars(ParsedVcfStats[] pvss) {
		int num = 0;
		for (ParsedVcfStats p: pvss) num+=p.numVars;
		return num;
	}

	private HashMap<String, Integer> mergeCountsPrintStats(ParsedVcfStats[] groupA) {
		HashMap<String, Integer> geneNameCounts = new HashMap<String, Integer>();
		double totalNumVars = 0;
		for (int i=0; i<groupA.length; i++) {
			ParsedVcfStats pvs = groupA[i];
			totalNumVars+= (double) pvs.numVars;
			IO.pl(Misc.removeExtension(pvs.vcfFile.getName())+"\t"+pvs.numVars+"\t"+pvs.geneNameCounts.size());
			for (String gn: pvs.geneNameCounts.keySet()) {
				int numHits = pvs.geneNameCounts.get(gn);
				Integer numTotalHits = geneNameCounts.get(gn);
				if (numTotalHits != null) numHits += numTotalHits;
				geneNameCounts.put(gn, new Integer(numHits));
				if (numHits > maxGeneCount) maxGeneCount = numHits;
			}
		}
		double aveVarsPerFile = totalNumVars/(double)groupA.length;
		IO.pl("\tMean # vars\t"+ Num.formatNumber(aveVarsPerFile, 2));
		return geneNameCounts;
	}

	private class ParsedVcfStats {
		File vcfFile = null;
		int numVars = 0;
		HashMap<String, Integer> geneNameCounts = null;
		
		ParsedVcfStats (File vcfFile, int numVars, HashMap<String, Integer> geneNameCounts) {
			this.vcfFile = vcfFile;
			this.numVars = numVars;
			this.geneNameCounts = geneNameCounts;
		}
	}
	
	private ParsedVcfStats[] loadVcfs(File[] vcfs) {
		ParsedVcfStats[] pvs = new ParsedVcfStats[vcfs.length];
		for (int i=0; i< vcfs.length; i++) pvs[i] = countVcfFile(vcfs[i]);
		return pvs;
	}

	private ParsedVcfStats countVcfFile (File vcf) {
		BufferedReader in = null;
		String line = null;
		HashMap<String, Integer> geneNameCounts = new HashMap<String, Integer>();
		int numVars = 0;
		try {
			in = IO.fetchBufferedReader(vcf);
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
					numVars++;
					String[] fields = Misc.TAB.split(line);
					
					//attempt to pull the clinvar gene for patho
					String gene = fetchClinvarPathoGene(fields[7]);
					if (gene == null) {
						//not found so find most freq affected
						String[] splitInfo = Misc.SEMI_COLON.split(fields[7]);
						LinkedHashMap<String, Integer> geneNameCount = fetchAffectedGenes (splitInfo);
						gene = findFirstLargest( geneNameCount);
					}
					Integer i = geneNameCounts.get(gene);
					if (i == null) i = new Integer(1);
					else i = new Integer (i.intValue()+1);
					geneNameCounts.put(gene,i);
//if (line.contains("|NF1|")) IO.pl(gene+" --- "+fields[7]);
				}
			}
			return new ParsedVcfStats(vcf, numVars, geneNameCounts);

		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing gene name from "+vcf+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
		return null;
	}
	
	private String fetchVcfSpliceJunctionEffect(String info) {
		if (info.contains("VCFSS=") == false) return null;
		
		String[] f = Misc.SEMI_COLON.split(info);
		String vcfss = null;
		for (String s: f) {
			if (s.startsWith("VCFSS=")) {
				vcfss= s.substring(6);
				break;
			}
		}
		if (vcfss == null) return null;
		
		String[] se = Misc.AMP.split(vcfss);
		//for each one parse the delta and gene name
		//SAMD11:D5S,939461,ACGGTGAGG,ACATGTACT,9.2,-27.1,9.2
		double maxDelta = -1;
		String maxGeneName = null;
		for (String effect: se){
			String[] tokens = Misc.COLON.split(effect);
			//for each of the user defined effects
			for (String code: passingVCFSS){
				if (tokens[1].startsWith(code)){
					String[] items = Misc.COMMA.split(tokens[1]);
					double delta = Double.parseDouble(items[6]);
					if (delta > maxDelta) {
						maxDelta = delta;
						maxGeneName = tokens[0];
					}
				}
			}
		}
		// check score and that there aren't two gene names!
		if (maxDelta >= minimumVCFSSDiff && maxGeneName.contains(",") == false) return maxGeneName;
		return null;
	}
	
	private String fetchClinvarPathoGene(String info) {
		//is there a CLNSIG=Pathogenic or CLNSIG=Likely_pathogenic? Then return the associated gene
		if (info.contains("CLNSIG=Pathogenic") || info.contains("CLNSIG=Likely_pathogenic")) {
			Matcher mat = geneInfo.matcher(info);
			if (mat.matches()) return mat.group(1);
		}
		return null;
	}

	/**Walks the LHM finding the key with the max counts. For those with the same # counts, the first is returned.*/
	public static String findFirstLargest(LinkedHashMap<String, Integer> nameCount) {
		int max = 0;
		String name = null;
		for (String s: nameCount.keySet()) {
			int count = nameCount.get(s);
			if (count > max) {
				max = count;
				name = s;
			}
		}
		return name;
	}
	
	/**
	 * ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
	 * 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
	 *    0          1                2               3          4            5
	 * ANN=GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20),
	 * GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20)
	 * 
	 * */
	public static LinkedHashMap<String, Integer> fetchAffectedGenes(String[] splitInfo) throws IOException {
		//find ANN=
		String annValue = null;
		for (String i: splitInfo){
			if (i.startsWith("ANN=")) {
				annValue = i;
				break;
			}
		}
		if (annValue == null) throw new IOException("Failed to find the ANN= in "+Misc.stringArrayToString(splitInfo, " "));

		//drop ANN=
		annValue = annValue.substring(4);

		//split ANN by comma
		String[] splitAnns = Misc.COMMA.split(annValue);

		//Save observed genes, but don't double count!
		LinkedHashMap<String, Integer> geneNameHits = new LinkedHashMap<String, Integer>();
		for (String ann: splitAnns){
			//split ANN and save gene name
			String[] splitAnn = Misc.PIPE.split(ann);
			Integer i = geneNameHits.get(splitAnn[3]);
			if (i == null) i = new Integer(1);
			else i = new Integer (i.intValue()+1);
			geneNameHits.put(splitAnn[3],i);
		}
		
		return geneNameHits;
		
	}
	
	/**
	 * ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
	 * 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
	 *    0          1                2               3          4            5
	 * ANN=GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20),
	 * GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20)
	 * 
	 * */
	public static String fetchFirstAffectedGene(String[] splitInfo) throws IOException {
		//find ANN=
		String annValue = null;
		for (String i: splitInfo){
			if (i.startsWith("ANN=")) {
				annValue = i;
				break;
			}
		}
		if (annValue == null) throw new IOException("Failed to find the ANN= in "+Misc.stringArrayToString(splitInfo, " "));

		//drop ANN=
		annValue = annValue.substring(4);

		//split ANN by comma
		String[] splitAnns = Misc.COMMA.split(annValue);
			
		//split first ANN and return it's gene name
		String[] splitAnn = Misc.PIPE.split(splitAnns[0]);
		return splitAnn[3];
		
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFGeneFreqJointGenoComp(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcfFile = new File(args[++i]); break;
					case 'd': minReadDepth = Double.parseDouble(args[++i]); break;
					//case 'a': alleleFreq = Double.parseDouble(args[++i]); break;
					case 'a': groupAPattern = Pattern.compile(args[++i]); break;
					case 'b': groupBPattern = Pattern.compile(args[++i]); break;
					case 'q': qual = Double.parseDouble(args[++i]); break;
					case 'g': genotypeQuality = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (vcfFile.canRead() == false) Misc.printErrAndExit("\nPlease provide the path to a joint genotyped, multi sample vcf file, "+vcfFile);
		if (groupAPattern == null || groupBPattern == null) Misc.printErrAndExit("\nPlease both a group A '"+groupAPattern+"' and group B '"+groupBPattern+"' java pattern to match sample names");
	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         VCF Gene Freq Joint Geno Comp : Aug 2020                 **\n" +
				"**************************************************************************************\n" +
				"For each affected gene in a multi sample vcf file, this tool contrasts the number of \n"+
				"mutant and wild type samples using a 2x2 two-tailed fisher's exact test. Run the \n"+
				"https://tinyurl.com/yygmj4oq workflow to annotate a joint genotyped, multi sample vcf\n"+
				"with SnpEff, Clinvar, and VCFSpliceScanner info. Then filter the vcf file with the\n"+
				"USeq AnnotatedVcfParser app to select variants of interest before running this app.\n" +

				"\nRequired Params:\n"+
				"-v Path to a multi-sample, joint genotyped, filtered vcf file (xxx.vcf(.gz/.zip OK))\n"+
				"-a Group A sample name pattern, must match entire name, see Java Pattern class\n"+
				"-a Group B sample name pattern, must match entire name, see Java Pattern class\n"+
				"\n"+
				"-q Minimum VCF QUAL, defaults to 20\n"+
				"-d Minimum sample read depth, defaults to 20\n"+
                "-g Minimum sample genotype quality, defaults to 20\n"+
				
				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VCFGeneFreqJointGenoComp -v jg.vcf.gz\n" +
				"       -d 20 -a '^M_.+' -b '^1.+' \n"+

		"\n**************************************************************************************\n");

	}
}
