package edu.utah.seq.vcf.json;

import java.io.*;
import java.util.regex.*;
import org.json.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;

/**
 * Takes one or more patient json reports from Tempus tests and converts the variants into vcf format. This is a mix of somatic, inherited, snv/indel, and cnvs.
 * @author david.nix@hci.utah.edu 
 **/
public class TempusJson2Vcf {

	//user defined fields
	private File[] jsonFiles = null;
	private File indexedFasta = null;
	private File saveDirectory = null;
	private String acceptedSchema = "1.3";
	
	//internal fields
	private IndexedFastaSequenceFile fasta = null;
	private String source = null;
	private HashMap<String, Bed> cnvGeneNameBed = null;
	
	//counters across all datasets
	TreeMap<String, Integer> bioInfPipelines = new TreeMap<String, Integer>();
	TreeMap<String, Integer> reportStatus = new TreeMap<String, Integer>();
	TreeMap<String, Integer> physicians = new TreeMap<String, Integer>();
	TreeMap<String, Integer> testCodes = new TreeMap<String, Integer>();
	TreeMap<String, Integer> testDescriptions = new TreeMap<String, Integer>();
	TreeMap<String, Integer> diagnosis = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleCategories = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleSites = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleTypes = new TreeMap<String, Integer>();
	
	//counters from TempusVariants
	TreeMap<String, Integer> genes = new TreeMap<String, Integer>();
	TreeMap<String, Integer> referenceGenome = new TreeMap<String, Integer>();
	TreeMap<String, Integer> variantType = new TreeMap<String, Integer>();
	TreeMap<String, Integer> variantDescription = new TreeMap<String, Integer>();
	
	Histogram tumorPercentages = new Histogram(0,100, 20);
	Histogram tumorAF = new Histogram(0,100, 20);
	Histogram tumorDP = new Histogram(0,5000, 50);
	Histogram ageAtDiagnosis = new Histogram(0,100, 20);
	
	private int numSomaticPotentiallyActionableMutations = 0;
	private int numSomaticVariantsOfUnknownSignificance = 0;
	private int numSomaticBiologicallyRelevantVariants = 0;
	private int numSomaticPotentiallyActionableCopyNumberVariants = 0;
	private int numInheritedRelevantVariants = 0;
	private int numInheritedIncidentalFindings = 0;
	private int numInheritedVariantsOfUnknownSignificance = 0;
	private int numFusionVariants = 0;
	
	//working data for a particular report
	private File workingJsonFile = null;
	private TempusReport workingReport = null;
	private TempusPatient workingPatient = null;
	private TempusOrder workingOrder = null;
	private TempusSpecimen[] workingSpecimens = null;
    private TempusResults workingResults = null;
	private int workingNumSomaticPotentiallyActionableMutations = 0;
	private int workingNumSomaticVariantsOfUnknownSignificance = 0;
	private int workingNumSomaticBiologicallyRelevantVariants = 0;
	private int workingNumSomaticPotentiallyActionableCopyNumberVariants = 0;
	private int workingNumInheritedRelevantVariants = 0;
	private int workingNumInheritedIncidentalFindings = 0;
	private int workingNumInheritedVariantsOfUnknownSignificance = 0;
	private int workingNumFusionVariants = 0;
	
	//for the vcf header
	private static String altCNV = "##ALT=<ID=CNV,Description=\"Copy number variable region\">";
	private static String infoEG = "##INFO=<ID=EG,Number=1,Type=String,Description=\"Effected gene\">";
	private static String infoCL = "##INFO=<ID=CL,Number=1,Type=String,Description=\"Tempus classification\">";
	private static String infoFE = "##INFO=<ID=FE,Number=1,Type=String,Description=\"Functional effect on gene\">";
	private static String infoDP = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
	private static String infoAF = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
	private static String infoIMPRECISE = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">";
	private static String infoSVTYPE = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">";
	private static String infoSVLEN = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles. Negative values for deletions, positive for amplifications.\">";
	private static String infoEND = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">";
	private static String infoCN = "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy Number\">";


	//constructors
	public TempusJson2Vcf(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();
			
			printStats();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running Tempus2Vcf app!");
		}
	}
	
	public void doWork() throws Exception{
		
		//Create fasta fetcher
		if (indexedFasta != null) {
			fasta = new IndexedFastaSequenceFile(indexedFasta);
			if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);
		}

		IO.pl("Parsing and coverting...\n");
		IO.pl("Name\t"
				+ "NumSomPotActMuts\t"
				+ "NumSomVUS\t"
				+ "NumSomBioRelVars\t"
				+ "NumSomPotActCNVs\t"
				+ "NumInhRelVars\t"
				+ "NumInhIncFinds\t"
				+ "NumInhVUS\t"
				+ "NumFusVars");
		
		for (int i=0; i< jsonFiles.length; i++){
			
			//process file
			workingJsonFile = jsonFiles[i];
			IO.p(workingJsonFile.getName());
			convert();
			
			//build vcf
			writeVcf();
			
			//finish line and output counters for this json file
			IO.pl("\t"+
					workingNumSomaticPotentiallyActionableMutations+"\t"+
					workingNumSomaticVariantsOfUnknownSignificance+"\t"+
					workingNumSomaticBiologicallyRelevantVariants+"\t"+
					workingNumSomaticPotentiallyActionableCopyNumberVariants+"\t"+
					workingNumInheritedRelevantVariants+"\t"+
					workingNumInheritedIncidentalFindings+"\t"+
					workingNumInheritedVariantsOfUnknownSignificance+"\t"+
					workingNumFusionVariants);
			
			resetWorkingCounters();
		}
		
		//close the fasta lookup fetcher
		if (indexedFasta != null) fasta.close();

	}

	private void resetWorkingCounters() {
		//increment master counters
		numSomaticPotentiallyActionableMutations+= workingNumSomaticPotentiallyActionableMutations;
		numSomaticVariantsOfUnknownSignificance+= workingNumSomaticVariantsOfUnknownSignificance;
		numSomaticBiologicallyRelevantVariants+= workingNumSomaticBiologicallyRelevantVariants;
		numSomaticPotentiallyActionableCopyNumberVariants+= workingNumSomaticPotentiallyActionableCopyNumberVariants;
		numInheritedRelevantVariants+= workingNumInheritedRelevantVariants;
		numInheritedIncidentalFindings+= workingNumInheritedIncidentalFindings;
		numInheritedVariantsOfUnknownSignificance+= workingNumInheritedVariantsOfUnknownSignificance;
		numFusionVariants += workingNumFusionVariants;
		
		//reset working counters
		workingNumSomaticPotentiallyActionableMutations = 0;
		workingNumSomaticVariantsOfUnknownSignificance = 0;
		workingNumSomaticBiologicallyRelevantVariants = 0;
		workingNumSomaticPotentiallyActionableCopyNumberVariants = 0;
		workingNumInheritedRelevantVariants = 0;
		workingNumInheritedIncidentalFindings = 0;
		workingNumInheritedVariantsOfUnknownSignificance = 0;
		workingNumFusionVariants = 0; 
		
		workingJsonFile = null;
		workingReport = null;
		workingPatient = null;
		workingOrder = null;
		workingSpecimens = null;
	    workingResults = null;
	}

	private void writeVcf() {
		
		try {
			String name = Misc.removeExtension(workingJsonFile.getName());
			Gzipper out = new Gzipper(new File (saveDirectory, name+".vcf.gz"));
			//Gzipper outTxt = new Gzipper(new File (saveDirectory, name+".txt.gz"));
			
			//add header
			out.println(buildVcfHeader());
			
			//add sorted vcfs
			ArrayList<Bed> bedVcfs = new ArrayList<Bed>();
			int counter = 0;
			for (TempusVariant tv: workingResults.getVariants()) {
				//outTxt.println(tv.toString());
				String vcf = tv.toVcf(counter);
				if (vcf != null) {
					counter++;
					bedVcfs.add(fetchBed(vcf));
				}
			}
			Bed[] b = new Bed[bedVcfs.size()];
			bedVcfs.toArray(b);
			Arrays.sort(b);
			for (Bed t: b) out.println(t.getName());
			
			out.close();
			//outTxt.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: issue writing out vcf for "+workingJsonFile);
		}
	}
	
	private Bed fetchBed(String vcf){
		String[] t = Misc.TAB.split(vcf);
		int pos = Integer.parseInt(t[1]);
		//String chromosome, int start, int stop, String name, double score, char strand
		Bed b = new Bed(t[0], pos, pos+1, vcf, 0, '.' );
		return b;
	}

	private String buildVcfHeader(){
		StringBuilder sb = new StringBuilder();
		
		//add in standard header
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("##source=\""+source+"\"\n");
		sb.append("##file-path="+workingJsonFile+"\n");
		sb.append("##parse-date="+Misc.getDateNoSpaces()+"\n");
		
		//add in meta data
		LinkedHashMap<String,String> meta = new LinkedHashMap<String, String>();
		workingReport.addMetaData(meta);
		workingPatient.addMetaData(meta);
		workingOrder.addMetaData(meta);
		for (int i=0; i<workingSpecimens.length; i++)workingSpecimens[i].addMetaData(meta, i);
		workingResults.addMetaData(meta);
		for (TempusVariant tv: workingResults.getVariants()) {
			if (tv.getReferenceGenome() != null) meta.put("referenceGenome", tv.getReferenceGenome());
		}
		for (String key: meta.keySet()){
			String value = meta.get(key).trim();
			sb.append("##");
			sb.append(key);
			//any whitespace?
			Matcher mat = Misc.WHITESPACE.matcher(value);
			if (mat.find()){
				sb.append("=\"");
				sb.append(value);
				sb.append("\"\n");
			}
			else {
				sb.append("=");
				sb.append(value);
				sb.append("\n");
			}
		}
		
		//add info
		sb.append(altCNV+"\n");
		sb.append(infoEG+"\n");
		sb.append(infoCL+"\n");
		sb.append(infoFE+"\n");
		sb.append(infoDP+"\n");
		sb.append(infoAF+"\n");
		sb.append(infoIMPRECISE+"\n");
		sb.append(infoSVTYPE+"\n");
		sb.append(infoSVLEN+"\n");
		sb.append(infoEND+"\n");
		sb.append(infoCN+"\n");

		//chrom line
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		return sb.toString();
	}
	
	private void convert() {
		try {
			
			String jString = IO.loadFile(workingJsonFile, " ", true);
	        JSONObject object = new JSONObject(jString);
	        
	        //schema
	        JSONObject meta = object.getJSONObject("metadata");
	        String jsonSchema = Json.getStringAttribute(meta, "schemaVersion");
	        if (jsonSchema == null || jsonSchema.equals(acceptedSchema) == false) Misc.printErrAndExit("\nIncorrect schema! Aborting. Only works with "+ acceptedSchema+" found "+jsonSchema);
	        
	        //report
	        workingReport = new TempusReport(object, this);
	        
	        //patient
	        workingPatient = new TempusPatient(object, this);
	        
	        //order
	        workingOrder = new TempusOrder(object, this);

	        //specimens, should be just two!
	        workingSpecimens = TempusSpecimen.getSpecimens(object, this);
	        
	        //results
	        workingResults = new TempusResults(object, this); 
	        
	        //check variant ref bps
	        if (indexedFasta != null) checkRefSeqs();
	        
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing this json file "+workingJsonFile);
		}
	}	
	
	private void checkRefSeqs() throws Exception {
		for (TempusVariant tv: workingResults.getVariants()) {
			if (tv.getRef()!=null && tv.getChromosome()!=null && tv.getPos()!=null) {
				int pos = Integer.parseInt(tv.getPos());
				ReferenceSequence rs = fasta.getSubsequenceAt(tv.getChromosome(), pos, pos+tv.getRef().length()-1);
				String test = new String(rs.getBases());
				if (test.equals(tv.getRef()) == false) throw new Exception("\tWARNING: variant ref does not match seq from fasta ('"+test+"')?\n"+tv);
			}
		}
	}

	public void printStats(){
		IO.pl("\nSummary Stats:\n");
		
		IO.pl(jsonFiles.length+ "\tNum Json files parsed");
		
		IO.pl(numSomaticPotentiallyActionableMutations+  			"\tNum Somatic Potentially Actionable Mutations");
		IO.pl(numSomaticVariantsOfUnknownSignificance+   			"\tNum Somatic Variants Of Unknown Significance");
		IO.pl(numSomaticBiologicallyRelevantVariants+   			"\tNum Somatic Biologically Relevant Variants");
		IO.pl(numSomaticPotentiallyActionableCopyNumberVariants+  	"\tNum Somatic Potentially Actionable Copy Number Variants");
		
		IO.pl(numInheritedRelevantVariants+              "\tNum Inherited Relevant Variants");
		IO.pl(numInheritedIncidentalFindings+            "\tNum Inherited Incidental Findings");
		IO.pl(numInheritedVariantsOfUnknownSignificance+ "\tNum Inherited Variants Of Unknown Significance");
		IO.pl(numFusionVariants+						 "\tNum Fusion Variants");
		
		IO.pl("\nBioInf Pipelines");
		Misc.printTreeMap(bioInfPipelines, "\t", "\t");
		IO.pl("\nPhysicians:");
		Misc.printTreeMap(physicians, "\t", "\t");
		IO.pl("\nReport Status (qns is failed):");
		Misc.printTreeMap(reportStatus, "\t", "\t");
		IO.pl("\nTest Codes:");
		Misc.printTreeMap(testCodes, "\t", "\t");
		IO.pl("\nTest Descriptions:");
		Misc.printTreeMap(testDescriptions, "\t", "\t");
		IO.pl("\nDiagnosis:");
		Misc.printTreeMap(diagnosis, "\t", "\t");
		IO.pl("\nSample Categories:");
		Misc.printTreeMap(sampleCategories, "\t", "\t");
		IO.pl("\nSample Sites:");
		Misc.printTreeMap(sampleSites, "\t", "\t");
		IO.pl("\nSample Types:");
		Misc.printTreeMap(sampleTypes, "\t", "\t");
		
		IO.pl("\nGene Mutations:");
		Misc.printTreeMap(genes, "\t", "\t");
		IO.pl("\nReference Genome:");
		Misc.printTreeMap(referenceGenome, "\t", "\t");
		IO.pl("\nVariant Type:");
		Misc.printTreeMap(variantType, "\t", "\t");
		IO.pl("\nVariant Description:");
		Misc.printTreeMap(variantDescription, "\t", "\t");
		
		IO.pl("\nPercent Tumor in Sample");
		tumorPercentages.printScaledHistogram();
		IO.pl("\nTumor Variant Allele Percentages");
		tumorAF.printScaledHistogram();
		IO.pl("\nTumor Variant Read Depth");
		tumorDP.printScaledHistogram();
		IO.pl("\nAge at Diagnosis");
		ageAtDiagnosis.printScaledHistogram();
		
	}
	
	public static void add(String key, TreeMap<String, Integer> map){
		if (key == null) return;
		int count = 0;
		if (map.containsKey(key)) count = map.get(key);
		map.put(key, ++count);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TempusJson2Vcf(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		source = useqVersion+" Args: USeq/TempusJson2Vcf "+ Misc.stringArrayToString(args, " ");
		System.out.println("\n"+ source +"\n");
		File forExtraction = null;
		File bed = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'j': forExtraction = new File(args[++i]); break;
					case 'f': indexedFasta = new File(args[++i]); break;
					case 'b': bed = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull json files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a Tempus json file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".json");
		tot[1] = IO.extractFiles(forExtraction,".json.gz");
		tot[2] = IO.extractFiles(forExtraction,".json.zip");
		jsonFiles = IO.collapseFileArray(tot);
		if (jsonFiles == null || jsonFiles.length ==0 || jsonFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.json(.zip/.gz OK) file(s)!\n");

		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save directory!\n"+saveDirectory);
		saveDirectory.mkdirs();
		if (saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: your save directory does not appear to be a directory?\n");

		if (jsonFiles == null || jsonFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any json files to parse?!\n");
		
		//load hash of gene name and Bed to use in coordinates for CNVs
		if (bed != null) {
			if (indexedFasta == null) Misc.printErrAndExit("\nError: addition of CNV info requires an indexed fasta\n");
			cnvGeneNameBed = Bed.parseBedForNames(bed);
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Tempus Json 2 Vcf: March 2019                        **\n" +
				"**************************************************************************************\n" +
				"Parses json Tempus reports to vcf. Leave in PHI to enable calculating age at\n"+
				"diagnosis. Summary statistics calculated for all reports. Vcfs will contain a mix of \n"+
				"somatic and inherited snvs, indels, and cnvs.\n"+

				"\nOptions:\n"+
				"-j Path to Tempus json report or directory containing such, xxx.json(.gz/.zip OK)\n"+
				"-s Path to a directory for saving the results.\n"+
				"-b Path to a bed file for converting CNV gene names to coordinates where the name\n"+
				"     column contains just the gene name.\n"+
				"-f Path to the reference fasta with xxx.fai index. Required for CNV conversions.\n"+
				
				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/TempusJson2Vcf -j /F1/TempusJsons\n" +
				"     -f /Ref/human_g1k_v37.fasta -s /F1/VCF/ -b /Ref/b37TempusGeneRegions.bed.gz \n\n" +

				"**************************************************************************************\n");

	}

	public void setWorkingNumSomaticPotentiallyActionableMutations(int workingNumSomaticPotentiallyActionableMutations) {
		this.workingNumSomaticPotentiallyActionableMutations = workingNumSomaticPotentiallyActionableMutations;
	}

	public void setWorkingNumSomaticVariantsOfUnknownSignificance(int workingNumSomaticVariantsOfUnknownSignificance) {
		this.workingNumSomaticVariantsOfUnknownSignificance = workingNumSomaticVariantsOfUnknownSignificance;
	}

	public void setWorkingNumInheritedRelevantVariants(int workingNumInheritedRelevantVariants) {
		this.workingNumInheritedRelevantVariants = workingNumInheritedRelevantVariants;
	}

	public void setWorkingNumInheritedIncidentalFindings(int workingNumInheritedIncidentalFindings) {
		this.workingNumInheritedIncidentalFindings = workingNumInheritedIncidentalFindings;
	}

	public void setWorkingNumInheritedVariantsOfUnknownSignificance(int workingNumInheritedVariantsOfUnknownSignificance) {
		this.workingNumInheritedVariantsOfUnknownSignificance = workingNumInheritedVariantsOfUnknownSignificance;
	}

	public void setWorkingNumSomaticBiologicallyRelevantVariants(int workingNumSomaticBiologicallyRelevantVariants) {
		this.workingNumSomaticBiologicallyRelevantVariants = workingNumSomaticBiologicallyRelevantVariants;
	}

	public void setWorkingNumSomaticPotentiallyActionableCopyNumberVariants(
			int workingNumSomaticPotentiallyActionableCopyNumberVariants) {
		this.workingNumSomaticPotentiallyActionableCopyNumberVariants = workingNumSomaticPotentiallyActionableCopyNumberVariants;
	}

	public void setWorkingNumFusionVariants(int workingNumFusionVariants) {
		this.workingNumFusionVariants = workingNumFusionVariants;
	}

	public Histogram getAgeAtDiagnosis() {
		return ageAtDiagnosis;
	}

	public TempusReport getWorkingReport() {
		return workingReport;
	}

	public TempusPatient getWorkingPatient() {
		return workingPatient;
	}

	public TempusOrder getWorkingOrder() {
		return workingOrder;
	}

	public TempusSpecimen[] getWorkingSpecimens() {
		return workingSpecimens;
	}

	public TempusResults getWorkingResults() {
		return workingResults;
	}

	public HashMap<String, Bed> getCnvGeneNameBed() {
		return cnvGeneNameBed;
	}

	public IndexedFastaSequenceFile getFasta() {
		return fasta;
	}
	
}
