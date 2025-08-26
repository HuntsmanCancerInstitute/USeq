package edu.utah.seq.vcf.json.tempusv3;

import java.io.*;
import java.util.regex.*;
import org.json.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;

/**
 * Takes one or more patient json v3.3+ reports from Tempus tests and converts the variants into vcf format. This is a mix of somatic, inherited, snv/indel, and cnvs.
 * @author david.nix@hci.utah.edu 
 **/
public class TempusV3Json2Vcf {

	//user defined fields
	private File[] jsonFiles = null;
	private File[] vcfFiles = null;
	private File geneAliasFile = null;
	private File indexedFasta = null;
	private File saveDirectory = null;
	private boolean includePHI = false;
	private LinkedHashSet<String> keysToExport = null;
	
	//internal fields
	private static String acceptedSchema = "v3."; 
	private IndexedFastaSequenceFile fasta = null;
	private String source = null;
	private HashMap<String, Bed> cnvGeneNameBed = null;
	private String jsonSchema = null;
	private ArrayList<TempusV3JsonSummary> summaries = new ArrayList<TempusV3JsonSummary>();
	private LinkedHashMap<String,String>[] allReportAttributes = null;
	private HashMap<String, String[]> geneAliases = null;
	private ArrayList<TempusV3Variant> failedToFindCooridinates = new ArrayList<TempusV3Variant>();
	private HashMap<String, TempusV3JsonCollection> orderCollections = null;
	
	//counters across all datasets
	TreeMap<String, Integer> bioInfoPipeline = new TreeMap<String, Integer>();
	TreeMap<String, Integer> reportStatus = new TreeMap<String, Integer>();
	TreeMap<String, Integer> physicians = new TreeMap<String, Integer>();
	TreeMap<String, Integer> testCodes = new TreeMap<String, Integer>();
	TreeMap<String, Integer> testDescriptions = new TreeMap<String, Integer>();
	TreeMap<String, Integer> diagnosis = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleCategories = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleSites = new TreeMap<String, Integer>();
	
	//counters from TempusVariants
	TreeMap<String, Integer> genes = new TreeMap<String, Integer>();
	TreeMap<String, Integer> genomicSourceClass = new TreeMap<String, Integer>();
	TreeMap<String, Integer> variantType = new TreeMap<String, Integer>();
	TreeMap<String, Integer> variantDescription = new TreeMap<String, Integer>();
	
	Histogram tumorPercentages = new Histogram(0,100, 20);
	Histogram tumorAF = new Histogram(0,100, 20);
	Histogram ageAtDiagnosis = new Histogram(0,100, 20);
	
	private int numPotentiallyActionable = 0;
	private int numBiologicallyRelevant = 0;
	private int numPathogenic = 0;
	private int numLikelyPathogenic = 0;
	private int numRiskAllele = 0;
	private int numUnknownSignificance = 0;
	
	//working data for a particular report
	private LinkedHashMap<String,String> reportAttributes = null;
	private File workingJsonFile = null;
	private JSONObject workingJsonObject = null;
	private TempusV3Report workingReport = null;
	private TempusV3Patient workingPatient = null;
	private TempusV3Order workingOrder = null;
	private TempusV3Specimen[] workingSpecimens = null;
    private TempusV3GenomicVariants workingResults = null;

	private int workingNumPotentiallyActionable = 0;
	private int workingNumBiologicallyRelevant = 0;
	private int workingNumPathogenic = 0;
	private int workingNumLikelyPathogenic = 0;
	private int workingNumRiskAllele = 0;
	private int workingNumUnknownSignificance = 0;
    
	private String[] workingSomVcfLines = null;
	private String[] workingGermVcfLines = null;
	
	//for the vcf header
	private static String altCNV = "##ALT=<ID=CNV,Description=\"Copy number variable region\">";
	private static String infoEG = "##INFO=<ID=EG,Number=1,Type=String,Description=\"Effected gene(s)\">";
	private static String infoCL = "##INFO=<ID=CL,Number=1,Type=String,Description=\"Tempus classifications: when a normal is provided, germline variants are 'pathogenic', 'likelyPathogenic', 'riskAllele'; somatic variants are 'potentiallyActionable', 'biologicallyRelevant', 'unknownSignificance' \">";
	private static String infoFG = "##INFO=<ID=FG,Number=1,Type=String,Description=\"Fused gene from a structural variant\">";
	private static String infoFE = "##INFO=<ID=FE,Number=1,Type=String,Description=\"Functional effect on gene\">";
	private static String infoAF = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
	private static String infoIMPRECISE = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">";
	private static String infoSVTYPE = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">";
	private static String infoSVLEN = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles. Negative values for deletions, positive for amplifications.\">";
	private static String infoEND = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">";
	private static String infoDESC = "##INFO=<ID=DESC,Number=1,Type=String,Description=\"Description of the type of rearrangement observed\">";
	private static String infoMATEID = "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of the mate vcf record\">";


	//constructors
	public TempusV3Json2Vcf(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			printStats();
			
			printSpreadsheet();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running TempusV3Json2Vcf app!");
		}
	}
	
	public void doWork() throws Exception{

		//Create fasta fetcher
		if (indexedFasta != null) {
			fasta = new IndexedFastaSequenceFile(indexedFasta);
			if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);
		}

		//create container for each report attributes
		allReportAttributes = new LinkedHashMap[jsonFiles.length];

		IO.pl("Parsing jsons...");
		

		//parse all the single json files
		for (int i=0; i< jsonFiles.length; i++){

			//process file
			workingJsonFile = jsonFiles[i];
			IO.p("\t"+workingJsonFile.getName());

			if (convert()) {
				summaries.add(new TempusV3JsonSummary(workingJsonObject, workingOrder, workingPatient, workingReport, workingSpecimens, workingResults));
				resetWorkingCounters();
				allReportAttributes[i] = reportAttributes;
			}
			else IO.pl();
		}
		
		//group summaries by tempusOrderId, cannot do by tempus patient id, these are not system wide
		groupSummaries();
		
		//save vcfs
		IO.pl("\nWriting merged json VCFs... ");
		for (TempusV3JsonCollection c: orderCollections.values()) writeVcfs(c);

		//close the fasta lookup fetcher
		if (indexedFasta != null) fasta.close();

	}
	
	private void groupSummaries() {
		IO.pl("\nGrouping summaries by 'tempusOrderId' all of which are focused on one tumor sample...");
		orderCollections = new HashMap<String, TempusV3JsonCollection>();
		for (TempusV3JsonSummary s: summaries) {
			String tempusOrderId = s.getTempusOrder().getTempusOrderId();
			TempusV3JsonCollection c = orderCollections.get(tempusOrderId);
			if (c==null) {
				c = new TempusV3JsonCollection(tempusOrderId);
				orderCollections.put(tempusOrderId, c);
			}
			c.getJsonSummaries().add(s);
		}
		IO.pl("\t"+orderCollections.size()+"\tOrders");
		
		IO.pl("\nDepreciate duplicate summaries with an appended report...");
		for (String tempusOrderId: orderCollections.keySet()) {
			IO.pl("TempusOrderId: "+tempusOrderId);

			TempusV3JsonCollection collection = orderCollections.get(tempusOrderId);
			HashMap<String, ArrayList<TempusV3JsonSummary>> codeSummary = new HashMap<String, ArrayList<TempusV3JsonSummary>>();
			//split by test code
			for (TempusV3JsonSummary sum: collection.getJsonSummaries()) {
				String code = sum.getTempusOrder().getTestCode();
				ArrayList<TempusV3JsonSummary> al = codeSummary.get(code);
				if (al==null) {
					al = new ArrayList<TempusV3JsonSummary>();
					codeSummary.put(code, al);
				}
				al.add(sum);
			}

			//look to see if any with the same test code have a report status of 'amendment'
			for (String code: codeSummary.keySet()) {
				ArrayList<TempusV3JsonSummary> al = codeSummary.get(code);
				//duplicates
				if (al.size()>1) {
					//any amendment?
					TempusV3JsonSummary toKeep = null;
					for (TempusV3JsonSummary sum: al) {
						if (sum.getTempusReport().getReportStatus().contains("amendment")) {
							toKeep = sum;
							break;
						}
					}
					//toKeep assigned? nuke others
					if (toKeep != null) {
						for (TempusV3JsonSummary sum: al) {
							if (sum.getTempusReport().getReportStatus().contains("amendment") == false) {
								IO.pl("\tExcluding "+sum.getTempusOrder().getAccessionId());
								sum.setDepreciated(true);
							}
						}
					}
				}
			}
		}

		//print stats
		for (String tempusOrderId: orderCollections.keySet()) {
			IO.pl("TempusOrderId: "+tempusOrderId);
			//get the collection
			TempusV3JsonCollection collection = orderCollections.get(tempusOrderId);

			ArrayList<TempusV3JsonSummary> al = collection.getJsonSummaries();
			IO.p("\t"+al.size()+"\tReports: ");
			for (TempusV3JsonSummary sum: al) {
				TempusV3Order order = sum.getTempusOrder();
				IO.p(order.getAccessionId()+ "_"+ order.getTempusOrderId());
				IO.p(" ");
			}
			IO.pl();

		}
	}
	
	private void loadGeneAliases() throws IOException {
		geneAliases = new HashMap<String, String[]>();
		BufferedReader in = IO.fetchBufferedReader(geneAliasFile);
		String line = null;
		String[] fields = null;
		HashSet<String> aliases = new HashSet<String>();
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.startsWith("#") || line.length()==0) continue;
			fields = Misc.TAB.split(line);
			if (fields.length==1) {
				geneAliases.put(line, new String[] {line});
			}
			else {
				aliases.clear();
				fields[0] = fields[0].trim();
				//only add it if a length of 2 or more
				if (fields[0].length()>1) aliases.add(fields[0]);
				for (int i=1; i< fields.length; i++) {
					String[] subFields = Misc.PIPE.split(fields[i]);
					for (String sf: subFields) {
						sf = sf.trim();
						if (sf.length()>1) aliases.add(sf);
					}
				}
				String[] allAliases = Misc.hashSetToStringArray(aliases);
				for (String s: allAliases) {
					//only add if not present
					if (geneAliases.containsKey(s) == false) geneAliases.put(s, allAliases);
				}
			}
			
		}
		in.close();
	}

	/**PMR ID, Patient Molecular Repo ID, HCI patient identifier, will NOT be present in a new json from Tempus*/
	private void addMDPID() throws IOException {
		String mdpid = Misc.fetchMDPID(workingJsonFile);
		if (mdpid != null) reportAttributes.put("MolecularDataPatientId", mdpid);
	}
	
	private void printSpreadsheet() throws IOException {
		IO.pl("\nExporting patient summary spreadsheet... ");
		//merge all the keys
		LinkedHashSet<String> allKeys = new LinkedHashSet<String>();
		for (LinkedHashMap<String,String> s : allReportAttributes) if (s!=null) allKeys.addAll(s.keySet());
		IO.pl("\tAvailable attributes: "+allKeys);
		
		if (keysToExport == null) keysToExport = allKeys;
		else {
			//check that the keys they want are part of the allKeys
			boolean failed = false;
			for (String k: keysToExport) {
				if (allKeys.contains(k) == false) {
					failed = true;
					IO.pl("\tFailed to find "+k);
				}
			}
			if (failed) {
				IO.pl("\tAvailable attributes: "+allKeys);
				throw new IOException("\nError printing requested attributes to spreadsheet.");
			}
		}
		
		//make writer add header
		String phi = "";
		if (includePHI) phi = "_PHI";
		PrintWriter txtOut = new PrintWriter (new FileWriter(new File (saveDirectory, "aggregatePatientInfo"+phi+".xls")));
		for (String k: keysToExport) {
			txtOut.print(k);
			txtOut.print("\t");
		}
		txtOut.println();
		
		//write a line per sample
		for (LinkedHashMap<String,String> s : allReportAttributes) {
			if (s== null) continue;
			//for each key
			for (String k: keysToExport) {
				String v = s.get(k);
				if (v!=null) txtOut.print(v);
				txtOut.print("\t");
			}
			txtOut.println();
		}
		
		txtOut.close();
	}

	private void resetWorkingCounters() {
		//increment master counters
		numPotentiallyActionable += workingNumPotentiallyActionable;
		numBiologicallyRelevant += workingNumBiologicallyRelevant;
		numPathogenic += workingNumPathogenic;
		numLikelyPathogenic += workingNumLikelyPathogenic;
		numRiskAllele += workingNumRiskAllele;
		numUnknownSignificance += workingNumUnknownSignificance; 
		
		//reset working counters
		workingNumPotentiallyActionable = 0;
		workingNumBiologicallyRelevant = 0;
		workingNumPathogenic = 0;
		workingNumLikelyPathogenic = 0;
		workingNumRiskAllele = 0;
		workingNumUnknownSignificance = 0; 
		
		workingJsonFile = null;
		workingReport = null;
		workingPatient = null;
		workingOrder = null;
		workingSpecimens = null;
	    workingResults = null;
	}

	private Bed fetchBed(String vcf){
		String[] t = Misc.TAB.split(vcf);
		int pos = Integer.parseInt(t[1]);
		//String chromosome, int start, int stop, String name, double score, char strand
		Bed b = new Bed(t[0], pos, pos+1, vcf, 0, '.' );
		return b;
	}
	
	private void writeVcfs(TempusV3JsonCollection collection) throws Exception {

		ArrayList<TempusV3JsonSummary> al = collection.getJsonSummaries();
		//collect accessionIDs
		String[] ids = new String[al.size()];
		int i=0;
		for (TempusV3JsonSummary sum: al) ids[i++] = sum.getTempusOrder().getAccessionId();
		String combIds = Misc.stringArrayToString(ids, "_");

		Gzipper out = new Gzipper(new File (saveDirectory, combIds+".vcf.gz"));
		//add header
		out.println(fetchVcfHeader(al));
		ArrayList<Bed> bedVcfs = new ArrayList<Bed>();
		int counter = 0;

		for (TempusV3JsonSummary sum: al) { 
			//skip depreciated tests due to the existance of an amendment
			if (sum.isDepreciated()) continue;
			for (TempusV3Variant tv: sum.getTempusV3Results().getVariants()) {
				String vcf = tv.toVcf(counter);
				if (vcf != null) {
					counter++;
					bedVcfs.add(fetchBed(vcf));
				}
			}
		}
		//sort the variants by position and print to file, might be none so an empty vcf
		Bed[] b = new Bed[bedVcfs.size()];
		bedVcfs.toArray(b);
		Arrays.sort(b);
		for (Bed t: b) out.println(t.getName());
		out.close();
	}

	private String fetchVcfHeader(ArrayList<TempusV3JsonSummary> al){
		StringBuilder sb = new StringBuilder();
		ArrayList<String> jsonFilePaths = new ArrayList<String>();
		LinkedHashMap<String,String> meta = new LinkedHashMap<String, String>();

		for (TempusV3JsonSummary sum: al) {
			//depreciated?
			if (sum.isDepreciated()) continue;
			//json file path
			jsonFilePaths.add(sum.getTempusReport().getJsonFile().getName());
			//meta data
			sum.getTempusReport().addMetaData(meta); 
			sum.getTempusPatient().addMetaData(meta); 
			sum.getTempusOrder().addMetaData(meta); 
			TempusV3GenomicVariants varRes = sum.getTempusV3Results();
			//any variant results? If so add in the sample meta data
			String tmb = varRes.getTumorMutationalBurden();
			int numVars = varRes.getVariants().size();
			if (tmb != null || numVars >=0) {
				TempusV3Specimen[] spec = sum.getTempusSpecimens();
				for (int i=0; i<spec.length; i++) spec[i].addMetaData(meta, i); 
			}
			sum.getTempusV3Results().addMetaData(meta); //check
		}

		//build header
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("##tempusParseCmd=\""+source+"\"\n");
		sb.append("##tempusJsonFiles=\""+Misc.stringArrayListToString(jsonFilePaths, "; ")+"\"\n");
		StringBuilder vcfNames = new StringBuilder (vcfFiles[0].getName());
		for (int i=1; i< vcfFiles.length; i++) {
			vcfNames.append("; ");
			vcfNames.append(vcfFiles[i].getName());
		}
		sb.append("##tempusVcfFiles=\""+vcfNames+"\"\n");
		sb.append("##tempusParseDate=\""+Misc.getDateNoSpaces()+"\"\n");
	
		for (String key: meta.keySet()){
			String value = meta.get(key);
			if (value == null) continue;
			value = value.trim();
			sb.append("##");
			sb.append(key);
			sb.append("=\"");
			sb.append(value);
			sb.append("\"\n");
		}
		
		//add info
		sb.append(altCNV+"\n");
		sb.append(infoEG+"\n");
		sb.append(infoCL+"\n");
		sb.append(infoFE+"\n");
		sb.append(infoFG+"\n");
		sb.append(infoAF+"\n");
		sb.append(infoIMPRECISE+"\n");
		sb.append(infoSVTYPE+"\n");
		sb.append(infoSVLEN+"\n");
		sb.append(infoEND+"\n");
		sb.append(infoDESC+"\n");
		sb.append(infoMATEID+"\n");

		//chrom line
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		return sb.toString();
	}
	 
	private boolean convert() {
		try {
			
			String jString = IO.loadFile(workingJsonFile, " ", true);
	        workingJsonObject = new JSONObject(jString);

	        //schema
	        if (checkSchema(workingJsonObject)== false) return false;
	        
	        reportAttributes = new LinkedHashMap<String,String>();
			addMDPID();
	        reportAttributes.put("jsonFile", workingJsonFile.getCanonicalPath());
	        
	        //report
	        workingReport = new TempusV3Report(workingJsonFile, workingJsonObject, this);
	        workingReport.addAttributes(reportAttributes);
	        
	        //patient
	        workingPatient = new TempusV3Patient(workingJsonObject, this, true);
	        workingPatient.addAttributes(reportAttributes, includePHI);
	        
	        //order
	        workingOrder = new TempusV3Order(workingJsonObject, this);
	        workingOrder.addAttributes(reportAttributes);

	        //specimens, should be a maximum of two!
	        workingSpecimens = TempusV3Specimen.getSpecimens(workingJsonObject, this);
	        TempusV3Specimen.addAttributes(reportAttributes, workingSpecimens);
	        
	        //load any somatic and germline vcf files, might be null
	        workingSomVcfLines = loadVcfLines (workingOrder.getAccessionId(), true);
	        workingGermVcfLines = loadVcfLines (workingOrder.getAccessionId(), false);
	        if (workingSomVcfLines == null && workingGermVcfLines == null) IO.pl("\tWARNING: no vcf lines for "+workingOrder.getAccessionId());
	        else IO.pl();
	        //results
	        workingResults = new TempusV3GenomicVariants(workingJsonObject, this);
	        workingResults.addAttributes(reportAttributes);

	        //check variant ref bps
	        if (indexedFasta != null) checkRefSeqs();
	        return true;
	        
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing this json file "+workingJsonFile);
			return false;
		}
	}	
	
	public static TempusV3JsonSummary parseJsonNoVariants(File v3JsonFile) throws Exception {
		String jString = IO.loadFile(v3JsonFile, " ", true);
        JSONObject object = new JSONObject(jString);
        
        //schema
        if (checkSchema(object)== false) return null;
        
        //parse components
        TempusV3Report report = new TempusV3Report(v3JsonFile, object, null);
        //check status, skip quality failing reports, "biomarker_qns" and "qns"
        if (report.getReportStatus().contains("qns")) return null;
        
        TempusV3Patient patient = new TempusV3Patient(object, null, true);
        TempusV3Order order = new TempusV3Order(object, null);
        TempusV3Specimen[] specimens = TempusV3Specimen.getSpecimens(object, null);
        TempusV3GenomicVariants results = null;

        return new TempusV3JsonSummary(object, order, patient, report, specimens, results);
	}
	
	private String[] loadVcfLines(String accessionId, boolean loadSomaticVcfs) throws IOException {
		//TL-25-1CGRW7DEBN_20250227.soma.freebayes.vcf
		//TL-25-1CGRW7DEBN_20250227.germ.pindel.vcf
		ArrayList<File> vcfsToLoad = new ArrayList<File>();
		for (File v: vcfFiles) {
			String fileName = v.getName();
			if (fileName.contains(accessionId)) {
				if (loadSomaticVcfs == true && fileName.contains(".soma")) vcfsToLoad.add(v);
				else if (loadSomaticVcfs == false && fileName.contains(".germ")) vcfsToLoad.add(v);
			}
		}
		if (vcfsToLoad.size()==0) return null;
		return parseVcfs(vcfsToLoad);
	}

	private String[] parseVcfs(ArrayList<File> vcfsToLoad) throws IOException {
		ArrayList<String> vcfLines = new ArrayList<String>();
		for (File f: vcfsToLoad) {
			BufferedReader in = IO.fetchBufferedReader(f);
			String line = null;
			while ((line = in.readLine())!=null) {
				if (line.startsWith("#") == false) vcfLines.add(line);
			}
		}
		if (vcfLines.size()==0) return null;
		String[] lines = new String[vcfLines.size()];
		vcfLines.toArray(lines);
		return lines;
	}

	private static boolean checkSchema(JSONObject object) {
		if (object.has("metadata") == false) return false;
        JSONObject meta = object.getJSONObject("metadata");
        String jsonSchema = Json.getStringAttribute(meta, "schemaVersion");
        if (jsonSchema == null) Misc.printErrAndExit("\nschemaVersion not found! Aborting. ");
        if (jsonSchema.startsWith(acceptedSchema)) return true;
        return false;
	}

	private void checkRefSeqs() throws Exception {
		for (TempusV3Variant tv: workingResults.getVariants()) {
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
			IO.pl(numPotentiallyActionable+  			"\tNum Potentially Actionable");
			IO.pl(numBiologicallyRelevant+   			"\tNum Biologically Relevant");
			IO.pl(numPathogenic+   			"\tNum Pathogenic");
			IO.pl(numLikelyPathogenic+   			"\tNum Likely Pathogenic");
			IO.pl(numRiskAllele+  	"\tNum Risk Allele");
			IO.pl(numUnknownSignificance+  	"\tNum Unknown Significance");
			IO.pl(failedToFindCooridinates.size()+						 "\tNum Variants skipped for failing to match genomic coordinates.");
		

		IO.pl("\nBioInf Pipelines");
		Misc.printTreeMap(bioInfoPipeline, "\t", "\t");
		IO.pl("\nPhysicians:");
		Misc.printTreeMap(physicians, "\t", "\t");
		IO.pl("\nReport Status (qns is failed):");
		Misc.printTreeMap(reportStatus, "\t", "\t");
		IO.pl("\nTest Codes:");
		Misc.printTreeMap(testCodes, "\t", "\t");
		IO.pl("\nTest Descriptions:");
		Misc.printTreeMap(testDescriptions, "\t", "\t");
		IO.pl("\nOrigin Path Lab Diagnosis:");
		Misc.printTreeMap(diagnosis, "\t", "\t");
		IO.pl("\nSample Categories:");
		Misc.printTreeMap(sampleCategories, "\t", "\t");
		IO.pl("\nSample Sites:");
		Misc.printTreeMap(sampleSites, "\t", "\t");
		IO.pl("\nPercent Tumor in Sample");
		tumorPercentages.printScaledHistogram();
		
			IO.pl("\nGene Mutations:");
			Misc.printTreeMap(genes, "\t", "\t");
			IO.pl("\nGenomic Source Class:");
			Misc.printTreeMap(genomicSourceClass, "\t", "\t");
			IO.pl("\nVariant Type:");
			Misc.printTreeMap(variantType, "\t", "\t");
			IO.pl("\nVariant Description:");
			Misc.printTreeMap(variantDescription, "\t", "\t");
			IO.pl("\nTumor Variant Allele Percentages:");
			tumorAF.printScaledHistogram();
			//IO.pl("\nAge at Diagnosis");		Not used?
			//ageAtDiagnosis.printScaledHistogram();
		
		if (failedToFindCooridinates.size()!=0) {
			IO.pl("\n"+failedToFindCooridinates.size()+"\tTempus variants could not be matched with genomic coordinates:");
			for (TempusV3Variant tv: failedToFindCooridinates) IO.pl(tv.toString());
		}

	}
	
	public void printStatsMinimal(){
		IO.pl("\tTempus BioInf Pipelines: "+bioInfoPipeline);
		IO.pl("\tTempus Physicians: "+physicians);
		IO.pl("\tTempus Report Status (qns is failed): "+reportStatus);
		IO.pl("\tTempus Test Codes: "+testCodes);
		IO.pl("\tTempus Test Descriptions: "+testDescriptions);
		IO.pl("\tTempus Diagnosis: "+diagnosis);
		IO.pl("\tTempus Sample Categories: "+sampleCategories);
		IO.pl("\tTempus Sample Sites: "+sampleSites);
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
		new TempusV3Json2Vcf(args);
	}		

	/**This method will process each argument and assign new varibles
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		source = useqVersion+" Args: USeq/TempusV3Json2Vcf "+ Misc.stringArrayToString(args, " ");
		System.out.println("\n"+ source +"\n");
		File forExtraction = null;
		String attributes = null;
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
					case 'g': geneAliasFile = new File(args[++i]); break;
					case 'b': bed = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'i': includePHI = true; break;
					case 'a': attributes = args[++i]; break;
					
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

		//pull vcf files
		File[][] vcfs = new File[3][];
		vcfs[0] = IO.extractFiles(forExtraction, ".vcf");
		vcfs[1] = IO.extractFiles(forExtraction,".vcf.gz");
		vcfs[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(vcfs);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your Tempus xxx.vcf(.zip/.gz OK) file(s)!\n");

		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save directory!\n"+saveDirectory);
		saveDirectory.mkdirs();
		if (saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: your save directory does not appear to be a directory?\n");
		
		//load hash of gene name and Bed to use in coordinates for CNVs
		if (bed != null) {
			if (indexedFasta == null) Misc.printErrAndExit("\nError: addition of CNV info requires an indexed fasta\n");
			cnvGeneNameBed = Bed.parseBedForNames(bed);
		}
		
		//Load the gene aliases
		if (geneAliasFile == null || geneAliasFile.canRead()== false) Misc.printErrAndExit("\nError: failed to find your gene aliases file : "+geneAliasFile);
		loadGeneAliases();
		
		if (attributes !=null) {
			String[] at = Misc.COMMA.split(attributes);
			keysToExport = new LinkedHashSet<String>();
			for (String s: at) keysToExport.add(s);
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Tempus V3 Json 2 Vcf: May 2025                        **\n" +
				"**************************************************************************************\n" +
				"This app groups json v3+ Tempus reports by patient and tumor. It then merges the \n"+
				"variant, IHC, MMR, RNA results and generates a combine vcf with both germline and \n"+
				"somatic snvs, indels, cnvs, and gene fusions.  Clinical attributes are included in the\n"+
				"vcf header as well as written out to a spreadsheet. Vcfs from tumor only submissions\n"+
				"may contain a mix of somatic and inherited variants. Vt normalize the exported vcfs,\n"+
				"https://github.com/atks/vt and CrossMap https://github.com/liguowang/CrossMap the\n"+
				"coordinates to hg38.\n"+

				"\nOptions:\n"+
				"-j Path to a directory containing json and vcf files downloaded from Tempus for one\n"+
				"     or more patients, (.json/.vcf)(.gz/.zip OK)\n"+
				"-s Path to a directory for saving the parsed results.\n"+
				"-b Path to a bed file for converting CNV and fusion gene names to coordinates where the\n"+
				"     bed name column contains just the gene name.\n"+
				"-f Path to the reference hg19 fasta with xxx.fai index for gene conversions.\n"+
				"-g Path to a gene alias file for genomic coordiante conversions.\n"+
				"-i Include PHI in spreadsheet output, defaults to excluding it. \n"+
				"-a Print this list of attributes in spreadsheet, comma delimited, case sensitive, no\n"+
				"     spaces. Defaults to all. Run without to get a list of available attributes.\n"+
				
				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/TempusV3Json2Vcf -j /F1/TempusJsonsVcfs\n" +
				"     -f /Ref/human_g1k_v37.fasta -s /F1/ParsedResults/ -b /Ref/b37TempusGeneRegions.bed\n"+
				"     -a 'accessionId,tumorMutationalBurden,msiStatus,diagnosis' -g /Ref/geneAli.txt.gz\n\n" +

				"**************************************************************************************\n");

	}

	public Histogram getAgeAtDiagnosis() {
		return ageAtDiagnosis;
	}
	public TempusV3Report getWorkingReport() {
		return workingReport;
	}
	public TempusV3Patient getWorkingPatient() {
		return workingPatient;
	}
	public TempusV3Order getWorkingOrder() {
		return workingOrder;
	}
	public TempusV3Specimen[] getWorkingSpecimens() {
		return workingSpecimens;
	}
	public TempusV3GenomicVariants getWorkingResults() {
		return workingResults;
	}
	public HashMap<String, Bed> getCnvGeneNameBed() {
		return cnvGeneNameBed;
	}
	public IndexedFastaSequenceFile getFasta() {
		return fasta;
	}
	public String getJsonSchema() {
		return jsonSchema;
	}
	public ArrayList<TempusV3JsonSummary> getSummaries() {
		return summaries;
	}
	public String[] getWorkingSomVcfLines() {
		return workingSomVcfLines;
	}
	public String[] getWorkingGermVcfLines() {
		return workingGermVcfLines;
	}
	public HashMap<String, String[]> getGeneAliases() {
		return geneAliases;
	}
	public ArrayList<TempusV3Variant> getFailedToFindCooridinates() {
		return failedToFindCooridinates;
	}

	public void setWorkingNumPotentiallyActionable(int workingNumPotentiallyActionable) {
		this.workingNumPotentiallyActionable = workingNumPotentiallyActionable;
	}

	public void setWorkingNumBiologicallyRelevant(int workingNumBiologicallyRelevant) {
		this.workingNumBiologicallyRelevant = workingNumBiologicallyRelevant;
	}

	public void setWorkingNumLikelyPathogenic(int workingNumLikelyPathogenic) {
		this.workingNumLikelyPathogenic = workingNumLikelyPathogenic;
	}
	
	public void setWorkingNumPathogenic(int workingNumPathogenic) {
		this.workingNumPathogenic = workingNumPathogenic;
	}

	public void setWorkingNumRiskAllele(int workingNumRiskAllele) {
		this.workingNumRiskAllele = workingNumRiskAllele;
	}

	public void setWorkingNumUnknownSignificance(int workingNumUnknownSignificance) {
		this.workingNumUnknownSignificance = workingNumUnknownSignificance;
	}

	public File getWorkingJsonFile() {
		return workingJsonFile;
	}
	
}
