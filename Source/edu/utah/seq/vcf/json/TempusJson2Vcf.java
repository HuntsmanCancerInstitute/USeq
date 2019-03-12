package edu.utah.seq.vcf.json;

import java.io.*;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.regex.*;
import org.json.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import java.util.concurrent.TimeUnit;


/**
 * Takes one or more patient xml reports from FoundationOne tests and converts most of the variants into vcf format.
 * 
 * @author david.nix@hci.utah.edu 
 **/
public class TempusJson2Vcf {

	//user defined fields
	private File[] jsonFiles;
	private File indexedFasta = null;
	private File saveDirectory;
	private String acceptedSchema = "1.3";
	
	//internal fields
	private IndexedFastaSequenceFile fasta; 
	private long numSnvs = 0;
	private long numIndels = 0;
	private long numCopyNumber = 0;
	private long numRearrangements = 0;
	private String source = null;
	
	//counters
	TreeMap<String, Integer> bioInfPipelines = new TreeMap<String, Integer>();
	TreeMap<String, Integer> physicians = new TreeMap<String, Integer>();
	TreeMap<String, Integer> testCodes = new TreeMap<String, Integer>();
	TreeMap<String, Integer> testDescriptions = new TreeMap<String, Integer>();
	TreeMap<String, Integer> diagnosis = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleCategories = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleSites = new TreeMap<String, Integer>();
	TreeMap<String, Integer> sampleTypes = new TreeMap<String, Integer>();
	TreeMap<String, Integer> somaticGenes = new TreeMap<String, Integer>();
	Histogram tumorPercentages = new Histogram(0,100, 20);
	Histogram tumorAF = new Histogram(0,100, 20);
	Histogram tumorDP = new Histogram(0,5000, 50);
	
	//working data for a particular report
	private File workingJsonFile;

	//constructors
	public TempusJson2Vcf(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();
			
			//stats
			/*System.out.println("\nParsing stats for "+jsonFiles.length+" files:");
			System.out.println(numSnvs+ "\t# Short Pass");
			System.out.println(numIndels+ "\t# Short Fail");
			System.out.println(numCopyNumber+ "\t# CNV Pass");
			System.out.println(numRearrangements+ "\t# CNV Fail");*/
			
			//sampleInfo
			//printSampleInfo();

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
		//fasta = new IndexedFastaSequenceFile(indexedFasta);
		//if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);

		System.out.println("Parsing and coverting...\n");
		for (int i=0; i< jsonFiles.length; i++){
			
			//process file
			workingJsonFile = jsonFiles[i];
			System.out.println(workingJsonFile.getName());
			convert();
			
			//build vcf
			//writeVcf();
			
			//create an SI
			/*
			SampleInfo si = new SampleInfo(reportAttributes);
			if (si.allLoaded() == false) {
				System.err.println("WARNING: missing sample info for:");
				System.err.println(si.toString());
				problemParsing = true;
			}
			sampleInfo.add(si);
			
			if (problemParsing) failingFiles.add(workingXmlFile.getName());*/
		}
		
		//close the fasta lookup fetcher
		//fasta.close();
	
		printStats();

	}

	private void writeVcf() {
		/*
		try {
			String name = Misc.removeExtension(workingJsonFile.getName());
			File vcf = new File (saveDirectory, name+".vcf");
			PrintWriter out = new PrintWriter( new FileWriter(vcf));
			
			//add header
			out.print(buildVcfHeader());
			
			//Create a bed record for each and sort
			Bed[] b = new Bed[shortVariants.size() + copyVariants.size()+ (2*rearrangeVariants.size())];
			int counter =0;
			
			for (FoundationShortVariant sv : shortVariants) b[counter++] = fetchBed(sv.toVcf(counter));
			for (FoundationCopyVariant cnv : copyVariants) b[counter++] = fetchBed(cnv.toVcf(counter));
			for (FoundationRearrangeVariant r : rearrangeVariants) {
				String[] v = r.toVcf(counter);
				 b[counter++] = fetchBed(v[0]);
				 b[counter++] = fetchBed(v[1]);
			}
			Arrays.sort(b);
			for (Bed t: b) out.println(t.getName());
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: issue writing out vcf for "+workingXmlFile);
		}*/

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
		
		//add in meta info
		/*
		for (String key: reportAttributes.keySet()){
			String value = reportAttributes.get(key).trim();
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
		}*/
		//add FILTER?
		//if (problemParsing) sb.append("##FILTER=<ID=ci,Description=\"Converting from xml to vcf issue, treat skeptically.\"\n");
		
		//add ALT and INFO lines
		TreeSet<String> alt = new TreeSet<String>();
		TreeSet<String> info = new TreeSet<String>();
		
		//if (shortVariants.size()!=0) FoundationShortVariant.appendInfoLines(info);
		//if (copyVariants.size()!=0) FoundationCopyVariant.appendInfoAltLines(alt, info);
		//if (rearrangeVariants.size()!=0) FoundationRearrangeVariant.appendInfoAltLines(alt, info);
		
		if (alt.size()!=0) {
			sb.append(Misc.stringArrayToString(Misc.setToStringArray(alt), "\n")); 
			sb.append("\n");
		}
		if (info.size()!=0) {
			sb.append(Misc.stringArrayToString(Misc.setToStringArray(info), "\n")); 
			sb.append("\n");
		}
		
		//chrom line
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
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
	        TempusReport report = new TempusReport(object, this);
	        
	        //patient
	        TempusPatient patient = new TempusPatient(object, this);
	        
	        //order
	        TempusOrder order = new TempusOrder(object, this);

	        
	        //specimens, should be just two!
	        TempusSpecimen[] specimens = TempusSpecimen.getSpecimens(object, this);
	        
	        //results
	        TempusResults results = new TempusResults(object, this); 
	        
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing this json file "+workingJsonFile);
		}
	}	
	
	public void printStats(){
		IO.pl("\nSummary Stats:\n");
		IO.pl("BioInfPipelines");
		Misc.printTreeMap(bioInfPipelines, "\t", "\t");
		IO.pl("\nPhysicians:");
		Misc.printTreeMap(physicians, "\t", "\t");
		IO.pl("\nTestCodes:");
		Misc.printTreeMap(testCodes, "\t", "\t");
		IO.pl("\nTestDescriptions:");
		Misc.printTreeMap(testDescriptions, "\t", "\t");
		IO.pl("\nDiagnosis:");
		Misc.printTreeMap(diagnosis, "\t", "\t");
		IO.pl("\nSampleCategories:");
		Misc.printTreeMap(sampleCategories, "\t", "\t");
		IO.pl("\nSampleSites:");
		Misc.printTreeMap(sampleSites, "\t", "\t");
		IO.pl("\nSampleTypes:");
		Misc.printTreeMap(sampleTypes, "\t", "\t");
		IO.pl("\nSomaticGeneMutations:");
		Misc.printTreeMap(somaticGenes, "\t", "\t");
		
		IO.pl("\nPercent Tumor in Sample");
		tumorPercentages.printScaledHistogram();
		IO.pl("\nTumor Variant Allele Percentages");
		tumorAF.printScaledHistogram();
		IO.pl("\nTumor Variant Read Depth");
		tumorDP.printScaledHistogram();
		
		
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
		source = useqVersion+" Args: "+ Misc.stringArrayToString(args, " ");
		System.out.println("\n"+ source +"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'j': forExtraction = new File(args[++i]); break;
					case 'f': indexedFasta = new File(args[++i]); break;
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


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Tempus Json 2 Vcf: March 2019                        **\n" +
				"**************************************************************************************\n" +
				"Attempts to parse json Tempus reports to vcf. Consider removing PHI elements: \n"+
				"grep -vwE '(firstName|lastName|emr_id|DoB)' TL18.json > TL18.clean.json\n"+

				"\nOptions:\n"+
				"-j Path to Tempus json reports or directory containing such.\n"+
				"-s Path to a directory for saving the results.\n"+
				"-f Path to the reference fasta with xxx.fai index\n"+
				"-o Skip variants that clearly fail to convert, e.g. var seq doesn't match fasta.\n"+
				"     Defaults to marking 'ci' in FILTER field.\n"+
				
				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/TempusJson2Vcf -j /F1/TemJsons\n" +
				"     -f /Ref/human_g1k_v37.fasta -s /F1/VCF/ \n\n" +

				"**************************************************************************************\n");

	}
	
	public IndexedFastaSequenceFile getFasta() {
		return fasta;
	}
	
}
