package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;
import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.UsernamePasswordCredentials;
import org.apache.commons.httpclient.auth.AuthScope;
import org.apache.commons.httpclient.methods.GetMethod;
import org.json.JSONArray;
import org.json.JSONObject;
import edu.utah.hci.query.LocalQuery;
import edu.utah.hci.query.UserQuery;

/**Calculates call frequencies for each vcf record from a db of vcf files and callable region files. 
 * Most mutated single site in p53 is seen at 0.039 across all cancers in TCGA.  0.57 for BRAF V600E across all cancers.*/
public class VCFCallFrequency {

	private File[] vcfFiles;
	private File saveDirectory;
	private String queryURL = null;
	private String host = null;
	private String realm = null;
	private String queryKey = null;
	private String userName = null; 
	private String password = null; 
	private String vcfFileFilter = null;
	private String bedFileFilter = null;
	private int numberRecordsPerQuery = 1000;
	private String vcfSearchUrl = null;
	private String bedSearchUrl = null;
	private String fetchOptionsUrl = null;
	private double[] idBedCount = null;
	private double[] idVcfCount;
	private double maxCallFreq = 1.0;
	private int minBedCount = 8;
	private boolean appendFilter = true;
	//local searching
	private File dataDir = null;
	private File indexDir = null;
	private LocalQuery localQuery = null;
	private UserQuery userQueryVcf = null;
	private UserQuery userQueryBed = null;

	private boolean debug = false;
	private HttpClient client = new HttpClient();
	private Histogram hist = new Histogram(0,1,20);
	private int totalQueries = 0;
	private int totalQueriesWithHits = 0;
	private int numFailingCallFreq = 0;
	private int numPassingCallFreq = 0;
	private int numWithLowBedCount = 0;


	public VCFCallFrequency (String[] args) {
		
		
		
		long startTime = System.currentTimeMillis();
		processArgs(args);

		if (queryURL != null) buildSearchUrls(); 
		else localQuery = new LocalQuery(dataDir, indexDir, false);

		for (File vcf: vcfFiles) {
			IO.p(vcf.getName());
			annotate(vcf);
		}
		//print stats
		IO.pl("\nTotal Vcf Queries: "+ totalQueries +", with matches: "+totalQueriesWithHits+", passing callFreq: "+numPassingCallFreq+", failing callFreq: "+numFailingCallFreq+ ", too low bed count to score: "+numWithLowBedCount);
		IO.pl("\nHistogram of call frequencies meeting the minimum bed count:");
		hist.printScaledHistogram();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");

	}

	private void buildSearchUrls() {
		//fetch key
		fetchEncodedQueryKey();

		//set search url
		vcfSearchUrl = queryURL+"search?key="+ queryKey+"&regExAll=.vcf.gz;"+vcfFileFilter+"&matchVcf=true&vcf=";
		bedSearchUrl = queryURL+"search?key="+ queryKey+"&regExAll=.bed.gz;"+bedFileFilter+"&bed=";
		fetchOptionsUrl = queryURL+"search?key="+ queryKey+"&fetchOptions=true";
		IO.pl("fetchOptions:\t"+fetchOptionsUrl);
		IO.pl("\tCall this in a browser to pick an appropriate file filter for your analysis. Only those data files with the file filter in their path will be queried.");
		IO.pl("VcfSearchUrl:\t"+vcfSearchUrl);
		IO.pl("BedSearchUrl:\t"+bedSearchUrl);
		IO.pl();
	}

	public void fetchEncodedQueryKey() {
		try {
			String urlStr = queryURL+"fetchKey";
			GetMethod getMethod = new GetMethod(urlStr);
			UsernamePasswordCredentials upc = new UsernamePasswordCredentials(userName, password);
			AuthScope as = new AuthScope(host, 8080, realm);
			client.getState().setCredentials(as, upc);
			int status = client.executeMethod(getMethod);
			String responseBody = getMethod.getResponseBodyAsString();
			getMethod.releaseConnection();

			if (status == 200) queryKey= encodeString(responseBody);
			else throw new Exception("\nFaild to fetch an authentication key? \n"+responseBody);

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private String encodeString(String value) throws UnsupportedEncodingException {
		return URLEncoder.encode(value, StandardCharsets.UTF_8.toString());
	}

	public void annotate(File vcfFile){

		File modVcf = new File (saveDirectory, Misc.removeExtension(vcfFile.getName())+".callFreq.vcf.gz");
		String line = null;
		Gzipper out = null;
		BufferedReader in = null;
		boolean addInfo = true;
		boolean addFilter = true;
		if (maxCallFreq == 1.0) addFilter = false;
		try {
			out = new Gzipper(modVcf);
			in = IO.fetchBufferedReader(vcfFile);
			int counter = 0;
			String[][] vcfRecords = new String[numberRecordsPerQuery][];
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0)continue;
				if (line.startsWith("#")) {
					out.println(line);
					if (addInfo && line.startsWith("##INFO")){
						out.println("##INFO=<ID=CF,Number=1,Type=Float,Description=\"Frequency the variant has been observed in '"+vcfFileFilter
								+ "' variant datasets, number vcf matches, number bed matches.\">");
						addInfo = false;
					}
					else if (addFilter && line.startsWith("##FILTER")){
						out.println("##FILTER=<ID=CallFreq,Description=\"Variant exceeded a call frequency of "+maxCallFreq+
								" from querying vcf and callable region bed files in '"+vcfFileFilter+"'\">");
						addFilter = false;
					}
				}
				else {
					vcfRecords[counter] = Misc.TAB.split(line);
					totalQueries++;
					//check length
					if (vcfRecords[counter].length < 8) throw new Exception ("This vcf record is missing 8 manditory fields?! "+line);
					counter++;
					if (counter == numberRecordsPerQuery) {
						processRecordBlock(vcfRecords, out);
						counter = 0;
						vcfRecords = new String[numberRecordsPerQuery][];
					}
				}
			}
			//process last
			processRecordBlock(vcfRecords, out);
			IO.pl();

			out.close();
			in.close();

		} catch (Exception e) {
			out.closeNoException();
			IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing this record "+line);
		} 
	}

	private void processRecordBlock(String[][] vcfRecords, Gzipper out) throws Exception {

		//run queries
		JSONObject jsonResultsVcf = null;
		JSONObject jsonResultsBed = null;

		//remote service?
		if (queryURL != null){
			//create search strings, any vcfs to query?
			String[] searchStrings = makeSearchUrls(vcfRecords);
			if (searchStrings != null) {
				//execute vcf query
				jsonResultsVcf = new JSONObject(query(searchStrings[0]));
				jsonResultsBed =  new JSONObject(query(searchStrings[1]));
			}
		}
		//local file system query?
		else {
			if (loadVcfsBedsToQuery(vcfRecords)) {
				localQuery.query(userQueryVcf);
				localQuery.query(userQueryBed);
				jsonResultsVcf  = userQueryVcf.getResults();
				jsonResultsBed = userQueryBed.getResults();
			}
		}

		//parse and annotate
		if (jsonResultsVcf != null && jsonResultsBed != null) {

			if (debug) IO.p("VCF RETURN\n"+jsonResultsVcf.toString(1));			
			loadVcfResults(jsonResultsVcf);

			if (debug) IO.p("BED RETURN\n"+jsonResultsBed.toString(1));
			loadBedResults(jsonResultsBed);

			//for each record
			for (int i=0; i< vcfRecords.length; i++){
				//pull counts and calc ratio
				String[] vcf = vcfRecords[i];
				if (vcf == null) break;
				double vcfCount = idVcfCount[i];
				double bedCount = idBedCount[i];
				double ratio = 0;
				if (bedCount > 0) {
					ratio = vcfCount/bedCount;
					if (ratio > 1.0) ratio = 1.0;
				}
				String cf = "CF="+Num.formatNumber(ratio, 3)+","+(int)vcfCount+","+(int)bedCount;

				//modify INFO field of each record
				vcfRecords[i][7] = cf+";"+vcfRecords[i][7];
				boolean printMe = true;

				//check bed count, must pass this before any filtering is applied.
				if (bedCount >= minBedCount) {
					hist.count(ratio);
					if (ratio > maxCallFreq){
						//just append or skip
						if (appendFilter) vcfRecords[i][6] = modifyFilter(vcfRecords[i][6]);
						else printMe = false;
						numFailingCallFreq++;
					}
					else numPassingCallFreq++;
				}
				else numWithLowBedCount++;

				//print it?
				if (printMe){
					out.print(vcf[0]);
					for (int x=1; x< vcf.length; x++){
						out.print("\t");
						out.print(vcf[x]);
					}
					out.println();
				}
			}
		}
		//else throw new IOException("Error! Failed to generate search urls?! Aborting.");
		System.out.print(".");

	}

	private String modifyFilter(String filter) {
		if (filter.equals(".") || filter.toUpperCase().equals("PASS")) return "CallFreq";
		return "CallFreq;"+filter;
	}

	/**Loads the number of hits returned to each vcf query. Some won't be returned and thus will be zero.*/
	private void loadVcfResults(JSONObject results) {
		//clear prior, all set to zero
		idVcfCount = new double[numberRecordsPerQuery];

		//any results?
		if (results.has("queryResults")){
			JSONArray ja = results.getJSONArray("queryResults");

			int numWithHits = ja.length();
			totalQueriesWithHits+= numWithHits;

			for (int i=0; i< numWithHits; i++){
				JSONObject jo = ja.getJSONObject(i);
				double numHits = jo.getDouble("numberHits");
				String vcfRecord = jo.getString("input");
				String[] splitVcf = Misc.TAB.split(vcfRecord);
				int index = Integer.parseInt(splitVcf[2]);
				//set in array
				idVcfCount[index] = numHits;
			}
		}
	}


	/**Loads the number of hits returned to each bed query. Some may be returned and thus will be zero.*/
	private void loadBedResults(JSONObject results) {
		//clear prior, set to zero
		idBedCount = new double[numberRecordsPerQuery];

		//any results?
		if (results.has("queryResults")){
			JSONArray ja = results.getJSONArray("queryResults");
			int numWithHits = ja.length();
			for (int i=0; i< numWithHits; i++){
				JSONObject jo = ja.getJSONObject(i);
				double numHits = jo.getDouble("numberHits");
				String bedRecord = jo.getString("input");
				String[] splitbed = Misc.TAB.split(bedRecord);
				int index = Integer.parseInt(splitbed[3]);
				//set in array
				idBedCount[index] = numHits;
			}
		}
	}

	private String query(String searchString) {
		try {
			GetMethod getMethod = new GetMethod(searchString);
			int status = client.executeMethod(getMethod);
			String responseBody = getMethod.getResponseBodyAsString();
			getMethod.releaseConnection();
			if (status == 200) return responseBody;
			else throw new Exception("\nFaild to execute the search?\n"+searchString+"\n"+responseBody);

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		return null;
	}

	private String[] makeSearchUrls(String[][] vcfRecords) {
		StringBuilder[] vcfBed = pullVcfsBedsToQuery(vcfRecords);
		if (vcfBed[0].length() == 0) return null;
		if (debug) IO.pl("VcfBedToQuery:\n\t"+vcfBed[0]+"\n\t"+vcfBed[1]);

		//make cmd with key
		String urlVcfStr = vcfBed[0].insert(0, vcfSearchUrl).toString();
		String urlBedStr = vcfBed[1].insert(0, bedSearchUrl).toString();
		if (debug) IO.pl("SearchURL:\n\t"+urlVcfStr+"\n\t"+urlBedStr);
		return new String[]{urlVcfStr,urlBedStr};
	}

	private StringBuilder[] pullVcfsBedsToQuery(String[][] vcfRecords) {
		//pull vcfs to query
		StringBuilder sbVcf = new StringBuilder();
		StringBuilder sbBed = new StringBuilder();
		for (int i=0; i< vcfRecords.length; i++){
			String[] vcf = vcfRecords[i];
			if (vcf == null) break;
			//build vcf query
			//pull chr, pos, ref, alt
			sbVcf.append(vcf[0]); sbVcf.append("%09");
			//pos
			sbVcf.append(vcf[1]); sbVcf.append("%09");	
			//replace id with index number of vcfRecords[index][]
			sbVcf.append(i); sbVcf.append("%09");
			//ref
			sbVcf.append(vcf[3]); sbVcf.append("%09");
			//alt
			sbVcf.append(vcf[4]);
			//remainder
			sbVcf.append("%09.%09.%09.;");

			//build bed region query
			//chrom
			sbBed.append(vcf[0]); sbBed.append("%09");
			//start
			int zeroPos = Integer.parseInt(vcf[1])-1;
			sbBed.append(zeroPos); sbBed.append("%09");
			//max stop
			sbBed.append(zeroPos+maxSize(vcf[4])); sbBed.append("%09");
			//name, score, strand
			sbBed.append(i); sbBed.append("%090%09.;");
		}
		return new StringBuilder[]{sbVcf, sbBed};
	}

	private boolean loadVcfsBedsToQuery(String[][] vcfRecords) {
		boolean loaded = false;
		userQueryVcf.clearResultsRegionsRecords();
		userQueryBed.clearResultsRegionsRecords();

		for (int i=0; i< vcfRecords.length; i++){
			String[] vcf = vcfRecords[i];
			if (vcf == null) break;
			StringBuilder sbVcf = new StringBuilder();
			StringBuilder sbBed = new StringBuilder();
			//build vcf query
			//pull chr, pos, ref, alt
			sbVcf.append(vcf[0]); sbVcf.append("\t");
			//pos
			sbVcf.append(vcf[1]); sbVcf.append("\t");	
			//replace id with index number of vcfRecords[index][]
			sbVcf.append(i); sbVcf.append("\t");
			//ref
			sbVcf.append(vcf[3]); sbVcf.append("\t");
			//alt
			sbVcf.append(vcf[4]);
			//remainder
			sbVcf.append("\t.\t.\t.");
			userQueryVcf.addVcfRecord(sbVcf.toString());

			//build bed region query
			//chrom
			sbBed.append(vcf[0]); sbBed.append("\t");
			//start
			int zeroPos = Integer.parseInt(vcf[1])-1;
			sbBed.append(zeroPos); sbBed.append("\t");
			//max stop
			sbBed.append(zeroPos+maxSize(vcf[4])); sbBed.append("\t");
			//name, score, strand
			sbBed.append(i); sbBed.append("\t0\t.");
			userQueryBed.addBedRegion(sbBed.toString());
			loaded = true;
		}
		return loaded;
	}

	public static int maxSize(String alts){
		String[] a = Misc.COMMA.split(alts);
		int max = a[0].length();
		for (int i=1; i< a.length; i++) {
			if (a[i].length() > max) max = a[i].length();
		}
		return max;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFCallFrequency(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		File configFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[++i]); break;
					case 'v': vcfFileFilter = args[++i]; break;
					case 'b': bedFileFilter = args[++i]; break;
					case 'x': appendFilter = false; break;
					case 'c': configFile = new File(args[++i]); break;
					case 'i': indexDir = new File(args[++i]); break;
					case 'd': dataDir = new File(args[++i]); break;
					case 'm': maxCallFreq = Double.parseDouble(args[++i]); break;
					case 'o': minBedCount = Integer.parseInt(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'e': debug = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//config file? or local
		if (configFile != null){
			if (configFile.canRead() == false) Misc.printErrAndExit("\nSorry, cannot read your config file: "+configFile+"\n");
			HashMap<String, String> config = IO.loadFileIntoHashMapLowerCaseKey(configFile);
			if (config.containsKey("queryurl")) queryURL = config.get("queryurl");
			if (config.containsKey("host")) host = config.get("host");
			if (config.containsKey("realm")) realm = config.get("realm");
			if (config.containsKey("username")) userName = config.get("username");
			if (config.containsKey("password")) password = config.get("password");
			if (config.containsKey("maxcallfreq")) maxCallFreq = Double.parseDouble(config.get("maxcallfreq"));
			if (config.containsKey("vcffilefilter")) vcfFileFilter = config.get("vcffilefilter");
			if (config.containsKey("bedfilefilter")) bedFileFilter = config.get("bedfilefilter");
		}
		//local search?
		if (queryURL == null){
			if (dataDir == null || indexDir == null || dataDir.isDirectory()== false || indexDir.isDirectory() == false) {
				Misc.printErrAndExit("\nProvide either a configuration file for remotely accessing a genomic query service or "
						+ "two directory paths to the Data and Index directories uses by the GQueryIndexer app.\n");;
			}
		}

		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments:");
		IO.pl("\t-f Vcfs "+forExtraction);
		IO.pl("\t-s SaveDir "+saveDirectory);
		IO.pl("\t-v Vcf File Filter "+vcfFileFilter);
		IO.pl("\t-b Bed File Filter "+bedFileFilter);
		IO.pl("\t-m MaxCallFreq "+maxCallFreq);
		IO.pl("\t-o MinBedCount "+minBedCount);
		IO.pl("\t-x Remove failing "+(appendFilter==false));
		IO.pl("\t-e Verbose "+debug);
		if (queryURL != null){
			IO.pl("\tQueryUrl "+queryURL);
			IO.pl("\tHost "+host);
			IO.pl("\tRealm "+realm);
			IO.pl("\tUserName "+userName);
			//check config params
			if (queryURL == null) Misc.printErrAndExit("\nError: failed to find a queryUrl in the config file, e.g. queryUrl http://hci-clingen1.hci.utah.edu:8080/Query/");
			if (queryURL.endsWith("/") == false) queryURL = queryURL+"/";
			if (host == null) Misc.printErrAndExit("\nError: failed to find a host in the config file, e.g. host hci-clingen1.hci.utah.edu");
			if (realm == null) Misc.printErrAndExit("\nError: failed to find a realm in the config file, e.g. realm QueryAPI");
			if (userName == null) Misc.printErrAndExit("\nError: failed to find a userName in the config file, e.g. userName FCollins");
			if (password == null) Misc.printErrAndExit("\nError: failed to find a password in the config file, e.g. password g0QueryAP1");

		}
		else {
			IO.pl("\tDataDir "+dataDir);
			IO.pl("\tIndexDir "+indexDir);
		}
		IO.pl();

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		//check params
		if (vcfFileFilter == null) Misc.printErrAndExit("\nError: provide a vcf file filter, e.g. Hg38/Somatic/Avatar/Vcf ");
		if (vcfFileFilter == null) Misc.printErrAndExit("\nError: provide a bed file filter, e.g. Hg38/Somatic/Avatar/Bed ");
		if (saveDirectory == null) Misc.printErrAndExit("\nError: provide a directory to save the annotated vcf files.");
		else saveDirectory.mkdirs();
		if (saveDirectory.exists() == false) Misc.printErrAndExit("\nError: could not find your save directory? "+saveDirectory);

		userQueryVcf = new UserQuery().addRegexAll(".vcf.gz").addRegexAll(vcfFileFilter).matchVcf();
		userQueryBed = new UserQuery().addRegexAll(".bed.gz").addRegexAll(bedFileFilter);
	}	

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                            VCF Call Frequency: Feb 2020                          **\n" +
				"**************************************************************************************\n" +
				"Calculates a vcf call frequency for each variant when pointed at a genomic Query\n"+
				"service (https://github.com/HuntsmanCancerInstitute/GQuery) or the Data and Index\n"+
				"directories the service is accessing. CallFreq's are calculated \n"+
				"by first counting the number of exact vcf matches present and dividing it by the\n"+
				"number of intersecting callable bed files. Use this tool to flag variants with high call rates that\n"+
				"are potential false positives. Some are not, e.g. BRAF V600E occurs in 58% of cancers\n"+
				"with a BRAF mutation. So treate the call freq as an annotation requiring context level\n"+
				"interpretation. Use the file filter to limit which files are included in the call freq\n"+
				"calculations. Only file paths containing the file filter will be included. Be sure to\n"+
				"place both the sample vcf and associated callable region bed files in the same Query\n"+
				"index folder path as defined by -f. \n"+

				"\nRequired Options:\n"+
				"-f Full path to a file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				"-s Directory to save the annotated vcf files\n"+
				"-v GQuery service vcf file filter, e.g. Hg38/Somatic/Avatar/Vcf\n"+
				"-b GQuery service callable region bed file filter, e.g. Hg38/Somatic/Avatar/Bed\n"+
				"-c Config txt file containing two tab delimited columns with host, queryUrl, realm, \n"+
				"     userName, password, and (optionally) fileFilter and or maxCallFreq. 'chmod 600'\n"+
				"     the file! e.g.: \n"+
				"     host hci-clingen1.hci.utah.edu\n"+
				"     queryUrl http://hci-clingen1.hci.utah.edu:8080/GQuery/\n"+
				"     realm GQuery\n"+
				"     userName FColins\n"+
				"     password g0QueryAP1\n"+
				"-i (Alternative to -c), provide a path to the Index directory generated by the USeq \n"+
				"     QueryIndexer app.\n"+
				"-d (Alternative to -c), provide a path to the Data directory used in creating the Index.\n"+

				"\nOptions:\n"+
				"-m Maximum call freq, defaults to 1, before appending 'CallFreq' to the FILTER field.\n"+
				"-x Remove failing max call freq records, not recommended.\n"+
				"-o Minimum bed call count before applying a max call freq filter, defaults to 8.\n"+
				"-e Print verbose debugging output.\n"+


				"\nExample: java -jar pathToUSeq/Apps/VCFCallFrequency -f Vcf/ -s CFVcfs -v \n"+
				"    Hg38/Somatic/Avatar/Vcf -b Hg38/Somatic/Avatar/Bed -m 0.05 -c vcfCFConfig.txt \n\n" +
				"**************************************************************************************\n");

	}

}
