package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
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
import com.eclipsesource.json.JsonArray;
import com.eclipsesource.json.JsonObject;
import com.eclipsesource.json.JsonValue;  

/**Calculates call frequencies for each vcf record from a db of vcf files and callable region files.*/
public class VCFCallFrequency {

	private File[] vcfFiles;
	private File saveDirectory;
	private String queryURL = null;
	private String host = null;
	private String realm = null;
	private String queryKey = null;
	private String userName = null; 
	private String password = null; 
	private String fileFilter = null;
	private int numberRecordsPerQuery = 1000;
	private String vcfSearchUrl = null;
	private String bedSearchUrl = null;
	private String fetchOptionsUrl = null;
	private double[] idBedCount = null;
	private double maxCallFreq = 1.0;

	private boolean debug = false;
	private HttpClient client = new HttpClient();
	private Histogram hist = new Histogram(0,1,20);
	private int numToQueries;
	private int totalQueries = 0;
	private int totalQueriesWithHits = 0;
	private int numFailingCallFreq = 0;
	private int numPassingCallFreq = 0;

	public VCFCallFrequency (String[] args) {
		long startTime = System.currentTimeMillis();
		processArgs(args);

		buildSearchUrls(); 

		for (File vcf: vcfFiles) {
			IO.p(vcf.getName());
			annotate(vcf);
		}
		//print stats
		IO.p("\nTotal Vcf Queries: "+ totalQueries +", with matches: "+totalQueriesWithHits+", passing callFreq: "+numPassingCallFreq+", failing callFreq: "+numFailingCallFreq);
		IO.p("\nHistogram of call frequencies");
		hist.printScaledHistogram();

		//finish and calc run time

		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.p("\nDone! "+Math.round(diffTime)+" Sec\n");

	}

	private void buildSearchUrls() {
		//fetch key
		fetchEncodedQueryKey();

		//set search url
		vcfSearchUrl = queryURL+"search?key="+ queryKey+"&regExAll=.vcf.gz;"+fileFilter+"&matchVcf=true&vcf=";
		bedSearchUrl = queryURL+"search?key="+ queryKey+"&regExAll=.bed.gz;"+fileFilter+"&bed=";
		fetchOptionsUrl = queryURL+"search?key="+ queryKey+"&fetchOptions=true";
		IO.p("fetchOptions:\t"+fetchOptionsUrl);
		IO.p("VcfSearchUrl:\t"+vcfSearchUrl);
		IO.p("BedSearchUrl:\t"+bedSearchUrl);
		IO.p();
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
		System.out.print("\t");
		File modVcf = new File (saveDirectory, Misc.removeExtension(vcfFile.getName())+".callFreq.vcf.gz");
		String line = null;
		Gzipper out = null;
		BufferedReader in = null;
		boolean addInfo = true;
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
						out.println("##INFO=<ID=CallFreq,Number=1,Type=Float,Description=\"Frequency the variant has been observed in '"+fileFilter
								+ "' variant datasets.\">");
						addInfo = false;
					}
				}
				else {
					vcfRecords[counter] = Misc.TAB.split(line);
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
			IO.p();

			out.close();
			in.close();

		} catch (Exception e) {
			out.closeNoException();
			IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing this record "+line);
		} 
	}

	private void processRecordBlock(String[][] vcfRecords, Gzipper out) throws IOException {
		//create search strings, any vcfs to query?
		String[] searchStrings = makeSearchUrls(vcfRecords);

		if (searchStrings != null) {

			//execute vcf query
			String jsonResultsVcf = query(searchStrings[0]);
			//IO.p("VCF RETURN\n"+jsonResultsVcf);			
			JsonObject joVcf = JsonObject.readFrom(jsonResultsVcf);
			JsonValue jvVcf = joVcf.get("queryResults");

			//any results?
			if (jvVcf != null){
				JsonArray ja = jvVcf.asArray();
				int numWithHits = ja.size();
				System.out.print(numToQueries+"->"+numWithHits+", ");
				totalQueries+= numToQueries;
				totalQueriesWithHits+= numWithHits;


				//execute bed query, only do with vcf queries that return hits
				String jsonBedResults = query(searchStrings[1]);
				//IO.p("BED RETURN\n"+jsonBedResults);
				JsonObject joBed = JsonObject.readFrom(jsonBedResults);
				JsonValue jvBed = joBed.get("queryResults");
				loadBedResults(jvBed);

				//for each vcf query hit
				for (int i=0; i< numWithHits; i++){
					JsonObject joq = ja.get(i).asObject();
					// get vcf's index
					JsonValue input = joq.get("input");
					String vcfRecord = input.asString();

					String[] splitVcf = Misc.TAB.split(vcfRecord);
					int vcfIndex = Integer.parseInt(splitVcf[2]);
					// get number of intersecting vcf files
					double numVcfHits = (double)joq.get("numberHits").asInt();

					//get bed count and calc ratio
					double callRatio = 0;
					double bedCount = idBedCount[vcfIndex];
					if (bedCount > 0) callRatio = numVcfHits/bedCount;

					//modify INFO field of each record
					vcfRecords[vcfIndex][7] = "CallFreq="+Num.formatNumber(callRatio, 3)+";"+vcfRecords[vcfIndex][7];
					if (debug) IO.p("\nVCFRecord:"+vcfRecord+"\n\tNumVcfHits:"+numVcfHits+" NumBedHits:"+bedCount+" VcfIndex:"+ vcfIndex +
							"\n\tModVCFRec:"+Misc.stringArrayToString(vcfRecords[vcfIndex], "_"));

					//add info to histogram
					hist.count(callRatio);

					//fail it?
					if (callRatio > maxCallFreq) vcfRecords[vcfIndex][0] = "Failed";
				}
			}

			//add CallFreq=0 to those without any hits
			for (int i=0; i< vcfRecords.length; i++){
				String[] vcf = vcfRecords[i];
				if (vcf == null) break;
				if (vcf[7].startsWith("CallFreq=") == false) {
					if (vcf[7].equals("Failed")){
						vcf[7] = "CallFreq=0;"+vcf[7];
						hist.count(0);
						if (debug) IO.p("\nVCFRecordNoHit:"+Misc.stringArrayToString(vcf, "_"));
					}
				}

				//print it to file?
				if (vcf[0].equals("Failed")) numFailingCallFreq++;
				else {
					numPassingCallFreq++;
					out.print(vcf[0]);
					for (int x=1; x< vcf.length; x++){
						out.print("\t");
						out.print(vcf[x]);
					}
					out.println();
				}

			}
		}

	}

	private void loadBedResults(JsonValue jvBed) {
		//clear prior
		idBedCount = new double[numberRecordsPerQuery];

		// TODO pull bed input line and the corresponding vcf id from the bed name field
		JsonArray ja = jvBed.asArray();
		int numWithHits = ja.size();

		for (int i=0; i< numWithHits; i++){
			JsonObject joq = ja.get(i).asObject();
			double numHits = joq.get("numberHits").asInt();
			JsonValue input = joq.get("input");
			String bedRecord = input.asString();
			String[] splitbed = Misc.TAB.split(bedRecord);
			int index = Integer.parseInt(splitbed[3]);

			//set in array
			idBedCount[index] = numHits;
			//if (debug) IO.p("Bed "+input +" vcfIndex:"+index+ " numFileHits:"+ numHits);
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
		if (debug) IO.p("VcfBedToQuery:\n\t"+vcfBed[0]+"\n\t"+vcfBed[1]);

		//make cmd with key
		String urlVcfStr = vcfBed[0].insert(0, vcfSearchUrl).toString();
		String urlBedStr = vcfBed[1].insert(0, bedSearchUrl).toString();
		if (debug) IO.p("SearchURL:\n\t"+urlVcfStr+"\n\t"+urlBedStr);
		return new String[]{urlVcfStr,urlBedStr};
	}

	private StringBuilder[] pullVcfsBedsToQuery(String[][] vcfRecords) {
		numToQueries = 0;
		//pull vcfs to query
		StringBuilder sbVcf = new StringBuilder();
		StringBuilder sbBed = new StringBuilder();
		for (int i=0; i< vcfRecords.length; i++){
			String[] vcf = vcfRecords[i];
			if (vcf == null) break;
			numToQueries++;
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
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'q': queryURL = args[++i]; break;
					case 'h': host = args[++i]; break;
					case 'r': realm = args[++i]; break;
					case 'u': userName = args[++i]; break;
					case 'p': password = args[++i]; break;
					case 'f': fileFilter = args[++i]; break;
					case 'm': maxCallFreq = Double.parseDouble(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'd': debug = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		IO.p("\n"+IO.fetchUSeqVersion()+" Arguments:");
		IO.p("\t-v "+forExtraction);
		IO.p("\t-s "+saveDirectory);
		IO.p("\t-q "+queryURL);
		IO.p("\t-h "+host);
		IO.p("\t-r "+realm);
		IO.p("\t-f "+fileFilter);
		IO.p("\t-m "+maxCallFreq);
		IO.p("\t-d "+forExtraction);
		IO.p();


		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");



		//check params

		if (queryURL == null) Misc.printErrAndExit("\nError: provide a query URL, e.g. http://hci-clingen1.hci.utah.edu:8080/Query/");
		if (queryURL.endsWith("/") == false) queryURL = queryURL+"/";
		if (host == null) Misc.printErrAndExit("\nError: provide a host, e.g. hci-clingen1.hci.utah.edu");
		if (realm == null) Misc.printErrAndExit("\nError: provide a realm, e.g. QueryAPI");
		if (userName == null) Misc.printErrAndExit("\nError: provide a user name");
		if (password == null) Misc.printErrAndExit("\nError: provide a password");
		if (fileFilter == null) Misc.printErrAndExit("\nError: provide a fileFilter, e.g. /B37/Somatic/Avatar/ ");


		if (saveDirectory == null) Misc.printErrAndExit("\nError: provide a directory to save the annotated vcf files.");
		else saveDirectory.mkdirs();
		if (saveDirectory.exists() == false) Misc.printErrAndExit("\nError: could not find your save directory? "+saveDirectory);
	}	

	public static void printDocs(){
		IO.p("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Call Frequency: March 2018                      **\n" +
				"**************************************************************************************\n" +
				"Calculates a vcf call frequency for each variant when pointed at a genomic Query\n"+
				"sevice (https://github.com/HuntsmanCancerInstitute/Query). CallFreq's are calculated \n"+
				"by first counting the number of exact vcf matches present and dividing it by the\n"+
				"number of intersecting bed files. Use this to flag variants with high call rates that\n"+
				"are potential false positives. Use the file filter to limit which files to use in the\n"+
				"calculations. Only file paths containing the file filter will be included. Be sure to\n"+
				"place both the sample vcf and associated callable region bed files in the same Query\n"+
				"index folder.\n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				"-s Directory to save the annotated vcf files\n"+
				"-h Query service host, e.g. hci-clingen1.hci.utah.edu\n"+
				"-q Query service url, e.g. http://hci-clingen1.hci.utah.edu:8080/Query/\n"+
				"-r Query service realm, e.g. QueryAPI\n"+
				"-u Query service username\n"+
				"-p Query service password\n"+
				"-f Query service file filter, e.g. /B37/Somatic/Avatar/\n"+
				"-d Debug output\n"+


				"\nExample: java -jar pathToUSeq/Apps/VCFCallFrequency -v VCFss/ -s Annotated \n\n" +
				"**************************************************************************************\n");

	}

}
