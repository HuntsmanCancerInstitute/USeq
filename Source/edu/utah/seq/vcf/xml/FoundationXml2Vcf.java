package edu.utah.seq.vcf.xml;

import java.io.*;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.regex.*;
import javax.xml.parsers.*;
import org.w3c.dom.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import java.util.concurrent.TimeUnit;


/**
 * Takes one or more patient xml reports from FoundationOne tests and converts most of the variants into vcf format.
 * Outputs a summary clinical data spreadsheet
 * @author david.nix@hci.utah.edu 
 **/
public class FoundationXml2Vcf {

	//user defined fields
	private File[] xmlFiles;
	private File indexedFasta = null;
	private File saveDirectory;
	private boolean skipBadVariants = false;
	private boolean includePHI = false;
	private LinkedHashSet<String> keysToExport = null;
	
	//internal fields
	private String source;
	private DocumentBuilderFactory factory;
	private DocumentBuilder builder;
	private IndexedFastaSequenceFile fasta; 
	private long numShortPass = 0;
	private long numShortFail = 0;
	private long numCopyPass = 0;
	private long numCopyFail = 0;
	private long numRearrangePass = 0;
	private long numRearrangeFail = 0;
	private ArrayList<String> failingFiles = new ArrayList<String>();
	private String genomeVersion = "hg19";
	private ArrayList<SampleInfo> sampleInfo = new ArrayList<SampleInfo>();
	private HashMap<String, String> geneNameStrand = null;
	
	//working data for a particular report
	private File workingXmlFile;
	private ArrayList<FoundationShortVariant> shortVariants = new ArrayList<FoundationShortVariant>();
	private ArrayList<FoundationCopyVariant> copyVariants = new ArrayList<FoundationCopyVariant>();
	private ArrayList<FoundationRearrangeVariant> rearrangeVariants = new ArrayList<FoundationRearrangeVariant>();
	private LinkedHashMap<String,String> reportAttributes = null;
	private LinkedHashMap<String,String>[] allReportAttributes = null;
	private boolean problemParsing;

	//constructors
	public FoundationXml2Vcf(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();
			
			//stats
			System.out.println("\nParsing stats for "+xmlFiles.length+" files:");
			System.out.println(numShortPass+ "\t# Short Pass");
			System.out.println(numShortFail+ "\t# Short Fail");
			System.out.println(numCopyPass+ "\t# CNV Pass");
			System.out.println(numCopyFail+ "\t# CNV Fail");
			System.out.println(numRearrangePass+ "\t# Rearrange Pass");
			System.out.println(numRearrangeFail+ "\t# Rearrange Fail");

			printSampleInfo();
			
			printSpreadsheet();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running FoundationXml2Vcf app!");
		}
	}

	private void printSpreadsheet() throws IOException {
		IO.pl("Exporting patient summary spreadsheet... ");
		//merge all the keys
		LinkedHashSet<String> allKeys = new LinkedHashSet<String>();
		for (LinkedHashMap<String,String> s : allReportAttributes) allKeys.addAll(s.keySet());
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
		PrintWriter txtOut = new PrintWriter (new FileWriter(new File (saveDirectory, "aggregatePatientInfo.xls")));
		for (String k: keysToExport) {
			txtOut.print(k);
			txtOut.print("\t");
		}
		txtOut.println();
		
		//write a line per sample
		for (LinkedHashMap<String,String> s : allReportAttributes) {
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

	private void loadGeneNameStrand() {
		geneNameStrand = new HashMap<String, String>();
		for (int i=0; i< geneStrand.length; i++) {
			String name = geneStrand[i];
			i++;
			String strand = geneStrand[i];
			geneNameStrand.put(name, strand);
		}
		
	}

	public void doWork() throws Exception{
		
		//Create fasta fetcher
		fasta = new IndexedFastaSequenceFile(indexedFasta);
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);
		
		loadGeneNameStrand();
		
		//Get the DOM Builder Factory and builder
		factory = DocumentBuilderFactory.newInstance();
		builder = factory.newDocumentBuilder();
		allReportAttributes = new LinkedHashMap[xmlFiles.length];

		System.out.println("Parsing and coverting...");
		for (int i=0; i< xmlFiles.length; i++){
			
			//clear any prior data
			reportAttributes = new LinkedHashMap<String,String>();
			shortVariants.clear();
			copyVariants.clear();
			rearrangeVariants.clear();
			problemParsing = false;
			
			//process file
			workingXmlFile = xmlFiles[i];
			System.out.println(workingXmlFile.getName());
			convert();
			
			//build vcf
			writeVcf();
			
			//create an SI
			SampleInfo si = new SampleInfo(reportAttributes);
			if (si.allLoaded() == false) {
				System.err.println("WARNING: missing sample info for:");
				System.err.println(si.toString());
				problemParsing = true;
			}
			sampleInfo.add(si);
			
			if (problemParsing) failingFiles.add(workingXmlFile.getName());
			allReportAttributes[i] = reportAttributes;
		}
		
		//close the IO
		fasta.close();

	}

	private void printSampleInfo() {
		System.out.println();
		System.out.println( SampleInfo.fetchSummaryInfo(sampleInfo) );
	}

	private void writeVcf() {
		try {
			String name = Misc.removeExtension(workingXmlFile.getName());
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
		sb.append("##file-path="+workingXmlFile+"\n");
		sb.append("##parse-date="+Misc.getDateNoSpaces()+"\n");
		
		//add in meta info
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
		}
		//add FILTER?
		if (problemParsing) sb.append("##FILTER=<ID=ci,Description=\"Converting from xml to vcf issue, treat skeptically.\"\n");
		
		//add ALT and INFO lines
		TreeSet<String> alt = new TreeSet<String>();
		TreeSet<String> info = new TreeSet<String>();
		
		if (shortVariants.size()!=0) FoundationShortVariant.appendInfoLines(info);
		if (copyVariants.size()!=0) FoundationCopyVariant.appendInfoAltLines(alt, info);
		if (rearrangeVariants.size()!=0) FoundationRearrangeVariant.appendInfoAltLines(alt, info);
		
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
			//Load and Parse the XML document
			Document document = builder.parse(workingXmlFile);

			//Iterating through the nodes and extract the data
			NodeList nodeList = document.getDocumentElement().getChildNodes();
			boolean foundFinalReport = false;
			boolean foundVariantReport = false;
			for (int i = 0; i < nodeList.getLength(); i++) {
				Node node = nodeList.item(i);
				if (node instanceof Element) {	
					//Is it the FinalReport?
					if (node.getNodeName().equals("FinalReport")) {
						foundFinalReport = true;
						parseFinalReport(node);
					}
					//Is it the variant-report?
					else if (node.getNodeName().equals("variant-report")) {
						foundVariantReport = true;
						parseVariantReport(node);
					}
					//must watch out for new style reports
					else if (node.getNodeName().equals("rr:ResultsPayload")) {
						NodeList nodeList2 = node.getChildNodes();
						for (int x = 0; x < nodeList2.getLength(); x++) {
							Node node2 = nodeList2.item(x);
							if (node2 instanceof Element) {	
								//Is it the FinalReport?
								if (node2.getNodeName().equals("FinalReport")) {
									foundFinalReport = true;
									parseFinalReport(node2);
								}
								//Is it the variant-report?
								else if (node2.getNodeName().equals("variant-report")) {
									foundVariantReport = true;
									parseVariantReport(node2);
								}
							}
						}
						break;
					}
				}
			}
			if (foundFinalReport == false || foundVariantReport == false) throw new Exception("Failed to find the FinalReport("+foundFinalReport+") or the variant-report("+foundVariantReport+")");
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing this xml file "+workingXmlFile);
		}
	}
	
	/**Parses the variant-report section
	 * @throws ParseException 
	 * @throws DOMException */
	private void parseFinalReport(Node node) throws DOMException, ParseException {
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//PMI?
				if (name.equals("PMI")) parsePmiInfo(cNode);
				//Genes?
				else if (name.equals("Genes")) parseGenes(cNode);
			}
		}
	}
	
	private void parseGenes(Node node) throws DOMException, ParseException {
		NodeList sampleList = node.getChildNodes();
		
		for (int j = 0; j < sampleList.getLength(); j++) {
			Node cNode = sampleList.item(j);
			String name = cNode.getNodeName();

			if (cNode instanceof Element && name.equals("Gene")) {
				NodeList geneNodes = cNode.getChildNodes();
				//for each node under a particular Gene
				for (int i=0; i<geneNodes.getLength(); i++) {
					Node gNode = geneNodes.item(i);
					String gName = gNode.getNodeName();
					//check the Name
					if (gName.equals("Name")) {
						Node last = gNode.getLastChild();
						if (last != null){
							String x = last.getTextContent().trim();
							//is it Microsatellite status?
							if (x.equals("Microsatellite status")) parseMSIOrTMB(geneNodes, "MSI");
							//is it Tumor Mutation Burden?
							else if (x.equals("Tumor Mutation Burden")) parseMSIOrTMB(geneNodes, "TMB");
						}
					}
				}
			}
		}	
	}

	
	private void parseMSIOrTMB(NodeList nodes, String key) {
		//for each node under a particular Gene
		for (int i=0; i<nodes.getLength(); i++) {
			Node mNode = nodes.item(i);
			String name = mNode.getNodeName();
			//is it the Alterations?
			if (name.equals("Alterations")) {
				NodeList altList = mNode.getChildNodes();
				for (int j=0; j<altList.getLength(); j++) {
					Node aNode = altList.item(j);
					if (aNode.getNodeName().equals("Alteration")) {
						NodeList aList = aNode.getChildNodes();
						for (int x=0; x<aList.getLength(); x++) {
							Node last= aList.item(x);
							if (last.getNodeName().equals("Name")) {
								addLastChild(key, last);
								return;
							}
						}
					}
				}
			}
		}
	}

	private void parsePmiInfo(Node node) throws DOMException, ParseException {
		NodeList sampleList = node.getChildNodes();
		Date collDate = null;
		Date dob = null;
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd");
		
		for (int j = 0; j < sampleList.getLength(); j++) {
			Node cNode = sampleList.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				if (name.equals("Gender")) addLastChild("Gender", cNode);
				else if (name.equals("OrderingMD")) addLastChild("OrderingMD", cNode);
				else if (name.equals("SpecSite")) addLastChild("SpecSite", cNode);
				else if (name.equals("SubmittedDiagnosis")) addLastChild("SubmittedDiagnosis", cNode);
				else if (name.equals("CollDate")) {
					collDate = formatter.parse(cNode.getLastChild().getTextContent());
					if (includePHI) addLastChild("CollDate", cNode);
				}
				else if (name.equals("DOB")) {
					dob = formatter.parse(cNode.getLastChild().getTextContent());
					if (includePHI) addLastChild("DOB", cNode);
				}
				else if (name.equals("ReportId")) addLastChild("ReportId", cNode);
				else if (includePHI) {
					if (name.equals("MRN")) addLastChild("MRN", cNode);
					else if (name.equals("FirstName")) addLastChild("FirstName", cNode);
					else if (name.equals("LastName")) addLastChild("LastName", cNode);
					else if (name.equals("ReceivedDate")) addLastChild("ReceivedDate", cNode);
				}
			}
		}	
		
		if (collDate != null && dob != null){
			long duration  = collDate.getTime() - dob.getTime();
			double days = TimeUnit.MILLISECONDS.toDays(duration);
			String years = Num.formatNumberOneFraction(days/365.0);
			reportAttributes.put("AgeAtCollection", years);
		}
	}
	
	/**Adds the last child to the reportAttributes*/
	private void addLastChild(String key, Node node){
		Node last = node.getLastChild();
		if (last != null){
			String x = last.getTextContent().trim();
			if (x.length()!=0) reportAttributes.put(key, x);
			//if (key.equals("MSI") || key.equals("TMB"))IO.pl("Adding "+key+" : "+x);
		}
	}

	/**Parses the variant-report section*/
	private void parseVariantReport(Node node) {
		
		//parse attributes in the variant-report line, e.g. disease, disease-ontology, gender, diagnosis, ....
		LinkedHashMap<String,String> att = parseNodeAttributes(node);
		
		//remove xml attributes
		att.remove("xmlns");
		att.remove("xmlns:xsi");
		att.remove("xsi:schemaLocation");
		att.remove("gender");
		att.remove("test-request");
		reportAttributes.putAll(att);
		
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);

			//Identifying the child tag of employee encountered. 
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				
				//samples? parse and add the attributes
				if (name.equals("samples")){
					LinkedHashMap<String,String> sampleAtt = parseSamples(cNode);
					if (sampleAtt.size() == 0) Misc.printErrAndExit("\nError: failed to parse sample info from variant-report node.\n");
					reportAttributes.putAll(sampleAtt);
				}
				else if (name.equals("short-variants")) parseShortVariants(cNode);
				else if (name.equals("copy-number-alterations")) parseCopyVariants(cNode);
				else if (name.equals("rearrangements")) parseRearrangementVariants(cNode);
				else if (name.equals("quality-control")) reportAttributes.put("quality-control", cNode.getAttributes().getNamedItem("status").getNodeValue());
			}
		}
	}
	
	/**Parses out rearrangement variants in the form, <rearrangement description="fusion" in-frame="Yes" other-gene="ERG" 
	 * pos1="chr21:42868152-42868305" pos2="chr21:39872810-39873056" status="known" supporting-read-pairs="16" targeted-gene="TMPRSS2">*/
	private void parseRearrangementVariants(Node node) {
		NodeList sampleList = node.getChildNodes();
		int idCounter = 0;
		for (int j = 0; j < sampleList.getLength(); j++) {
			Node cNode = sampleList.item(j);
			if (cNode instanceof Element) {
				if (cNode.getNodeName().equals("rearrangement")) {
					LinkedHashMap<String, String> var = parseNodeAttributes(cNode);
					FoundationRearrangeVariant fc = new FoundationRearrangeVariant(var, this, idCounter++);
					if (fc.isFailedParsing()) numRearrangeFail++;
					else {
						rearrangeVariants.add(fc);
						numRearrangePass++;
					}
				}
			}
		}	
	}
	
	/**Parses out multiple copy number variants in the form, <copy-number-alteration copy-number="7" equivocal="true" gene="KRAS" 
	 * number-of-exons="5 of 5" position="chr12:25362722-25398327" ratio="1.47" status="known" type="amplification"/>*/
	private void parseCopyVariants(Node node) {
		NodeList sampleList = node.getChildNodes();
		for (int j = 0; j < sampleList.getLength(); j++) {
			Node cNode = sampleList.item(j);
			if (cNode instanceof Element) {
				if (cNode.getNodeName().equals("copy-number-alteration")) {
					LinkedHashMap<String, String> var = parseNodeAttributes(cNode);
					FoundationCopyVariant fc = new FoundationCopyVariant(var, this);
					if (fc.isFailedParsing()) numCopyFail++;
					else {
						copyVariants.add(fc);
						numCopyPass++;
					}
				}
			}
		}	
	}


	/**Parses out multiple short variants in the form, <short-variant cds-effect="394C&gt;T" depth="872" functional-effect="missense" 
	 * gene="IDH1" percent-reads="42.0" position="chr2:209113113" protein-effect="R132C" status="known" strand="-" transcript="NM_005896">*/
	private void parseShortVariants(Node node) {
		NodeList sampleList = node.getChildNodes();
		for (int j = 0; j < sampleList.getLength(); j++) {
			Node cNode = sampleList.item(j);
			if (cNode instanceof Element) {
				if (cNode.getNodeName().equals("short-variant")) {
					LinkedHashMap<String, String> var = parseNodeAttributes(cNode);
					FoundationShortVariant fv = new FoundationShortVariant(var, this);
					if (fv.isFailedParsing()) {
						numShortFail++;
						problemParsing = true;
						if (skipBadVariants == false) shortVariants.add(fv);
					}
					else {
						shortVariants.add(fv);
						numShortPass++;
					}
				}
			}
		}	
	}

	/**Parses <samples>  <sample bait-set="T7" mean-exon-depth="963.88" name="SA-1352726" nucleic-acid-type="DNA"/> in the variant-report */
	private LinkedHashMap<String, String> parseSamples(Node node) {
		LinkedHashMap<String, String> att = new LinkedHashMap<String, String>();
		NodeList sampleList = node.getChildNodes();
		for (int j = 0; j < sampleList.getLength(); j++) {
			Node cNode = sampleList.item(j); 
			if (cNode instanceof Element) {
				//sample? parse and add the attributes, should be only one, not so! now seeing RNA
				if (cNode.getNodeName().equals("sample")) {
					LinkedHashMap<String, String> kv = parseNodeAttributes(cNode);
					String nt = kv.remove("nucleic-acid-type").toLowerCase();
					for (String k: kv.keySet()) {
						att.put(nt+ "-"+ k, kv.get(k));
					}
				}
			}
		}		
		return att;
	}

	/**Loads a hashmap with all of the key: value attribute pairs of the node, not its children though.*/
	private LinkedHashMap<String, String> parseNodeAttributes(Node node) {
		LinkedHashMap<String, String> att = new LinkedHashMap<String, String>();
		NamedNodeMap map = node.getAttributes();
		int numAtt = map.getLength();
		for (int i=0; i< numAtt; i++){
			Node n =  map.item(i);
			String key = n.getNodeName();
			String value = n.getNodeValue();
			att.put(key, value);
		}
		return att;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FoundationXml2Vcf(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		source = useqVersion+" Args: "+ Misc.stringArrayToString(args, " ");
		System.out.println("\n"+ source +"\n");
		String attributes = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'x': xmlFiles = IO.extractFiles(new File(args[++i]), "xml"); break;
					case 'f': indexedFasta = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'o': skipBadVariants = true; break;
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

		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save directory!\n"+saveDirectory);
		saveDirectory.mkdirs();
		if (saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: your save directory does not appear to be a directory?\n");
		if (xmlFiles == null || xmlFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any xml files to parse?!\n");
		if (indexedFasta == null || indexedFasta.canRead() == false) Misc.printErrAndExit("\nError: cannot find your fasta reference?\n");
		
		if (attributes !=null) {
			String[] at = Misc.COMMA.split(attributes);
			keysToExport = new LinkedHashSet<String>();
			for (String s: at) keysToExport.add(s);
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Foundation Xml 2 Vcf: April 2020                       **\n" +
				"**************************************************************************************\n" +
				"Parses xml foundation reports to vcf. This is an inprecise process with some\n"+
				"insertions, multi snv, and multi vars. Output vcf variants have not been normalized.\n"+
				"Recommend running vt normalize, see https://genome.sph.umich.edu/wiki/Vt\n"+

				"\nOptions:\n"+
				"-x Path to a FoundationOne xml report or directory containing such.\n"+
				"-s Path to a directory for saving the results.\n"+
				"-f Path to the reference fasta with xxx.fai index used by Foundation for coordinates.\n"+
				"-o Skip variants that clearly fail to convert, e.g. var seq doesn't match fasta.\n"+
				"     Defaults to marking 'ci' in FILTER field.\n"+
				"-i Include PHI in spreadsheet output.\n"+
				"-a Print this list of attributes in spreadsheet, comma delimited, case sensitive, no\n"+
				"     spaces. Defaults to all.\n"+
				
				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/FoundationXml2Vcf -x FXmlReports/\n" +
				"     -f /Ref/b37.fasta -s FVcfs/ -a ReportId,Gender,SubmittedDiagnosis,MSI,TMB\n\n" +

				"**************************************************************************************\n");

	}
	
	public IndexedFastaSequenceFile getFasta() {
		return fasta;
	}

	public boolean isSkipBadVariants() {
		return skipBadVariants;
	}

	public String getGenomeVersion() {
		return genomeVersion;
	}
	
	/*Table of gene name and strand, needed for a few cases where Foundation doesn't report the strand for an effected gene. Arrrrggggg! */
	private String[] geneStrand  = new String[]{"CDC14A","+", "KIF1B","+", "PRMT6","+", "TNFRSF14","+", "VAV3","-", "CSF1","+", "RBM15","+", "MTOR","-", "WNT2B","+", "PTPN22","-",
			"TRIM33","-", "NRAS","-", "NGF","-", "MAD2L2","-", "VTCN1","-", "FAM46C","+", "MTHFR","-", "NPPB","-", "HSD3B1","+", "NOTCH2","-", "DVL1","-",
			"PRDM2","+", "NOTCH2NL","-", "CHD1L","+", "BCL9","+", "APH1A","-", "MCL1","-", "CTSS","-", "ARNT","-", "SETDB1","+", "MLLT11","+", "FLG","-",
			"NPR1","+", "CRTC2","-", "CREB3L4","+", "HAX1","+", "CTRC","+", "IL6R","+", "SHC1","-", "CKS1B","+", "ZBTB7B","+", "FDPS","+", "ASH1L","-",
			"RIT1","-", "RAB25","+", "LMNA","+", "VHLL","-", "IQGAP3","-", "HDGF","-", "PRCC","+", "NTRK1","+", "INSRR","-", "PEAR1","+", "ETV3L","-",
			"ETV3","-", "CD1D","+", "SPEN","+", "NCSTN","+", "EPHA2","-", "SDHC","+", "FCGR2A","+", "FCGR3A","-", "DDR2","+", "CDK11B","-", "PBX1","+",
			"SDHB","-", "CDK11A","-", "PRRX1","+", "FASLG","+", "ABL2","-", "PAX7","+", "PTGS2","-", "CDC73","+", "NEK7","+", "PLA2G2A","-", "CACNA1S","-",
			"ASCL5","-", "TNNT2","-", "ELF3","+", "LGR6","+", "UBE2T","-", "KDM5B","-", "BTG2","+", "PIK3C2B","-", "MDM4","+", "CDK18","+", "ELK4","-",
			"IKBKE","+", "TRAF3IP3","+", "IRF6","-", "NEK2","-", "SMYD2","+", "TGFB2","+", "CDC42","+", "WNT4","-", "SKI","+", "TLR5","-", "EPHA8","+",
			"H3F3A","+", "PARP1","-", "ITPKB","-", "PSEN2","+", "EPHB2","+", "JMJD4","-", "WNT9A","-", "WNT3A","+", "KDM1A","+", "EGLN1","-", "ID3","-",
			"MDS2","+", "RYR2","+", "IL22RA1","-", "FH","-", "IFNLR1","-", "SMYD3","-", "RUNX3","-", "EXTL1","+", "CNKSR1","+", "ARID1A","+", "MAP3K6","-",
			"FGR","-", "PTPRC","+", "PRDM16","+", "LCK","+", "HDAC1","+", "SFPQ","-", "CSF3R","-", "ZC3H12A","+", "MYCL","-", "TIE1","+", "MPL","+",
			"CDC20","+", "PTPRF","+", "KDM4A","+", "ARTN","+", "PLK3","+", "PTCH2","-", "MUTYH","-", "MAST2","+", "PIK3R3","-", "RAD54L","+", "TAL1","-",
			"CMPK1","+", "AKT3","-", "CDKN2C","+", "GLIS1","-", "PCSK9","+", "JUN","-", "CYP2J2","-", "CHD5","-", "NFIA","+", "ROR1","+", "HES2","-",
			"JAK1","-", "LEPR","+", "IL23R","+", "IL12RB2","+", "CAMTA1","+", "DIRAS3","-", "MSH4","+", "FUBP1","-", "TNFRSF9","-", "ERRFI1","-", "ADGRL2","+",
			"PRKACB","+", "BCL10","-", "RBMXL1","-", "BRDT","+", "GFI1","-", "RPL5","+", "BCAR3","-", "PIK3CD","+", "DPYD","-", "HES4","-", "CHUK","-",
			"WNT8B","+", "HIF1AN","+", "PAX2","+", "TLX1","+", "BTRC","+", "FGF8","-", "LDB1","-", "NFKB2","+", "SUFU","+", "CYP17A1","-", "NT5C2","-",
			"SMC3","+", "SHOC2","+", "TCF7L2","+", "ADRB1","+", "GRK5","+", "FGFR2","-", "BUB3","+", "MKI67","-", "MGMT","+", "ZMYND11","+", "SUV39H2","+",
			"CACNB2","+", "ALOX5","+", "MLLT10","+", "BMI1","+", "ABI1","-", "MAP3K8","+", "CREM","+", "KLF6","-", "RET","+", "NCOA4","-", "MAPK8","+",
			"ERCC6","-", "UBE2D1","+", "IL15RA","-", "CCDC6","-", "IL2RA","-", "CDK1","+", "ARID5B","+", "EGR2","-", "JMJD1C","-", "CTNNA3","-", "TET1","+",
			"PRF1","-", "C10orf54","-", "KAT6B","+", "NUTM2B","+", "GATA3","+", "NRG3","+", "WAPL","-", "BMPR1A","+", "NUTM2A","+", "KLLN","-", "PTEN","+",
			"ACTA2","-", "FAS","+", "IFIT2","+", "IFIT3","+", "IFIT1","+", "TNKS2","+", "TBC1D12","+", "HELLS","+", "CYP2C19","+", "CYP2C9","+", "CYP2C8","-",
			"BLNK","-", "GTPBP4","+", "GOT1","-", "NKX2-3","+", "ABCC2","+", "PGR","-", "YAP1","+", "BIRC3","+", "BIRC2","+", "KCNQ1","+", "DYNC2H1","+",
			"PDGFD","-", "GUCY1A2","-", "ATM","+", "C11orf65","-", "POU2AF1","-", "PPP2R1B","-", "SDHD","+", "DRD2","-", "ZBTB16","+", "CDKN1C","-",
			"APOA1","-", "IL10RA","+", "UBE4A","+", "KMT2A","+", "DDX6","-", "CBL","+", "CHEK1","+", "TEAD1","+", "ETS1","-", "FLI1","+", "KCNJ5","+",
			"PRDM10","-", "CYP2R1","-", "PIK3C2A","-", "CTSD","-", "MYOD1","+", "H19","-", "PRMT3","+", "IGF2","-", "FANCF","-", "ASCL2","-", "LGR4","-",
			"BDNF","-", "PAX6","-", "WT1","-", "LMO2","-", "ELF5","-", "EHF","+", "CD44","+", "TRAF6","-", "NUP98","-", "RRM1","+", "EXT2","+", "PRDM11","+",
			"CREB3L1","+", "DDB2","+", "MYBPC3","-", "SPI1","-", "PTPRJ","+", "FOLH1","-", "CTNND1","+", "MS4A1","+", "SDHAF2","+", "FEN1","+", "HRAS","-",
			"SLC22A6","-", "VEGFB","+", "ESRRA","+", "SF1","-", "MAP4K2","-", "MEN1","-", "MAP3K11","-", "RELA","-", "KAT5","+", "FOSL1","-", "RBM14","+",
			"KDM2A","+", "RPS6KB2","+", "AIP","+", "GSTP1","+", "KMT5B","-", "LRP5","+", "CCND1","+", "FGF19","-", "FGF4","-", "FGF3","-", "FADD","+",
			"PPFIA1","+", "RBMXL2","+", "PHOX2A","-", "SLCO2B1","+", "WNT11","-", "EMSY","+", "PAK1","-", "RSF1","-", "GAB2","-", "LMO1","-", "PICALM","-",
			"TRIM66","-", "EED","+", "ASCL3","-", "FAT3","+", "MRE11A","-", "KDM4D","+", "WEE1","+", "CEP57","+", "MAML2","-", "SPIC","+", "IGF1","-",
			"ASCL1","+", "STYK1","-", "PRDM4","-", "ASCL4","+", "MYL2","-", "SH2B3","+", "PTPN11","+", "TBX3","-", "ETV6","+", "TAOK3","-", "HNF1A","+",
			"LRP6","-", "KDM2B","-", "SETD1B","+", "BCL7A","+", "CLIP1","-", "KMT5A","+", "NCOR2","-", "CDKN1B","+", "POLE","-", "WNT5B","+", "PIK3C2G","+",
			"CACNA1C","+", "SLCO1B3","+", "SLCO1B1","+", "SLCO1A2","-", "KRAS","-", "KDM5A","-", "FOXM1","-", "TEAD4","+", "CAPRIN2","-", "PKP2","-",
			"PRMT8","+", "LRRK2","+", "CCND2","+", "ADAMTS20","-", "FGF23","-", "FGF6","-", "RAD51AP1","+", "ARID2","+", "HDAC7","-", "VDR","-", "WNT10B","-",
			"WNT1","+", "KMT2D","-", "DHH","-", "PRPF40B","+", "SMARCD1","+", "ATF1","+", "POU6F1","-", "ACVR1B","+", "NR4A1","+", "RARG","-", "ESPL1","+",
			"MAP3K12","-", "HOXC13","+", "HOXC11","+", "HOXC10","+", "CBX5","-", "NTF3","+", "CDK2","+", "ERBB3","+", "STAT2","-", "BAZ2A","-", "NAB2","+",
			"STAT6","-", "GLI1","+", "DDIT3","-", "CDK4","-", "WIF1","-", "CHD4","-", "HMGA2","+", "ING4","-", "ZNF384","-", "DYRK2","+", "LAG3","+",
			"MDM2","+", "YEATS4","+", "PTPN6","+", "FRS2","+", "PTPRB","-", "PTPRR","-", "LGR5","+", "WNK1","+", "E2F7","-", "PTPRQ","+", "KITLG","-",
			"RAD52","-", "BTG1","-", "ELK3","+", "CDK17","-", "IKBIP","-", "FGF14","-", "ERCC5","+", "IRS2","-", "ING1","+", "SOX1","+", "CUL4A","+",
			"DCUN1D2","-", "TPTE2","-", "LATS2","-", "FGF9","+", "PARP4","-", "CDK8","+", "CDX2","-", "FLT3","-", "FLT1","-", "HSPH1","-", "BRCA2","+",
			"PDS5B","+", "SMAD9","-", "FOXO1","-", "ELF1","-", "HTR2A","-", "NUDT15","+", "RB1","+", "SETDB2","+", "ATP7B","-", "NEK5","-", "NEK3","-",
			"DIAPH3","-", "DACH1","-", "DIS3","-", "KLF5","+", "KLF12","-", "LMO7","+", "SPRY2","-", "GPC5","+", "SOX21","-", "YY1","+", "HSP90AA1","-",
			"TRAF3","+", "XRCC3","-", "AKT1","-", "PARP2","+", "APEX1","+", "PRMT5","-", "AJUBA","-", "CEBPE","-", "BCL2L2","+", "MYH7","-", "REC8","+",
			"TINF2","-", "RIPK3","-", "NFATC4","+", "FOXG1","+", "PRKD1","-", "BAZ1A","-", "NFKBIA","-", "NKX2-1","-", "NKX2-8","-", "PAX9","+", "MIPOL1","+",
			"FOXA1","-", "SSTR1","+", "FANCM","+", "SOS2","-", "MAP4K5","-", "SAV1","-", "CDKN3","+", "HIF1A","+", "ESR2","-", "FNTB","+", "MAX","-",
			"RAD51B","+", "MAP3K9","-", "PSEN1","+", "NUMB","-", "PGF","-", "MLH3","-", "NEK9","-", "FOS","+", "TSHR","+", "PTPN21","-", "FOXN3","-",
			"DICER1","-", "TCL1B","+", "TCL1A","-", "BCL11B","-", "SETD3","-", "NOP10","-", "NUTM1","+", "ACTC1","-", "B2M","+", "SPRED1","+", "BUB1B","+",
			"PAK6","+", "RAD51","+", "LTK","-", "TYRO3","+", "MGA","+", "JMJD7","+", "TP53BP1","-", "FBN1","-", "SCG5","+", "SHC4","-", "GREM1","+",
			"FGF7","+", "MAPK6","+", "TCF12","+", "TPM1","+", "MAP2K1","+", "SMAD6","+", "SMAD3","+", "MAP2K5","+", "SKOR1","+", "CD276","+", "PML","+",
			"CYP1A2","+", "CSK","+", "NRG4","-", "BCL2A1","-", "ARNT2","+", "NTRK3","-", "FANCI","+", "IDH2","-", "IQGAP1","+", "CRTC3","+", "BLM","+",
			"FES","+", "CHD2","+", "IGF1R","+", "SSTR5","+", "CIITA","+", "SOCS1","-", "TNFRSF17","+", "MYH11","-", "ERCC4","+", "MKL2","+", "GLIS2","+",
			"ABCC1","+", "NTHL1","-", "TSC2","+", "TRAF7","+", "MLST8","+", "PALB2","-", "PLK1","+", "PDPK1","+", "CES1","-", "KDM8","+", "IL4R","+",
			"IL21R","+", "SULT1A1","-", "AXIN1","-", "CD19","+", "TAOK2","+", "MAPK3","-", "RNF40","+", "SETD1A","+", "ZNF668","-", "VKORC1","-", "KAT8","+",
			"PRSS8","-", "FUS","+", "SLX4","-", "CREBBP","-", "ZNF423","-", "BRD7","-", "NOD2","+", "CYLD","+", "CHD9","+", "NUP93","+", "USB1","+",
			"SETD6","+", "CDH11","-", "CDH5","+", "CES2","+", "CBFB","+", "HSD11B2","+", "CTCF","+", "NFATC3","+", "JMJD8","-", "PRMT7","+", "CDH3","+",
			"CDH1","+", "HAS3","+", "NQO1","-", "PHLPP2","-", "ZFHX3","-", "MPG","+", "MAF","-", "PLCG2","+", "FOXL1","+", "CBFA2T3","-", "CDK10","+",
			"FANCA","-", "MC1R","+", "PRDM7","-", "GRIN2A","-", "SOX8","+", "MAP2K4","+", "WNT3","-", "WNT9B","+", "NCOR1","-", "FLCN","-", "COPS3","-",
			"SMYD4","-", "GID4","+", "RPA1","+", "HNF1B","-", "MAPK7","+", "SLC47A1","+", "SLC47A2","-", "MAP2K3","+", "MLLT6","+", "YWHAE","-", "NEK8","+",
			"TAOK1","+", "SLC6A4","-", "NF1","+", "SUZ12","+", "RHOT1","+", "RAD51D","-", "CRK","-", "CDK12","+", "STARD3","+", "ERBB2","+", "GRB7","+",
			"IKZF3","-", "CDC6","+", "RARA","+", "TOP2A","-", "SMARCE1","-", "JUP","-", "FKBP10","+", "KAT2A","-", "STAT5B","-", "STAT5A","+", "STAT3","-",
			"EZH1","-", "WNK4","+", "BRCA1","-", "ETV4","-", "TAF15","+", "G6PC3","+", "HDAC5","-", "MAP3K14","-", "CBX1","-", "MINK1","+", "HOXB3","-",
			"CRHR1","+", "HOXB13","-", "PHB","-", "SPOP","-", "KAT7","+", "ABCC3","+", "RABEP1","+", "NLRP1","-", "HLF","+", "MSI2","+", "RNF43","-",
			"RAD51C","+", "RPS6KB1","+", "PPM1D","+", "TBX2","+", "BRIP1","-", "TLK2","+", "ACE","+", "MAP3K3","+", "SMARCD2","-", "CD79B","-", "DDX5","-",
			"SMURF2","-", "GNA13","-", "AXIN2","-", "BPTF","+", "PRKAR1A","+", "MAP2K6","+", "SOX9","+", "DVL2","-", "GPS2","-", "SSTR2","+", "TNK1","+",
			"FGF11","+", "GRB2","-", "CDK3","+", "TP53","-", "JMJD6","-", "SRSF2","-", "TMC6","-", "TMC8","+", "BIRC5","+", "KDM6B","+", "CHD3","+",
			"CBX2","+", "CBX8","-", "CBX4","-", "RNF213","+", "RPTOR","+", "ALOX12B","-", "AATK","-", "ASPSCR1","+", "AURKB","-", "CSNK1D","-", "CTC1","-",
			"ZNF750","-", "PTPN2","-", "ROCK1","-", "ESCO1","-", "MIB1","+", "GATA6","+", "ZNF521","-", "SS18","-", "SMCHD1","+", "CDH2","-", "DSC2","-",
			"DSG2","+", "ASXL3","+", "PIK3C3","+", "SETBP1","+", "SMAD2","-", "SMAD7","-", "MBD1","-", "MAPK4","+", "SMAD4","+", "DCC","+", "TCF4","-",
			"MALT1","+", "PMAIP1","+", "CDH20","+", "PHLPP1","+", "BCL2","-", "KDSR","-", "CDH7","+", "TYMS","+", "YES1","-", "PTPRM","+", "NFATC1","+",
			"DNMT1","-", "TYK2","-", "KEAP1","-", "DNM2","+", "CARM1","+", "SMARCA4","+", "LDLR","+", "CNOT3","+", "EPOR","-", "STK11","+", "JUNB","+",
			"MAST1","+", "CALR","+", "NFIX","+", "LYL1","-", "PRKACA","-", "RPS15","+", "NOTCH3","-", "BRD4","-", "MBD3","-", "CYP4F2","-", "TCF3","-",
			"BABAM1","+", "JAK3","-", "IL12RB1","-", "PIK3R2","+", "JUND","-", "CRTC1","+", "MEF2B","-", "MAU2","+", "PBX4","-", "DOT1L","+", "CCNE1","+",
			"GNA11","+", "ELANE","+", "TSHZ3","-", "CEBPA","-", "CEBPG","+", "NFIC","+", "FZR1","+", "CD22","+", "ETV2","+", "KMT2B","+", "PSENEN","+",
			"NFKBID","-", "ALKBH6","-", "ZNF607","-", "MATK","-", "SPRED3","+", "RYR1","+", "MAP4K1","-", "NFKBIB","+", "PAK4","+", "IFNL3","-", "MED29","+",
			"PIAS4","+", "MAP3K10","+", "AKT2","-", "NUMBL","-", "CYP2A6","-", "MAP2K2","-", "CYP2B6","+", "AXL","+", "TGFB1","-", "SHC2","-", "CD79A","+",
			"POU2F2","-", "GSK3A","-", "ERF","-", "CIC","+", "XRCC1","-", "SH3GL1","-", "BCL3","+", "CBLC","+", "RELB","+", "ERCC2","-", "ERCC1","-",
			"FOSB","+", "FOXA3","+", "ARHGAP35","+", "BBC3","-", "CARD8","-", "LMTK3","-", "FGF21","+", "PPP1R15A","+", "BAX","+", "UHRF1","+", "NTF4","-",
			"TEAD2","-", "FLT3LG","+", "PRMT1","+", "KDM4B","+", "POLD1","+", "SPIB","+", "PPP2R1A","+", "BIRC8","-", "TNNI3","-", "KMT5C","+", "U2AF2","+",
			"ZNF444","+", "ZNF471","+", "PEG3","-", "AURKC","+", "NRTN","+", "TRIM28","+", "MLLT1","-", "PSPN","-", "FGF22","+", "KHSRP","-", "CD70","-",
			"VAV1","+", "INSR","-", "MAP2K7","+", "MAP4K4","+", "IL1R2","+", "IL1R1","+", "IL18R1","+", "IL18RAP","+", "ODC1","-", "NCK2","+", "RANBP2","+",
			"ABCB11","-", "BUB1","-", "BCL2L11","+", "ROCK2","-", "MERTK","+", "PAX8","-", "E2F6","-", "GLI2","+", "TRIB2","+", "ERCC3","-", "MAP3K2","-",
			"MAP3K19","-", "CXCR4","-", "SPOPL","+", "LRP1B","-", "ZEB2","-", "ACVR2A","+", "NR4A2","-", "ACVR1","-", "TANC1","+", "BAZ2B","-", "MYCN","+",
			"XIRP2","+", "TLK1","-", "PDK1","+", "SP3","-", "HOXD13","+", "HOXD11","+", "HOXD10","+", "HOXD3","+", "HOXD4","+", "SMC6","-", "HNRNPA3","+",
			"NFE2L2","-", "GEN1","+", "PPP1R1C","+", "COL3A1","+", "PMS1","+", "NAB1","+", "STAT1","-", "STAT4","-", "SF3B1","-", "SGO2","+", "CASP8","+",
			"CDK15","+", "CD28","+", "CTLA4","+", "ICOS","+", "RHOB","+", "CREB1","+", "IDH1","-", "APOB","-", "ERBB4","-", "IKZF2","-", "BARD1","-",
			"ATIC","+", "STK36","+", "WNT6","+", "WNT10A","+", "FEV","-", "IHH","-", "EPHA4","-", "PAX3","-", "CUL3","-", "IRS1","-", "SP110","-", "SP140","+",
			"SP140L","+", "SP100","+", "DIS3L2","+", "UGT1A9","+", "UGT1A4","+", "UGT1A1","+", "ATAD2B","-", "TRAF3IP1","+", "TWIST2","+", "HDAC4","-",
			"NCOA1","+", "DNMT3A","-", "ASXL2","-", "FOSL2","+", "ALK","-", "CEBPZ","-", "CYP1B1","-", "SOS1","-", "MAP4K3","-", "EML4","+", "EPCAM","+",
			"MSH2","+", "MSH6","+", "FBXO11","-", "LHCGR","-", "FSHR","-", "FANCL","-", "BCL11A","-", "REL","+", "XPO1","-", "PDCD1","-", "SPRED2","-",
			"MXD1","+", "PCBP1","+", "TGFA","-", "SMYD5","+", "TET3","+", "MOB1A","-", "TLX2","+", "CTNNA2","+", "TCF7L1","+", "GGCX","-", "KDM3A","+",
			"ID2","+", "SMYD1","+", "ADAM17","-", "YWHAQ","-", "TMEM127","-", "ZAP70","+", "AFF3","-", "SEC23B","+", "NKX2-4","-", "NKX2-2","-", "PAX1","+",
			"FOXA2","-", "SSTR4","+", "ID1","+", "BCL2L1","-", "HCK","+", "PLAGL2","-", "ASXL1","+", "DNMT3B","+", "CBFA2T2","+", "E2F1","-", "GFRA4","-",
			"SRC","+", "CDC25B","+", "TRIB3","+", "MAFB","-", "TOP1","+", "PLCG1","+", "CHD6","-", "PTPRT","-", "MYBL2","+", "YWHAB","+", "STK4","+",
			"CD40","+", "ZMYND8","-", "NCOA3","+", "PTGIS","-", "CEBPB","+", "NFATC2","-", "ZNF217","-", "AURKA","-", "CTCFL","-", "GNAS","+", "SS18L1","+",
			"GATA5","-", "PTK6","-", "SRMS","-", "ARFRP1","-", "PRPF6","+", "TPTE","+", "NRIP1","-", "GABPA","+", "BACH1","+", "OLIG2","+", "IFNAR2","+",
			"IL10RB","+", "IFNAR1","+", "RUNX1","-", "SETD4","-", "CBR3","+", "ERG","-", "ETS2","+", "BRWD1","-", "TMPRSS2","-", "RIPK4","-", "PRDM15","-",
			"ABCG1","+", "U2AF1","-", "ICOSLG","-", "DNMT3L","-", "BTG3","-", "SLC19A1","-", "PRMT2","+", "IFNGR2","+", "IL17RA","+", "CECR2","+", "BID","-",
			"COMT","+", "CRKL","+", "LZTR1","+", "MAPK1","-", "BCR","+", "SMARCB1","+", "MN1","-", "CHEK2","-", "XBP1","-", "ZNRF3","+", "EWSR1","+",
			"NF2","+", "SF3A1","-", "PATZ1","-", "YWHAH","+", "CSF2RB","+", "IL2RB","-", "SSTR3","-", "RAC2","-", "CARD10","-", "SOX10","-", "CSNK1E","-",
			"CBX6","-", "CBX7","-", "PDGFB","-", "MKL1","-", "EP300","+", "ZC3H7B","+", "TEF","+", "SMC1B","-", "WNT7B","-", "PPARA","+", "CYP2D6","-",
			"BRD1","-", "PIM3","+", "HDAC10","-", "MAPK12","-", "MAPK11","-", "GSTT1","-", "FANCD2","+", "TFG","+", "VHL","+", "NFKBIZ","+", "CBLB","-",
			"TIGIT","+", "ZBTB20","-", "VGLL4","-", "CD80","-", "GSK3B","-", "SLC15A2","+", "CD86","+", "CASR","+", "HSPBAP1","-", "DIRC2","+", "PPARG","+",
			"MYLK","-", "UMPS","+", "RAF1","-", "RUVBL1","-", "GATA2","-", "RPN1","-", "PIK3R4","-", "NEK11","+", "RYK","-", "EPHB1","+", "HDAC11","+",
			"STAG1","-", "NCK1","+", "IL20RB","+", "WNT7A","-", "PIK3CB","-", "FOXL2","-", "TMEM43","+", "XPC","-", "ATR","-", "HLTF","-", "WWTR1","-",
			"MED12L","+", "CCNL1","-", "MLF1","+", "SMC4","+", "MECOM","-", "TERC","-", "PRKCI","+", "SKIL","+", "TBL1XR1","-", "PIK3CA","+", "ZNF639","+",
			"SOX2","+", "DCUN1D1","-", "KLHL6","-", "DVL3","+", "THPO","-", "EPHB3","+", "C3orf70","-", "MAP3K13","+", "ETV5","-", "BCL6","-", "LPP","+",
			"TP63","+", "IL1RAP","+", "FGF12","-", "HES1","+", "TNK2","-", "PAK2","+", "KAT2B","+", "SGO1","-", "RARB","+", "TOP2B","-", "NEK10","-",
			"TGFBR2","+", "IL5RA","-", "CRBN","-", "MLH1","+", "MYD88","+", "ACVR2B","+", "SCN5A","-", "ZNF620","+", "CTNNB1","+", "SETMAR","+", "MYL3","-",
			"SETD2","-", "SMARCC1","-", "MAP4","-", "CDC25A","-", "RHOA","-", "MST1","-", "MST1R","-", "TLR9","-", "BAP1","-", "PBRM1","-", "NEK4","-",
			"IL17RB","+", "WNT5A","-", "FHIT","-", "PTPRG","+", "MITF","+", "FOXP1","-", "ROBO2","+", "CADM2","+", "VGLL3","-", "EPHA3","+", "SRGAP3","-",
			"SETD5","+", "EPHA6","+", "BRPF1","+", "IL17RC","+", "NFKB1","+", "UBE2D3","-", "CENPE","-", "TET2","+", "INTS12","-", "LEF1","-", "EGF","+",
			"MAD2L1","-", "PRDM5","-", "FGF2","+", "FAT4","+", "PLK4","+", "JADE1","+", "NKX3-2","-", "ELF2","-", "SETD7","-", "MAML3","-", "IL15","+",
			"INPP4B","-", "GAB1","+", "SMARCA5","+", "SMAD1","+", "ARHGAP10","+", "NR3C2","-", "FBXW7","-", "TLR2","+", "PDGFC","-", "PALLD","+", "NEK1","-",
			"FBXO8","-", "VEGFC","-", "FGFR3","+", "IRF2","-", "FAT1","-", "WHSC1","+", "NSD2","+", "SLIT2","+", "RBPJ","+", "GRK4","+", "TLR10","-",
			"TLR1","-", "TLR6","-", "PDS5A","-", "RHOH","+", "PHOX2B","-", "TXK","-", "TEC","-", "CHIC2","-", "PDGFRA","+", "KIT","+", "KDR","-", "ADGRL3","+",
			"EPHA5","-", "GNRHR","-", "MOB1B","+", "EPGN","+", "EREG","+", "AREG","+", "BTC","-", "PRDM8","+", "FGF5","+", "FAM175A","-", "MAPK10","-",
			"PTPN13","+", "AFF1","+", "ABCG2","-", "BMPR1B","+", "RAP1GDS1","+", "FER","+", "APC","+", "DMXL1","+", "PRDM6","+", "TERT","-", "ACSL6","-",
			"IL3","+", "IRF1","-", "RAD50","+", "TCF7","+", "SMAD5","+", "WNT8A","+", "BRD8","-", "CDC25C","-", "KDM3B","+", "EGR1","+", "CTNNA1","+",
			"TMEM173","-", "UBE2D2","+", "NRG2","-", "HBEGF","-", "TRIO","+", "DIAPH1","-", "HDAC3","-", "FGF1","-", "ARHGAP26","+", "NR3C1","-",
			"PPP2R2B","-", "SPINK1","-", "ADRB2","+", "CSF1R","-", "PDGFRB","-", "CDX1","+", "FAT2","-", "HAVCR2","-", "ITK","+", "EBF1","-", "PTTG1","+",
			"TENM2","+", "DOCK2","+", "TLX3","+", "NPM1","+", "FGF18","+", "FBXW11","-", "NKX2-5","-", "DRD1","-", "FGFR4","+", "NSD1","+", "NHP2","-",
			"MAML1","+", "CDK7","+", "MAPK9","-", "FLT4","-", "SDHA","+", "PRDM9","+", "CDH10","-", "PRLR","-", "IL7R","+", "SKP2","+", "NIPBL","+",
			"GDNF","-", "LIFR","-", "OSMR","+", "RICTOR","-", "CARD6","+", "GHR","+", "FGF10","-", "IL6ST","-", "MAP3K1","+", "SETD9","+", "PLK2","-",
			"HTR1A","-", "SLC6A3","-", "PIK3R1","+", "HMGCR","+", "IQGAP2","+", "MTRR","+", "DHFR","-", "MSH3","+", "BRD9","-", "RASA1","+", "POU5F2","-",
			"CHD1","-", "PRDM1","+", "MAK","-", "FOXO3","+", "HLA-G","+", "CDK19","-", "TRAF3IP2","-", "FYN","-", "HDAC2","-", "FRK","-", "VGLL2","+",
			"ROS1","-", "HLA-A","+", "HEY2","+", "RSPO3","+", "FOXQ1","+", "SGK1","-", "MYB","+", "BCLAF1","-", "MAP3K5","-", "IL20RA","-", "IL22RA2","-",
			"IFNGR1","-", "TNFAIP3","+", "PERP","-", "ECT2L","+", "PLAGL1","-", "ZC3H12D","-", "LATS1","-", "ESR1","+", "SYNE1","-", "JARID2","+",
			"ARID1B","+", "SOD2","-", "IGF2R","+", "SLC22A1","+", "SLC22A2","-", "SLC22A3","+", "MAP3K4","+", "PARK2","-", "QKI","+", "TBP","+", "HLA-E","+",
			"TPMT","-", "KDM1B","+", "DEK","-", "ID4","+", "E2F3","+", "DDR1","+", "POU5F1","-", "HLA-C","-", "HIST1H3B","-", "HLA-B","-", "HIST1H1E","+",
			"HIST1H4E","+", "TNF","+", "DUSP22","+", "RIPK1","+", "EHMT2","-", "PTPRK","-", "STK19","+", "HLA-DRB5","-", "HLA-DRB6","-", "CYP21A2","+",
			"PBX2","-", "NOTCH4","-", "HMGA1","+", "SPDEF","-", "PPARD","+", "FANCE","+", "TEAD3","-", "FKBP5","-", "MAPK14","+", "MAPK13","+", "BRPF3","+",
			"ETV7","-", "CDKN1A","+", "PIM1","+", "HLA-DRB1","-", "HLA-DQA1","+", "HLA-DQB1","-", "IRF4","+", "HLA-DQA2","+", "HLA-DQB2","-", "TAP2","-",
			"TAP1","-", "HLA-DOB","-", "HLA-DMB","-", "HLA-DMA","-", "FOXP4","+", "BRD2","+", "TFEB","-", "FRS3","-", "CCND3","-", "HLA-DOA","-",
			"HLA-DPA1","-", "HLA-DPB1","+", "HLA-DPB2","+", "PTK7","+", "POLH","+", "VEGFA","+", "HSP90AB1","+", "NFKBIE","-", "DAXX","-", "RUNX2","+",
			"PHF1","+", "PKHD1","-", "ICK","-", "RAB23","-", "ADGRB3","+", "DSP","+", "PHIP","-", "TBX18","-", "PNRC1","+", "BACH2","-", "MAP3K7","-",
			"EPHA7","-", "HLA-DRA","+", "HLA-F","+", "PRDM13","+", "EPHB4","-", "CUX1","+", "KMT2E","+", "RINT1","+", "PIK3CG","+", "SLC26A3","-", "FOXP2","+",
			"TFEC","-", "MET","+", "WNT2","-", "CFTR","+", "WNT16","+", "WASL","-", "POT1","-", "GRM8","-", "PAX4","-", "IRF5","+", "SMO","+", "CREB3L2","-",
			"TRIM24","+", "ETV1","-", "TBXAS1","+", "KDM7A","-", "BRAF","-", "EPHA1","-", "EZH2","-", "KCNH2","-", "CDK5","-", "SMARCD3","-", "RHEB","-",
			"PRKAG2","-", "KMT2C","-", "XRCC2","-", "PAXIP1","-", "SHH","-", "MNX1","-", "AHR","+", "HDAC9","+", "MAD1L1","-", "TWIST1","-", "CBX3","+",
			"HOXA3","-", "HOXA9","-", "HOXA10","-", "HOXA11","-", "HOXA13","-", "JAZF1","-", "CARD11","-", "FKBP9","+", "POU6F2","+", "CDK13","+", "INHBA","-",
			"GLI3","-", "UBE2D4","+", "PDGFA","-", "IKZF1","+", "GRB10","-", "PRKAR1B","-", "EGFR","+", "PMS2","-", "RAC1","+", "SBDS","-", "BAZ1B","-",
			"PRSS1","+", "POR","+", "MAGI2","-", "PRSS2","+", "GLCCI1","+", "SEMA3C","-", "HGF","-", "GRM3","+", "ABCB4","-", "ABCB1","-", "EPHB6","+",
			"CDK14","+", "AKAP9","+", "CDK6","-", "SAMD9","-", "KEL","-", "SEM1","-", "ASNS","-", "LMTK2","+", "TRRAP","+", "SMURF1","-", "ARPC1A","+",
			"ARPC1B","+", "CYP3A5","-", "CYP3A4","-", "YWHAZ","-", "UBR5","-", "RSPO2","-", "BLK","+", "GATA4","+", "RAD21","-", "EXT1","-", "ATAD2","-",
			"RNF139","+", "TRIB1","+", "POU5F1B","+", "CASC11","-", "MYC","+", "PVT1","+", "WISP1","+", "NDRG1","-", "PTK2","-", "RECQL4","-", "TUSC3","+",
			"FGF20","-", "NAT2","+", "MAPK15","+", "FGF17","+", "HR","-", "NKX3-1","-", "NKX2-6","-", "PPP2R2A","+", "PTK2B","+", "ESCO2","+", "ELP3","+",
			"WRN","+", "NRG1","+", "ZNF703","+", "ADGRA2","+", "ASH2L","+", "LSM1","-", "BAG4","+", "NSD3","-", "FGFR1","-", "IDO1","+", "SFRP1","-",
			"KAT6A","-", "IKBKB","+", "CEBPD","-", "PRKDC","-", "SOX17","+", "LYN","+", "MOS","-", "PLAG1","-", "CHD7","+", "MCPH1","+", "MYBL1","-",
			"PREX2","+", "C8orf34","+", "PRDM14","-", "NCOA2","-", "TERF1","+", "ELOC","-", "HEY1","-", "ZNF704","-", "E2F5","+", "RIPK2","+", "NBN","-",
			"RUNX1T1","-", "RAD54B","-", "CCNE2","-", "TNKS","+", "MTDH","+", "STK3","-", "SMC2","+", "ABCA1","-", "TAL2","+", "KLF4","-", "MUSK","+",
			"TLR4","+", "TRAF1","-", "PTGS1","+", "NEK6","+", "PPP6C","-", "PBX3","+", "CDK9","+", "ENG","-", "SET","+", "PRRX2","+", "PRDM12","+", "ABL1","+",
			"NUP214","+", "TSC1","-", "GFI1B","+", "DBH","+", "VAV2","-", "BRD3","-", "RXRA","+", "NOTCH1","-", "TRAF2","+", "EHMT1","+", "NFIB","-",
			"PSIP1","-", "SMARCA2","+", "MLLT3","-", "MTAP","+", "CDKN2A","-", "CDKN2B","-", "TEK","+", "TAF1L","-", "PRSS3","+", "CNTFR","-", "IL11RA","+",
			"FANCG","-", "PAX5","-", "ZCCHC7","+", "ZBTB5","-", "SHB","-", "GLIS3","-", "JAK2","+", "CD274","+", "PDCD1LG2","+", "UHRF2","+", "KDM4C","+",
			"SMC5","+", "GNAQ","-", "PTPRD","-", "NTRK2","+", "CTSL","+", "CDK20","-", "SHC3","-", "CKS2","+", "SYK","+", "ROR2","-", "WNK2","+", "PHF2","+",
			"NUTM2F","-", "FANCC","-", "PTCH1","-", "NUTM2G","+", "XPA","-", "GALNT12","+", "TGFBR1","+", "NR4A3","+", "BTK","-", "GLA","-", "MID1","-",
			"IRS4","-", "PAK3","+", "IL13RA1","+", "CRLF2","-", "ZBTB33","+", "CUL4B","-", "XIAP","+", "STAG2","+", "SH2D1A","+", "CSF2RA","+", "TLR7","+",
			"TLR8","+", "SMARCA1","-", "BCORL1","+", "ELF4","-", "GPC3","-", "IL3RA","+", "PHF6","+", "VGLL1","+", "CD40LG","+", "RBMX","-", "FGF13","-",
			"SOX3","-", "P2RY8","-", "FANCB","-", "AFF2","+", "MAMLD1","+", "IRAK1","-", "G6PD","-", "DKC1","+", "BMX","+", "MTCP1","-", "ZRSR2","+",
			"MAP3K15","-", "EIF1AX","-", "ARX","-", "NR0B1","-", "BCOR","-", "USP9X","+", "DDX3X","+", "MAOA","+", "KDM6A","+", "RBM10","+", "CDK16","+",
			"ARAF","+", "ELK1","-", "SSX1","+", "SSX3","-", "SSX4","+", "WAS","+", "SUV39H1","+", "GATA1","+", "HDAC6","+", "PIM2","-", "TFE3","-",
			"FOXP3","-", "CCNB3","+", "MAGED1","+", "SSX2","-", "KDM5C","-", "SMC1A","-", "PHF8","-", "WNK3","-", "KLF8","+", "SPRY3","+", "IL9R","+",
			"AMER1","-", "TBX22","+", "AR","+", "FOXO4","+", "IL2RG","-", "MED12","+", "ZMYM3","-", "NONO","+", "TAF1","+", "HDAC8","-", "CHIC1","+",
			"FGF16","+", "ATRX","-", "BRWD3","-", "DACH2","+", "TBL1X","+", "DIAPH2","+", "USP9Y","+", "UTY","-", "KDM5D","-"};

	public HashMap<String, String> getGeneNameStrand() {
		return geneNameStrand;
	}
	
}
