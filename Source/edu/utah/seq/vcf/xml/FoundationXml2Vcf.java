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
 * 
 * @author david.nix@hci.utah.edu 
 **/
public class FoundationXml2Vcf {

	//user defined fields
	private File[] xmlFiles;
	private File indexedFasta = null;
	private File saveDirectory;
	private boolean skipBadVariants = false;
	
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
	
	//working data for a particular report
	private File workingXmlFile;
	private ArrayList<FoundationShortVariant> shortVariants = new ArrayList<FoundationShortVariant>();
	private ArrayList<FoundationCopyVariant> copyVariants = new ArrayList<FoundationCopyVariant>();
	private ArrayList<FoundationRearrangeVariant> rearrangeVariants = new ArrayList<FoundationRearrangeVariant>();
	private LinkedHashMap<String,String> reportAttributes = new LinkedHashMap<String,String>();
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
			
			//sampleInfo
			printSampleInfo();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running FoundationXml2Vcf app!");
		}
	}

	public void doWork() throws Exception{
		
		//Create fasta fetcher
		fasta = new IndexedFastaSequenceFile(indexedFasta);
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);
		
		//Get the DOM Builder Factory and builder
		factory = DocumentBuilderFactory.newInstance();
		builder = factory.newDocumentBuilder();

		System.out.println("Parsing and coverting...\n");
		for (int i=0; i< xmlFiles.length; i++){
			
			//clear any prior data
			reportAttributes.clear();
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
		}
		
		//close the fasta lookup fetcher
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
		sb.append(Misc.stringArrayToString(Misc.setToStringArray(info), "\n")); 
		sb.append("\n");
		
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
			for (int i = 0; i < nodeList.getLength(); i++) {
				Node node = nodeList.item(i);
				if (node instanceof Element) {
					//Is it the FinalReport?
					if (node.getNodeName().equals("FinalReport")) parseFinalReport(node);
					//Is it the variant-report?
					if (node.getNodeName().equals("variant-report")) parseVariantReport(node);
				}
			}
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
				if (name.equals("PMI")){
					parsePmiInfo(cNode);
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
				else if (name.equals("DOB")) dob = formatter.parse(cNode.getLastChild().getTextContent());
				else if (name.equals("CollDate")) collDate = formatter.parse(cNode.getLastChild().getTextContent());
			}
		}	
		
		if (collDate != null && dob != null){
			long duration  = collDate.getTime() - dob.getTime();
			double days = TimeUnit.MILLISECONDS.toDays(duration);
			String years = Num.formatNumberOneFraction(days/365.0);
			reportAttributes.put("age-at-collection", years);
		}
	}
	
	/**Adds the last child to the reportAttributes*/
	private void addLastChild(String key, Node node){
		Node last = node.getLastChild();
		if (last != null){
			String x = last.getTextContent().trim();
			if (x.length()!=0) reportAttributes.put(key, x);
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
					if (sampleAtt == null) Misc.printErrAndExit("\nError: failed to parse sample info from variant-report node.\n");
					reportAttributes.putAll(sampleAtt);
				}
				else if (name.equals("short-variants")) parseShortVariants(cNode);
				else if (name.equals("copy-number-alterations")) parseCopyVariants(cNode);
				else if (name.equals("rearrangements")) parseRearrangementVariants(cNode);
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
		NodeList sampleList = node.getChildNodes();
		for (int j = 0; j < sampleList.getLength(); j++) {
			Node cNode = sampleList.item(j); 
			if (cNode instanceof Element) {
				//sample? parse and add the attributes, should be only one
				if (cNode.getNodeName().equals("sample")) return parseNodeAttributes(cNode);
			}
		}		
		return null;
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


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Foundation Xml 2 Vcf: June 2016                      **\n" +
				"**************************************************************************************\n" +
				"Attempts to parse xml foundation reports to vcf. This is an inprecise process with\n"+
				"some insertions, multi snv, and multi vars. VCF variants have not been normalized.\n"+
				"Consider left aligning and demultiplexing.\n"+

				"\nOptions:\n"+
				"-x Path to a FoundationOne xml report or directory containing such.\n"+
				"-s Path to a directory for saving the results.\n"+
				"-f Path to the reference fasta with xxx.fai index\n"+
				"-o Skip variants that clearly fail to convert, e.g. var seq doesn't match fasta.\n"+
				"     Defaults to marking 'ci' in FILTER field.\n"+
				
				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/FoundationXml2Vcf -x /F1/TRF145179.xml\n" +
				"     -f /Ref/human_g1k_v37.fasta -s /F1/VCF/ \n\n" +

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
	
}
