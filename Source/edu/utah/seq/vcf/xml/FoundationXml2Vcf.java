package edu.utah.seq.vcf.xml;

import java.io.*;
import java.util.regex.*;
import javax.xml.parsers.*;
import org.w3c.dom.*;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import util.gen.*;
import java.util.*;


/**
 * Takes one or more patient xml reports from FoundationOne tests and converts the variants into vcf format.
 * 
 * @author david.nix@hci.utah.edu 
 **/
public class FoundationXml2Vcf {

	//user defined fields
	private File[] xmlFiles;
	private File indexedFasta = null;
	private File saveDirectory;
	
	//internal fields
	private String source;
	private DocumentBuilderFactory factory;
	private DocumentBuilder builder;
	private HashMap<String, String> seqRef = new HashMap<String,String>();
	private IndexedFastaSequenceFile fasta; 
	
	//working data for a particular report
	private File workingXmlFile;
	private ArrayList<FoundationShortVariant> shortVariants = new ArrayList<FoundationShortVariant>();
	private LinkedHashMap<String,String> reportAttributes = new LinkedHashMap<String,String>();

	//constructors
	public FoundationXml2Vcf(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running FoundationXml2Vcf app!");
		}
	}

	public void doWork() throws Exception{
		
		loadSeqRef();
		
		//Create fasta fetcher
		fasta = new IndexedFastaSequenceFile(indexedFasta);
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);
		
		//Get the DOM Builder Factory and builder
		factory = DocumentBuilderFactory.newInstance();
		builder = factory.newDocumentBuilder();

		System.out.println("Parsing and coverting...");
		for (int i=0; i< xmlFiles.length; i++){
			
			//clear any prior data
			reportAttributes.clear();
			shortVariants.clear();
			
			//process file
			workingXmlFile = xmlFiles[i];
			System.out.println("\t"+ workingXmlFile.getName());
			convert();
			
			//build vcf
			writeVcf();
			
		}
		
		//close the fasta lookup fetcher
		fasta.close();

	}
	
	private void writeVcf() {
		try {
			String name = Misc.removeExtension(workingXmlFile.getName());
			File vcf = new File (saveDirectory, name+".vcf");
			PrintWriter out = new PrintWriter( new FileWriter(vcf));
			
			//add header
			out.print(buildVcfHeader());
			
			//add short variants
//TODO: need to merge and sort
			for (FoundationShortVariant sv : shortVariants) out.println(sv.toVcf());
			
			//add others?
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: issue writing out vcf for "+workingXmlFile);
		}

	}

	private String buildVcfHeader(){
		StringBuilder sb = new StringBuilder();
		
		//add in standard header
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("##source=\""+source+"\"\n");
		sb.append("##file-path=\""+workingXmlFile+"\"\n");
		sb.append("##parse-date="+Misc.getDateNoSpaces()+"\n");
		
		//add in meta info
		for (String key: reportAttributes.keySet()){
			String value = reportAttributes.get(key);
			sb.append("##");
			sb.append(key);
			sb.append("=\"");
			sb.append(value);
			sb.append("\"\n");
		}
		
		//add in INFO
		if (shortVariants.size()!=0) FoundationShortVariant.appendInfoLines(sb);
		//others?
		
		//chrom line
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		return sb.toString();
	}

	/**For reversing the strand.*/
	private void loadSeqRef() {
		seqRef.put("G","C");
		seqRef.put("A","T");
		seqRef.put("T","A");
		seqRef.put("C","G");
	}

	public HashMap<String, String> getSeqRef() {
		return seqRef;
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
	
	/**Parses the variant-report section*/
	private void parseFinalReport(Node node) {
		
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);

			//Identifying the child tag of employee encountered. 
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				
				//samples? parse and add the attributes
				if (name.equals("PMI")){
					parseShortVariants(cNode);
				}
				
			}
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
				
				else if (name.equals("short-variants")){
					parseShortVariants(cNode);
				}
				
				else if (name.equals("copy-number-alterations")){
					
				}
				
				else if (name.equals("rearrangements")){
					
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
//System.out.println("VAR: "+var);
					FoundationShortVariant fv = new FoundationShortVariant(var, this);
					shortVariants.add(fv);
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
			//System.out.println("\t"+key+" : "+value);
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
				"**                              Foundation Xml 2 Vcf: May 2016                      **\n" +
				"**************************************************************************************\n" +
				"Note, variants have not been normalized.  Consider left aligning and demultiplexing.\n"+

				"\nRequired Options:\n"+
				"-r A regions bed file (chr, start, stop,...) to intersect, see\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1 , gz/zip OK.\n"+
				"-s Path to a directory for saving the results.\n"+
				"-t Number concurrent threads to run, will need > 4G RAM each.\n"+
				"-c GATK command to execute, see the example below, modify to match your enviroment.\n"+
				"     Most resources require full paths. Don't set -o or -L\n"+

				"\nExample: java -Xmx24G -jar pathToUSeq/Apps/GatkRunner -r /SS/targets.bed -t 8 -s\n" +
				"     /SS/HC/ -c 'java -Xmx4G -jar /SS/GenomeAnalysisTK.jar -T MuTect2 \n"+
				"    -R /SS/human_g1k_v37.fasta --dbsnp /SS/dbsnp_138.b37.vcf \n"+
				"    --cosmic /SS/v76_GRCh37_CosmicCodingMuts.vcf.gz -I:tumor /SS/sarc.bam -I:normal \n"+
				"    /SS/norm.bam'\n\n" +

				"**************************************************************************************\n");

	}
	
	public IndexedFastaSequenceFile getFasta() {
		return fasta;
	}
	
}
