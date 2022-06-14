package edu.utah.seq.vcf.xml.caris;

import java.io.*;
import java.text.ParseException;
import java.util.regex.*;
import javax.xml.parsers.*;
import org.w3c.dom.*;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.gen.*;
import java.util.*;

/**
 * Takes one or more patient xml reports from Caris tests and converts most of the variants into vcf format by using the Caris vcf to pull genomic coordinates which are not provided in the xml
 * Outputs a summary clinical data spreadsheet
 * @author david.nix@hci.utah.edu 
 **/
public class CarisXmlVcfParser {

	//user defined fields
	private HashMap<String, File[]> xmlVcfFiles = null;
	private File saveDirectory;
	private boolean includePHI = false;
	private boolean saveAllVcfRecords = false;
	private boolean skipUnPaired = false;

	//internal fields
	private String source;
	private DocumentBuilderFactory factory;
	private DocumentBuilder builder; 
	private String genomeVersion = "##reference=hg38";
	private HashSet<String> ihcTestNames = null;
	private LinkedHashMap<String,String> workingReportAttributes = null;
	private LinkedHashMap<String,String>[] allReportAttributes = null;
	private String[] statLines = null;
	private HashMap<String, UCSCGeneLine[]> name2GeneModels = null;

	//working data for a particular report
	private File workingXmlFile;
	private String workingReportName = null;
	private int workingNumberFailedGenomicMatches = 0;

	//vcf info
	private File workingVcfFile;
	private StringBuilder workingVcfHeader;
	private LinkedHashMap<String, SimpleVcf> workingVcfs = new LinkedHashMap<String, SimpleVcf>();

	//xml objects
	private ArrayList<GenomicAlteration> workingGenomicAlterations = new ArrayList<GenomicAlteration>();
	private ArrayList<CNVAlteration> workingCNVAlterations = new ArrayList<CNVAlteration>();
	private ArrayList<Translocation> workingTranslocations = new ArrayList<Translocation>();
	private ArrayList<ExpressionAlteration> workingExpressionAlterations = new ArrayList<ExpressionAlteration>();
	private ArrayList<MethylationAlteration> workingMethylationAlterations = new ArrayList<MethylationAlteration>();



	//constructors
	public CarisXmlVcfParser(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			parseDatasets();
			
			printSpreadsheet();
			
			//print the stats lines
			IO.pl("\nParsing Statistics:");
			IO.pl("Name\tVcfRecords\tGenomicAlterationReports\tCNVs\tFusionReports\tIHCTests\tFailedGenomicAltMatches");
			for (String s: statLines) IO.pl(s);

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Sec\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running the CarisXmlVcfParser app.");
		}
	}

	private void printSpreadsheet() throws IOException {
		IO.pl("\nExporting patient summary spreadsheet... ");
		//merge all the keys
		LinkedHashSet<String> allKeys = new LinkedHashSet<String>();
		for (LinkedHashMap<String,String> s : allReportAttributes) allKeys.addAll(s.keySet());

		//make writer add header
		PrintWriter txtOut = new PrintWriter (new FileWriter(new File (saveDirectory, "aggregatePatientInfo.xls")));
		for (String k: allKeys) {
			txtOut.print(k);
			txtOut.print("\t");
		}
		txtOut.println();

		//write a line per sample
		for (LinkedHashMap<String,String> s : allReportAttributes) {
			//for each key
			for (String k: allKeys) {
				String v = s.get(k);
				if (v!=null) txtOut.print(v);
				txtOut.print("\t");
			}
			txtOut.println();
		}
		txtOut.close();
	}
	
//Add it MGMT-Me methylation status like an IHC test to the spreadsheet

	public void parseDatasets() throws Exception{

		//build hash for ihc expression tests
		buildIHCHash();

		//Get the DOM Builder Factory and builder
		factory = DocumentBuilderFactory.newInstance();
		builder = factory.newDocumentBuilder();
		allReportAttributes = new LinkedHashMap[xmlVcfFiles.size()];
		statLines = new String[allReportAttributes.length];

		System.out.println("Parsing and coverting...");
		Iterator<String> it = xmlVcfFiles.keySet().iterator();
		int index = 0;
		while (it.hasNext()) {
			workingReportName = it.next();
			File[] xv = xmlVcfFiles.get(workingReportName);
			workingXmlFile = xv[0];
			workingVcfFile = xv[1];
//IO.pl("\tXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\t"+workingReportName);

			//clear any prior data
			workingReportAttributes = new LinkedHashMap<String,String>();
			workingVcfs.clear();
			workingGenomicAlterations.clear();
			workingCNVAlterations.clear();
			workingTranslocations.clear();
			workingExpressionAlterations.clear();
			workingMethylationAlterations.clear();
			workingNumberFailedGenomicMatches = 0;
			
			//add xml name, often contains the patient id
			workingReportAttributes.put("XmlFileName", Misc.removeExtension(workingXmlFile.getName()));

			loadVcf();

			loadXml();

			matchVcfWithGenomicAlts();

			saveVcfFile();
			
			saveCnvFile();
			
			saveFusionFile();

			statLines[index] = statParsing();
			allReportAttributes[index++] = workingReportAttributes;
		}
	}

	private void matchVcfWithGenomicAlts() {
		//for each GermlineAlteration find the corresponding vcf record
		for (GenomicAlteration ga: workingGenomicAlterations) {
			
			//check if present, HLA- tests don't have this nor any vcf like record so just skip
			String patho = ga.getPathogenicity();
			if (patho == null) continue;
			
			patho = patho.toLowerCase();
			//skip wt
			if (patho.contains("wild") || patho.contains("indeterminate")) continue;
			String gaKey = ga.fetchKey();
			SimpleVcf vcf = workingVcfs.get(gaKey);
			if (vcf == null) {
				//IO.el(ga.toString());
				IO.el("Failed to find the vcf for the xml variant "+gaKey+" in the vcf file "+ workingVcfFile+ 
						"\nConsider manually matching by editing the xml file. Skipping xml variant.");
				workingNumberFailedGenomicMatches++;
			}
			else {
				patho = Misc.WHITESPACE.matcher(patho).replaceAll("_");
				vcf.appendInfo("XRV;XRP="+patho+";");
			}
		}
	}

	private String statParsing() {
		ArrayList<String> al = new ArrayList<String>();
		//name
		al.add(workingReportName);
		//# vcf records
		al.add(new Integer(workingVcfs.size()).toString());
		//# genomic alterations
		al.add(new Integer(workingGenomicAlterations.size()).toString());
		//# cnv alterations
		al.add(new Integer(workingCNVAlterations.size()).toString());		
		//# translocation alterations
		al.add(new Integer(workingTranslocations.size()).toString());	
		//# IHC alterations
		al.add(new Integer(workingExpressionAlterations.size()).toString());
		//# Failed to match genomic alterations
		al.add(workingNumberFailedGenomicMatches+"");
		return Misc.stringArrayListToString(al, "\t");
	}

	private void buildIHCHash() {
		String[] ihcNames = {"p16","PD-L1 (SP263)", "PD-L1 (22c3)", "PD-L1 (SP142)", "PD-L1 FDA(SP142)", "PD-L1 FDA (28-8)", 
				"MLH1", "PMS2", "MSH2", "MSH6", "ALK", "PTEN", "Mismatch Repair Status", "Her2/Neu", 
				"TrkA/B/C", "Androgen Receptor", "Folfox Responder Similarity Score", "ROS1", "ER" };
		ihcTestNames = new HashSet<String>();
		for (String n: ihcNames) ihcTestNames.add(n);

	}

	private void loadVcf() throws IOException {
		workingVcfHeader = new StringBuilder();
		BufferedReader in = IO.fetchBufferedReader(workingVcfFile);
		String line;
		while ((line = in.readLine())!=null) {
			if (line.startsWith("#")) {
				//reference? check it
				if (line.startsWith("##reference=")) {
					if (line.toLowerCase().equals(genomeVersion) == false) Misc.printErrAndExit("ERROR: the reference genome does not appear to be "+genomeVersion);
				}
				//CHROM line
				if (line.startsWith("#CHROM")) {
					workingVcfHeader.append(SimpleVcf.xrv);
					workingVcfHeader.append(SimpleVcf.xrp);
				}
				//fix their malformed header line that is causing a vt to error
				if (line.startsWith("##FORMAT=<ID=SA,")) line = "##FORMAT=<ID=SA,Number=1,Type=String,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling\">";
				workingVcfHeader.append(line); 
				workingVcfHeader.append("\n");
			}
			else {
				SimpleVcf vcf = new SimpleVcf(line,0);
				String key = vcf.fetchCarisKey();
				if (key == null) key = Misc.getRandomString(10);
				workingVcfs.put(key, vcf);
			}


		}
		in.close();
	}

	private void saveVcfFile() {
		if (workingVcfs.size() == 0) return;
		Gzipper out = null;
		try {
			//for each vcf record
			ArrayList<String> toSave = new ArrayList<String>();
			int counter = 1;
			for (SimpleVcf sv : workingVcfs.values()) {
				if (saveAllVcfRecords) {
					sv.setId("Caris_"+counter);
					counter++;
					toSave.add(sv.toString());
				}
				else if (sv.getInfo().startsWith("XRV;XRP=")) {
					sv.setId("Caris_"+counter);
					counter++;
					toSave.add(sv.toString());
				}
			}
			//if (toSave.size() == 0) return; don't do this always want a vcf file output even if empty
			
			File vcf = new File (saveDirectory, workingReportName+".vcf.gz");
			out = new Gzipper(vcf);

			//add header
			out.print(workingVcfHeader);
			for (String s: toSave) out.println(s);

			out.close();
		} catch (IOException e) {
			if (out != null) {
				out.closeNoException();
				out.getGzipFile().deleteOnExit();
			}
			e.printStackTrace();
			Misc.printErrAndExit("\nError: issue writing out vcf for "+workingXmlFile);
		}
	}
	
	private void saveCnvFile() {
		if (workingCNVAlterations.size() == 0) return;
		Gzipper out = null;
		try {
			
			
			ArrayList<SimpleBed> bedToSave = new ArrayList<SimpleBed>();
			for (int i=0; i< workingCNVAlterations.size(); i++) {
				CNVAlteration cnv = workingCNVAlterations.get(i);
				String res = cnv.getResult().toLowerCase();
				if (res.contains("not detected") || res.contains("indeterminate")) continue;
				bedToSave.add(cnv.toBed(name2GeneModels));
			}
			if (bedToSave.size() == 0) return;
			
			SimpleBed[] bedRegions = new SimpleBed[bedToSave.size()];
			bedToSave.toArray(bedRegions);
			Arrays.sort(bedRegions);
			
			File bedFile = new File (saveDirectory, workingReportName+".cnv.bed.gz");
			out = new Gzipper(bedFile);

			//add header
			out.println("#ChrHg38\tGeneStart\tGeneStop\tGeneName:Type\tScore\tGeneStrand");

			//for each bed
			for (SimpleBed cnv: bedRegions) {
				if (cnv!= null) out.println(cnv);
			}

			out.close();
		} catch (IOException e) {
			if (out != null) {
				out.closeNoException();
				out.getGzipFile().deleteOnExit();
			}
			e.printStackTrace();
			Misc.printErrAndExit("\nError: issue writing out vcf for "+workingXmlFile);
		}
	}
	
	private void saveFusionFile() {
		if (workingTranslocations.size() == 0) return;
		Gzipper out = null;
		try {
			
			ArrayList<SimpleBed> beds = new ArrayList<SimpleBed>();
			int num = workingTranslocations.size();
			for (int i=0; i< num; i++) {
				Translocation trans = workingTranslocations.get(i);
				String res = trans.getResult().toLowerCase();
				if (res.contains("not detected") || res.contains("indeterminate")) continue;
				
				SimpleBed[] sbs = workingTranslocations.get(i).toBed(name2GeneModels);
				for (SimpleBed sb: sbs) if (sb!=null) beds.add(sb);
			}
			if (beds.size()==0) return;
			
			SimpleBed[] bedRegions = new SimpleBed[beds.size()];
			beds.toArray(bedRegions);
			Arrays.sort(bedRegions);
			
			
			File bedFile = new File (saveDirectory, workingReportName+".fusions.bed.gz");
			out = new Gzipper(bedFile);

			//add header
			out.println("#ChrHg38\tGeneStart\tGeneStop\tGeneName1:GeneName2:Type\tScore\tGeneStrand");

			//for each bed
			for (SimpleBed fus: bedRegions) out.println(fus);

			out.close();
		} catch (IOException e) {
			if (out != null) {
				out.closeNoException();
				out.getGzipFile().deleteOnExit();
			}
			e.printStackTrace();
			Misc.printErrAndExit("\nError: issue writing out vcf for "+workingXmlFile);
		}
	}


	private void loadXml() {
		try {
			//Load and Parse the XML document
			Document document = builder.parse(workingXmlFile);
			//Iterating through the nodes and extract the data
			NodeList nodeList = document.getDocumentElement().getChildNodes();
			for (int i = 0; i < nodeList.getLength(); i++) {
				Node node = nodeList.item(i);
				if (node instanceof Element) {
					String nodeName = node.getNodeName();
					if (nodeName.equals("testDetails")) parseTestDetails(node);
					else if (nodeName.equals("patientInformation")) parsePatientInformation(node);
					else if (nodeName.equals("physicianInformation")) parsePhysicianInformation(node);
					else if (nodeName.equals("specimenInformation")) parseSpecimenInformation(node);
					else if (nodeName.equals("tests")) parseTests(node);
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing this xml file "+workingXmlFile);
		}
	}

	/**Parses the testDetails section
	 * @throws ParseException 
	 * @throws DOMException 
	 * @throws IOException */
	private void parseTests(Node node) throws DOMException, ParseException, IOException {
		//parse the tests
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//testName
				if (name.equals("testName")) {
					String testName = cNode.getTextContent().trim();
					if (testName.equals("Exome Panel - Clinical Genes")) parseGenesTest(childNodes);
					else if (testName.equals("Exome Panel - Additional Genes")) parseGenesTest(childNodes);
					else if (testName.equals("Exome CND Panel - Clinical Genes")) parseCNVPanel(childNodes);
					else if (testName.equals("Exome CNA Panel - Clinical Genes")) parseCNVPanel(childNodes);
					else if (testName.equals("Exome CNA Panel - Additional Genes")) parseCNVPanel(childNodes);
					else if (testName.equals("Transcriptome Detection_v1 Panel")) parseTranslocationPanel(childNodes);
					else if (testName.equals("Transcriptome Detection_v1 Variant Panel")) parseTranslocationPanel(childNodes);
					else if (testName.equals("MGMT-Me")) parseMethyl(childNodes);
					else if (ihcTestNames.contains(testName)) parseIHC(childNodes);
					else if (testName.equals("Her2 CISH")) parseCNVPanel(childNodes);
					else if (testName.equals("EBER ISH")){
						IO.el("\nSkipping the parsing of 'EBER ISH' test for Epstein-Barr virus genome copy number\n");
					}
					else if (testName.equals("ROS1 FISH")){
						IO.el("\nSkipping the parsing of 'ROS1 FISH' test for ROS1 cytogenic alteration.");
					}
					else if (testName.equals("RNA Expression")){
						IO.el("\nSkipping the parsing of 'RNA Expression' test for a small set of genes.");
					}
					else if (testName.equals("eKarotypeGraph")){
						//silently skip, this is a binary result
					}
					else throw new IOException("Found an unknown test! "+testName);
				}
			}
		}
	}

	private void parseIHC(NodeList childNodes) throws IOException {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//find testResults
				if (name.equals("testResults")) {
					NodeList subNodes = cNode.getChildNodes();
					for (int i=0; i< subNodes.getLength(); i++) {
						Node n = subNodes.item(i);
						if (n instanceof Element) {
							String subName = n.getNodeName();
							if (subName.equals("expressionAlteration")) this.workingExpressionAlterations.add( new ExpressionAlteration(this.workingReportAttributes, n.getChildNodes()));
							else throw new IOException("Found something other than 'expressionAlteration' "+subName);
						}
					}
				}
			}
		}
	}
	
	private void parseMethyl(NodeList childNodes) throws IOException {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//find testResults
				if (name.equals("testResults")) {
					NodeList subNodes = cNode.getChildNodes();
					for (int i=0; i< subNodes.getLength(); i++) {
						Node n = subNodes.item(i);
						if (n instanceof Element) {
							String subName = n.getNodeName();
							if (subName.equals("epigeneticAlteration")) workingMethylationAlterations.add( new MethylationAlteration(workingReportAttributes, n.getChildNodes()));
							else throw new IOException("Found something other than 'epigeneticAlteration' "+subName);
						}
					}
				}
			}
		}
	}

	private void parseTranslocationPanel(NodeList childNodes) throws IOException {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//find testResults
				if (name.equals("testResults")) {
					NodeList subNodes = cNode.getChildNodes();
					for (int i=0; i< subNodes.getLength(); i++) {
						Node n = subNodes.item(i);
						if (n instanceof Element) {
							String subName = n.getNodeName();
							if (subName.equals("translocation")) workingTranslocations.add(new Translocation(n.getChildNodes()));
							else throw new IOException("Found something other than 'translocation' "+subName);
						}
					}
				}
			}
		}
	}

	private void parseCNVPanel(NodeList childNodes) throws IOException {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//find testResults
				if (name.equals("testResults")) {
					NodeList subNodes = cNode.getChildNodes();
					for (int i=0; i< subNodes.getLength(); i++) {
						Node n = subNodes.item(i);
						if (n instanceof Element) {
							String subName = n.getNodeName();
							if (subName.equals("copyNumberAlteration")) workingCNVAlterations.add(new CNVAlteration(n.getChildNodes()));
							else throw new IOException("Found something other than 'copyNumberAlteration' "+subName);
						}
					}
				}
			}
		}
	}

	//Exome Panel - Clinical Genes
	private void parseGenesTest(NodeList childNodes) throws IOException {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//find testResults
				if (name.equals("testResults")) {
					NodeList subNodes = cNode.getChildNodes();
					for (int i=0; i< subNodes.getLength(); i++) {
						Node n = subNodes.item(i);
						if (n instanceof Element) {
							String subName = n.getNodeName();
							if (subName.equals("tumorMutationBurden")) parseTMB(n.getChildNodes());
							else if (subName.equals("microsatelliteInstability")) parseMSI(n.getChildNodes());
							else if (subName.equals("genomicLevelHeterozygosity")) parseLOH(n.getChildNodes());
							else if (subName.equals("genomicAlteration")) workingGenomicAlterations.add(new GenomicAlteration(n.getChildNodes()));
						}
					}
				}

			}
		}

	}

	private void parseLOH(NodeList childNodes) {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node n = childNodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("result")) addLastChild("result", n, "LOHCall");
				else if (name.equals("LOHpercentage")) addLastChild("LOHpercentage", n, "LOHScore");
			}
		}
	}

	private void parseMSI(NodeList childNodes) {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node n = childNodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("msiCall")) {
					addLastChild("msiCall", n, "MSICall");
					return;
				}
			}
		}
	}

	private void parseTMB(NodeList childNodes) {
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node n = childNodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("mutationBurdenCall")) addLastChild("mutationBurdenCall", n, "TMBCall");
				else if (name.equals("mutationBurdenScore")) addLastChild("mutationBurdenScore", n, "TMBScore");
			}
		}
	}

	/**Parses the specimenInformation section
	 * @throws ParseException 
	 * @throws DOMException 
	 * @throws IOException */
	private void parseSpecimenInformation(Node node) throws DOMException, ParseException, IOException {
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//right now there's only one
				if (name.equals("tumorSpecimenInformation")){
					NodeList ccNodes = cNode.getChildNodes();
					for (int x=0; x< ccNodes.getLength(); x++) {
						Node n = ccNodes.item(x);
						if (n instanceof Element) {
							String nm = n.getNodeName();
							//specimenID
							if (nm.equals("specimenID")) addLastChild("specimenID",n);
							//specimenType
							else if (nm.equals("specimenType")) addLastChild("specimenType",n);
							//specimenAccessionID
							else if (nm.equals("specimenAccessionID")) addLastChild("specimenAccessionID",n);
							//specimenSite
							else if (nm.equals("specimenSite")) addLastChild("specimenSite",n);
							//specimenCollectionDate
							else if (nm.equals("specimenCollectionDate")) addLastChild("specimenCollectionDate",n);
						}
					}
				}
				else throw new IOException("Found something other than tumorSpecimenInformation! "+name);
			}
		}
	}

	/**Parses the testDetails section
	 * @throws ParseException 
	 * @throws DOMException */
	private void parseTestDetails(Node node) throws DOMException, ParseException {
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//testName
				if (name.equals("testName")) addLastChild("testName",cNode);
				//testCode
				else if (name.equals("testCode")) addLastChild("testCode",cNode);
				//labReportID
				else if (name.equals("labReportID")) addLastChild("labReportID",cNode);
				//orderedDate
				else if (name.equals("orderedDate")) addLastChild("orderedDate",cNode);

			}
		}
	}

	/**Parses the patientInformation section
	 * @throws ParseException 
	 * @throws DOMException */
	private void parsePatientInformation(Node node) throws DOMException, ParseException {
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				if (includePHI) {
					//lastName
					if (name.equals("lastName")) addLastChild("lastName",cNode);
					//firstName
					else if (name.equals("firstName")) addLastChild("firstName",cNode);
					//dob
					else if (name.equals("dob")) addLastChild("dob",cNode);
					//mrn
					else if (name.equals("mrn")) addLastChild("mrn",cNode);
				}
				//gender
				if (name.equals("gender")) addLastChild("gender",cNode);
				//icd_code
				else if (name.equals("icd_code")) addLastChild("icd_code",cNode);
				//diagnosis
				else if (name.equals("diagnosis")) addLastChild("diagnosis",cNode);
				//pathologicDiagnosis
				else if (name.equals("pathologicDiagnosis")) addLastChild("pathologicDiagnosis",cNode);
				//primarySite
				else if (name.equals("primarySite")) addLastChild("primarySite",cNode);
				//lineage
				else if (name.equals("lineage")) addLastChild("lineage",cNode);
				//subLineage
				else if (name.equals("subLineage")) addLastChild("subLineage",cNode);
				//zipcode
				else if (name.equals("zipcode")) addLastChild("zipcode",cNode);
			}
		}
	}

	/**Parses the physicianInformation section
	 * @throws ParseException 
	 * @throws DOMException */
	private void parsePhysicianInformation(Node node) throws DOMException, ParseException {
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//testName
				if (name.equals("fullName")) addLastChild("fullName",cNode, "physicianName");
			}
		}
	}

	/**Adds the last child to the reportAttributes*/
	private void addLastChild(String key, Node node){
		Node last = node.getLastChild();
		if (last != null){
			String x = last.getTextContent().trim();
			//x = x.replaceAll("&#13", " ");
			x = Misc.WHITESPACE.matcher(x).replaceAll(" ");
			if (x.length()!=0) workingReportAttributes.put(key, x);
		}
	}

	/**Adds the last child to the reportAttributes*/
	public void addLastChild(String key, Node node, String newKey){
		Node last = node.getLastChild();
		if (last != null){
			String x = last.getTextContent().trim();
			//x = x.replaceAll("&#13", " ");
			x = Misc.WHITESPACE.matcher(x).replaceAll(" ");
			if (x.length()!=0) workingReportAttributes.put(newKey, x);
		}
	}

	public static String fetch(String key, Node node){
		Node last = node.getLastChild();
		if (last != null){
			String x = last.getTextContent().trim();
			x = Misc.RETURN.matcher(x).replaceAll("");
			if (x.length()!=0) return x;
		}
		return null;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CarisXmlVcfParser(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		source = useqVersion+" Args: "+ Misc.stringArrayToString(args, " ");
		System.out.println("\n"+ source +"\n");
		File ucscFile = null;
		File vcfXmlDir = null;
		
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': vcfXmlDir = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'u': ucscFile = new File(args[++i]); break;
					case 'i': includePHI = true; break;
					case 'k': skipUnPaired = true; break;
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

		if (vcfXmlDir != null && vcfXmlDir.exists()) parseXmlVcfFiles(vcfXmlDir);
		else Misc.printErrAndExit("\nError: cannot find the directory -d containing your xml and vcf files?! ");

		if (ucscFile == null || ucscFile.exists() == false) Misc.printErrAndExit("\nError: cannot find or read your hg38 merged UCSC refFlat gene table? "+ucscFile);


		//main transcripts
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(ucscFile, 0);
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your transcript's coordinates are reversed. Check that each start is less than the stop.\n");
		IO.pl("Loading gene models...");
		name2GeneModels = reader.getGeneNameTranscripts();



	}	

	private void parseXmlVcfFiles(File vcfXmlDir) {
		HashMap<String, File> vcfs = new HashMap<String, File>();
		HashMap<String, File> xmls = new HashMap<String, File>();
		File[] files = IO.extractFiles(vcfXmlDir);
		for (int i=0; i< files.length; i++) {
			String fullName = files[i].getName();
			if (fullName.endsWith(".gz"))fullName = fullName.substring(0, fullName.length()-3);
			String[] nameParts = Misc.UNDERSCORE.split(fullName);
			if (fullName.endsWith(".vcf")) {
				boolean itIsACaris = carisVcf(files[i]);
				if (itIsACaris) vcfs.put(nameParts[0], files[i]);
			}
			if (fullName.endsWith(".xml")) xmls.put(nameParts[0], files[i]);
		}
		//same number?
		if (vcfs.size()==0 || (vcfs.size()!=xmls.size())) {
			if (skipUnPaired == false) Misc.printErrAndExit("\nERROR: the # vcf and xml files differ in "+vcfXmlDir);
		}
		xmlVcfFiles = new HashMap<String, File[]>();
		for (String name: xmls.keySet()) {
			File xml = xmls.get(name);
			File vcf = vcfs.get(name);
			if (vcf == null) {
				if (skipUnPaired == false) Misc.printErrAndExit("\nERROR: failed to find a matching vcf file for "+xml);
			}
			else xmlVcfFiles.put(name, new File[] {xml,vcf});
		}
	}

	private boolean carisVcf(File file) {
		try {
			BufferedReader in = IO.fetchBufferedReader(file);
			String line = null;
			while ((line = in.readLine())!=null) {
				line = line.trim();
				if (line.startsWith("##INFO=<ID=CI,Number=.,Type=String,Description=")) {
					in.close();
					return true;
				}
				if (line.startsWith("#CHROM")) {
					in.close();
					return false;
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return false;
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Caris Xml Vcf Parser: May 2022                         **\n" +
				"**************************************************************************************\n" +
				"This tool parses Caris paired xml and vcf report files to generate: new vcfs where xml\n"+
				"reported genomic alternations are annotated, bed files of copy number changes and gene\n"+
				"fusions as well as a summary spreadsheet of the non NGS test results and report\n"+
				"attributes. Use Vt to normalize the vcf output files. \n"+

				"\nOptions:\n"+
				"-d Path to a directory containing paired xml and vcf files from Caris.\n"+
				"-s Path to a directory for saving the results. Only impacting vcf, cnv, fusions are\n"+
				"     saved, if none present, no file(s) are created.\n"+
				"-u Path to a Hg38 UCSC RefFlat or RefSeq merged gene file for CNV and gene fusion\n"+
				"   coordinate extraction. See: http://www.genome.ucsc.edu/FAQ/FAQformat.html#format9\n"+
				"   (refSeqGeneName name2(optional) chr strand txStart txEnd cdsStart cdsEnd exonCount\n"+
				"   exonStarts exonEnds). Be sure to merge transcripts with the USeq MergeUCSCGeneTable\n"+
				"   app and remove non standard chromosomes.\n"+
				"-i Include PHI in spreadsheet output, defaults to excluding.\n"+
				"-a Include all vcf records in output, defaults to just those with an xml match.\n"+
				"-k Skip datasets that are missing either an xml or vcf report.\n"+

				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/CarisXmlVcfParser -d CarisReports/\n" +
				"     -s ParsedCarisReports/ -u ~/GRCh38/hg38RefSeq13Dec2020_MergedStdChr.ucsc.gz \n\n" +

				"**************************************************************************************\n");
	}

}
