package edu.utah.seq.run.caris;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.TimeUnit;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import util.gen.IO;
import util.gen.Misc;

public class CarisPatient {
	
	//fields
	private String testID = null;
	private String hciID = null;
	private CarisDataWrangler cdw = null;
	private ArrayList<String[]> objectInfo = new ArrayList<String[]>();
	private ArrayList<String> vcfNames = new ArrayList<String>();
	private ArrayList<String> xmlNames = new ArrayList<String>();
	private ArrayList<String> pdfNames = new ArrayList<String>();
	private ArrayList<String> dnaNames = new ArrayList<String>();
	private ArrayList<String> rnaNames = new ArrayList<String>();
	
	//for retrieving the hci patient id
	private String firstName = null;
	private String lastName = null;
	private String dob = null;
	private String gender = null;
	private String mrn = null;
	
	private File testDir = null;
	private File clinReportDir = null;
	private boolean ready = true;
	
	public static final SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	public CarisPatient(String testID, CarisDataWrangler cdw) {
		this.testID = testID;
		this.cdw = cdw;
	}
	
	public void addObjectLine(String[] tokens) {
		objectInfo.add(tokens);
	}

	public void parseFileLines(int minHours) throws Exception {
		if (checkAges(minHours) == false) return;
		if (loadArrayLists() == false) return;
		checkDatasets();
	}
	
	private void checkDatasets() {
		boolean oneVcf = vcfNames.size() == 1;
		boolean oneXml = xmlNames.size() == 1;
		boolean twoDna = dnaNames.size() == 2;
		boolean oneRna = rnaNames.size() == 1;
		//must have a vcf, xml, and two dna fastq files; if rna present, must have two files
		if (oneVcf == false || oneXml == false || twoDna == false || oneRna == true) ready = false;
		IO.pl("\tXml\t"+oneXml+"\t"+ xmlNames);
		IO.pl("\tVcf\t"+oneVcf+"\t"+ vcfNames);
		IO.pl("\tDNA\t"+twoDna+"\t"+ dnaNames);
		IO.pl("\tRNA\t"+(oneRna==false) +"\t"+ rnaNames);
		IO.pl("\tOK\t"+ready);
	}

	private boolean loadArrayLists() {
		//2022-03-12 11:33:49 2077226966 RNA_TN22-109832_S29_L001_R2_001.fastq.gz
		//    0         1          2                         3
		for (String[] tokens : objectInfo) {
			String fileName = tokens[3];
			if (fileName.endsWith(".fastq.gz")) {
				if (fileName.startsWith("RNA_")) rnaNames.add(fileName);
				else if (fileName.startsWith("DNA_")) dnaNames.add(fileName);
				else {
					cdw.getErrorMessages().add("Couldn't source the fastq type from "+tokens[3]);
					ready = false;
					return false;
				}
			}
			else if (fileName.endsWith(".xml")) xmlNames.add(fileName);
			else if (fileName.endsWith(".vcf")) vcfNames.add(fileName);
			else if (fileName.endsWith(".pdf")) pdfNames.add(fileName);
			else {
				cdw.getErrorMessages().add("Unknown type type "+tokens[3]);
				ready = false;
				return false;
			}
		}
		return true;
	}

	public boolean checkAges(int minHours) throws ParseException {
		//check ages, must all be > min hours
		//2022-03-12 11:33:49 2077226966 RNA_TN22-109832_S29_L001_R2_001.fastq.gz
		//    0         1          2                         3
		for (String[] tokens : objectInfo) {
			Date objectDate = sdf.parse(tokens[0]+" "+tokens[1]);
			long diff = System.currentTimeMillis() - objectDate.getTime();
			long diffHours = TimeUnit.HOURS.convert(diff, TimeUnit.MILLISECONDS);
			if (diffHours< minHours) {
				IO.pl("\tTooYoung\t"+tokens[3]);
				ready = false;
				return false;
			}
		}
		return true;
	}

	public boolean isReady() {
		return ready;
	}

	public void downloadDatasets() throws Exception {
		File fastqDir = new File (testDir, "Fastq");
		fastqDir.mkdir();

		//cp in vcf
		cdw.cp(vcfNames.get(0), new File (clinReportDir, vcfNames.get(0)));
		//cp in xml, already done so skip
		//cdw.cp(xmlNames.get(0), new File (clinDir, xmlNames.get(0)));
		//cp in DNA0
		cdw.cp(dnaNames.get(0), new File (fastqDir, "/TumorDNA/"+dnaNames.get(0)));
		//cp in DNA1
		cdw.cp(dnaNames.get(1), new File (fastqDir, "/TumorDNA/"+dnaNames.get(1)));
		//any RNA?
		if (rnaNames.size() == 2) {
			//cp in RNA0
			cdw.cp(rnaNames.get(0), new File (fastqDir, "/TumorRNA/"+rnaNames.get(0)));
			//cp in RNA1
			cdw.cp(rnaNames.get(1), new File (fastqDir, "/TumorRNA/"+rnaNames.get(1)));
		}

		
		
	}

	public void fetchHCIId() throws Exception {

		//cp in xml
		File xml = new File (cdw.getTmpDir(), xmlNames.get(0));
		cdw.cp(xmlNames.get(0), xml);
		
		//load patient info
		loadPatientInfo(xml);
		
//TODO: fetch HCI patient ID given the available PHI
hciID = "PH_"+Misc.getRandomString(10);

		//make the patient dir
		testDir = new File (cdw.getJobsDirectory(), hciID+"/Caris/"+testID+"/");
		testDir.mkdirs();
		clinReportDir = new File (testDir, "/ClinicalReport/");
		clinReportDir.mkdirs();
		
		//move the report into the ClinicalReport folder
		xml.renameTo(new File(clinReportDir, xml.getName()));

	}

	private void loadPatientInfo(File xml) throws Exception {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document document = builder.parse(xml);
		NodeList nodeList = document.getDocumentElement().getChildNodes();
		for (int i = 0; i < nodeList.getLength(); i++) {
			Node node = nodeList.item(i);
			if (node instanceof Element) {
				String nodeName = node.getNodeName();
				if (nodeName.equals("patientInformation")) {
					parsePatientInformation(node);
					break;
				}
			}
		}
	}

	private void parsePatientInformation(Node node) {

		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//lastName
				if (name.equals("lastName")) lastName = getLastChild("lastName",cNode);
				//firstName
				else if (name.equals("firstName")) firstName = getLastChild("firstName",cNode);
				//dob
				else if (name.equals("dob")) dob = getLastChild("dob",cNode);
				//mrn
				else if (name.equals("mrn")) mrn = getLastChild("mrn",cNode);
				//gender
				if (name.equals("gender")) gender = getLastChild("gender",cNode);
			}
		}
	}

	/**Returns the last child if it exits*/
	private String getLastChild(String key, Node node){
		Node last = node.getLastChild();
		if (last != null){
			String x = last.getTextContent().trim();
			return x;
		}
		return null;
	}

}
