package edu.utah.seq.run.caris;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import util.gen.IO;
import util.gen.Misc;

public class CarisXmlParser {
	
	private File xmlFile = null;
	private File deidentifiedXmlFile = null;
	
	//patient info PHI
	private String firstName = null;
	private String lastName = null;
	private String dob = null; 	//1965-06-15, 1950-07-26, don't use this use dobMonth, dobDay, dobYear
	private int dobMonth = -1;
	private int dobDay = -1;
	private int dobYear = -1;
	private String gender = null; //M or F
	private String mrn = null;
	
	//physician info
	private String firstNamePhysician = null;
	private String lastNamePhysician = null;

	
	public CarisXmlParser(File xmlFile) throws Exception {
		this.xmlFile = xmlFile;
		loadPatientInfo();
		//printInfo();
	}
	
	/**lastName firstName dobMonth(1-12) dobDay(1-31) dobYear(1900-2050) gender(M|F) mrn */
	public String fetchSubjectMatchMakerLine() {
		StringBuilder sb = new StringBuilder();
		if (lastName != null) sb.append(lastName);
		else sb.append(".");
		sb.append("\t");
		
		if (firstName != null) sb.append(firstName);
		else sb.append(".");
		sb.append("\t");
		
		if (dobMonth != -1) sb.append(new Integer(dobMonth).toString());
		else sb.append(".");
		sb.append("\t");
		
		if (dobDay != -1) sb.append(new Integer(dobDay).toString());
		else sb.append(".");
		sb.append("\t");
		
		if (dobYear != -1) sb.append(new Integer(dobYear).toString());
		else sb.append(".");
		sb.append("\t");
		
		if (gender != null) sb.append(gender);
		else sb.append(".");
		sb.append("\t");
		
		if (mrn != null) sb.append(mrn);
		else sb.append(".");
		return sb.toString();
	}
	
	private void printInfo() {
		IO.pl("\n"+ xmlFile.getName()+" -> "+deidentifiedXmlFile.getName());
		IO.pl(lastName+", "+firstName);
		IO.pl(dob+" -> "+dobYear+" "+dobMonth+" "+dobDay);
		IO.pl(gender+" "+mrn);
		
	}

	private void loadPatientInfo() throws Exception {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document document = builder.parse(xmlFile);
		NodeList nodeList = document.getDocumentElement().getChildNodes();
		for (int i = 0; i < nodeList.getLength(); i++) {
			Node node = nodeList.item(i);
			if (node instanceof Element) {
				String nodeName = node.getNodeName();
				if (nodeName.equals("patientInformation")) {
//IO.pl("Parsing patient info...");
					parseAndCleanPatientInformation(node);
				}
				else if (nodeName.equals("physicianInformation")) {
//IO.pl("Parsing physician info...");
					parsePhysicianInformation(node);
				}
			}
		}
//IO.pl("Printing deidentified doc...");
		printDocument(document);
		
	}
	
	public void printDocument(Document doc) throws Exception {
		String name =  Misc.removeExtension(xmlFile.getCanonicalPath());
		if (firstNamePhysician == null || lastNamePhysician == null) deidentifiedXmlFile = new File(name+"_deid.xml");
		else deidentifiedXmlFile = new File(name+"_deid_"+firstNamePhysician+"_"+lastNamePhysician+".xml");
		
		FileOutputStream out = new FileOutputStream(deidentifiedXmlFile);
	    TransformerFactory tf = TransformerFactory.newInstance();
	    Transformer transformer = tf.newTransformer();
	    transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "no");
	    transformer.setOutputProperty(OutputKeys.METHOD, "xml");
	    transformer.setOutputProperty(OutputKeys.INDENT, "yes");
	    transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
	    transformer.transform(new DOMSource(doc), new StreamResult(new OutputStreamWriter(out, "UTF-8")));
	    out.close();
	}
	
	private boolean canceledTests(File xml) throws IOException {
		BufferedReader in = IO.fetchBufferedReader(xml);
		String line = null;
		while ((line = in.readLine())!=null) {
			if (line.contains("<test_cancellation_reason>")) return true;
		}
		return false;
	}

	private void parseAndCleanPatientInformation(Node node) throws Exception {
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//lastName
				if (name.equals("lastName")) lastName = getLastChildAndDelete("lastName",cNode);
				//firstName
				else if (name.equals("firstName")) firstName = getLastChildAndDelete("firstName",cNode);
				//dob
				else if (name.equals("dob")) {
					dob = getLastChildAndDelete("dob",cNode);
					parseDoB();
				}
				//mrn
				else if (name.equals("mrn")) mrn = getLastChildAndDelete("mrn",cNode);
				//gender, don't delete, needed for SampleConcordance
				else if (name.equals("gender")) {
					gender = getLastChild("gender",cNode);
					parseGender();
				}
				//fullName
				else if (name.equals("fullName")) cNode.getParentNode().removeChild(cNode);
				//middleName
				else if (name.equals("middleName")) cNode.getParentNode().removeChild(cNode);
				//contactInformation
				else if (name.equals("contactInformation")) cNode.getParentNode().removeChild(cNode);
				//middleName
				else if (name.equals("addressInformation")) cNode.getParentNode().removeChild(cNode);
			}
		}
	}
	
	private void parseGender() throws Exception {
		if (gender == null) return;
		//Male or Female
		if (gender.equals("Male")) gender = "M";
		else if (gender.equals("Female")) gender = "F";
		else throw new Exception("ERROR: parsing gender, must be Male or Female, see <gender> in "+xmlFile.getCanonicalPath());
	}

	private void parseDoB() throws IOException, Exception {
		if (dob == null || dob.length()==0) return;
		//1965-06-15, 1950-07-26
		String[] t = Misc.DASH.split(dob);
		if (t.length!=3) throw new Exception("ERROR: parsing date of birth, not 3 fields after splitting on dash, see <dob> in "+xmlFile.getCanonicalPath());
		if (t[0].length()!=0) {
			dobYear = Integer.parseInt(t[0]);
			if (dobYear< 1900 || dobYear > 2050) throw new IOException("ERROR: dob year field '"+t[0]+"' is malformed, must be 1900-2050, see <dob> in "+xmlFile.getCanonicalPath());
		}
		
		if (t[1].length()!=0) {
			dobMonth = Integer.parseInt(t[1]);
			if (dobMonth< 1 || dobMonth > 12) throw new IOException("ERROR: dob month field '"+t[1]+"' is malformed, must be 1-12, see <dob> in "+xmlFile.getCanonicalPath());
		}
		
		if (t[2].length()!=0) {
			dobDay = Integer.parseInt(t[2]);
			if (dobDay< 1 || dobDay > 31) throw new IOException("ERROR: dob day field '"+t[2]+"' is malformed, must be 1-31, see <dob> in "+xmlFile.getCanonicalPath());
		}
	}

	private void parsePhysicianInformation(Node node) {
		//parse its children
		NodeList childNodes = node.getChildNodes();
		for (int j = 0; j < childNodes.getLength(); j++) {
			Node cNode = childNodes.item(j);
			if (cNode instanceof Element) {
				String name = cNode.getNodeName();
				//remove any ' in the name since this will be used in the file name
				//lastName
				if (name.equals("lastName")) {
					lastNamePhysician = getLastChild("lastName",cNode);
					lastNamePhysician = Misc.SINGLE_QUOTE.matcher(lastNamePhysician).replaceAll("");
					
				}
				//firstName
				else if (name.equals("firstName")) {
					firstNamePhysician = getLastChild("firstName",cNode);
					firstNamePhysician = Misc.SINGLE_QUOTE.matcher(firstNamePhysician).replaceAll("");
				}
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
	
	/**Returns the last child if it exits*/
	private String getLastChildAndDelete(String key, Node node){
		Node last = node.getLastChild();
		if (last != null){
			String x = last.getTextContent().trim();
			node.getParentNode().removeChild(node);
			return x;
		}
		return null;
	}

	/*For testing **/
	public static void main(String[] args) throws Exception {
		//File x = new File("/Users/u0028003/Downloads/Caris/TN22-126203_2022-04-11_11_30.xml");
		//new CarisXmlParser(x);
		File dir = new File("/Users/u0028003/Downloads/Caris/Delme_PHI");
		File[] xmls = IO.extractFiles(dir, ".xml");
		for (File x : xmls) new CarisXmlParser(x);
	}

	public File getXmlFile() {
		return xmlFile;
	}
	public File getDeidentifiedXmlFile() {
		return deidentifiedXmlFile;
	}
	public String getFirstName() {
		return firstName;
	}
	public String getLastName() {
		return lastName;
	}
	public int getDobMonth() {
		return dobMonth;
	}
	public int getDobDay() {
		return dobDay;
	}
	public int getDobYear() {
		return dobYear;
	}
	/**M or F*/
	public String getGender() {
		return gender;
	}
	public String getMrn() {
		return mrn;
	}
	public String getFirstNamePhysician() {
		return firstNamePhysician;
	}
	public String getLastNamePhysician() {
		return lastNamePhysician;
	}

}
