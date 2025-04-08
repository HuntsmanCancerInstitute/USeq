package edu.utah.seq.run.tempus;

import java.io.File;
import java.io.IOException;
import org.json.JSONObject;
import util.gen.IO;
import util.gen.Misc;

public class TempusJsonParser {
	
	private File jsonFile = null;
	private String s3Path = null;
	private JSONObject mainJsonObject = null;
	private File deidentifiedJsonFile = null;
	
	//patient info PHI
	private String firstName = null;
	private String lastName = null;
	private String dob = null; 	//1965-06-15, 1950-07-26, don't use this use dobMonth, dobDay, dobYear
	private int dobMonth = -1;
	private int dobDay = -1;
	private int dobYear = -1;
	private String gender = null; //M or F
	private String mrn = null;
	
	//order info
	private String firstNamePhysician = null;
	private String lastNamePhysician = null;
	private String testId = null;
	private String testCode = null;
	private String testDate = null;
	private boolean parsed = false;
	private boolean isDNATest = false;

	
	public TempusJsonParser(File jsonFile, String s3Path) throws Exception {
		this.jsonFile = jsonFile;
		this.s3Path = s3Path;
		loadJsonInfo();
		checkTestOrder();
		//printInfo();
	}
	
	/*
	 	XE.V1 Complete exome (germline, tumor) + mRNA (tumor)
		XO Panel of 1714 genes (germline, tumor) + mRNA (tumor)
		XT.V1 Panel of 595 genes (germline, tumor)
		XT.V1 Panel of 595 genes (germline, tumor) + mRNA (tumor)
		XT.V1 Panel of 595 genes (tumor) + mRNA (tumor)
		XT.V2 Panel of 596 genes (germline, tumor) + mRNA (tumor)
		XT.V2 Panel of 596 genes (tumor) + mRNA (tumor)
		XT.V3 Panel of 648 genes (germline, tumor) + mRNA (tumor)
		XT.V3 Panel of 648 genes (tumor) + mRNA (tumor)
		XT.V4 Panel of 648 genes (germline, tumor) + mRNA (tumor)
		
		HRD - not supported
		MMR - not supported
		PD-L1 - not supported
		RS.v2 - not supported
		XE
		XF - not supported
		XO
		XT
	 */
	private void checkTestOrder() {
		if (parsed == false) return;
		if (testCode.startsWith("XE") || testCode.startsWith("XO") || testCode.startsWith("XT")) isDNATest = true;
		else isDNATest = false;
	}

	/**lastName firstName dobMonth(1-12) dobDay(1-31) dobYear(1900-2050) gender(M|F) mrn */
	public String fetchSubjectMatchMakerLine() {
		StringBuilder sb = new StringBuilder();
		if (lastName != null) sb.append(lastName);
		sb.append("\t");
		
		if (firstName != null) sb.append(firstName);
		sb.append("\t");
		
		if (dobMonth != -1) sb.append(new Integer(dobMonth).toString());
		sb.append("\t");
		
		if (dobDay != -1) sb.append(new Integer(dobDay).toString());
		sb.append("\t");
		
		if (dobYear != -1) sb.append(new Integer(dobYear).toString());
		sb.append("\t");
		
		if (gender != null) sb.append(gender);
		sb.append("\t");
		
		if (mrn != null) sb.append(mrn);
		return sb.toString();
	}
	
	private void printInfo() {
		IO.pl("\n"+ jsonFile.getName()+" -> "+deidentifiedJsonFile.getName());
		IO.pl(lastName+", "+firstName);
		IO.pl(dob+" -> "+dobYear+" "+dobMonth+" "+dobDay);
		IO.pl(gender+" "+mrn);
		
	}

	private void loadJsonInfo() {
		try {
			String jString = IO.loadFile(jsonFile, " ", true);
			mainJsonObject = new JSONObject(jString);

			//report
			TempusJsonReport report = new TempusJsonReport(mainJsonObject);
			//2018-10-12T02:50:44+00:00
			String sd = report.getSignout_date();
			testDate = sd.substring(0, sd.indexOf('T')); 

			//patient
			TempusJsonPatient patient = new TempusJsonPatient(mainJsonObject);
			firstName = patient.getFirstName();
			lastName = patient.getLastName();
			gender = patient.getSex();
			dob = patient.getDateOfBirth();
			mrn = patient.getEmr_id();
			parseGender();
			parseDoB();

			//order
			TempusJsonOrder order = new TempusJsonOrder(mainJsonObject);
			testId = order.getAccessionId();
			testCode = order.getTestCode();
			String firstLast = order.getPhysician();
			firstLast = Misc.SINGLE_QUOTE.matcher(firstLast).replaceAll("");
			String[] split = Misc.WHITESPACE.split(firstLast);
			firstNamePhysician = split[0];
			StringBuilder sb = new StringBuilder (split[1]);
			for (int i=2; i< split.length; i++) {
				sb.append("-");
				sb.append(split[i]);
			}
			lastNamePhysician = sb.toString();
			parsed = true;

		} catch (Exception e) {
				IO.el("\tBroken json File! Skipping -> "+jsonFile+" "+e.getMessage());
				parsed = false;
			}
		}
	
	public void printDeIDDocument() throws Exception {
		//testId_testCode_deid_firstNamePhy_lastNamePhy_gender.json
		String name = testId+"_"+testCode+"_"+testDate;
		if (firstNamePhysician == null || lastNamePhysician == null) deidentifiedJsonFile = new File(jsonFile.getParentFile(), name+"_deid.json");
		else deidentifiedJsonFile = new File(jsonFile.getParentFile(), name+"_deid_"+firstNamePhysician+"_"+lastNamePhysician+"_"+gender+".json");
		IO.writeString(mainJsonObject.toString(5), deidentifiedJsonFile);
	}
	
	private void parseGender() throws Exception {
		if (gender == null) return;
		//Male or Female
		if (gender.equals("MALE")) gender = "M";
		else if (gender.equals("FEMALE")) gender = "F";
		else {
			gender = "NA";
			throw new Exception("ERROR: parsing gender, must be Male or Female, see <sex> in "+jsonFile.getCanonicalPath());
		}
	}

	private void parseDoB() throws Exception {
		if (dob == null || dob.length()==0) return;
		//1965-06-15, 1950-07-26
		String[] t = Misc.DASH.split(dob);
		if (t.length!=3) throw new Exception("ERROR: parsing date of birth, not 3 fields after splitting on dash, see <dob> in "+jsonFile.getCanonicalPath());
		if (t[0].length()!=0) {
			dobYear = Integer.parseInt(t[0]);
			if (dobYear< 1900 || dobYear > 2050) throw new IOException("ERROR: dob year field '"+t[0]+"' is malformed, must be 1900-2050, see <dob> in "+jsonFile.getCanonicalPath());
		}
		
		if (t[1].length()!=0) {
			dobMonth = Integer.parseInt(t[1]);
			if (dobMonth< 1 || dobMonth > 12) throw new IOException("ERROR: dob month field '"+t[1]+"' is malformed, must be 1-12, see <dob> in "+jsonFile.getCanonicalPath());
		}
		
		if (t[2].length()!=0) {
			dobDay = Integer.parseInt(t[2]);
			if (dobDay< 1 || dobDay > 31) throw new IOException("ERROR: dob day field '"+t[2]+"' is malformed, must be 1-31, see <dob> in "+jsonFile.getCanonicalPath());
		}
	}


	/*For testing **/
	public static void main(String[] args) throws Exception {
		File dir = new File("/Users/u0028003/Downloads/Tempus/");
		File[] jsons = IO.extractFiles(dir, ".json");
		for (File x : jsons) new TempusJsonParser(x, "MockS3Path/ddd/ddd");
	}

	public File getJsonFile() {
		return jsonFile;
	}
	public File getDeidentifiedJsonFile() {
		return deidentifiedJsonFile;
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

	public String getTestCode() {
		return testCode;
	}

	public String getTestId() {
		return testId;
	}

	public String getTestDate() {
		return testDate;
	}

	public String getS3Path() {
		return s3Path;
	}

	public boolean isParsed() {
		return parsed;
	}

	public boolean isDNATest() {
		return isDNATest;
	}

}
