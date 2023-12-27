package edu.utah.seq.run.invitae;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import org.json.JSONObject;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class InvitaeReport {
	
	//fields needed for processing
	private String reportId = null;
	private String requisitionId = null;
	private String pmrId = null;
	
	//only one patient line
	private String[] patientInfo = null;
	private String subjectMatchMakerLine = null;
	private File reportDir = null;
	
	//zero or more variants
	private ArrayList<String[]> variantInfo = new ArrayList<String[]>();
	InvitaeVariant[] invitaeVariants = null;
	Bed[] testedGeneCoorHg38 = null;
	
	//constructor
	public InvitaeReport(String reportId, String requisitionId, String genes, HashMap<String, Bed> iGenCoor, String[] patientInfo) throws IOException {
		this.reportId = reportId;
		this.requisitionId = requisitionId;
		this.patientInfo = patientInfo;
		
		//check genes: ATM,BRCA1,BRCA2,CDH1,CHEK2,PALB2,PTEN,STK11,TP53; many duplicates so need to hash
		boolean ok = true;
		String[] genesTested = Misc.COMMA.split(genes);
		LinkedHashSet<String> uniqueGenes = new LinkedHashSet<String>();
		for (String g: genesTested) uniqueGenes.add(g);
		testedGeneCoorHg38 = new Bed[uniqueGenes.size()];
		int counter = 0;
		for (String g: uniqueGenes) {
			String geneName = Misc.WHITESPACE.matcher(g).replaceAll("_").toUpperCase();
			Bed coor = iGenCoor.get(geneName);
			if (coor == null) {
				ok = false;
				IO.el("ERROR: missing hg38 genomic coordinates for '"+g+"' '"+geneName+"'");
			}
			else testedGeneCoorHg38[counter++] = coor;
		}
		if (ok == false) throw new IOException("Add missing Hg38 gene coordinates to the interrogated gene file.");
	}
	
	public String getPatientInfo(int index) {
		String f = patientInfo[index];
		if (f!=null) return f.trim();
		return null;
	}

	public String getReportId() {
		return reportId;
	}
	
	public void parseVariants(InvitaeBulkCsvParser idw) throws IOException {
		if (variantInfo.size() == 0) return;
		
		//for each info line, some need to be skipped
		ArrayList<InvitaeVariant> passingVars = new ArrayList<InvitaeVariant>();
		for (int i=0; i< variantInfo.size(); i++) {
			InvitaeVariant iv = new InvitaeVariant(variantInfo.get(i), idw);
			if (iv.isPassed()) passingVars.add(iv);
		}
		//any passing variants?
		if (passingVars.size()!=0) {
			invitaeVariants = new InvitaeVariant[passingVars.size()];
			passingVars.toArray(invitaeVariants);
			Arrays.sort(invitaeVariants);			
		}
	}
	
	/**lastName firstName dobMonth(1-12) dobDay(1-31) dobYear(1900-2050) gender(M|F) mrn 
	 * @throws IOException */
	public void setSubjectMatchMakerLine(int lastName, int firstName, int dob, int gender, int mrn) throws IOException {
		StringBuilder sb = new StringBuilder();
		
		//LastName
		String ln = patientInfo[lastName].trim();
		if (ln.length()!=0) sb.append(ln);
		sb.append("\t");

		//FirstName, sometimes    Kathleen "Kathy", uggg!
		String fn = patientInfo[firstName].trim();
		if (fn.contains("\"")) {
			String[] split = Misc.DOUBLE_QUOTE.split(fn);
			fn = split[0].trim();
		}
		if (fn.length()!=0) sb.append(fn);
		sb.append("\t");
		
		// 4/25/1961
		
		int dobMonth = -1;
		int dobDay = -1;
		int dobYear = -1;
		String dobAll = patientInfo[dob].trim();
		if (dobAll.length()!=0) {
			String[] splitDoB = Misc.FORWARD_SLASH.split(dobAll);
			if (splitDoB.length!=3) throw new IOException("ERROR parsing dob from "+reportId+" see "+dobAll);
			dobMonth = Integer.parseInt(splitDoB[0]);
			dobDay = Integer.parseInt(splitDoB[1]);
			dobYear = Integer.parseInt(splitDoB[2]);
		}
		
		//Month
		if (dobMonth != -1) sb.append(new Integer(dobMonth).toString());
		sb.append("\t");
		
		//Day
		if (dobDay != -1) sb.append(new Integer(dobDay).toString());
		sb.append("\t");
		
		//Year
		if (dobYear != -1) sb.append(new Integer(dobYear).toString());
		sb.append("\t");
		
		//Gender
		String g = patientInfo[gender].toLowerCase().trim();
		if (g.length()!=0) {
			if (g.equals("female")) sb.append("F");
			else if (g.equals("male")) sb.append("M");
			else throw new IOException("ERROR parsing gender from "+reportId+" see "+g);
		}
		sb.append("\t");
		
		String m = patientInfo[mrn];
		if (m.length()!=0) sb.append(m);
		
		subjectMatchMakerLine = sb.toString();
	}

	public String getSubjectMatchMakerLine() {
		return subjectMatchMakerLine;
	}

	public String getPmrId() {
		return pmrId;
	}

	public void setPmrId(String pmrId) {
		this.pmrId = pmrId;
	}

	public void addClinicalReport(File jobsDirectory, int indexOrderingPhysician, int indexReportDate, String[] headerNames) throws IOException {
		reportDir = new File(jobsDirectory, pmrId+"/Invitae/"+requisitionId+"_"+reportId);
		reportDir.mkdirs();
		File crDir = new File(reportDir, "ClinicalReport");
		crDir.mkdir();
		if (crDir.exists()==false) throw new IOException("ERROR: failed to make the ClinicalReport dir in the job dir for "+reportDir);
		
		//A042126_SL600535_SL600473_SL567649_TWSv2_GU_F.json - Avatar
		//Id   capture   diseaseGroup  gender
		//TL-23-Z6XQ3HVI_XT.V4_2023-08-24_deid_Sonam_Puri_M.json - Tempus
		//Id  test  date  ordering physician  gender
		
		//OP
		String oc = patientInfo[indexOrderingPhysician];
		oc = Misc.WHITESPACE.matcher(oc).replaceAll("_");
		//ReportDate 2/23/2023 19:02
		String ymd = parseReportDate(patientInfo[indexReportDate]);
		//ID reportDate orderingClinician - Tempus, don't need gender since we're not recalling CNVs
		File json = new File(crDir, requisitionId+"_"+reportId+"_"+ymd+"_"+oc+".json");
		
		JSONObject main = new JSONObject();
		// watch out and exclude patient_first_name	patient_middle_initial	patient_last_name	mrn	dob
		for (int i=0; i< patientInfo.length; i++) {
			if (headerNames[i].startsWith("patient_") || headerNames[i].equals("mrn") || headerNames[i].equals("dob")) continue;
			if (patientInfo[i].length()!=0) main.put(headerNames[i], patientInfo[i]);
		}
		IO.writeString(main.toString(3), json);
		
	}
	
	public void addDataFiles(String vcfHeader, String bedHeader) throws IOException {
	/*AJobs/xjk2RU7Ybx/Avatar/
	 	A043846_FT-SA238602_FT-SA238627D_FT-SA238627R/
	 		GermlineVariantCalling/
	 			A043846_FT-SA238602_FT-SA238627D_FT-SA238627R_Anno/
	 				Vcfs/
	 					A043846_FT-SA238602_FT-SA238627D_FT-SA238627R_Anno_Hg38.anno.vcf.gz
	 IJobs/AVw5Cf7VHA/Invitae/
	 	RQ322930_AR302357/
	 		GermlineVariantCalling/
	 			RQ322930_AR302357_Hg19/
	 				RQ322930_AR302357_Hg19.txt				- from their spreadsheet, don't use xxx.vcf since tabix will read the END and error out when before the pos
	 				 
	 			RQ322930_AR302357_Hg38/					
	 				RQ322930_AR302357_Hg38_TestedGenes.bed.gz 	- from the spreadsheet genes
	 				RQ322930_AR302357_Hg38.vcf.gz 				- from CrossMap
	 			RQ322930_AR302357_Hg38_Anno/					- to do with the Annotator workflow		
	*/
		//make data folders
		File gDir = new File(reportDir, "GermlineVariantCalling/");
		gDir.mkdirs();
		
		// for the interrogated gene file
		File hg38Dir = new File(gDir, reportDir.getName()+"_Hg38");
		hg38Dir.mkdir();
		addBed(new File(hg38Dir, reportDir.getName()+"_Hg38.bed"), bedHeader);
		
		// any variants?  lots of tests have no findings, thus "No Pathogenic sequence variants or deletions/duplications identified."
		if (invitaeVariants!= null) {
			File hg19Dir = new File(gDir, reportDir.getName()+"_Hg19");
			hg19Dir.mkdir();
			//add vcf
			addVcf(vcfHeader, new File (hg19Dir, reportDir.getName()+"_Hg19.vcf"));
		}
	}
	
	private void addBed(File file, String bedHeader) throws IOException {
		//sort the regions
		Arrays.sort(testedGeneCoorHg38);
		//write them out with the header
		PrintWriter out = new PrintWriter (new FileWriter( file));
		out.println(bedHeader);
		for (Bed b: testedGeneCoorHg38) out.println(b.toString());
		out.close();
	}

	private void addVcf(String vcfHeader, File file) throws IOException {
		PrintWriter out = new PrintWriter (new FileWriter( file));
		out.print(vcfHeader);
		if (invitaeVariants != null) {
			for (int i=0; i< invitaeVariants.length; i++) out.println(invitaeVariants[i].toVcf(i)); 
		}
		out.close();
	}

	public static String parseReportDate(String rd) throws IOException {
		//   2/23/2023 19:02
		String slashDate = Misc.WHITESPACE.split(rd)[0];
		String[] dmy = Misc.FORWARD_SLASH.split(slashDate);
		if (dmy.length != 3) throw new IOException("ERROR: failed to parse the date from "+rd);
		return dmy[2]+"-"+dmy[1]+"-"+dmy[0];
	}

	public ArrayList<String[]> getVariantInfo() {
		return variantInfo;
	}



}
