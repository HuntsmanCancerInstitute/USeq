package edu.utah.seq.vcf.json.tempusv3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;

import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusV3Report {

	private File jsonFile = null;
	private String reportId = null;
	private String signoutDate = null;
	private String bioInfoPipeline = null;
	private String notes = null;
	private String reportStatus = null;
	private ArrayList<String> warningMessages = new ArrayList<String>();
	private static final Pattern patT = Pattern.compile("T");
	
	/*
     "report": {
          "notes": null,
          "reportId": "0e18a291-9f50-4daa-b027-e171861f8edb",
          "workflow": {
               "reportType": "DNA",
               "reportStatus": "standard",
               "details": "No potentially actionable variants and no reportable treatment options found."
          },
          "signingPathologist": "Rongqin Ren, MD, PhD",
          "modifiesReportId": "N/A",
          "bioInfoPipeline": "5.2.5",
          "signoutDate": "2025-03-02T19:46:09"
     },
	 */
	public TempusV3Report(File jsonFile, JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException {
		this.jsonFile = jsonFile;
		JSONObject report = object.getJSONObject("report");
		reportId = Json.getStringAttribute(report, "reportId");
		signoutDate = Json.getStringAttribute(report, "signoutDate");
		bioInfoPipeline = Json.getStringAttribute(report, "bioInfoPipeline");
		notes = Json.getStringAttribute(report, "notes");
		if (notes != null) notes = Misc.WHITESPACE.matcher(notes).replaceAll(" ");
		if (tempusJson2Vcf!=null) TempusV3Json2Vcf.add(bioInfoPipeline, tempusJson2Vcf.bioInfoPipeline);

		if (report.has("workflow")) {
			JSONObject workflow = report.getJSONObject("workflow");
			//report status, watch this one, if qns then expect junk results
			reportStatus = Json.getStringAttribute(workflow, "reportStatus");
			if (reportStatus.contains("standard")==false) warningMessages.add("Report Status indicates it is not standard, possible quality control problem. View results with caution.");
			if (tempusJson2Vcf!=null) TempusV3Json2Vcf.add(reportStatus, tempusJson2Vcf.reportStatus);
		}
		
	}
	
	public String getSignOutDateNoTime() throws IOException {
		//"signoutDate": "2025-03-02T19:46:09"
		String[] split = patT.split(signoutDate);
		if (split.length !=2 || split[0].startsWith("20")==false) throw new IOException("Failed to split the sign out data "+signoutDate);
		return split[0];
	}

	public String getReportId() {
		return reportId;
	}

	public String getBioInfPipeline() {
		return bioInfoPipeline;
	}

	public String getNotes() {
		return notes;
	}

	public ArrayList<String> getWarningMessages() {
		return warningMessages;
	}

	public void setWarningMessages(ArrayList<String> warningMessages) {
		this.warningMessages = warningMessages;
	}

	/**Added to Vcf header*/
	public void addMetaData(LinkedHashMap<String, String> meta) {
		Misc.addConcatinatedValue (meta,"tempusReportId", reportId, "; ");
		Misc.addConcatinatedValue (meta,"tempusReportStatus", reportStatus, "; ");
		Misc.addConcatinatedValue (meta,"tempusSignoutDate", signoutDate, "; ");
		Misc.addConcatinatedValue (meta,"tempusBioInfPipeline", bioInfoPipeline, "; ");
		Misc.addConcatinatedValue (meta, "tempusReportNotes", notes, "; ");
	}
	/**Added to spreadsheet*/
	public void addAttributes(LinkedHashMap<String, String> meta) {
		meta.put("jsonFile", jsonFile.toString());
		meta.put("reportId", reportId);
		meta.put("reportStatus", reportStatus);
		meta.put("signoutDate", signoutDate);
		meta.put("bioInfPipeline", bioInfoPipeline);
		if (notes != null) meta.put("notes", notes);
	}

	public File getJsonFile() {
		return jsonFile;
	}

	public String getSignoutDate() {
		return signoutDate;
	}

	public String getBioInfoPipeline() {
		return bioInfoPipeline;
	}

	public String getReportStatus() {
		return reportStatus;
	}
	
}
