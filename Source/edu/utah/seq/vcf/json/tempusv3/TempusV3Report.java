package edu.utah.seq.vcf.json.tempusv3;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusV3Report {

	private String reportId = null;
	private String signoutDate = null;
	private String bioInfoPipeline = null;
	private String notes = null;
	private String reportStatus = null;
	private ArrayList<String> warningMessages = new ArrayList<String>();
	
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
	public TempusV3Report(JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException {
		JSONObject report = object.getJSONObject("report");
		reportId = Json.getStringAttribute(report, "reportId");
		signoutDate = Json.getStringAttribute(report, "signoutDate");
		if (signoutDate == null) signoutDate = Json.getStringAttribute(report, "signoutDate");
		bioInfoPipeline = Json.getStringAttribute(report, "bioInfoPipeline");
		notes = Json.getStringAttribute(report, "notes");
		if (notes != null) notes = Misc.WHITESPACE.matcher(notes).replaceAll(" ");
		TempusV3Json2Vcf.add(bioInfoPipeline, tempusJson2Vcf.bioInfoPipeline);

		if (report.has("workflow")) {
			JSONObject workflow = report.getJSONObject("workflow");
			//report status, watch this one, if qns then expect junk results
			reportStatus = Json.getStringAttribute(workflow, "reportStatus");
			if (reportStatus.contains("standard")==false) warningMessages.add("Report Status indicates it is not standard, possible quality control problem. View results with caution.");
			TempusV3Json2Vcf.add(reportStatus, tempusJson2Vcf.reportStatus);
		}
		
	}

	public String getReportId() {
		return reportId;
	}

	public String getsignoutDate() {
		return signoutDate;
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
		meta.put("tempusReportId", reportId);
		meta.put("tempusReportStatus", reportStatus);
		meta.put("tempusSignoutDate", signoutDate);
		meta.put("tempusBioInfPipeline", bioInfoPipeline);
		if (notes != null) meta.put("tempusReportNotes", notes);
	}
	/**Addes to spreadsheet*/
	public void addAttributes(LinkedHashMap<String, String> meta) {
		meta.put("reportId", reportId);
		meta.put("reportStatus", reportStatus);
		meta.put("signoutDate", signoutDate);
		meta.put("bioInfPipeline", bioInfoPipeline);
		if (notes != null) meta.put("notes", notes);
	}
	
}
