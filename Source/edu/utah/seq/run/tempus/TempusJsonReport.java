package edu.utah.seq.run.tempus;

import java.util.ArrayList;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusJsonReport {

	private String reportId = null;
	private String signout_date = null;
	private String bioInfPipeline = null;
	private String notes = null;
	private String reportStatus = null;
	private ArrayList<String> warningMessages = new ArrayList<String>();
	
	/*
    "report": {
        "reportId": "03cfd69c-5ea0-4ba6-92f6-bc43175c1222",
        "workflow": {
            "reportStatus": "standard",
            "details": null,
            "reportType": "DNA"
        },
        "signing_pathologist": "Timothy Taxter, M.D.",
        "signout_date": "2018-10-12T02:50:44+00:00",
        "bioInfPipeline": "1.3.5",
        "notes": "The tumor shows a loss of heterozygosity in TP53.RNA analysis is being performed and will be reported in the Tempus online portal when complete."
     },

	 */
	public TempusJsonReport(JSONObject main) throws JSONException {
		JSONObject report = main.getJSONObject("report");
		reportId = Json.getStringAttribute(report, "reportId");
		signout_date = Json.getStringAttribute(report, "signout_date");
		if (signout_date == null) signout_date = Json.getStringAttribute(report, "signoutDate");
		bioInfPipeline = Json.getStringAttribute(report, "bioInfPipeline");
		notes = Json.getStringAttribute(report, "notes");
		if (notes != null) notes = Misc.WHITESPACE.matcher(notes).replaceAll(" ");
		if (report.has("workflow")) {
			JSONObject workflow = report.getJSONObject("workflow");
			//report status, watch this one, if qns then expect junk results
			reportStatus = Json.getStringAttribute(workflow, "reportStatus");
			if (reportStatus.contains("qns")) warningMessages.add("Report Status indicates quality control problem. View results with caution.");
			
		}
		
	}

	public String getReportId() {
		return reportId;
	}

	public String getSignout_date() {
		return signout_date;
	}

	public String getBioInfPipeline() {
		return bioInfPipeline;
	}

	public String getNotes() {
		return notes;
	}

	public ArrayList<String> getWarningMessages() {
		return warningMessages;
	}
}
