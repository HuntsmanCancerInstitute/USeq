package edu.utah.seq.vcf.json.tempusv3;

import java.io.File;

import org.json.JSONObject;

public class TempusV3JsonSummary {

	private JSONObject mainObject = null;
	private TempusV3Order tempusV3Order = null;
	private TempusV3Patient tempusV3Patient = null;
	private TempusV3Report tempusV3Report = null;
	private TempusV3Specimen[] tempusV3Specimens = null;
	private TempusV3GenomicVariants tempusV3Results = null;
	
	TempusV3JsonSummary (JSONObject mainObject, TempusV3Order tempusV3Order, TempusV3Patient tempusV3Patient, TempusV3Report tempusV3Report, TempusV3Specimen[] tempusSpecimens, TempusV3GenomicVariants tempusV3Results){
		this.mainObject = mainObject;
		this.tempusV3Order = tempusV3Order;
		this.tempusV3Patient = tempusV3Patient;
		this.tempusV3Report = tempusV3Report;
		this.tempusV3Specimens = tempusSpecimens;
		this.tempusV3Results = tempusV3Results;
	}

	public TempusV3Order getTempusOrder() {
		return tempusV3Order;
	}

	public TempusV3Patient getTempusPatient() {
		return tempusV3Patient;
	}

	public TempusV3Report getTempusReport() {
		return tempusV3Report;
	}

	public TempusV3Specimen[] getTempusSpecimens() {
		return tempusV3Specimens;
	}

	public TempusV3GenomicVariants getTempusV3Results() {
		return tempusV3Results;
	}

	public JSONObject getMainObject() {
		return mainObject;
	}
}
