package edu.utah.seq.vcf.json.tempusv3;

public class TempusV3JsonSummary {

	private TempusV3Order tempusV3Order = null;
	private TempusV3Patient tempusV3Patient = null;
	private TempusV3Report tempusV3Report = null;
	private TempusV3Specimen[] tempusV3Specimens = null;
	
	TempusV3JsonSummary (TempusV3Order tempusV3Order, TempusV3Patient tempusV3Patient, TempusV3Report tempusV3Report, TempusV3Specimen[] tempusSpecimens){
		this.tempusV3Order = tempusV3Order;
		this.tempusV3Patient = tempusV3Patient;
		this.tempusV3Report = tempusV3Report;
		this.tempusV3Specimens = tempusSpecimens;
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
}
