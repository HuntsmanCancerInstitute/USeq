package edu.utah.seq.vcf.json;

public class TempusJsonSummary {

	private TempusOrder tempusOrder = null;
	private TempusPatient tempusPatient = null;
	private TempusReport tempusReport = null;
	private TempusSpecimen[] tempusSpecimens = null;
	
	TempusJsonSummary (TempusOrder tempusOrder, TempusPatient tempusPatient, TempusReport tempusReport, TempusSpecimen[] tempusSpecimens){
		this.tempusOrder = tempusOrder;
		this.tempusPatient = tempusPatient;
		this.tempusReport = tempusReport;
		this.tempusSpecimens = tempusSpecimens;
	}

	public TempusOrder getTempusOrder() {
		return tempusOrder;
	}

	public TempusPatient getTempusPatient() {
		return tempusPatient;
	}

	public TempusReport getTempusReport() {
		return tempusReport;
	}

	public TempusSpecimen[] getTempusSpecimens() {
		return tempusSpecimens;
	}
}
