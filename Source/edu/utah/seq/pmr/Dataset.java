package edu.utah.seq.pmr;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import edu.utah.seq.vcf.json.TempusJsonSummary;
import edu.utah.seq.vcf.json.TempusOrder;
import edu.utah.seq.vcf.json.TempusPatient;
import edu.utah.seq.vcf.json.TempusReport;
import edu.utah.seq.vcf.json.TempusSpecimen;

public class Dataset {
	
	private String source;		//Avatar,                             Tempus,       Caris
	private String datasetId;	//A032049_SL419345_SL419548_SL420681, TL-20-B70ACE, TN21-123058_2022-01-27
	private ArrayList<String> partialPaths = new ArrayList<String>();
	private ArrayList<File> clinicalInfoFiles = null;
	private AvatarClinicalInfo avatarClinicalInfo = null;  //only for Avatar datasets
	private TempusJsonSummary tempusJsonReportInfo = null; //only for Tempus datasets
	private LinkedHashMap<String, String> carisClinicalInfo = null; //only for Caris datasets

	public Dataset(String source, String datasetId) {
		this.source = source;
		this.datasetId = datasetId;
	}
	
	public String toString(String patientMolecularRepoId) {
		StringBuilder sb = new StringBuilder();
		sb.append(patientMolecularRepoId+"\n");
		//Avatar?
		if (avatarClinicalInfo != null) sb.append(avatarClinicalInfo.toString());
		//Caris
		else if (carisClinicalInfo != null) {
			for (String key: carisClinicalInfo.keySet()) {
				sb.append("  ");
				sb.append(key); 
				sb.append(" : "); 
				sb.append(carisClinicalInfo.get(key)); 
				sb.append("\n");
			}
		}
		//Tempus
		else {
			LinkedHashMap<String, String> meta = new LinkedHashMap<String, String>();
			tempusJsonReportInfo.getTempusPatient().addAttributes(meta, false);
			tempusJsonReportInfo.getTempusOrder().addAttributes(meta);
			tempusJsonReportInfo.getTempusReport().addAttributes(meta);
			TempusSpecimen.addAttributes(meta, tempusJsonReportInfo.getTempusSpecimens());
			tempusJsonReportInfo.getTempusReport().addAttributes(meta);
			for (String key: meta.keySet()) {
				sb.append("  ");
				sb.append(key); 
				sb.append(" : "); 
				sb.append(meta.get(key)); 
				sb.append("\n");
			}
		}
		return sb.toString();
	}
	
	public String checkForMissing(String val) {
		if (val == null || val.length()==0 || val.toLowerCase().contains("null")) return null;
		return val;
	}

	public String getSource() {
		return source;
	}

	public String getDatasetId() {
		return datasetId;
	}

	public ArrayList<String> getPartialPaths() {
		return partialPaths;
	}

	public ArrayList<File> getClinicalInfoFiles() {
		return clinicalInfoFiles;
	}

	public void setClinicalInfoFile(File clinicalInfo) {
		if (clinicalInfoFiles == null) clinicalInfoFiles = new ArrayList<File>();
		clinicalInfoFiles.add(clinicalInfo);
	}

	public AvatarClinicalInfo getAvatarClinicalInfo() {
		return avatarClinicalInfo;
	}

	public void setAvatarClinicalInfo(AvatarClinicalInfo avatarClinicalInfo) {
		this.avatarClinicalInfo = avatarClinicalInfo;
	}

	public TempusJsonSummary getTempusJsonReportInfo() {
		return tempusJsonReportInfo;
	}

	public void setTempusJsonReportInfo(TempusJsonSummary tempusJsonReportInfo) {
		this.tempusJsonReportInfo = tempusJsonReportInfo;
	}

	public void setCarisXmlReportInfo(LinkedHashMap<String, String> linkedHashMap) {
		carisClinicalInfo = linkedHashMap;
		
	}

	public LinkedHashMap<String, String> getCarisClinicalInfo() {
		return carisClinicalInfo;
	}

}
