package edu.utah.seq.pmr;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import edu.utah.seq.vcf.json.TempusJsonSummary;

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
