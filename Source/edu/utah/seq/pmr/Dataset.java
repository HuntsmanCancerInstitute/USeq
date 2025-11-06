package edu.utah.seq.pmr;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import edu.utah.seq.vcf.json.TempusJsonSummary;
import edu.utah.seq.vcf.json.TempusOrder;
import edu.utah.seq.vcf.json.TempusPatient;
import edu.utah.seq.vcf.json.TempusReport;
import edu.utah.seq.vcf.json.TempusSpecimen;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonSummary;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Specimen;
import util.gen.IO;

public class Dataset {
	
	private String source;		//Avatar,                             Tempus,       Caris
	private String datasetId;	//A032049_SL419345_SL419548_SL420681, TL-20-B70ACE, TN21-123058_2022-01-27
	private ArrayList<String> partialPaths = new ArrayList<String>();
	private ArrayList<File> clinicalInfoFiles = null;
	private AvatarClinicalInfo avatarClinicalInfo = null;  //only for Avatar datasets
	private TempusJsonSummary tempusJsonReportInfoPreV3 = null; //only for Tempus datasets before v3
	private TempusV3JsonSummary tempusJsonReportInfoV3 = null; //for Tempus datasets v3+
	private LinkedHashMap<String, String> carisClinicalInfo = null; //only for Caris datasets

	public Dataset(String source, String datasetId) {
		this.source = source;
		this.datasetId = datasetId;
	}
	
	
	
	public String toString(String patientMolecularRepoId) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append(patientMolecularRepoId+"\n");
		sb.append("  File(s) : "+IO.fileArrayListToString(clinicalInfoFiles, ", ")+"\n");  
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
		//Tempus pre v3
		else if (tempusJsonReportInfoPreV3 != null) {
			LinkedHashMap<String, String> meta = new LinkedHashMap<String, String>();
			tempusJsonReportInfoPreV3.getTempusPatient().addAttributes(meta, false);
			tempusJsonReportInfoPreV3.getTempusOrder().addAttributes(meta);
			tempusJsonReportInfoPreV3.getTempusReport().addAttributes(meta);
			TempusSpecimen.addAttributes(meta, tempusJsonReportInfoPreV3.getTempusSpecimens());
			for (String key: meta.keySet()) {
				sb.append("  ");
				sb.append(key); 
				sb.append(" : "); 
				sb.append(meta.get(key)); 
				sb.append("\n");
			}
		}
		//Tempus v3
		else if (tempusJsonReportInfoV3 != null) {
			LinkedHashMap<String, String> meta = new LinkedHashMap<String, String>();
			tempusJsonReportInfoV3.getTempusPatient().addAttributes(meta, false);
			tempusJsonReportInfoV3.getTempusOrder().addAttributes(meta);
			tempusJsonReportInfoV3.getTempusReport().addAttributes(meta);
			TempusV3Specimen.addAttributes(meta, tempusJsonReportInfoV3.getTempusSpecimens());
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

	public TempusJsonSummary getTempusJsonReportInfoPreV3() {
		return tempusJsonReportInfoPreV3;
	}

	public void setTempusJsonReportInfoPreV3(TempusJsonSummary tempusJsonReportInfo) {
		this.tempusJsonReportInfoPreV3 = tempusJsonReportInfo;
	}
	
	public TempusV3JsonSummary getTempusJsonReportInfoV3() {
		return tempusJsonReportInfoV3;
	}

	public void setTempusJsonReportInfoV3(TempusV3JsonSummary tempusJsonReportInfo) {
		this.tempusJsonReportInfoV3 = tempusJsonReportInfo;
	}

	public void setCarisXmlReportInfo(LinkedHashMap<String, String> linkedHashMap) {
		carisClinicalInfo = linkedHashMap;
		
	}

	public LinkedHashMap<String, String> getCarisClinicalInfo() {
		return carisClinicalInfo;
	}

}
