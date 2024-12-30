package edu.utah.seq.run.avproj.adw;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;

import org.json.JSONArray;
import org.json.JSONObject;

import edu.utah.hci.misc.Util;
import util.gen.IO;
import util.gen.Misc;

/**Represents one patients tumor and normal datasets for TNRunner2 processing.  Some patients will have several of these due to multiple tumor samples and or multiple platforms. 
 * The tumor exome and tumor RNA are merged into the same Tumor Sample*/
public class AvatarAnalysisJob {

	//fields
	private TumorSampleADW tumorSample = null;
	private ArrayList<NormalSampleADW> normalSamples = new ArrayList<NormalSampleADW>();
	private boolean matchedPlatform = true;
	private String submittingGroup = null;
	private PatientADW patient = null;
	private File normalDNAFastqDir = null;

	public AvatarAnalysisJob(PatientADW patient, TumorSampleADW tumorSample, ArrayList<NormalSampleADW> normalSamples, boolean matchedPlatform) {
		this.patient = patient;
		this.tumorSample = tumorSample;
		this.normalSamples = normalSamples;
		this.matchedPlatform = matchedPlatform;
	}

	/**Returns patientId_normalId1-normalId2_tumorExomeId_tumorRnaId most will have just one normal*/
	public String getComparisonKey(String patientId) {
		StringBuilder sb = new StringBuilder();
		//patientId
		sb.append(patientId); 

		//normalIds
		if (normalSamples.size()==0) sb.append("_NA");
		else {
			sb.append("_");
			sb.append(normalSamples.get(0).getNormalDnaName());
			for (int i=1; i< normalSamples.size(); i++) {
				sb.append("-");
				sb.append(normalSamples.get(i).getNormalDnaName());
			}
		}

		//tumorExomeId
		sb.append("_");
		if (tumorSample.getTumorDnaName() == null) sb.append("NA");
		else sb.append(tumorSample.getTumorDnaName());

		//tumorRnaId
		sb.append("_");
		if (tumorSample.getTumorRnaName() == null) sb.append("NA");
		else sb.append(tumorSample.getTumorRnaName());

		return sb.toString();
	}

	public TumorSampleADW getTumorSample() {
		return tumorSample;
	}

	public ArrayList<NormalSampleADW> getNormalSamples() {
		return normalSamples;
	}

	public boolean isMatchedPlatform() {
		return matchedPlatform;
	}

	public void makeAnalysisJobAsterOld(File fastqDownloadRoot, String nameAJ, File testDir, ArrayList<String> cmds, ClinicalMolLinkage linkage, HashSet<String> dirPathsToDownload) throws IOException {
		String downloadRoot = fastqDownloadRoot.getCanonicalPath();

		//Fastq
		File fastq = new File (testDir, "Fastq");
		fastq.mkdir();
		//TumorDNA
		if (tumorSample.getTumorDnaName() != null) {
			tumorSample.getTumorDnaFastqFiles();
			File tumorDNA = new File (fastq, "TumorDNA");
			tumorDNA.mkdir();
			String tumorFastqDir = tumorDNA.getCanonicalPath()+"/";

			//make link commands, downloads with Aster are a major issue
			for (String path: tumorSample.getTumorWesFastqPathsToFetch()) {
				String ln = "ln -s "+ downloadRoot+ path+ " "+ tumorFastqDir;
				cmds.add(ln);
			}
			//update the dirPathsToDownload
			dirPathsToDownload.add(tumorSample.getTumorWesFastqPathsToFetch().toString());
			dirPathsToDownload.add(tumorSample.getTumorWesFastqPathsToFetch().toString());
		}
		//TumorRNA
		if (tumorSample.getTumorRnaName()!= null) {
			File tumorRNA = new File (fastq, "TumorRNA");
			tumorRNA.mkdir();
			String tumorFastqDir = tumorRNA.getCanonicalPath()+"/";

			//make link commands, downloads with Aster are a major issue
			for (String path: tumorSample.getTumorRnaFastqPathsToFetch()) {
				String ln = "ln -s "+ downloadRoot+ path+ " "+ tumorFastqDir;
				cmds.add(ln);
			}
			
			//update the dirPathsToDownload
			dirPathsToDownload.add(tumorSample.getTumorRnaFastqPathsToFetch().toString());
			
		}
		//NormalDNAs
		if (normalSamples.size()!=0) {
			File normalDNA = new File (fastq, "NormalDNA");
			normalDNA.mkdir();
			String normalFastqDir = normalDNA.getCanonicalPath()+"/";
			for (NormalSampleADW ns: normalSamples) {

				//make link commands, downloads with Aster are a major issue
				for (String path: ns.getNormalWesFastqPathsToFetch()) {
					String ln = "ln -s "+ downloadRoot+ path+ " "+ normalFastqDir;
					cmds.add(ln);
				}
				
				//update the dirPathsToDownload
				dirPathsToDownload.add(ns.getNormalWesFastqPathsToFetch().toString());
			} 
		}
		//ClinicalReport
		File clinicalReport = new File (testDir, "ClinicalReport");
		clinicalReport.mkdir();
		writeJson(nameAJ, clinicalReport, linkage);
	}

	public void makeAnalysisJobAsterNew(File fastqDownloadRoot, String nameAJ, File testDir, ArrayList<String> cmds, ClinicalMolLinkage linkage, LinkedHashSet<String> filesToDownload) throws IOException {
		String downloadRoot = fastqDownloadRoot.getCanonicalPath();
		//Fastq
		File fastq = new File (testDir, "Fastq");
		fastq.mkdir();
		//TumorDNA
		if (tumorSample.getTumorDnaName() != null) {
			tumorSample.getTumorDnaFastqFiles();
			File tumorDNA = new File (fastq, "TumorDNA");
			tumorDNA.mkdir();
			String tumorFastqDir = tumorDNA.getCanonicalPath()+"/";
			//for each file
			for (String path: tumorSample.getTumorWesFastqPathsToFetch()) {
				String fileName = path.substring(path.lastIndexOf('/'));
				String ln = "ln -s "+ downloadRoot + fileName+ " " +tumorFastqDir;
				cmds.add(ln);
				filesToDownload.add(path);
			}
		}
		//TumorRNA
		if (tumorSample.getTumorRnaName()!= null) {
			File tumorRNA = new File (fastq, "TumorRNA");
			tumorRNA.mkdir();
			String tumorFastqDir = tumorRNA.getCanonicalPath()+"/";

			for (String path: tumorSample.getTumorRnaFastqPathsToFetch()) {
				String fileName = path.substring(path.lastIndexOf('/'));
				String ln = "ln -s "+ downloadRoot + fileName+ " " +tumorFastqDir;
				cmds.add(ln);
				filesToDownload.add(path);
			}
			
		}
		//NormalDNAs
		if (normalSamples.size()!=0) {
			normalDNAFastqDir = new File (fastq, "NormalDNA");
			normalDNAFastqDir.mkdir();
			String normalFastqDir = normalDNAFastqDir.getCanonicalPath()+"/";
			for (NormalSampleADW ns: normalSamples) {

				for (String path: ns.getNormalWesFastqPathsToFetch()) {
					String fileName = path.substring(path.lastIndexOf('/'));
					String ln = "ln -s "+ downloadRoot + fileName+ " " +normalFastqDir;
					cmds.add(ln);
					filesToDownload.add(path);
				}
			} 
		}
		//ClinicalReport
		File clinicalReport = new File (testDir, "ClinicalReport");
		clinicalReport.mkdir();
		writeJson(nameAJ, clinicalReport, linkage);
	}


	private void writeJson(String nameAJ, File clinicalReport, ClinicalMolLinkage linkage) throws IOException {

		//platform
		String platform = "NA";
		if (isMatchedPlatform() == false) platform = "Mixed";
		else if (tumorSample.getPlatformName() != null) platform = tumorSample.getPlatformName();
		else if (normalSamples.size()!=0) {
			NormalSampleADW ns = normalSamples.get(0);
			platform = ns.getPlatformName();
		}

		//group from tumor DNA sample then try from tumor RNA sample
		String group = "NA";
		if (tumorSample.getTumorDnaName() != null) {
			group = linkage.getKeyDiseaseType().get(patient.getPatientId()+"_"+tumorSample.getTumorDnaName());
			if (group == null) throw new IOException("\nFailed to find a disease type for the tumor sample "+tumorSample.getTumorDnaName());
		}
		else if (tumorSample.getTumorRnaName() != null) {
			group = linkage.getKeyDiseaseType().get(patient.getPatientId()+"_"+tumorSample.getTumorRnaName());
			if (group == null) throw new IOException("\nFailed to find a disease type for the tumor sample "+tumorSample.getTumorDnaName());
		}

		//testId_platform_groupId_gender.json
		File json = new File (clinicalReport, nameAJ+"_"+platform+"_"+group+"_"+patient.getSubjectMatchMaker().getGender()+".json"); 

		if (platform == null) {
			Misc.printErrAndExit("Null plat "+this.getComparisonKey(patient.getPatientId()));
		}

		//non matched?
		boolean mixed = false;
		if (platform.equals("Mixed")) {
			File f = new File (clinicalReport, "NON_MATCHED_CAPTURE_PLATFORMS");
			f.createNewFile();
			mixed = true;
		}

		JSONObject main = new JSONObject();
		main.put("mixedCapturePlatforms", mixed);
		main.put("diseaseGroup", group);

		//Patient json
		JSONObject p = patient.fetchJson();
		main.put("Patient", p);

		//TumorDNA json
		if (tumorSample.getTumorDnaName() != null) {
			JSONObject tj = tumorSample.fetchTumorDnaJson(linkage);
			main.put("TumorDNA", tj);
		}

		//TumorRNA json
		if (tumorSample.getTumorRnaName()!= null) {
			JSONObject tj = tumorSample.fetchTumorRnaJson(linkage);
			main.put("TumorRNA", tj);
		}

		//NormalDNAs json
		if (normalSamples.size()!=0) {
			JSONArray ja = new JSONArray();
			for (NormalSampleADW ns: normalSamples) {
				JSONObject nj = ns.fetchJson(linkage);
				ja.put(nj);
			}
			main.put("NormalDNA", ja);

		}

		IO.writeString(main.toString(3), json);
	}

	public File getNormalDNAFastqDir() {
		return normalDNAFastqDir;
	}

}
