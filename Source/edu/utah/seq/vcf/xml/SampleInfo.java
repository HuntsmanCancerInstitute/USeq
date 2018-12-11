package edu.utah.seq.vcf.xml;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.TreeMap;

import util.gen.Misc;

public class SampleInfo {
	
	//fields
	private String sampleName;
	private String submittedDiagnosis;
	private String specSite;
	private String disease;
	private String diseaseOntology;
	private String pathologyDiagnosis;
	private String tissueOfOrigin;
	private String qualityControl;
	
	public SampleInfo(LinkedHashMap<String,String> ra){
		sampleName= ra.get("test-request");
		submittedDiagnosis= ra.get("SubmittedDiagnosis");
		specSite= ra.get("SpecSite");
		disease= ra.get("disease");
		diseaseOntology= ra.get("disease-ontology");
		pathologyDiagnosis= ra.get("pathology-diagnosis");
		tissueOfOrigin= ra.get("tissue-of-origin");
		qualityControl= ra.get("quality-control");
		
		if (submittedDiagnosis != null) submittedDiagnosis = submittedDiagnosis.toLowerCase();
		if (specSite != null) specSite = specSite.toLowerCase();
		if (disease != null) disease = disease.toLowerCase();
		else disease = "unknown";
		if (diseaseOntology != null) diseaseOntology = diseaseOntology.toLowerCase();
		if (pathologyDiagnosis != null) pathologyDiagnosis = pathologyDiagnosis.toLowerCase();
		if (tissueOfOrigin != null) tissueOfOrigin = tissueOfOrigin.toLowerCase();
		if (qualityControl != null) qualityControl = qualityControl.toLowerCase();
	}
	
	public boolean allLoaded(){
		if (sampleName == null || submittedDiagnosis == null || specSite == null || 
				disease == null || diseaseOntology == null | pathologyDiagnosis == null 
				|| tissueOfOrigin == null || qualityControl == null) return false;
		return true;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("Name\t"+sampleName); sb.append("\n");
		sb.append("QualityControl\t"+qualityControl); sb.append("\n");
		sb.append("SubmittedDiagnosis\t"+submittedDiagnosis); sb.append("\n");
		sb.append("SpecSite\t"+specSite); sb.append("\n");
		sb.append("Disease\t"+disease); sb.append("\n");
		sb.append("DiseaseOntology\t"+diseaseOntology); sb.append("\n");
		sb.append("PathologyDiagnosis\t"+pathologyDiagnosis); sb.append("\n");
		sb.append("TissueOfOrigin\t"+tissueOfOrigin); sb.append("\n");
		return sb.toString();
	}
	
	public static String fetchSummaryInfo(ArrayList<SampleInfo> sampleInfo) {
		//aggregate items
		TreeMap<String, ArrayList<String>> submittedDiagnosisTM = new TreeMap<String, ArrayList<String>>();
		TreeMap<String, ArrayList<String>> specSiteTM = new TreeMap<String, ArrayList<String>>();
		TreeMap<String, ArrayList<String>> diseaseTM = new TreeMap<String, ArrayList<String>>();
		TreeMap<String, ArrayList<String>> diseaseOntologyTM = new TreeMap<String, ArrayList<String>>();
		TreeMap<String, ArrayList<String>> pathologyDiagnosisTM = new TreeMap<String, ArrayList<String>>();
		TreeMap<String, ArrayList<String>> tissueOfOriginTM = new TreeMap<String, ArrayList<String>>();
		TreeMap<String, ArrayList<String>> qualityControlTM = new TreeMap<String, ArrayList<String>>();
		
		for (SampleInfo si: sampleInfo){
			if (si.getSubmittedDiagnosis() != null) {
				ArrayList<String> samples = submittedDiagnosisTM.get(si.getSubmittedDiagnosis());
				if (samples == null) {
					samples = new ArrayList<String>();
					submittedDiagnosisTM.put(si.getSubmittedDiagnosis(), samples);
				}
				samples.add(si.getSampleName());
			}
			
			if (si.getSpecSite() != null) {
				ArrayList<String> samples = specSiteTM.get(si.getSpecSite());
				if (samples == null) {
					samples = new ArrayList<String>();
					specSiteTM.put(si.getSpecSite(), samples);
				}
				samples.add(si.getSampleName());
			}
			
			if (si.getDisease() != null) {
				ArrayList<String> samples = diseaseTM.get(si.getDisease());
				if (samples == null) {
					samples = new ArrayList<String>();
					diseaseTM.put(si.getDisease(), samples);
				}
				samples.add(si.getSampleName());
			}
			
			if (si.getDiseaseOntology() != null) {
				ArrayList<String> samples = diseaseOntologyTM.get(si.getDiseaseOntology());
				if (samples == null) {
					samples = new ArrayList<String>();
					diseaseOntologyTM.put(si.getDiseaseOntology(), samples);
				}
				samples.add(si.getSampleName());
			}
			
			if (si.getPathologyDiagnosis() != null) {
				ArrayList<String> samples = pathologyDiagnosisTM.get(si.getPathologyDiagnosis());
				if (samples == null) {
					samples = new ArrayList<String>();
					pathologyDiagnosisTM.put(si.getPathologyDiagnosis(), samples);
				}
				samples.add(si.getSampleName());
			}
			
			if (si.getTissueOfOrigin() != null) {
				ArrayList<String> samples = tissueOfOriginTM.get(si.getTissueOfOrigin());
				if (samples == null) {
					samples = new ArrayList<String>();
					tissueOfOriginTM.put(si.getTissueOfOrigin(), samples);
				}
				samples.add(si.getSampleName());
			}
			
			if (si.getQualityControl() != null) {
				ArrayList<String> samples = qualityControlTM.get(si.getQualityControl());
				if (samples == null) {
					samples = new ArrayList<String>();
					qualityControlTM.put(si.getQualityControl(), samples);
				}
				samples.add(si.getSampleName());
			}
		}
		
		
		//summarize
		StringBuilder sb = new StringBuilder();
		addInfo(sb, "Disease", diseaseTM);
		addInfo(sb, "DiseaseOntology", diseaseOntologyTM);
		addInfo(sb, "SubmittedDiagnosis", submittedDiagnosisTM);
		addInfo(sb, "PathologyDiagnosis", pathologyDiagnosisTM);
		addInfo(sb, "TissueOfOrigin", tissueOfOriginTM);
		addInfo(sb, "SpecSite", specSiteTM);
		addInfo(sb, "QualityControl", qualityControlTM);
		
		return sb.toString();
	}

	
	private static void addInfo(StringBuilder sb, String name, TreeMap<String, ArrayList<String>> type) {
		sb.append(name); sb.append(":\n");
		for (String typeName: type.keySet()){
			ArrayList<String> samples = type.get(typeName);
			sb.append("\t");
			sb.append(typeName); 
			sb.append("\t");
			sb.append(samples.size());
			sb.append("\t");
			for (String sampleName: samples){
				sb.append(sampleName);
				sb.append(" ");
			}
			sb.append("\n");
		}
		sb.append("\n");
	}

	public String getSampleName() {
		return sampleName;
	}
	public String getSubmittedDiagnosis() {
		return submittedDiagnosis;
	}
	public String getSpecSite() {
		return specSite;
	}
	public String getDisease() {
		return disease;
	}
	public String getDiseaseOntology() {
		return diseaseOntology;
	}
	public String getPathologyDiagnosis() {
		return pathologyDiagnosis;
	}
	public String getTissueOfOrigin() {
		return tissueOfOrigin;
	}

	public String getQualityControl() {
		return qualityControl;
	}
	
}
