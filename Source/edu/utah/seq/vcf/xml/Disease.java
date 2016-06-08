package edu.utah.seq.vcf.xml;

import java.util.LinkedHashMap;

public class Disease {
	String diseaseOntology;
	String pathologyDiagnosis;
	
	public Disease (String diseaseOntology, String pathologyDiagnosis){
		this.diseaseOntology = diseaseOntology;
		this.pathologyDiagnosis = pathologyDiagnosis;
	}

	public Disease(LinkedHashMap<String, String> reportAttributes) {
		diseaseOntology = reportAttributes.get("disease-ontology");
		if (diseaseOntology != null) diseaseOntology = diseaseOntology.toLowerCase();
		pathologyDiagnosis = reportAttributes.get("pathology-diagnosis");
	}

	public String getDiseaseOntology() {
		return diseaseOntology;
	}

	public String getPathologyDiagnosis() {
		return pathologyDiagnosis;
	}
}
