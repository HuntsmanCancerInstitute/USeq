package edu.utah.seq.vcf.xml.caris;

import java.io.File;
import java.io.IOException;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class GenomicAlteration {

	//fields
	private String biomarkerName = null;
	private String pathogenicity = null;
	private String geneName = null;
	private String hgvsCodingChange = null;
	private String chr = null;
	private String transcriptId = null;
	private String molecularConsequence = null;
	private String alleleFrequency = null;
	private String readDepth = null;
	private String ref = null;
	private String alt = null;
	
	public GenomicAlteration(NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("biomarkerName")) biomarkerName = CarisXmlVcfParser.fetch("biomarkerName", n);
				else if (name.equals("result")) pathogenicity = CarisXmlVcfParser.fetch("result", n);
				else if (name.equals("gene")) geneName = CarisXmlVcfParser.fetch("gene", n);
				else if (name.equals("hgvsCodingChange")) hgvsCodingChange = CarisXmlVcfParser.fetch("hgvsCodingChange", n).replaceFirst("&gt;", ">");
				else if (name.equals("chromosome")) chr = CarisXmlVcfParser.fetch("chromosome", n);
				else if (name.equals("molecularConsequence")) molecularConsequence = CarisXmlVcfParser.fetch("molecularConsequence", n);
				else if (name.equals("alterationDetails")) parseAlterationDetails(n.getChildNodes());
				else if (name.equals("alleleFrequencyInformation")) parseAlleleFrequencyInformation(n.getChildNodes());
				else if (name.equals("readInformation")) parseReadInformation(n.getChildNodes());
			}
		}
//IO.pl(toString());
	}
	
	private void parseReadInformation(NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("readDepth")) readDepth = CarisXmlVcfParser.fetch("readDepth", n);
				else throw new IOException("ERROR: found something other than 'readDepth' "+name);
			}
		}
	}
	
	private void parseAlleleFrequencyInformation(NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("alleleFrequency")) alleleFrequency = CarisXmlVcfParser.fetch("alleleFrequency", n);
				else throw new IOException("ERROR: found something other than 'alleleFrequency' "+name);
			}
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (biomarkerName!=null) sb.append("\nbiomarkerName\t"+biomarkerName);
		if (pathogenicity!=null) sb.append("\npathogenicity\t"+pathogenicity);
		if (geneName!=null) sb.append("\ngeneName\t"+geneName);
		if (hgvsCodingChange!=null) sb.append("\nhgvsCodingChange\t"+hgvsCodingChange);
		if (chr!=null) sb.append("\nchr\t"+chr);
		if (ref!=null) sb.append("\nref\t"+ref);
		if (alt!=null) sb.append("\nalt\t"+alt);
		if (transcriptId!=null) sb.append("\ntranscriptId\t"+transcriptId);
		if (molecularConsequence!=null) sb.append("\nmolecularConsequence\t"+molecularConsequence);
		if (alleleFrequency!=null) sb.append("\nalleleFrequency\t"+alleleFrequency);
		if (readDepth!=null) sb.append("\nreadDepth\t"+readDepth);
		return sb.toString();
	}

	private void parseAlterationDetails(NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("transcriptAlterationDetails")) parseTranscriptAlterationDetails(n.getChildNodes());
				else throw new IOException("ERROR: found something other than 'transcriptAlterationDetails' "+name);
			}
		}
	}



	private void parseTranscriptAlterationDetails(NodeList nodes) {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("transcriptID")) transcriptId = CarisXmlVcfParser.fetch("transcriptID", n);
				else if (name.equals("referenceNucleotide")) ref = CarisXmlVcfParser.fetch("referenceNucleotide", n);
				else if (name.equals("observedNucleotide")) alt = CarisXmlVcfParser.fetch("observedNucleotide", n);

			}
		}
	}


	public String fetchKey() {
		//chromosome_ref_alt_readDepth_geneName_hgvsCodingChange
		return chr +"_"+ref+"_"+alt+"_"+readDepth+"_"+geneName+ "_" +hgvsCodingChange; 
	}

	public String getBiomarkerName() {
		return biomarkerName;
	}

	public String getPathogenicity() {
		return pathogenicity;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getHgvsCodingChange() {
		return hgvsCodingChange;
	}

	public String getChr() {
		return chr;
	}

	public String getTranscriptId() {
		return transcriptId;
	}

	public String getMolecularConsequence() {
		return molecularConsequence;
	}

	public String getAlleleFrequency() {
		return alleleFrequency;
	}

	public String getReadDepth() {
		return readDepth;
	}
}
