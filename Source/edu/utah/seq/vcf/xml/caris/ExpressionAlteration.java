package edu.utah.seq.vcf.xml.caris;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import util.gen.IO;
import util.gen.Misc;

public class ExpressionAlteration {

	//fields
	private String expressionType = null;
	private String biomarkerName = null;
	private String result = null;
	private ArrayList<String> genes = new ArrayList<String>();
	private String stainPercent = null;
	private String threshold = null;
	private String isExpressed = null;
	private String equivocal = null;
	private String intensity = null;
	
	public ExpressionAlteration(LinkedHashMap<String, String> report, NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("biomarkerName")) biomarkerName = CarisXmlVcfParser.fetch("biomarkerName", n);
				else if (name.equals("expressionType")) expressionType = CarisXmlVcfParser.fetch("expressionType", n);
				else if (name.equals("result")) result = CarisXmlVcfParser.fetch("result", n);
				else if (name.equals("gene")) genes.add(CarisXmlVcfParser.fetch("gene", n));
				else if (name.equals("stainPercent")) stainPercent = CarisXmlVcfParser.fetch("stainPercent", n);
				else if (name.equals("threshold")) {
					threshold = CarisXmlVcfParser.fetch("threshold", n);
					if (threshold != null) {
						threshold = threshold.replace("&lt;", "<");
						threshold = threshold.replace("&gt;", ">");
					}
				}
				else if (name.equals("isExpressed")) isExpressed = CarisXmlVcfParser.fetch("isExpressed", n);
				else if (name.equals("equivocal")) equivocal = CarisXmlVcfParser.fetch("equivocal", n);
				else if (name.equals("intensity")) intensity = CarisXmlVcfParser.fetch("intensity", n);
			}
		}
		//IO.pl(toString());
		addReportInfo(report);
	}
	

	private void addReportInfo(LinkedHashMap<String, String> report) {
		String bn = biomarkerName.replaceAll(" ", "_");
		String res = result.replaceAll(" ", "_");
		report.put(bn+":Result", res);
		if (stainPercent != null) {
			String sp = stainPercent.replaceAll(" ", "_");
			report.put(bn+":StainPercent", sp);
		}
		if (intensity != null) {
			String in = intensity.replaceAll(" ", "_");
			report.put(bn+":Intensity", in);
		}
	}


	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (expressionType!=null) sb.append("\nexpressionType\t"+expressionType);
		if (biomarkerName!=null) sb.append("\nbiomarkerName\t"+biomarkerName);
		if (result!=null) sb.append("\nresult\t"+result);
		if (genes.size()!=0) sb.append("\ngene\t"+Misc.stringArrayListToString(genes, ", "));
		if (stainPercent!=null) sb.append("\nstainPercent\t"+stainPercent);
		if (intensity!=null) sb.append("\nintensity\t"+intensity);
		if (threshold!=null) sb.append("\nthreshold\t"+threshold);
		if (isExpressed!=null) sb.append("\nisExpressed\t"+isExpressed);
		if (equivocal!=null) sb.append("\nequivocal\t"+equivocal);
		return sb.toString();
	}
}
