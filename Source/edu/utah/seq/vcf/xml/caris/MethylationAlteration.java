package edu.utah.seq.vcf.xml.caris;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import util.gen.IO;
import util.gen.Misc;

public class MethylationAlteration {

	//fields
	private String biomarkerName = null;
	private String gene = null;
	private String result = null;
	
	public MethylationAlteration(LinkedHashMap<String, String> report, NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("biomarkerName")) biomarkerName = CarisXmlVcfParser.fetch("biomarkerName", n);
				else if (name.equals("gene")) gene = CarisXmlVcfParser.fetch("gene", n);
				else if (name.equals("result")) result = CarisXmlVcfParser.fetch("result", n);
			}
		}
		//IO.pl(toString());
		addReportInfo(report);
	}

	private void addReportInfo(LinkedHashMap<String, String> report) {
		String bn = biomarkerName.replaceAll(" ", "_");
		String res = result.replaceAll(" ", "_");
		report.put(bn+":Result", res);
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (biomarkerName!=null) sb.append("\nbiomarkerName\t"+biomarkerName);
		if (gene!=null) sb.append("\ngene\t"+gene);
		if (result!=null) sb.append("\nresult\t"+result);
		return sb.toString();
	}
}
