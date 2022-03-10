package edu.utah.seq.vcf.xml.caris;

import java.io.IOException;
import java.util.HashMap;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import util.bio.annotation.Bed;
import util.bio.parsers.UCSCGeneLine;
import util.gen.IO;
import util.gen.Misc;

public class CNVAlteration {

	//fields
	private String biomarkerName = null;
	private String result = null;
	private String geneName = null;
	//hg38 and hg19
	private String genomeVersion = null;
	private String chr = null;
	private int start = -1;
	private int end = -1;
	private String copyNumberType = null;
	private String copyNumberScore = null;
	
	//for CISH
	private String copyNumberControl = null;
	private String copyNumberCounted = null;
	private String copyNumberRatio = null;
	private String threshold = null;

	public CNVAlteration(NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("biomarkerName")) biomarkerName = CarisXmlVcfParser.fetch("biomarkerName", n);
				else if (name.equals("gene")) geneName = CarisXmlVcfParser.fetch("gene", n);
				else if (name.equals("result")) result = CarisXmlVcfParser.fetch("result", n);
				else if (name.equals("chromosome")) chr = CarisXmlVcfParser.fetch("chromosome", n);
				else if (name.equals("genomicCoordinates")) {
					String coor = CarisXmlVcfParser.fetch("genomicCoordinates", n);
					String[] f = Misc.COLON.split(coor);
					if (f.length!=3) throw new IOException("ERROR: failed to split three : delimited fields from "+coor);
					String[] ss = Misc.DASH.split(f[2]);
					start = Integer.parseInt(ss[0]);
					end = Integer.parseInt(ss[1]);
				}
				else if (name.equals("genomeBuild")) genomeVersion = CarisXmlVcfParser.fetch("genomeBuild", n);
				else if (name.equals("copyNumberType")) copyNumberType = CarisXmlVcfParser.fetch("copyNumberType", n);
				else if (name.equals("copyNumber")) copyNumberScore = CarisXmlVcfParser.fetch("copyNumber", n);
				//for CISH
				else if (name.equals("copyNumberControl")) copyNumberControl = CarisXmlVcfParser.fetch("copyNumberControl", n);
				else if (name.equals("copyNumberCounted")) copyNumberCounted = CarisXmlVcfParser.fetch("copyNumberCounted", n);
				else if (name.equals("copyNumberRatio")) copyNumberRatio = CarisXmlVcfParser.fetch("copyNumberRatio", n);
				else if (name.equals("threshold")) {
					threshold = CarisXmlVcfParser.fetch("threshold", n);
					if (threshold != null) {
						threshold = threshold.replace("&lt;", "<");
						threshold = threshold.replace("&gt;", ">");
					}
				}
			}
		}
//IO.pl(toString());
	}
	

	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (biomarkerName!=null) sb.append("\nbiomarkerName\t"+biomarkerName);
		if (geneName!=null) sb.append("\ngeneName\t"+geneName);
		if (result!=null) sb.append("\nresult\t"+result);
		if (genomeVersion!=null) sb.append("\ngenomeVersion\t"+genomeVersion);
		if (chr!=null) {
			sb.append("\nchr\t"+chr);
			sb.append("\nstart\t"+start);
			sb.append("\nend\t"+end);
		}
		if (copyNumberScore!=null) sb.append("\ncopyNumberScore\t"+copyNumberScore);
		if (copyNumberType!=null) sb.append("\ncopyNumberType\t"+copyNumberType);
		if (copyNumberControl!=null) sb.append("\ncopyNumberControl\t"+copyNumberControl);
		if (copyNumberCounted!=null) sb.append("\ncopyNumberCounted\t"+copyNumberCounted);
		if (copyNumberRatio!=null) sb.append("\ncopyNumberRatio\t"+copyNumberRatio);
		if (threshold!=null) sb.append("\nthreshold\t"+threshold);
		return sb.toString();
	}

	/**Using gene model lookup instead of parsed coordinates since these are mixed hg19/hg38.
	 * @throws IOException */
	public SimpleBed toBed(HashMap<String, UCSCGeneLine[]> name2GeneModels) throws IOException {
		// Chr\tGeneStart\tGeneStop\tGeneName:Type\tScore\tGeneStrand
		UCSCGeneLine[] genes = name2GeneModels.get(geneName);
		if (genes == null) throw new IOException("ERROR: failed to find a gene model for cnv coordinate extraction for "+geneName+"\n"+toString());
		UCSCGeneLine model = null;
		if (chr != null) {
			for (UCSCGeneLine line : genes) {
				if (line.getChrom().equals(chr)) {
					model = line;
					break;
				}
			}
			if (model == null) throw new IOException("ERROR: failed to find a chrom specific gene model for cnv coordinate extraction for "+geneName);
		}
		else model = genes[0];
		String score = copyNumberScore;
		if (score == null) score = ".";
		return new SimpleBed(model.getChrom(), model.getTxStart(), model.getTxEnd(), geneName+":"+result.toLowerCase().replace(" ", "_"), score, model.getStrand());
	}


	public String getResult() {
		return result;
	}
}
