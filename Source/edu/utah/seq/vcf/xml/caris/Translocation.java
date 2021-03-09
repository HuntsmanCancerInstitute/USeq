package edu.utah.seq.vcf.xml.caris;

import java.io.IOException;
import java.util.HashMap;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import util.bio.parsers.UCSCGeneLine;
import util.gen.IO;
import util.gen.Misc;

public class Translocation {

	//fields
	private String biomarkerName = null;
	private String result = null;
	private String fusionIsoform = null;
	private String geneName1 = null;
	private String geneName2 = null;
	private String genomeVersion = null;
	private String chrGene1 = null;
	private String chrGene2 = null;
	private int startGene1 = -1;
	private int startGene2 = -1;


	public Translocation(NodeList nodes) throws IOException {
		for (int j = 0; j < nodes.getLength(); j++) {
			Node n = nodes.item(j);
			if (n instanceof Element) {
				String name = n.getNodeName();
				if (name.equals("biomarkerName")) biomarkerName = CarisXmlVcfParser.fetch("biomarkerName", n);
				else if (name.equals("result")) result = CarisXmlVcfParser.fetch("result", n);
				else if (name.equals("genomeBuild")) {
					//these are mixed, undefined is likely hg19, also seeing hg38
					genomeVersion = CarisXmlVcfParser.fetch("genomeBuild", n);
				}
				else if (name.equals("fusionISOForm")) fusionIsoform = CarisXmlVcfParser.fetch("fusionISOForm", n);
				else if (name.equals("gene1")) geneName1 = CarisXmlVcfParser.fetch("gene1", n);
				else if (name.equals("gene2")) geneName2 = CarisXmlVcfParser.fetch("gene2", n);
				else if (name.equals("genomicBreakpoint")) {
					//   chr2:42522656:+/chr2:29446394:-   hg19 coordinates!! others hg38
					String bp = CarisXmlVcfParser.fetch("genomicBreakpoint", n);
					String[] f = Misc.FORWARD_SLASH.split(bp);
					if (f.length!=2) throw new IOException("ERROR: failed to split two / delimited fields from "+bp);
					String[] ccg1 = Misc.COLON.split(f[0]);
					chrGene1 = ccg1[0];
					startGene1 = Integer.parseInt(ccg1[1]);
					String[] ccg2 = Misc.COLON.split(f[1]);
					chrGene2 = ccg2[0];
					startGene2 = Integer.parseInt(ccg2[1]);
				}
			}
		}
//IO.pl(toString());
	}
	

	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (biomarkerName!=null) sb.append("\nbiomarkerName\t"+biomarkerName);
		if (result!=null) sb.append("\nresult\t"+result);
		if (fusionIsoform!=null) sb.append("\nfusionIsoform\t"+fusionIsoform);
		if (geneName1!=null) sb.append("\ngeneName1\t"+geneName1);
		if (geneName2!=null) sb.append("\ngeneName2\t"+geneName2);
		if (chrGene1!=null) sb.append("\nchrGene1\t"+chrGene1);
		if (chrGene2!=null) sb.append("\nchrGene2\t"+chrGene2);
		if (genomeVersion!=null) sb.append("\ngenomeVersion\t"+genomeVersion);
		if (startGene1!=-1) sb.append("\nstartGene1\t"+startGene1);
		if (startGene2!=-1) sb.append("\nstartGene2\t"+startGene2);
		return sb.toString();
	}


	/**Using gene model lookup instead of parsed coordinates since these are mixed hg19/hg38.
	 * @throws IOException */
	public SimpleBed[] toBed(HashMap<String, UCSCGeneLine[]> name2GeneModels) throws IOException {
		SimpleBed bed1 = null;
		SimpleBed bed2 = null;
		
		//Any gene names?
		if (geneName1!=null || geneName2!=null || chrGene1!=null || chrGene2!=null) {
			UCSCGeneLine gene1 = fetchModel(geneName1, chrGene1, name2GeneModels);
			UCSCGeneLine gene2 = fetchModel(geneName2, chrGene2, name2GeneModels);
			String name = gene1.getDisplayName()+":"+ gene2.getDisplayName()+ ":"+ result.toLowerCase().replace(" ", "_");
			bed1 = new SimpleBed(gene1.getChrom(), gene1.getTxStart(), gene1.getTxEnd(), name, null, gene1.getStrand());
			bed2 = new SimpleBed(gene2.getChrom(), gene2.getTxStart(), gene2.getTxEnd(), name, null, gene2.getStrand());
			
		}
		else if (biomarkerName != null) {
			String sGeneName = biomarkerName;
			//cMet = MET, ARv7 = AR
			if (sGeneName.endsWith("cMET")) sGeneName = "MET";
			else if (sGeneName.endsWith("ARv7")) sGeneName = "AR";
			UCSCGeneLine gene1 = fetchModel(sGeneName, null, name2GeneModels);
			String name = gene1.getDisplayName()+ ":none:"+ result.toLowerCase().replace(" ", "_");
			bed1 = new SimpleBed(gene1.getChrom(), gene1.getTxStart(), gene1.getTxEnd(), name, null, gene1.getStrand());
		}
		
		return new SimpleBed[] {bed1, bed2};
	}


	private UCSCGeneLine fetchModel(String geneName, String chr, HashMap<String, UCSCGeneLine[]> name2GeneModels) throws IOException {
		UCSCGeneLine[] genes = name2GeneModels.get(geneName);
		if (genes == null) throw new IOException("ERROR: failed to find a gene model for fusion coordinate extraction for "+geneName);
		UCSCGeneLine model = null;
		if (chr != null) {
			for (UCSCGeneLine line : genes) {
				if (line.getChrom().equals(chr)) {
					model = line;
					break;
				}
			}
			if (model == null) throw new IOException("ERROR: failed to find a chrom specific gene model for fusion coordinate extraction for "+geneName);
		}
		else model = genes[0];
		return model;
	}
	


}
