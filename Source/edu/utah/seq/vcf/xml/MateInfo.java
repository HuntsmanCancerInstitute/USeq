package edu.utah.seq.vcf.xml;

import java.util.regex.Matcher;
import htsjdk.samtools.reference.ReferenceSequence;

public class MateInfo {

	//fields
	private String coordinates;
	private String chromosome;
	private int position;
	private String reference;
	private int end;
	private String id; 
	private String gene;
	private boolean failedParsing = false;
	private FoundationXml2Vcf parser;
	
	public MateInfo(String coordinates, String id, String gene, FoundationXml2Vcf parser){
		this.coordinates = coordinates;
		this.id = id;
		this.gene = gene;
		this.parser = parser;
		
		parseChromPos();
		if (failedParsing == false) parseRef();
	}
	
	private void parseRef() {
		//fetch the first base
		ReferenceSequence p = parser.getFasta().getSubsequenceAt(chromosome, position, position);
		reference = new String(p.getBases());
	}
	
	private void parseChromPos() {
		//parse  coordinates
		Matcher mat = FoundationRearrangeVariant.coorPat.matcher(coordinates);
		if (mat.matches()){
			chromosome = mat.group(1);
			position = Integer.parseInt(mat.group(2));
			end = Integer.parseInt(mat.group(3));
		}
		else {
			System.err.println("\tError in MateInfo: failed to match the chr:pos-end for '"+coordinates+ "'");
			failedParsing = true;
		}
	}

	public String getPartialVcf(MateInfo other) {
		StringBuilder sb = new StringBuilder();
		//CHROM
		sb.append(chromosome); sb.append("\t");
		//POS
		sb.append(position); sb.append("\t");
		//ID
		sb.append(id); sb.append("\t");
		//REF
		sb.append(reference); 
		//ALT
		sb.append("\t<BND>\t");
		//QUAL
		sb.append(".\t");
		//FILTER
		sb.append(".\t");
		//INFO
		sb.append("EG=");
		sb.append(gene);
		sb.append(";MG=");
		sb.append(other.gene);
		sb.append(";MATEID=");
		sb.append(other.id);
		sb.append(";MC=");
		sb.append(other.coordinates);
		sb.append(";END=");
		sb.append(end);
		return sb.toString();
	}

	public boolean isFailedParsing() {
		return failedParsing;
	}
}
