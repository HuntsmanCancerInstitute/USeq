package edu.utah.seq.vcf;

import java.util.ArrayList;


/**For parsing a VCFRecord with multiple samples*/
public class VCFRecord implements Comparable<VCFRecord> {
	
	//fields
	private String chromosome;
	private int position; //interbase coordinates! not 1 based
	private String rsNumber;
	private String reference;
	private String[] alternate;
	private float quality;
	private String filter;
	private VCFInfo info;
	private String format;
	private VCFSample[] sample;
	private String originalRecord;
	public static final String PASS = "PASS";
	public static final String FAIL = "FAIL";
	private float score = 0;
	
	/**Only extracts some of the fields from a record*/
	public VCFRecord(String record, VCFParser vcfParser, boolean loadSamples, boolean loadInfo) throws Exception{
		originalRecord = record;
		String[] fields = vcfParser.TAB.split(record);
		if (vcfParser.numberFields !=0){
			if (fields.length != vcfParser.numberFields) throw new Exception("\nIncorrect number of fields in -> "+record+"\nTry uncompressing the vcf file?");
		}
		else if (fields.length < vcfParser.minimumNumberFields ) throw new Exception("\nIncorrect number of fields in -> "+record+"\nTry uncompressing the vcf file?");
		else if (vcfParser.numberFields == 0) vcfParser.numberFields = fields.length;
		
		//must subtract 1 from position to put it into interbase coordinates
		chromosome = fields[vcfParser.chromosomeIndex] ;
		position = Integer.parseInt(fields[vcfParser.positionIndex]) - 1;
		reference = fields[vcfParser.referenceIndex];
		alternate = VCFParser.COMMA.split(fields[vcfParser.alternateIndex]);
		rsNumber = fields[vcfParser.rsIndex];
		if (fields[vcfParser.qualityIndex].equals(".")) quality = 0;
		else quality = Float.parseFloat(fields[vcfParser.qualityIndex]);
		filter = fields[vcfParser.filterIndex];
		if (loadInfo){
			info = new VCFInfo();
			info.parseInfoGatk(fields[vcfParser.infoIndex]);
		}
		format = fields[vcfParser.formatIndex];
		
		if (loadSamples){
			ArrayList<VCFSample> al = new ArrayList<VCFSample>();
			//must watch out for blank samples
			for (int i=vcfParser.firstSampleIndex; i< fields.length; i++) {
				if (fields[i].equals(".") == false) al.add(new VCFSample(fields[i], format));
			}
			if (al.size()!= 0){
				sample = new VCFSample[al.size()];
				al.toArray(sample);
			}
			else sample = null;
		}
	}
	
	/**Return modified record line.*/
	public String getModifiedRecord(ArrayList<String> infoToUse, String style) {
		String infoLine;
		
		infoLine = this.getModifiedInfoString(infoToUse,style);
		
		String altString = "";
		for (String alt: this.alternate) {
			altString += "," + alt;
		}
		altString = altString.substring(1);
		
		String modifiedRecord = String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",this.chromosome,String.valueOf(this.position+1),
				this.rsNumber,this.reference,altString,String.valueOf(this.quality),this.filter,infoLine,this.format);
		for (VCFSample s: this.sample) {
			modifiedRecord += "\t" + s.getUnmodifiedSampleString();
		}
		return modifiedRecord;
	}
	
	
	
	/**Return original unmodified record line.*/
	public String toString(){
		return originalRecord;
	}
	
	/**Checks if any of the alternate alleles match and optionally that the genotype of the first sample are identical.*/
	public boolean matchesAlternateAlleleGenotype(VCFRecord vcfRecord, boolean requireGenotypeMatch) {
		//check alternate allele
		String[] otherAlts = vcfRecord.getAlternate();
		boolean match = false;
		for (String o : otherAlts){
			for (String t : alternate){
				if (o.equals(t)){
					match = true;
					break;
				}
			}
		}
		if (match == false) return false;
		
		//check genotype of first sample
		if (requireGenotypeMatch){
			if (vcfRecord.getSample()[0].getGenotypeGT().equals(sample[0].getGenotypeGT()) == false) return false;
		}
		return true;
	}

	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public String getReference() {
		return reference;
	}

	public void setReference(String reference) {
		this.reference = reference;
	}

	public String[] getAlternate() {
		return alternate;
	}

	public void setAlternate(String[] alternate) {
		this.alternate = alternate;
	}

	public String getFilter() {
		return filter;
	}

	public void setFilter(String filter) {
		this.filter = filter;
	}
	
	/** get subset of modified info fields */
	public String getModifiedInfoString(ArrayList<String> infoToUse) {
		return info.buildInfoString(infoToUse, VCFInfo.UNMODIFIED);
	}
	
	/** get subset of modified info fields in a specific style */
	public String getModifiedInfoString(ArrayList<String> infoToUse, String style) {
		return info.buildInfoString(infoToUse, style);
	}
	
	/** get raw, unmodified info fields */
	public String getUnmodifiedInfoString() {
		return info.getInfoString();
	}

	public void setInfoString(String info) {
		this.info.overwriteInfoString(info);
	}
	public float getQuality() {
		return quality;
	}
	
	public VCFInfo getInfoObject() {
		return this.info;
	}

	public VCFSample[] getSample() {
		return sample;
	}

	public void setSample(VCFSample[] sample) {
		this.sample = sample;
	}

	public float getScore() {
		return score;
	}

	public void setScore(float score) {
		this.score = score;
	}

	public String getFormat() {
		return format;
	}

	public void setFormat(String format) {
		this.format = format;
	}

	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public String getOriginalRecord() {
		return originalRecord;
	}

	public void setOriginalRecord(String originalRecord) {
		this.originalRecord = originalRecord;
	}

	public void setQuality(float quality) {
		this.quality = quality;
	}

	public boolean isSNP() {
		if (alternate[0].length() == 1 && reference.length() == 1 && alternate.equals(".") == false) return true;
		return false;
	}
	
	public int compareTo(VCFRecord other) {
		//sort by chromosome
		int x = this.chromosome.compareTo(other.chromosome);
		if (x !=0) return x;
		//sort by position
		if (this.position < other.position) return -1;
		if (this.position > other.position) return 1;
		return 0;
	}
	
	private Integer maxInt(Integer val1, Integer val2) {
		if (val1 >= val2) {
			return val1;
		} else {
			return val2;
		}
	}
	
	public String getSpreadsheetOutput(ArrayList<String> infoToAdd, String style, VCFComments comments) {
		String endPos = String.valueOf(this.getPosition() + 1 + this.maxInt(this.getReference().length()-1,0));
		StringBuffer full = new StringBuffer(this.getChromosome() + "\t" + this.getPosition() + "\t" + endPos + "\t" + this.getReference() + "\t" + 
							this.getAlternate()[0] + "\t" + String.valueOf(this.getQuality()));
		full.append("\t" + this.info.buildInfoForTable(infoToAdd, style, comments));
		for (VCFSample sample: this.sample) {
			String genotype = sample.getGenotypeGT();
			if (sample.isNoCall()) {
				full.append("\tNA");
			} else {
				full.append("\tGT:" + genotype + ":" + sample.getAlleleCount());
			}
		}
		full.append("\n");
		return full.toString();
	}

	/**Looks to see if any of the alternates are an A <-> G or T <-> C change. 
	 * If not returns false.  Note all non SNPs will be false. 
	 * So you might want to check it first.*/
	public boolean isTransition() {
		if (reference.equals("A")){
			for (String base : alternate){
				if (base.equals("G")) return true;
			}
		}
		else if (reference.equals("T")){
			for (String base : alternate){
				if (base.equals("C")) return true;
			}
		}
		else if (reference.equals("G")){
			for (String base : alternate){
				if (base.equals("A")) return true;
			}
		}
		else if (reference.equals("C")){
			for (String base : alternate){
				if (base.equals("T")) return true;
			}
		}
		return false;
	}	


}
