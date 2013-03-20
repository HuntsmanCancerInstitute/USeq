package edu.utah.seq.vcf;


/**For parsing a VCFRecord with multiple samples*/
public class VCFRecord {
	
	//fields
	private String chromosome;
	private int position; //interbase coordinates! not 1 based
	private String reference;
	private String alternate;
	private float quality;
	private String filter;
	private String info;
	private String format;
	private VCFSample[] sample;
	private String originalRecord;
	public static final String PASS = "PASS";
	public static final String FAIL = "FAIL";
	public float score = 0;
	
	/**Only extracts some of the fields from a record*/
	public VCFRecord(String record, VCFParser vcfParser, boolean loadSamples) throws Exception{
		originalRecord = record;
		String[] fields = vcfParser.TAB.split(record);
		if (vcfParser.numberFields !=0){
			if (fields.length != vcfParser.numberFields) throw new Exception("\nIncorrect number of fields in -> "+record);
		}
		else if (fields.length < vcfParser.minimumNumberFields ) throw new Exception("\nIncorrect number of fields in -> "+record);
		else if (vcfParser.numberFields == 0) vcfParser.numberFields = fields.length;
		
		//must subtract 1 from position to put it into interbase coordinates
		chromosome = fields[vcfParser.chromosomeIndex] ;
		position = Integer.parseInt(fields[vcfParser.positionIndex]) - 1;
		reference = fields[vcfParser.referenceIndex];
		alternate = fields[vcfParser.alternateIndex];
		quality = Float.parseFloat(fields[vcfParser.qualityIndex]);
		filter = fields[vcfParser.filterIndex];
		info = fields[vcfParser.infoIndex];
		format = fields[vcfParser.formatIndex];
		
		if (loadSamples){
			//check format
			if (format.equals(vcfParser.expectedSampleFormat) == false) throw new Exception("\nSample format does not match expected format ('"+vcfParser.expectedSampleFormat+ "') for record -> "+record);
			//make samples
			sample = new VCFSample[fields.length - vcfParser.firstSampleIndex];
			int index = 0;
			for (int i=vcfParser.firstSampleIndex; i< fields.length; i++) sample[index++] = new VCFSample(fields[i], vcfParser);
		}
	}
	
	/**Return original unmodified record line.*/
	public String toString(){
		return originalRecord;
	}
	
	/**Checks that alternate allele and optionally that the genotype of the first sample are identical.*/
	public boolean matchesAlternateAlleleGenotype(VCFRecord vcfRecord, boolean requireGenotypeMatch) {
		//check alternate allele
		if (vcfRecord.getAlternate().equals(alternate) == false) return false;
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

	public String getAlternate() {
		return alternate;
	}

	public void setAlternate(String alternate) {
		this.alternate = alternate;
	}

	public String getFilter() {
		return filter;
	}

	public void setFilter(String filter) {
		this.filter = filter;
	}

	public String getInfo() {
		return info;
	}

	public void setInfo(String info) {
		this.info = info;
	}
	public float getQuality() {
		return quality;
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


}
