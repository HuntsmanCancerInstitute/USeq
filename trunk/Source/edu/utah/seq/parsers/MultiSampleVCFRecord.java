package edu.utah.seq.parsers;

/**For parsing a VCFRecord with multiple samples*/
public class MultiSampleVCFRecord {
	
	//fields
	private int position; //interbase coordinates! not 1 based
	private String reference;
	private String alternate;
	private String quality;
	private String filter;
	private String info;
	private String format;
	private VCFSample[] sample;
	private String originalRecord;
	
	/**Only extracts some of the fields from a record*/
	public MultiSampleVCFRecord(String record, MultiSampleVCFParser vcfParser) throws Exception{
		originalRecord = record;
		String[] fields = vcfParser.TAB.split(record);
		if (vcfParser.numberFields !=0){
			if (fields.length != vcfParser.numberFields) throw new Exception("\nIncorrect number of fields in -> "+record);
		}
		else if (fields.length < vcfParser.minimumNumberFields ) throw new Exception("\nIncorrect number of fields in -> "+record);
		else if (vcfParser.numberFields == 0) vcfParser.numberFields = fields.length;
		
		//must subtract 1 from position to put it into interbase coordinates
		position = Integer.parseInt(fields[vcfParser.positionIndex]) - 1;
		reference = fields[vcfParser.referenceIndex];
		alternate = fields[vcfParser.alternateIndex];
		quality = fields[vcfParser.qualityIndex];
		filter = fields[vcfParser.filterIndex];
		info = fields[vcfParser.infoIndex];
		format = fields[vcfParser.formatIndex];
		
		//check format
		if (format.equals(vcfParser.expectedSampleFormat) == false) throw new Exception("\nSample format does not match expected format ('"+vcfParser.expectedSampleFormat+ "') for record -> "+record);
		
		//make samples
		sample = new VCFSample[fields.length - vcfParser.firstSampleIndex];
		int index = 0;
		for (int i=vcfParser.firstSampleIndex; i< fields.length; i++) sample[index++] = new VCFSample(fields[i], vcfParser);
	}
	
	public String toString(){
		return originalRecord;
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
	public String getQuality() {
		return quality;
	}

	public VCFSample[] getSample() {
		return sample;
	}

	public void setSample(VCFSample[] sample) {
		this.sample = sample;
	}


}
