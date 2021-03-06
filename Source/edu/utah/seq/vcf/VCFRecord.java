package edu.utah.seq.vcf;

import java.util.ArrayList;

import util.gen.Misc;


/**For parsing a VCFRecord with multiple samples*/
public class VCFRecord implements Comparable<VCFRecord> {
	
	//fields
	private String chromosome;
	private int position; //interbase coordinates! not 1 based
	private String rsNumber; //id
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
	
	private boolean missingQual = false;
	

	
	/**Only extracts some of the fields from a record*/
	public VCFRecord(String record, VCFParser vcfParser, boolean loadSamples, boolean loadInfo) throws Exception{
		originalRecord = record;
		String[] fields = vcfParser.TAB.split(record);
		if (vcfParser.numberFields !=0){
			if (fields.length != vcfParser.numberFields) throw new Exception("\nIncorrect number of fields in -> "+record+"\nTry uncompressing the vcf file?");
		}
		else if (fields.length < vcfParser.minimumNumberFields ) throw new Exception("\nIncorrect number of fields ("+vcfParser.minimumNumberFields+") in -> "+record+"\nTry uncompressing the vcf file?");
		else if (vcfParser.numberFields == 0) vcfParser.numberFields = fields.length;
		
		//must subtract 1 from position to put it into interbase coordinates
		chromosome = fields[vcfParser.chromosomeIndex] ;
		position = Integer.parseInt(fields[vcfParser.positionIndex]) - 1;
		
		String[] refByComma = VCFParser.COMMA.split(fields[vcfParser.referenceIndex]);
		String[] refBySlash = VCFParser.SLASH.split(fields[vcfParser.referenceIndex]);
		if (refBySlash.length > 1) {
			reference = refBySlash[0];
		} else if (refByComma.length > 1) {
			reference = refByComma[0];
		} else {
			reference = fields[vcfParser.referenceIndex];
		}
		
		String[] altByComma = VCFParser.COMMA.split(fields[vcfParser.alternateIndex]);
		String[] altBySlash = VCFParser.SLASH.split(fields[vcfParser.alternateIndex]);
		if (altByComma.length > altBySlash.length) {
			alternate = altByComma;
		} else {
			alternate = altBySlash;
		}
		
		rsNumber = fields[vcfParser.rsIndex];
		if (fields[vcfParser.qualityIndex].equals(".")) quality = 0;
		else quality = Float.parseFloat(fields[vcfParser.qualityIndex]);
		filter = fields[vcfParser.filterIndex];
		if (loadInfo){
			info = new VCFInfo();
			info.parseInfoGatk(fields[vcfParser.infoIndex]);
		}
		
		if (loadSamples){
			//watch out for missing sample info
			if (fields.length == 8) format = ".";
			else format = fields[vcfParser.formatIndex];
			ArrayList<VCFSample> al = new ArrayList<VCFSample>();
			//must watch out for blank samples
			for (int i=vcfParser.firstSampleIndex; i< fields.length; i++) {
//				if (fields[i].equals(".") == false) 
				al.add(new VCFSample(fields[i], format));
			}
			if (al.size()!= 0){
				sample = new VCFSample[al.size()];
				al.toArray(sample);
				
				//Check for missing qualities.  This is valid for things like varscan, so simply store the fact that there were missing quals.
				for (VCFSample s: sample) {
					if (s.isMissingQual()) {
						missingQual = true;
					}
				}
			}
			else sample = null;
		}
	}
	
	public VCFRecord() {}
	
	/**Returns the position + longest alt or ref, interbase coor.*/
	public int getMaxEndPosition(){
		int max = reference.length();
		for (String a: alternate) if (a.length() > max) max = a.length();
		return max+ position;
	}

	/**Return modified record line.*/
	public String getModifiedRecord(ArrayList<String> infoToUse, String style) {
		String infoLine;
		
		infoLine = getModifiedInfoString(infoToUse,style);
		
		String altString = "";
		for (String alt: alternate) {
			altString += "," + alt;
		}
		altString = altString.substring(1);
		
		String modifiedRecord = String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",chromosome,String.valueOf(position+1),
				rsNumber,reference,altString,String.valueOf(quality),filter,infoLine,format);
		for (VCFSample s: sample) {
			modifiedRecord += "\t" + s.getUnmodifiedSampleString();
		}
		return modifiedRecord;
	}
	
	
	
	/**Return original unmodified record line.*/
	public String toString(){
		return originalRecord;
	}
	
	/**Returns Chrom thru Info fields*/
	public String getTruncatedRecord(){
		StringBuilder sb = new StringBuilder();
		sb.append(chromosome); sb.append("\t");
		sb.append((position+1)); sb.append("\t");
		sb.append(rsNumber); sb.append("\t");
		sb.append(reference); sb.append("\t");
		sb.append(alternate[0]);
		for (int i=1; i<alternate.length; i++){
			sb.append(",");
			sb.append(alternate[i]);
		}
		sb.append("\t");
		sb.append(quality); sb.append("\t");
		sb.append(filter); sb.append("\t");
		sb.append(info.getInfoString());
		return sb.toString();
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
		if (match == false) {
			return false;
		}
		
		//check genotype of first sample
		if (requireGenotypeMatch){
			if (vcfRecord.getSample()[0].getGenotypeGT().equals(sample[0].getGenotypeGT()) == false) {
				return false;
			}
		}
		return true;
	}
	
	/**Checks if the position, ref, and any of the alternate alleles match.*/
	public boolean matchesPosRefAlt(VCFRecord vcfRecord) {
		if (this.position != vcfRecord.position) return false;
		if (this.reference.equals(vcfRecord.reference) == false) return false;
		
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
		return info;
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
	
	public boolean isMissingQual() {
		return this.missingQual;
	}
	
	/**Checks that all alternates are a snv.  TT->CC is an insertion, not a double snv.*/
	public boolean isSNP() {
		if (reference.length() != 1) return false;
		//check each alternate for non snp
		for (int i=0; i< alternate.length; i++){
			if (alternate[i].length() != 1 || alternate[i].equals(".")) return false;
		}
		return true;
	}
	/**Checks that each alternate is an insertion.
	 * TT -> CC is considered a balanced deletion/insertion and returns true here.*/
	public boolean isInsertion(){
		//check each alternate for non insertion
		for (int i=0; i< alternate.length; i++) if (alternate[i].length() < reference.length()) return false;
		return true;
	}
	/**Checks that each alternate is a deletion, smaller than ref.
	 * TT -> CC is considered an insertion and returns false here.*/
	public boolean isDeletion(){
		for (int i=0; i< alternate.length; i++) if (alternate[i].length() >= reference.length()) return false;
		return true;
	}
	
	/**Checks to see if any of the alternates matches the type (0 SNV, 1 INS, 2 DEL)*/
	public boolean matchesVariantType(int type){
		//snp
		if (type == 0){
			for (int i=0; i< alternate.length; i++){
				if (alternate[i].length() == 1 || alternate[i].equals(".") == false) return true;
			}
		}
		//insertion
		else if (type == 1){
			for (int i=0; i< alternate.length; i++){
				if (alternate[i].length() > reference.length()) return true;
			}
		}
		else {
			for (int i=0; i< alternate.length; i++){
				if (alternate[i].length() < reference.length()) return true;
			}
		}
		return false;
	}
	
	public int compareTo(VCFRecord other) {
		//sort by chromosome
		int x = chromosome.compareTo(other.chromosome);
		if (x !=0) return x;
		//sort by position
		if (position < other.position) return -1;
		if (position > other.position) return 1;
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
		String endPos = String.valueOf(getPosition() + 1 + maxInt(getReference().length()-1,0));
		StringBuffer altSb = new StringBuffer("");
		for (String alt: alternate) {
			altSb.append(","+alt);
		}
		String altS = altSb.toString().substring(1, altSb.length());
		StringBuffer full = new StringBuffer(getChromosome() + "\t" + getPosition() + "\t" + endPos + "\t" + getReference() + "\t" + 
							altS + "\t" + String.valueOf(getQuality()));
		full.append("\t" + info.buildInfoForTable(infoToAdd, style, comments));
		for (VCFSample sample: this.sample) {
			String genotype = sample.getGenotypeGT();
			if (sample.isNoCall()) {
				full.append("\tNA");
			} else {
				//Create an array list of all observed alleles.  This is used to create
				//a string-based genotype
				ArrayList<String> bases = new ArrayList<String>();
				bases.add(reference);
				for (String alt: alternate) {
					bases.add(alt);
				}
				
				//use the number-based genotype to look up the bases
				String[] genotypes = genotype.split("/");
				StringBuilder gtString = new StringBuilder("");
				for (String g: genotypes) {
					int index = Integer.parseInt(g);
					gtString.append("/"+bases.get(index));
				}
				
				//Create count string.  This contains the reference count + all alternate allele counts.
				//Some VCFs have RD and AD (alternate counts) columns, some just AD (all allele counts). 
				//If AD+1 = ref base + obs alt bases, add ref counts to string.
				String[] counts = sample.getAlleleCount().split(",");
				StringBuilder sb = new StringBuilder("");
				if (counts.length+1 == bases.size()) {
					sb.append("," + sample.getReferenceCount());
				}
				
				for (String c: counts) {
					sb.append(","+c);
				}
				
				full.append("\t" + gtString.toString().substring(1, gtString.length()) + ":" + genotype + ":" + sb.toString().substring(1, sb.length()));	
							
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

	/**Adds a chr onto the chromosome field if it doesn't start with such.*/
	public void appendChr() {
		if (chromosome.startsWith("chr") == false) chromosome = "chr"+chromosome;
		
	}
	public void correctChrMTs(){
		if (chromosome.equals("chrMT")) chromosome = "chrM";
		else if (chromosome.equals("MT")) chromosome = "M";
	}
	public String getChrPosRefAlt(boolean zeroPosition){
		int pos = position;
		if (zeroPosition == false) pos++;
		return chromosome+"_"+ pos +"_"+reference+"_"+Misc.stringArrayToString(alternate, ",");
	}

	public String getRsNumber() {
		return rsNumber;
	}

	public void setRsNumber(String rsNumber) {
		this.rsNumber = rsNumber;
	}

	public void appendId(String primaryName) {
		if (rsNumber.contains(primaryName)==false) {
			if (rsNumber.length() == 0) rsNumber = primaryName;
			else rsNumber = rsNumber+";"+primaryName;
		}
	}
	public void appendFilter(String otherFilter) {
		if (otherFilter.equals(".")) return;
		if (filter.contains(otherFilter)==false) {
			if (filter.length() == 0) filter = otherFilter;
			else filter = filter+";"+otherFilter;
		}
	}

	public int getSizeIndel() {
		int max = reference.length();
		for (String a: alternate) if (a.length() > max) max = a.length();
		return max;
	}

	/**For non snvs, checks to see that the first base in the ref and alt are the same.*/
	public boolean checkLeadingRef() {
		if (isSNP()) return true;
		char ref = reference.charAt(0);
		for (String alt: alternate){
			if (alt.charAt(0) != ref) return false;
		}
		return true;
	}

}
