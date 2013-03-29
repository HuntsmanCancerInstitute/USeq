package edu.utah.seq.vcf;


/* GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> 0/0, 0/1, 1/1
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"> comma delimited
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality"> bigger numbers better!
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
 */
public class VCFSample {

	//fields, default to null or -1
	private String genotypeGT = null;
	private int readDepthDP = -1;
	private int genotypeQualityGQ = -1;
	private boolean noCall = false;
	private String originalRecord = null;
	private String originalFormat = null;
	private String alleleCounts = null;

	/**Finding vcf files with mixed sample formats so must determine each record by record :( 
	 * Add rippers as needed.*/
	public VCFSample(String sample, String sampleFormat) throws Exception{
		this.originalRecord = sample;
		this.originalFormat = sampleFormat;
		//is it a no call?
		if (sample.equals("./.")){
			noCall = true;
		}
		else {
			String[] data = VCFParser.COLON.split(sample);
			String[] format = VCFParser.COLON.split(sampleFormat);
			if (data.length != format.length) throw new Exception("Incorrect number of fields in sample -> "+sample+" for indicated format -> "+sampleFormat);
			//attempt to parse GT, DP, GQ
			for (int i=0; i< format.length; i++){
				if (format[i].equals("GT")) genotypeGT = data[i];
				else if (format[i].equals("DP")) readDepthDP = Integer.parseInt(data[i]);
				else if (format[i].equals("GQ")) genotypeQualityGQ = Integer.parseInt(data[i]);
				else if (format[i].equals("AD")) {
					alleleCounts = data[i];
				}
			}
			
		}
	}
	
	public String getAlleleCount() {
		return alleleCounts;
	}


	public String getGenotypeGT() {
		return genotypeGT;
	}

	public int getReadDepthDP() {
		return readDepthDP;
	}

	public int getGenotypeQualityGQ() {
		return genotypeQualityGQ;
	}

	public boolean isNoCall() {
		return noCall;
	}

	public void setNoCall(boolean noCall) {
		this.noCall = noCall;
	}
	
	public String getUnmodifiedSampleString() {
		return this.originalRecord;
	}
}
