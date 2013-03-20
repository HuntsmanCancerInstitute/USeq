package edu.utah.seq.vcf;


/* GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> 0/0, 0/1, 1/1
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"> comma delimited
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality"> bigger numbers better!
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
 */
public class VCFSample {

	//fields
	private String genotypeGT;
	private int readDepthDP;
	private int genotypeQualityGQ;
	private boolean noCall = false;

	public VCFSample(String sample, VCFParser parser) throws Exception{
		//is it a no call?
		if (sample.equals("./.")){
			noCall = true;
		}
		else {
			String[] fields = parser.COLON.split(sample);
			if (fields.length != parser.numberFieldsInSample) throw new Exception("Incorrect number of fields in sample -> "+sample);
			if (parser.sampleGenotypeGTIndex !=-1) genotypeGT = fields[parser.sampleGenotypeGTIndex];
			if (parser.sampleReadDepthDPIndex !=-1) readDepthDP = Integer.parseInt(fields[parser.sampleReadDepthDPIndex]);
			if (parser.sampleGenotypeQualityGQIndex !=-1) genotypeQualityGQ = Integer.parseInt(fields[parser.sampleGenotypeQualityGQIndex]);
		}
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

}
