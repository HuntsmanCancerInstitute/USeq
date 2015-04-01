package edu.utah.seq.vcf;

import java.util.regex.Pattern;


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
	private float genotypeQualityGQ = -1;
	private boolean noCall = false;
	private String originalRecord = null;
	private String originalFormat = null;
	private String alleleCounts = null;
	private String referenceCounts = null;
    private String alternateCounts = null;
    private String[] data = null;
    private String[] format = null;
    private boolean missingQual = false;
	
	private static final Pattern PIPE = Pattern.compile("\\|");
	

	/**Finding vcf files with mixed sample formats so must determine each record by record :( 
	 * Add rippers as needed.*/
	public VCFSample(String sample, String sampleFormat) throws Exception{
		this.originalRecord = sample;
		this.originalFormat = sampleFormat;
		//is it a no call?
		if (sample.equals(".") || sample.substring(0,3).startsWith("./.") ){
			noCall = true;
		}
		else {
			data = VCFParser.COLON.split(sample);
			format = VCFParser.COLON.split(sampleFormat);
			if (data.length != format.length) throw new Exception("Incorrect number of fields in sample -> "+sample+" for indicated format -> "+sampleFormat);
			//attempt to parse GT, DP, GQ
			for (int i=0; i< format.length; i++){
				if (format[i].equals("GT")) {
					genotypeGT = data[i];
					//replace any | with /
					genotypeGT = PIPE.matcher(genotypeGT).replaceAll("/");
					//replace 1/0 with 0/1
					if (genotypeGT.equals("1/0")) genotypeGT = "0/1";
				}
				else if (format[i].equals("DP")) {
					if (data[i].equals(".")) {
						noCall = true;
						break;
					} else {
						readDepthDP = Integer.parseInt(data[i]);
					}
				}
				else if (format[i].equals("GQ")) {
					if (data[i].equals(".")) {
						genotypeQualityGQ = 0;
						missingQual = true;
					} else {
						genotypeQualityGQ = Float.parseFloat(data[i]);
					}
					
				}
				else if (format[i].equals("AD")) {
					alleleCounts = data[i];
					String[] multiple = data[i].split(",");
					if (multiple.length > 1) {
						if (referenceCounts == null) {
							this.referenceCounts = multiple[0];
						}
						this.alternateCounts = multiple[1];
					} else if (multiple.length ==  1) {
						this.alternateCounts = multiple[0];
					} else {
						throw new Exception("There is no data in the 'AD' field, please make sure your file is formed correctly.  If it is formed correctly, "
								+ "please harass the HCI core to fix this problem");
					}
				} else if (format[i].equals("RD")) {
					this.referenceCounts = data[i];
				}
				else if (format[i].equals("DP4")) {
					//# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases
					String[] multiple = data[i].split(",");
					int rf = Integer.parseInt(multiple[0]);
					int rr = Integer.parseInt(multiple[1]);
					int af = Integer.parseInt(multiple[2]);
					int ar = Integer.parseInt(multiple[3]);
					referenceCounts = (rf+rr)+"";
					alternateCounts = (af+ar)+"";
				}
			}
			//not all vcf records have counts! so don't require it
			/*if (referenceCounts == null || alternateCounts == null) {
				throw new Exception("Could not parse counts for either the reference (" + referenceCounts + ") or alternate (" + alternateCounts + ") allele.");
			}*/
			
		}
	}
	
	/**Looks in the Format for the label and returns the Sample value.*/
	public String getFormatData(String label){
		for (int i=0; i< format.length; i++){
			if (format[i].equals(label)) return data[i];
		}
		return null;
	}
	
	/**Returns -1 if either ref or alt counts are null.*/
	public double getAltRatio(){
		if (referenceCounts == null || alternateCounts == null) return -1;
		double ref = Double.parseDouble(referenceCounts);
		double alt = Double.parseDouble(alternateCounts);
		return alt/(alt+ref);
	}
	
	public String getReferenceCount() {
		return referenceCounts;
	}
	
	public String getAlternateCounts() {
		return alternateCounts;
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

	public float getGenotypeQualityGQ() {
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
	
	public boolean isMissingQual() {
		return this.missingQual;
	}

	public String[] getData() {
		return data;
	}

	public String[] getFormat() {
		return format;
	}

	public void setReferenceCounts(String referenceCounts) {
		this.referenceCounts = referenceCounts;
	}

	public void setAlternateCounts(String alternateCounts) {
		this.alternateCounts = alternateCounts;
	}
}
