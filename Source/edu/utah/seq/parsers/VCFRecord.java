package edu.utah.seq.parsers;

import java.util.regex.Pattern;

/**For parsing a VCFRecord
 * 
##fileformat=VCFv4.1
##source=VarScan2
##INFO=<ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 15">
##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">
##INFO=<ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">
##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">
##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">
##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">
##FILTER=<ID=indelError,Description="Likely artifact due to indel reads at this position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	14907	.	A	G	.	PASS	ADP=59;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:22:60:59:10:49:83.05%:1.3961E-23:29:29:6:4:32:17
chr1	14930	.	A	G	.	PASS	ADP=60;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:21:60:60:12:48:80%:1.5902E-22:29:29:10:2:32:16
chr1	15118	.	A	G	.	PASS	ADP=15;WT=0;HET=1;HOM=0;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	0/1:2:16:15:9:6:40%:8.4291E-3:26:29:0:9:0:6
chr1	15211	.	T	G	.	PASS	ADP=13;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:7:13:13:0:13:100%:9.6148E-8:0:29:0:0:1:12
chr1	16378	.	T	C	.	PASS	ADP=11;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:5:11:11:0:11:100%:1.4176E-6:0:28:0:0:9:2
chr1	16534	.	C	T	.	PASS	ADP=26;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:11:26:26:3:23:88.46%:7.3681E-12:29:29:3:0:14:9
chr1	16571	.	G	A	.	PASS	ADP=21;WT=0;HET=1;HOM=0;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	0/1:5:21:21:8:13:61.9%:7.9741E-6:28:28:5:3:4:9
chr1	714427	.	G	A	.	PASS	ADP=8;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:4:8:8:0:8:100%:7.77E-5:0:28:0:0:7:1
chr1	718386	.	A	G	.	str10	ADP=25;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:14:25:25:0:25:100%:7.9107E-15:0:28:0:0:23:2
chr1	718555	.	T	C	.	PASS	ADP=113;WT=0;HET=0;HOM=1;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:64:113:113:1:112:99.12%:1.994E-65:30:29:1:0:90:22
...
*/
public class VCFRecord {
	
	//fields
	private int position;
	private char reference;
	private char alternate;
	private String filter;
	private String info;
	private String[] sample;
	
	public static final int chromosomeIndex= 0;
	public static final int positionIndex = 1;
	public static final int referenceIndex = 3;
	public static final int alternateIndex = 4;
	public static final int filterIndex = 6;
	public static final int infoIndex = 7;
	public static final int sampleIndex = 9;
	public static final int numberColumnsInVCFRecord = 10;
	public static final Pattern colon = Pattern.compile(":");
	
	/**Only extracts some of the fields from a record*/
	public VCFRecord(String[] fields) throws Exception{
		position = Integer.parseInt(fields[positionIndex]);
		reference = fields[referenceIndex].charAt(0);
		alternate = fields[alternateIndex].charAt(0);
		filter = fields[filterIndex];
		info = fields[infoIndex];
		sample = colon.split(fields[sampleIndex]);
	}
	
	public String getSampleField(int index){
		return sample[index];
	}

	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public char getReference() {
		return reference;
	}

	public void setReference(char reference) {
		this.reference = reference;
	}

	public char getAlternate() {
		return alternate;
	}

	public void setAlternate(char alternate) {
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

	public String[] getSample() {
		return sample;
	}

	public void setSample(String[] sample) {
		this.sample = sample;
	}
}
