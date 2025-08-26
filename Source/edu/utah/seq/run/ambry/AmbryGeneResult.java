package edu.utah.seq.run.ambry;

public class AmbryGeneResult {

	//fields, using Ambry header names
	private AmbryResult ambryResult;
	private String geneResult;			// BARD1
	private String nucleotide_id;		// NM_000465.2
	private String zygostiy;			// Heterozygous or Mosaic
	private String c_variant;			// c.1935_1954dupTGAACAGGAAGAAAAGTATG (bad 5'UTR_EX8dup, 5'UTR_EX8dup, EX9_11del, EX3del, 5'UTR_EX1del, EX12del, 5'UTR_3'UTRdel, EX2_3del, 5'UTR_EX1dup)
	private String p_variant;			// p.E652Vfs*69
	private String method;				// Sequencing & Del/Dup
	private String categoryDescription;// Pathogenic Mutation, Variant of Unknown Significance, Variant Likely Pathogenic, Variant Likely Benign, Benign
	private String chromosome;			// 16
	private String genomicStart;		// 68855966
	private String ref;					// G
	private String alt;					// T
	
	//for the vcf header
	static String infoGene = "##INFO=<ID=aGene,Number=1,Type=String,Description=\"Ambry gene name\">\n";
	static String infoTranscript = "##INFO=<ID=aTrans,Number=1,Type=String,Description=\"Ambry gene transcript\">\n";
	static String infoZygosity = "##INFO=<ID=aZyg,Number=1,Type=String,Description=\"Ambry variant zygosity\">\n";
	static String infoCDot = "##INFO=<ID=aCDot,Number=1,Type=String,Description=\"Ambry cDot\">\n";
	static String infoPDot = "##INFO=<ID=aPDot,Number=1,Type=String,Description=\"Ambry pDot\">\n";
	static String infoInterpretation = "##INFO=<ID=aInterp,Number=1,Type=String,Description=\"Ambry variant effect\">\n";
	static String infoMethod = "##INFO=<ID=aMeth,Number=1,Type=String,Description=\"Ambry variant detection method\">\n";
	
	//cDot info
	private boolean cDotPass = false;
	private String transcriptCDot = null; 	// NM_000465.2:c.1935_1954dupTGAACAGGAAGAAAAGTATG
	private String[] cDotHg38Vcf = null;	// from Jannovar chr16 68855966 G T
	
	//vcf info
	private String vcfHg19Key = null;
	private String[] vcfHg38ChromPosRefAlt = null;
	
	//already has all chrom info
	private boolean chromPass = false;
	
	public AmbryGeneResult(AmbryResult ambryResult, String geneResult, String nucleotide_id, String zygostiy, String c_variant, String p_variant, 
			String method, String categoryDescription, String chromosome, String genomicStart, String ref, String alt) {
		this.ambryResult = ambryResult;
		this.geneResult = geneResult;
		this.nucleotide_id = nucleotide_id;
		this.zygostiy = zygostiy;
		this.c_variant = c_variant;
		this.p_variant = p_variant;
		this.method = method;
		this.categoryDescription = categoryDescription;
		this.chromosome = chromosome;
		// these are hg19
		this.genomicStart = genomicStart;
		this.ref = ref;
		this.alt = alt;
		
		checkResults();
	}
	
	public static String fetchVcfHeader() {
		StringBuilder sb = new StringBuilder();
		sb.append("##fileformat=VCFv4.1\n");
		sb.append("##genomeBuild=Hg38\n");
		sb.append(infoGene);
		sb.append(infoTranscript);
		sb.append(infoZygosity);
		sb.append(infoCDot);
		sb.append(infoPDot);
		sb.append(infoInterpretation);
		sb.append(infoMethod);
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		return sb.toString();
	}
	
	private void checkResults() {
		// do they provide chrom pos ref alt? 
		if (chromosome !=null && genomicStart != null && ref != null && alt != null) chromPass = true;
		// cDot info OK?
		if (c_variant==null || c_variant.length()==0 || c_variant.startsWith("c.")==false || c_variant.contains("SVA") || nucleotide_id==null || nucleotide_id.length()==0) cDotPass = false;
		else {
			cDotPass = true;
			transcriptCDot = nucleotide_id+ ":"+ c_variant;
		}
	}

	public String getTranscriptCDot() {
		return transcriptCDot;
	}

	public String[] getCDotHg38Vcf() {
		return cDotHg38Vcf;
	}

	public void setCDotHg38Vcf(String[] cDotVcf) {
		this.cDotHg38Vcf = cDotVcf;
	}

	public String getGeneResult() {
		return geneResult;
	}

	public String getNucleotide_id() {
		return nucleotide_id;
	}

	public String getZygostiy() {
		return zygostiy;
	}

	public String getC_variant() {
		return c_variant;
	}

	public String getP_variant() {
		return p_variant;
	}

	public String getMethod() {
		return method;
	}

	public String getCategoryDescription() {
		return categoryDescription;
	}

	public String getChromosome() {
		return chromosome;
	}

	public String getGenomicStart() {
		return genomicStart;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt() {
		return alt;
	}

	public boolean iscDotPass() {
		return cDotPass;
	}

	public boolean isChromPass() {
		return chromPass;
	}

	public String getVcfHg19Key() {
		if (vcfHg19Key !=null) return vcfHg19Key;
		if (chromosome==null || genomicStart==null ||  ref==null || alt==null) return null;
		vcfHg19Key = chromosome+"_"+ genomicStart+"_"+ ref+"_"+ alt;
		return vcfHg19Key;
	}

	public String[] getVcfHg38ChromPosRefAlt() {
		return vcfHg38ChromPosRefAlt;
	}

	public void setVcfHg38ChromPosRefAlt(String[] vcfHg38ChromPosRefAlt) {
		this.vcfHg38ChromPosRefAlt = vcfHg38ChromPosRefAlt;
	}

	public String[] getHg38Coordinates() {
		// trust the crossmapped Ambry coordinates first then cDot if necessary
		if (vcfHg38ChromPosRefAlt != null) return vcfHg38ChromPosRefAlt;
		if (cDotHg38Vcf != null) return cDotHg38Vcf;
		return null;
	}
	
	public String getVcfLine(String id) {
		//chrom pos ref alt
		String[] coor = getHg38Coordinates();
		
		// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
		StringBuilder sb = new StringBuilder();
		
		// chrom
		sb.append(coor[0]); sb.append("\t");
		// pos
		sb.append(coor[1]); sb.append("\t");
		// id
		sb.append(id); sb.append("\t");
		// ref
		sb.append(coor[2]); sb.append("\t");
		// alt
		sb.append(coor[3]);
		//qual filter
		sb.append("\t.\t.\t");
		//info
		sb.append("aGene=");
		sb.append(geneResult);
		sb.append(";");
		if (nucleotide_id !=null) {
			sb.append("aTrans=");
			sb.append(nucleotide_id);
			sb.append(";");
		}		
		if (zygostiy !=null) {
			sb.append("aZyg=");
			sb.append(zygostiy);
			sb.append(";");
		}
		if (c_variant !=null) {
			sb.append("aCDot=");
			sb.append(c_variant.replaceAll("'", "_"));
			sb.append(";");
		}
		if (p_variant !=null) {
			sb.append("aPDot=");
			sb.append(p_variant);
			sb.append(";");
		}
		if (method !=null) {
			sb.append("aMeth=");
			sb.append(method.replaceAll(" ", ""));
			sb.append(";");
		}
		if (categoryDescription !=null) {
			sb.append("aInterp=");
			sb.append(categoryDescription.replaceAll(" ", ""));
			sb.append(";");
		}

		return sb.toString();
	}

	public AmbryResult getAmbryResult() {
		return ambryResult;
	}


}
