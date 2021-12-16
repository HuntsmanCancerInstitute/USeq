package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.RegionScoreText;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/*
##fileformat=VCFv4.1
##ApplyRecalibration="analysis_type=ApplyRecalibration input_file=[] read_buffer_size=null phone_home=STANDARD gatk_key=null tag=NA read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=hg19.fasta nonDeterministicRandomSeed=false disableRandomization=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=1000 baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false performanceLog=null useOriginalQualities=false BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false defaultBaseQualities=-1 validation_strictness=SILENT remove_program_records=false keep_program_records=false unsafe=null num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false version=false input=[(RodBinding name=input source=rawSNP.vcf)] recal_file=(RodBinding name=recal_file source=snp.recal) tranches_file=snp.tranches out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub ts_filter_level=99.0 ignore_filter=null mode=SNP filter_mismatching_base_and_quals=false"
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -7.0875 <= x < 1.6023">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model at VQS Lod < -12535.2851">
##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -12535.2851 <= x < -7.0875">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model">
##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">
##SelectVariants="analysis_type=SelectVariants input_file=[] read_buffer_size=null phone_home=STANDARD gatk_key=null tag=NA read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=hg19.fasta nonDeterministicRandomSeed=false disableRandomization=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=1000 baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false performanceLog=null useOriginalQualities=false BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false defaultBaseQualities=-1 validation_strictness=SILENT remove_program_records=false keep_program_records=false unsafe=null num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false version=false variant=(RodBinding name=variant source=snp.vcf) discordance=(RodBinding name= source=UNBOUND) concordance=(RodBinding name= source=UNBOUND) out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sample_name=[] sample_expressions=null sample_file=null exclude_sample_name=[] exclude_sample_file=[] select_expressions=[VQSLOD>4.0] excludeNonVariants=false excludeFiltered=false restrictAllelesTo=ALL keepOriginalAC=false mendelianViolation=false mendelianViolationQualThreshold=0.0 select_random_fraction=0.0 remove_fraction_genotypes=0.0 selectTypeToInclude=[] keepIDs=null fullyDecode=false forceGenotypesDecode=false justRead=false maxIndelSize=2147483647 ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES=false filter_mismatching_base_and_quals=false"
##UnifiedGenotyper="analysis_type=UnifiedGenotyper input_file=[carcinoid5.bam, carcinoid6.bam, carcinoid7.bam, colon1.bam, colon2.bam, colon3.bam, colon4.bam] read_buffer_size=null phone_home=STANDARD gatk_key=null tag=NA read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=hg19.fasta nonDeterministicRandomSeed=false disableRandomization=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=250 baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false performanceLog=null useOriginalQualities=false BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false defaultBaseQualities=-1 validation_strictness=SILENT remove_program_records=false keep_program_records=false unsafe=null num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false version=false genotype_likelihoods_model=SNP pcr_error_rate=1.0E-4 computeSLOD=false annotateNDA=false pair_hmm_implementation=ORIGINAL min_base_quality_score=17 max_deletion_fraction=0.05 min_indel_count_for_genotyping=5 min_indel_fraction_per_sample=0.25 indel_heterozygosity=1.25E-4 indelGapContinuationPenalty=10 indelGapOpenPenalty=45 indelHaplotypeSize=80 indelDebug=false ignoreSNPAlleles=false allReadsSP=false ignoreLaneInfo=false reference_sample_calls=(RodBinding name= source=UNBOUND) reference_sample_name=null sample_ploidy=2 min_quality_score=1 max_quality_score=40 site_quality_prior=20 min_power_threshold_for_calling=0.95 min_reference_depth=100 exclude_filtered_reference_sites=false heterozygosity=0.001 genotyping_mode=DISCOVERY output_mode=EMIT_VARIANTS_ONLY standard_min_confidence_threshold_for_calling=30.0 standard_min_confidence_threshold_for_emitting=30.0 alleles=(RodBinding name= source=UNBOUND) max_alternate_alleles=6 contamination_fraction_to_filter=0.05 contamination_fraction_per_sample_file=null p_nonref_model=EXACT_INDEPENDENT logRemovedReadsFromContaminationFiltering=null exactcallslog=null dbsnp=(RodBinding name=dbsnp source=hg19.dbsnp.vcf) comp=[] out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub debug_file=null metrics_file=null annotation=[] excludeAnnotation=[] filter_mismatching_base_and_quals=false"
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##contig=<ID=chrM,length=16571,assembly=hg19>
##reference=file:///home/u0028003/PIs/Neklason/Alignments/hg19.fasta
##source=SelectVariants
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	carcinoid5	carcinoid6	carcinoid7	colon1	colon2	colon3	colon4
chr1	54708	.	G	C	51.35	PASS	AC=3;AF=0.750;AN=4;BaseQRankSum=0.736;DP=3;Dels=0.00;FS=4.771;HaplotypeScore=0.4999;MLEAC=3;MLEAF=0.750;MQ=70.00;MQ0=0;MQRankSum=0.736;QD=17.12;ReadPosRankSum=0.736;VQSLOD=10.70;culprit=MQ	GT:AD:DP:GQ:PL	./.	0/1:1,1:2:34:36,0,34	./.	1/1:0,1:1:3:41,3,0	./.	./.	./.
chr1	714427	rs12028261	G	A	57.36	PASS	AC=4;AF=1.00;AN=4;DB;DP=2;Dels=0.00;FS=0.000;HaplotypeScore=0.0000;MLEAC=4;MLEAF=1.00;MQ=70.00;MQ0=0;QD=28.68;VQSLOD=22.08;culprit=HaplotypeScore	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:39,3,0	./.	./.	./.	./.	1/1:0,1:1:3:41,3,0
chr1	725822	rs199845677	G	A	45.59	PASS	AC=1;AF=0.167;AN=6;BaseQRankSum=-1.231;DB;DP=5;Dels=0.00;FS=3.979;HaplotypeScore=0.3333;MLEAC=1;MLEAF=0.167;MQ=70.00;MQ0=0;MQRankSum=0.358;QD=15.20;ReadPosRankSum=0.358;VQSLOD=11.31;culprit=HaplotypeScore	GT:AD:DP:GQ:PL	0/0:1,0:1:3:0,3,42	0/1:1,2:3:32:72,0,32	./.	./.	./.	0/0:1,0:1:3:0,3,43	./.

 */
public class VCFParser {

	//user defined fields
	private File vcfFile;

	//internal fields
	private VCFRecord[] vcfRecords;
	private VCFComments comments;
	public static final Pattern TAB = Pattern.compile("\\t");
	public static final Pattern COLON = Pattern.compile(":");
	public static final Pattern COMMA = Pattern.compile(",");
	public static final Pattern SLASH = Pattern.compile("/");
	public static final Pattern ID = Pattern.compile(".+=<ID=(.+)");
	private ArrayList<String> badVcfRecords = new ArrayList<String>();
	private HashMap<String, VCFLookUp> chromosomeVCFRecords = null;
	private boolean loadRecords = true;
	private boolean loadSamples = true;
	private boolean loadInfo = true;
	private boolean missingQual = false;

	//indexs for ripping vcf records
	int chromosomeIndex= 0;
	int positionIndex = 1;
	int rsIndex = 2;
	int referenceIndex = 3;
	int alternateIndex = 4;
	int qualityIndex = 5;
	int filterIndex = 6;
	int infoIndex = 7;
	int formatIndex = 8;
	int firstSampleIndex = 9;
	int minimumNumberFields = 8;
	int numberFields = 0;



	private String[] sampleNames = null;

	//Chunking arguments
	int chunkNumber = 0;
	int chunkSize = Integer.MAX_VALUE; //hopefully a vcf file will never have more than 2.1 Billion lines



	//Constructors
	/**If loadRecords is false then just the header and #CHROM line with sample names is parsed.
	 * If loadSamples is false then none of the sample info is loaded, just the required vcf fields. */
	public VCFParser(File vcfFile, boolean loadRecords, boolean loadSamples, boolean loadInfo) {
		this. vcfFile = vcfFile;
		this.loadRecords = loadRecords;
		this.loadSamples = loadSamples;
		this.loadInfo = loadInfo;
		parseVCF();
	}

	public VCFParser(File vcfFile, boolean loadRecords, boolean loadSamples, boolean loadInfo, int chunkNumber, int chunkSize) {
		this. vcfFile = vcfFile;
		this.loadRecords = loadRecords;
		this.loadSamples = loadSamples;
		this.loadInfo = loadInfo;
		this.chunkNumber = chunkNumber;
		this.chunkSize = chunkSize;
		parseVCF();
	}
	
	/**If loadRecords is false then just the header and #CHROM line with sample names is parsed.
	 * If loadSamples is false then none of the sample info is loaded, just the required vcf fields. */
	public VCFParser(File vcfFile) {
		this. vcfFile = vcfFile;
		this.loadRecords = true;
		this.loadSamples = true;
		this.loadInfo = true;
		missingQual = true;
		parseVCF();
	}

	public VCFParser() {}

	//Methods

	/**Adds a chr onto chromosome names that lack it. Converts chrMT to chrM.*/
	public void appendChrFixMT(){
		for (VCFRecord r : vcfRecords){
			r.appendChr();
			r.correctChrMTs();
		}
	}
	
	/**Call after initializeParser(). Returns a record or null.  If null also closes the BufferedReader.
	 * No checking of position or chrom!*/
	public VCFRecord fetchNext(BufferedReader in) throws Exception{
		String line = in.readLine(); 
		if (line == null) {
			in.close();
			return null;
		}
		VCFRecord vcf = new VCFRecord(line, this, loadSamples, loadInfo);
		if (!missingQual && vcf.isMissingQual()) missingQual = true;
		return vcf;
	}
	
	/**Call first if you then want to load a vcf file record by record. Don't forget to close the BufferedReader when done!*/
	public BufferedReader initializeParser() {
		BufferedReader in = null;
		String line = null;
		try {
			in  = IO.fetchBufferedReader(vcfFile);
			ArrayList<String> commentsAL = new ArrayList<String>();

			//find "#CHROM" line and parse sample names
			boolean notFound = true;
			while ((line=in.readLine()) != null){
				//comments
				if (line.startsWith("#")){
					commentsAL.add(line);
					if (line.startsWith("#CHROM")) {
						notFound = false;
						String[] header = TAB.split(line);
						int numSampleNames = header.length - firstSampleIndex;
						if (numSampleNames < 1) sampleNames = null;
						else {
							sampleNames = new String[numSampleNames];
							int index =0;
							for (int i=firstSampleIndex; i< header.length; i++) sampleNames[index++] = header[i];
						}
						break;
					}
				}
			}
			if (notFound) throw new Exception("\nFailed to find the #CHROM header line.");
			comments = new VCFComments(commentsAL);
			return in;
		}catch (Exception e) {
			System.err.println("\nAborting, problem parsing vcf file -> "+vcfFile);
			e.printStackTrace();
			System.exit(1);
		} 
		return null;
	}


	/**Parses a sorted vcf file. Indicate whether you want to load the data and not just parse the header and #CHROM line.*/
	public void parseVCF() {
		BufferedReader in = null;
		String line = null;
		int badCounter = 0;
		try {
			in  = IO.fetchBufferedReader(vcfFile);
			ArrayList<VCFRecord> records = new ArrayList<VCFRecord>(10000);
			ArrayList<String> commentsAL = new ArrayList<String>();

			//find "#CHROM" line and parse sample names
			boolean notFound = true;
			while ((line=in.readLine()) != null){
				//comments
				if (line.startsWith("#")){
					commentsAL.add(line);
					if (line.startsWith("#CHROM")) {
						notFound = false;
						String[] header = TAB.split(line);
						int numSampleNames = header.length - firstSampleIndex;
						if (numSampleNames < 1) sampleNames = null;
						else {
							sampleNames = new String[numSampleNames];
							int index =0;
							for (int i=firstSampleIndex; i< header.length; i++) sampleNames[index++] = header[i];
						}
						break;
					}
				}
			}
			if (notFound) throw new Exception("\nFailed to find the #CHROM header line.");
			comments = new VCFComments(commentsAL);

			//load data?
			if (loadRecords == false) return;

			//Initilize variables
			String oldChrom = "";
			int oldPosition = -1;
			HashSet<String> chromNames = new HashSet<String>();

			//Setup chunking variables
			int minRecord = this.chunkSize * this.chunkNumber;
			int maxRecord = this.chunkSize * this.chunkNumber + this.chunkSize;
			int recordCount = 0;			
			
			while ((line=in.readLine()) != null){
				if (recordCount < minRecord) { //Ignore records less than starting point
					//pass
				} else if (recordCount >= maxRecord) { //Quit looking for records beyond ending point
					break;
				} else { //Juuust right.
					try {					
						VCFRecord vcf = new VCFRecord(line, this, loadSamples, loadInfo);
						
						if (!missingQual && vcf.isMissingQual()) {
							//this is a bothersome warning that disrupts alot of apps that work with NIST and VarScan vcf files so I'm Nixing it. 
							//System.out.println("\n\nWARNING: Sample variant quality (GQ) is missing from record: " + line + ".\n\nMissing sample qualities are set to zero.  Missing quality information is normal "
									//+ "when parsing VCF files from tumor/normal varscan, but not VCF files from GATK.  This error is only reported once to reduce "
									//+ "on-screen messages, so there may be other variants throughout the file with missing qualities.\n\n");
							missingQual = true;
						}
						//old chrom
						if (vcf.getChromosome().equals(oldChrom)){
							//check position
							if (vcf.getPosition() < oldPosition) throw new Exception("New vcf record position is < prior position!  Is this file sorted?");
						}
						//nope new
						else {
							//chrom seen before?
							if (chromNames.contains(vcf.getChromosome())) throw new Exception("New vcf record chromosome has been seen before!  Is this file sorted?");
							oldChrom = vcf.getChromosome();
							chromNames.add(oldChrom);
						}
						//add record and set position
						records.add(vcf);	
						oldPosition = vcf.getPosition();

					} catch (Exception e) {
						System.err.println("Skipping malformed VCF Record-> "+line);
						System.err.println("Error-> "+e.getMessage());
						if (badCounter++ > 10) {
							throw new Exception("\nToo many malformed VCF Records.\n");
						}
						badVcfRecords.add(line);
					}					
				}
				recordCount++;
			}

			//save array
			vcfRecords = new VCFRecord[records.size()];
			records.toArray(vcfRecords);

		}catch (Exception e) {
			System.err.println("\nAborting, problem parsing vcf file -> "+vcfFile);
			e.printStackTrace();
		} finally{
			try {
				in.close();
			} catch (IOException e) {}
		}
	}
	
	/**Merges headers eliminating duplicate lines.  Does a bad ID name collision checking, silently keeps first one. Source's are concatenated with a _.
	 * Returns null if CHROM lines differ. */
	public static String[] mergeHeaders(VCFParser[] parsers, boolean tossSampleInfo, boolean useFirstCHROM) {
		LinkedHashSet<String> other = new LinkedHashSet<String>();
		LinkedHashSet<String> contig = new LinkedHashSet<String>();
		LinkedHashSet<String> info = new LinkedHashSet<String>();
		LinkedHashSet<String> filter = new LinkedHashSet<String>();
		LinkedHashSet<String> format = new LinkedHashSet<String>();
		TreeSet<String> source = new TreeSet<String>();
		String chromLine = null;
		
		//for each header line, load hash
		for (VCFParser p : parsers){
			String[] header = p.getStringComments();
			for (String h: header){
				h=h.trim();
				if (h.startsWith("##contig")){
					if (contig.contains(h) == false) contig.add(h);
				}
				else if (h.startsWith("##INFO")){
					if (info.contains(h) == false) info.add(h);
				}
				else if (h.startsWith("##FILTER")){
					if (filter.contains(h) == false) filter.add(h);
				}
				else if (h.startsWith("##FORMAT")){
					if (format.contains(h) == false) format.add(h);
				}
				else if (h.startsWith("##source=")){
					source.add(parseSourceName(h));
				}
				else if (h.startsWith("#CHROM") && tossSampleInfo == false) {
					if (chromLine == null) chromLine = h;
					else if (useFirstCHROM==false && chromLine.equals(h) == false) return null;
				}
				else if (other.contains(h) == false) {
					other.add(h);
				}
			}
		}
		//remove ID dups from contig, filter, format, info
		ArrayList<String> contigAL = mergeHeaderIds(contig);
		ArrayList<String> filterAL = mergeHeaderIds(filter);
		ArrayList<String> formatAL = mergeHeaderIds(format);
		ArrayList<String> infoAL = mergeHeaderIds(info);
		
		if (tossSampleInfo) {
			formatAL.clear();
			chromLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
		}
		
		//Concatenate sources
		String concatSource = null;
		int numS = source.size();	
		if (numS>0){
			StringBuilder sb = new StringBuilder("##source=");
			Iterator<String> it = source.iterator();
			while (it.hasNext()){
				sb.append(it.next());
				if (numS>1) sb.append("_");
			}
			if (numS>1) sb.append("USeqMerged");
			concatSource = sb.toString();
		}
		
		int num = other.size();
		num+= contigAL.size();
		num+= infoAL.size();
		num+= filterAL.size();
		num+= formatAL.size();
		if (chromLine != null) num++;
		if (concatSource != null) num++;
		
		String[] merge = new String[num];
		int counter = 0;
		for (String s : other) merge[counter++] = s;
		if (concatSource != null) merge[counter++] = concatSource;
		for (String s : contigAL) merge[counter++] = s;
		for (String s : filterAL) merge[counter++] = s;
		for (String s : infoAL) merge[counter++] = s;
		for (String s : formatAL) merge[counter++] = s;
		if (chromLine != null) merge[num-1] = chromLine;
		
		//Misc.printArray(merge);
		
		return merge;
	}
	
	/**Parses the name of the app used to create the vcf from it's ##source= header, tested with strelka and lofreq.*/
	public static String parseSourceName(String h){
		//examples
		//##source=strelka
		//##source=/uufs/chpc.utah.edu/common/home/arup-storage3/Reference/Apps/LoFreq/lofreq_star-2.1.2/bin/lofreq call -m 20 -q 13 -Q 13 -f /uufs/chpc.utah.edu/common/home/arup-storage3/Reference/Data/B37/BWA/APCompCanPseudoMaskedv3/mod_3_human_g1k_v37_decoy.fasta -o 11054N_Lofreq/lofreq.snv.vcf -s -S /uufs/chpc.utah.edu/common/home/arup-storage3/Reference/Data/B37/DbSNP/dbSNP_146_20151104/All_20151104.vcf.gz -l /uufs/chpc.utah.edu/common/home/arup-storage3/Reference/Bed/B37/Somatic/APCompCancer/V3/A3370/compcancer_321transcripts_118gene_cds_v2_bed3.bed 11054N_final.bam

		//split on = sign
		String[] s = Misc.PATTERN_EQUALS.split(h);
		
		//split on whitespace to pull actual source
		String[] sourceToken = Misc.WHITESPACE.split(s[1]);
		
		//split on /, to pull lofreq name
		String[] name = Misc.FORWARD_SLASH.split(sourceToken[0]);
		return name[name.length-1];
	}
	
	/*Parses ID, keeps first one if dups found.*/
	public static ArrayList<String> mergeHeaderIds(LinkedHashSet<String> comments){
		ArrayList<String> toReturn = new ArrayList<String>();
		HashSet<String> ids = new HashSet<String>();
		//hash on =<ID=xxxx, if collision keeps first one
		for (String s : comments){
			String[] f = COMMA.split(s);
			Matcher mat = ID.matcher(f[0]);
			if (mat.matches()){
				//System.out.println(f[0]+" -> "+mat.group(1)+" ori "+s);
				if (ids.contains(mat.group(1)) == false){
					toReturn.add(s);
					ids.add(mat.group(1));
				}
			}
			else Misc.printErrAndExit("\nError: failed to find a match to '<ID=xxxx,' in this  VCF header line:\n"+s+"\n");
		}
		return toReturn;
	}

	public double calculateTiTvRatio(){
		double numTransitions = 0;
		double numTransversions = 0;
		for (VCFRecord r: vcfRecords) {
			if (r.isSNP() == false) continue;
			if (r.isTransition()) numTransitions++;
			else numTransversions++;
		}
		return numTransitions/numTransversions;
	}

	/**This reloads the chromosomeVCFRecords HashMap<String (chromosome), VCFLookUp>() with the current vcfRecords. */
	public void splitVCFRecordsByChromosome() {
		try {
			chromosomeVCFRecords = new HashMap<String, VCFLookUp>();

			if (vcfRecords == null || vcfRecords.length == 0) return;

			String oldChrom = vcfRecords[0].getChromosome();
			int oldPosition = vcfRecords[0].getPosition();
			ArrayList<VCFRecord> recordsAL = new ArrayList<VCFRecord>();
			recordsAL.add(vcfRecords[0]);
			ArrayList<Integer> positionsAL = new ArrayList<Integer>();
			positionsAL.add(oldPosition);
			//for each record
			for (int i=1; i< vcfRecords.length; i++){
				String testChrom = vcfRecords[i].getChromosome();
				int testPosition = vcfRecords[i].getPosition();
				//is it the same chrom?
				if (oldChrom.equals(testChrom) ){
					//check position
					if (testPosition < oldPosition) throw new Exception("\nError: subsequent vcf record's position is less than prior position!  Is this file sorted? See -> "+vcfRecords[i]);
				}
				//nope so close old and start new
				else {
					//close old
					VCFRecord[] v = new VCFRecord[recordsAL.size()];
					recordsAL.toArray(v);
					recordsAL.clear();
					chromosomeVCFRecords.put(oldChrom, new VCFLookUp(Num.arrayListOfIntegerToInts(positionsAL), v));
					positionsAL.clear();
					//create new
					//has new chrom been seen before?
					if (chromosomeVCFRecords.containsKey(testChrom)) throw new Exception ("\nError: subsequent vcf record's reference chromosome has been seen before! Is this file sorted? See -> "+vcfRecords[i]);
					oldChrom = testChrom;
				}
				//save info
				recordsAL.add(vcfRecords[i]);
				positionsAL.add(testPosition);
				oldPosition = testPosition;
			}
			//set last
			VCFRecord[] v = new VCFRecord[recordsAL.size()];
			recordsAL.toArray(v);
			chromosomeVCFRecords.put(oldChrom, new VCFLookUp(Num.arrayListOfIntegerToInts(positionsAL), v));
		} catch (Exception e){
			System.err.println("\nError spliting records by chromosome!");
			e.printStackTrace();
		}
	}

	public int[] filterVCFRecords(HashMap<String,RegionScoreText[]> goodRegions){
		int starting = countMatchingVCFRecords(VCFRecord.PASS);

		//call chrom splitter
		getChromosomeVCFRecords();
		
		//set scores to -1
		for (VCFRecord r: vcfRecords) r.setScore(-1f);

		//for each interrogated region
		for (String chr: goodRegions.keySet()){		
			VCFLookUp vcf = chromosomeVCFRecords.get(chr);
			if (vcf == null) continue;
			//make boolean array representing covered bases
			RegionScoreText[] r = goodRegions.get(chr);
			if (r == null || r.length == 0) continue;
			int lastBase = RegionScoreText.findLastBase(r);
			boolean[] coveredBases = new boolean[lastBase];
			for (int i=0; i< r.length; i++){
				int stop = r[i].getStop();
				for (int j=r[i].getStart(); j< stop; j++){
					coveredBases[j] = true;
				}
			}
			VCFRecord[] vcfRecords = vcf.getVcfRecord();
			//for each vcf position, is it covered?
			for (int i=0; i< vcfRecords.length; i++){
				if (vcfRecords[i].getPosition() < lastBase && coveredBases[vcfRecords[i].getPosition()] == true) {
					vcfRecords[i].setFilter(VCFRecord.PASS);
				}
				else {
					vcfRecords[i].setFilter(VCFRecord.FAIL);
				}
				vcfRecords[i].setScore(0f);
			}
		}
		//fail any records with a -1 score, these weren't scored thus couldn't have intersected a region
		for (VCFRecord r: vcfRecords) {
			if (r.getScore() == -1f) {
				r.setScore(0f);
				r.setFilter(VCFRecord.FAIL);
			}
		}
		int ending = countMatchingVCFRecords(VCFRecord.PASS);
		
		//remove non passing
		filterVCFRecords(VCFRecord.PASS);
		
		return new int[]{starting,ending};
	}

	/**This is a hard filter that only keeps VCFRecords whose Filter field matches the matchFilterText2Keep.
	 * Will reload the chromosomeVCFRecords if needed.*/
	public void filterVCFRecords(String matchFilterText2Keep){
		ArrayList<VCFRecord> keep = new ArrayList<VCFRecord>();
		for (VCFRecord r : vcfRecords){
			if (r.getFilter().equals(matchFilterText2Keep)) keep.add(r);
		}
		vcfRecords = new VCFRecord[keep.size()];
		keep.toArray(vcfRecords);

		//remake split chrom data?
		if (chromosomeVCFRecords != null) splitVCFRecordsByChromosome();
	}
	
	public ArrayList<VCFRecord> filterVCFRecordsOnFILTER(String filterTxtToExclude) {
		ArrayList<VCFRecord> keep = new ArrayList<VCFRecord>();
		ArrayList<VCFRecord> fail = new ArrayList<VCFRecord>();
		for (VCFRecord r : vcfRecords){
			if (r.getFilter().contains(filterTxtToExclude) == false) keep.add(r);
			else fail.add(r);
		}
		vcfRecords = new VCFRecord[keep.size()];
		keep.toArray(vcfRecords);

		//remake split chrom data?
		if (chromosomeVCFRecords != null) splitVCFRecordsByChromosome();
		return fail;
	}
	
	public ArrayList<VCFRecord> filterVCFRecordsOnQUAL(float minQualFilt) {
		ArrayList<VCFRecord> keep = new ArrayList<VCFRecord>();
		ArrayList<VCFRecord> fail = new ArrayList<VCFRecord>();
		for (VCFRecord r : vcfRecords){
			if (r.getQuality() >= minQualFilt) keep.add(r);
			else fail.add(r);
		}
		vcfRecords = new VCFRecord[keep.size()];
		keep.toArray(vcfRecords);

		//remake split chrom data?
		if (chromosomeVCFRecords != null) splitVCFRecordsByChromosome();
		return fail;
	}

	public int countMatchingVCFRecords(String matchFilterText){
		int numMatches = 0;
		for (VCFRecord r : vcfRecords){
			if (r.getFilter().equals(matchFilterText)) numMatches++;
		}
		return numMatches;
	}

	public int countNonMatchingVCFRecords(String matchFilterText){
		int numNonMatches = 0;
		for (VCFRecord r : vcfRecords){
			if (r.getFilter().equals(matchFilterText) == false) numNonMatches++;
		}
		return numNonMatches;
	}

	/**Sets the filter field in each record to the indicated text.*/
	public void setFilterFieldOnAllRecords (String text){
		for (VCFRecord r : vcfRecords) r.setFilter(text);
	}

	/**Sets the filter field in each record to the indicated text if it is '.'*/
	public void setFilterFieldPeriodToTextOnAllRecords (String text){
		for (VCFRecord r : vcfRecords) {
			if (r.getFilter().equals(".")) r.setFilter(text);
		}
	}

	/**Sets the GQ sample score as the thresholding score in each VCFRecord.
	 * Returns the min and max scores found.*/
	public float[] setGenotypeQualityGQScore(int sampleIndex) {
		float minScore = vcfRecords[0].getSample()[sampleIndex].getGenotypeQualityGQ();
		float maxScore = minScore;
		for (VCFRecord r : vcfRecords) {
			float score = r.getSample()[sampleIndex].getGenotypeQualityGQ();
			if (score < minScore) minScore = score;
			if (score > maxScore) maxScore = score;
			r.setScore((float)score);
		}
		return new float[]{minScore, maxScore};
	}

	/**Sets the QUAL score as the thresholding score in each VCFRecord.
	 * Returns the min and max scores found.*/
	public float[] setRecordQUALAsScore() {
		float minScore = vcfRecords[0].getQuality();
		float maxScore = minScore;		
		for (VCFRecord r : vcfRecords) {
			float score = r.getQuality();
			if (score < minScore) minScore = score;
			if (score > maxScore) maxScore = score;
			r.setScore((float)score);		
		}
		return new float[]{minScore, maxScore};
	}

	/**Sets the infoNamed score as the thresholding score in each VCFRecord.
	 * Returns the min and max scores found.
	 * Throw an exception if it can't be found or parsed as a double.
	 * e.g. VQSLOD or AVR*/
	public float[] setRecordScore(String infoName) throws Exception{
		float minScore = vcfRecords[0].getInfoObject().getInfoFloat(infoName);
		float maxScore = minScore;
		for (VCFRecord r : vcfRecords) {
			float score = r.getInfoObject().getInfoFloat(infoName);
			if (score < minScore) minScore = score;
			if (score > maxScore) maxScore = score;
			r.setScore((float)score);
		}
		return new float[]{minScore, maxScore};
	}

	/**Prints out all unmodified records to a file with the specified suffix
	 * 
	 * @param suffix
	 */
	public void printRecords(File outfile) {
		printRecords(outfile,false);
	}

	/**Prints out all modified records to a file with the given suffix 
	 * 
	 * @param suffix
	 * @param modified
	 */
	public void printRecords(File outfile,boolean modified) {
		printRecords(outfile,modified,this.getVcfComments().getInfoOrder());
	}

	/**Prints out all modified records to a file with the given suffix.  Only specified INFO fields are used.
	 * 
	 * @param suffix
	 * @param modified
	 * @param infoToUse
	 */
	public void printRecords(File outfile,boolean modified, ArrayList<String> infoToUse) {
		printRecords(outfile,modified,infoToUse,VCFInfo.UNMODIFIED);
	}


	/** The method that actually writes the records to file. This version of the method writes to a single file*/
	public void printRecords(File records, boolean modified, ArrayList<String> infoToUse, String style) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(records));

			//Write header
			this.getVcfComments().setAnnotationStyle(style);
			String[] comments;
			if (infoToUse != null) {
				comments = this.getStringComments(infoToUse);
			} else {
				comments = this.getStringComments();
			}
			for (String comment: comments) {
				bw.write(comment + "\n");
			}

			//Write body
			if (modified) {
				for (VCFRecord r: vcfRecords) {
					bw.write(r.getModifiedRecord(infoToUse, style) + "\n");
				}
			} else {
				for (VCFRecord r: vcfRecords) {
					bw.write(r + "\n");
				}
			}
			bw.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void printRecords(VCFRecord[] sortedRecords, File file) {
		try {
			Gzipper out = new Gzipper(file);
			//print header
			for (String c: comments.getOriginalComments()) out.println(c);
			//print records
			for (VCFRecord r: sortedRecords) out.println(r);
			out.close();
		} catch (Exception e) {
			System.err.println("\nProblem printing out vcf records for "+file);
			e.printStackTrace();
			System.exit(1);
		}

	}


	/**The method that actually writes the records to file.  This version writes Passed/Failed records to separate files*/
	public void printFilteredRecords(File records, String filter) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(records));

			for (String comment: this.comments.getOriginalComments()) {
				bw.write(comment + "\n");
			}

			for (VCFRecord r: vcfRecords) {
				if (r.getFilter().equals(filter)) {
					bw.write(r + "\n");
				}
			}
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void removeSNPs() {
		//set all to pass
		setFilterFieldOnAllRecords(VCFRecord.PASS);
		//set snps to fail
		for (VCFRecord r : vcfRecords){
			if (r.isSNP()) r.setFilter(VCFRecord.FAIL);
		}
		//remove fails
		filterVCFRecords(VCFRecord.PASS);
	}

	public void removeNonSNPs() {
		//set all to fail
		setFilterFieldOnAllRecords(VCFRecord.FAIL);
		//set snps to Pass
		for (VCFRecord r : vcfRecords){
			if (r.isSNP()) r.setFilter(VCFRecord.PASS);
		}
		//keep pass
		filterVCFRecords(VCFRecord.PASS);
	}
	
	/**Returns snvs, indels, other.*/
	public VCFRecord[][] splitVCFRecordsBySnvIndelOther() {
		ArrayList<VCFRecord> snvAL = new ArrayList<VCFRecord>();
		ArrayList<VCFRecord> indelAL = new ArrayList<VCFRecord>();
		ArrayList<VCFRecord> otherAL = new ArrayList<VCFRecord>();
		for (VCFRecord rec: vcfRecords){
			if (rec.isSNP()) snvAL.add(rec);
			else if (rec.isDeletion() || rec.isInsertion()) indelAL.add(rec);
			else otherAL.add(rec);
		}
		VCFRecord[] snvs = new VCFRecord[snvAL.size()];
		VCFRecord[] indels = new VCFRecord[indelAL.size()];
		VCFRecord[] other = new VCFRecord[otherAL.size()];
		snvAL.toArray(snvs);
		indelAL.toArray(indels);
		otherAL.toArray(other);
		return new VCFRecord[][]{snvs, indels, other};
	}

	public VCFRecord[] getVcfRecords() {
		return vcfRecords;
	}

	public void setVcfRecords(VCFRecord[] vcfRecords) {
		this.vcfRecords = vcfRecords;
	}

	public String[] getStringComments() {
		String[] c= new String[comments.getOriginalComments().size()];
		comments.getOriginalComments().toArray(c);
		return c;
	}

	public VCFComments getVcfComments() {
		return comments;
	}

	public String[] getStringComments(ArrayList<String> infoFieldsToUse) {
		return comments.getComments(infoFieldsToUse);
	}

	public String[] getSampleNames() {
		return sampleNames;
	}

	public void setSampleNames(String[] sampleNames) {
		this.sampleNames = sampleNames;
	}

	public ArrayList<String> getBadVcfRecords() {
		return badVcfRecords;
	}

	public HashMap<String, VCFLookUp> getChromosomeVCFRecords() {
		if (chromosomeVCFRecords == null) splitVCFRecordsByChromosome();
		return chromosomeVCFRecords;
	}

	public File getVcfFile() {
		return vcfFile;
	}

}
