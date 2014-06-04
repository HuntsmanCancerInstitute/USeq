package edu.utah.tomato.util;

import java.io.File;
import java.util.HashSet;

public class TFConstants {
	public static String ANALYSIS_EXOME_NOVO_RAW = "exome_novo_raw";
	public static String ANALYSIS_EXOME_NOVO_VQSR = "exome_novo_vqsr";
	public static String ANALYSIS_EXOME_BWA_RAW = "exome_bwa_raw";
	public static String ANALYSIS_EXOME_BWA_VQSR = "exome_bwa_vqsr";
	public static String ANALYSIS_EXOME_ALIGN_NOVO = "exome_align_novo";
	public static String ANALYSIS_EXOME_ALIGN_BWA = "exome_align_bwa";
	public static String ANALYSIS_EXOME_ALIGN_BEST = "exome_align_best";
	public static String ANALYSIS_EXOME_METRICS = "exome_metrics";
	public static String ANALYSIS_EXOME_VARIANT_RAW = "exome_variant_raw";
	public static String ANALYSIS_EXOME_VARIANT_VQSR = "exome_variant_vqsr";
	public static String ANALYSIS_EXOME_VARIANT_BEST = "exome_variant_best";
	public static String ANALYSIS_EXOME_BEST = "exome_best";
	public static String ANALYSIS_CUSTOM = "custom";
	
	public static String FILE_FASTQ1 = "fastq1";
	public static String FILE_FASTQ2 = "fastq2";
	
	public static String FILE_BAM = "bam";
	public static String FILE_BAI = "bai";
	public static String FILE_LANE_BAM = "lane_bam";
	public static String FILE_LANE_BAI = "lane_bai";
	public static String FILE_SAMPLE_BAM = "sample_bam";
	public static String FILE_SAMPLE_BAI = "sample_bai";
	public static String FILE_SPLIT_BAM = "split_bam";
	public static String FILE_SPLIT_BAI = "split_bai";
	public static String FILE_FINAL_BAM = "final_bam";
	public static String FILE_FINAL_BAI = "final_bai";
	public static String FILE_REDUCE_BAM = "reduce_bam";
	public static String FILE_REDUCE_BAI = "reduce_bai";
	
	public static String FILE_ID_FINAL_BAM = "idFinalBam";
	public static String FILE_ID_FINAL_BAI = "idFinalBai";
	public static String FILE_ID_REDUCE_BAM = "idReduceBam";
	public static String FILE_ID_REDUCE_BAI = "idReduceBai";
	
	public static String FILE_VCF_RAW = "vcf_raw";
	public static String FILE_VCF_RAW_IDX = "vcf_raw_idx";
	public static String FILE_VCF_FILTER = "vcf_filter";
	public static String FILE_VCF_FILTER_IDX = "vcf_filter_idx";
	public static String FILE_VCF_PASSING = "vcf_passing";
	public static String FILE_VCF_PASSING_IDX = "vcf_passing_idx";
	
	public static String FILE_SAM = "sam";
	public static String FILE_METRICS = "pdf";
	
	public static String TARGET_AGILENTV5UTR = "AgilentAllExonV5UTR";
	public static String TARGET_AGILENTV5 = "AgilentAllExonV5";
	public static String TARGET_AGILENTV4 = "AgilentAllExonV4";
	public static String TARGET_AGILENT50MB = "AgilentAllExon50MB";
	public static String TARGET_NIMBLEGENV2 = "NimbleGenEZCapV2";
	public static String TARGET_NIMBLEGENV3 = "NimbleGenEZCapV3";
	public static String TARGET_TRUSEQ = "TruSeq";
	public static String TARGET_CUSTOM = "custom";
	public static String TARGET_ALL = "all";
	
	public static String PREFIX_SAMPLEID = "sampleid";
	public static String PREFIX_SAMPLENAME = "samplename";
	public static String PREFIX_PUID = "puid";
	
	public static HashSet<String> prefixTypes;
	static
	{
		prefixTypes = new HashSet<String>();
		prefixTypes.add(TFConstants.PREFIX_SAMPLENAME);
		prefixTypes.add(TFConstants.PREFIX_SAMPLEID);
		prefixTypes.add(TFConstants.PREFIX_PUID);
	}

	
	public static HashSet<String> validTypes = new HashSet<String>() {
		private static final long serialVersionUID = 1L;

	{
		add(TFConstants.ANALYSIS_CUSTOM);
		add(TFConstants.ANALYSIS_EXOME_BWA_RAW);
		add(TFConstants.ANALYSIS_EXOME_BWA_VQSR);
		add(TFConstants.ANALYSIS_EXOME_NOVO_RAW);
		add(TFConstants.ANALYSIS_EXOME_NOVO_VQSR);
		add(TFConstants.ANALYSIS_EXOME_ALIGN_BWA);
		add(TFConstants.ANALYSIS_EXOME_ALIGN_NOVO);
		add(TFConstants.ANALYSIS_EXOME_ALIGN_BEST);
		add(TFConstants.ANALYSIS_EXOME_METRICS);
		add(TFConstants.ANALYSIS_EXOME_VARIANT_RAW);
		add(TFConstants.ANALYSIS_EXOME_VARIANT_VQSR);
		add(TFConstants.ANALYSIS_EXOME_VARIANT_BEST);
		add(TFConstants.ANALYSIS_EXOME_BEST);
	}};
	
	public static HashSet<String> validTargets = new HashSet<String>() {
		private static final long serialVersionUID = 1L;

	{
		add(TFConstants.TARGET_AGILENTV5UTR);
		add(TFConstants.TARGET_AGILENTV5);
		add(TFConstants.TARGET_AGILENTV4);
		add(TFConstants.TARGET_AGILENT50MB);
		add(TFConstants.TARGET_NIMBLEGENV2);
		add(TFConstants.TARGET_NIMBLEGENV3);
		add(TFConstants.TARGET_TRUSEQ);
		add(TFConstants.TARGET_CUSTOM);
		add(TFConstants.TARGET_ALL);
	}};
	
}
