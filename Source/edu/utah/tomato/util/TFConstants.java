package edu.utah.tomato.util;

import java.io.File;
import java.util.HashSet;

public class TFConstants {
	public static String MENU_UGP_FULL = "ugp_full";
	public static String MENU_UGP_ALIGN = "ugp_align";
	public static String MENU_UGP_VARIANT = "ugp_variant";
	public static String MENU_CORE_FULL = "core_full";
	public static String MENU_CORE_ALIGN = "core_align";
	public static String MENU_CORE_VARIANT = "core_variant";
	public static String MENU_METRICS = "metrics";
	
	
	public static String ANALYSIS_CORE_FULL_NOV = "core_full_nov";
	public static String ANALYSIS_CORE_FULL_BWA = "core_full_bwa";
	public static String ANALYSIS_CORE_ALIGN_NOV = "core_align_nov";
	public static String ANALYSIS_CORE_ALIGN_BWA = "core_align_bwa";
	public static String ANALYSIS_CORE_VARIANT = "core_variant";
	public static String ANALYSIS_METRICS = "metrics";
	public static String ANALYSIS_UGP_FULL_NOV_RAW = "ugp_full_nov_raw";
	public static String ANALYSIS_UGP_FULL_NOV_VQSR = "ugp_full_nov_vqsr";
	public static String ANALYSIS_UGP_FULL_BWA_RAW = "ugp_full_bwa_raw";
	public static String ANALYSIS_UGP_FULL_BWA_VQSR = "ugp_full_bwa_vqsr";
	public static String ANALYSIS_UGP_ALIGN_NOV = "ugp_align_nov";
	public static String ANALYSIS_UGP_ALIGN_BWA = "ugp_align_bwa";
	public static String ANALYSIS_UGP_VARIANT_RAW = "ugp_variant_raw";
	public static String ANALYSIS_UGP_VARIANT_VQSR = "ugp_variant_vqsr";

		
	public static String FILE_FASTQ1 = "fastq1";
	public static String FILE_FASTQ2 = "fastq2";
	
	public static String FILE_BAM = "bam";
	public static String FILE_BAI = "bai";
	public static String FILE_SAMPLE_BAM = "sample_bam";
	public static String FILE_SAMPLE_BAI = "sample_bai";
	public static String FILE_FINAL_BAM = "final_bam";
	public static String FILE_FINAL_BAI = "final_bai";
	
	public static String FILE_GVCF = "gvcf";
	public static String FILE_GVCF_IDX = "gvcf_dx";
	public static String FILE_GVCF_STUDY = "gvcf_study";
	public static String FILE_GVCF_STUDY_IDX = "gvcf_study_idx";

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

	public static HashSet<String> validMenu = new HashSet<String>() {
	{
		add(TFConstants.MENU_METRICS);
		add(TFConstants.MENU_UGP_FULL);
		add(TFConstants.MENU_UGP_ALIGN);
		add(TFConstants.MENU_UGP_VARIANT);
		add(TFConstants.MENU_CORE_FULL);
		add(TFConstants.MENU_CORE_ALIGN);
		add(TFConstants.MENU_CORE_VARIANT);
	}};
	
	public static HashSet<String> alignOptions = new HashSet<String>() {
	{
		add(TFConstants.MENU_CORE_FULL);
		add(TFConstants.MENU_CORE_ALIGN);
		add(TFConstants.MENU_UGP_FULL);
		add(TFConstants.MENU_UGP_ALIGN);
	}};
	
	public static HashSet<String> variantOptions = new HashSet<String>() {
	{
		add(TFConstants.MENU_UGP_FULL);
		add(TFConstants.MENU_UGP_VARIANT);
	}};
	
	public static HashSet<String> alignTypes = new HashSet<String>() {
	{
		add("core_bwa");
		add("core_nov");
		add("ugp_nov");
		add("ugp_bwa");
	}};
	
	public static HashSet<String> metricsTypes = new HashSet<String>() {
	{
		add("metrics");
	}};
	
	public static HashSet<String> variantTypes = new HashSet<String>() {
	{
		add("core_variant");
		add("ugp_raw");
		add("ugp_vqsr");
	}};
	
	
	public static HashSet<String> validTypes = new HashSet<String>() {
	{
		add(TFConstants.ANALYSIS_CORE_FULL_NOV);
		add(TFConstants.ANALYSIS_CORE_FULL_BWA);
		add(TFConstants.ANALYSIS_CORE_ALIGN_NOV);
		add(TFConstants.ANALYSIS_CORE_ALIGN_BWA);
		add(TFConstants.ANALYSIS_CORE_VARIANT);
		add(TFConstants.ANALYSIS_METRICS);
		add(TFConstants.ANALYSIS_UGP_FULL_BWA_VQSR);
		add(TFConstants.ANALYSIS_UGP_FULL_BWA_RAW);
		add(TFConstants.ANALYSIS_UGP_FULL_NOV_VQSR);
		add(TFConstants.ANALYSIS_UGP_FULL_BWA_RAW);
		add(TFConstants.ANALYSIS_UGP_ALIGN_NOV);
		add(TFConstants.ANALYSIS_UGP_ALIGN_BWA);
		add(TFConstants.ANALYSIS_UGP_VARIANT_RAW);
		add(TFConstants.ANALYSIS_UGP_VARIANT_VQSR);
	}};
	
	
	
	public static HashSet<String> validTargets = new HashSet<String>() {
	
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
