package edu.utah.tomato;

import java.io.File;
import java.util.HashSet;

public class TFConstants {
	public static String ANALYSIS_EXOME_NOVO_RAW = "exome_novo_raw";
	public static String ANALYSIS_EXOME_NOVO_VQSR = "exome_novo_vqsr";
	public static String ANALYSIS_EXOME_BWA_RAW = "exome_bwa_raw";
	public static String ANALYSIS_EXOME_BWA_VQSR = "exome_bwa_vqsr";
	public static String ANALYSIS_EXOME_ALIGN_NOVO = "exome_align_novoalign";
	public static String ANALYSIS_EXOME_ALIGN_BWA = "exome_align_bwa";
	public static String ANALYSIS_EXOME_METRICS = "exome_metrics";
	public static String ANALYSIS_EXOME_VARIANT_RAW = "exome_variant_raw";
	public static String ANALYSIS_EXOME_VARIANT_VQSR = "exome_variant_vqsr";
	public static String ANALYSIS_CUSTOM = "custom";
	
	public static String FILE_FASTQ1 = "fastq1";
	public static String FILE_FASTQ2 = "fastq2";
	public static String FILE_BAM = "bam";
	public static String FILE_BAI = "bai";
	public static String FILE_LANE_BAM = "lane_bam";
	public static String FILE_LANE_BAI = "lane_bai";
	public static String FILE_SAMPLE_BAM = "sample_bam";
	public static String FILE_SAMPLE_BAI = "sample_bai";
	public static String FILE_SPLIT_LANE_BAM = "split_lane_bam";
	public static String FILE_SPLIT_LANE_BAI = "split_lane_bai";
	public static String FILE_REALIGN_SAMPLE_BAM = "realign_sample_bam";
	public static String FILE_REALIGN_SAMPLE_BAI = "realign_sample_bai";
	public static String FILE_REDUCE_BAM = "reduce_bam";
	public static String FILE_REDUCE_BAI = "reduce_bai";
	
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
		add(TFConstants.ANALYSIS_EXOME_METRICS);
		add(TFConstants.ANALYSIS_EXOME_VARIANT_RAW);
		add(TFConstants.ANALYSIS_EXOME_VARIANT_VQSR);
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
