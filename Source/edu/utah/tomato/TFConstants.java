package edu.utah.tomato;

import java.io.File;
import java.util.HashSet;

public class TFConstants {
	public static String ANALYSIS_EXOME_NOVOALIGN = "exome_novoalign";
	public static String ANALYSIS_EXOME_BWA = "exome_bwa";
	public static String ANALYSIS_EXOME_ALIGN_NOVOALIGN = "exome_align_novoalign";
	public static String ANALYSIS_EXOME_ALIGN_BWA = "exome_align_bwa";
	public static String ANALYSIS_EXOME_METRICS = "exome_metrics";
	public static String ANALYSIS_EXOME_VARIANT_RAW = "exome_variant_raw";
	public static String ANALYSIS_CUSTOM = "custom";
	
	public static String FILE_FASTQ1 = "fastq1";
	public static String FILE_FASTQ2 = "fastq2";
	public static String FILE_BAM = "bam";
	public static String FILE_BAI = "bai";
	public static String FILE_SAM = "sam";
	public static String FILE_METRICS = "pdf";
	
	public static String TARGET_AGILENT = "agilent";
	public static String TARGET_NIMBLEGEN = "nimblegen";
	public static String TARGET_TRUSEQ = "truseq";
	public static String TARGET_CUSTOM = "custom";
	public static String TARGET_ALL = "all";
	
	public static File templateDir = new File("/home/BioApps/tomatoFarmer/");
	
	public static HashSet<String> validTypes = new HashSet<String>() {
		private static final long serialVersionUID = 1L;

	{
		add(TFConstants.ANALYSIS_CUSTOM);
		add(TFConstants.ANALYSIS_EXOME_BWA);
		add(TFConstants.ANALYSIS_EXOME_NOVOALIGN);
		add(TFConstants.ANALYSIS_EXOME_ALIGN_BWA);
		add(TFConstants.ANALYSIS_EXOME_ALIGN_NOVOALIGN);
		add(TFConstants.ANALYSIS_EXOME_METRICS);
		add(TFConstants.ANALYSIS_EXOME_VARIANT_RAW);
	}};
	
	public static HashSet<String> validTargets = new HashSet<String>() {
		private static final long serialVersionUID = 1L;

	{
		add(TFConstants.TARGET_AGILENT);
		add(TFConstants.TARGET_NIMBLEGEN);
		add(TFConstants.TARGET_TRUSEQ);
		add(TFConstants.TARGET_CUSTOM);
		add(TFConstants.TARGET_ALL);
	}};
	
}
