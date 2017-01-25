package edu.utah.pysano.utils;

import java.util.HashSet;


public class Constants {
	public static HashSet<String> availClusters = new HashSet<String>() {
	{
		add("ember");
		add("kingspeak");
		add("kingspeak_20");
		add("kingspeak_24");
		add("timbuktu");
	}};
	
	public static HashSet<String> validMenu = new HashSet<String>() {
	{
		add("ugp_full");
		add("ugp_align");
		add("ugp_variant");
		add("metrics");
		add("ugp_align_metrics");
		add("ugp_haplotype_only");
	}};
	
	public static HashSet<String> validTargets = new HashSet<String>() {
	{
		add("AgilentAllExonV5UTR");
		add("AgilentAllExonV5");
		add("AgilentAllExonV4");
		add("AgilentAllExon50MB");
		add("NimbleGenEZCapV2");
		add("NimbleGenEZCapV3");
		add("TruSeq");
	}};
}
