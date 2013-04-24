package edu.utah.seq.vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
//##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
//##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
//##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
//##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
//##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
//##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
//##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
//##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
//##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
//##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
//##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
//##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
//##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
//##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
//##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
//##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
//##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
//##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
//##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
//##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
//##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
//##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model">
//##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">

public class VCFInfo {	
	private HashMap<String, String> hashInfoString = new HashMap<String,String>();
	
	public static String UNMODIFIED = "UNMODIFIED";
	public static String CLEAN = "CLEAN";
	
	
	private String infoString = null;
	

	public VCFInfo() {
		// TODO Auto-generated constructor stub
	}
	
	public void parseInfoGatk(String infoString) {
		String[] entries = infoString.split(";");
		for (String e: entries) {
			//create simpleInfoString
			this.infoString = infoString;
					
			//Get key-par, account for flag values
			String[] keyValue = e.split("=",2);
			String key = keyValue[0];
			String value = null;
			
			if (keyValue.length == 1) {
				value = "true";
			} else {
				value = keyValue[1];
			}
			
			hashInfoString.put(key,value);
		}
	}
	
	/** This method replaces the raw info string, which is built simply by appending added info fields */
	public void overwriteInfoString(String infoString) {
		this.infoString = infoString;
	}
	
	/** This method returns the raw info string, which is built simply by appending added info fields */
	public String getInfoString() {
		return this.infoString;
	}
	
	/** This method adds a new piece of data to the info object. The data will be added to the end of the raw
	 * info string and can be used when buiding custom info strings.
	 * @param identifier  name of info field
	 * @param info        value of info field
	 */
	public void addInfo(String identifier,String info) {
		this.hashInfoString.put(identifier, info);
		this.infoString += ";" + identifier + "=" + info;
	}
	
	
	/** Checks for the existence of a peice of inforamtion.  Returns true if the info object has data for the field
	 * and false if it doesn't
	 * @param entry
	 * @return
	 */
	public boolean doesInfoEntryExist(String entry) {
		if (hashInfoString.containsKey(entry)) {
			return true;
		} else {
			return false;
		}
	}
	
	/** Returns value of key.  If key doesn't exist, then an empty string is returned. If a style is 
	 * specified, certain strings will be reformatted */
	public String getInfo(String key, String style) {
		if (doesInfoEntryExist(key)) {
			String unmodified = hashInfoString.get(key);
			return checkForMods(key,unmodified,style);
		} else {
			return "";
		}
	}
	
	public float getInfoFloat(String key) throws Exception{
		String val = hashInfoString.get(key);
		if (val == null) throw new Exception ("Info key not found? -> "+key);
		return Float.parseFloat(val);
	}
	
	private String checkForMods(String key, String value, String style) {
		String moddedValue;
		if (style.equals(VCFInfo.CLEAN) && key.equals("SIFT")) {
			moddedValue = String.valueOf(1-Float.parseFloat(value));
		} else if (style.equals(VCFInfo.CLEAN) && key.equals("VarDesc")) {
			moddedValue = value.split(",")[0];
		} else {
			moddedValue = value;
		}
		return moddedValue;
	}

	/** This method builds a custom info string by adding values from the info fields listed in infoToAdd.
	 * 
	 * @param infoToAdd    Info fields that will make up the infostring
	 * @return             custom infostring
	 */
	public String buildInfoString(ArrayList<String> infoToAdd, String style) {
		StringBuilder infoString = new StringBuilder(""); 
		for (String info: infoToAdd) {
			if (hashInfoString.containsKey(info)) {
				String value = hashInfoString.get(info);
				if (value.equals("true")) {
					infoString.append(";" + info);
				} else {
					infoString.append(";" + info + "=" + getInfo(info,style));
				}
			}
		}
		return infoString.toString().replaceFirst(";","");
	}
	
	/**This method builds an annotation string for an output table.
	 * 
	 */
	public String buildInfoForTable(ArrayList<String> infoToAdd, String style, VCFComments comments) {
		StringBuilder infoString = new StringBuilder("");
		for (String info: infoToAdd) {
			if (hashInfoString.containsKey(info)) {
				if (comments.getFormat().get(info).equals("Flag")) {
					infoString.append("\t1");
				} else {
					infoString.append("\t" + getInfo(info,style));
				}
			}else {
				if (comments.getFormat().get(info).equals("Flag")) {
					infoString.append("\t0");
				} else {
					infoString.append("\tNA");
				}
			}
		}
		return infoString.toString().trim();
	}
	
	/** This method creates a to add list from a to skip list.  Its not efficient to parse the header for each record, so this run just once.
	 * 
	 * @param infoToSkip   Info fields that will not be part of the infostring
	 * @param vcfComments  Comments section of the VCF file
	 * @return             custom infostring
	 */
	public static ArrayList<String> buildToAddFromToSkip(ArrayList<String> infoToSkip,String[] vcfComments) {
		HashSet<String> skipSet = new HashSet<String>(infoToSkip);
		ArrayList<String> infoToAdd = new ArrayList<String>();
		Pattern p = Pattern.compile("##INFO=<ID=(.+?),.+");
		for (String comment: vcfComments) {
			Matcher m = p.matcher(comment);
			if (m.matches() && !skipSet.contains(m.group(1))) {
				infoToAdd.add(m.group(1));
			}
		}
		return infoToAdd;
	}


}
