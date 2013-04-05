package edu.utah.seq.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VCFComments {
	//This class was made to handle info data, it be expanded to handle other fields.
	private ArrayList<String> nonStandardComments = new ArrayList<String>();
	
	//INFO related data
	private ArrayList<String> infoOrder = new ArrayList<String>();
	private HashMap<String,String> infoData = new HashMap<String,String>();
	private HashMap<String,String> infoDesc = new HashMap<String,String>();
	private HashMap<String,String> infoFormat = new HashMap<String,String>();
	private Pattern infoPattern = Pattern.compile("##INFO=<ID=(.+?),.+?Type=(.+?),Description=\"(.+?)\">");
	
	//FILTER related data
	private ArrayList<String> filterComments = new ArrayList<String>();
	private Pattern filterPattern = Pattern.compile("##FILTER.+");
	
	//FORMAT related data
	private ArrayList<String> formatComments = new ArrayList<String>();
	private Pattern formatPattern = Pattern.compile("##FORMAT.+");
	
	//COLUMN related data
	private String sampleString = null;
	private ArrayList<String> sampleList = new ArrayList<String>();
	private Pattern samplePattern = Pattern.compile("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t(.+?)");
	
	//VCF style
	private String styleString = null;
	private String style = null;
	private Pattern stylePattern = Pattern.compile("##annotationStyle=(.+?)");
	
	public VCFComments(ArrayList<String> comments) {
		setComments(comments.toArray(new String[0]));
	}
	
	/** Get comment array, restrict info lines to those in the list */
	public String[] getComments(ArrayList<String> infoFieldsToUse) {
		ArrayList<String> fullComments = new ArrayList<String>();
		fullComments.addAll(nonStandardComments);
		if (styleString != null) {
			fullComments.add(styleString);
		}
		for (String info: infoOrder) {
			if (infoFieldsToUse.contains(info)) {
				fullComments.add(infoData.get(info));
			}
		}
		fullComments.addAll(filterComments);
		fullComments.addAll(formatComments);
		fullComments.add(sampleString);
		
		return fullComments.toArray(new String[0]);
	}
	
	/** Get comment array, no arguments returns all info fields */
	public String[] getComments() {
		return getComments(infoOrder);
	}
	
	public void addInfo(String infofield) {
		Matcher m = infoPattern.matcher(infofield);
		if (!m.matches()) {
			System.out.println("Info header line is not properly formatted");
			System.exit(1);
		}
		infoOrder.add(m.group(1));
		infoData.put(m.group(1), infofield);	
	}
	
	public void setComments(String[] comments) {
		for (String comment: comments) {
			Matcher infoMatcher = infoPattern.matcher(comment);
			Matcher filterMatcher = filterPattern.matcher(comment);
			Matcher formatMatcher = formatPattern.matcher(comment);
			Matcher sampleMatcher = samplePattern.matcher(comment);
			Matcher styleMatcher = stylePattern.matcher(comment);
			if (infoMatcher.matches()) {
				String infoName = infoMatcher.group(1);
				infoOrder.add(infoName);
				infoData.put(infoName, comment);
				infoFormat.put(infoName, infoMatcher.group(2));
				infoDesc.put(infoName,infoMatcher.group(3));
			} else if (filterMatcher.matches()) {
				filterComments.add(comment);
			} else if (formatMatcher.matches()) {
				formatComments.add(comment);
			} else if (sampleMatcher.matches()) {
				sampleString = comment;
				String[] sampleNames = sampleMatcher.group(1).split("\t");
				for (String sn: sampleNames) {
					sampleList.add(sn);
				}
			} else if (styleMatcher.matches()) {
				styleString = comment;
				style = styleMatcher.group(1);
			} else {
				nonStandardComments.add(comment);
			}
		}	
	}
	
	public String getAnnotationStyle() {
		return style;
	}
	
	public void setAnnotationStyle(String style) {
		this.styleString = "##annotationStyle=" + style;
	}
	
	public ArrayList<String> getInfoDesc(ArrayList<String> infoToUse) {
		ArrayList<String> result = new ArrayList<String>();
		for (String info: infoToUse) {
			if (infoDesc.containsKey(info)) {
				result.add(infoDesc.get(info));
			}
		}
		
		return result;
	}
	
	public HashMap<String,String> getInfo() {
		return infoData;
	}
	
	public ArrayList<String> getInfoOrder() {
		return infoOrder;
	}
	
	public HashMap<String,String> getFormat() {
		return infoFormat;
	}
	
	public ArrayList<String> getSampleList() {
		return this.sampleList;
	}
	
	

}
 