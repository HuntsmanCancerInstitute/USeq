package edu.utah.seq.vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VCFComments {
	//This class only handles the info section right now

	private ArrayList<String> preInfo = new ArrayList<String>();
	private ArrayList<String> postInfo = new ArrayList<String>();
	private ArrayList<String> infoOrder = new ArrayList<String>();
	private HashMap<String,String> infoData = new HashMap<String,String>();
	private Pattern p = Pattern.compile("##INFO=<ID=(.+?),.+");
	
	public VCFComments(ArrayList<String> comments) {
		setComments(comments.toArray(new String[0]));
	}
	
	/** Get comment array, restrict info lines to those in the list */
	public String[] getComments(ArrayList<String> infoFieldsToUse) {
		ArrayList<String> fullComments = new ArrayList<String>();
		fullComments.addAll(preInfo);
		for (String info: infoOrder) {
			if (infoFieldsToUse.contains(info)) {
				fullComments.add(infoData.get(info));
			}
			
		}
		fullComments.addAll(postInfo);
		
		return fullComments.toArray(new String[0]);
	}
	
	/** Get comment array, no arguments returns all info fields */
	public String[] getComments() {
		return getComments(infoOrder);
	}
	
	public void addInfo(String infofield) {
		Matcher m = p.matcher(infofield);
		if (!m.matches()) {
			System.out.println("Info header line is not properly formatted");
			System.exit(1);
		}
		
		infoOrder.add(m.group(1));
		infoData.put(m.group(1), infofield);	
	}
	
	public void setComments(String[] comments) {
		boolean found = false;
		for (String comment: comments) {
			Matcher m = p.matcher(comment);
			if (m.matches()) {
				found = true;
				infoOrder.add(m.group(1));
				infoData.put(m.group(1), comment);
			} else if (found) {
				postInfo.add(comment);
			} else {
				preInfo.add(comment);
			}
		}
	}
	
	public HashMap<String,String> getInfo() {
		return infoData;
	}
	
	public ArrayList<String> getInfoOrder() {
		return infoOrder;
	}
	

}
 