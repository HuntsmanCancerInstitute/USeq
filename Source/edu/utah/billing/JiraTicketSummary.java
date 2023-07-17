package edu.utah.billing;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import util.gen.IO;
import util.gen.Num;

public class JiraTicketSummary {
	
	//pull Work Type, Hourly | Infrastructure | FTE - cannot be blank
	private String workType = null;
	//pull Account Name, Judson-Torres, Robert Lab - might be blank for Infrastructure
	private String groupToBill = null;
	//pull Issue Key, BSD-642 - cannot be blank
	private String jiraTicket = null;
	//pull Full name, Timothy Parnell - cannot be blank
	private String analystName = null;
	//pull Billed Hours, 2 - cannot be blank
	private String hoursString = null;
	//pull Issue summary, ChIP Seq data analysis - can be blank
	private String issueSummary = null;
	//pull Work Description, Evaluate new results. Generate new plots. Annotate peaks. Post to GNomEx and Email. - can be blank
	private String workPerformed = null;
	//should be empty!
	private ArrayList<String> errors = new ArrayList<String>();
	
	public JiraTicketSummary (String[] cells, String line, LinkedHashMap<String, Integer> headerKeyIndex) {
		
		
		//pull Work Type, Hourly | Infrastructure | FTE - cannot be blank
		workType = cells[headerKeyIndex.get("Work Type")].trim();
		if (workType.length()==0) errors.add("Error: missing 'Work Type' in -> "+line);
		
		//pull Account Name, Judson-Torres, Robert Lab - might be blank for Infrastructure
		groupToBill = cells[headerKeyIndex.get("Account Name")].trim();
		if (workType.startsWith("Infrastructure")== false && groupToBill.length()==0) errors.add("Error: missing 'Account Name' for a non-infrastructure job in -> "+line);
		if (workType.startsWith("Infrastructure") && groupToBill.length()!=0) errors.add("Error: Infrastructure job has an Account Name -> "+line);
	
		//pull Issue Key, BSD-642 - cannot be blank
		jiraTicket = cells[headerKeyIndex.get("Issue Key")].trim();
		if (jiraTicket.length()==0) errors.add("Error: missing 'Issue Key' in -> "+line);
		
		//pull Full name, Timothy Parnell - cannot be blank
		analystName = cells[headerKeyIndex.get("Full name")].trim();
		if (analystName.length()==0) errors.add("Error: missing 'Full name' in -> "+line);
		
		//pull Billed Hours, 2 - cannot be blank
		hoursString = cells[headerKeyIndex.get("Billed Hours")].trim();
		if (hoursString.length()==0) errors.add("Error: missing 'Billed Hours' in -> "+line);
		
		//pull Issue summary, ChIP Seq data analysis - can be blank
		issueSummary = cells[headerKeyIndex.get("Issue summary")].trim();
		
		//pull Work Description, Evaluate new results. Generate new plots. Annotate peaks. Post to GNomEx and Email. - can be blank
		workPerformed = cells[headerKeyIndex.get("Work Description")].trim();
		
	}
	
	public static final String getToStringHeader() {
		return "JiraTicketID\tGroup\tAnalyst\tHours\tTicketSummary\tWorkPerformed";
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(jiraTicket); sb.append("\t");
		sb.append(groupToBill); sb.append("\t");
		sb.append(analystName); sb.append("\t");
		double hours = Double.parseDouble(hoursString);
		sb.append(Num.formatNumber(hours, 2)); sb.append("\t");
		sb.append(issueSummary); sb.append("\t");
		sb.append(workPerformed);
		return sb.toString();
	}

	public String getWorkType() {
		return workType;
	}

	public String getGroupToBill() {
		return groupToBill;
	}

	public String getJiraTicket() {
		return jiraTicket;
	}

	public String getAnalystName() {
		return analystName;
	}

	public String getHoursString() {
		return hoursString;
	}

	public String getIssueSummary() {
		return issueSummary;
	}

	public String getWorkPerformed() {
		return workPerformed;
	}

	public ArrayList<String> getErrors() {
		return errors;
	}

}
