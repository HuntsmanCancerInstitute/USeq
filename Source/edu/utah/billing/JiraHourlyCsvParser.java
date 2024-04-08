package edu.utah.billing;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.TreeMap;

import util.gen.IO;
import util.gen.Misc;

public class JiraHourlyCsvParser {

	private File jiraCsv = null;
	private boolean debug = false;
	private LinkedHashMap<String, Integer> headerKeyIndex = null;
	private TreeMap<String, ArrayList<JiraTicketSummary>> groupNameTickets = new TreeMap<String, ArrayList<JiraTicketSummary>>();
	private ArrayList<JiraTicketSummary> errorTickets = new ArrayList<JiraTicketSummary>();
	
	public JiraHourlyCsvParser (File jiraCsv, boolean debug) throws IOException {
		this.jiraCsv = jiraCsv;
		this.debug = debug;
		parseIt();
	}
	
	private void parseIt() throws IOException {
		IO.pl("\nParsing the Jira ticket report...");
		BufferedReader in = IO.fetchBufferedReader(jiraCsv);
		String line;
		String [] cells;

		while ((line = in.readLine()) != null) {
			//empty?
			line = line.trim();
			if (line.length()== 0) continue;
			cells = Misc.splitCsvLine(line);
			if (cells[0].contains("Issue Key")) parseHeader(cells);
			else {
				//data line, add to group
				if (headerKeyIndex == null) throw new IOException ("Failed to parse a header line from "+jiraCsv);
				JiraTicketSummary jts = new JiraTicketSummary(cells, line, headerKeyIndex);
				//any errors?
				if (jts.getErrors().size()!=0) errorTickets.add(jts);
				else {
					String groupName = jts.getGroupToBill();
					if (groupName == null || groupName.length()==0) groupName = "NA";
					ArrayList<JiraTicketSummary> al = groupNameTickets.get(groupName);
					if (al == null) {
						al = new ArrayList<JiraTicketSummary>();
						groupNameTickets.put(groupName, al);
					}
					al.add(jts);					
				}
				if (jts.getJiraTicket()==null) {
					Misc.printErrAndExit("Missing jt#: ");
				}
			}
		}
		in.close();
		
		//any errors?
		if (errorTickets.size()!=0) {
			IO.pl("\nErrors with the following Jira tickets, correct and restart:");
			for (JiraTicketSummary jts: errorTickets) {
				for (String e: jts.getErrors()) IO.pl("\t"+e);
			}
			System.exit(1);
		}
		
		//check groups
		boolean errorsFound = false;
		for (String groupName: groupNameTickets.keySet()) {
			if (groupName.equals("NA")) continue;
			int numFTEFound = 0;
			int numHourlyFound = 0;
			int numInfrastructureFound = 0;
			for (JiraTicketSummary jts: groupNameTickets.get(groupName)) {
				if (jts.getWorkType().equals("FTE")) numFTEFound++;
				else if (jts.getWorkType().equals("Hourly")) numHourlyFound++;
				else if (jts.getWorkType().equals("Infrastructure")) numHourlyFound++;
			}
			//any infrastructure? no need to check this, should have thrown an error above
			if (numInfrastructureFound !=0) Misc.printErrAndExit("Non zero infrastructure count for "+groupName);
			
			//both FTE and Hourly?
			int total = numFTEFound + numHourlyFound;
			if (numFTEFound != 0 && numFTEFound != total) {
				errorsFound = true;
				IO.pl("\tMixed FTE and Hourly tickets were found for "+groupName);
				for (JiraTicketSummary jts: groupNameTickets.get(groupName)) {
					IO.pl("\t\t"+jts.getWorkType()+"\t"+jts.toString());
				}
			}
		}
		if (errorsFound) System.exit(1);
		
		//any NA?
		if (groupNameTickets.containsKey("NA")) {
			IO.pl("\n\tThe following are missing an Account and were labled Infrastructure, checkem!");
			for (JiraTicketSummary jts: groupNameTickets.get("NA")) {
				IO.pl("\t"+jts.toString());
				
			}
		}
		
		//print out the FTE
		IO.pl("\n\tThe following are FTE and will be excluded from hourly billing, checkem!");
		ArrayList<String> groupsToRemove = new ArrayList<String>();
		for (String groupName: groupNameTickets.keySet()) {
			if (groupName.equals("NA")) {
				groupsToRemove.add("NA");
				continue;
			}
			boolean fteGroup = false;

			for (JiraTicketSummary jts: groupNameTickets.get(groupName)) {
				if (jts.getWorkType().equals("FTE")) {
					IO.pl("\t"+jts.toString());
					fteGroup = true;
				}
			}
			if (fteGroup) groupsToRemove.add(groupName);
		}
		
		//sum the hours for FTE and remove them from the jira ticket map
		for (String gn: groupsToRemove) {
			if (gn.equals("NA")== false) {
				float totalHours = 0;
				for (JiraTicketSummary jts: groupNameTickets.get(gn)) totalHours+= Float.parseFloat(jts.getHoursString());
				IO.pl("\t\t"+gn+"\t"+totalHours);
			}
			groupNameTickets.remove(gn);
		}
	}
	
	
	private void parseHeader(String[] cells) throws IOException {
		if (debug) IO.pl("\nParsing Jira Report Header:");
		headerKeyIndex = new LinkedHashMap<String, Integer>();
		for (int i=0; i< cells.length; i++) {
			if (headerKeyIndex.containsKey(cells[i]) == false) {
				headerKeyIndex.put(cells[i], i);
				if (debug) IO.pl("\t"+cells[i]+"\t"+i);
			}
			else if (debug) IO.pl("\t"+cells[i]+"\tDuplicate header key skipping"+i);
		}

		//check it contains the required cells
		String[] toFind = {"CBI - Work Type", "Account Name", "Issue Key", "Full name", "Billed Hours", "Issue summary", "Work Description"};

		for (String tf: toFind) {
			if (headerKeyIndex.containsKey(tf) == false) throw new IOException("Failed to find the '"+tf+"' header key in "+headerKeyIndex);
		}
	}

	public TreeMap<String, ArrayList<JiraTicketSummary>> getGroupNameTickets() {
		return groupNameTickets;
	}
	
	
}
