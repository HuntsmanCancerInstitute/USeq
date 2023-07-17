package edu.utah.billing;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**
 * Generates a billing report for Cancer Bioinformatics from hourly and cloud expenses.
 */
public class CBiBilling {

	//fields
	private File jiraReportFile;
	private File wafDirectory; 
	private File outputDirectory;
	private File masterAcountInfo;
	private File cloudReportsDirectory;
	private File awsAccountsFile;
	
	//internal fields
	private CarahsoftXlsxParser carahsoftParser;
	private WafXlsxParser hourly = null;
	private WafXlsxParser cloud = null;
	private MasterAccountInfoParser masterAccountInfo = null;
	private AwsXlsxAccountParser awsXlsxAccountParser;
	private boolean debug = false;
	private LinkedHashMap<String, Integer> headerKeyIndex = null;
	private TreeMap<String, ArrayList<JiraTicketSummary>> groupNameTickets = new TreeMap<String, ArrayList<JiraTicketSummary>>();
	private ArrayList<JiraTicketSummary> errorTickets = new ArrayList<JiraTicketSummary>();
	private ArrayList<String> accountNumbersCloudBilled = new ArrayList<String>();
	private ArrayList<String> nonCancerBilling = new ArrayList<String>();
	private float totalHours = 0;
	private float totalExpenses = 0;
	private float firstTierPricePerHour = 70f;
	private float secondTierPricePerHour = 100f;
	private int maxHoursForFirstTier = 70;
	private float hscCostSharing = 356f;
	private float totalHscBilling = 0f;
	private String date = null;

	//constructor
	public CBiBilling (String[] args) throws Exception{
		processArgs (args);
		
		parseWafs();

		parseJiraHours();
		
		parseMasterAccountInfo();
		
		parseAWSAccounts();
		
		printInvoices();
		
		masterAccountInfo.saveUpdatedInfo();
		
		IO.pl("\nComplete!");
		
	}


	private void printInvoices() {
		IO.pl("\nPrinting Invoices...");
		HashMap<String, Float> groupNameHoursBilled = this.masterAccountInfo.getGroupNameHoursBilled();
		
		//two sources of billing, AWS and Hourly
		//Hourly
		for (String groupName: groupNameTickets.keySet()) {
			ArrayList<String> txtOut = new ArrayList<String>();
			IO.pl("GroupName:\t"+groupName);
			txtOut.add("HCI Cancer Bioinformatics Invoice\t\t"+date+"\n\n"+groupName+"\n");
			//pull Aliases
			HashSet<String> nameAliases = masterAccountInfo.getGroupNameAliases().get(groupName);
			//pull Hourly WAF
			ArrayList<String[]> wafLines = null;
			for (String a: nameAliases) {
				wafLines = hourly.getGroupNameWafLines().get(a);
				if (wafLines != null) break;
			}
			if (wafLines == null) IO.pl("HourlyWAFLine:\tNo Hourly WAF");
			else {
				IO.pl("HourlyWAFHeader:\t"+hourly.getHeaderTabbed());
				for (String[] wl: wafLines) IO.pl("HourlyWAFLine:\t"+Misc.stringArrayToString(wl, "\t"));
			}
			float totalHoursToBill = 0;
			IO.pl("HourlyJiraHeader:\t"+JiraTicketSummary.getToStringHeader());
			txtOut.add("Hourly Billing:\n\n"+JiraTicketSummary.getToStringHeader());
			for (JiraTicketSummary jt: groupNameTickets.get(groupName)) {
				IO.pl("HourlyJiraLine:\t"+jt.toString());
				totalHoursToBill += Float.parseFloat(jt.getHoursString());
				txtOut.add(jt.toString());
			}
			
			//pull the prior # of hours billed
			Float pastTotal = null;
			String groupNameInMaster = null;
			for (String a: nameAliases) {
				pastTotal = groupNameHoursBilled.get(a);
				if (pastTotal!=null) {
					groupNameInMaster = a;
					break;
				}
			}
			
			//calculate the expense, either $70/hr or $100/hr
			float pricePerHour = firstTierPricePerHour;
			if (pastTotal > maxHoursForFirstTier) pricePerHour = secondTierPricePerHour;
			float hourlyExpenses = totalHoursToBill* pricePerHour;
			IO.pl("HourlyTotal:\t$"+Num.formatNumber(hourlyExpenses, 2)+ "\t("+ totalHoursToBill+ " x $"+pricePerHour+"/hr)");
			txtOut.add("\t$"+Num.formatNumber(hourlyExpenses, 2)+ "\t("+ totalHoursToBill+ " x $"+pricePerHour+"/hr)\tHourly Expenses\n");
			totalHours += totalHoursToBill;
			
			//update the master
			groupNameHoursBilled.put(groupNameInMaster, new Float (pastTotal.floatValue()+ totalHoursToBill));
			
			//pull non cancer
			String cancerStatus = masterAccountInfo.getGroupNameCancerStatus().get(groupName);
			if (cancerStatus.contains("Non")) {
				float additionalBilling = hscCostSharing * totalHoursToBill;
				nonCancerBilling.add(groupName+"\t"+Num.formatNumber(totalHoursToBill,3)+"\t$"+ Num.formatNumber(additionalBilling, 2));
				totalHscBilling+= additionalBilling;
			}
				
			
			//look for Aws billing
			float cloudExpenses = printCloudInvoice(nameAliases,txtOut);
			IO.pl("TotalExpenses:\t$"+Num.formatNumber(hourlyExpenses+cloudExpenses, 2));
			IO.pl();
			totalExpenses += (hourlyExpenses+cloudExpenses);
			txtOut.add("\n$"+Num.formatNumber(hourlyExpenses+cloudExpenses, 2)+"\tTotal Expenses");
			
			//write out the txt file details
			String fileName = Misc.COMMA_WHITESPACE.matcher(groupName).replaceAll("_")+".txt";
			IO.writeString(Misc.stringArrayListToString(txtOut, "\n"), new File (outputDirectory, fileName));
			//if (groupName.contains("Ulrich")) Misc.printErrAndExit(Misc.stringArrayListToString(txtOut, "\n"));
		}
		
		//Cloud for those not printed with the hourly
		printJustCloudInvoices();
		
		IO.pl("TotalHoursBilled:\t"+ Num.formatNumber(totalHours, 2));
		IO.pl("TotalExpensesBilled:\t$"+ Num.formatNumber(totalExpenses, 2));
		
		printHSCBilling();
		
		
	}
	
	private void printHSCBilling() {
		IO.pl("\nHSC Invoice...");
		IO.pl("GroupName\tHoursBilled\tHSCShare($"+(int)hscCostSharing+"/hour)");
		for (String gl: nonCancerBilling) IO.pl(gl);
		IO.pl("\nTotalHSCCostSharing:\t$"+ Num.formatNumber(totalHscBilling, 2));
	}


	private void printJustCloudInvoices() {
		//remove any account numbers that were already printed with the hourly invoices
		TreeMap<String, Float> accountNumberTotals = carahsoftParser.getAwsAccountNumberTotals();
		for (String ac: accountNumbersCloudBilled) accountNumberTotals.remove(ac);

		//any remaining?
		if (accountNumberTotals.size()!=0) {
			TreeMap<String, String> awsAccountGroupName = awsXlsxAccountParser.getAwsAccountGroupName();
			TreeMap<String, ArrayList<String>> groupNameBillingInfo = new TreeMap<String, ArrayList<String>>();

			//for each numbered account with a bill
			for (String accountNumber: accountNumberTotals.keySet()) {

				//pull the groupName for the numberedAccount
				String groupName = awsAccountGroupName.get(accountNumber);

				ArrayList<String> billInfo = groupNameBillingInfo.get(groupName);
				if (billInfo == null) {
					billInfo = new ArrayList<String>();
					groupNameBillingInfo.put(groupName, billInfo);
				}
				float cost = accountNumberTotals.get(accountNumber);
				billInfo.add("AWS Acc# "+accountNumber+"\t$"+Num.formatNumber(cost, 2));
			}

			//for each group 
			TreeMap<String, ArrayList<String[]>> groupNameCloudWafs = cloud.getGroupNameWafLines();
			for (String groupName: groupNameBillingInfo.keySet()) {
				ArrayList<String> txtOut = new ArrayList<String>();
				IO.pl("GroupName:\t"+groupName);
				txtOut.add("HCI Cancer Bioinformatics Invoice\t\t"+date+"\n\n"+groupName+"\n\nCloud Billing:\n");
				
				if (groupNameCloudWafs.containsKey(groupName) == false) IO.pl("CloudWAFLine:\tNo Cloud WAF");
				else {
					IO.pl("CloudWAFHeader:\t"+cloud.getHeaderTabbed());
					for (String[] cw: groupNameCloudWafs.get(groupName)) {
						IO.pl("CloudWAFLine:\t"+Misc.stringArrayToString(cw, "\t"));
					}
				}
				float totalAwsCost = 0f;
				for (String line : groupNameBillingInfo.get(groupName)) {
					String[] tokens = Misc.TAB.split(line);
					float cost = Float.parseFloat(tokens[1].substring(1));
					totalAwsCost+= cost;
					IO.pl("CloudBilling:\t"+ line);
					txtOut.add(line);
				}
				String totalString = Num.formatNumber(totalAwsCost, 2);
				IO.p("CloudTotal:\t$"+totalString+"\n");
				IO.pl("TotalExpenses:\t$"+totalString+"\n");
				totalExpenses += totalAwsCost;
				
				txtOut.add("\t$"+totalString+"\tCloud Expenses\n\n");
				txtOut.add("$"+totalString+"\tTotal Expenses");
				
				//write out the txt file details
				String fileName = Misc.COMMA_WHITESPACE.matcher(groupName).replaceAll("_")+".txt";
				IO.writeString(Misc.stringArrayListToString(txtOut, "\n"), new File (outputDirectory, fileName));
				
			}
		}
	}

	
	private float printCloudInvoice(HashSet<String> groupAliases, ArrayList<String> txtOut) {
		float totalAwsCost = 0f;
		
		//any cloud reports? often these aren't available for a given month
		if (carahsoftParser.getAwsAccountNumberTotals().size()!=0) {
			TreeMap<String, Float> awsAccountNumberTotals = carahsoftParser.getAwsAccountNumberTotals();
			TreeMap<String, String> awsAccountGroupName = awsXlsxAccountParser.getAwsAccountGroupName();

			
			//for each numbered account with a bill
			StringBuilder out = new StringBuilder();
			txtOut.add("Cloud Billing:\n");
			for (String accountNumber: awsAccountNumberTotals.keySet()) {
				
				//pull the groupName for the numberedAccount
				String groupName = awsAccountGroupName.get(accountNumber);
				
				//is this in the groupAliases?
				if (groupAliases.contains(groupName)) {
					float cost = awsAccountNumberTotals.get(accountNumber);
					totalAwsCost += cost;
					out.append("CloudBilling:\tAWS Acc# "+accountNumber+"\t$"+Num.formatNumber(cost, 2)+"\n");
					txtOut.add("AWS Acc# "+accountNumber+"\t$"+Num.formatNumber(cost, 2));
					accountNumbersCloudBilled.add(accountNumber);
				}
			}
			if (totalAwsCost != 0f) {
				out.append("CloudTotal:\t$"+Num.formatNumber(totalAwsCost, 2)+"\n");
				txtOut.add("\t$"+Num.formatNumber(totalAwsCost, 2)+"\tCloud Expenses\n");
				//check for the WAFs
				//for each groupName alias, look for cloudWAFlines
				TreeMap<String, ArrayList<String[]>> cloudWafs = cloud.getGroupNameWafLines();
				
				boolean wafFound = false;
				for (String gn: groupAliases) {
					if (cloudWafs.containsKey(gn)) {
						wafFound = true;
						IO.pl("CloudWAFHeader:\t"+cloud.getHeaderTabbed());
						for (String[] cw: cloudWafs.get(gn)) {
							IO.pl("CloudWAFLine:\t"+Misc.stringArrayToString(cw, "\t"));
						}
					}
				}
				if (wafFound ==false) IO.pl("CloudWAFLine:\tNo Cloud WAF");
				IO.p(out);
			}
		}
		return totalAwsCost;
	}


	private void parseAWSAccounts() {
		IO.pl("\nParsing AWS accounts...");
		//any cloud reports?
		if (cloudReportsDirectory != null) {
			carahsoftParser = new CarahsoftXlsxParser(cloudReportsDirectory, debug);
			//anything parsed
			if (carahsoftParser.getAwsAccountNumberTotals().size()!=0) {
				//parse the Aws Accounts
				awsXlsxAccountParser = new AwsXlsxAccountParser(awsAccountsFile, debug);
				//for each carasoft account charge, look to see if it's in the aws accounts
				TreeMap<String, String> awsAccountGroupName = awsXlsxAccountParser.getAwsAccountGroupName();
				TreeMap<String, Float> awsAccountNumberTotals = carahsoftParser.getAwsAccountNumberTotals();
				ArrayList<String> missingAccountNumber = new ArrayList<String>();
				for (String awsAccountNumber: awsAccountNumberTotals.keySet()) {
					if (awsAccountGroupName.containsKey(awsAccountNumber)==false) missingAccountNumber.add(awsAccountNumber);
				}
				
				//any missing account numbers
				if (missingAccountNumber.size()!=0) {
					for (String acc: missingAccountNumber) IO.el("\tMissing "+acc+" in "+awsAccountsFile+", correct and restart.");
					System.exit(1);
				}
				
				//check for the WAFs
				//for each account number from Carahsoft with charges
				TreeMap<String, ArrayList<String[]>> cloudWafs = cloud.getGroupNameWafLines();
				for (String awsAccountNumber: awsAccountNumberTotals.keySet()) {
					//fetch the name from the aws accounts
					String groupName = awsAccountGroupName.get(awsAccountNumber);
					//look for the name in the cloud WAFs
					if (cloudWafs.containsKey(groupName) == false) {
						IO.pl("\tMissing Cloud WAF entry for -> "+groupName);
					}
				}
			}
		}
	}
	
	private void parseMasterAccountInfo() {
		
		masterAccountInfo = new MasterAccountInfoParser(masterAcountInfo, debug);
		
		ArrayList<String> missingGroups = new ArrayList<String>();
		HashMap<String, HashSet<String>> aliasesMap = masterAccountInfo.getGroupNameAliases();
		TreeMap<String, ArrayList<String[]>> userWafs = hourly.getGroupNameWafLines();
		
		//check jira tickets
		IO.pl("\nChecking Aliases and if Jira tickets have an Hourly WAF... ");
		
		//the only groups left here are Hourly ticket holders
		boolean missingAlias = false;
		for (String jg: groupNameTickets.keySet()) {
			if (aliasesMap.containsKey(jg) == false) {
				missingGroups.add("Missing alias entry for -> "+jg);
				missingAlias = true;
			}
			else {
				//check for WAF entry
				boolean found = false;
				for (String a: aliasesMap.get(jg)) {
					if (userWafs.containsKey(a)) {
						found = true;
						break;
					}
				}
				if (found == false) missingGroups.add("Missing Hourly WAF entry for -> "+jg);
			}
		}
		for (String me: missingGroups) IO.pl("\t"+me);
		if (missingAlias) Misc.printErrAndExit("\t\tCorrect missing alias and restart");
	}


/*
	private void parseAliases() {
		
		aliases = new AliasXlsxParser(groupNameAliasesFile, debug);
		
		ArrayList<String> missingGroups = new ArrayList<String>();
		TreeMap<String, HashSet<String>> aliasesMap = aliases.getGroupNameAliases();
		TreeMap<String, ArrayList<String[]>> userWafs = hourly.getGroupNameWafLines();
		
		//check jira tickets
		IO.pl("\nChecking Aliases and if Jira tickets have an Hourly WAF... ");
		
		//the only groups left here are Hourly ticket holders
		boolean missingAlias = false;
		for (String jg: groupNameTickets.keySet()) {
			if (aliasesMap.containsKey(jg) == false) {
				missingGroups.add("Missing alias entry for -> "+jg);
				missingAlias = true;
			}
			else {
				//check for WAF entry
				boolean found = false;
				for (String a: aliasesMap.get(jg)) {
					if (userWafs.containsKey(a)) {
						found = true;
						break;
					}
				}
				if (found == false) missingGroups.add("Missing Hourly WAF entry for -> "+jg);
			}
		}
		for (String me: missingGroups) IO.pl("\t"+me);
		if (missingAlias) Misc.printErrAndExit("\t\tCorrect missing alias and restart");
	}
	*/


	private void parseJiraHours() throws IOException {
		IO.pl("\nParsing the Jira ticket report...");
		BufferedReader in = IO.fetchBufferedReader(jiraReportFile);
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
				if (headerKeyIndex == null) throw new IOException ("Failed to parse a header line from "+jiraReportFile);
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
					IO.pl("\t\t"+jts.toString());
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
		String[] toFind = {"Work Type", "Account Name", "Issue Key", "Full name", "Billed Hours", "Issue summary", "Work Description"};

		for (String tf: toFind) {
			if (headerKeyIndex.containsKey(tf) == false) throw new IOException("Failed to find the '"+tf+"' header key in "+headerKeyIndex);
		}
	}


	private void parseWafs() {
		IO.pl("\nParsing the WAF tracking spreadsheets...");
		File[] xlsFiles = IO.extractFiles(wafDirectory, ".xlsx");
		for (File f: xlsFiles) {
			if (f.getName().contains("WAF") && f.getName().startsWith("~")== false) {
				if (f.getName().contains("Cloud")) cloud = new WafXlsxParser(f, debug);
				else if (f.getName().contains("Cloud") == false) hourly = new WafXlsxParser(f, debug);
			}
		}
		if (cloud == null || hourly == null) Misc.printErrAndExit("\nFailed to parse both an hourly and cloud WAF tracking schedule xlsx file.");

		
	}







	public static void main(String[] args) throws Exception {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CBiBilling(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'j': jiraReportFile = new File (args[i+1]); i++; break;
					case 'w': wafDirectory = new File (args[i+1]); i++; break;
					case 'o': outputDirectory = new File (args[i+1]); i++; break;
					case 'm': masterAcountInfo = new File (args[i+1]); i++; break;
					case 'c': cloudReportsDirectory = new File (args[i+1]); i++; break;
					case 'a': awsAccountsFile = new File (args[i+1]); i++; break;
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (jiraReportFile == null || jiraReportFile.exists() == false || jiraReportFile.getName().endsWith(".csv") == false) Misc.printExit("\nError: cannot find your Jira report ending in .csv ! "+ jiraReportFile);
		if (masterAcountInfo == null || masterAcountInfo.exists() == false) Misc.printExit("\nError: cannot find your group name aliases xlsx file! "+ masterAcountInfo);
		if (wafDirectory == null || wafDirectory.exists() == false) Misc.printExit("\nError: cannot find your user account info file ! "+ wafDirectory);
		if (awsAccountsFile == null || awsAccountsFile.exists() == false) Misc.printExit("\nError: cannot find your AWS account xlsx file ! "+ awsAccountsFile);
		if (outputDirectory == null ) Misc.printExit("\nError: cannot find your output directory ! "+ outputDirectory);
		outputDirectory.mkdirs();
		date = Misc.getDateTime();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Correlate:    Nov 2008                            **\n" +
				"**************************************************************************************\n" +
				"Calculates all pair-wise Pearson correlation coefficients (r) and if indicated will\n" +
				"perform a hierarchical clustering on the files.\n\n"+				

				"Parameters:\n" +
				"-d The full path directory text containing serialized java float[] files (xxx.celp\n"+
				"      see CelProcessor app).\n"+
				"-a Files provided are float[][] files (xxx.cela) and need to be collapsed to float[]\n"+
				"-c Cluster files.\n\n" +

				"Example: java -Xmx256M -jar pathTo/T2/Apps/Correlate -d /Mango/PCels/ -c -a\n\n" +


		"**************************************************************************************\n");		
	}	
}
