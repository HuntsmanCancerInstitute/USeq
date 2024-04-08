package edu.utah.billing;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**
 * Generates a billing report for Cancer Bioinformatics from hourly and cloud expenses.
 */
public class CBiBilling2 {

	//fields
	private File jiraReportFile;
	private File wafDirectory; 
	private File outputDirectory;
	private File masterAcountInfo;
	private File cloudReportsDirectory;
	private File awsAccountsFile;
	private File expenseFile;
	private boolean debug = false;
	
	//internal fields
	private MasterAccountInfoParser2 masterAccountInfo = null;
	private HashMap<String, BillingGroup> aliasesBillingGroups = null;
	
	private WafXlsxParser2 hourlyWafParser = null;
	private WafXlsxParser2 cloudWafParser = null;
	
	private AwsXlsxAccountParser2 awsXlsxAccountParser = null;
	private TDSynnexXlsxParser tDSynnexXlsxParser = null;
	private MiscExpenseXlsxParser miscExpenseParser = null;

	private ArrayList<String> accountNumbersCloudBilled = new ArrayList<String>();
	private ArrayList<String> nonCancerBilling = new ArrayList<String>();
	
	static float firstTierPricePerHour = 70f;
	static float secondTierPricePerHour = 100f;
	static int maxHoursForFirstTier = 70;
	static float hscCostSharing = 356f;
	static float minimalBillingExpense = 5f;
	static String contactEmail = "hci-srbilling@hci.utah.edu and cbi@hci.utah.edu";
	static String cbiPolicyUrl = "https://uofuhealth.utah.edu/huntsman/shared-resources/gcb/cbi/cost";
	
	private int clients = 0;
	private int ccClientsBilled = 0;
	private int ncClientsBilled = 0;
	private float ccTotalHours = 0;
	private float ncTotalHours = 0;
	private float totalHourlyExpenses = 0;
	private float totalHSCExpenses = 0f;
	private float totalComputeExpenses = 0f;
	private ArrayList<BillingGroup> bgsNoHourlyWaf = new ArrayList<BillingGroup>();
	private ArrayList<BillingGroup> bgsNoCloudWaf = new ArrayList<BillingGroup>();
	private String date = null;

	//constructor
	public CBiBilling2 (String[] args) throws Exception{
		processArgs (args);
		
		createBillingGroups();
		
		parseWafs();
		
		//any AWS cloud compute stuff to parse
		if (cloudReportsDirectory != null) {
			parseAwsAccounts();
			parseAWSCloudAccountInvoices();
		}
		
		//any misc expenses, one offs like annual software licenses
		if (expenseFile != null) parseMiscExpenses();

		parseJiraHours();
		
		printInvoices();
		
		calculateSummaryStatistics();
		
		printSummaryStatsMissingWafs();
		
		masterAccountInfo.saveUpdatedInfo();
		
		IO.pl("\nComplete!");
		
	}


	private void printSummaryStatsMissingWafs() {
		
		IO.pl("Summary Statistics (CC CancerCenter, NC NonCC) ...");
		IO.pl("\tClients             :\t"+ clients);
		IO.pl("\tCCClientsBilled       :\t"+ ccClientsBilled);
		IO.pl("\tNCClientsBilled       :\t"+ ncClientsBilled);
		IO.pl("\tTotalCCClientHours    :\t"+ Num.formatNumber(ccTotalHours, 2));
		IO.pl("\tTotalNCClientHours    :\t"+ Num.formatNumber(ncTotalHours, 2));
		IO.pl("\tTotalHourlyExpenses :\t$"+ Num.formatNumber(totalHourlyExpenses, 2));
		IO.pl("\tTotalComputeExpenses:\t$"+ Num.formatNumber(totalComputeExpenses, 2));
		IO.pl("\tTotalHSCCostSharing :\t$"+ Num.formatNumber(totalHSCExpenses, 2));
		
		IO.pl("\nTotalBilling        :\t$"+ Num.formatNumber((totalHourlyExpenses+ totalHSCExpenses+ totalComputeExpenses), 2));
		
		IO.pl("\nClients Missing Current WAFs in the Tracking Spreadsheet! Check if these are in the UBox before adding to the contact list!");
		IO.pl("\nClients Missing Hourly WAFs (GroupName CancerMembership UBoxStatus):");
		for (BillingGroup bg: bgsNoHourlyWaf) IO.pl("\t"+bg.getGroupName()+"\t"+bg.isCancerCenterMember()+"\t?");
		IO.pl("\nClients Missing Cloud WAFs (GroupName CancerMember):");
		for (BillingGroup bg: bgsNoCloudWaf) IO.pl("\t"+bg.getGroupName()+"\t"+bg.isCancerCenterMember()+"\t?");
	}


	private void calculateSummaryStatistics() {
		
		for (BillingGroup bg: masterAccountInfo.getBillingGroups()) {
			if (bg.getTotalExpenses() > 0) {
				clients++;
				
				if (bg.getTotalExpenses() >= minimalBillingExpense) {
					if (bg.isCancerCenterMember()) {
						ccClientsBilled++;
						ccTotalHours+= bg.getTotalHoursToBill();
					}
					else {
						ncClientsBilled++;
						ncTotalHours+= bg.getTotalHoursToBill();
					}
					
					totalHourlyExpenses+= bg.getTotalHourlyExpenses();
					totalHSCExpenses+= bg.getAdditionalHourlyExpenses();
					totalComputeExpenses+= bg.getTotalComputeExpenses();
					// missing wafs?
					if (bg.getTotalHourlyExpenses()>0 && bg.getHourlyWafs().size()==0) bgsNoHourlyWaf.add(bg);
					if (bg.getTotalComputeExpenses()>0 && bg.getComputeWafs().size()==0) bgsNoCloudWaf.add(bg);
				}
			}
			
		}
	}


	private void parseJiraHours() throws IOException {
		JiraHourlyCsvParser jiraParser = new JiraHourlyCsvParser(jiraReportFile, debug);
		
		//add to BillingGroups and id those missing
		TreeMap<String, ArrayList<JiraTicketSummary>> groupNameTickets = jiraParser.getGroupNameTickets();
		
		ArrayList<String> missingGroupNames = new ArrayList<String>();
		for (String gn: groupNameTickets.keySet()) {
			BillingGroup bg = this.aliasesBillingGroups.get(gn);
			if (bg == null) missingGroupNames.add(gn);
			else bg.getJiraTickets().addAll(groupNameTickets.get(gn));
		}
		
		//any missing? if so exit.
		if (missingGroupNames.size()!=0) {
			IO.el("The following group names from the Jira Hourly parsing are missing from the aliases in the MasterAccountsInfo sheet, add them and restart:");
			for (String gn: missingGroupNames) IO.el("\t"+gn);
			System.exit(1);
		}
	}

	private void parseMiscExpenses() {
		IO.pl("\nParsing the miscellaneous compute expense spreadsheet...");
		//parse the xlsx spreadsheet
		miscExpenseParser = new MiscExpenseXlsxParser(expenseFile, debug);
		
		//add to BillingGroups and id those missing
		TreeMap<String, ArrayList<MiscExpense>> groupNameExpense = miscExpenseParser.getGroupNameExpense();
		ArrayList<String> missingGroupNames = new ArrayList<String>();
		for (String gn: groupNameExpense.keySet()) {
			BillingGroup bg = this.aliasesBillingGroups.get(gn);
			if (bg == null) missingGroupNames.add(gn);
			else bg.getMiscExpenses().addAll(groupNameExpense.get(gn));
		}
		
		//any missing? if so exit.
		if (missingGroupNames.size()!=0) {
			IO.el("The following group names from the MiscExpenseXlsx parsing are missing from the aliases in the MasterAccountsInfo sheet, add them and restart:");
			for (String gn: missingGroupNames) IO.el("\t"+gn);
			System.exit(1);
		}
		
	}

	private void parseAwsAccounts() {
		
		IO.pl("\nParsing the AWS Account Info spreadsheet...");
		awsXlsxAccountParser = new AwsXlsxAccountParser2(awsAccountsFile, debug);
		
		//check that all of the groupNames are in the billing groups
		ArrayList<String> missingNames = new ArrayList<String>();
		for (String groupName: awsXlsxAccountParser.getAwsAccountGroupName().values()) {
			BillingGroup bg = aliasesBillingGroups.get(groupName);
			if (bg == null) missingNames.add(groupName);
		}
		
		if (missingNames.size()!=0) {
			IO.el("The following group names from the AwsXlsxAccount parsing are missing from the aliases in the MasterAccountsInfo sheet, add them and restart:");
			for (String gn: missingNames) IO.el("\t"+gn);
			System.exit(1);
		}
	}

	private void createBillingGroups() throws IOException {
		IO.pl("\nParsing the Master Account Info spreadsheet...");
		
		// Parse the master account xlsx sheet
		masterAccountInfo = new MasterAccountInfoParser2(masterAcountInfo, debug);
		aliasesBillingGroups = masterAccountInfo.getAliasesBillingGroups();
	}

	private void printInvoices() throws IOException {
		IO.pl("\nPrinting Invoices...");
		
		//find BillingGroups with > minimumExpense, not billing tiny users, < $5
		ArrayList<BillingGroup> bgsToBillCancer = new ArrayList<BillingGroup>();
		ArrayList<BillingGroup> bgsToBillNonCancerWithHSCCostSharing = new ArrayList<BillingGroup>();
		ArrayList<BillingGroup> bgsToBillNonCancerNoHSCCostSharing = new ArrayList<BillingGroup>();
		for (BillingGroup bg: masterAccountInfo.getBillingGroups()) {
			bg.calculateTotalExpenses();
			if (bg.getTotalExpenses() >= minimalBillingExpense) {
				if (bg.isCancerCenterMember()) bgsToBillCancer.add(bg);
				else if (bg.getAdditionalHourlyExpenses()> 0) bgsToBillNonCancerWithHSCCostSharing.add(bg);
				else bgsToBillNonCancerNoHSCCostSharing.add(bg);
			}
		}
		
		
		//for HCI and HSC Billing
		IO.pl("\n####################### Cancer Center Members #########################\n");
		for (BillingGroup bg: bgsToBillCancer) {
			IO.pl("----------------------------------------------------------------\n");
			String invoice = bg.generateInvoice(date, hourlyWafParser.getHeaderTabbed(), cloudWafParser.getHeaderTabbed(), true, false);
			IO.pl(invoice);	
		}
		
		IO.pl("####################### Non Cancer Center Members No HSC Cost Sharing #########################\n");
		for (BillingGroup bg: bgsToBillNonCancerNoHSCCostSharing) {
			IO.pl("----------------------------------------------------------------\n");
			String invoice = bg.generateInvoice(date, hourlyWafParser.getHeaderTabbed(), cloudWafParser.getHeaderTabbed(), true, false);
			IO.pl(invoice);	
		}
		
		IO.pl("####################### Non Cancer Center Members With HSC Cost Sharing #########################\n");
		for (BillingGroup bg: bgsToBillNonCancerWithHSCCostSharing) {
			IO.pl("----------------------------------------------------------------\n");
			String invoice = bg.generateInvoice(date, hourlyWafParser.getHeaderTabbed(), cloudWafParser.getHeaderTabbed(), true, false);
			IO.pl(invoice);	
		}
		
		//print individual invoices
		bgsToBillCancer.addAll(bgsToBillNonCancerNoHSCCostSharing);
		bgsToBillCancer.addAll(bgsToBillNonCancerWithHSCCostSharing);
		for (BillingGroup bg: bgsToBillCancer) {
			String fileName = Misc.COMMA_WHITESPACE.matcher(bg.getGroupName()).replaceAll("_")+"_CBIInvoice_"+date+".txt";
			String invoice = bg.generateInvoice(date, hourlyWafParser.getHeaderTabbed(), cloudWafParser.getHeaderTabbed(), false, true);
			IO.writeString(invoice, new File (outputDirectory, fileName));
		}		
	}

	private void parseAWSCloudAccountInvoices() throws IOException {
		IO.pl("\nParsing TDSynnex AWS account invoices...");
		
		TreeMap<String, String> accountNumberGroupName = awsXlsxAccountParser.getAwsAccountGroupName();
		

		//for each billing group charged, add to it their AWS expenses
		ArrayList<String> missingNames = new ArrayList<String>();
		tDSynnexXlsxParser = new TDSynnexXlsxParser(cloudReportsDirectory, debug);
		TreeMap<String, Float> accountNumberTotalExpense= tDSynnexXlsxParser.getAwsAccountNumberTotals();
		for (String awsAccountNumber: accountNumberTotalExpense.keySet()) {
			//find the billing group name 
			String billingGroupName = accountNumberGroupName.get(awsAccountNumber);
			if (billingGroupName == null) missingNames.add(awsAccountNumber);
			
			else {
				BillingGroup bg = aliasesBillingGroups.get(billingGroupName);
				bg.getAwsAccountExpenses().add(new AwsAccountExpense(awsAccountNumber, accountNumberTotalExpense.get(awsAccountNumber)));
			}
		}
		if (missingNames.size()!=0) {
			IO.el("The following group names from the TDSynnex Xlsx billing are missing from the AWS Account Xlsx info sheet, add them and restart:");
			for (String gn: missingNames) IO.el("\t"+gn);
			System.exit(1);
		}
	}
	
	private void parseWafs() {
		IO.pl("\nParsing the WAF tracking spreadsheets...");
		File[] xlsFiles = IO.extractFiles(wafDirectory, ".xlsx");
		for (File f: xlsFiles) {
			String name = f.getName().toLowerCase();
			if (name.contains("waf") && f.getName().startsWith("~")== false) {
				if (name.contains("cloud")) cloudWafParser = new WafXlsxParser2(f, debug);
				else hourlyWafParser = new WafXlsxParser2(f, debug);

			}
		}
		
		if (cloudWafParser == null || hourlyWafParser == null) Misc.printErrAndExit("\nFailed to parse both an hourly and cloud WAF tracking schedule xlsx file.");
		
		//add cloud compute WAF lines to each Billing Group
		ArrayList<String> missingAliasCloud = new ArrayList<String>();
		for (String groupName: cloudWafParser.getGroupNameWafLines().keySet()) {
			BillingGroup bg = aliasesBillingGroups.get(groupName);
			if (bg == null) missingAliasCloud.add(groupName);
			else bg.getComputeWafs().addAll(cloudWafParser.getGroupNameWafLines().get(groupName));
		};
		if (missingAliasCloud.size()!=0) {
			IO.el("The following cloud WAF group names are missing from the MasterAccountsInfo sheet, add them and restart:");
			for (String gn: missingAliasCloud) IO.el("\t"+gn);
			System.exit(1);
		}
		
		//add hourly  WAF lines to each Billing Group
		ArrayList<String> missingAliasHourly = new ArrayList<String>();
		for (String groupName: hourlyWafParser.getGroupNameWafLines().keySet()) {
			BillingGroup bg = aliasesBillingGroups.get(groupName);
			
			if (bg == null) missingAliasHourly.add(groupName);
			else bg.getHourlyWafs().addAll(hourlyWafParser.getGroupNameWafLines().get(groupName));
		}
		if (missingAliasHourly.size()!=0) {
			IO.el("The following hourly WAF group names are missing from the MasterAccountsInfo sheet, add them and restart:");
			for (String gn: missingAliasHourly) IO.el("\t"+gn);
			System.exit(1);
		}
		
	}

	public static void main(String[] args) throws Exception {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CBiBilling2(args);
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
					case 'e': expenseFile = new File (args[i+1]); i++; break;
					case 'v': debug = true; break;
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
		date = Misc.getDateNoSpaces();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              CBI Billing2:   April 2024                          **\n" +
				"**************************************************************************************\n" +
				"Generates billing reports for the Cancer Bioinformatics Shared Resource. Only groups\n"+
				"with > $5 of hourly+ compute costs are billed.\n\n"+

				"Required Parameters:\n" +
				"-j Path to the exported cvs Jira 'Logged Time' report.\n"+
				"-m Path to the masterAccountInfo.xlsx spreadsheet updated from the prior month.\n"+
				"-a Path to the awsAccounts.xlsx spreadsheet.\n"+
				"-w Path to a dir containing the hourly and cloud 'WAF Tracking Schedule' xlsx files.\n" +
				"-c If available, path to a dir with the cloud AWS TDSynnex xlsx expense reports.\n"+
				"-e If available, path to a miscellaneous compute usage expense xlsx spreadsheet.\n"+
				"-o Path to write the Invoices.\n"+

				"\nExample: java -Xmx256M -jar pathTo/USeq/Apps/CBiBilling -j jiraTime.cvs -m \n"+
				"   masterAccountInfo.xlsx -w WAFs/ -c TDSynnex/ -o Invoices\n" +


		"**************************************************************************************\n");		
	}	
}
