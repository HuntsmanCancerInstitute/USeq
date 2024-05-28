package edu.utah.billing;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import util.gen.Misc;
import util.gen.Num;

/**Object to represent a particular lab or group we support, e.g. AC Tan Laboratory.  It is loaded with all of the info related to billing and generates a group specific invoice.
 * Add new billing info here. */
public class BillingGroup {
	
	//Master Account Info
	private boolean cancerCenterMember = true;
	private double totalHoursBilled = 0d;
	private LinkedHashSet<String> aliases = null;
	private String groupName = null;
	
	//Hourly WAF
	private ArrayList<String[]> hourlyWafs = new ArrayList<String[]>();
	private boolean missingHourlyWafs = false;
	
	//Compute WAF
	private ArrayList<String[]> computeWafs = new ArrayList<String[]>();
	private boolean missingComuteWafs = false;
	
	//Monthly+ AWS Expenses, TDSynnex; where doing these quarterly since some are < $10; thus might be null;
	private ArrayList<AwsAccountExpense> awsAccountExpenses = new ArrayList<AwsAccountExpense>();
	
	//Misc Compute Expenses, one off stuff. IPA license, AWS downloads
	private ArrayList<MiscExpense> miscExpenses = new ArrayList<MiscExpense>();
	
	//Jira Hourly Tracking
	private ArrayList<JiraTicketSummary> jiraTickets = new ArrayList<JiraTicketSummary>();
	
	//Expenses
	private float totalHoursToBill = 0;
	private float totalHourlyExpenses = 0f;
	private String totalHourlySummary = null;
	private float additionalHourlyExpenses = 0f;	//for HSC billing
	private String additionalHourlySummary = null;
	private float totalComputeExpenses = 0f;
	
	public BillingGroup (boolean isCancerMember, double totalHoursBilled, LinkedHashSet<String> aliases) {
		cancerCenterMember = isCancerMember;
		this.totalHoursBilled = totalHoursBilled;
		this.aliases = aliases;
		//take first alias as groupName
		groupName = aliases.iterator().next();
	}
	
	public void calculateTotalExpenses() {
		//compute?
		if (awsAccountExpenses.size()!=0) totalComputeExpenses+= AwsAccountExpense.fetchTotalExpense(awsAccountExpenses);
		if (miscExpenses.size()!=0) totalComputeExpenses+= MiscExpense.fetchTotalExpense(miscExpenses);
		
		//hourly?
		if (jiraTickets.size()!=0) {
			for (JiraTicketSummary jts: jiraTickets) totalHoursToBill+= Float.parseFloat(jts.getHoursString());
			
			//calculate the expense, either $70/hr or $100/hr
			float pricePerHour = CBiBilling2.firstTierPricePerHour;
			if (totalHoursBilled > CBiBilling2.maxHoursForFirstTier) pricePerHour = CBiBilling2.secondTierPricePerHour;
			totalHourlyExpenses = totalHoursToBill* pricePerHour;
			totalHourlySummary = "Hourly Total:\t$"+Num.formatNumber(totalHourlyExpenses, 2)+ "\t("+ totalHoursToBill+ " x $"+pricePerHour+"/hr)";
			
			totalHoursBilled += totalHoursToBill;
			
			//non cancer?
			if (cancerCenterMember==false) {
				additionalHourlyExpenses = CBiBilling2.hscCostSharing * totalHoursToBill;
				additionalHourlySummary = "\t(Hourly HSC Cost Sharing Total:\t$"+Num.formatNumber(additionalHourlyExpenses, 2)+ "\t("+ totalHoursToBill+ " x $"+CBiBilling2.hscCostSharing+"/hr))";
			}
		}
	}
	
	public float getTotalExpenses() {
		return totalHourlyExpenses+ additionalHourlyExpenses+ totalComputeExpenses;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (cancerCenterMember) sb.append("Cancer\t");
		else sb.append("Non-Cancer\t");
		
		sb.append(totalHoursBilled);
		
		for (String a: aliases) {
			sb.append("\t");
			sb.append(a);
		}
		
		return sb.toString();
	}
	
	public String generateInvoice(String date, String hourlyWafHeader, String cloudWafHeader, boolean includeAdditionalHourlyExpenses, boolean pretty) {
		ArrayList<String> txtOut = new ArrayList<String>();

		txtOut.add("HCI Cancer Bioinformatics Invoice - "+date+" - "+groupName);
		if (pretty)txtOut.add("");
		txtOut.add("Aliases:\t"+Misc.linkedSetToString(aliases, "; "));
		txtOut.add("YTD Hourly Usage:\t"+Num.formatNumber(totalHoursBilled, 1)+" hrs");
		if (pretty)txtOut.add("");

		//Hourly Expenses
		if (totalHourlyExpenses !=0) {
			txtOut.add("Hourly Billing:");
			if (hourlyWafs.size()==0) {
				txtOut.add("\tHourly WAF:\tNo Hourly WAF - contact PI");
				missingHourlyWafs = true;
			}
			else {
				txtOut.add("\tHourly WAF:\t"+hourlyWafHeader);
				for (String[] wl: hourlyWafs) txtOut.add("\tHourly WAF:\t"+Misc.stringArrayToString(wl, "\t"));
			}
			if (pretty)txtOut.add("");

			txtOut.add("\tHourly Tickets:\t"+JiraTicketSummary.getToStringHeader());
			for (JiraTicketSummary jt: jiraTickets) txtOut.add("\tHourly Tickets:\t"+jt.toString());
			if (pretty)txtOut.add("");

			txtOut.add("\t\t"+totalHourlySummary);


			//any HSC billing?
			if (additionalHourlyExpenses > 0 && includeAdditionalHourlyExpenses) {
				txtOut.add("\t"+additionalHourlySummary);
			}
			if (pretty)txtOut.add("");
		}
		
		
		//Compute Usage billing
		if (totalComputeExpenses > 0) {
//if (totalHourlyExpenses ==0 && pretty) txtOut.add("");
			
			txtOut.add("Compute Billing:");
			if (computeWafs.size()==0) {
				txtOut.add("\tCloud WAF:\tNo Cloud WAF - contact PI");
				missingComuteWafs = true;
			}
			else {
				txtOut.add("\tCloud WAF:\t"+cloudWafHeader);
				for (String[] wl: computeWafs) txtOut.add("\tCloud WAF:\t"+Misc.stringArrayToString(wl, "\t"));
			}
			if (pretty)txtOut.add("");
			
			//Any AWS compute expenses?
			for (AwsAccountExpense aae: awsAccountExpenses) {
				String ae = "\tCloud Amazon Web Services (AWS) Acc:\t"+aae.getAwsAccountNumber()+"\t$"+Num.formatNumber(aae.getTotalExpense(), 2);
				txtOut.add(ae);
			}
			//Any Misc compute expenses?
			if (pretty && awsAccountExpenses.size()!=0 && miscExpenses.size()!=0) txtOut.add("");
			for (MiscExpense aae: miscExpenses) {
				String me = "\tMisc Compute:\t$"+Num.formatNumber(aae.getCost(), 2)+"\t"+aae.getDescription();
				txtOut.add(me);
			}
			if (pretty)txtOut.add("");
			txtOut.add("\t\tCompute Total:\t$"+Num.formatNumber(totalComputeExpenses, 2));
			if (pretty)txtOut.add("");
		}
		
		txtOut.add("Total Billing:\t$"+Num.formatNumber(totalComputeExpenses+ totalHourlyExpenses, 2)+"\n");
		if (pretty) {
			txtOut.add("Questions?\n\tEmail: "+CBiBilling2.contactEmail);
			txtOut.add("\tOperating Policies and WAF forms: "+CBiBilling2.cbiPolicyUrl);
			txtOut.add("");
		}
		
		return Misc.stringArrayListToString(txtOut, "\n");
	}
	
	public String generateHscBillingLine(int wafIndex, Integer[] indexesToPull) {
		String[] fields = hourlyWafs.get(wafIndex);
		StringBuilder sb = new StringBuilder(fields[indexesToPull[0]]);
		for (int i=1; i< indexesToPull.length; i++) {
			sb.append("\t");
			sb.append(fields[indexesToPull[i]]);
		}
		return sb.toString();
	}


	public boolean isCancerCenterMember() {
		return cancerCenterMember;
	}

	public double getTotalHoursBilled() {
		return totalHoursBilled;
	}

	public LinkedHashSet<String> getAliases() {
		return aliases;
	}

	public ArrayList<String[]> getHourlyWafs() {
		return hourlyWafs;
	}

	public ArrayList<String[]> getComputeWafs() {
		return computeWafs;
	}

	public double getTotalHourlyExpenses() {
		return totalHourlyExpenses;
	}

	public double getTotalComputeExpenses() {
		return totalComputeExpenses;
	}

	public ArrayList<AwsAccountExpense> getAwsAccountExpenses() {
		return awsAccountExpenses;
	}

	public ArrayList<MiscExpense> getMiscExpenses() {
		return miscExpenses;
	}

	public ArrayList<JiraTicketSummary> getJiraTickets() {
		return jiraTickets;
	}

	public String getTotalHourlySummary() {
		return totalHourlySummary;
	}

	public float getAdditionalHourlyExpenses() {
		return additionalHourlyExpenses;
	}

	public String getAdditionalHourlySummary() {
		return additionalHourlySummary;
	}

	public boolean isMissingHourlyWafs() {
		return missingHourlyWafs;
	}

	public boolean isMissingComuteWafs() {
		return missingComuteWafs;
	}

	public float getTotalHoursToBill() {
		return totalHoursToBill;
	}

	public String getGroupName() {
		return groupName;
	}
}
