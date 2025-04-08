package edu.utah.billing;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import util.gen.IO;
import util.gen.Misc;

/**Parses a generic AWS bill, duplicate accounts will be summed. # comments and empty lines skipped.

#Account	Expense
8525629597	91.06280995
8525629597	3.461092659
8971640797	0.005837915
...

 * */
public class AwsAccountBilling {

	private HashMap<String, Float> accountTotal = new HashMap<String, Float>();
	
	public AwsAccountBilling (File parsedBill) throws IOException {
		BufferedReader in = IO.fetchBufferedReader(parsedBill);
		String line = null;
		String[] fields = null;
		
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#")) continue;
			fields = Misc.TAB.split(line);
			if (fields.length!=2) throw new IOException("\nERROR: failed to find two colums in the AWS bill for -> "+line);
			Float cost = accountTotal.get(fields[0]);
			if (cost == null) accountTotal.put(fields[0], Float.parseFloat(fields[1]));
			else {
				Float total = cost + Float.parseFloat(fields[1]);
				accountTotal.put(fields[0], total);
			}
		}
		in.close();
	}

	public HashMap<String, Float> getAccountTotal() {
		return accountTotal;
	}

	public static void main (String[] args) throws IOException {
		File test = new File("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2025/3_CBI_Mar_2025_WithAWS/AWS/fourPointsAWSBillParsing.txt");
	
		AwsAccountBilling aab = new AwsAccountBilling(test);
		IO.pl(aab.getAccountTotal());
	}
	
}
