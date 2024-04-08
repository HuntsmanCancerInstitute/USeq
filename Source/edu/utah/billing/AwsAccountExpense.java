package edu.utah.billing;

import java.util.ArrayList;

public class AwsAccountExpense {

	private String awsAccountNumber = null;
	private float totalExpense = 0f;
	
	public AwsAccountExpense (String number, float total) {
		awsAccountNumber = number;
		totalExpense = total;
	}

	public String getAwsAccountNumber() {
		return awsAccountNumber;
	}

	public float getTotalExpense() {
		return totalExpense;
	}
	
	public static float fetchTotalExpense(ArrayList<AwsAccountExpense> accounts) {
		float total = 0f;
		for (AwsAccountExpense aae: accounts) total+= aae.getTotalExpense();
		return total;
	}
}
