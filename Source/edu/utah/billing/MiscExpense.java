package edu.utah.billing;

import java.util.ArrayList;

public class MiscExpense {
	private double cost = 0;
	private String description = null;

	public MiscExpense (double cost, String description){
		this.cost = cost;
		this.description = description;
	}
	
	public static float fetchTotalExpense(ArrayList<MiscExpense> accounts) {
		float total = 0f;
		for (MiscExpense aae: accounts) total+= aae.getCost();
		return total;
	}

	public double getCost() {
		return cost;
	}

	public String getDescription() {
		return description;
	}
	
}