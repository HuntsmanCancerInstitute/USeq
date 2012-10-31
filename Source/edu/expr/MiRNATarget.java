package edu.expr;

import java.util.*;

public class MiRNATarget {

	//fields
	private String name;
	private String alias;
	private double score;
	private double pvalue;
	
	//constructor
	public MiRNATarget (String name, String alias, String score, String pvalue){
		this.name = name;
		this.alias = alias;
		this.score = Double.parseDouble(score);
		this.pvalue = Double.parseDouble(pvalue);
	}
	
	public String toString(){
		return name+"\t"+alias+"\t"+score+"\t"+pvalue;
	}

	public String getAlias() {
		return alias;
	}

	public String getName() {
		return name;
	}

	public double getPvalue() {
		return pvalue;
	}

	public double getScore() {
		return score;
	}
}
