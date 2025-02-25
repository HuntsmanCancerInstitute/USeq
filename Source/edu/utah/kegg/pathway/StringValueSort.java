package edu.utah.kegg.pathway;

public class StringValueSort implements Comparable<StringValueSort> {
	
	private StringBuilder cargo = null;
	private double sortValue = 0;
	
	public StringValueSort (StringBuilder cargo, double sortValue){
		this.cargo = cargo;
		this.sortValue = sortValue;
	}

	/**Sorts smallest to largest*/
	public int compareTo(StringValueSort o) {
		if (this.sortValue< o.sortValue) return -1;
		if (this.sortValue> o.sortValue) return 1;
		return 0;
	}

	public StringBuilder getCargo() {
		return cargo;
	}

	public double getSortValue() {
		return sortValue;
	}

}
