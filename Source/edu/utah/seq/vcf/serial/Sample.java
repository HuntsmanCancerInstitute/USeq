package edu.utah.seq.vcf.serial;

import java.util.HashMap;

public class Sample {
	
	private Integer number = null;
	private String name = null;
	
	private HashMap<String, Variant> variants = new HashMap<String, Variant>();

	public Sample(Integer number, String name) {
		this.number = number;
		this.name = name;
	}

	public Integer getNumber() {
		return number;
	}

	public String getName() {
		return name;
	}

	public HashMap<String, Variant> getVariants() {
		return variants;
	}

}
