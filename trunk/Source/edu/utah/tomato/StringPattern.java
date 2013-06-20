package edu.utah.tomato;

import java.util.regex.Pattern;

public class StringPattern {

	private String name;
	private Pattern pattern;
	
	public StringPattern(String name,Pattern pattern) {
		this.name = name;
		this.pattern = pattern;
	}
	
	public String getName() {
		return name;
	}
	
	public Pattern getPattern() {
		return pattern;
	}

}
