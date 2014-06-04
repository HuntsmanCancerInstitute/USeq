package edu.utah.tomato.model;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TFMatchObject {

	private String name;
	private Pattern pattern;
	private String matchType;
	private Matcher match;
	
	public TFMatchObject(String matchName,Pattern matchingPattern, String prefixType) {
		this.name = matchName;
		this.pattern = matchingPattern;
		
		this.matchType = prefixType;
		this.match = null;
	}
	
	public String getMatchName() {
		return name;
	}
	
	public Pattern getMatchingPattern() {
		return pattern;
	}
	
	public String getPrefixType() {
		return matchType;
	}
	
	public boolean testMatch(String potential) {
		Matcher m = this.pattern.matcher(potential);
		if (m.matches()) {
			this.match = m;
			return true;
		} else {
			return false;
		}
	}
	
	public String getMatchPrefix() {
		return this.match.group(1);
	}
	
	public String getMatchFull() {
		return this.match.group(0);
	}
	
	//public String getSuffix() {
	//	return this.match.group(2);
	//}

}
