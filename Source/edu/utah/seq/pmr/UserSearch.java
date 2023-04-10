package edu.utah.seq.pmr;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.hci.misc.Util;
import util.gen.IO;
import util.gen.Misc;

public class UserSearch {
	
	private boolean caseInsensitive = true;
	private boolean allTermsMustMatch = false;
	private boolean partialMatchesAllowed = true;
	private boolean printMatchedDatasetURIs = false;
	private boolean clearPriorResults = false;
	private String printDataTable = null;
	private String dataTableToSearch = null;
	private String[] searchTerms = null; 
	private boolean goodToSearch = false;
	private boolean verbose = false;

	public UserSearch(String input) throws IOException {
		processArgs(input);
	}
	
	/**This method will process each argument and assign new variables
	 * @throws IOException 
	*/
	public void processArgs(String input) throws IOException {
		//remove any ' or " from the input
		input = Util.QUOTE_SINGLE.matcher(input).replaceAll("");
		input = Util.QUOTE_DOUBLE.matcher(input).replaceAll("");
		//IO.pl("\tINPUT: "+input);
		String[] args = Util.WHITESPACE.split(input);
		
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': dataTableToSearch = args[++i]; break;
					case 's': searchTerms = Misc.COMMA.split(args[++i]); break;
					case 'c': caseInsensitive = false; break;
					case 'e': partialMatchesAllowed = false; break;
					case 'a': allTermsMustMatch = true; break;
					case 'n': clearPriorResults = true; break;
					case 'v': verbose = true; break;
					case 'f': printMatchedDatasetURIs = true; break;
					case 'p': printDataTable = args[++i]; break;
					default: Misc.printErrAndExit("\nProblem, unknown search option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this search parameter: -"+test+"\n");
				}
			}
		}
		
		if (dataTableToSearch != null && searchTerms != null) goodToSearch = true;
		

	}

	public boolean isCaseInsensitive() {
		return caseInsensitive;
	}

	public boolean isAllTermsMustMatch() {
		return allTermsMustMatch;
	}

	public boolean isPartialMatchesAllowed() {
		return partialMatchesAllowed;
	}

	public boolean isPrintMatchedDatasetURIs() {
		return printMatchedDatasetURIs;
	}

	public boolean isClearPriorResults() {
		return clearPriorResults;
	}

	public String getPrintDataTable() {
		return printDataTable;
	}

	public String getDataTableToSearch() {
		return dataTableToSearch;
	}

	public String[] getSearchTerms() {
		return searchTerms;
	}

	public boolean isGoodToSearch() {
		return goodToSearch;
	}

	public boolean isVerbose() {
		return verbose;
	}

}
