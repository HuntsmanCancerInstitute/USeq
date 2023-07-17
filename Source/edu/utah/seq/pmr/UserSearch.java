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
	private boolean addMatches = false;
	private boolean printDatasetInfo = false;
	private boolean printDatasetNames = false;
	private boolean verbose = false;

	public UserSearch(String input) throws IOException {
		processArgs(input);
	}
	
	public static void printInteractiveMenu() {
		StringBuilder sb = new StringBuilder("\nInteractive Search Options:\n");
		sb.append("   -d Data table to search, choose from: Diagnosis, PhysicianDiseaseGroups, SpecimenSites, SpecimenIds, Sex, DatasetIds\n");
		sb.append("   -s Search terms, comma delimited, no spaces\n");
		sb.append("   -c Terms are case sensitive.\n");
		sb.append("   -e Terms are exact, no partial matching.\n");
		sb.append("   -m All terms must match.\n");
		sb.append("   -a Add matches to the prior result set.\n");
		sb.append("   -x New search, clear any prior results.\n");
		sb.append("   -p Print the contents of the named data table, see -d\n");
		sb.append("   -f Print the available file AWS URIs for the matched datasets.\n");
		sb.append("   -i Print the clinical and test details for the matched datasets.\n");
		sb.append("   -n Print the matched dataset names.\n");
		sb.append("   -v Verbose output.\n");
		sb.append("      e.g. '-d Diagnosis -s Liver,Metastasis -m' then '-d Sex -s F' then '-i -f'\n");
		IO.p(sb.toString());
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
					case 'm': allTermsMustMatch = true; break;
					case 'a': addMatches = true; break;
					case 'x': clearPriorResults = true; break;
					case 'v': verbose = true; break;
					case 'f': printMatchedDatasetURIs = true; break;
					case 'i': printDatasetInfo = true; break;
					case 'n': printDatasetNames = true; break;
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

	public boolean isAddMatches() {
		return addMatches;
	}

	public boolean isPrintDatasetInfo() {
		return printDatasetInfo;
	}

	public boolean isPrintDatasetNames() {
		return printDatasetNames;
	}

}
