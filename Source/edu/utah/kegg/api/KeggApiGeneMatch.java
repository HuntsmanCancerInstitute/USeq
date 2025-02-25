package edu.utah.kegg.api;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Pattern;

import org.biojavax.bio.seq.Position;

import util.gen.IO;
import util.gen.Misc;

public class KeggApiGeneMatch {
	
	//fields
	private String geneSymbol;
	private String keggGeneId = null;
	private String[] keggGeneNameAliases;
	private String keggDescription;
	private boolean ok = false;
	private String[] results = null;
	private String matchingResult = null;
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(geneSymbol);sb.append("\n");
		sb.append(ok);
		for (String r: results) {
			sb.append("\n");
			sb.append(r);
		}
		if (ok) {
			sb.append("\n");
			sb.append(keggGeneId);
			sb.append("\n");
			sb.append(keggDescription);
		}
		return sb.toString();
	}
	
	public String toStringSummary() {
		//GeneSymbol KeggGeneId MatchingKeggQueryResponse
		StringBuilder sb = new StringBuilder();
		sb.append(geneSymbol);sb.append("\t");
		sb.append(keggGeneId); sb.append("\t");
		sb.append(matchingResult); 
		return sb.toString();
	}
	
	public static final Pattern COMMA_SPACE = Pattern.compile(",\\s");
	public static final Pattern SEMI_COLON_SPACE = Pattern.compile(";\\s");
	
	public KeggApiGeneMatch(File apiResults, String geneSymbol) {
		this.geneSymbol = geneSymbol;
		
		results = IO.loadFileIntoStringArray(apiResults);
		
		if (results.length == 1) process(results[0]);
		
		else if (results.length > 1) {
			//look at the first alias in each result
			for (int i=0; i< results.length; i++) {
				//does the line start with hsa: ?
				if (results[i].startsWith("hsa:") == false) continue;
		
				//split on tab then ;
				String[] idInfo = Misc.TAB.split(results[i]);
				if (idInfo.length !=2) continue;
									
				String[] aliasesDescription = SEMI_COLON_SPACE.split(idInfo[1]);
				keggGeneNameAliases = COMMA_SPACE.split(aliasesDescription[0]);
				if (keggGeneNameAliases[0].equals(geneSymbol)) {
					process (results[i]);
					if (ok) return;
				}
				
			}
			//none found, process all
			for (int i=0; i< results.length; i++) {
				process (results[i]);
				if (ok) return;
			}
		}
	}
	
	public void process (String l) {
		//hsa:25	ABL1, ABL, BCR-ABL, CHDSKM, JTK7, bcr/abl, c-ABL, c-ABL1, p150, v-abl; ABL proto-oncogene 1, non-receptor tyrosine kinase

		//does the line start with hsa: ?
		if (l.startsWith("hsa:") == false) return;
		
		//split on tab then ;
		String[] idInfo = Misc.TAB.split(l);

		if (idInfo.length !=2) return;
							
		String[] aliasesDescription = SEMI_COLON_SPACE.split(idInfo[1]);
		if (aliasesDescription.length !=2) return;
		
		//is the geneSymbol in the aliases?
		keggGeneNameAliases = COMMA_SPACE.split(aliasesDescription[0]);
		for (String s: keggGeneNameAliases) {
			if (s.equals(geneSymbol)) {
				keggGeneId = idInfo[0].substring(4);
				keggDescription = aliasesDescription[1];
				matchingResult = l;
				ok = true;
				return;
			}
		}
	}

	public String getGeneSymbol() {
		return geneSymbol;
	}
	public String getKeggGeneId() {
		return keggGeneId;
	}
	public String[] getKeggGeneNameAliases() {
		return keggGeneNameAliases;
	}
	public String getKeggDescription() {
		return keggDescription;
	}
	public boolean isOk() {
		return ok;
	}
	public String[] getResults() {
		return results;
	}
	public static Pattern getCommaSpace() {
		return COMMA_SPACE;
	}
	public static Pattern getSemiColonSpace() {
		return SEMI_COLON_SPACE;
	}

}
