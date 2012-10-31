package util.bio.parsers;
import java.io.*;

import util.apps.ScoreChromosomes;
import util.gen.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SynonymMatching {
	
	//fields
	private File standardGeneNameFile;
	//list of synonyms from Agilent data, all syns on same line, TAB delimited
	private File totalAgilentMixedGeneNameFile;
	//pick list of Agilent gene names user want's to convert to standard
	private File agilentPickGeneList;
	private String[] agilentPicksToConvert;
	//hashes
	private HashMap <String,String> standardVsRefSeqNames;
	private HashMap <String,String> refSeqVsStandardNames;
	private HashMap <String,String> agilentToStandardGeneNames;
	
	
	//constructors
	public SynonymMatching (String[] args){
		//process args
		processArgs(args);
		
		//load standardVsRefSeqNames
		standardVsRefSeqNames = IO.loadFileIntoHashMap(standardGeneNameFile);
		refSeqVsStandardNames = Misc.invert(standardVsRefSeqNames);
		
		//load agilentToStandardGeneNames HashMap
		loadAgilentToStandardGeneNames();
		
		//load agilent picks
		agilentPicksToConvert = IO.loadFile(agilentPickGeneList);
		
		//match
		matchPicks();
	}
	
	public void matchPicks(){
		System.out.println("\nUserPick\tMatchingStandardName");
		for (int i=0; i< agilentPicksToConvert.length; i++){
			//any match
			if (agilentToStandardGeneNames.containsKey(agilentPicksToConvert[i])){
				System.out.println(agilentPicksToConvert[i] +"\t"+ agilentToStandardGeneNames.get(agilentPicksToConvert[i]));
			}
			else System.out.println(agilentPicksToConvert[i] +"\t");
		}
	}
	

	public void loadAgilentToStandardGeneNames(){
		agilentToStandardGeneNames = new HashMap();
		
		try {
			BufferedReader in = new BufferedReader (new FileReader (totalAgilentMixedGeneNameFile));
			String line;
			String[] syms;
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0 || line.startsWith("#")) continue;
				//split line on tabs
				syms = line.split("\\t+");
				//for each synonym look for match in standardNames
				String standardGene = null;
				for (int i=0; i< syms.length; i++){
					syms[i] = syms[i].trim();
					if (standardVsRefSeqNames.containsKey(syms[i])){
						standardGene = syms[i];
						break;
					}
				}
				//was there a match? if not then look for in refSeqNames
				if (standardGene == null){
					for (int i=0; i< syms.length; i++){
						syms[i] = syms[i].trim();
						if (refSeqVsStandardNames.containsKey(syms[i])){
							standardGene = refSeqVsStandardNames.get(syms[i]);
							break;
						}
					}
				}
				//were any matches made to standard genes?
				if (standardGene != null) for (int i=0; i< syms.length; i++) agilentToStandardGeneNames.put(syms[i], standardGene);
		
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': standardGeneNameFile = new File(args[i+1]); i++; break;
					case 'a': totalAgilentMixedGeneNameFile = new File(args[i+1]); i++; break;
					case 'p': agilentPickGeneList = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
	}
	
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SynonymMatching(args);
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Synonym Matching: Jan 2008                           **\n" +
				"**************************************************************************************\n" +
				"SM attempts to assign a standard text to each pick using synonym tables. For each text\n" +
				"in the picks file, it is used to fetch the associated synonyms and then SM attempts to\n" +
				"match it to a standard text or the alternative text.  If a match is found, the original\n" +
				"pick and it's associated standard text are printed to screen.\n\n" +
				
				"-s The full path file text for a two column tab delimited text file containing the\n" +
				"      standard text and an alternative text.\n" +
				"-a The full path file text for a multi column tab delimited text file containing\n" +
				"      synonyms. All names on a line are considered synonyms.\n" +
				"-p The full path file text for a one column tab delimited text file containing select\n" +
				"      names from the synonyms file.\n" +
				"\n" +
				"Example: java -Xmx500M -jar pathTo/T2/Apps/SynonymMatching -s /Anno/zv7Stnd_GbNames.txt\n" +
				"      -a /Anno/zv7AgilentSynonyms.txt -p /Data/upRegPicksAgilent.txt\n\n" +

		"**************************************************************************************\n");
	}

	

}
