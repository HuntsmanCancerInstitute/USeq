package edu.expr;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class FileCrossFilter {

	private File matcherFile;
	private String[] uniqueIDs;
	private File[] filesToParse;
	private int[] indexMatcher;
	private int[] indexFilesToParse;
	private HashMap<String,String> toParseMap = new HashMap<String, String>();
	private PrintWriter out = null;
	

	public FileCrossFilter (String[] args){
		processArgs(args);
		//load matcher file
		uniqueIDs = IO.parseColumns(matcherFile, indexMatcher);
		//for each file parse and print matching lines
		parseAndPrintMatches();
		System.out.println("\nDone!\n");
	}

	public void parseAndPrintMatches() {
		try {
			System.out.println("Parsing...\nFileName\tNumMatches\tNumFailedMatches");
			for (int i=0; i< filesToParse.length; i++){
				System.out.print(filesToParse[i].getName());
				
				//make writer for the parsed file
				File file = new File (filesToParse[i].getParentFile(), Misc.removeExtension(filesToParse[i].getName())+"_Matched.xls");
				out = new PrintWriter ( new FileWriter (file));
				
				//load the file into a hashmap where the key is the columns they defined, the value the line of the file
				loadHashMap(filesToParse[i]);
				
				//for each matcher key
				int numMatches = 0;
				int numFails = 0;
				for (int j=0; j< uniqueIDs.length; j++){
					if (toParseMap.containsKey(uniqueIDs[j])) {
						out.println(toParseMap.get(uniqueIDs[j]));
						numMatches++;
					}
					else {
						out.println();
						numFails++;
					}
				}
				out.close();
				System.out.println("\t"+numMatches+"\t"+numFails);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	private void loadHashMap(File toParse) throws IOException {
		toParseMap.clear();
		BufferedReader in = IO.fetchBufferedReader(toParse);
		String line;
		String[] columns;
		while ((line = in.readLine())!= null){
			if (line.startsWith("#")) {
				out.println(line);
				continue;
			}
			columns = Misc.TAB.split(line);
			StringBuilder sb = new StringBuilder(columns[indexFilesToParse[0]]);
			for (int i=1; i< indexFilesToParse.length; i++){
				sb.append("\t");
				sb.append(columns[indexFilesToParse[i]]);
			}
			//does it exist?
			String key = sb.toString();			
			if (toParseMap.containsKey(key)) Misc.printErrAndExit("\nERROR: Duplicate key found in file to parse, see '"+key+"' in "+toParse);
			toParseMap.put(key, line);
		}
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': filesToParse = IO.extractFiles(new File(args[++i])); break;
					case 'm': matcherFile = new File(args[++i]); break;
					case 'a': indexMatcher = Num.parseInts(args[++i], Misc.COMMA); break;
					case 'b': indexFilesToParse = Num.parseInts(args[++i], Misc.COMMA); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					e.printStackTrace();
				}
			}
		}
		//check agilentFiles
		if (matcherFile == null ) Misc.printExit("\nCannot find your matcher file?\n");
		if (filesToParse == null) Misc.printExit("\nCannot find your file(s) to parse??\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            File Cross Filter: Sept 2017                         **\n" +
				"**************************************************************************************\n" +
				"FCF takes one or more columns in the matcher file and uses these as a key to parse and\n"+
				"save matching keys in the to parse files. Use this to parse lines in files that match\n"+
				"those in another. Keys must be unique. The order and number of the rows in the matcher\n"+
				"file is preserved, if a match is not found in the parsed file, a blank line is inserted\n"+
				"instead.\n\n" +

				"-m Path to a tab delimited txt (.gz/.zip OK) file to use in matching.\n" +
				"-a One or more column indexs in the matcher file to use as the key.\n"+
				"-p Path to a file or directory of files to parse (.gz/.zip OK).\n" +
				"-b One or more column indexes in the parse file(s) to use as a key.\n"+

				"\nExample: java -jar pathTo/USeq/Apps/FileCrossFilter -m intRegions.bed -a 0,1,2\n" +
				"     -p SpreadSheetData/ -b 0,1,2 \n\n" +

				"**************************************************************************************\n");		
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new FileCrossFilter(args);
	}
}
