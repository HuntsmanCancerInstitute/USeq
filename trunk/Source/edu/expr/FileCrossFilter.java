package edu.expr;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.IO;
import util.gen.Misc;

public class FileCrossFilter {
	
	private File matcherFile;
	private String[] uniqueIDs;
	private File[] filesToParse;
	private int indexMatcher = 0;
	private int indexFilesToParse = 0;
	private boolean ignoreDuplicateKeys = false;
	
	public FileCrossFilter (String[] args){
		processArgs(args);
		//load matcher file
		uniqueIDs = IO.parseColumn(matcherFile, indexMatcher);
		//for each file parse and print matching lines
		parseAndPrintMatches();
		System.out.println("\nDone!\n");
	}
	
	public void parseAndPrintMatches(){
		for (int i=0; i< filesToParse.length; i++){
			System.out.println("\tParsing "+filesToParse[i].getName());
			HashMap map = IO.parseFile(filesToParse[i], indexFilesToParse,ignoreDuplicateKeys );
			if (map == null) Misc.printExit("Error: duplicate key found in "+filesToParse[i].getName() );
			try {
				File file = new File (filesToParse[i].getParentFile(), Misc.removeExtension(filesToParse[i].getName())+"_Matched.xls");
				PrintWriter out = new PrintWriter ( new FileWriter (file));
				//for each uniqueID
				for (int j=0; j< uniqueIDs.length; j++){
					if (map.containsKey(uniqueIDs[j])) out.println(map.get(uniqueIDs[j]));
					else {
						out.println();
						System.err.println("\tKey not found! "+uniqueIDs[j]+" in "+filesToParse[i].getName());
					}
				}
				out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': filesToParse = IO.extractFiles(new File(args[i+1]));  i++; break;
					case 'm': matcherFile = new File(args[i+1]); i++; break;
					case 'i': ignoreDuplicateKeys = true; break;
					case 'a': indexMatcher = Integer.parseInt(args[i+1]); i++; break;
					case 'b': indexFilesToParse = Integer.parseInt(args[i+1]); i++; break;
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
				"**                            File Cross Filter: March 2008                         **\n" +
				"**************************************************************************************\n" +
				"FCF take a column in the matcher file and uses it to parse the rows from other files.\n" +
				"Useful for pulling out and printing in order the rows that match the first file.\n\n" +
				
				"-m Full path file text for a tab delimited txt file to use in matching.\n" +
				"-f Full path file text to parse, can specify a directory too.\n" +
				"-i Ignore duplicate keys.\n"+
				"-a Column index containing the unique IDs in the matcher, defaults to 0.\n"+
				"-b Column index containing the unique IDs in the parsers, defaults to 0.\n"+
				
				"Example: java -jar pathTo/T2/Apps/FileCrossFilter -f /extendedArrayData/ -m /old/\n" +
				"     originalArray.txt -a 2 -b 2\n\n" +
				
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
