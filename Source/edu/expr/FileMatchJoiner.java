package edu.expr;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.IO;
import util.gen.Misc;

public class FileMatchJoiner {
	
	private File keyFile;
	private HashMap<String,String> keys;
	private File[] filesToParse;
	private int keyIndex = 0;
	private int indexFilesToParse = 0;
	private boolean ignoreDuplicateKeys = false;
	private boolean printOnlyMatches = false;
	private boolean skipDuplicateKeys = false;
	
	public FileMatchJoiner (String[] args){
		processArgs(args);
		//load key file
		if (skipDuplicateKeys) keys = IO.parseUniqueKeyFile(keyFile, keyIndex);
		else keys = IO.parseFile(keyFile, keyIndex, ignoreDuplicateKeys);
		if (keys == null) Misc.printExit("Problem parsing key? Duplicates?");
		System.out.println("Number unique keys "+keys.size());
		
		//for each file parse and print matching lines
		parseAndPrintMatches();
		System.out.println("\nDone!\n");
	}
	
	public void parseAndPrintMatches(){
		for (int i=0; i< filesToParse.length; i++){
			System.out.println("\tParsing "+filesToParse[i].getName());
			try {
				File file = new File (filesToParse[i].getParentFile(), Misc.removeExtension(filesToParse[i].getName())+"_Joined.xls");
				PrintWriter out = new PrintWriter ( new FileWriter (file));
				BufferedReader in = IO.fetchBufferedReader(filesToParse[i]);
				String[] bits;
				Pattern pat = Pattern.compile("\\s");
				String line;
				//for each line
				int numMatches = 0;
				int numNoMatches = 0;
				int numBadColumns = 0;
				while ((line = in.readLine())!=null){
					bits = pat.split(line);
					if (indexFilesToParse >= bits.length) {
						System.out.println("Not enought columns, can't match -> "+line);
						if (printOnlyMatches == false) out.println(line);
						numBadColumns++;
					}
					else if (keys.containsKey(bits[indexFilesToParse]) == false){
						System.out.println("No match for -> "+bits[indexFilesToParse]+" from line -> "+line);
						if (printOnlyMatches == false) out.println(line);
						numNoMatches++;
					}
					else {
						out.println(line+"\t"+keys.get(bits[indexFilesToParse]));
						numMatches++;
					}
				}
				System.out.println("\n"+numMatches+"\t# Matches");
				System.out.println("\n"+numNoMatches+"\t# No Matches");
				System.out.println("\n"+numBadColumns+"\t# Bad Column lines");
				in.close();
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
					case 'k': keyFile = new File(args[i+1]); i++; break;
					case 'i': ignoreDuplicateKeys = true; break;
					case 'j': skipDuplicateKeys = true; break;
					case 'p': printOnlyMatches = true; break;
					case 'a': keyIndex = Integer.parseInt(args[i+1]); i++; break;
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
		if (keyFile == null ) Misc.printExit("\nCannot find your key file?\n");
		if (filesToParse == null) Misc.printExit("\nCannot find your file(s) to parse??\n");
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            File Match Joiner:  July 2008                         **\n" +
				"**************************************************************************************\n" +
				"FMJ loads a file and a particular column containing unique entries, a key, and then\n" +
				"appends the key line to lines in the parsed file that match a particular column.\n" +
				"Usefull for appending say chromosome coordinates to snp ids data, etc.\n\n" +
				
				"-k Full path file text for a tab delimited txt file (key) containing unique entries.\n" +
				"-f Ditto but for the file to parse, can specify a directory too.\n" +
				"-i Collapse duplicate keys.\n"+
				"-j Skip duplicate keys.\n"+
				"-a Column index containing the unique IDs in key, defaults to 0.\n"+
				"-b Column index containing the unique IDs in parsers, defaults to 0.\n"+
				"-p Print only matches.\n\n"+
				
				"Example: java -jar pathTo/Apps/FileMatchJoiner -k /snpChromMap.txt -m /SNPData/\n" +
				"     --b 2 -p\n\n" +
				
		"**************************************************************************************\n");		
	}
	
	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new FileMatchJoiner(args);
	}
}
