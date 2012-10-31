package util.apps;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;



/**
 * Removes empty lines from text files in a given directory.
 * */
public class FileTrimmer {
	//fields
	private File directory;
	private boolean stripDuplicates = false;
	
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
					case 'f': directory = new File(args[i+1]); i++; break;
					case 's': stripDuplicates = true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     File Trimmer: March 2005                                      **\n" +
				"**************************************************************************************\n" +
				"Removes any empty lines from any text files in a given directory. Can also strip\n" +
				"duplicate lines.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path text for the directory containing the text files.\n" +
				"-s Strip duplicate lines.\n"+
	
				"\n" +
				"Example: java FileJoiner -f /affy/SplitFiles/ -s\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		System.out.println("\nLaunching FileJoiner...");
		new FileTrimmer(args);
	}
	
	public FileTrimmer(String[] args){
		processArgs(args);
		File[] files = IO.extractFiles(directory);
		String line;
		try{
			for (int i=0; i<files.length; i++){
				StringBuffer sb = new StringBuffer();
				HashSet hash = new HashSet();
				BufferedReader in = new BufferedReader(new FileReader(files[i]));
				System.out.println("\tProcessing "+files[i]+"...");
				while ((line = in.readLine()) !=null) { 
					line = line.trim();
					//only print non empty lines
					if (line.length()!=0){
						if (hash.contains(line)) continue;
						else {
							sb.append(line+"\n");
							hash.add(line);
						}
					}
				}
				in.close();
				//overwrite original file
				IO.writeString(sb.toString(), new File(files[i]+".trm"));
			}
		}
		catch (IOException e){
			e.printStackTrace();
		}
	}
}

