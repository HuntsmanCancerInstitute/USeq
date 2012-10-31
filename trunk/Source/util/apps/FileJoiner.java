package util.apps;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**
 * Joins text files into a single file, avoiding line concatenations. 
 * This is a problem with using 'cat * > combine.txt'.  
 * Removes empty lines.
 */
public class FileJoiner {
	//fields
	private File directory;
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		String[] celDirs = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': directory = new File(args[i+1]); i++; break;
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
		//check to see if they entered required params
		if (directory==null || directory.isDirectory() == false){
			System.out.println("\nCannot find your directory!\n");
			System.exit(0);
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             File Joiner: Feb 2005                                **\n" +
				"**************************************************************************************\n" +
				"Joins text files into a single file, avoiding line concatenations. This is a problem\n" +
				"with using 'cat * >> combine.txt'.  Removes empty lines.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path text for the directory containing the text files.\n" +
	
				"\n" +
				"Example: java -jar pathTo/T2/Apps/FileJoiner -f /affy/SplitFiles/\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		System.out.println("\nLaunching FileJoiner...");
		new FileJoiner(args);
	}
	
	public FileJoiner(String[] args){
		processArgs(args);
		File[] files = IO.extractFiles(directory);
		File combineFile = new File(directory,"combineFile.txt");
		
		int counter = 0;
		String line;
		try{
			PrintWriter out = new PrintWriter(new FileWriter(combineFile));
			for (int i=0; i<files.length; i++){
				BufferedReader in = IO.fetchBufferedReader(files[i]);
				System.out.println("\tReading "+files[i]+"...");
				while ((line = in.readLine()) !=null) { 
					line = line.trim();
					//only print non empty lines
					if (line.length()!=0){
						out.println(line);
						counter++;
					}
				}
				in.close();
			}
			out.close();
		}
		catch (IOException e){
			e.printStackTrace();
		}
		System.out.println(counter+" Lines written to "+combineFile+"\n");
	}
}

