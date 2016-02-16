package util.apps;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
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
	private String[] orderedFileNames;
	private File output;
	
	
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
					case 'o': orderedFileNames = args[++i].split(","); break;
					case 'c': output = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
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
				"**                             File Joiner: Feb 2016                                **\n" +
				"**************************************************************************************\n" +
				"Joins text files into a single file, avoiding line concatenations. This is a problem\n" +
				"with using 'cat * >> combine.txt'. Removes empty lines. Option to follow custom order.\n\n"+
				
				"Parameters:\n"+
				"-f Full path text for the directory containing the text files.\n" +
				"-o (Optional) Order the files using this comma delimited list, no spaces. Not all\n"+
				"         need to exist.\n"+
				"-c (Optional) Concatinated results file.\n"+
	
				"\n" +
				"Example: java -jar pathTo/T2/Apps/FileJoiner -f /affy/SplitFiles/\n" +
				"    -o 1.fasta,2.fasta,3.fasta,4.fasta\n" +
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
		
		//fetch the files to concat
		File[] files = IO.extractFiles(directory);
		
		//possibly reorder 
		files = reorder(files);

		//make output file if needed
		if (output == null) output = new File(directory,"combineFile.txt.gz");
		
		int counter = 0;
		String line;
		try{
			Gzipper out = new Gzipper(output);
			for (int i=0; i<files.length; i++){
				BufferedReader in = IO.fetchBufferedReader(files[i]);
				System.out.println("\tWriting "+files[i]+"...");
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
		System.out.println(counter+" Lines written to "+output+"\n");
	}
	private File[] reorder(File[] files) {
		//anything to order?
		if (orderedFileNames == null) return files;
		
		//add order to a hash
		HashMap<String, Integer> orderedNames = new HashMap<String, Integer>();
		for (int i=0; i< orderedFileNames.length; i++) orderedNames.put(orderedFileNames[i], new Integer(i));
		
		//split into ordered and not found
		ArrayList<File> notFoundFiles = new ArrayList<File>();
		File[] orderedSubSet = new File[orderedFileNames.length];
		
		//for each incoming file
		for (int i=0; i< files.length; i++){
			String testName = files[i].getName();
			if (orderedNames.containsKey(testName) == false) notFoundFiles.add(files[i]);
			else {
				int index = orderedNames.get(testName);
				orderedSubSet[index] = files[i];
			}
		}
		
		//combine
		File[] ordered = new File[files.length];
		int index = 0;
		for (int i=0; i< orderedSubSet.length; i++){
			if (orderedSubSet[i] !=null) ordered[index++] = orderedSubSet[i];
		}
		int num = notFoundFiles.size();
		for (int i=0; i< num; i++){
			ordered[index++] = notFoundFiles.get(i);
		}
		return ordered;
	}
}

