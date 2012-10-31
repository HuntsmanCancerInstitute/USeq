package util.apps;
import java.io.*;

import util.gen.*;
import util.bio.annotation.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**Randomized the lines of a text.*/
public class RandomizeTextFile {
	
	private File textFile;
	private int numberLinesToPrint = 0;

	public RandomizeTextFile(String[] args){
		//process args
		processArgs(args);
		
		//load regions 
		System.out.println("Loading "+textFile);
		String[] lines = IO.loadFile(textFile);
		
		//randomize
		System.out.println("Randomizing...");
		Misc.randomize(lines, System.currentTimeMillis());
		
		//print to file
		File results = new File(textFile.getParent(), Misc.removeExtension(textFile.getName())+"_Randomized.txt");
		if (numberLinesToPrint == 0 || numberLinesToPrint > lines.length) numberLinesToPrint = lines.length;
		System.out.println("Saving "+numberLinesToPrint +" lines to "+results);
		try {
			PrintWriter out = new PrintWriter ( new FileWriter( results));
			for (int i=0; i< numberLinesToPrint; i++) out.println(lines[i]);
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		System.out.println("\nDone!\n");
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new RandomizeTextFile(args);
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
					case 'f': textFile = new File (args[++i]); break;
					case 'n': numberLinesToPrint = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (textFile == null || textFile.exists()==false){
			Misc.printErrAndExit("\nPlease enter a text file to randomize.\n");
		}
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Randomize Text File: May 2009                         **\n" +
				"**************************************************************************************\n" +
				"Randomizes the lines of a text file.\n" +

				"\nOptions:\n"+
				"-f Text file to randomize.\n"+
				"-n Number of lines to print, defaults to all.\n"+

				"\nExample: java -Xmx4000M -jar pathTo/Apps/RandomizeTextFile -n 24560 -t\n" +
				"       /TilingDesign/oligos.txt\n\n"+

		"************************************************************************************\n");
	}
}
