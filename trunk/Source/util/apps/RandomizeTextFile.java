package util.apps;
import java.io.*;

import util.gen.*;
import util.bio.annotation.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**Randomized the lines of a text.*/
public class RandomizeTextFile {

	private File[] textFiles;
	private int numberLinesToPrint = 0;

	public RandomizeTextFile(String[] args){
		//process args
		processArgs(args);

		System.out.println("FileName\tLinesInFile\tLinesPrinted");
		for (int i=0; i< textFiles.length; i++){
			System.out.print(textFiles[i].getName());

			//load regions 
			String[] lines = IO.loadFile(textFiles[i]);
			System.out.print("\t"+lines.length);

			//randomize
			Misc.randomize(lines, System.currentTimeMillis());

			//print to file
			File results = new File(textFiles[i].getParent(), Misc.removeExtension(textFiles[i].getName())+"_Randomized.txt.gz");
			int num = numberLinesToPrint;
			if (num == 0 || num > lines.length) num = lines.length;

			try {
				Gzipper out = new Gzipper(results);
				for (int j=0; j< num; j++) out.println(lines[j]);
				out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
			System.out.println("\t"+num);

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
					case 'f': textFiles = IO.extractFiles(new File (args[++i])); break;
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
		if (textFiles == null || textFiles[0].exists()==false){
			Misc.printErrAndExit("\nPlease enter a text file to randomize.\n");
		}
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Randomize Text File: May 2013                         **\n" +
				"**************************************************************************************\n" +
				"Randomizes the lines of a text file(s).\n" +

				"\nOptions:\n"+
				"-f Full path to a text file or directory containing such to randomize. Gzip/zip OK.\n"+
				"-n Number of lines to print, defaults to all.\n"+

				"\nExample: java -Xmx4G -jar pathTo/Apps/RandomizeTextFile -n 24560 -f\n" +
				"       /TilingDesign/oligos.txt.gz\n\n"+

		"************************************************************************************\n");
	}
}
