package util.apps;
import java.io.*;
import util.gen.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;

/**
 
 * @author nix
 *
 */
public class FileColumnSplitter {
	//fields
	private File file;
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		System.out.println("Launching FileSplitter...");
		new FileColumnSplitter(args);
	}

	public FileColumnSplitter(String[] args){
		processArgs(args);

		try{
			BufferedReader in = IO.fetchBufferedReader(file);
			//read first line
			String line = in.readLine();
			//split it and make print writers
			Pattern tab = Pattern.compile("\\t");
			String[] cells = tab.split(line);
			int numberColumns = cells.length;
			PrintWriter[] outs = new PrintWriter[numberColumns];
			for (int i=0; i< numberColumns; i++) outs[i] = new PrintWriter( new File(file.getParentFile(), cells[i]));

			//for each line
			while ((line = in.readLine())!= null){
				cells = tab.split(line);
				if (cells.length != numberColumns) Misc.printErrAndExit("\nThe following line doesn't contain the correct number of columns:\n"+line);
				for (int i=0; i< numberColumns; i++) outs[i].println(cells[i]);
			}

			//close writers
			for (int i=0; i< outs.length; i++) outs[i].close();
			in.close();
		}
		catch (IOException e){
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
					case 'f': file = new File(args[++i]); break;
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
				"**                          File Column Splitter: July 2011                         **\n" +
				"**************************************************************************************\n" +
				"Splits a text file by column to using the header row as the file names.\n\n"+

				"Required Parameters:\n"+
				"-f Full path file text to split by column, first row should contain file names.\n"+

				"\n" +
				"Example: java -Xmx256M -jar pathTo/T2/FileSplitter -f /affy/geneIntensities.txt\n" +
				"\n" +
		"**************************************************************************************\n");
	}	


}

