package util.apps;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.*;

/**Prints out select columns from a tab delimited file, skips lines, appends words, etc. */
public class PrintSelectColumns {
	private int numberToSkip = 0;
	private String columnWord = "";
	private int[] columnIndexes;
	private int maxIndex;
	private File[] filesToParse = null;
	private int printLastNumberOfLines = 0;
	private boolean appendFileName = false;
	private boolean appendRowNumberColumn = false;
	private boolean skipBlanksAndOddNum = false;
	private boolean printAvailableColumns = false;


	public PrintSelectColumns (String[] args) {
		//process args
		processArgs(args);
		//calc max column index
		maxIndex = columnIndexes[Num.findMaxIntIndex(columnIndexes)];
		String line = null;
		try {
			//for each file
			for (int z =0; z<filesToParse.length; z++){
				if (filesToParse[z].isDirectory()) {
					System.out.println("\tSkipping directory file "+filesToParse[z].getName());
					continue;
				}
				File parsedFile = new File (Misc.removeExtension(filesToParse[z].getCanonicalPath())+".PSC.xls");
				BufferedReader in = IO.fetchBufferedReader(filesToParse[z]);
				PrintWriter out = new PrintWriter (new FileWriter(parsedFile));
				System.out.println("Parsing "+filesToParse[z].getName());

				//set file text as columnWord?
				if (appendFileName){
					columnWord = filesToParse[z].getName()+"\t";
				}
				//skip header lines?
				if (printLastNumberOfLines !=0) {
					int numLines = (int)IO.countNumberOfLines(filesToParse[z]);
					numberToSkip = numLines-printLastNumberOfLines;
				}
				for (int i=0; i< numberToSkip; i++) in.readLine();
				int rowNumber = 0;
				while ((line = in.readLine()) != null){
					line = line.trim();
					//is it empty?
					if (line.length() ==0) {
						if (skipBlanksAndOddNum) {
							System.out.println("\tSkipping empty line # "+rowNumber);
							continue;
						}
						else {
							out.println();
							System.out.println("\tPrinting empty line # "+rowNumber);
						}
					}

					//correct number of columns?
					String[] columns = line.split("\t");
					StringBuffer sb = new StringBuffer (columnWord);
					if (appendRowNumberColumn)  {
						sb.append(rowNumber);
						sb.append("\t");
					}
					if (printAvailableColumns){
						for (int i=0; i< columnIndexes.length; i++){
							if ( columnIndexes[i] < columns.length){
								sb.append(columns[columnIndexes[i]]);
								if (i != (columns.length -1))sb.append("\t");
							}
						}
					}
					else {
						if (columns.length <= maxIndex ){
							if (skipBlanksAndOddNum) {
								System.out.println("\tSkipping odd num column line # "+rowNumber+" line ->"+line);
								continue;
							}
							else Misc.printExit("\nError: line found with odd number of columns in "+filesToParse[z].getName()+" line number "+rowNumber+"\n");
						}
						sb.append(columns[columnIndexes[0]]);
						for (int i=1; i< columnIndexes.length; i++){
							sb.append("\t");
							sb.append(columns[columnIndexes[i]]);
						}
					}
					out.println(sb);
					rowNumber++;
				}
				in.close();
				out.close();
			}
		}
		catch (Exception e){
			System.out.println("Line -> "+line);
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
					case 'f': filesToParse = IO.extractFiles(new File (args[i+1])); i++; break;
					case 'i': columnIndexes = Num.parseInts(args[i+1].split(",")); i++; break;
					case 'n': numberToSkip = Integer.parseInt(args[i+1]); i++; break;
					case 'l': printLastNumberOfLines = Integer.parseInt(args[i+1]); i++; break;
					case 'c': columnWord = args[i+1]+"\t"; i++; break;
					case 'd': appendFileName = true; break;
					case 'r': appendRowNumberColumn = true; break;
					case 's': skipBlanksAndOddNum = true; break;
					case 'a': printAvailableColumns = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}

	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Print Select Columns: Sept 2010                        **\n" +
				"**************************************************************************************\n" +
				"Spread sheet manipulation.\n\n"+

				"Required Parameters:\n"+
				"-f Full path file or directory text for tab delimited text file(s)\n" +
				"-i Column indexs to print, comma delimited, no spaces\n"+
				"-n Number of initial lines to skip\n"+
				"-l Print only this last number of lines\n"+
				"-c Column word to append onto the start of each line\n"+
				"-r Append a row number column as the first column in the output\n"+
				"-d Append f ile text onto the start of each line\n"+
				"-s Skip blank lines and those with less than the indicated number of columns.\n"+
				"-a Print all available columns.\n"+

				"\n" +
				"Example: java -jar pathTo/T2/PrintSelectColumns -f /TabFiles/ -i 0,3,9 -n 1 -c chr\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new PrintSelectColumns(args);
	}

}


