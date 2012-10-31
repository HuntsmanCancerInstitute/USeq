package edu.expr;
import java.util.*;
import java.io.*;
import java.util.regex.*;
import util.gen.*;

/**Parses Agilent data based on a list of column names. */
public class ParseAgilentData {

		//fields
		private File[] agilentFiles;
		private boolean removeControls = false;
		private boolean removeDuplicates = false;
		private String[] columns2Parse;
		private String controlTypeHeadingName = "ControlType";
		private String uniqueIDHeadingName = "ProbeName";

		
		//constructor
		public ParseAgilentData( String[] args){
			//process arguments
			processArgs(args);
			
			//for each file parse
			for (int i =0; i< agilentFiles.length; i++){
				System.out.println("\nParsing "+agilentFiles[i]+"...");
				
				//parse it
				parseAgilentFile(agilentFiles[i]);
				
				
			}
			
			System.out.println("Done!");
		}
		
		public boolean parseAgilentFile (File dataFile){
			String line = "";
			try {
				File parsedFile = new File(dataFile.getParentFile(), Misc.removeExtension(dataFile.getName())+"_Parsed.xls");
				BufferedReader in = new BufferedReader (new FileReader (dataFile));
				LinkedHashMap columnHeadings = null;
			
				//find and parse column heading from FEATURES row
				boolean featuresFound = false;
				while ((line = in.readLine()) != null){
					if (line.startsWith("FEATURES")){
						columnHeadings = new LinkedHashMap();
						String[] features = line.split("\\t");
						for (int i=0; i< features.length; i++){
							columnHeadings.put(features[i], new Integer(i));
						}
						featuresFound = true;
						break;
					}
				}
				if (featuresFound == false){
					System.err.println("\nError: FEATURES line not found in file?! "+dataFile);
					return false;
				}
				
				//find ControlType index?
				int controlTypeIndex = 0;
				if (removeControls){
					if (columnHeadings.containsKey(controlTypeHeadingName) == false) Misc.printExit("\nError: cannot find the '"+ controlTypeHeadingName +"' heading, aborting.\n");
					controlTypeIndex = ((Integer)columnHeadings.get(controlTypeHeadingName)).intValue();
				}
				
				//find uniqueProbeName index?
				int uniqueProbeNameIndex = 0;
				HashSet uniqueIDs = null;
				if (removeDuplicates){
					if (columnHeadings.containsKey(uniqueIDHeadingName) == false) Misc.printExit("\nError: cannot find the '"+ uniqueIDHeadingName +"' heading, aborting.\n");
					uniqueProbeNameIndex = ((Integer)columnHeadings.get(uniqueIDHeadingName)).intValue();
					uniqueIDs = new HashSet();
				}
				
				//find indexes to parse
				int[] indexes2Parse = new int[columns2Parse.length];
				for (int i=0; i< columns2Parse.length; i++){
					if (columnHeadings.containsKey(columns2Parse[i]) == false) {
						Misc.printExit("\nError: cannot find '"+columns2Parse[i]+ "' in available headings "+columnHeadings.keySet()+" aborting.\n");
					}
					indexes2Parse[i] = ((Integer)columnHeadings.get(columns2Parse[i])).intValue();
				}
				
				//print header
				PrintWriter out = new PrintWriter (new FileWriter (parsedFile));
				out.println(Misc.stringArrayToString(columns2Parse, "\t"));
				
				//run through remaining
				int numberDataLines = 0;
				int numberDuplicates = 0;
				int numberControls = 0;
				while ((line = in.readLine()) != null){
					if (line.startsWith("DATA") == false) continue;
					String[] items = line.split("\\t");
					//remove controls?
					if (removeControls && items[controlTypeIndex].equals("0") == false ) {
						numberControls ++;
						continue;
					}
					//check for duplicates?
					if (removeDuplicates){
						if (uniqueIDs.contains(items[uniqueProbeNameIndex])) {
							numberDuplicates++;
							continue;
						}
						else {
							uniqueIDs.add(items[uniqueProbeNameIndex]);
						}
					}
					//print line
					StringBuffer sb = new StringBuffer(items[indexes2Parse[0]]);
					for (int i=1; i<indexes2Parse.length; i++){
						sb.append("\t");
						sb.append(items[indexes2Parse[i]]);
					}
					out.println(sb);
					numberDataLines++;
				}
				in.close();
				out.close();
				if (removeControls) System.out.println("\tNumber controls "+numberControls);
				if (removeDuplicates) System.out.println("\tNumber duplicates "+numberDuplicates);
				System.out.println("\tNumber data lines "+numberDataLines);
				return true;
			} catch (Exception e){
				System.out.println("\nError parsing "+dataFile.getName()+" line "+line+"\n");
				e.printStackTrace();
				return false;
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
						case 'f': agilentFiles = IO.extractFiles(new File(args[i+1]));  i++; break;
						case 'r': removeControls = true; break;
						case 'd': removeDuplicates = true; break;
						case 'e': uniqueIDHeadingName = args[i+1]; i++; removeDuplicates = true; break;
						case 'c': columns2Parse = args[i+1].split(","); i++; break;
						case 'n': controlTypeHeadingName = args[i+1]; i++; removeControls = true; break;
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
			if (agilentFiles == null || agilentFiles.length ==0) Misc.printExit("\nCannot find your Agilent data file(s)?\n");
			if (columns2Parse == null) Misc.printExit("\nWhich columns do you want to parse?\n");
		}	
		
		public static void printDocs(){ 
			System.out.println("\n" +
					"**************************************************************************************\n" +
					"**                         Parse Agilent Data: June 2007                         **\n" +
					"**************************************************************************************\n" +
					"Parses Agilent microarray data files based on particular column names.\n\n" +
					
					"-f Full path file text or directory containing text Agilent data files.\n" +
					"-c Comma delimited, no spaces, list of columns to parse, see FEATURES row.\n"+
					"-r Remove control probe data.\n"+
					"-d Remove duplicates based on a unique probe text column.\n"+
					"-e Name of unique probe text column, defaults to 'ProbeName'.\n"+
					"-n Name of control probe column, defaults to 'ControlType', non controls are '0'.\n\n" +
					
					"Example: java -jar pathTo/T2/Apps/ParseAgilentData -f /badData/ -r -c\n" +
					"      ProbeName,GeneName,Description,LogRatio,IsManualFlag\n\n" +
					
			"**************************************************************************************\n");		
		}
		
		public static void main(String[] args) {
			if (args.length == 0) {
				printDocs();
				System.exit(0);
			}
			else new ParseAgilentData(args);
		}
		
	}


