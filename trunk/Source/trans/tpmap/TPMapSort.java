package trans.tpmap;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;

import util.gen.*;

/**
 * This sorts a tpmap first by chromosome, second by base position.
 */
public class TPMapSort {
	//fields
	private File tpmapFile;

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
					case 'f': tpmapFile = new File(args[i+1]); i++; break;
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
		if (tpmapFile==null || tpmapFile.canRead() == false){
			System.out.println("\nCannot find your bpmapFile file!\n");
			System.exit(0);
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            TPMap Sort: Jan 2005                                  **\n" +
				"**************************************************************************************\n" +
				"Sorts a text bpmap file based on chromosome and oligo start positions.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path file text for the text tpmap file.\n" +

				"\n" +
				"Example: java -Xmx500M -jar pathTo/T2/Apps/TPMapSort -f /affy/tpmap.txt\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length !=2){
			printDocs();
			System.exit(0);
		}
		System.out.println("Launching TPMapSort...");
		new TPMapSort(args);
	}

	
	public TPMapSort(String[] args){
		processArgs(args);
		
		try {
			long start = System.currentTimeMillis();
			//count lines
			System.out.println("Counting # lines in file...");
			int numLines = (int)IO.countNumberOfLines(tpmapFile);
			
			//read in and make array
			System.out.println("Making objects...");
			TPMapLine[] lines = new TPMapLine[numLines];
			String line;
			BufferedReader in = new BufferedReader(new FileReader(tpmapFile));
			int counter =0;
			while ((line = in.readLine()) !=null) { 
				lines[counter++] = new TPMapLine(line);
			}			
			
			in.close();
			
			//sort array
			System.out.println("Sorting...");
			Arrays.sort(lines);
			
			//print out sorted array
			System.out.println("Writing new file...");
			PrintWriter outRes = new PrintWriter(new FileWriter(tpmapFile.getCanonicalPath()+"Sorted"));
			for (int i=0; i<numLines; i++){
				outRes.println(lines[i].getLine());
			}

			int elapse = (int)(System.currentTimeMillis()-start)/1000;
			System.out.println("Finished "+elapse+" seconds");
			outRes.close();
			
		} catch (Exception e) {e.printStackTrace();}
		
	}
	

	
}
