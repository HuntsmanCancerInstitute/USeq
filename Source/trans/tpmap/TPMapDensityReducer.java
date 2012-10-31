package trans.tpmap;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;

import util.gen.*;

/**
 * Randomly drops a certain fraction of tpmap lines.
 */
public class TPMapDensityReducer {
	//fields
	private File tpmapFile;
	private double fractionToDrop;
	
	public TPMapDensityReducer(String[] args){
		processArgs(args);
		
		try {
			//create Random
			Random randomGenerator = new Random();
			//read in 
			String line;
			BufferedReader in = new BufferedReader(new FileReader(tpmapFile));
			PrintWriter out = new PrintWriter(new FileWriter(tpmapFile.getCanonicalPath()+".reduced"));
			double total =0;
			double numberSaved =0;
			while ((line = in.readLine()) !=null) { 
				total++;
				double test = randomGenerator.nextDouble();
				if (test >= fractionToDrop){
					out.println(line);
					numberSaved++;
				}
				
			}			
			
			in.close();
			out.close();
			
			double fractionDropped = 1- numberSaved/total;
			System.out.println("Fraction dropped = "+Num.formatNumber(fractionDropped, 4)+" ("+(int)(total-numberSaved)+"/"+(int)total+")\n");
			
		} catch (Exception e) {e.printStackTrace();}
		
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
					case 't': tpmapFile = new File(args[i+1]); i++; break;
					case 'f': fractionToDrop =Double.parseDouble(args[i+1]);i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test+"\n");
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
				"**                            TPMap Density Reducer: Feb 2007                       **\n" +
				"**************************************************************************************\n" +
				"Randomly drops a certain percentage of the total TPMapLines. Fraction to drop =\n" +
				"      1- (current bp tiling density/ desired density)\n\n"+
				
				"Required Parameters:\n"+
				"-t Full path file text for the text tpmap file.\n" +
				"-f Fraction of tpmap lines to drop."+

				"\n" +
				"Example: java -Xmx500M -jar pathTo/T2/Apps/TPMapDensityReducer -f /bpmap.txt -f 0.25\n" +
				"\n" +
			
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		System.out.println("\nLaunching TPMapDensityReducer...");
		new TPMapDensityReducer(args);
	}

	


	
}
