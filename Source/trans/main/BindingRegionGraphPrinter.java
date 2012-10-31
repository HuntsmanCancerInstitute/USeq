package trans.main;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import trans.anno.*;
import util.gen.*;


/**
 * For printing chromosome specific .sgr files from IntervalPloter .txt picks files.
 * Tab delimited: rank subWinMedianScore trimmedMeanPickScore chrom start stop seq.
 */
public class BindingRegionGraphPrinter {
	//fields
	private File[] files;
	private String chromosomePrefix = "chr";
	
	public BindingRegionGraphPrinter(String[] args){
		processArgs(args);
		
		//for each picks.txt file
		for (int i=0; i<files.length; i++){
			//parse picks
			BindingRegion[] bindingRegions = AnnotateRegions.parsePicksFile(files[i],0);
			//find chromosomes
			Arrays.sort(bindingRegions);		
			//make array of PrintWriter to print out files, add first entry
			PrintWriter out = null;
			try{
				out = new PrintWriter( new FileWriter( files[i].getCanonicalPath()+".sgr" ) );
				
			}catch (IOException e){
				e.printStackTrace();
				System.exit(1);
			}
			//run thru array of BindingRegion and add lines,
			int numBindingRegions = bindingRegions.length;
			for (int j=0; j<numBindingRegions; j++){
				double score = bindingRegions[j].getScore() * 10;
				//print 
				StringBuffer sb = new StringBuffer();
				String chromName = chromosomePrefix+bindingRegions[j].getChromosome() + "\t";
				//make zero entry stop
				sb.append(chromName);
				sb.append((bindingRegions[j].getEnd()+1));
				sb.append("\t");
				sb.append(0);
				sb.append("\n");
				//make stop entry?
				if (bindingRegions[j].getEnd() != bindingRegions[j].getStart() ){
					sb.append(chromName);
					sb.append(bindingRegions[j].getEnd());
					sb.append("\t");
					sb.append(score);
					sb.append("\n");
				}
				//make start
				sb.append(chromName);
				sb.append(bindingRegions[j].getStart());
				sb.append("\t");
				sb.append(score);
				sb.append("\n");
				//make zero
				sb.append(chromName);
				sb.append((bindingRegions[j].getStart()-1));
				sb.append("\t");
				sb.append(0);
				sb.append("\n");
				out.print(sb);
			}
			//close print writers
			out.close();
		}
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		File directory = null;
		
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
		if (directory==null || directory.exists() == false){
			System.out.println("\nEnter a binding region text file or directory!\n");
			System.exit(0);
		}
		
		if (directory.isDirectory()) files = IO.extractFiles(directory , ".txt");	
		else {
			files = new File[1];
			files[0] = directory;
		}
		
	}	
	
	//main
	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}		
		new BindingRegionGraphPrinter(args);
	}	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**          Binding GenomicRegion Graph Printer: March 2005                                **\n" +
				"**************************************************************************************\n" +
				"BRGP converts a Binding GenomicRegion Picks txt file (from the IntensityPlotter: rank\n" +
				"subWinMedianScore trimmedMeanPickScore chrom start stop seq) to an .sgr file for\n" +
				"import into Affymetrix's IGB.\n"+
				"\n" +
				"Use the following options when running BRGP:\n" +
				"-f Full path file text for the .txt picks file, if a directory is specified, all '.txt'\n" +
				"      files found within will be processed. (Required)\n" +
				"\n" +
				"Example: java BindingRegionGraphPrinter -f /my/affy/res/hbPicks.txt \n" +
				"\n" +
		"**************************************************************************************\n");
	}
}
