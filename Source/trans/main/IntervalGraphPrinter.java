package trans.main;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import util.gen.*;

/**
 * Prints an array of Interval as an '.sgr' file for import into IGB.
 */
public class IntervalGraphPrinter {
	//fields
	private File[] files;
	private Interval[] intervals;
	private double scoreCutOff =0;
	private int numberCutOff =0;
	private String chromosomePrefix = "";
	private boolean printRatios = true;
	private boolean printScores = true;
	private int defaultScoreIndex = 1;
	
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
					case 's': scoreCutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'r': numberCutOff =Integer.parseInt(args[i+1]); i++; break;
					case 'c': chromosomePrefix =args[i+1]; i++; break;
					case 'g': printScores=false; break;
					case 'd': printRatios=false; break;
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
			System.out.println("\nEnter a serialized Interval[] array file or directory!\n");
			System.exit(0);
		}
		
		//get files to process
		if (directory.isDirectory()){
			files = directory.listFiles();
			if 	(files == null || files.length==0){
				System.out.println("Cannot find the directory or files in the directory?!\n\n");
				System.exit(0);	
			}
		}
		else files = new File[]{directory};	
		
	}
	
	public IntervalGraphPrinter (String[] args){
		processArgs(args);
		printSaveIntervalScoreGraphFiles();
	}
	
	
	/**Prints and saves a Affy .sgr file*/
	public void printSaveIntervalScoreGraphFiles(){
		//for each file
		for (int i=0; i<files.length; i++){
			//get intervals, watch for class cast exceptions
			try {
				intervals = (Interval[])IO.fetchObject(files[i]);
				sortIntervals();
				System.out.println("Creating an .sgr and .bed files for "+files[i].getName());
				int numIntervals = intervals.length;
				
				//run thru intervals saving an ArrayList of BindingRegion s
				ArrayList goodIntervals = new ArrayList();
				double score;
				if (numberCutOff!=0 && numberCutOff<numIntervals) numIntervals = numberCutOff;
				for (int j=0; j<numIntervals; j++){
					//make entry line
					score = intervals[j].getSortBy();
					if (score>= scoreCutOff) goodIntervals.add(intervals[j]);
				}
				//recycle intervals to include only good intervals
				intervals = new Interval[goodIntervals.size()];
				goodIntervals.toArray(intervals);
				//print good intervals
				File sgrFile = new File (files[i].getCanonicalPath()+".sgr");
				if (intervals[0].getOligos() != null && printScores){
					printScores(intervals, sgrFile, printRatios);
				}
				else IO.writeString(convertToSGR(intervals), sgrFile);
				//zip it
				IO.zipAndDelete(sgrFile);
				//print bed file
				IO.writeString(convertToBED(intervals), files[i].getCanonicalPath()+".bed");
				
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		System.out.println("\nDone!\n");
	}
	
	/**Prints to a File a line for every Oligo in each Interval.
	 *Tab delimited: Chromosome, Start of the oligo, the log2(aveT/ aveC ratio).
	 **/
	public static void printScores(Interval[] intervals, File file, boolean printRatios){
		double log2 = Math.log(2);
		try{
			PrintWriter out = new PrintWriter (new FileWriter(file));
			int num = intervals.length;
			int numTs = intervals[0].getNumberTreatmentIntensities();
			int numCs = intervals[0].getNumberControlIntensities();
			for (int i=0; i<num; i++){
				Oligo[] oligos = intervals[i].getOligos();
				int numOligos = oligos.length;
				for (int j=0; j<numOligos; j++){
					StringBuffer sb = new StringBuffer();
					sb.append(intervals[i].getChromosome());
					sb.append("\t");
					sb.append(oligos[j].getStart());
					sb.append("\t");
					double ratio;
					if (printRatios) ratio = Num.mean(oligos[j].getTreatmentIntensities(numTs)) / Num.mean(oligos[j].getControlIntensities(numCs));
					else {
						double t = Num.mean(oligos[j].getTreatmentIntensities(numTs));
						double c = Num.mean(oligos[j].getControlIntensities(numCs));
						ratio = Num.relativeDifference(t,c);
					}
					sb.append(Math.log(ratio)/log2);
					out.println(sb);
				}
				
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Converts an array of Interval to a .sgr String format.*/
	public String convertToSGR(Interval[] intervals){
		int num = intervals.length;
		StringBuffer b = new StringBuffer();
		for (int i=0; i< num; i++){
			String chrom = chromosomePrefix+intervals[i].getChromosome() +"\t";
			int end = intervals[i].getStartLastOligo()+intervals[i].getSizeOfOligoMinusOne();
			//make zero entry stop
			b.append(chrom);
			b.append((end+1));
			b.append("\t");
			b.append(0);
			b.append("\n");
			//make stop 
			b.append(chrom);
			b.append(end);
			b.append("\t");
			b.append(intervals[i].getSortBy());
			b.append("\n");
			//make start
			b.append(chrom);
			b.append(intervals[i].getStart1stOligo());
			b.append("\t");
			b.append(intervals[i].getSortBy());
			b.append("\n");
			//make zero
			b.append(chrom);
			b.append((intervals[i].getStart1stOligo()-1));
			b.append("\t");
			b.append(0);
			b.append("\n");
		}
		return b.toString();
	}
	
	/**Converts an array of Interval to a .bed String format.*/
	public String convertToBED(Interval[] intervals){
		int num = intervals.length;
		StringBuffer b = new StringBuffer();
		for (int i=0; i< num; i++){
			b.append(chromosomePrefix);
			b.append(intervals[i].getChromosome());
			b.append("\t");
			b.append(intervals[i].getStart1stOligo());
			b.append("\t");
			int end = intervals[i].getStartLastOligo()+intervals[i].getSizeOfOligoMinusOne();
			b.append(end);
			b.append("\n");
		}
		return b.toString();
	}
	
	/**Sorts Interval[] based on subWindow if present, if not, then uses best window*/
	public void sortIntervals(){
		int numIntervals = intervals.length;
		SubWindow sub;
		for (int j=0; j< numIntervals; j++){
			sub = intervals[j].getBestSubWindow();
			if (sub != null)intervals[j].setSortBy(sub.getMedianRatio());
			else {
				if (j==0) System.out.println("Warning: no sub window! Setting sort by to best window score, index "+defaultScoreIndex);
				double score = intervals[j].getBestWindow().getScores()[defaultScoreIndex];
				intervals[j].setSortBy(score);
			}
		}	
		Arrays.sort(intervals);
	}
	
	//main
	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);
		}		
		new IntervalGraphPrinter(args);
	}	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Interval Graph Printer: Jun 2006                          **\n" +
				"**************************************************************************************\n" +
				"IGP converts serialized Interval[] arrays to graph files (.sgr & .bed) for import into\n" +
				"Affymetrix' IGB. Files will be written to the data directory. Intervals are sorted by\n" +
				"the median ratio of the best sub window, if present, or the best window. \n" +
				"\n" +
				"Use the following options when running IGP:\n\n" +
				"-f Full path file text for the Interval[], if a directory is specified, all files\n" +
				"      within will be processed. (Required)\n" +
				"-s Score cut off, print everything above this score, defaults to all.\n" +
				"-r Rank cut off, prints everything above this rank, ie the top 200, defaults to all.\n" +
				"-g Print goal post Interval representations using the sub window score, if present, or\n" +
				"      the best window score, default is to print individual oligo log2(aveT/aveC)\n" +
				"      scores.\n" +
				"-c Chromosome prefix to prepend onto each .sgr line. (This is sometimes needed to match\n" +
				"      IGB's quickload seqIDs.)\n"+
				"\n" +
				"Example: java -jar pathTo/T2/Apps/IntervalGraphPrinter -f /my/affy/res/ -s 50 -r 200 \n" +
				"\n" +

		"**************************************************************************************\n");
	}
}
