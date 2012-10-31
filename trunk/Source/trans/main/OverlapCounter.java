package trans.main;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.*;
/**
 * Provides an intersection analysis on multiple Interval[]s.
 *
 */
public class OverlapCounter {
	//fields
	private int suggestedNumber = 10000000;
	private boolean writeFiles = false; //default is to not write files
	private File[] files;
	private int maxGap = 0;
	private boolean comparePeaks = false;
	
	public OverlapCounter(String[] args){
		processArgs(args);
		
		//open each file, sort and save, also check for class cast exceptions
		System.out.println("\nLaunching OverlapCounter:\n\tMaxGap "+maxGap);
		
		//compare peaks, find closest peak print distance
		if (comparePeaks){
			System.out.println("\nFinding closest Interval peak...");
			for (int i=0; i<files.length; i++){
				Interval[] one = (Interval[])IO.fetchObject(files[i]);
				int j = i+1;
				for (; j<files.length; j++){
					Interval[] two = (Interval[])IO.fetchObject((files[j]));
					findClosestPeaks(one, two, files[i], files[j]);
					findClosestPeaks(two, one, files[j], files[i]);
				}
			}
		}
		//compare intervals
		else {
			System.out.println("\tComparing Intervals...");
			//make table to hold info
			String[][] table = new String[files.length][files.length];
			for (int i=0; i<files.length; i++){
				Interval[] one = (Interval[])IO.fetchObject(files[i]);
				int j = i+1;
				for (; j<files.length; j++){
					Interval[] two = (Interval[])IO.fetchObject((files[j]));
					table[j][i] = testOverlap(one, two, files[i], files[j]);
				}
			}	
			//print table headder, 1st row names
			System.out.println();
			String[] names = new String[files.length];
			for (int i=0; i<files.length; i++){
				names[i] = files[i].getName();
			}
			names = Misc.trimCommon(names);
			
			//print table
			printTable(table, names);	
		}
	}
	
	public static void printTable(Object[][] t, String[] names){
		for (int x=0; x< names.length-1; x++){
			System.out.print("\t"+names[x]);
		}
		System.out.println();
		for (int x=1; x<t.length; x++){
			System.out.print(names[x]+"\t");
			for (int z =0; z<t[x].length; z++){
				if (t[x][z] != null) System.out.print(t[x][z]+"\t");
				else System.out.print("\t");
			}
			System.out.println();
		}
	}
	
	public void findClosestPeaks(Interval[] one, Interval[] two, File fileOne, File fileTwo){
		int oneLen = one.length;
		int twoLen = two.length;
		if (oneLen> suggestedNumber) oneLen = suggestedNumber;
		if (twoLen> suggestedNumber) twoLen = suggestedNumber;
		String nameOne = fileOne.getName().replaceFirst(".tmpFile\\d+","");
		String nameTwo = fileTwo.getName().replaceFirst(".tmpFile\\d+","");
		System.out.println("Comparing "+nameOne+ " to "+nameTwo + " (coordinates, distance)");
		for (int i=0; i<oneLen; i++){
			int gap = 1000000000;
			Interval closest = null;
			for (int j=0; j< twoLen; j++){
				int distance = overlapPeaks(one[i],two[j]);
				if (distance !=-1 && distance< gap) {
					gap = distance;
					closest = two[j];
				}
				if (gap == 0) break;
			}
			String id = "";
			if (closest !=null ) id = closest.getChromosome()+":"+closest.getStart1stOligo()+"-"+closest.getStartLastOligo()+closest.getSizeOfOligoMinusOne();
			System.out.println(one[i].getChromosome()+":"+one[i].getStart1stOligo()+"-"+one[i].getStartLastOligo()+one[i].getSizeOfOligoMinusOne()+"\t"+gap +"\t"+id);
		}
	}
	
	/**Returns distance between top peaks. 
	 * if no peak uses middle of sub window
	 * if no sub window uses middle of interval
	 * returns -1 if different chromosomes, 0 if direct match, 1 if immediately adjacent, 2 if 1bp gap...*/
	public int overlapPeaks (Interval one, Interval two){
		//do chromosomes match
		if (one.getChromosome().equals(two.getChromosome()) == false ) return -1;
		int peakOne = findPeak(one);
		int peakTwo = findPeak(two);
		return Math.abs(peakOne-peakTwo);
	}
	
	/**Returns base position of highest scoring peak, 
	 * if no peak uses middle of sub window
	 * if no sub window uses middle of interval*/
	public static int findPeak(Interval i){
		//peak
		BindingPeak[] peaks = i.getBindingPeaks(); //this is a sorted list highest to lowest
		if (peaks!=null) return peaks[0].getPeakBP();
		//sub window
		SubWindow sub = i.getBestSubWindow();
		if (sub != null){
			Oligo[] oligos = sub.getOligos();
			int start = oligos[0].getStart();
			int stop = oligos[oligos.length-1].getStart();
			float diff = stop-start;
			return Math.round(diff/2) + start;
		}
		//interval
		float diff = (i.getStartLastOligo()+i.getSizeOfOligoMinusOne()) - i.getStart1stOligo();
		return Math.round(diff/2) + i.getStart1stOligo();
	}
	
	/**checks to see if intervals overlap by the minimum maxGap, can set maxGap negative to require an overlap.*/
	public boolean overlap (Interval one, Interval two){
		int oneLast = one.getStartLastOligo()+one.getSizeOfOligoMinusOne();
		int overlap;
		//to left
		if (oneLast< two.getStart1stOligo()){
			//check distance
			overlap = two.getStart1stOligo()-oneLast;
			if ( overlap <=maxGap )return true;
			else return false;
		}
		int twoLast = two.getStartLastOligo()+two.getSizeOfOligoMinusOne();
		//to right
		if (one.getStart1stOligo()> twoLast) {
			//check distance
			overlap = one.getStart1stOligo()-twoLast;
			if ( overlap <=maxGap ) return true;
			else return false;
		}
		//by default they overlap
		//entirely contained within?
		if ((one.getStart1stOligo()<=two.getStart1stOligo() && oneLast>=twoLast) || (two.getStart1stOligo()<=one.getStart1stOligo() && twoLast>= oneLast) ) {
			return true;
		}
		//partial overlap
		//two right of one
		if (one.getStart1stOligo()< two.getStart1stOligo()) overlap = two.getStart1stOligo()- oneLast; //want it negative
		//two left of one
		else overlap = one.getStart1stOligo()-twoLast; //want it negative
		if (overlap<= maxGap) return true;
		return false;
	}
	
	/**Returns the average % intersection.*/
	public String testOverlap(Interval[] one, Interval[] two, File fileOne, File fileTwo){
		int oneLen = one.length;
		int twoLen = two.length;
		if (oneLen> suggestedNumber) oneLen = suggestedNumber;
		if (twoLen> suggestedNumber) twoLen = suggestedNumber;
		
		String nameOne = fileOne.getName().replaceFirst(".tmpFile\\d+","");
		String nameTwo = fileTwo.getName().replaceFirst(".tmpFile\\d+","");
		
		//make HashMaps to hold Interval overlaps, need hash to kill any potential dups
		HashMap overlapOne = new HashMap();
		HashMap overlapTwo = new HashMap();
		//print overlap files?
		if (writeFiles){
			
			for (int i=0; i<oneLen; i++){
				for (int j=0; j< twoLen; j++){
					//check chromosome && overlap
					if (two[j].getChromosome().equals(one[i].getChromosome()) && (overlap(one[i], two[j]))){
						String oneSig = one[i].getChromosome()+one[i].getStart1stOligo()+"_"+one[i].getStartLastOligo();
						String twoSig = two[j].getChromosome()+two[j].getStart1stOligo()+"_"+two[j].getStartLastOligo();
						overlapOne.put(oneSig, one[i]);
						overlapTwo.put(twoSig, two[j]);
					}
				}
			}
			//Write arrays
			//write 1st overlap array
			Interval[] ints = new Interval[overlapOne.size()];
			Iterator it = overlapOne.keySet().iterator();
			int counter = 0;
			while (it.hasNext()){
				ints[counter++] = (Interval)overlapOne.get((String)it.next());
			}
			IO.saveObject(new File(fileOne.getParent()+File.separator+nameOne+"_Vs_"+nameTwo+".1stOverlap"), ints);
			//write 1st non overlap array
			ints = subtractIntervals(one, ints);
			IO.saveObject(new File(fileOne.getParent()+File.separator+nameOne+"_Vs_"+nameTwo+".1stNonOverlap"), ints);
			//write 2nd
			ints = new Interval[overlapTwo.size()];
			it = overlapTwo.keySet().iterator();
			counter = 0;
			while (it.hasNext()){
				ints[counter++] = (Interval)overlapTwo.get((String)it.next());
			}
			IO.saveObject(new File(fileTwo.getParent()+File.separator+nameOne+"_Vs_"+nameTwo+".2ndOverlap"), ints);
			ints = subtractIntervals(two, ints);
			IO.saveObject(new File(fileTwo.getParent()+File.separator+nameOne+"_Vs_"+nameTwo+".2ndNonOverlap"), ints);
		}
		else {
			for (int i=0; i<oneLen; i++){
				for (int j=0; j< twoLen; j++){
					//check chromosome
					if (two[j].getChromosome().equals(one[i].getChromosome())){
						//check overlap
						if (overlap(one[i], two[j])){
							String oneSig = one[i].getChromosome()+one[i].getStart1stOligo()+"_"+one[i].getStartLastOligo();
							String twoSig = two[j].getChromosome()+two[j].getStart1stOligo()+"_"+two[j].getStartLastOligo();
							overlapOne.put(oneSig, one[i]);
							overlapTwo.put(twoSig, two[j]);
						}
					}
				}
			}
			
		}
		System.out.println(nameOne+" Vs "+nameTwo+"\t");
		System.out.println("\t"+percent(overlapOne.size(), oneLen));
		System.out.println("\t"+percent(overlapTwo.size(), twoLen));
		double oneF = (double)overlapOne.size()/(double)oneLen;
		double twoF = (double)overlapTwo.size()/(double)twoLen;
		double ave = (oneF + twoF)/2.0;
		if (ave == 0) return "0";
		return Num.formatNumber(100*ave, 1);
	}
	
	public static String percent(int overlaps, int total){
		return 100*(double)overlaps/(double)total+"% ("+overlaps+"/"+total+")";
	}
	
	/**Returns the difference between two interval arrays, be sure to get bigger and smaller correct.*/
	public static Interval[] subtractIntervals(Interval[] bigger, Interval[] smaller){
		int lenBigger = bigger.length;
		int lenSmaller = smaller.length;
		ArrayList al = new ArrayList();
		String compBigger;
		String compSmaller;
		boolean found;
		//for each of the bigger intervals
		for (int i=0; i<lenBigger; i++) {
			//make comp phrase
			compBigger = bigger[i].getChromosome()+bigger[i].getStart1stOligo()+"_"+bigger[i].getStartLastOligo();
			found = false;
			//for each of the smaller intervals look for the bigger comp phrase
			for (int j=0; j<lenSmaller; j++){
				compSmaller = smaller[j].getChromosome()+smaller[j].getStart1stOligo()+"_"+smaller[j].getStartLastOligo();
				//if found set boolean and break looping
				if (compBigger.equals(compSmaller)){
					found = true;
					break;
				}
			}
			//if not found then add to subtracted
			if (found == false) {
				al.add(bigger[i]);
			}
		}
		Interval[] sub = new Interval[al.size()];
		al.toArray(sub);
		return sub;
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
					case 'f': files = new File(args[i+1]).listFiles(); i++; break;
					case 'n': suggestedNumber =Integer.parseInt(args[i+1]); i++; break;
					case 'm': maxGap =Integer.parseInt(args[i+1]); i++; break;
					case 'p': comparePeaks = true; break;
					case 'w': writeFiles = true; break;
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
		if (files==null || files.length==0){
			System.out.println("Cannot find the directory or files in the directory?!\n\n");
			System.exit(0);
		}
	}
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Overlap Counter: March 2005                              **\n" +
				"**************************************************************************************\n" +
				"\nOC performs a pairwise intersection analysis between each Interval[]  file within a\n" +
				"directory.  Sets of overlapping and non-overlapping Intervals, for each file, can be\n" +
				"written to disk. \n\n" +
				
				"Parameters:\n"+
				"-f Full path file text for a directory containing Interval files, required.\n" +
				"-n The number of Intervals to compare. Enter say 200 to compare the top 200 after\n" +
				"       sorting. If they don't exist, all will be used.  Default is all.\n" +
				"-m Max gap between boarders, default is 0bp. Negative values force an overlap.\n"+
				"-p Print distance to closest interval peak.\n"+
				"-w Write the overlapping and non-overlapping Interval files to disk, the default is no.\n\n" +
				
				"Example: java -Xmx128M -jar pathTo/T2/Apps/OverlapCounter -f /affy/res/Intervals/\n" +
				"      -n 250 -p -w -m -100\n\n" +
				
		"**************************************************************************************\n");
	}
	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new OverlapCounter(args);
	}
	
}