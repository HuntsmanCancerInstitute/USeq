package trans.tpmap;
import java.util.*;
import java.io.*;

import util.gen.*;

/**
 * Creates Windows (start position, stop position) of a given length, containing a given number of start positions and maximum number of features.
 * Remember for any given position it may have multiple features, ie spotted multiple times on the array due to multiple hits.
 * Thus a window with 3 positions, each with 10 spottings has 30 features for the window. Can limit the amount added by setting the
 * maxNumFeatures field, it defaults to 1,000,000,000.
 * Used to assemble subsets of the array data for performing statistical tests.  Each window contains a unique set
 * of features, tries to maximize its size to the max size, and is only saved if it contains the minimum number of 
 * features.  
 */
public class WindowMaker {
	//fields
	private int windowSize; 
	private int minNumFeatures;
	private int maxNumFeatures = 1000000000;
	int[][] runs; //final list of start stop indexes
	
	public WindowMaker (int windowSize, int minNumOligos){
		this.windowSize = windowSize;
		this.minNumFeatures = minNumOligos;
	}
	
	public WindowMaker (int windowSize, int minNumOligos, int maxNumFeatures){
		this.windowSize = windowSize;
		this.minNumFeatures = minNumOligos;
		this.maxNumFeatures = maxNumFeatures;
	}
	
	/**Takes int[] files containing bp start positions for each feature, ordered. Saves the results as a 
	 * int[][] where each window has a start index and stop index that correspond to the indexes in the original 
	 * bp start positions.  Indexes of indexes...*/
	public void makeWindowsFromFiles(String[] files){	
		for (int i=0; i<files.length; i++){
			int[] vals = (int[])IO.fetchObject(new File(files[i]));
			//Misc.printArray(vals);
			//Watch out for single oligo hits on a chromosome
			if (vals.length<2) continue;
			makeWindows(vals);
			//System.out.println ("\t"+runs.length + " windows found with at least "+(minNumFeatures)+" features and a max size of "+
				//	windowSize+"bp.");
			//save int[][] to disk
			System.out.println("\tSaving int[][] ... "+files[i]+"Win");
			IO.saveObject(new File(files[i]+"Win"),runs);
		}	
	}
	/**Given a sorted set of values, a start and stop index (inclusive) returns the number of unique values.*/
	public static int numberUniqueValues(int[] values, int startIndex, int stopIndex){
		int totalUnique = 1;
		int currentValue = values[startIndex];
		for (int i=startIndex+1; i<=stopIndex; i++){
			if (currentValue != values[i]){
				//diff value
				totalUnique++;
				currentValue = values[i];
			}
		}
		return totalUnique;
	}
	
	/**Takes an int[] of sorted bp start positions and finds unique windows of 
	 * a maximized width with a minimal # of features.  Returns
	 * an int[window number][start index, stop index]. Stop is inclusive!
	 */
	public int[][] makeWindows(int[] bpStartPositions){
		int length = bpStartPositions.length;
		//only one position?
		if (length == 1 ){
			if (1 >= minNumFeatures) return new int[][] { {0, 0} };
			return new int[0][0];
		}
		int startIndex = 0;
		int endIndex = 1;
		int oldEndIndex =-1;
		int adv = 0;
		ArrayList windows = new ArrayList(length);
		boolean go = true;
		//run through array
		while (go) {
			//grow and make
			// if a properly sized window is found
			if (bpStartPositions[endIndex] - bpStartPositions[startIndex] > windowSize) {
				//back off one
				endIndex--;
				//advance endIndex til bpStartPositions changes (there might be duplicate bpStartPositions[]! hell yes there are!)
				while (true) {
					adv = endIndex + 1;
					//check if can advance stop
					if (endIndex == length) {
						go = false;
						break;
					} 
					//check to see if bpStartPositions are the same
					else {
						if (bpStartPositions[endIndex] != bpStartPositions[adv]) break;
						endIndex++;
					}
				}
				
				//save window if min number of features reached and less than max
				int numOligos = numberUniqueValues(bpStartPositions, startIndex, endIndex);
				int numFeatures = endIndex- startIndex;
				if (numOligos >= minNumFeatures && oldEndIndex!=endIndex && numFeatures <= maxNumFeatures) {
					oldEndIndex = endIndex;
					windows.add(new int[] { startIndex, endIndex });					
				}
				//advance startIndex
				while (bpStartPositions[startIndex] == bpStartPositions[++startIndex]) {}
			}
		
			endIndex++;
		
			//check if at stop
			if (endIndex == length) {				
				endIndex--;
				//save window if min number of features reached
				int numOligos = numberUniqueValues(bpStartPositions, startIndex, endIndex);
				int numFeatures = endIndex- startIndex;
				if (numOligos >= minNumFeatures && numFeatures <= maxNumFeatures) {
					windows.add(new int[] { startIndex, endIndex });
					//System.out.println("start: " + startIndex + " stop: " + endIndex+ " features: "+(1+endIndex-startIndex));
				}
				go = false;
			}			
		}
		//convert ArrayList to int[][]
		int size = windows.size();
		runs = new int[size][];
		windows.toArray(runs);		
		return runs;
	}
	
	public static String[] extractChrFiles (File directory){		
		String[] files;			
		files = directory.list();
		int num = files.length;
		ArrayList chromFiles = new ArrayList();
		try{
			String path = directory.getCanonicalPath() + File.separator;
			for (int i=0; i< num; i++)  {
				if (files[i].indexOf("chr")!=-1 && files[i].endsWith("Matches") == false) chromFiles.add(path+files[i]);
			}
			files = new String[chromFiles.size()];
			chromFiles.toArray(files);
		}catch(IOException e){e.printStackTrace();}
		return files;
	}

	public int getMaxNumFeatures() {
		return maxNumFeatures;
	}

	public void setMaxNumFeatures(int maxNumFeatures) {
		this.maxNumFeatures = maxNumFeatures;
	}
	public static void main(String[] args){
		WindowMaker wm = new WindowMaker(100, 10);
		String[] files = {"/Users/nix/Desktop/Delme/chr12"};
		wm.makeWindowsFromFiles(files);
		
		
	}
}
