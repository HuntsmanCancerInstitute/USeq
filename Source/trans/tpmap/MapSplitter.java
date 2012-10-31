package trans.tpmap;
import java.io.*;
import java.util.*;

import util.gen.*;

/**
 * Splits a '.bpmap' text file into different chromosomes, helper script for {@link TPMapProcessor}.
 *
 */
public class MapSplitter {
	private ArrayList info;
	private int[] startPositions;
	private int[] matches;
	private int numberFeatures;
	private MapFeature[] features;
	private File saveDirectory;
	
	/**For running in combination with BPMapProcessor*/
	public MapSplitter (MapFeature[] features, File saveDirectory){
		this.features = features;
		this.saveDirectory = saveDirectory;
		numberFeatures = features.length;
		splitTPMap();
	}
	public void splitTPMap(){
		//get information about .bpmap file 
		//chromosome text, start index, stop index, length;
		System.out.println("\tGathering information about tpmap MapFeature[] file...");
		gatherInfo();	
		//System.out.println(info+"\n"+info.size());
		
		//save bpmapInfo to disk for later use
		IO.saveObject(new File(saveDirectory, "tpmap.faInfo"), info);
		
		System.out.println("\tSplitting features into different chromosomes...");
		//break up startPositions and duplicates into individual int[], one for each chromosome
		breakSavePositions();
	}
	public void breakSavePositions(){
		//initialize int[] arrays for each chromosome
		int size = info.size();	
		int[] starts;
		int[] hits;
		String chromName;
		int startIndex;
		int length;
		//start at 1 since 1st num in the number of data lines
		for (int i=1; i< size; i+=4){
			chromName = (String)info.get(i);
			startIndex = ((Integer)info.get(i+1)).intValue();
			length = ((Integer)info.get(i+3)).intValue();
			starts = new int[length];
			System.arraycopy(startPositions,startIndex,starts,0,length);
			//save starts to disk using the chrom text
			IO.saveObject(new File(saveDirectory,"tpmap.fa"+chromName), starts);
			//save hits/ matches to genome?
			if (matches[0] != -1){
				hits = new int[length];
				System.arraycopy(matches,startIndex,hits,0,length);
				//save starts to disk using the chrom text
				IO.saveObject(new File(saveDirectory,"tpmap.fa"+chromName+"Matches"), hits);
			}
		}
	}
	
	
	
	/**Use this method to break a int[] array into different int[] by chromosome.
	 * Saves the int[] arrays to disk using the fullPathBaseName+ the chromosome text.
	 * bpmapInfo comes from running the BPMapSplitter on a .bpmap file.*/
	public static String breakSaveIntensityValues(ArrayList bpmapInfo, int[] intensities, String fullPathBaseName){
		//initialize int[] arrays for each chromosome
		int size = bpmapInfo.size();	
		int[] values;	
		String chromName;
		int startIndex;
		int length;
		String uniqueId = ".tmpFile"+new Random().nextInt(10000);
		//start at 1 since 1st num in the number of data lines
		for (int i=1; i< size; i+=4){
			chromName = (String)bpmapInfo.get(i);
			startIndex = ((Integer)bpmapInfo.get(i+1)).intValue();
			length = ((Integer)bpmapInfo.get(i+3)).intValue();
			values = new int[length];
			System.arraycopy(intensities,startIndex,values,0,length);
			//save baseFullPathName to disk appending the chrom text
			IO.saveObject(new File(fullPathBaseName+chromName+uniqueId), values);	
		}
		return uniqueId;
	}	
	/**Use this method to break a float[] array into different float[] by chromosome.
	 * Saves the float[] arrays to disk using the fullPathBaseName+ the chromosome text.
	 * bpmapInfo comes from running the BPMapSplitter on a .bpmap file.*/
	public static String breakSaveIntensityValues(ArrayList bpmapInfo, float[] intensities, String fullPathBaseName){
		//initialize int[] arrays for each chromosome
		int size = bpmapInfo.size();	
		float[] values;	
		String chromName;
		int startIndex;
		int length;
		String uniqueId = ".tmpFile"+new Random().nextInt(10000);
		//start at 1 since 1st num in the number of data lines
		for (int i=1; i< size; i+=4){
			chromName = (String)bpmapInfo.get(i);
			startIndex = ((Integer)bpmapInfo.get(i+1)).intValue();
			length = ((Integer)bpmapInfo.get(i+3)).intValue();
			values = new float[length];
			System.arraycopy(intensities,startIndex,values,0,length);
			//save baseFullPathName to disk appending the chrom text
			IO.saveObject(new File(fullPathBaseName+chromName+uniqueId), values);	
		}
		return uniqueId;
	}	
	/**Use this method to break a int[][] array into different int[][] by chromosome.
	 * Saves the int[][] arrays to disk using the fullPathBaseName+ the chromosome text.
	 * bpmapInfo comes from running the BPMapSplitter on a .bpmap file.*/
	public static String breakSaveIntensityValues(ArrayList bpmapInfo, int[][] intensities, String fullPathDir){
		//initialize int[] arrays for each chromosome
		int size = bpmapInfo.size();	
		int[][] values;	
		String chromName;
		int startIndex;
		int length;
		int numInt = intensities.length;
		String uniqueId = ".tmpFile"+new Random().nextInt(10000);
		//start at 1 since 1st num in the number of data lines
		for (int i=1; i< size; i+=4){
			chromName = (String)bpmapInfo.get(i);
			startIndex = ((Integer)bpmapInfo.get(i+1)).intValue();
			length = ((Integer)bpmapInfo.get(i+3)).intValue();
			values = new int[numInt][length];
			for (int j=0; j<numInt; j++) System.arraycopy(intensities[j],startIndex,values[j],0,length);
			//save baseFullPathName to disk appending the chrom text
			IO.saveObject(new File(fullPathDir,chromName+uniqueId), values);	
		}
		return uniqueId;
	}
	/**Use this method to break a float[][] array into different float[][] by chromosome.
	 * Saves the float[][] arrays to disk using the fullPathBaseName+ the chromosome text.
	 * bpmapInfo comes from running the TPMapSplitter on a .tpmap file.*/
	public static String breakSaveIntensityValues(ArrayList bpmapInfo, float[][] intensities, String fullPathDir){
		//initialize int[] arrays for each chromosome
		int size = bpmapInfo.size();	
		float[][] values;	
		String chromName;
		int startIndex;
		int length;
		int numInt = intensities.length;
		String uniqueId = ".tmpFile"+new Random().nextInt(100000);
		//start at 1 since 1st num in the number of data lines
		for (int i=1; i< size; i+=4){
			chromName = (String)bpmapInfo.get(i);
			startIndex = ((Integer)bpmapInfo.get(i+1)).intValue();
			length = ((Integer)bpmapInfo.get(i+3)).intValue();
			values = new float[numInt][length];
			for (int j=0; j<numInt; j++) System.arraycopy(intensities[j],startIndex,values[j],0,length);
			//save baseFullPathName to disk appending the chrom text
			IO.saveObject(new File(fullPathDir,chromName+uniqueId), values);	
		}
		return uniqueId;
	}	
	
	public void gatherInfo(){
		startPositions = new int[numberFeatures];
		matches = new int[numberFeatures];
		info = new ArrayList();
		//add number of lines in bpmap file
		info.add(new Integer(numberFeatures));
		String chrom = "";
		int counter =0;
		int start = 0;
		int stop = 0;
		for (int x=0; x< numberFeatures; x++){
			//check for chromosome change
			if (features[x].chromosome.equals(chrom)==false){
				if (counter!=0){
					stop = counter-1;
					info.add(new Integer(stop));
					info.add(new Integer(stop-start+1));
				}
				chrom = features[x].chromosome;
				start = counter;
				info.add(chrom);
				info.add(new Integer(start));
			}
			//save startPosition and matches to genome
			startPositions[counter] = features[x].start;
			matches[counter] = features[x].matches;
			counter++;
		}
		//calculate last info
		stop = counter-1;
		info.add(new Integer(stop));
		info.add(new Integer(stop-start+1));			
	}
}
