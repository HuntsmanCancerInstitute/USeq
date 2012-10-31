package trans.tpmap;
import java.io.*;
import java.util.*;

import util.gen.*;

/**
 * For extracting duplicate oligo/ features from Affy's bpmap file, can be used as a 
 * stand alone but better to use {@link TPMapProcessor}. 
 * 
 * The TPMapDuplicateFilter first creates a hashMap of unique oligo sequences paired with an 
 * 		ArrayList of chip coordinates where a particular oligo was spotted.  In most cases, 
 * 		each oligo was spotted only once. In some cases the same oligo was placed on the array 
 * 		many times.  This repeat information will be used by the CelMapper program to average 
 * 		the intensities of all the repeats, for a particular oligo, and assign this value to 
 * 		each spot.  This repeat information is saved as an unordered array of TPMapFeature[] to disk.  
 * 
 * The TPMapDuplicate Filter then reads line by line thru the bpmap txt file and removes any 
 * 		adjacent repeated sequences.  (After averaging they all will have the same value.)
 * 		It then removes any oligo that maps more than once to the genome.  This "dupFree" 
 * 		bpmap text file is written to disk.
 * 
 * Recommend rewriting this to use first hashMap to eliminate dups and skip the second step.
 */
public class TPMapDuplicateFilter {
	// fields
	private HashMap hashMap;
	private TPMapFeature[][] duplicateFeatures;
	private int lineNumber; //3153283 number of lines in bpmap file, might get reset by one of the methods
	private File bpmapFile;
	private File saveDirectory;
	private String newFileName;
	
	/**For running this alone.*/	
	public TPMapDuplicateFilter(String fileName){
		bpmapFile = new File(fileName);
		saveDirectory = new File(bpmapFile.getParent());
		newFileName = bpmapFile.getName();
		filterTPMap();
	}
	
	/**For running in combination with TPMapProcessor*/
	public TPMapDuplicateFilter (File bpmapFile, File saveDirectory, String newFileName){
		this.bpmapFile = bpmapFile;
		this.saveDirectory = saveDirectory;
		this.newFileName = newFileName;
		lineNumber = (int)IO.countNumberOfLines(bpmapFile);
		System.out.println("\tNumber of lines in bpmap file...: "+lineNumber);
		filterTPMap();
	}
	
	public void filterTPMap(){	
		System.out.println("\tLoading hash with bpmap...");
		//load hash
		loadBPHashMap (lineNumber, bpmapFile);
			
		//run thru hash counting repeats
		ArrayList dupFeats = new ArrayList(30000);
		int[] dups = new int[101];
		Set keys = hashMap.keySet();
		Iterator it = keys.iterator();
		ArrayList al;
		int size;
		Object obj;
		while (it.hasNext()){
			obj = it.next();
			al = (ArrayList)hashMap.get(obj);			
			size = al.size();
			if (size>=100)dups[100] +=1;
			else dups[size] +=1;
			//add to TPMapFeature[][]
			if (size!=1){
				dupFeats.add(al);
			}
		}
		System.out.println("\nA histogram of oligo frequence on the chip (occurances, # oligos)");
		//Histogram.printHistogram(dups);

		//convert ArrayList into a TPMapFeature[][] to use in averaging oligo intensities
		int numDupFeats = dupFeats.size();
		duplicateFeatures = new TPMapFeature[numDupFeats][];
		ArrayList alFeats;
		for (int i=0; i< numDupFeats; i++){
			alFeats = (ArrayList)dupFeats.get(i);		
			TPMapFeature[] feats = new TPMapFeature[alFeats.size()];
			alFeats.toArray(feats);
			duplicateFeatures[i] = feats;
		}
		dupFeats = null;
		
		//save DupFeatures
		IO.saveObject(new File (saveDirectory, newFileName+"Dups"), duplicateFeatures);
		System.out.println("\n\tSaving DuplicateFeature[][] file...");
		duplicateFeatures = null;
		
		System.out.println ("\tRemoving multiple listings...");
		
		//make new bpmap file
		removeAdjacentDuplicates();	
		File adjDupRes = new File (saveDirectory, newFileName+"Dup");
		
		System.out.println ("\tRemoving oligos that map to multiple locations...");
		removeDuplicateMappingOligos(adjDupRes);
		adjDupRes.delete();		
	}
	
	/**Removes any oligo that is found more than once in the bpmap file.  Run
	 * removeAdjacentDuplicates() first!*/
	public void removeDuplicateMappingOligos(File bpmap){
		//build a Linked hash map of the bpmap file
		LinkedHashMap map = new LinkedHashMap(lineNumber);
		//read in line by line
		String line;
		String[] tokens;
		try {
			BufferedReader in = new BufferedReader(new FileReader(bpmap));
			while ((line = in.readLine()) !=null) {
				line = line.trim();              
				tokens = line.split("\\s+");
				if (tokens.length == 8){
					//check if sequence is present in hash
					if (map.containsKey(tokens[0])){
						//duplicate so set value as null
						map.put(tokens[0], null);
					}
					else map.put(tokens[0],line);
				}
			 }
			 in.close();
			 PrintWriter out = new PrintWriter(new FileWriter(bpmap.toString()+"Free"));
			 //iterate thru hash saving non null value lines
			 int numDups=0;
			 int numLines =0;
			 Iterator it = map.keySet().iterator();
			 while (it.hasNext()){
				Object obj = map.get(it.next());
				if (obj!=null){
					numLines++;
					out.println(obj);
				}
				else numDups++;
			 }
			 out.close();
			 //reset global number of lines field
			lineNumber = numLines;
			 System.out.println("\tSaving DupFree bpmap: "+bpmap.toString()+"Free");
			 System.out.println ("\t"+numDups+ " Duplicates were removed, "+numLines+" unique oligos were saved.\n");
		 }
		 catch (IOException e) {
			 e.printStackTrace();
		 }	
	}	
	/**Throws out data lines that are immediately adjacent and have the same sequence.*/
	public void removeAdjacentDuplicates(){
		//read in line by line
		String line;
		String[] oldTokens = null;
		String oldLine;
		String[] tokens;
		lineNumber =0;
		try {
			File outFile = new File (saveDirectory, newFileName+"Dup");
			 PrintWriter out = new PrintWriter(new FileWriter(outFile));
			 BufferedReader in = new BufferedReader(new FileReader(bpmapFile));
			
			//load old params
			while ((oldLine = in.readLine()) !=null) {
				oldLine = oldLine.trim();              
				oldTokens = oldLine.split("\\s+");
				if (oldTokens.length == 8){
					lineNumber++;
					break;
				}
			}
			while ((line = in.readLine()) !=null) {
				line = line.trim();              
				tokens = line.split("\\s+");
				
				if (tokens.length == 8){
					lineNumber++;
					//check if they are the same sequence
					if (tokens[0].equals(oldTokens[0])==false){
						//print old tokens
						out.println(oldLine);
						oldLine = line;
						oldTokens = tokens;
					}
				}
			 }
			 //print last line, close IO
			 out.println(oldLine);
			 out.close();
			 in.close();
		 }
		 catch (IOException e) {
			 e.printStackTrace();
		 }
	}
	
	public void loadBPHashMap(int initialValue, File file){
		//initialize hash
		hashMap = new HashMap(initialValue);
		
		//read in line by line
		String line;
		String[] tokens;
		lineNumber =0;
		TPMapFeature test;
		int num;
		ArrayList al;
		boolean add;
		try {
			 BufferedReader in = new BufferedReader(new FileReader(file));
			 while ((line = in.readLine()) !=null) {
				line = line.trim();              
				tokens = line.split("\\s+");
				
				if (tokens.length == 8){
					lineNumber++;
					//check if sequence is in hashMap
					if (hashMap.containsKey(tokens[0])){
						add = true;
						test = new TPMapFeature(tokens);
						al = (ArrayList)hashMap.get(tokens[0]);
						num = al.size();
						for (int i=0; i<num; i++){							
							if (matches(test, (TPMapFeature)al.get(i))){
			 					add = false;
								break;
							}
						}
						//add new TPMapCoordinate to ArrayList
						if (add) {
							al.add(test);
						} 
					}
					//if not, create new entry
					else {
						ArrayList al2 = new ArrayList(1);
						test = new TPMapFeature(tokens);
						al2.add(test);						
						hashMap.put(tokens[0], al2);
					}
				}
			 }
			 in.close();
		 }
		 catch (IOException e) {
			 e.printStackTrace();
		 }
	}	
	
	/**Only checks coordinates!*/
	public static boolean matches (TPMapFeature a, TPMapFeature b){
		if (a.getPMX() != b.getPMX()) return false;
		if (a.getPMY() != b.getPMY()) return false;
		if (a.getMMX() != b.getMMX()) return false;
		if (a.getMMY() != b.getMMY()) return false;
		return true;
	}
	public static void main(String args[]) {
		new TPMapDuplicateFilter (args[0]);
	}
	public int getLineNumber() {
		return lineNumber;
	}

}