package trans.tpmap;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;
import util.gen.*;

import util.gen.IO;

/**
 * App for processing the tpmap file for use by downstream TiMAT applications. 
 * 
 * @author nix
 */
public class TPMapProcessor {
	
	//fields
	private File tpmapFile;
	private File saveDirectory;
	private int windowSize = 675;
	private MapFeature[] features;
	private int[][] controlIndexes;
	private int minNumFeatures = 1;   //setting to 1, otherwise funky probe distributions can occur.
	private WindowMaker windowMaker;
	
	//constructor
	public TPMapProcessor(String[] args){
		processArgs(args);
		windowMaker = new WindowMaker(windowSize, minNumFeatures);
		
		//make a save directory, will over write any existing TPMapFiles dir
		String truncatedName = Misc.removeExtension(tpmapFile.getName());
		saveDirectory = new File(tpmapFile.getParent(), truncatedName+"_"+minNumFeatures+"oligos_"+windowSize+"bp_TPMapFiles");
		if (saveDirectory.exists()) IO.deleteDirectory(saveDirectory);
		saveDirectory.mkdir();
		
		//make MapFeature[] 
		convertTPMapToFeatures(tpmapFile);
		
		//filter MapFeature[] by windows
		if (minNumFeatures !=1) filterMapFeatures();
		
		//save MapFeature[]
		IO.saveObject(new File(saveDirectory, "tpmap.fa"), features);
		
		//save control Indexes, int[oligo #][indexes in MapFeature[], may be more than one if duplicates exist]
		if (controlIndexes != null){
			IO.saveObject(new File(saveDirectory, "tpmap.controlIndexes"), controlIndexes);
		}
		
		//fire MapSplitter
		System.out.println("Launching MapSplitter:");
		new MapSplitter(features, saveDirectory);
		
		//fire WindowMaker
		System.out.println("\nLaunching WindowMaker:");
		String[] chrFiles = WindowMaker.extractChrFiles(saveDirectory);
		if (chrFiles.length==0){
			System.out.println("\nError: No 'chr' files were found. Be sure your chromosome descriptions in the tpmap file begin " +
			"with 'chr' (e.g. chr1, chr2, chrX, chr2L; not 1, 2, X).\n");
			System.exit(0);
		}
		windowMaker.makeWindowsFromFiles(chrFiles);
		
		System.out.println("\nDone!\n");	
	}
	public void filterMapFeatures(){
		System.out.println("Filtering MapFeature[]:");
		ArrayList filteredFeatures = new ArrayList();
		ArrayList chromStarts = new ArrayList();
		ArrayList chromFeatures = new ArrayList();
		String currentChrom = features[0].chromosome;
		chromStarts.add(new Integer(features[0].start));
		chromFeatures.add(features[0]);
		//collect chromosome specific features and their starts
		for (int i=1; i< features.length; i++){
			String testChrom = features[i].chromosome;
			if (testChrom.equals(currentChrom)) {
				chromStarts.add(new Integer(features[i].start));
				chromFeatures.add(features[i]);
			}
			else {
				//filter chrom features
				//watch out for single oligo hits on a chromosome
				if (chromStarts.size() > 1) {
					ArrayList chromFilteredFeatures = filterChromosomeSpecificFeatures(chromStarts, chromFeatures);
					filteredFeatures.addAll(chromFilteredFeatures);
				}
				//reset to new chrom
				currentChrom = testChrom;
				chromStarts.clear();
				chromFeatures.clear();
				chromStarts.add(new Integer(features[i].start));
				chromFeatures.add(features[i]);
			}
		}
		//add last
		ArrayList chromFilteredFeatures = filterChromosomeSpecificFeatures(chromStarts, chromFeatures);
		filteredFeatures.addAll(chromFilteredFeatures);
		
		//convert to MapFeature[]
		System.out.println("\t"+(features.length- filteredFeatures.size()) + " oligo features were removed for falling outside a window.\n");
		features = new MapFeature[filteredFeatures.size()];
		filteredFeatures.toArray(features);
	}
	
	public ArrayList filterChromosomeSpecificFeatures (ArrayList startsAL, ArrayList chromMapFeaturesAL){
		int numFeatures = chromMapFeaturesAL.size();
		//convert starts to int[]
		int[] starts = Num.arrayListOfIntegerToInts(startsAL);
		int[][] windows = windowMaker.makeWindows(starts);
		boolean[] saveIndexes = new boolean[numFeatures];
		//for each window flip boolean it it should be saved, windows are overlapping
		for (int i=0; i< windows.length; i++){
			int[] startStop = windows[i];
			for (int j=startStop[0]; j<= startStop[1]; j++){
				saveIndexes[j] = true;
			}
		}
		//filter
		ArrayList filtered = new ArrayList();
		for (int i=0; i< numFeatures; i++) if (saveIndexes[i]) filtered.add(chromMapFeaturesAL.get(i));
		return filtered;
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
					case 'f': tpmapFile = new File(args[i+1]); i++; break;
					case 'w': windowSize=Integer.parseInt(args[i+1]); i++; break;
					case 'n': minNumFeatures=Integer.parseInt(args[i+1]); i++; break;
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
			System.out.println("\nCannot find your tpmap file!\n");
			System.exit(0);
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            TPMap Processor: Aug 2006                             **\n" +
				"**************************************************************************************\n" +
				"\nThe TPMapProcessor takes a tpmap file, creates sets of oligos (Windows) to be used\n" +
				"in statistical analysis, and splits the tpmap file into separate chromosomes. The\n" +
				"entire output folder is needed by several other TiMAT2 applications. Save it!\n\n" +
				
				"Parameters:\n"+
				"-w The maximum length for a Window in bp, default 675.\n" +
				"-f Full path file text for the text tpmap file.\n" +
				"-n Minimum number of oligos, defaults to 1.  Don't change if using the q-value\n" +
				"      estimation in ScanChip! Minimal oligo requirements can be set later.\n\n" +
				
				"Example: java -Xmx1000M -jar pathTo/T2/Apps/TPMapProcessor -f /affy/tpmap.txt -w 500 \n\n" +
				
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length <2){
			printDocs();
			System.exit(0);
		}
		new TPMapProcessor(args);
	}
	
	/**Converts a text tpmap to a MapFeature[] and tracks duplicates*/
	public void convertTPMapToFeatures(File tpmapFile){
		features = null;
		try {
			ArrayList al = new ArrayList(1000000);
			//read in and make array
			System.out.println("\nMaking MapFeature[]:");
			String line;
			BufferedReader in = new BufferedReader(new FileReader(tpmapFile));
			HashMap hash = new HashMap(2000);
			int index=0;
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() !=0 && line.startsWith("#") == false) {
					MapFeature mf = new MapFeature(line);
					al.add(mf);
					//is it a control probe?
					if (mf.chromosome.equals(MummerMapper.controlChromosomeName)){
						//get sequence
						String seq = line.split("\\s+")[0];
						//seq is already in hash, add to ArrayList of index numbers, it's a duplicate
						if (hash.containsKey(seq)){
							ArrayList indexAL = (ArrayList)hash.get(seq);
							indexAL.add(new Integer(index));
						}
						//seq is not in hash make new entry
						else{
							ArrayList indexAL = new ArrayList();
							indexAL.add(new Integer(index));
							hash.put(seq, indexAL);
						}
					}
					index++;
				}
			}
			in.close();
			//convert to MapFeature[]
			features = new MapFeature[al.size()];
			al.toArray(features);
			System.out.println("\t"+features.length+" Oligo features found.");
			
			//if applicable convert Hash to int[][], control oligos may not be mapped.
			if (hash.size()!=0){
				int total = 0;
				int counter = 0;
				controlIndexes = new int[hash.size()][];
				Iterator it = hash.keySet().iterator();
				while (it.hasNext()){
					String seq = (String)it.next();
					ArrayList indexDupsAL = (ArrayList)hash.get(seq);
					//ArrayList indexDupsAL = (ArrayList)hash.get(it.next());
					int[] indexDups = Num.arrayListOfIntegerToInts(indexDupsAL);
					controlIndexes[counter++] = indexDups;
					total += indexDups.length;					
				}
				System.out.println("\t"+hash.size()+" Unique and "+total+" Total control oligos found.\n");
			}
			else System.out.println("\tNo control oligos found.\n");
			
		} catch (Exception e) {e.printStackTrace();}
	}
}

