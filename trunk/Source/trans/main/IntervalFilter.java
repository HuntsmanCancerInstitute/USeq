package trans.main;
import java.io.*;
import java.util.regex.*;
import java.util.*;
import trans.misc.*;
import trans.anno.*;
import util.gen.*;


/**
 * Sorts an Interval[] into two groups, those that pass or fail a variety of filters.
 */ 
public class IntervalFilter {
	//fields
	private File[] files;
	//which filter to apply, boolean flags
	private boolean windowFilter = false;
	private boolean subWindowFilter = false;
	private boolean stndDev = false;
	private boolean deleteSpecifics = false;
	private boolean deleteIntersecting = false;
	private double subWindowValue = 0;
	private double windowValue = 0;
	private double stndDevValue = 0;			//really coefficient of the variation (stdn/mean)
	private File specificsToDelete = null;
	private File[] intersectingRegionFilesToFlag = null;
	private int scoreIndex = 1;
	private double minFractionIntersection = 0.25;
	
	private Interval[] intervals;
	private ArrayList bad = new ArrayList();
	private ArrayList good = new ArrayList();
	
	public IntervalFilter (String[] args){
		processArgs(args);
		//for each file attempt to get Interval[]
		int numFiles = files.length;
		for (int i=0; i<numFiles; i++){
			try{
				System.out.println("\nProcessing "+files[i]);
				
				//initialize interval groups
				intervals = (Interval[])IO.fetchObject(files[i]);
				good = Misc.objectArrayToArrayList(intervals);
				bad = new ArrayList();
				
				//go thru each boolean and filter
				if (deleteSpecifics) deleteSpecificIntervals();
				if (deleteIntersecting) deleteIntersectingIntervals();
				if (windowFilter) filterWindowScore();
				if (subWindowFilter) filterSubWindow();
				if (stndDev) filterByStndDev();
				
				//print final
				System.out.println("\n"+good.size() +" Good Intervals");
				System.out.println(bad.size() +" Bad Intervals");
				
				//save
				if (bad.size()!=0){
					if (deleteIntersecting){
						System.out.println("Saving intersected Interval[]...");
						IO.saveObject(new File(files[i].getCanonicalPath()+"Intr"), intervals);
					}
					System.out.println("Saving filtered Interval[]s:");
					intervals = new Interval[good.size()];
					good.toArray(intervals);
					IO.saveObject(new File(files[i].getCanonicalPath()+"Good"), intervals);
					System.out.println("\t"+files[i].getCanonicalPath()+"Good");
					intervals = new Interval[bad.size()];
					bad.toArray(intervals);
					IO.saveObject(new File(files[i].getCanonicalPath()+"Bad"), intervals);
					System.out.println("\t"+files[i].getCanonicalPath()+"Bad\n");
				}
				//save flagged if so filtered
				else if (intersectingRegionFilesToFlag != null && bad.size()==0) {
					IO.saveObject(new File(files[i].getCanonicalPath()+"Good"), intervals);
				}
				else {
					System.out.println("No bad intervals, nothing to save.\n");
				}
				
				
			}catch (Exception e){
				System.out.println("Skipping "+files[i]);
				e.printStackTrace();
			}
		}
	}
	

	
	public void deleteSpecificIntervals(){
		//load list to delete
		GenomicRegion[] toDelete = GenomicRegion.parseRegions(specificsToDelete);
		Interval interval;
		int counter = 0;
		int numToDelete = toDelete.length;
		for (int i=0; i<good.size(); i++){
			interval = (Interval)good.get(i);
			//run thru list toDelete looking for a match, assumes chrom-start are unique
			for (int j=0; j<numToDelete; j++){
				if ( (interval.getChromosome().equals(toDelete[j].getChromosome())) && (interval.getStart1stOligo() == toDelete[j].getStart()) ){
					bad.add(interval);
					good.remove(i);
					i--;
					counter++;
					break;
				}
			}
		}
		System.out.println(counter + " additional intervals removed by filtering for particular bad apples.");
	}
	
	public void deleteIntersectingIntervals(){
		//load regions to use in intersection, sort, split by chromosome
		HashMap[] flagRegions = new HashMap[intersectingRegionFilesToFlag.length];
		String[] fileNames = new String[flagRegions.length];
		RegionComparator comp= new RegionComparator();
		for (int i=0; i< intersectingRegionFilesToFlag.length; i++){
			GenomicRegion[] toSort = GenomicRegion.loadWriteBinaryRegions(intersectingRegionFilesToFlag[i]); 
			Arrays.sort(toSort, comp);
			flagRegions[i]= GenomicRegion.splitByChromosome(toSort);
			fileNames[i] = intersectingRegionFilesToFlag[i].getName();
		}
		if (fileNames.length>1) fileNames = Misc.trimCommon(fileNames);
		else fileNames[0] = Misc.removeExtension(fileNames[0]);
		
		Interval interval;
		int counter = 0;
		
		//break good Intervals list by chromosome
		HashMap chromAL = Util.splitIntervalArrayListByChromosome(good);
		good.clear();
		
		//for each interval set split by chromosome
		Iterator it = chromAL.keySet().iterator();
		while (it.hasNext()){
			
			//get ArrayList of chromosome specific Intervals
			ArrayList chromSpecificIntervals = (ArrayList) chromAL.get(it.next());
			String currentChromosome = ((Interval)chromSpecificIntervals.get(0)).getChromosome();
			//System.out.println("CurrChrom "+currentChromosome +"  NumIntervals "+chromSpecificIntervals.size());		
			
			//make GenomicRegion[][] of chromosome specific regions
			GenomicRegion[][] chromSpecificRegions = new GenomicRegion[flagRegions.length][];
			for (int x=0; x< flagRegions.length; x++){
				chromSpecificRegions[x] = (GenomicRegion[]) flagRegions[x].get(currentChromosome); //may return null!
				//System.out.println(x+" NumRegions "+chromSpecificRegions[x].length+" "+chromSpecificRegions[x][0].getChromosome());
			}
			
			//for each Interval
			for (int i=0; i<chromSpecificIntervals.size(); i++){
				interval = (Interval)chromSpecificIntervals.get(i);
				double intervalLength = interval.getStartLastOligo()+interval.getSizeOfOligoMinusOne()-interval.getStart1stOligo() +1;
				double[] fractionIntersect = new double[flagRegions.length];
				boolean badInterval = false;
				
				//for each set of regions
				for (int k=0; k< chromSpecificRegions.length; k++){
					int bpIntersection = 0;
					//any regions?
					if (chromSpecificRegions[k] == null) continue;
					
					//run thru chromSpecificRegions counting total intesection, these are sorted so break after missing after a hit
					boolean hit = false;
					for (int j=0; j<chromSpecificRegions[k].length; j++){
						int numBpInt = chromSpecificRegions[k][j].overlap(interval);
						//overlap
						if (numBpInt > 0 ){
							//System.out.println("Overlap: "+numBpInt);
							//System.out.println("\t Int: "+interval.toString()+ "  Length "+intervalLength);
							//System.out.println("\t Reg: "+chromSpecificRegions[k][j].toString());
							if (numBpInt == intervalLength){
								bpIntersection = numBpInt;
								break;
							}
							else {
								bpIntersection += numBpInt;
								hit = true;
							}
						}
						else if (hit) break;
					}
					//calculate fraction of intersection
					fractionIntersect[k] = ((double)bpIntersection)/intervalLength;
					if (fractionIntersect[k] >= minFractionIntersection) badInterval = true;
				}
				
				//load results
				interval.setRegionNames(fileNames);
				interval.setFractionIntersections(fractionIntersect);
				
				//flag it?
				if(badInterval){
					bad.add(interval);
					chromSpecificIntervals.remove(i);
					i--;
					counter++;
				}
				
			}
			good.addAll(chromSpecificIntervals);
		}
		System.out.println(counter + " intervals removed by filtering for intersection with particular regions.");
	}
	
	public void filterWindowScore(){
		Interval interval;
		int counter = 0;
		for (int i=0; i<good.size(); i++){
			interval = (Interval)good.get(i);
			if (windowValue > interval.getBestWindow().getScores()[scoreIndex]){
				bad.add(interval);
				good.remove(i);
				i--;
				counter++;
			}
		}
		System.out.println(counter + " additional intervals removed by filtering for a window score of less than "+windowValue);
	}
	
	
	
	/**Filters by the median ratio of the sub window*/
	public void filterSubWindow(){
		Interval interval;
		int counter = 0;
		for (int i=0; i<good.size(); i++){
			interval = (Interval)good.get(i);		
			if (interval.getBestSubWindow() == null){
				System.out.println("\tNo sub window!, deleting");
				bad.add(interval);
				good.remove(i);
				i--;
				counter++;
			}
			else if (subWindowValue > interval.getBestSubWindow().getMedianRatio()){
				bad.add(interval);
				good.remove(i);
				i--;
				counter++;
			}
		}
		System.out.println(counter + " additional intervals removed by filtering for a fold difference in the sub window of less than "+subWindowValue);
	}
	
	//not really standard deviation but coefficient of the variation, std/mean
	public void filterByStndDev(){
		Interval interval;
		Oligo[] oligos;
		int numOligos;
		float treatments[][];
		float controls[][];
		double stndT;
		double stndC;
		int counter = 0;
		double meanT;
		double meanC;
		float[] collapsedT;
		float[] collapsedC;
		
		for (int i=0; i<good.size(); i++){
			interval = (Interval)good.get(i);
			
			oligos = interval.extractSubRegionOligos(interval.getBestWindow().getStart1stOligo(),
					interval.getBestWindow().getStartLastOligo());	
			numOligos = oligos.length;
			treatments = new float[numOligos][];
			controls = new float[numOligos][];
			for (int j=0; j<numOligos; j++){
				treatments[j] = oligos[j].getTreatmentIntensities(interval.getNumberTreatmentIntensities());
				controls[j] = oligos[j].getControlIntensities(interval.getNumberControlIntensities());
			}
			collapsedT = Num.collapseFloatArray(treatments);
			collapsedC = Num.collapseFloatArray(controls);
			meanT = Num.mean(collapsedT);
			meanC = Num.mean(collapsedC);
			
			stndT = Num.standardDeviation(collapsedT, meanT)/meanT;
			stndC = Num.standardDeviation(collapsedC, meanC)/meanC;
			if (stndDevValue < stndT || stndDevValue < stndC){
				bad.add(interval);
				good.remove(i);
				i--;
				counter++;
			}
		}
		System.out.println(counter + " additional intervals removed by filtering for standard of deviations, for either the treatment values or control values, greater than "+stndDevValue);
	}
	
	
	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Interval Filter: April 2007                           **\n" +
				"**************************************************************************************\n" +
				"IF filters Interval[] arrays dividing and saving Intervals that pass and fail.\n"+
				"The reported number of cut intervals is sequential and may not reflect the total\n"+
				"number that would fail each particular test when using a single filter.\n\n" +
				
				"-i Score index to use in window filtering, see ScanChromosome or ScanChipNoPerm\n"+
				"-a Filter by the best window score based on the score index (minimum).\n"+
				"-b Filter by the best median ratio sub window (minimum).\n"+
				"-f Filter by the coefficient of the variation (stnd dev window oligos/ mean) for\n" +
				"      the intensities, either treatment or control (singleton,\n" +
				"      one chip driver) (maximum)\n"+
				"-e Flag intervals that intersect particular regions, enter full path names, comma\n" +
				"      delimited, containing tab delimited: chrom start stop.\n"+
				"-m Minimal fraction of intersection with particular regions for removal, defaults to\n"+
				"      0.25.  Measured as cumulative bases covered, relative to the interval.\n"+
				"-g Remove particular intervals, enter a full path text for a tab delimited text\n"+
				"      file containing rows with chromosome start stop.\n"+
				"-k Full path file text for the Interval[], if a directory is specified, all files\n" +
				"      within will be processed. (Required)\n\n" +

				"Example: java -Xmx256M -jar pathTo/T2/Apps/IntervalFilter -k /affy/res/ -a 1.5 -i 1\n" +
				"      -e /repeats/hg17_simpleSegDupRepMask.txt -m 0.5\n" +
				"\n" +
		"**************************************************************************************\n");
	}

	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File directory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		String intersectingFileNames = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': scoreIndex = Integer.parseInt(args[i+1]); i++; break;
					case 'k': directory = new File(args[i+1]); i++; break;
					case 'a': windowFilter = true; windowValue = Double.parseDouble(args[i+1]); i++; break;
					case 'b': subWindowFilter = true; subWindowValue = Double.parseDouble(args[i+1]); i++; break;
					case 'f': stndDev = true; stndDevValue = Double.parseDouble(args[i+1]); i++; break;
					case 'g': specificsToDelete = new File(args[i+1]); i++; deleteSpecifics = true; break;
					case 'e': intersectingFileNames = args[i+1]; i++; deleteIntersecting = true; break;
					case 'm': minFractionIntersection = Double.parseDouble(args[i+1]); i++; break;
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
		
		if (deleteSpecifics == true && (specificsToDelete == null || specificsToDelete.exists() == false) ){
			System.out.println("\nEnter a tab delimited text file (chromosome start stop) file containing specific intervals to remove.\n");
			System.exit(0);
		}

		if (intersectingFileNames != null){
			intersectingRegionFilesToFlag = IO.extractFiles(intersectingFileNames);
			if (intersectingRegionFilesToFlag == null || intersectingRegionFilesToFlag.length ==0) Misc.printExit("\nCannot find any files with which to intersect your intervals?! "+intersectingFileNames);
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
	
	public static void main(String[] args) {
		if (args.length <2){
			printDocs();
			System.exit(0);
		}
		new IntervalFilter(args);
	}
}
