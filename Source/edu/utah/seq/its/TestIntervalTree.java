package edu.utah.seq.its;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;


/**Test app to work with Jannovars interval tree search classes.  Nix modified to use interbase coordinates.*/
public class TestIntervalTree {

	public static void main(String[] args) {
		
		//make set of intervals
		//List<Interval<Integer>> intervals = getRegions();
		//List<Interval<String>> intervals = getIntervalListString();
		//make Integer Tree disabling searchingForNeighbors, this can be set to true later on if desired.
		//IntervalTree<Integer> tree = new IntervalTree<Integer>(intervals, false);
		//IntervalTree<String> tree = new IntervalTree<String>(intervals, false);
		//List<String> qy = tree.search(0, 3);
		//List<Integer> qy = tree.search(0, 3);
		//System.out.println(qy);
		//search it for 20K times
		/*for (int x=0; x<10; x++){
			long start = System.currentTimeMillis();
			for (int i=0; i< 20000; i++){
				List<Integer> qy = tree.search(i, i+10);
			}
			System.out.println(x+ "\t"+ (System.currentTimeMillis() - start));
		}*/
		
		//make 100 chIPSeq datasets
		new TestIntervalTree();
		
		
	}
	
	public TestIntervalTree (){
		
		//fetch simulated datasets in IntervalTrees
		HashMap<Integer, IntervalTree<Peak>>[] chIPData = getChIPPeaks();
		
		//fetch a simulated test dataset
		HashMap<Integer, Peak[]> testPeaks = getPeaks();
		
		//intersect
		intersectPeaks (testPeaks, chIPData);
		
	}
	
	private void intersectPeaks(HashMap<Integer, Peak[]> testPeaks, HashMap<Integer, IntervalTree<Peak>>[] chIPData) {
		//for each test dataset
		long start = System.currentTimeMillis();
		for (int i=0; i< chIPData.length; i++){
			long begin = System.currentTimeMillis();
			HashMap<Integer, IntervalTree<Peak>> chIPDataset = chIPData[i];
			ArrayList<Peak> peaks = search(testPeaks, chIPDataset);
			System.out.println (peaks.size()+ "\tintersect\t"+ (System.currentTimeMillis()-begin));
		}
		long stop = System.currentTimeMillis();
		System.out.println("Total time for search "+ (stop-start));
		
	}

	public class Peak{
		int start;
		int stop;
		float transFDR;
		float log2Rto;
		public Peak(int start, int stop, float transFDR, float log2Rto) {
			this.start = start;
			this.stop = stop;
			this.transFDR = transFDR;
			this.log2Rto = log2Rto;
		}
		public int getStart() {
			return start;
		}
		public int getStop() {
			return stop;
		}
		public float getTransFDR() {
			return transFDR;
		}
		public float getLog2Rto() {
			return log2Rto;
		}
	}
	
	public HashMap<Integer, IntervalTree<Peak>>[] getChIPPeaks(){
		HashMap<Integer, IntervalTree<Peak>>[] datasets = new HashMap[100];
		for (int i=0; i< 100; i++) {
			datasets[i] = getIntervalTreePeaks();
		}
		return datasets;
	}
	
	public HashMap<Integer, Peak[]> getPeaks(){
		HashMap<Integer, Peak[]> chromIdTrees = new HashMap<Integer, Peak[]>();
		Random rand = new Random();
		//for 23 chroms
		int numPeaks = 0;
		for (int i=0; i< 23; i++){
			//for 100+ peaks in each chrom
			int num = 100 + rand.nextInt(1000);
			Peak[] p = new Peak[num];
			for (int j=0; j< num; j++){
				int start = rand.nextInt(100000);
				int stop = start + 100 + rand.nextInt(500);
				p[j] = new Peak(start, stop, 20f, 2.3f);
			}
			numPeaks+= num;
			chromIdTrees.put(new Integer(i), p);
		}
		System.out.println(numPeaks+ "\tMock peaks for intersection");
		return chromIdTrees;
	}
	
	public HashMap<Integer, IntervalTree<Peak>> getIntervalTreePeaks(){
		HashMap<Integer, IntervalTree<Peak>> chromIdTrees = new HashMap<Integer, IntervalTree<Peak>>();
		Random rand = new Random();
		//for 23 chroms
		int numPeaks = 0;
		for (int i=0; i< 23; i++){
			//for 100+ peaks in each chrom
			int num = 100 + rand.nextInt(1000);
			List<Interval<Peak>> intervals = new ArrayList<Interval<Peak>>(num);
			for (int j=0; j< num; j++){
				int start = rand.nextInt(100000);
				int stop = start + 100 + rand.nextInt(500);
				Peak p = new Peak(start, stop, 20f, 2.3f);
				intervals.add( new Interval<Peak>(start, stop, p) );
			}
			numPeaks+= intervals.size();
			IntervalTree<Peak> tree = new IntervalTree<Peak>(intervals, false);
			chromIdTrees.put(new Integer(i), tree);
		}
		System.out.println(numPeaks+ "\tMock peaks");
		return chromIdTrees;
	}
	
	/**Takes a HashMap of chromosomeIDs and their associated Peaks and searches them against the chromosomeID matched IntervalTrees.
	 * @return ArrayList of Peaks from the dataSet that intersect.*/
	public static ArrayList<Peak> search(HashMap<Integer, Peak[]> chromIdQuery, HashMap<Integer, IntervalTree<Peak>> chromIdDataSet){
		ArrayList<Peak> dataSetIds = new ArrayList<Peak>();
		//for each chromosome ID
		for (Integer chrId: chromIdQuery.keySet()){
			//does the chrom exist in the dataSet?
			IntervalTree<Peak> tree = chromIdDataSet.get(chrId);
			if (tree == null) continue;
			//for each queryRegion
			Peak[] queryRegions = chromIdQuery.get(chrId);
			//System.out.println("Num p in chrom "+queryRegions.length);
			//int counter = 0;
			for (Peak r: queryRegions){
				ArrayList<Peak> intersectingPeaks = tree.search(r.getStart(), r.getStop());
				//System.out.println("\tInts "+intersectingPeaks.size());
				dataSetIds.addAll(intersectingPeaks);
				//if (counter++ == 10) System.exit(0);
			}
		}
		return dataSetIds;
	}


	public static List<Interval<Integer>> getRegions(){
		List<Interval<Integer>> intervals = new ArrayList<Interval<Integer>>(20000);
		for (int i=0; i< 20000; i++){
			Integer x = new Integer(i);
			intervals.add( new Interval<Integer>(i, i+100, x) );
		}
		return intervals;
	}
	
	public static List<Interval<String>> getIntervalListString() {
	 	List<Interval<String>> ilist = new ArrayList<Interval<String>>();
	 	
	 	ilist.add(new Interval<String>(0,8,"a"));
	 	ilist.add(new Interval<String>(2,3,"b"));
	 	ilist.add(new Interval<String>(4,14,"c"));
	 	ilist.add(new Interval<String>(5,7,"d"));
	 	ilist.add(new Interval<String>(7,12,"e"));
	 	ilist.add(new Interval<String>(9,22,"f"));
	 	ilist.add(new Interval<String>(16,19,"g"));
	 	ilist.add(new Interval<String>(17,20,"h"));
	 	ilist.add(new Interval<String>(29,33,"i"));
	 	ilist.add(new Interval<String>(30,32,"j"));
	 
	 	return ilist;
	 
	}
}
