package edu.utah.seq.analysis;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

import trans.roc.ParseSgrsForParticularRegions;
import trans.roc.PositiveComparator;
import util.bio.seq.Seq;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

import edu.utah.seq.analysis.DefinedRegionDifferentialSeq;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.analysis.multi.Condition;
import edu.utah.seq.analysis.multi.Replica;
import edu.utah.seq.useq.data.Region;
import trans.roc.Positive;

import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;


public class ScoreEnrichedRegions {
	//primitives
	private boolean verbose = true;
	private int maxAlignmentsDepth=50000;
	private int pseudocounts =10;
	private double fractionGCTolerance=0.1;
	private int numberRandom=1000;
	private double[] readCounts = new double[2];
	
	//complex
	private Condition[] conditions;
	private HashMap<String, File> chromGCFileMap;
	private HashMap<String, ArrayList<Positive>> regions;
	private HashMap<String, Region[]> interrogatedRegions;
	private ArrayList<String> chromOrder;
	
	//main data containers
	private ArrayList<GeneCountList> regionData;
	private ArrayList<ArrayList<GeneCountList>> randomData;
	HashSet<String> badRegions;
	
	//Files
	private File outputFile ;
	
	/**Simple class that holds a string and a index.  Similar to R's named vector */
	private class NameInteger {
		private String name;
		private Integer index;
		
		public NameInteger(String name, Integer index) {
			this.name = name;
			this.index = index;
		}
	}
	
	/**Fetches an old or makes a new Integer to represent the sam read name (e.g. fragment name)*/
	private class NameTracker {
		private final Pattern BAD_NAME = Pattern.compile("(.+)/[12]$");
		private int posIndex;
		private int negIndex;
		private HashMap<String,Integer> fragNameIndex;
		
		public NameTracker() {
			this.posIndex = 1;
			this.negIndex = -1;
			this.fragNameIndex = new HashMap<String, Integer>(10000);
		}
		
		public NameInteger fetchFragmentNameIndex(SAMRecord sam) {
			String samReadName = sam.getReadName();
			Integer index;
			Matcher mat = BAD_NAME.matcher(samReadName);
			if (mat.matches()) {
				samReadName = mat.group(1);
			}
			
			if (fragNameIndex.containsKey(samReadName)) {
				index = fragNameIndex.get(samReadName);
			} else {
				if (sam.getReadNegativeStrandFlag()) {
					index = new Integer(negIndex--);
				} else {
					index = new Integer(posIndex++);
				}
				fragNameIndex.put(samReadName, index);
			}
			return new NameInteger(samReadName,index);
		}
		
		public void removeIndex(String name) {
			fragNameIndex.remove(name);
		}
	}
	
	private class IRegions {
		private Region[] currentInterrogatedRegions;
		private int totalBPInterrogatedRegions;
		private int[] startsBPInterrogatedRegions;
		
		public IRegions(String chrom, HashMap<String,Region[]> regions) {
			currentInterrogatedRegions = regions.get(chrom);
			if (currentInterrogatedRegions == null) Misc.printExit("\nError: cannot find interrogated regions for "+ chrom);
			totalBPInterrogatedRegions = Region.totalBP(currentInterrogatedRegions);
			startsBPInterrogatedRegions = Region.startsInBases(currentInterrogatedRegions);
		}
		
	}
	
	private class GeneCountList {
		private int totalCounts;
		private HashMap<String,Integer> countMap;
		
		public GeneCountList() {
			countMap = new HashMap<String,Integer>();
			totalCounts = 0;
		}
		
		public void addGene(String name, Integer count) {
			totalCounts += count;
			countMap.put(name, count);
		}
		
		public void remove(ArrayList<String> badnames) {
			for (String name: badnames) {
				countMap.remove(name);
			}
			
		}
		
	}

	public ScoreEnrichedRegions(String[] args) {
		//process arguments
		processArgs(args);
		
		//initialize containers
		regionData = new ArrayList<GeneCountList>();
		randomData = new ArrayList<ArrayList<GeneCountList>>();
		badRegions = new HashSet<String>();
		
		//Grab gene counts from data
		System.out.println("Load Conditions");
		for (int i=0; i< conditions.length; i++) {
			System.out.println("\nCondition " + (i+1));
			for (Replica r: conditions[i].getReplicas()){
				System.out.println("Replica " + (i+1));
				loadReplica(r,i);
			}
		}
		
		//Generate the per region data
		System.out.println("Calculate per region data");
		calculatePerRegion();
		
		//Clean counts
		System.out.println("Clean counts");
		cleanupCounts();
		
		//calculate aggregate score
		System.out.println("Calculate Aggregate");
		calculateAggregate();
		
	}
	
	
	
	private void calculatePerRegion() {
		//make sure there is actually data
		if (regionData.size() == 0) {
			System.out.println("There is not data loaded in the regionData array");
			System.exit(1);
		}
		
		//Open file handle
		
		try {
			BufferedWriter br = new BufferedWriter(new FileWriter(outputFile));
			br.write("#Index\tChrom\tStart\tStop\tLength\tCondition1\tCondition2\tRegion_Log2Ratio\tAverage_Random_Log2Ratio\tPvalue_Greater\tPvalue_Less\n");
			
			//calculate globals
			double scalarCT = readCounts[1] / (double)readCounts[0];
			double scalarTC = readCounts[0] / (double)readCounts[1];
			
			//Organize the region data
			ArrayList<GeneCountList>[] oRegionData = organizeData(regionData);
			
			//Iterate over each region
			int regionIndex = 0;
			for (String chrom: chromOrder) {
				for (Positive region: regions.get(chrom)) {
					//regionName
					String regionName = String.valueOf(regionIndex);
					
					//Calculate Log-Ratio of the region
					double realA = 0;
					double realB = 0;
					for (GeneCountList gcl: oRegionData[0]) {
						realA += gcl.countMap.get(regionName);
					}
					
					for (GeneCountList gcl: oRegionData[1]) {
						realB += gcl.countMap.get(regionName);
					}
					
					double RealLogRatio = this.calculateLog2Ratio(realA, realB, scalarTC, scalarCT);
					
					//Check if region is valid
					if (badRegions.contains(String.valueOf(regionIndex))) {
					  br.write(String.format("%s\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%s\t%s\t%s\n",regionName,region.getChromosome(),region.getStart(),region.getStop(),
							  region.getLength(),realA,realB,RealLogRatio,"Not Enough Random Regions","NA","NA"));
					} else {
						//Transpose Random data
						ArrayList<GeneCountList> gcla = new ArrayList<GeneCountList>();
						for (int sampleIndex=0;sampleIndex<regionData.size();sampleIndex++) {
							GeneCountList gcl = new GeneCountList();
							for (int randomIndex=0;randomIndex<numberRandom;randomIndex++) {
								String randomName = String.valueOf(randomIndex);
								gcl.addGene(randomName,randomData.get(sampleIndex).get(randomIndex).countMap.get(regionName));
							}
							gcla.add(gcl);
						}
						
						//Organize it by condition
						ArrayList<GeneCountList>[] oRandomData = organizeData(gcla);
						
						//variables
						double greater = 0;
						double less = 0;
						double average = 0;
						
						//Generate log2ratios for random regions
						for (int randomIndex=0;randomIndex<numberRandom;randomIndex++) {
							String randomName = String.valueOf(randomIndex);
							double randA = 0;
							double randB = 0;
							
							for (GeneCountList gcl: oRandomData[0]) {
								randA += gcl.countMap.get(randomName);
							}
							
							for (GeneCountList gcl: oRandomData[1]) {
								randB += gcl.countMap.get(randomName);
							}
							
							double RandLogRatio = this.calculateLog2Ratio(randA, randB, scalarTC, scalarCT);
							
							average += RandLogRatio;
							
							if (RandLogRatio > RealLogRatio) {
								greater += 1;
							} else if (RandLogRatio < RealLogRatio) {
								less += 1;
							} else {
								greater += 1;
								less += 1;
							}
						}
						average = average /(double)numberRandom;
						greater = greater / (double)numberRandom;
						less = less / (double)numberRandom;
						
						br.write(String.format("%s\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",regionName,region.getChromosome(),region.getStart(),region.getStop(),
								  region.getLength(),realA,realB,RealLogRatio,average,greater,less));
					}
					regionIndex += 1; 
				}
			} 
			br.close();
		} catch (IOException ioex) {
			System.out.println("Could write to file: " + ioex.getMessage());
			ioex.printStackTrace();
			System.exit(1);
		}
		
		
		
		
	}
	
	private void calculateAggregate() {
		if (regionData.size() == 0) {
			Misc.printErrAndExit("No data was loaded. Please check the input files");
		}
		ArrayList<GeneCountList>[] rd = this.organizeData(regionData);
		//int totalNumberRegions = regionData.get(0).countMap.size();
		
		double realMed = getLog2Ratio(rd[0],rd[1]);
		double[] fakeMed = new double[numberRandom];
		
	
		for (int i=0; i<fakeMed.length; i++) {
			ArrayList<GeneCountList> subList = new ArrayList<GeneCountList>();
			for (int j=0; j<randomData.size();j++) {
				subList.add(randomData.get(j).get(i));
			}
			ArrayList<GeneCountList>[] fd = this.organizeData(subList);
			fakeMed[i] = getLog2Ratio(fd[0],fd[1]);
		}
		
		double numGreaterThan = 0;
		double numLessThan = 0;
		double average = 0;
		
		for (double ratio: fakeMed) {
			average += ratio;
			
			if (ratio > realMed) {
				numGreaterThan += 1;
			} else if (ratio < realMed) {
				numLessThan += 1;
			} else {
				numGreaterThan += 1;
				numLessThan += 1;
			}
		}
		
		double scalarCT = readCounts[1] / (double)readCounts[0];
		double scalarTC = readCounts[0] / (double)readCounts[1];
		System.out.println(String.format("\nScaling factor Condition1: %.3f\nScaling factor Condition2: %.3f\n",scalarCT,scalarTC));
		
		
		Arrays.sort(fakeMed);
		double test = Num.median(fakeMed);
		double med = average / ((double)numberRandom);
		double numGreaterThanP = numGreaterThan/ ((double)numberRandom);
		if (numGreaterThanP > 1) numGreaterThanP = 1;
		double numLessThanP = numLessThan/ ((double)numberRandom);
		if (numLessThanP > 1) numLessThanP = 1;
		System.out.println(String.format("\nMedian Enrichment Targeted Regions: %.3f\n" +
										 "Median Median Enrichment Random Regions %.3f\n" +
										 "Mean Median Enrichment Random Regions: %.3f.\n" +
										 "P-value targeted more enriched than random regions: %.3f (%d/%d).\n" + 
										 "P-value targeted more enriched than random regions: %.3f (%d/%d).\n",
				(double)realMed,(double)test,(double)med,(double)numGreaterThanP,
				(int)numGreaterThan,numberRandom,(double)numLessThanP,(int)numLessThan,numberRandom));

	}
	
	private boolean[] loadGCContent(String chrom) {
		File file = chromGCFileMap.get(chrom);
		if (file == null) {
			Misc.printExit("\nError: cannot find a gc binary file for "+ chrom);
		}
		
		boolean[] currentGCContent = (boolean[])IO.fetchObject(file);
		
		return currentGCContent;
	}
	
	private void cleanupCounts() {
		//any genes with too many reads that were excluded?
		
		if (badRegions.size() !=0) {
			System.err.println("\nWARNING: The following genes/ regions were excluded from the analysis due to one or more bps exceeding the maximum read coverage of "+maxAlignmentsDepth+" . If these are genes/ regions you wish to interrogate, increase the maximum read coverage threshold. " +
					"Realize this will likely require additional memory. Often, these are contaminants (e.g. rRNA) or super abundant transcripts that should be dropped from the analysis.\n"+badRegions.toString());
			
			//remove flagged genes from all replicas
			//for each condition
			
			for (GeneCountList gcl: regionData) {
				gcl.remove(new ArrayList<String>(badRegions));
			}
			
			for (ArrayList<GeneCountList> gcla: randomData) {
				for (GeneCountList gcl: gcla) {
					gcl.remove(new ArrayList<String>(badRegions));
				}
			}
		}
		
	}
	
	
	private ArrayList<GeneCountList>[] organizeData(ArrayList<GeneCountList> data) {
		ArrayList<GeneCountList>[] countData = new ArrayList[conditions.length];
		int index=0;
		
		//split by condition and replica
		for (int i=0;i<conditions.length;i++) {
			countData[i] = new ArrayList<GeneCountList>();
			for (int j=0;j<conditions[i].getReplicas().length;j++) {
				countData[i].add(data.get(index));
				index += 1;
			}
		}
		return countData;
	}
	
	
	private void loadReplica(Replica replica, int condNumber){
		//the idea here is to take a sorted bam file and add the reads to an ArrayList, basically every base should have a ArrayList of reads that overlap that base
		//one can then count the number of overlapping fragments for an exon or a collection of exons by hashing the ArrayList
		//assumes each read (first or second) has the same name
		
		/**
		 * Intialize Stuff
		 */
		
		//make reader
		SAMFileReader reader = new SAMFileReader(replica.getBamFile());
		
		//initialize region list
		GeneCountList regionGCL = new GeneCountList();
		
		//initialize replica gene list
		ArrayList<GeneCountList> randomGCLA = new ArrayList<GeneCountList>();
		for (int i=0; i<numberRandom;i++) {
			GeneCountList randomGCL = new GeneCountList();
			randomGCLA.add(randomGCL);
		}
		
		//Genome variables
		HashMap<String, Integer> chromLength = new HashMap<String, Integer>();
		List<SAMSequenceRecord> seqs =  reader.getFileHeader().getSequenceDictionary().getSequences();
		for (SAMSequenceRecord sr: seqs) {
			chromLength.put(sr.getSequenceName(), sr.getSequenceLength());
		}
		
		/**
		 * Read through the SAM file by chromosome
		 */
		int regionIndex = 0;
		
		for (String chrom: chromOrder) {
			if (!(chromLength.containsKey(chrom))) {
				System.out.println("Alignment file does not contain chromosome: " + chrom);
				System.exit(1);
			}
			
			//SAM-specific variables
			System.out.println("Working on: " + chrom);
			
			SAMRecordIterator iterator = reader.query(chrom,0,chromLength.get(chrom),false);
			SAMRecord sam = null;
			
			
			//Chromosome-specific variables
			ArrayList<Integer>[] bpNames = new ArrayList[chromLength.get(chrom)];
			boolean[] badBases = new boolean[chromLength.get(chrom)];
			
			NameTracker nt = new NameTracker();
			
			while (iterator.hasNext()) {
 				sam = iterator.next();
				
				//unaligned? too many hits?
				if (alignmentFails(sam)) continue;
				
				loadBlocks(sam,0,bpNames,badBases,nt);
				readCounts[condNumber] += 1;
			}
			
			iterator.close();
			
			//Collect Data for original region
			System.out.println("Loading count data into specified regions");
			loadGeneCounts(regionGCL,regions.get(chrom),bpNames,badBases,chrom);
			
			//Get Random Regions for chromsome
			boolean[] chromSpecificGC = loadGCContent(chrom);
			IRegions ir = new IRegions(chrom,interrogatedRegions);
			
			//Get random Regions
			System.out.println("Creating random regions");
			ArrayList<ArrayList<Positive>> randomRegions = new ArrayList<ArrayList<Positive>>();
	
			
			for (Positive region: regions.get(chrom)) {
				ArrayList<Positive> rr = makeRandom(region,chromSpecificGC,ir,regionIndex);
				regionIndex += 1;
				if (rr.size() < numberRandom) {
					System.out.println(String.format("Could not find enough random regions that match this location: %s:%d-%d",chrom,region.getStart(),region.getStop()));
					badRegions.add(String.valueOf(region.getIndex()));
				}
				randomRegions.add(rr);
			}
			
			
			//Find the largest row in the data
			int maxSize = Integer.MIN_VALUE;
			for (ArrayList<Positive> row: randomRegions) {
				if (row.size() > maxSize) {
					maxSize = row.size();
				}
			}
			
			System.out.println("Transposing regions");
			//Transpose Regions
			ArrayList<ArrayList<Positive>> realRandomRegions = new ArrayList<ArrayList<Positive>>();
			for(int i=0; i<maxSize; i++) {
				ArrayList<Positive> col = new ArrayList<Positive>();
				for (ArrayList<Positive> row: randomRegions) {
					try {
						col.add(row.get(i));
					} catch (IndexOutOfBoundsException iobex) {
						//ArrayLists might not all have the same size, so this is here to make sure
						//errors aren't thrown.
					}
				}
				realRandomRegions.add(col);
			}
			randomRegions.clear();
			
			System.out.println("Loading gene counts into random regions");
			//Load counts into random regions
			for (int i=0;i<numberRandom;i++) {
				loadGeneCounts(randomGCLA.get(i),realRandomRegions.get(i),bpNames,badBases,chrom);
			}
			
			
		}	
		
		reader.close();
		
		//add gene list to final data
		randomData.add(randomGCLA);
		regionData.add(regionGCL);

		//Clean up objects
		seqs = null;
		chromLength = null;
	}
	
	private void loadBlocks(SAMRecord sam, int chrStartBp, ArrayList<Integer>[] bpNames, boolean[] badBases, NameTracker nt){
		NameInteger nameIndex = null;
		boolean addIt;
		ArrayList<int[]> blocks = DefinedRegionDifferentialSeq.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
		//add name to each bp
		for (int[] b : blocks){
			int start = b[0] - chrStartBp;
			int stop = b[1] - chrStartBp;
			for (int i=start; i < stop; i++){
				//bad base?
				if (badBases[i]) continue;
				addIt = true;
				
				//never seen before?
				if (bpNames[i] == null) {
					bpNames[i] = new ArrayList<Integer>();
				}
				
				//old so check size
				else if (bpNames[i].size() == maxAlignmentsDepth) {
					badBases[i] = true;
					addIt = false;
					bpNames[i] = null;	
					if (nameIndex !=null) {
						nt.removeIndex(nameIndex.name);
					}
				}
				
				// add it?
				if (addIt) {
					if (nameIndex == null) {
						nameIndex = nt.fetchFragmentNameIndex(sam);
					}
					bpNames[i].add(nameIndex.index);
				}
			}
		}
	}
	
	private void loadGeneCounts(GeneCountList geneCountList, ArrayList<Positive> cRegion, ArrayList<Integer>[] bpNames, boolean[] badBases, String chromosome){
		
		HashSet<Integer> allReads = new HashSet<Integer>();
		
		for (int i=0; i< cRegion.size(); i++){
			Positive region = cRegion.get(i);
			Integer index = region.getIndex();
			
			//skip region if deemed bad
			if (badRegions.contains(region.getIndex())) {
				continue;
			}
			
			//Initialize Stuff
			allReads.clear();
	
			for (int y=region.getStart(); y< region.getStop(); y++){
				//bad base? if so then flag entire region
				if (badBases[y]) {
					badRegions.add(index.toString());
					allReads.clear();
					break;
				}
				
				if (bpNames[y] != null) {
					allReads.addAll(bpNames[y]);
				} 
			}
				
			int numCounts = allReads.size() + pseudocounts;
			
			geneCountList.addGene(String.valueOf(index), numCounts);
		
		}	
		
		//clean up
		allReads = null;

	}
	
	/**Calculates the gc content.*/
	private double calculateFractionGCContent(int start, int stop, boolean[] currentGCContent){
		double ave = 0;
		int realStop = stop +1;
		if (realStop >= currentGCContent.length){
			realStop = currentGCContent.length -1;
		}
		for (int i = start; i< realStop; i++){
			if (currentGCContent[i]) ave++;
		}
		double length = realStop - start;
		return ave/ length;
	}
	
	
	/**Takes a Positive and finds 1000 chrom, length, and gc matched random regions.*/
	private ArrayList<Positive> makeRandom(Positive region, boolean[] currentGCContent, IRegions ir, int regionIndex){
		//calculate gc content of real region and set min max
		double realGC = calculateFractionGCContent(region.getStart(), region.getStop(), currentGCContent);		
		int sizeRealRegionMinOne = region.getStop() - region.getStart();
		double minGC = realGC - fractionGCTolerance;
		double maxGC = realGC + fractionGCTolerance;

		
		//try 100,000 times to find a gc and min num matched random region
		Random randomGenerator = new Random();
		ArrayList<Positive> randomRegions = new ArrayList<Positive>();
		
		for (int x=0; x<numberRandom; x++){
			for (int y=0; y< 100000; y++){
				//randomly pick an interrogatedRegion unbiased by size
				//pick a base from total
				int randomBase = randomGenerator.nextInt(ir.totalBPInterrogatedRegions);
				//find it's index
				int index = Arrays.binarySearch(ir.startsBPInterrogatedRegions, randomBase);
				//pull the region
				int indexInterrogatedRegion;
				if (index >=0) indexInterrogatedRegion = index;
				else indexInterrogatedRegion = -index -1;
				Region interrogatedRegion = ir.currentInterrogatedRegions[indexInterrogatedRegion];
				//get it's size
				int lengthIntReg = interrogatedRegion.getLength();
				//check size to be sure it's big enough
				if (lengthIntReg < sizeRealRegionMinOne) continue;
				//generate a random start stop matching the length of the region
				int start = randomGenerator.nextInt(lengthIntReg) + interrogatedRegion.getStart();
				//check to see if start is < 0
				if (start < 0) continue;
				int stop = start + sizeRealRegionMinOne;
				//check to see if stop goes past the stop of the interrogated region
				if (stop >= interrogatedRegion.getStop()) continue;
				//check gc
				double testGC = calculateFractionGCContent(start, stop,currentGCContent);
				if (testGC < minGC || testGC > maxGC) continue;
				//check num scores
				Positive testRegion = new Positive(regionIndex,region.getChromosome(), start, stop);
				//load region with scores
			
			
				randomRegions.add(testRegion);
				break;
			}
		}
		
		return randomRegions;
		
	}
	
	/**Calculate the median log2ratio given arrays of GeneCountLists for two Conditions
	 */
	private double getLog2Ratio(ArrayList<GeneCountList> cond1, ArrayList<GeneCountList> cond2) {
		//Combined region counts per condition
		float[] countsR1 = new float[cond1.get(0).countMap.size()];
		float[] countsR2 = new float[cond2.get(0).countMap.size()];
		
		//Collapse the data
		for (GeneCountList rep: cond1) { //treatment
			ArrayList<Integer> counts = new ArrayList<Integer>(rep.countMap.values());
			for (int i=0; i<counts.size(); i++) {
				countsR1[i] += counts.get(i);
			}
		}
		
		for (GeneCountList rep: cond2) { //control
			ArrayList<Integer> counts = new ArrayList<Integer>(rep.countMap.values());
			for (int i=0; i<counts.size(); i++) {
				countsR2[i] += counts.get(i);
			}
		}

		//Caculate scalar values
		//double scalarCT = countsA2 / (double)countsA1;
		//double scalarTC = countsA1 / (double)countsA2;
		double scalarCT = readCounts[1] / (double)readCounts[0];
		double scalarTC = readCounts[0] / (double)readCounts[1];
		
		//Create log2Ratio Container
		double[] log2Ratio = new double[countsR1.length];
		
		
		//Calculate ratio
		for (int i=0;i<countsR1.length;i++) {
			log2Ratio[i] = calculateLog2Ratio(countsR1[i], countsR2[i], scalarTC, scalarCT);
		}
		
		//Calculate median ratio
		double median = -1;
		
		//Sort array
		Arrays.sort(log2Ratio);

		if (log2Ratio.length == 1) {
			median = log2Ratio[0];
		}
		else if (log2Ratio.length == 2) {
			median = Num.mean(log2Ratio);
		}
		else if (log2Ratio.length > 2) {
			median = (double)Num.median(log2Ratio);
		}
		
		return median;
	}
	
	
	/**Calculates a log2( (tSum+1)/(cSum+1) ) on linearly scaled tSum and cSum based on the total observations.*/
	private float calculateLog2Ratio( double tSum, double cSum, double scalarTC, double scalarCT){
		double t;
		double c;
		if (tSum !=0 ) {
			t = tSum * scalarCT;
			c = cSum;
		}
		else {
			c = cSum * scalarTC;
			t = tSum;
		}
		double ratio = (t+1)/(c+1);
		return (float)Num.log2(ratio);
	}
	
	private boolean alignmentFails(SAMRecord sam){
		//aligned?
		if (sam.getReadUnmappedFlag()) {
			return true;
		} else {
			return false;
		}
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ScoreEnrichedRegions(args);

	}
	
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		
		File regionFile=null;
		File interrFile=null;
		File gcGenomeDir=null;
		File[] conditionDirectories=null;
		
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': conditionDirectories = IO.extractOnlyDirectories(new File (args[++i])); break;
					case 'r': regionFile = new File(args[++i]); break;
					case 'i': interrFile = new File(args[++i]); break;
					case 'g': gcGenomeDir = new File(args[++i]); break;
					case 'x': maxAlignmentsDepth = Integer.parseInt(args[++i]); break;
					case 'f': pseudocounts = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					case 'o': outputFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//look for necessary bam files
		if (conditionDirectories == null || conditionDirectories.length == 0) Misc.printErrAndExit("\nError: cannot find any condition directories?\n");
		if (conditionDirectories.length != 2) Misc.printErrAndExit("\nError: must provide only two Conditions for analysis.\n");
		if (outputFile == null) Misc.printErrAndExit("\nError: no output file specified\n");

		//if everything is kosher, load up condition information
		conditions = new Condition[conditionDirectories.length];
		for (int i=0; i<conditionDirectories.length;i++) {
			File dir = conditionDirectories[i];
			File[] bamFiles = IO.extractFiles(dir, ".bam");
			OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, true);
			conditions[i] = new Condition(dir);
		}
		
		//look for bed file
		if (regionFile == null) {
			Misc.printErrAndExit("\nPlease enter a region of interest file in Bed format.\n");
		} else if (interrFile == null) {
			Misc.printErrAndExit("\nPlease enter a interrogated region file in Bed format.\n");	
		}
		
		//load and sort Regions File
		Positive[] regionsA = ParseSgrsForParticularRegions.parseRegionFile(regionFile);
		Arrays.sort(regionsA, new PositiveComparator());
		
		//find smallest region
		int shortestLength = Integer.MAX_VALUE;
		for (int i=0; i< regionsA.length; i++){
			int size = regionsA[i].getLength();
			if (size < shortestLength) shortestLength = size;
		}
		
		regions = new HashMap<String,ArrayList<Positive>>();
		chromOrder = new ArrayList<String>();
		
		//Split regions into chromosomes
		for (Positive region: regionsA) {
			String chrom = region.getChromosome();
			if (!(chromOrder.contains(chrom))) {
				chromOrder.add(chrom);
			}
			if (!(regions.containsKey(chrom))) {
				regions.put(chrom, new ArrayList<Positive>());
			}
			regions.get(chrom).add(region);
		}


		//make hash of chromosome text and gc boolean file?
		File[] gcFiles = IO.extractFiles(gcGenomeDir, ".gc");
		chromGCFileMap = Seq.makeChromosomeNameFileHash( gcFiles );
		if (chromGCFileMap == null) Misc.printExit("\nError parsing xxx.gc files.\n");

		
		//Load in interrogated regions
		interrogatedRegions = Region.parseStartStops(interrFile, 0, 0, shortestLength);

	}	
	
	public static void printDocs(){ //fix in post
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Score Enriched Regions: December 2012                     **\n" +
				"**************************************************************************************\n" +
				"ScoreEnrichedRegions determines if the set of regions specified by the user is more or less\n" +
				"enriched than a randomly generated set of regions matched on chromosome, region\n" +
				"length and GC content.  The software determines if each individual region is more/less enriched\n" +
				"as well as the dataset as a whole.  Individual region p-values are caluculated by comparing\n" +
				"the region fold-enrichment to the fold-enrichment of each randomly generated region.  Aggregate\n" +
				"p-values are calculated by comparing the median fold-enrichment of the user-specified dataset\n " +
				"to median fold-enrichment each randomly generated datset. The software uses 1000 \n" +
				"randomly generated regions by default. Pseudocounts are added to moderate fold-change values" +
				

				"\nOptions:\n"+
				"-c Conditions directory containing one directory for each condition with one xxx.bam\n" +
				"       file per biological replica and their xxx.bai indexes. The BAM files should be \n" +
				"       sorted by coordinate using Picard's SortSam.\n" +
				"-r Regions of interest in Bed format (chr, start, stop,...), full path, See,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+
				"-i Interrogated regions in Bed format (chr, start, stop, ...), full path to use in drawing\n" +
				"       random regions.\n"+
				"-g A full path directory containing for chromosome specific gc content boolean arrays. See\n"+
				"       the ConvertFasta2GCBoolean app\n"+
				"-o Full path to the output file\n" +
				
				"\nAdvanced Options:\n"+
				"-x Max per base alignment depth, defaults to 50000. Genes containing such high\n"+
				"       density coverage are ignored. Warnings are thrown.\n"+
				"-f Psuedocounts to each region.  Defaults to 10\n" +
				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/ScoreEnrichedRegions -c\n" +
				"      /Data/TimeCourse/ESCells/ -r regionOfInterest.bed -i sequencedRegions.bed -g gcContent/ \n" +
				"     -o resultsGoHere.txt \n\n" +

		        "**************************************************************************************\n");

	}

}
