package trans.anno;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import trans.main.*;
import util.bio.calc.*;
import util.bio.parsers.MultiFastaParser;
import util.gen.*;
import edu.utah.seq.its.Interval;
import edu.utah.seq.its.IntervalTree;


/**
 * Performs an intersection analysis on binding region files (chr start stop (optional fraction gc)).
 *
 */
public class IntersectRegions {
	//fields
	private File[] firstRegionsFiles = null;
	private File[] secondRegionsFiles = null;
	//private BindingRegion[] one;
	private BindingRegion[][] twos;
	//private BindingRegion[] two;
	private BindingRegion[] oneOverlap;
	private BindingRegion[] twoOverlap;
	private BindingRegion[] oneNonOverlap;
	private BindingRegion[] twoNonOverlap;
	private int gap = 0;
	private int numberRandom = 1000;
	private int[] hitsPerTrial = null;
	private File chromosomeRegionsDirectory = null;
	private boolean intervalFiles = false;
	private boolean writeFiles = false;
	private Histogram histogram;
	private boolean printHistogram = false;
	private double[] minMaxBins = {-100,2400,100};
	private File genomicSequenceDirectory = null;
	private boolean entirelyContained = false; 
	private boolean saveIntersectionPairs = false;
	private ArrayList intersectedRegions = new ArrayList();
	private HashSet<String> performedComparisons = new HashSet<String>();
	private File intersectionFile = null;

	public IntersectRegions(String[] arguments){	
		processArgs(arguments);
		System.out.println("\nLaunching...");
		if (entirelyContained) System.out.println("\tIntersections are scored where 2nd regions are entirely contained in 1st regions.");
		else System.out.println("\tMax gap "+ (gap-1));

		//make StringBuffer to hold total bases and median lengths
		StringBuffer stats = new StringBuffer();
		stats.append("Name\t# Regions\tTotal BPs\tMedian Length\n");

		//make StringBuffer to hold intersection results, add header
		StringBuffer summaryLines = new StringBuffer();
		summaryLines.append("Max_Gap\tFirst\tSecond\tInt_1st_w/_2nd\t_\tInt_2nd_w/_1st\t_");
		if (chromosomeRegionsDirectory != null){
			summaryLines.append("\tLog2(Fold_Enrich)\t_\tCorr_P-Value\tAnti-Corr_P-Value\n");
		}
		else summaryLines.append("\n");
	
		HashMap<String, HashMap<String,IntervalTree<BindingRegion>>> it2List = new HashMap<String, HashMap<String,IntervalTree<BindingRegion>>>();

		//for each first file
		for (int x=0; x<firstRegionsFiles.length; x++){
			//parse BindingRegions for ones
			BindingRegion[] set1;
			if (intervalFiles) {
				set1 = parseIntervalFile(firstRegionsFiles[x]);
			} else {
				set1 = parseRegionsFile(firstRegionsFiles[x]);
			}
			
			String oneName = Misc.removeExtension(firstRegionsFiles[x].getName());
			if (set1.length == 0) {
				System.out.println("\tError: no regions in "+firstRegionsFiles[x]+"! Skipping!");
				continue;
			}
			
			double[] bpsLength = statBindingRegionArray(set1);
			stats.append(firstRegionsFiles[x].getName()+"\t"+set1.length+"\t"+(int)bpsLength[0]+ "\t"+Num.formatNumber(bpsLength[1],4)+"\n");
			
			//for each second file
			for (int i=0; i< secondRegionsFiles.length; i++){

				
				BindingRegion[] set2;
				if (twos[i] == null){
					if (intervalFiles) twos[i] = parseIntervalFile(secondRegionsFiles[i]);
					else twos[i] = parseRegionsFile(secondRegionsFiles[i]);
				}
				set2 = twos[i];				
				
				if (set2.length == 0) {
					System.out.println("\tError: no regions in "+secondRegionsFiles[i]+"! Skipping!");
					it2List.put(secondRegionsFiles[i].toString(), new HashMap<String,IntervalTree<BindingRegion>>());
					continue;
				}
				
				//are they different files?
				if (secondRegionsFiles[i].toString().equals(firstRegionsFiles[x].toString())) continue;
				String f_s = firstRegionsFiles[x].toString() + secondRegionsFiles[i].toString();
				String s_f = secondRegionsFiles[i].toString() + firstRegionsFiles[x].toString();
				if (performedComparisons.contains(f_s) || performedComparisons.contains(s_f)) continue;
				performedComparisons.add(f_s);
				performedComparisons.add(s_f);
				
				//Create or grab stored IntervalTree
				HashMap<String,IntervalTree<BindingRegion>> it2 = it2List.get(secondRegionsFiles[i].toString());
				if (it2 == null) {
					it2 = createIntevervalFromBindingRegion(set2);
					it2List.put(secondRegionsFiles[i].toString(), it2);
				}
				
				//count regions
				String twoName =  Misc.removeExtension(secondRegionsFiles[i].getName());
				summaryLines.append((gap-1)+"\t"+oneName+"\t"+twoName+"\t");

				//make histogram to hold length distribution?
				if (printHistogram) histogram = new Histogram(minMaxBins[0], minMaxBins[1], (int)minMaxBins[2]);

				//overlap
				intersectBindingRegions(set1, set2, it2);

				//save intersections?
				if (saveIntersectionPairs){
					intersectionFile = new File (firstRegionsFiles[x].getParentFile(), oneName+"_Int_"+twoName+"_Pairs.txt");
					String ints = Misc.stringArrayListToString(intersectedRegions, "\n");
					String rankName = "First_Rank";
					if (set1[0].getName()!=null) rankName = "First_Rank\tName";
					
					String rankName2 = "Second_Rank";
					if (set2[0].getName()!=null) rankName2 = "Second_Rank\tName";
					
					
					String headder = rankName+"\tChrom\tStart\tStop\t"+ rankName2+"\tChrom\tStart\tStop\n";
					IO.writeString(headder+ints, intersectionFile);
					intersectedRegions.clear();
				}

				//print to file
				if (writeFiles){
					if (oneOverlap.length !=0) writeBindingRegions(oneOverlap, new File(firstRegionsFiles[x].getParentFile(),oneName+"_Int_"+twoName+".txt"));
					if (twoOverlap.length !=0)writeBindingRegions(twoOverlap, new File(secondRegionsFiles[i].getParentFile(), twoName+"_Int_"+oneName+".txt"));
					if (oneNonOverlap.length !=0)writeBindingRegions(oneNonOverlap, new File(firstRegionsFiles[x].getParentFile(),oneName+"_Diff_"+twoName+".txt"));
					if (twoNonOverlap.length !=0)writeBindingRegions(twoNonOverlap, new File(secondRegionsFiles[i].getParentFile(), twoName+"_Diff_"+oneName+".txt"));
				}
				//calc overlap
				double intOne = (double)oneOverlap.length/(double)set1.length;
				double intTwo = (double)twoOverlap.length/(double)set2.length;
				summaryLines.append(Num.formatNumber(intOne, 4)+"\t("+oneOverlap.length+"/"+set1.length+")\t"+
						Num.formatNumber(intTwo, 4)+"\t("+twoOverlap.length+"/"+set2.length+")");

				//print histogram?
				if (printHistogram) histogram.printScaledHistogram();

				//randomize?
				if (chromosomeRegionsDirectory != null){
					hitsPerTrial = new int[numberRandom];
					intersectRandomRegions(set1,set2);
					int totalHits = Num.sumIntArray(hitsPerTrial);
					double aveHits = (double)totalHits/(double)numberRandom;
					double foldEnrichment = (double)twoOverlap.length/ aveHits;
					summaryLines.append("\t"+Num.formatNumber(Num.log2(foldEnrichment), 4)+ "\t("+twoOverlap.length+"/"+aveHits+")\t");
					double[] pValues = pValue();
					String pValueCorr = "<1/"+numberRandom;
					String pValueAntiCorr = "<1/"+numberRandom;
					if (pValues[0] != 0.0) pValueCorr = Num.formatNumber(pValues[0], 4);
					if (pValues[1] != 0.0) pValueAntiCorr = Num.formatNumber(pValues[1], 4);
					summaryLines.append(pValueCorr+"\t"+pValueAntiCorr);
				}
				summaryLines.append("\n");

			}
		}
		//print stats on lengths and total bps etc
		for (int i=0; i< secondRegionsFiles.length; i++){
			double[] bpsLength = statBindingRegionArray(twos[i]);
			stats.append(secondRegionsFiles[i].getName()+"\t"+twos[i].length+"\t"+(int)bpsLength[0]+ "\t"+Num.formatNumber(bpsLength[1],4)+"\n");
		}
		System.out.println(stats+"\n"+summaryLines);	
	}

	public double[] pValue(){
		double numHitsInReal = twoOverlap.length;
		double counterCorr = 0;
		double counterAntiCorr = 0;
		for (int i=0; i<hitsPerTrial.length; i++){
			if (hitsPerTrial[i] > numHitsInReal) counterCorr++;
			else if (hitsPerTrial[i] < numHitsInReal) counterAntiCorr++;
			else {
				counterCorr++;
				counterAntiCorr++;
			}
		}
		double pvalCorr = counterCorr/((double)numberRandom);
		double pvalAntiCorr = counterAntiCorr/((double)numberRandom);
		
		return new double[]{pvalCorr, pvalAntiCorr};
	}

	/**Increments the hitsPerTrial int[].*/
	public void intersectRandomRegions(BindingRegion[] set1, BindingRegion[] set2){
		//sort two by chromosome		
		Arrays.sort(set1);
		Arrays.sort(set2);
		String chrom = "";
		BindingRegion[] chromSpecificOne = null;
		RandomRegions rr = null;

		//for each binding region in two
		System.out.println("\tProcessing random...");
		
		Date startTime = new Date();
		float totalTime = 0;

		String chromosomeSequence = null;
		if (set2[0].getGcContent() != -1) System.out.println("\nTwo with GC Content ");
		for (int i=0; i< set2.length; i++){
			//different chromosome? make new RandomRegions
			if (set2[i].getChromosome().equals(chrom) == false){
				rr = null;
				chrom = set2[i].getChromosome();
				rr= makeRandomRegions(chrom);
				//fetch from one, those that are chromosome specific
				chromSpecificOne = fetchChromosomeSpecificBRs(set1, chrom);
				
				//set genomic sequence for GC calculation?
				if (genomicSequenceDirectory != null && set2[i].getGcContent() == -1) {
					//look for sequence file
					File chromFastaFile = new File (genomicSequenceDirectory,chrom+".fasta");
					if (chromFastaFile.exists() ==  false )  chromFastaFile = new File (genomicSequenceDirectory,chrom+".fasta.gz");
					if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,chrom+".fa");
					if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,chrom+".fa.gz");
					if (chromFastaFile.exists() ==  false ) Misc.printExit("Error: Cannot find genomic sequence file -> "+chromFastaFile);
					MultiFastaParser fastaParser = new MultiFastaParser();
					fastaParser.parseIt(chromFastaFile);
					chromosomeSequence = fastaParser.getSeqs()[0].toLowerCase();
				}
			}
			
			
			
			if (i % 1000 == 0 && i != 0) {
				Date currTime = new Date();
				float seconds = (currTime.getTime() - startTime.getTime()) / 1000;
				startTime = currTime;
				float minutes = seconds / 60;
				totalTime += minutes;
				float remaining1 = (set2.length - i) * (totalTime / i);
				//float remaining2 = (set2.length - i ) / 1000 * minutes;
				//System.out.println(seconds + " " + minutes + " " + remaining1 + " " + remaining2);
				//System.out.println(String.format("\tFinished testing random regions for %d of %d regions in file2. Time for last 1000 regions: %.2f min.  Estimated Remaining time: %.2f min.",i,set2.length,minutes,remaining1));
			}
			if (chromSpecificOne.length == 0) {
				continue;
			}
			
			//generate numberRandom set of of coordinates
			int[][] randomCoordinates; 

			double gc = 0;
			
			if (genomicSequenceDirectory == null) {
				randomCoordinates = rr.fetchRandomCoordinates(set2[i].getEnd()-set2[i].getStart()+1, numberRandom);
			}
			//match gc content
			else {
				int len = 0;
				int gcCount = 0;
				String sequence = chromosomeSequence.substring(set2[i].getStart(),set2[i].getEnd()).toLowerCase();
				for (char b: sequence.toCharArray()) {
					len++;
					if (b == 'g' || b == 'c') {
						gcCount++;
					}
				}
				gc = (double)gcCount / len ;
				randomCoordinates = rr.fetchRandomCoordinates(set2[i].getEnd()-set2[i].getStart()+1, numberRandom, gc);
			}
			//failed to make?
			if (randomCoordinates == null ) System.err.println(String.format("WARNING: Problem making random coordinates for %s.  GC content: %.2f. Length: %d.",set2[i].simpleSummaryLine(),gc,set2[i].getLength()));
			else {
				//convert to an array of BindingRegion
				BindingRegion[] randomBRs = makeBindingRegions(randomCoordinates, chrom);
//				for (BindingRegion br: randomBRs) {
//					System.out.println(br.getStart() + " " + br.getEnd() + " " + br.getGcContent() + "\n");
//				}
				//increment the hitsPerTrial counter
				intersectRandomBindingRegions(randomBRs, chromSpecificOne);

			}
			
		}
		System.out.println();
	}

	/**Given a sorted array of BindingRegion and a chromosome, returns the appropriate BindingRegions.*/
	public static BindingRegion[] fetchChromosomeSpecificBRs(BindingRegion[] sortedBRs, String chromosome){
		ArrayList al = new ArrayList();
		boolean found = false;
		for (int i=0; i< sortedBRs.length; i++){
			if (sortedBRs[i].getChromosome().equals(chromosome)){
				found = true;
				al.add(sortedBRs[i]);
			}
			else if (found) break;
		}
		BindingRegion[] brs = new BindingRegion[al.size()];
		al.toArray(brs);
		return brs;
	}

	/**Given an array of start stops and a chromosome, makes BindingRegions.*/
	public static BindingRegion[] makeBindingRegions(int[][] regions, String chromosome){
		BindingRegion[] brs = new BindingRegion[regions.length];
		for (int i=0; i< regions.length; i++){
			brs[i] = new BindingRegion(i, chromosome, regions[i][0], regions[i][1]);
		}
		return brs;
	}

	public RandomRegions makeRandomRegions(String chromosome){
		File chromFile = new File (chromosomeRegionsDirectory, chromosome);
		if (chromFile.exists() == false) Misc.printExit("\nCannot find the interrogated regions file for chromosome "+chromosome+" in "+chromosomeRegionsDirectory);
		if (genomicSequenceDirectory != null) return new RandomRegions(chromFile, genomicSequenceDirectory);
		return new RandomRegions(chromFile);
	}

	public static boolean writeBindingRegions(BindingRegion[] brs, File file){
		try{
			//first file
			PrintWriter out = new PrintWriter(new FileWriter(file));
			int num = brs.length;
			if (brs[0].getGcContent() != -1){
				for (int i=0; i<num; i++) out.println(brs[i].simpleSummaryLine()+"\t"+brs[i].getGcContent());
			}
			else for (int i=0; i<num; i++) out.println(brs[i].simpleSummaryLine());
			out.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}

	/**Assumes the same chromosome*/
	public void intersectRandomBindingRegions(BindingRegion[] allTrialsOfSingleBR, BindingRegion[] chrOne){
		//System.out.print("Constructing Tree...");
		IntervalTree<Integer> invTree = this.createIntevervalFromInteger(allTrialsOfSingleBR);
		//System.out.println("done");
		
		HashSet<Integer> found = new HashSet<Integer>();
		
		if (entirelyContained) {
			for (int i=0;i<chrOne.length;i++) {
			    BindingRegion br = chrOne[i];
				ArrayList<Integer> resList = invTree.searchContains(br.getStart(), br.getEnd());
			
				for (Integer idx: resList) {
					if (!found.contains(idx)) {
						hitsPerTrial[idx]++;
						found.add(idx);
					}
				}
			}
		} else {
			for (int i=0;i<chrOne.length;i++) {
			    BindingRegion br = chrOne[i];
				ArrayList<Integer> resList = invTree.search(br.getStart()-gap, br.getEnd()+gap);
				
				for (Integer idx: resList) {
					if (!found.contains(idx)) {
						hitsPerTrial[idx]++;
						found.add(idx);
					}	
				}
			}
		}
	}

	/**Returns true if two entirely contains one.*/
	public static boolean containedWithin(BindingRegion one, BindingRegion two){
		if (two.getStart()<= one.getStart() && two.getEnd()>= one.getEnd()) return true;
		return false;
	}

	public void intersectBindingRegions(BindingRegion[] set1, BindingRegion[] set2, HashMap<String,IntervalTree<BindingRegion>> it2){
		//find intersection
		//use hashmaps to hold non redundant binding regions
		HashSet<BindingRegion> oneOverlapHS = new HashSet<BindingRegion>();
		HashSet<BindingRegion> twoOverlapHS = new HashSet<BindingRegion>();
		
		if (entirelyContained) {
			for (BindingRegion br: set1) {
				String chrom = br.getChromosome();
				if (it2.containsKey(chrom)) {
					ArrayList<BindingRegion> results = it2.get(chrom).searchContains(br.getStart(), br.getEnd());
					for (BindingRegion r: results) {
						oneOverlapHS.add(br);
						twoOverlapHS.add(r);
						if (saveIntersectionPairs) {
							intersectedRegions.add(br.simpleSummaryLineWithRank()+"\t"+r.simpleSummaryLineWithRank()); 
						}
					}
				}
			}
		} else {
			for (BindingRegion br: set1) {
				String chrom = br.getChromosome();
				int minGap = 1000000000;
				if (it2.containsKey(chrom)) {
					ArrayList<BindingRegion> results = it2.get(chrom).search(br.getStart()-gap, br.getEnd()+gap);
					for (BindingRegion r: results) {
						int bpOverlap = br.bpIntersectionSameChromosome(r);				
						if (bpOverlap < minGap) {
							minGap = bpOverlap;
						}
						//direct overlap?  add to total
						if (bpOverlap <= gap ){
							oneOverlapHS.add(br);
							twoOverlapHS.add(r);
							if (saveIntersectionPairs) {
								intersectedRegions.add(br.simpleSummaryLineWithRank()+"\t"+r.simpleSummaryLineWithRank());
							}
						}
					}
					if (printHistogram && minGap < 1000000000) {
						histogram.count(minGap);
					}
				}
			}
		}
		
		//convert to arrays
		oneOverlap = new BindingRegion[oneOverlapHS.size()];
		int counter =0;
		Iterator it = oneOverlapHS.iterator();
		while (it.hasNext()){
			oneOverlap[counter++] = (BindingRegion)it.next();
		}
		twoOverlap = new BindingRegion[twoOverlapHS.size()];
		counter =0;
		it = twoOverlapHS.iterator();
		while (it.hasNext()){
			twoOverlap[counter++] = (BindingRegion)it.next();
		}

		//find difference
		ArrayList nonOverlapOne = new ArrayList();
		ArrayList nonOverlapTwo = new ArrayList();
		for (int i=0; i<set1.length; i++){
			if (oneOverlapHS.contains(set1[i]) == false) nonOverlapOne.add(set1[i]);
		}
		for (int i=0; i<set2.length; i++){
			if (twoOverlapHS.contains(set2[i]) == false) nonOverlapTwo.add(set2[i]);
		}
		//convert to arrays
		oneNonOverlap = new BindingRegion[nonOverlapOne.size()];
		nonOverlapOne.toArray(oneNonOverlap);
		twoNonOverlap = new BindingRegion[nonOverlapTwo.size()];
		nonOverlapTwo.toArray(twoNonOverlap);
	}


	public static void main(String[] args){
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new IntersectRegions(args);
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String firstFiles = null;
		String secondFiles = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': firstFiles = args[i+1]; i++; break;
					case 's': secondFiles = args[i+1]; i++; break;
					case 'n': numberRandom = Integer.parseInt(args[i+1]); i++; break;
					case 'g': gap = Integer.parseInt(args[i+1]); i++; break;
					case 'r': chromosomeRegionsDirectory = new File(args[i+1]); i++; break;
					case 'w': writeFiles = true; break;
					case 'x': saveIntersectionPairs = true; break;
					case 'e': entirelyContained = true; break;
					case 'c': genomicSequenceDirectory = new File(args[i+1]); i++; break;
					case 'p': printHistogram = true; break;
					case 'q': minMaxBins = Num.parseDoubles(args[i+1].split(",")); i++; break;
					case 'i': intervalFiles = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//Add one to gap!
		gap++;
		
		//look for required parameters
		if (firstFiles == null || secondFiles == null ){
			Misc.printExit("\nPlease complete one or more of the following required parameters: -f or -s .\n");
		}

		//parse first files
		String[] files = firstFiles.split(",");
		if (files.length == 1) firstRegionsFiles = IO.extractFiles(new File(firstFiles));
		else firstRegionsFiles = IO.extractFiles(firstFiles);

		//parse second files
		files = secondFiles.split(",");
		if (files.length == 1) secondRegionsFiles = IO.extractFiles(new File(files[0]));
		else secondRegionsFiles = IO.extractFiles(secondFiles);	
		twos = new BindingRegion[secondRegionsFiles.length][];

		//look for others
		if (genomicSequenceDirectory != null && genomicSequenceDirectory.exists() == false){
			Misc.printExit("\nError: cannot find your genomic sequences directory "+genomicSequenceDirectory.getName());
		}
		
		//split interrogated regions file?
		if (chromosomeRegionsDirectory != null){
			if (chromosomeRegionsDirectory.isDirectory() == false){
				//make it from bed file
				System.out.println("\nSplitting background interrogated regions file by chromosome...");
				splitBackgroundRegionsFile(chromosomeRegionsDirectory);
			}
			if (chromosomeRegionsDirectory.exists()== false) Misc.printExit("\nError: cannot find your chromosome split interrogated regions directory "+chromosomeRegionsDirectory.getName());
		}

	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Intersect Regions: May 2017                           **\n" +
				"**************************************************************************************\n" +
				"IR intersects lists of regions (tab delimited: chrom start stop(inclusive)). Random\n" +
				"regions can also be used to calculate a p-value and fold enrichment.\n\n"+

				"-f First regions files, a single file, or a directory of files.\n"+
				"-s Second regions files, a single file, or a directory of files.\n"+
				"-g Max gap, defaults to 0. A max gap of 0 = regions must at least abut or overlap,\n" +
				"      negative values force overlap (ie -1= 1bp overlap, be careful not to exceed the\n" +
				"      length of the smaller region), positive values enable gaps (ie 1=1bp gap).\n"+
				"-e Score intersections where second regions are entirely contained by first regions.\n"+
				"-r Make random regions matched to the second regions file(s) and intersect with the\n"+
				"      first.  Enter either a bed file or full path directory that contains chromosome\n"+
				"      specific interrogated regions files (ie named: chr1, chr2 ...: chrom start stop).\n"+
				"-c Match GC content of second regions file(s) when selecting random regions, rather\n" +
				"      slow. Provide a full path directory text containing chromosome specific genomic\n" +
				"      sequences.\n"+
				"-n Number of random region trials, defaults to 1000.\n" +
				"-w Write intersections and differences.\n" +
				"-x Write paired intersections.\n"+
				"-p Print length distribution histogram for gaps between first and closest second.\n" +
				"-q Parameters for histogram, comma delimited list, no spaces:\n"+
				"       minimum length, maximum length, number of bins.  Defaults to -100, 2400, 100.\n"+


				"\nExample: java -Xmx1500M -jar pathTo/Apps/IntersectRegions -f /data/miRNAs.txt\n" +
				"      -s /data/DroshaLists/ -g 500 -n 10000 -r /data/InterrogatedRegions/\n\n" +

				"\n" +
		"**************************************************************************************\n");		
	}

	public static BindingRegion[] parseRegionsFile(File picksFile){
		BindingRegion[] regions =null;
		try{
			BufferedReader in = IO.fetchBufferedReader(picksFile);
			String line;
			String[] tokens;			
			HashMap<String, BindingRegion> regionsHS = new HashMap<String,BindingRegion>();
			int rank = 1;
			String chromosome;
			int start;
			int stop;
			double gc;
			//chrom, start, stop, fractionGC?
			//  0      1      2      
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.equals("") || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length<3) continue;
				chromosome = tokens[0];
				start = Integer.parseInt(tokens[1]);
				stop = Integer.parseInt(tokens[2]);
				gc =-1;
				String name = null;
				if (tokens.length>3){
					try{
						gc = Double.parseDouble(tokens[3]);
					}catch(Exception e){}
					if (tokens.length == 6) name = tokens[3];
				}
				BindingRegion br = new BindingRegion(rank++, chromosome, start, stop);
				br.setGcContent(gc);
				br.setName(name);
				regionsHS.put(br.simpleSummaryLine(), br);
			}
			regions = new BindingRegion[regionsHS.size()];
			regionsHS.values().toArray(regions);

		}catch (IOException e){
			e.printStackTrace();
		}
		return regions;
	}
	
	public  void splitBackgroundRegionsFile(File bed){
		try{
			PrintWriter out = null;
			BufferedReader in = IO.fetchBufferedReader(bed);
			chromosomeRegionsDirectory = new File (bed.getParentFile(), "ChromSplit_"+Misc.removeExtension(bed.getName()));
			chromosomeRegionsDirectory.mkdir();
			String line;
			String[] tokens;			
			HashMap<String, PrintWriter> chromPW = new HashMap<String,PrintWriter>();
			int badLines = 0;
			//chrom, start, stop, .....
			//  0      1      2     
			while ((line = in.readLine()) !=null) {
				
				//OK line?
				line = line.trim();
				if (line.equals("") || line.startsWith("#")) continue;
				tokens = line.split("\\t");
				if (tokens.length<3) {
					if (badLines++ > 100) throw new IOException ("\nToo many bad/ malformed bed file lines from the interrogated regions file.");
					continue;
				}
				
				//fetch writer
				out = chromPW.get(tokens[0]);
				if (out == null){
					File f = new File (chromosomeRegionsDirectory, tokens[0]);
					out = new PrintWriter( new FileWriter(f));
					chromPW.put(tokens[0], out);
				}
				
				//print it
				out.println(tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]);
			}
			
			//close the print writers
			for (PrintWriter pw : chromPW.values()) pw.close();

		}catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing background bed file.\n");
		}
	}

	public static BindingRegion[] parseIntervalFile(File intervalFile){
		trans.main.Interval[] ints = (trans.main.Interval[])IO.fetchObject(intervalFile);
		BindingRegion[] regions = new BindingRegion[ints.length];
		for (int i=0; i< ints.length; i++) regions[i] = new BindingRegion(ints[i], i);
		return regions;
	}

	/**@return double[2]{total bases,median length}.*/
	public static double[] statBindingRegionArray(BindingRegion[] brs){
		//total bases
		double totalBases = 0;
		double[] lengths = new double[brs.length];
		for (int i=0; i< brs.length; i++){
			lengths[i] = brs[i].getLength();
			totalBases += lengths[i];
		}
		//calculate median length
		Arrays.sort(lengths);
		double medianLength = Num.median(lengths);
		return new double[] {totalBases, medianLength};
	}
	
	private HashMap<String, IntervalTree<BindingRegion>> createIntevervalFromBindingRegion(BindingRegion[] regions) {
		HashMap<String,ArrayList<edu.utah.seq.its.Interval<BindingRegion>>> invLists = new HashMap<String,ArrayList<edu.utah.seq.its.Interval<BindingRegion>>>();
		
		//Split up binding regions by chromosome
		for (BindingRegion br: regions) {
			String chrom = br.getChromosome();
			if (!invLists.containsKey(chrom)) {
				invLists.put(chrom, new ArrayList<edu.utah.seq.its.Interval<BindingRegion>>());
			}
			edu.utah.seq.its.Interval<BindingRegion> interval = new edu.utah.seq.its.Interval<BindingRegion>(br.getStart(), br.getEnd(), br);
			invLists.get(chrom).add(interval);
		}
		
		HashMap<String,IntervalTree<BindingRegion>> invTree = new HashMap<String,IntervalTree<BindingRegion>>();
		for (String chrom: invLists.keySet()) {
			invTree.put(chrom, new IntervalTree<BindingRegion>(invLists.get(chrom),false));
		}
		
		return invTree;
	}
	
	private IntervalTree<Integer> createIntevervalFromInteger(BindingRegion[] regions) {
		ArrayList<Interval<Integer>> invLists = new ArrayList<Interval<Integer>>();
		
		String chrom = null;
		//Split up binding regions by chromosome
		for (int i=0; i<regions.length; i++)  {
			BindingRegion br = regions[i];
			if (chrom == null) {
				chrom = br.getChromosome();
			} else if (!chrom.equals(br.getChromosome())) {
				System.out.println(String.format("More than one chromsome in list: %s %s",chrom,br.getChromosome()));
				System.exit(1);
			}
		
			Interval<Integer> interval = new edu.utah.seq.its.Interval<Integer>(br.getStart(), br.getEnd(), i);
			invLists.add(interval);
		}
		
		IntervalTree<Integer> invTree = new IntervalTree<Integer>(invLists,false);
	
		return invTree;
	}

}
