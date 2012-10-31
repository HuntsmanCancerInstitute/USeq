package trans.anno;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import trans.main.*;
import util.bio.calc.*;
import util.bio.parsers.MultiFastaParser;
import util.gen.*;

/**
 * Performs an intersection analysis on binding region files (chr start stop (optional fraction gc)).
 *
 */
public class IntersectRegions {
	//fields
	private File[] firstRegionsFiles = null;
	private File[] secondRegionsFiles = null;
	private BindingRegion[] one;
	private BindingRegion[][] twos;
	private BindingRegion[] two;
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
	private File intersectionFile = null;

	public IntersectRegions(String[] arguments){	
		processArgs(arguments);
		System.out.println("\nLaunching...");
		if (entirelyContained) System.out.println("\tIntersections are scored where 2nd regions are entirely contained in 1st regions.");
		else System.out.println("\tMax gap "+gap);

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

		//for each first file
		for (int x=0; x<firstRegionsFiles.length; x++){
			//parse BindingRegions for ones
			if (intervalFiles) one = parseIntervalFile(firstRegionsFiles[x]);
			else one = parseRegionsFile(firstRegionsFiles[x]);
			String oneName = Misc.removeExtension(firstRegionsFiles[x].getName());
			if (one.length == 0) {
				System.out.println("\tError: no regions in "+firstRegionsFiles[x]+"! Skipping!");
				continue;
			}
			double[] bpsLength = statBindingRegionArray(one);
			stats.append(firstRegionsFiles[x].getName()+"\t"+one.length+"\t"+(int)bpsLength[0]+ "\t"+Num.formatNumber(bpsLength[1],4)+"\n");

			//for each second file
			for (int i=0; i< secondRegionsFiles.length; i++){
				if (twos[i] == null){
					if (intervalFiles) twos[i] = parseIntervalFile(secondRegionsFiles[i]);
					else twos[i] = parseRegionsFile(secondRegionsFiles[i]);
				}
				two = twos[i];
				if (two.length == 0) {
					System.out.println("\tError: no regions in "+secondRegionsFiles[i]+"! Skipping!");
					continue;
				}
				//count regions
				String twoName =  Misc.removeExtension(secondRegionsFiles[i].getName());
				summaryLines.append(gap+"\t"+oneName+"\t"+twoName+"\t");

				//make histogram to hold length distribution?
				if (printHistogram) histogram = new Histogram(minMaxBins[0], minMaxBins[1], (int)minMaxBins[2]);

				//overlap
				intersectBindingRegions();

				//save intersections?
				if (saveIntersectionPairs){
					intersectionFile = new File (firstRegionsFiles[x].getParentFile(), oneName+"_Int_"+twoName+"_Pairs.txt");
					String ints = Misc.stringArrayListToString(intersectedRegions, "\n");
					String rankName = "First_Rank";
					if (one[0].getName()!=null) rankName = "First_Rank\tName";
					
					String rankName2 = "Second_Rank";
					if (two[0].getName()!=null) rankName2 = "Second_Rank\tName";
					
					
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
				double intOne = (double)oneOverlap.length/(double)one.length;
				double intTwo = (double)twoOverlap.length/(double)two.length;
				summaryLines.append(Num.formatNumber(intOne, 4)+"\t("+oneOverlap.length+"/"+one.length+")\t"+
						Num.formatNumber(intTwo, 4)+"\t("+twoOverlap.length+"/"+two.length+")");

				//print histogram?
				if (printHistogram) histogram.printScaledHistogram();

				//randomize?
				if (chromosomeRegionsDirectory != null){
					hitsPerTrial = new int[numberRandom];
					intersectRandomRegions();
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
	public void intersectRandomRegions(){
		//sort two by chromosome		
		Arrays.sort(two);
		Arrays.sort(one);
		String chrom = "";
		BindingRegion[] chromSpecificOne = null;
		RandomRegions rr = null;

		//for each binding region in two
		System.out.println("\tProcessing random...");

		String chromosomeSequence = null;
		if (two[0].getGcContent() != -1) System.out.println("\nTwo with GC Content ");
		for (int i=0; i< two.length; i++){
			//different chromosome? make new RandomRegions
			if (two[i].getChromosome().equals(chrom) == false){
				rr = null;
				chrom = two[i].getChromosome();
				rr= makeRandomRegions(chrom);
				//fetch from one, those that are chromosome specific
				chromSpecificOne = fetchChromosomeSpecificBRs(one, chrom);
				//set genomic sequence for GC calculation?
				if (genomicSequenceDirectory != null && two[i].getGcContent() == -1) {
					//look for sequence file
					File chromFastaFile = new File (genomicSequenceDirectory,chrom+".fasta");
					if (chromFastaFile.exists() ==  false )  chromFastaFile = new File (genomicSequenceDirectory,chrom+".fasta.gz");
					if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,chrom+".fa");
					if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,chrom+".fa.gz");
					if (chromFastaFile.exists() ==  false ) Misc.printExit("Error: Cannot find genomic sequence file -> "+chromFastaFile);
					MultiFastaParser fastaParser = new MultiFastaParser();
					fastaParser.parseIt(chromFastaFile);
					chromosomeSequence = fastaParser.getSeqs()[0];
				}
			}
			//generate numberRandom set of of coordinates
			int[][] randomCoordinates; 

			if (genomicSequenceDirectory == null) randomCoordinates = rr.fetchRandomCoordinates(two[i].getEnd()-two[i].getStart()+1, numberRandom);
			//match gc content
			else {
				double gc = two[i].getGcContent();
				randomCoordinates = rr.fetchRandomCoordinates(two[i].getEnd()-two[i].getStart()+1, numberRandom, gc);
			}
			//failed to make?
			if (randomCoordinates == null ) System.err.println("\nERROR: Problem making random coordinates for "+two[i].simpleSummaryLine()+" Skipping!\n");
			else {
				//convert to an array of BindingRegion
				BindingRegion[] randomBRs = makeBindingRegions(randomCoordinates, chrom);
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
		int numTwo = chrOne.length;
		if (entirelyContained){
			for (int i=0; i<numberRandom; i++){
				for (int j=0; j<numTwo; j++){
					if (containedWithin(allTrialsOfSingleBR[i], chrOne[j])){
						hitsPerTrial[i]++;
						break;
					}
				}
			}
		}
		else{
			for (int i=0; i<numberRandom; i++){
				for (int j=0; j<numTwo; j++){
					//get bp overlap; 0=abut, -1=1bp overlap, 1=1bp gap
					int bpOverlap = allTrialsOfSingleBR[i].bpIntersectionSameChromosome(chrOne[j]);
					//direct overlap?  add to total
					if (bpOverlap <= gap ){
						hitsPerTrial[i]++;
						break;
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

	public void intersectBindingRegions(){
		//find intersection
		//use hashmaps to hold non redundant binding regions
		HashSet oneOverlapHS = new HashSet();
		HashSet twoOverlapHS = new HashSet();
		//run through both sets
		int numOne = one.length;
		int numTwo = two.length;
		if (entirelyContained){
			//for each one
			for (int i=0; i<numOne; i++){
				//for each two
				for (int j=0; j<numTwo; j++){
					//check if on same chromosome
					if (one[i].getChromosome().equals(two[j].getChromosome())){
						if (containedWithin(two[j], one[i])){
							oneOverlapHS.add(one[i]);
							twoOverlapHS.add(two[j]);
							if (saveIntersectionPairs) intersectedRegions.add(one[i].simpleSummaryLineWithRank()+"\t"+two[j].simpleSummaryLineWithRank()); 
						}
					}
				}
			}
		}
		else{
			//for each one
			for (int i=0; i<numOne; i++){
				//for each two
				int minGap = 1000000000;
				for (int j=0; j<numTwo; j++){
					//check if on same chromosome
					if (one[i].getChromosome().equals(two[j].getChromosome())){
						//get bp overlap; 0=abut, -1=1bp overlap, 1=1bp gap
						int bpOverlap = one[i].bpIntersectionSameChromosome(two[j]);				
						if (bpOverlap < minGap) minGap = bpOverlap;
						//direct overlap?  add to total
						if (bpOverlap <= gap ){
							oneOverlapHS.add(one[i]);
							twoOverlapHS.add(two[j]);
							if (saveIntersectionPairs) intersectedRegions.add(one[i].simpleSummaryLineWithRank()+"\t"+two[j].simpleSummaryLineWithRank());
						}
					}
				}
				if (printHistogram && minGap < 1000000000) histogram.count(minGap);
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
		for (int i=0; i<numOne; i++){
			if (oneOverlapHS.contains(one[i]) == false) nonOverlapOne.add(one[i]);
		}
		for (int i=0; i<numTwo; i++){
			if (twoOverlapHS.contains(two[i]) == false) nonOverlapTwo.add(two[i]);
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
				"**                            Intersect Regions: May 2012                           **\n" +
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
		Interval[] ints = (Interval[])IO.fetchObject(intervalFile);
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

}
