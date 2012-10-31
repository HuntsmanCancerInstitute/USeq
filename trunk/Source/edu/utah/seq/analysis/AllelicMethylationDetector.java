package edu.utah.seq.analysis;

import java.io.*;

import java.util.regex.*;
import java.util.*;
import edu.utah.seq.useq.data.RegionScoreText;
import net.sf.samtools.*;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.bio.seq.Seq;
import util.gen.*;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.parsers.*;




/** Application for examining select regions for a bimodal distribution of methylation in reads.
 * @author Nix
 * */
public class AllelicMethylationDetector {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private File[] bamFiles;
	private HashMap<String, File> chromosomeFastaFiles = new HashMap<String, File>();
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private File bedFile;
	private int minimumCsInAlignment = 6;
	private double minimumReadsInRegion = 15;
	private float minimumFractionMethylation = 0.4f;
	private float maximumFractionMethylation = 0.6f;
	private int maximumNumberReads = 10000;

	//internal fields
	private AllelicRegionMaker irm;
	private RegionScoreText[] rst;
	private HashMap<String, RegionScoreText[]> chromRegions; 
	private String chromosome;
	private String[] chromosomes;
	private SAMFileReader[] samReaders;
	private String genomicSequence = null;
	private NovoalignBisulfiteParser novoalignBisulfiteParser;
	private ArrayList<PutativeImprintRegion> pirAL = new ArrayList<PutativeImprintRegion>();
	private PutativeImprintRegion[] pirs;
	
	//for cumulative stats
	private Histogram allHistogram = new Histogram(0,1, 20);
	private float allFractionSum = 0;
	private int[] allNonCon = {0,0};
	private ArrayList<int[]> allPairCounts = new ArrayList<int[]>();

	//constructor
	/**Stand alone.*/
	public AllelicMethylationDetector(String[] args){
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);
		
		allHistogram.setSkipZeroBins(false);

		//launch
		findRegions();

		//any regions?
		if (pirAL.size() == 0) System.out.println("\nNo regions found passing filters?!\n");
		else {
			//make regions
			pirs = new PutativeImprintRegion[pirAL.size()];
			pirAL.toArray(pirs);

			//score with Ken's method
			scoreRegions();

			//sort and print regions
			Arrays.sort(pirs);
			printRegions();

		}
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}
	
	public void scoreRegions(){
		System.out.println("\nScoring " + pirs.length + " regions for log likelihood of U distribution (it's R so patience required)...");
		ArrayList<int[][]> counts = new ArrayList<int[][]>();
		for (PutativeImprintRegion p: pirs){
			counts.add(p.nonConReadCounts);
		}
		double[][] stats = logLikeRatio(counts, saveDirectory, fullPathToR);
		for (int i=0; i< pirs.length; i++) {
			stats[i][3] = Num.minus10log10(stats[i][3]);
			pirs[i].kenStats = stats[i];
			pirs[i].sortBy = pirs[i].histogram.getTotalBinCounts();
		}
		
	}

	public void printRegions(){

		try {
			File bedFile = new File(saveDirectory, "putativeImprintedRegions.bed");
			PrintWriter bedOut = new PrintWriter (new FileWriter ( bedFile));
			File regionFile = new File(saveDirectory, "putativeImprintedRegions.txt");
			PrintWriter regionsOut = new PrintWriter (new FileWriter ( regionFile));
			
			PutativeImprintRegion composite = null;
			
			for (PutativeImprintRegion i : pirs){
				if (i.bedLine.startsWith("Composite")) composite = i;
				else {
				bedOut.println(i.bedLine);
				regionsOut.println(i);
				if (i.histogram.getTotalBinCounts() > 300)  i.histogram.printScaledHistogram(regionsOut);
				else i.histogram.printHistogram(regionsOut);
				regionsOut.println();
				}
			}
			
			//print composite
			bedOut.println(composite.bedLine);
			regionsOut.println(composite);
			if (composite.histogram.getTotalBinCounts() > 300)  composite.histogram.printScaledHistogram(regionsOut);
			else composite.histogram.printHistogram(regionsOut);

			bedOut.close();
			regionsOut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void makeSamReaders(){
		samReaders = new SAMFileReader[bamFiles.length];
		for (int i=0; i< samReaders.length; i++) samReaders[i] = new SAMFileReader(bamFiles[i]);
	}

	public void closeSamReaders(){
		for (int i=0; i< samReaders.length; i++) samReaders[i].close();
	}


	public void findRegions(){
		//load regions
		if (convertedPointDirs!= null) {
			irm = new AllelicRegionMaker(convertedPointDirs, nonConvertedPointDirs);
			chromosomes = irm.getChromosomes();
		}
		else chromosomes = Misc.setToStringArray(chromRegions.keySet());

		//make readers on each bam file
		makeSamReaders();

		//for each chromosome 
		System.out.println("\nScanning regions by chromosome");
		for (String c: chromosomes){
			chromosome = c;

			//fetch fasta file
			File seqFile = chromosomeFastaFiles.get(chromosome);
			if (seqFile == null) {
				System.err.println("\n\tWarning, could not find a fasta file for "+chromosome+", skipping!");
				continue;
			}				
			MultiFastaParser mfp = new MultiFastaParser(seqFile);
			genomicSequence = mfp.getSeqs()[0];

			//fetch regions
			if (irm != null) {
				rst = irm.scan(chromosome);
				if (rst == null){
					System.err.println("\tNo regions for "+chromosome);
					continue;
				}
			}
			else rst = chromRegions.get(chromosome);

			//make a bisulfite parser
			novoalignBisulfiteParser = new NovoalignBisulfiteParser(chromosome, genomicSequence);

			scanChromosome();

		}
		
		//composite
		float meanMethyl = allFractionSum/(float)allHistogram.getTotalBinCounts();
		if (allHistogram.getTotalBinCounts() > minimumReadsInRegion && meanMethyl>= minimumFractionMethylation && meanMethyl <= maximumFractionMethylation) {
			float[] binFractions = fetchBinHistogramFractions(allHistogram);
			if (binFractions != null) {
				PutativeImprintRegion pir = new PutativeImprintRegion(binFractions, meanMethyl, "Composite of all passing regions",  allHistogram, allPairCounts);
				pirAL.add(pir);
			}
		}

		//close streams
		closeSamReaders();

	}


	public float parseFloat(String f){
		if (f.equals("Inf") ) return Float.MIN_VALUE;
		else return Float.parseFloat(f);
	}

	/**Interbase coordinates!*/
	public ArrayList<SAMRecord> fetchAlignments (ExonIntron ei, SAMFileReader reader){
		ArrayList<SAMRecord> al = new ArrayList<SAMRecord>();
		SAMRecordIterator i = reader.queryOverlapping(chromosome, ei.getStart()+1, ei.getEnd());
		while (i.hasNext()) al.add(i.next());
		i.close();
		return al;
	}

	/**Interbase coordinates!*/
	public ArrayList<SAMRecord>[] fetchOverlappingAlignments (ExonIntron[] ei, SAMFileReader reader){
		ArrayList<SAMRecord>[] als = new ArrayList[ei.length];
		for (int i=0; i< ei.length; i++) als[i] = fetchAlignments( ei[i], reader);
		return als;
	}
	/**Interbase coordinates!*/
	public ArrayList<SAMRecord>[][] fetchOverlappingAlignments (ExonIntron[] ei, SAMFileReader[] readers){
		ArrayList<SAMRecord>[][] ohMy = new ArrayList[readers.length][];
		for (int i=0; i< readers.length; i++) ohMy[i] = fetchOverlappingAlignments(ei, readers[i]);
		return ohMy;
	}
	/**Interbase coordinates!*/
	public ArrayList<SAMRecord> fetchOverlappingAlignmentsCombine (ExonIntron[] ei, SAMFileReader[] readers){
		ArrayList<SAMRecord> ohMy = new ArrayList<SAMRecord>();
		HashSet<String> names = new HashSet<String>();
		for (int i=0; i< readers.length; i++) {
			ArrayList<SAMRecord>[] al = fetchOverlappingAlignments(ei, readers[i]);
			for (ArrayList<SAMRecord> a: al) {
				for (SAMRecord sam: a){
					if (names.contains(sam.getReadName()) == false){
						ohMy.add(sam);
						names.add(sam.getReadName());
					}
				}
			}
		}
		return ohMy;
	}



	/**Interbase coordinates!*/
	public void fetchAlignments (RegionScoreText region, SAMFileReader reader, ArrayList<SAMRecord> samRecordsAL){
		SAMRecordIterator i = reader.queryOverlapping(chromosome, region.getStart()+1, region.getStop());
		int readCount = 0;
		while (i.hasNext()) {
			if (readCount++ > maximumNumberReads) {
				System.err.println("WARNING: max read count reached for -> '"+region.toString()+"', skipping remaining reads.");
				i.close();
				return;
			}
			SAMRecord sam = i.next();
			if (sam.getReadUnmappedFlag()) continue;
			
			//fetch blocks of actual alignment
			ArrayList<int[]> blocks = SamAlignment.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
			//check to see if any intersects
			for (int[] b : blocks){
				if (region.intersects(b[0], b[1])){
					samRecordsAL.add(sam);
					break;
				}
			}
		}
		i.close();
	}

	/**Interbase coordinates!*/
	public ArrayList<SAMRecord> fetchOverlappingAlignments (RegionScoreText region, SAMFileReader[] readers){
		ArrayList<SAMRecord> ohMy = new ArrayList<SAMRecord>();
		for (int i=0; i< readers.length; i++)  fetchAlignments(region, readers[i], ohMy);
		return ohMy;
	}

	/**Collects number of observations under each gene's exons.*/
	private void scanChromosome(){
		System.out.println("\t"+chromosome+"\t"+rst.length+ " regions");
		//for each gene 
		for (RegionScoreText region : rst){
			Histogram histogram = new Histogram(0,1, 20);
			histogram.setSkipZeroBins(false);
			int start = region.getStart();
			int stop = region.getStop();

			//fetch the overlapping alignments [replicas][exons]
			ArrayList<SAMRecord> samAL = fetchOverlappingAlignments (region, samReaders);

			//fetch read pairs loaded with bisulfite calls
			ArrayList<AlignmentPair> pairs = novoalignBisulfiteParser.parseSamAlignments(samAL);

			//iterate through reads 
			float fractionSum = 0;
			ArrayList<int[]> pairCounts = new ArrayList<int[]>();
			for (AlignmentPair pa : pairs){
				//fetch counts for mCG contexts only
				int[] nonCon = pa.getmCGCounts(start, stop);
				//enough bases to score?
				if ((nonCon[0]+nonCon[1]) < minimumCsInAlignment) continue;
				float fraction = calculateFractionNon(nonCon);
				pairCounts.add(new int[]{nonCon[0], nonCon[1]});
				fractionSum += fraction;
				histogram.count(fraction);
			}
			//save it?
			float meanMethyl = fractionSum/(float)histogram.getTotalBinCounts();
			if (histogram.getTotalBinCounts() > minimumReadsInRegion && meanMethyl>= minimumFractionMethylation && meanMethyl <= maximumFractionMethylation) {
				float[] binFractions = fetchBinHistogramFractions(histogram);
				if (binFractions != null) {
					PutativeImprintRegion pir = new PutativeImprintRegion(binFractions, meanMethyl, region.getBedLine(chromosome),  histogram, pairCounts);
					pirAL.add(pir);
					//for all
					allPairCounts.addAll(pairCounts);
					allFractionSum += fractionSum;
					try {
						allHistogram.addCounts(histogram);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}

		}

	}

	private class PutativeImprintRegion implements Comparable<PutativeImprintRegion>{
		float[] binFractions;
		float meanMethylation;
		int numberReads;
		String bedLine;
		Histogram histogram;
		int[][] nonConReadCounts;
		//'logLR','p1','p2','p-value'
		double[] kenStats;
		double sortBy =0;


		private PutativeImprintRegion(float[] binFractions, float meanMethylation, String bedLine, Histogram histogram,ArrayList<int[]> counts){
			this.binFractions = binFractions;
			this.meanMethylation = meanMethylation;
			this.bedLine = bedLine;
			this.histogram = histogram;
			numberReads = (int)histogram.getTotalBinCounts();
			nonConReadCounts = new int[counts.size()][2];
			counts.toArray(nonConReadCounts);
		}
		
		public int compareTo(PutativeImprintRegion pir){
			if (pir.sortBy > this.sortBy) return 1;
			if (pir.sortBy < this.sortBy) return -1;
			return 0;
			
		}

		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(bedLine);
			sb.append("\n");

			sb.append(meanMethylation);
			sb.append("\tMean Alignment Methylation\n");
			
			sb.append(kenStats[0]+"\t"+kenStats[1]+"\t"+kenStats[2]+"\t"+kenStats[3]+"\t");
			sb.append("LogLR p1 p2 -10Log10(p-val) (Ken's test stats)\n");
			
			sb.append((int)histogram.getTotalBinCounts());
			sb.append("\tNumber Reads\n");

			sb.append(binFractions[0]);
			sb.append(" ");
			sb.append(binFractions[1]);
			sb.append(" ");
			sb.append(binFractions[2]);
			sb.append("\tBin Fractions");

			return sb.toString();
		}

	}

	/**Returns null if fraction bin failed or the first middle and last.*/
	public float[] fetchBinHistogramFractions(Histogram h){
		float[] fractions = calculateSplits(h);
		if (fractions[0] < .25 || fractions[2] < .25 || fractions[1] > .15) return null;
		return fractions;
	}

	/**Returns the fraction found*/
	public float[] calculateSplits(Histogram h){
		int[] counts = h.getBinCounts();
		//sum first 7
		float first = 0;
		for (int i=0; i<7; i++) first+= counts[i];
		//sum middle 5
		float middle = 0;
		for (int i=7; i<12; i++) middle+= counts[i];
		//sum last 8
		float last = 0;
		for (int i=12; i<counts.length; i++) last+= counts[i];
		float total = first + middle + last;
		return new float[]{first/total, middle/total, last/ total};
	}


	public static float calculateFractionNon (int[] nonCon){
		float non = nonCon[0] +1;
		float con = nonCon[1] +1;
		return non/(non+con);
	}

	private static String kensRFunction1 = 
		"CalcBin <- function(y,p){\n"+
		"# y is a pair of integers\n"+
		"# p is a proportion strictly between 0 and 1\n"+
		"# Note: A factor that cancels out in the likelihood ratio is ignored.\n"+
		"  return(p^y[1]*(1-p)^y[2])\n"+
		"}\n"+
		"\n"+
		"CalcBinLogLike <- function(x){\n"+
		"# x is a matrix of integer pairs. The maximum of the binomial likelihood \n"+
		"# occurs at the sample proportion of the first column divided by the total.\n"+
		"# Note: A factor that cancels out of the likelihood ratio is ignored.\n"+
		"  a <- apply(x,2,sum)\n"+
		"  p0 <- a[1]/sum(a)\n"+
		"  bin <- log(apply(x,1,CalcBin,p0))\n"+
		"  return(sum(bin))\n"+
		"}\n"+
		"\n"+
		"CalcMixLogLike <- function(x,del){\n"+
		"# x is a matrix of integer pairs. \n"+
		"# del is a difference from the overall proportion\n"+
		"# we are assuming a 50% mixture\n"+
		"  a <- apply(x,2,sum)\n"+
		"  p0 <- a[1]/sum(a)\n"+
		"  p1 <- p0 - del\n"+
		"  p2 <- p0 + del\n"+
		"  mix <- log(0.5) + log(apply(x,1,CalcBin,p1) + apply(x,1,CalcBin,p2))\n"+
		"  return(sum(mix))\n"+
		"}\n"+
		"\n"+
		"\n"+
		"CalcMaxLogLike <- function(x,inc =0.001){\n"+
		"# Calculates the final likelihood ratio by searching for a 2-component \n"+
		"# mixture centered at the mean proportion.\n"+
		"# inc is the increments for the search (0 < inc < min(p0,1-p0))\n"+
		"  CalcMixLogLike1 <- function(s){CalcMixLogLike(x,s)}\n"+
		"  a <- apply(x,2,sum)\n"+
		"  p0 <- a[1]/sum(a)\n"+
		"  p <- min(p0,1-p0)\n"+
		"  inc1 <- min(inc,p)\n"+
		"  n <- floor(p/inc)-1\n"+
		"  val <- as.array(inc*(0:n))\n"+
		"  y <- apply(val,1,CalcMixLogLike1)\n"+
		"  m <- max(y)\n"+
		"  v <- val[y == m]\n"+
		"  r <- c(m,p0-v,p0+v)\n"+
		"  return(r)\n"+
		"}\n"+
		"\n"+
		"CalculatePvalue <- function(loglikratio){\n"+
		"# Calculates a p-value using the assumption that twice the log likelihood rato \n"+
		"# has a central chi-square distribution with one df\n"+
		"return(pchisq(2*loglikratio,df = 1,lower.tail = FALSE))\n"+
		"}\n"+
		"\n"+
		"CalcLogLikeRatioStat <- function(x,c){\n"+
		"# Evaluates the difference between the a two-component mixture model and \n"+
		"# a single component mode.\n"+
		"# Returns loglikelihood ratio and the two proportions\n"+
		"a <- CalcBinLogLike(x)\n"+
		"b <- CalcMaxLogLike(x,c)\n"+
		"d <- CalculatePvalue(b[1]-a)\n"+
		"likrat <- b[1]-a\n"+
		"r <- cbind(likrat,b[2],b[3],d)\n"+
		"colnames(r) <- c('logLR','p1','p2','p-value')\n"+
		"return(r)\n"+
		"}\n\n";

	/**Runs Ken Boucher's logLikelihood ratio test for the U distribution in R. counts[regions][#NonCon, #Con].
	 * @return double[regions]['logLR','p1','p2','p-value'] best region is one with small p1 and large p2*/
	public static double[][] logLikeRatio (ArrayList<int[][]> counts, File tempDirectory, File fullPathToR) {
		double[][] values = null;
		try {
			//make File objects
			String randomWord = Passwords.createRandowWord(6)+".txt";
			File rScriptFile = new File(tempDirectory, "rScript_"+randomWord);
			File rOutFile = new File (tempDirectory, "rOut_"+randomWord);
			File rResultsFile = new File (tempDirectory, "rResults_"+randomWord);
			File matrixFile = new File (tempDirectory, "countMatrix_"+randomWord);

			//write matrix to file and keep track of stops
			PrintWriter dataOut = new PrintWriter( new FileWriter(matrixFile));
			int numberRegions = counts.size();
			int runningCount = 0;
			StringBuilder stopsSB = new StringBuilder("stops = c(");

			//for each region
			for (int[][] regionCounts: counts){
				//for each observation
				for (int[] obs: regionCounts){
					dataOut.print(obs[0]);
					dataOut.print("\t");
					dataOut.println(obs[1]);
					runningCount++;
				}
				stopsSB.append(runningCount);
				stopsSB.append(",");
			}
			dataOut.close();

			//close stops
			String stops = stopsSB.toString();
			stops = stops.substring(0, stops.length()-1) + ")\n";

			//make script
			StringBuilder script = new StringBuilder(kensRFunction1);
			script.append("numberStops = "+numberRegions +"\n");
			script.append(stops);
			script.append("bigMatrix = read.table('"+matrixFile +"') \n");
			script.append("results = matrix(nrow="+numberRegions+", ncol=4) \n");
			script.append("start = 1 \n");
			script.append("for (i in 1:numberStops){ \n");
			script.append("stop = stops[i] \n");
			script.append("subMatrix = bigMatrix[start : stops[i],] \n");
			script.append("results[i,] = CalcLogLikeRatioStat(subMatrix, 0.001) \n");
			script.append("start = stop +1 \n");
			script.append("} \n");
			script.append("write.table(results, file='"+rResultsFile.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE, sep = \"\t\") ");
			IO.writeString(script.toString(), rScriptFile);

			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					rScriptFile.getCanonicalPath(),
					rOutFile.getCanonicalPath()};			
			//execute
			IO.executeCommandLine(command);
			//load results
			values = Num.loadDoubleMatrix(rResultsFile);
			//check
			if (values == null || values.length != numberRegions) throw new Exception ("Number of results from R does not match the number of regions?! See tempFiles xxx+"+randomWord+" and try executing in R command shell.");

			//clean up
			matrixFile.deleteOnExit();
			rResultsFile.deleteOnExit();
			rOutFile.deleteOnExit();
			rScriptFile.deleteOnExit();
			rOutFile.deleteOnExit();

		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}

	public static void main(String[] args) {
		
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AllelicMethylationDetector(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File fastaDir = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': bamFiles = IO.extractFiles(args[++i], ".bam"); break;
					case 'f': fastaDir = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'c': convertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': nonConvertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'e': minimumCsInAlignment = Integer.parseInt(args[++i]); break;
					case 'a': minimumReadsInRegion = Double.parseDouble(args[++i]); break;
					case 'm': minimumFractionMethylation = Float.parseFloat(args[++i]); break;
					case 'x': maximumFractionMethylation = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//user defined regions?
		if (bedFile != null){
			convertedPointDirs = null;
			nonConvertedPointDirs = null;
			chromRegions = Bed.parseBedFile(bedFile, true);
		}
		//look for point directories?
		else if (convertedPointDirs != null || nonConvertedPointDirs != null){
			if (convertedPointDirs == null || convertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your converted PointData directories(s)!\n");
			//only one directory look deeper
			if (convertedPointDirs.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(convertedPointDirs[0]);
				if (otherDirs != null && otherDirs.length > 0) convertedPointDirs = otherDirs;
			}
			//nonConverted data
			if (nonConvertedPointDirs == null || nonConvertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your non converted PointData directories(s)!\n");
			//only one directory look deeper
			if (nonConvertedPointDirs.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(nonConvertedPointDirs[0]);
				if (otherDirs != null && otherDirs.length > 0) nonConvertedPointDirs = otherDirs;
			}
		}
		else Misc.printErrAndExit("\nPlease enter a regions file to use in scoring regions OR provide converted and non converted PointData to scan the genome.\n");
			
		

		//look for bam files
		if (bamFiles == null || bamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any treatment xxx.bam files?\n");

		//look for bai index files
		OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, false);

		//load fastaFiles into hash
		if (fastaDir == null) Misc.printErrAndExit("\nError: cannot find any fasta sequence files?\n");
		chromosomeFastaFiles = Seq.fetchChromosomeFastaFileHashMap(fastaDir);

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdir();

		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Allelic Methylation Detector:  September 2012                **\n" +
				"**************************************************************************************\n" +
				"AMD identifies regions displaying allelic methylation, e.g. ~50% average mCG\n" +
				"methylation yet individual read pairs show a bimodal fraction distribution of either\n" +
				"fully methylated or unmethylated.  \n\n"+

				"Options:\n"+
				"-s Save directory.\n"+
				"-f Fasta file directory.\n"+
				"-t BAM file directory containing one or more xxx.bam file with their associated xxx.bai\n" +
				"       index. The BAM files should be sorted by coordinate and have passed Picard\n" +
				"       validation.\n" +
				"-a Minimum number alignments per region, defaults to 15.\n"+
				"-e Minimum number Cs in each alignment, defaults to 6\n"+
				"-m Minimum region fraction methylation, defaults to 0.4\n"+
				"-x Maximum region fraction methylation, defaults to 0.6\n"+
				"-r Full path to R, defaults to /usr/bin/R\n"+
				"-c Converted CG context PointData directories, full path, comma delimited. These \n" +
				"       should contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. Use the ParsePointDataContexts on the output of the\n" +
				"       NovoalignBisulfiteParser to select CG contexts. \n" +
				"-n Non-converted PointData directories, ditto. \n" +
				"-b Provide a bed file (chr, start, stop,...), full path, to scan a list of regions\n" +
				"       instead of the genome.  See, http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+ 

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/ beta! \n\n" +

		"**************************************************************************************\n");

	}
}
