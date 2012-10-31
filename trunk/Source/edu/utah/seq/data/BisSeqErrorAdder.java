package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.cluster.HierarchicalClustering;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;
import edu.utah.seq.analysis.BaseContext;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.tpmap.*;

/**For processing bisulfite sequencing data
 * @author Nix
 * */
public class BisSeqErrorAdder {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private File saveDirectory;
	private File convertedSaveDirectory;
	private File nonConvertedSaveDirectory;
	private double fractionNonConverted;
	private String focusChromosome = "chrLambda";
	private double plusCorr = 0;
	private double minusCorr = 0;

	//internal fields
	private String[] chromosomes;
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;

	//by chromosome
	private String chromosome;


	//constructors
	/**Stand alone.*/
	public BisSeqErrorAdder(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		System.out.println("\nCalculating read count stats...");
		fetchDataHashes();
		fetchAllChromosomes();
		
		//score lambda
		System.out.println("\nScoring "+focusChromosome+ "...");
		scoreLambda();
		
		//convert it
		System.out.println("\nConverting...");
		convert();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods
	
	public void scoreLambda(){
		//find lambda and use to calculate multiplier
		for (String c: chromosomes){
			chromosome = c;
			if (chromosome.startsWith(focusChromosome) == false) continue;

			//plus
			PointData[] pd = convertedPlusPointData.get(chromosome);
			PointData converted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			pd = nonConvertedPlusPointData.get(chromosome); 
			PointData nonConverted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			double numCon = converted.getInfo().getScoreTotal();
			double numNonCon = nonConverted.getInfo().getScoreTotal();
			double total = numCon+numNonCon;
			//check fraction
			double frac = numNonCon/total;
			System.out.println ("\tCurrent fraction non converted + "+chromosome +" "+Num.formatNumber(frac, 4));
			if (fractionNonConverted < frac )Misc.printErrAndExit("\nTarget fraction non converted is less than current fraction! Aborting\n");
			double toFlip = fractionNonConverted * total - numNonCon;
			plusCorr = toFlip/numCon;
			System.out.println("\t\tFraction con obs to flip to nonCon "+Num.formatNumber(plusCorr, 4));

			//minus
			pd = convertedMinusPointData.get(chromosome);
			converted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			pd = nonConvertedMinusPointData.get(chromosome); 
			nonConverted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			numCon = converted.getInfo().getScoreTotal();
			numNonCon = nonConverted.getInfo().getScoreTotal();
			total = numCon+numNonCon;
			//check fraction
			frac = numNonCon/total;
			System.out.println ("\tCurrent fraction non converted - "+chromosome +" "+Num.formatNumber(frac, 4));
			if (fractionNonConverted < frac )Misc.printErrAndExit("\nTarget fraction non converted is less than current fraction! Aborting\n");
			toFlip = fractionNonConverted * total - numNonCon;
			minusCorr = toFlip/numCon;
			System.out.println("\t\tFraction con obs to flip to nonCon "+Num.formatNumber(minusCorr, 4));
			
			break;
			
		}
		if (plusCorr == 0 || minusCorr == 0) Misc.printErrAndExit("\nNo "+focusChromosome+" found to establish observations needing non conversion to meet your target fraction non conversion rate?");
		
	}

	public void convert() {

		//calc totals
		double totalCon = PointData.totalScoreMultiPointData(convertedPlusPointData);
		totalCon += PointData.totalScoreMultiPointData(convertedMinusPointData);
		double totalNonCon = PointData.totalScoreMultiPointData(nonConvertedPlusPointData);
		totalNonCon += PointData.totalScoreMultiPointData(nonConvertedMinusPointData);
		
		double postTotalNonCon =0;
		double postTotalCon =0;
		
		System.out.println("\nChromosome\tNonCon\tCon\tFracNonCon");
		//for each chromosome calc number to flip
		for (String c: chromosomes){
			chromosome = c;

			//plus
			PointData[] pd = convertedPlusPointData.remove(chromosome);
			PointData converted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			pd = nonConvertedPlusPointData.remove(chromosome); 
			PointData nonConverted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			double numCon = converted.getInfo().getScoreTotal();
			double numNonCon = nonConverted.getInfo().getScoreTotal();
			System.out.println("Pre  +"+chromosome+"\t"+calcFrac4Line(numNonCon, numCon));
			
			int toFlip = (int)Math.round(numCon * plusCorr);
			
			double[] nonCon = flip(nonConverted, converted, toFlip);
			System.out.println("Post +"+chromosome+"\t"+calcFrac4Line(nonCon[0], nonCon[1]));
			
			postTotalNonCon+= nonCon[0];
			postTotalCon+= nonCon[1];
			
			//minus
			pd = convertedMinusPointData.remove(chromosome);
			converted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			pd = nonConvertedMinusPointData.remove(chromosome); 
			nonConverted = PointData.mergePointData(PointData.convertArray2ArrayList(pd), false, true);
			numCon = converted.getInfo().getScoreTotal();
			numNonCon = nonConverted.getInfo().getScoreTotal();
			System.out.println("Pre  -"+chromosome+"\t"+calcFrac4Line(numNonCon, numCon));
			
			toFlip = (int)Math.round(numCon * minusCorr);
			nonCon = flip(nonConverted, converted, toFlip);
			System.out.println("Post -"+chromosome+"\t"+calcFrac4Line(nonCon[0], nonCon[1]));
			
			postTotalNonCon+= nonCon[0];
			postTotalCon+= nonCon[1];
		}
		
		System.out.println("Pre  All\t" + calcFrac4Line(totalNonCon, totalCon));
		System.out.println("Post All\t" + calcFrac4Line(postTotalNonCon, postTotalCon));
		
	}
	
	
	public  double[] flip(PointData nonConverted, PointData converted, int toFlip) {
		//expand counts to positions so random picking is fair
		int[] exPos = PointData.expandCounts(converted);
		converted.nullPositionScoreArrays();
		Num.randomize(exPos, System.currentTimeMillis());
		
		int[] randomPositions = new int[toFlip];
		System.arraycopy(exPos, 0, randomPositions, 0, toFlip);
		
		//rebuild converted minus extracted positions
		int[] positions = new int[exPos.length - toFlip];
		System.arraycopy(exPos, toFlip, positions, 0, positions.length);
		Arrays.sort(positions);
	
		//flatten positions
		PointData pd = PointData.flattenPositions(positions);
		converted.setPositions(pd.getPositions());
		converted.setScores(pd.getScores());
		
		//write it
		converted.writePointData(convertedSaveDirectory);
		double numCon = converted.getInfo().getScoreTotal();
		
		//clean up
		positions = null;
		exPos = null;
		converted = null;
		
		// for Non Converted
		
		//now make PointData from the randomPositions
		Arrays.sort(randomPositions);
		float[] counts = new float[toFlip];
		Arrays.fill(counts, 1f);
		pd = new PointData();
		pd.setInfo(nonConverted.getInfo().copy());
		pd.setPositions(randomPositions);
		pd.setScores(counts);
		
		//merge with existing nonCon
		pd = PointData.mergePairedPointData(nonConverted, pd, true);
		double numNonCon = pd.getInfo().getScoreTotal();
		
		//write it
		pd.writePointData(nonConvertedSaveDirectory);
		pd = null;
		
		return new double[]{numNonCon, numCon};
	}
	


	public static String calcFrac4Line(double nonCon, double con){
		double f = calcFrac(nonCon, con);
		return (int)nonCon+"\t"+(int)con+ "\t"+Num.formatNumber(f, 4);
	}
	
	public static double calcFrac(double nonCon, double con){
		return nonCon / (nonCon+con);
	}

	/**Fetches the names of all the chromosomes where data is found for all*/
	public void fetchAllChromosomes(){
		System.out.println("Scanning chromosomes for all four bisulfite datasets...");
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = convertedPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = convertedMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = nonConvertedPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = nonConvertedMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		
		//for each chrom
		ArrayList<String> chroms = new ArrayList<String>();
		for (String chrom: c){
			if (convertedPlusPointData.containsKey(chrom) && convertedMinusPointData.containsKey(chrom) && nonConvertedPlusPointData.containsKey(chrom) && 
					nonConvertedMinusPointData.containsKey(chrom)) chroms.add(chrom);
			else {
				System.err.println("\tFailed to find all for datasets for "+chrom+", skipping.");
				convertedPlusPointData.remove(chrom);
				convertedMinusPointData.remove(chrom);
				nonConvertedPlusPointData.remove(chrom);
				nonConvertedMinusPointData.remove(chrom);
			}
		}
		chromosomes=  Misc.stringArrayListToStringArray(chroms);
		Arrays.sort(chromosomes);
	}


	/**Collects and calculates a bunch of stats re the PointData.*/
	private void fetchDataHashes(){
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);	
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BisSeqErrorAdder(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': convertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': nonConvertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'f': fractionNonConverted = Double.parseDouble(args[++i]); break;
					case 't': focusChromosome = args[++i]; break;
					case 's': saveDirectory = new File(args[++i]); break;
					
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//look for point directories
		if (convertedPointDirs == null || convertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your converted PointData directories(s)!\n");
		//only one directory look deeper
		if (convertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(convertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) convertedPointDirs = otherDirs;
		}

		//nonConverted data
		if (nonConvertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your nonConverted PointData directories(s)!\n");
		//only one directory look deeper
		if (nonConvertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(nonConvertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) nonConvertedPointDirs = otherDirs;
		}

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		else if (saveDirectory.exists() == false) saveDirectory.mkdirs();
		convertedSaveDirectory = new File (saveDirectory, "Converted");
		convertedSaveDirectory.mkdir();
		nonConvertedSaveDirectory = new File (saveDirectory, "NonConverted");
		nonConvertedSaveDirectory.mkdir();
		


		
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                BisSeqErrorAdder: June 2012                       **\n" +
				"**************************************************************************************\n" +
				"Takes PointData from converted and non-converted C bisulfite sequencing data parsed\n" +
				"using the NovoalignBisulfiteParser and simulates a worse non-coversion rate by \n" +
				"randomly picking converted observations and making them non-converted. This is\n" +
				"accomplished by first measuring the non-conversion rate in the test chromosome (e.g.\n" +
				"chrLambda), calculating the fraction of converted C's need to flip to non-converted\n" +
				"to reach the target fraction non-converted and then using this flip fraction\n" +
				"to modify the other chromosome data. \n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-c Converted PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories.\n" +
				"-n Non-converted PointData directories, ditto. \n" +
				"-f Target fraction non-converted for test chromosome, this cannot be less than the\n" +
				"       current fraction.\n"+
				"-t Test chromosome, defaults to chrLambda* .\n"+


				"\n"+

				"Example: java -Xmx12G -jar pathTo/USeq/Apps/BisSeqErrorAdder -c /Data/Sperm/Converted\n" +
				"      -n /Data/Sperm/NonConverted -f 0.02 \n\n" +

		"**************************************************************************************\n");

	}
}
