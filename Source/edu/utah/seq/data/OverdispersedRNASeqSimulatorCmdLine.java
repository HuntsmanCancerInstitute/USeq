package edu.utah.seq.data;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;
import edu.utah.seq.useq.data.*;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;

/**Simulates over dispersed RNA-Seq datasets.*/
public class OverdispersedRNASeqSimulatorCmdLine {

	//fields
	private HashMap<String,UCSCGeneLine[]> geneModels;
	private int numberGenesToSkew = 500;
	private File[] treatmentPointDirs;
	private HashMap<String,PointData[]> treatmentPlusPointData;
	private HashMap<String,PointData[]> treatmentMinusPointData;
	private String chromosome;
	private PointData treatmentChromPlus = null;
	private PointData treatmentChromMinus = null;
	private int halfPeakShift = 0;
	private boolean useWeightedReads = false;
	private int minimumNumberReads = 20;
	private ArrayList<String> genesWithMinimumNumberReads = new ArrayList<String>();
	private ArrayList<String> genesWithReads = new ArrayList<String>();
	private HashMap<String, Double>genesToSkew = null;
	//initial skew values for differentially expressed genes
	private double minimumFraction = 0.2;
	private double maximumFraction = 0.8;
	private double minimumExcludeFraction = 0.45;
	private double maximumExcludeFraction = 0.55;

	//alignment file fields
	private File[] dataFiles;
	private File workingFile;
	private int chromosomeColumnIndex = 7;
	private int positionColumnIndex = 8;
	private int sequenceColumnIndex = 2;
	private int biggestIndex = positionColumnIndex;
	private Pattern spliceJunction = Pattern.compile("(.+)_(\\d+)_(\\d+)");
	private Pattern tab = Pattern.compile("\\t");
	private int radiusSpliceJunction = 34;
	private HashMap <String, DataOutputStream> chromOut;
	private HashMap <String, File> chromFiles; 
	private File tempDirectory;
	private PrintWriter tOut;
	private PrintWriter cOut;

	//constructor
	public OverdispersedRNASeqSimulatorCmdLine (String[] args){
		if (args.length == 0) Misc.printExit("UCSCGeneTable PointData DirectoryWithNovoAlignmentFiles");
		//load gene models
		System.out.print("Loading gene models... ");
		File geneFile = new File (args[0]);
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(geneFile, 0);
		reader.splitByChromosome();
		geneModels = reader.getChromSpecificGeneLines();

		//randomly order the gene names
		UCSCGeneLine[] lines = reader.getGeneLines();
		String[] uniqueGeneNames = new String[lines.length];
		for (int i=0; i< uniqueGeneNames.length; i++) uniqueGeneNames[i] = lines[i].getName();
		Misc.randomize(uniqueGeneNames, System.currentTimeMillis());


		//fetch point data, don't load
		treatmentPointDirs = IO.extractFiles(args[1]);
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (treatmentPointDirs);
		treatmentPlusPointData = PointData.convertArrayList2Array(combo[0]);
		treatmentMinusPointData = PointData.convertArrayList2Array(combo[1]);

		//scan gene models for those with minimum number of reads
		System.out.print("Scanning chromosomes for min reads...");
		Iterator<String> it = geneModels.keySet().iterator();
		while (it.hasNext()){
			chromosome = it.next();
			scanChromosome();
		}

		//randomly select the genes to skew and the skew value
		String[] loadedGenes = Misc.stringArrayListToStringArray(genesWithMinimumNumberReads);
		Misc.randomize(loadedGenes, System.currentTimeMillis());
		genesToSkew = new HashMap<String, Double>();
		for (int i=0; i< numberGenesToSkew; i++) {
			double skewFraction = fetchSkewFraction();
			genesToSkew.put(loadedGenes[i], new Double(skewFraction));
		}

		//reset skew values for modifying non differentially expressed genes
		minimumFraction = 0.45;
		maximumFraction = 0.55;
		minimumExcludeFraction = 1;
		maximumExcludeFraction = 0;

		//for each randomized alignment file
		dataFiles = IO.extractFiles(new File(args[2]));

		for (int i=0; i<dataFiles.length; i++){
			//make hashes
			chromOut = new HashMap <String, DataOutputStream>();
			chromFiles = new HashMap <String, File>();

			//make temp dir
			tempDirectory = new File (geneFile.getParentFile(), "TempDir");
			tempDirectory.mkdir();

			//parse file
			workingFile = dataFiles[i];
			parseFile();

			//close the binary writers
			closeWriters();

			//rescan gene models, this time skewing select genes according to the percents, others by a random amount.
			try {
				System.out.print("Scanning chromosomes to skew reads...");
				it = geneModels.keySet().iterator();
				File tFile = new File (geneFile.getParentFile(), "tReads"+i+".txt");
				File cFile = new File (geneFile.getParentFile(), "cReads"+i+".txt");
				tOut = new PrintWriter (new FileWriter(tFile));
				cOut = new PrintWriter (new FileWriter(cFile));
				while (it.hasNext()){
					chromosome = it.next();
					scanChromosomeBedFiles();
				}
				tOut.close();
				cOut.close();
			} catch (Exception e){
				e.printStackTrace();
			}
			//clean up
			IO.deleteDirectory(tempDirectory);
		}

	}

	/**Closes writers.*/
	public void closeWriters(){
		try{
			Iterator<DataOutputStream> it = chromOut.values().iterator();
			while (it.hasNext()) it.next().close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void parseFile(){
		System.out.print("\t"+workingFile);
		//split file to chromosome strand specific temp files
		boolean parsed = parseWorkingFile(); 
		if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
		System.out.println();
	}

	/**Splits a file by chromosome.*/
	public boolean parseWorkingFile(){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			String line;
			String[] tokens = null;
			int counter = 0;
			String currentChrom = "";
			DataOutputStream dos = null;
			while ((line = in.readLine()) !=null){
				try {
					if (line.startsWith("#")) continue;
					if (line.contains("chrAdapter")) continue;
					tokens = tab.split(line);
					if (tokens.length <= biggestIndex) continue;

					//parse chromosome, start, and stop
					String chr = tokens[chromosomeColumnIndex];
					int start = Integer.parseInt(tokens[positionColumnIndex]);
					int stop = start + tokens[sequenceColumnIndex].length();
					//check for splice junction
					Matcher mat = spliceJunction.matcher(chr);
					if (mat.matches()) {
						chr = mat.group(1);
						int startSplice = Integer.parseInt(mat.group(2));
						int stopSplice = Integer.parseInt(mat.group(3));
						start = startSplice - radiusSpliceJunction;
						stop = stopSplice + radiusSpliceJunction;
					}

					//get PrintWriter
					if (currentChrom.equals(chr) == false){
						currentChrom = chr;
						if (chromOut.containsKey(currentChrom)) dos = chromOut.get(currentChrom);
						else {
							//watch out for leading >
							String filtChr = currentChrom;
							if (filtChr.startsWith(">")) filtChr = filtChr.substring(1);
							//watch out for trailing .fasta and .fa
							if (filtChr.endsWith(".fasta")) filtChr = filtChr.substring(0, filtChr.length()-6);
							if (filtChr.endsWith(".fa")) filtChr = filtChr.substring(0, filtChr.length()-3);

							File f = new File(tempDirectory, filtChr);
							chromFiles.put(filtChr, f);
							dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f)));
							chromOut.put(currentChrom, dos);
						}
					}
					//save data
					//position
					dos.writeInt(start);
					dos.writeInt(stop);
					dos.writeFloat(0f);
					//length of line
					dos.writeInt(line.length());
					//line
					dos.writeBytes(line);
				} catch (Exception e){
					System.out.println("\nProblem parsing line -> "+line +" Skipping!\n");
					e.printStackTrace();
				}
			}
			in.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void scanChromosome(){
		//fetch data
		if (fetchData() == false) return;
		//fetch, shift, and merge all positions from the treatment and control
		shiftStripPointData();
		//scan
		scanGenes();
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void scanChromosomeBedFiles(){
		//get UCSCGeneLine[], these are sorted by start position and length, shortest first
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//load split RegionScoreText data, sorted by position and length, shortest first
		if (chromFiles.get(chromosome) == null) return;
		RegionScoreText[] alignments = RegionScoreText.oldLoadBinary_DEPRECIATED(chromFiles.get(chromosome), true);

		//for each gene 
		for (int i=0; i< genes.length; i++){

			//is this a gene to differentially skew reads?
			if (genesToSkew.containsKey(genes[i].getName())){
				//get skew fraction
				double skewFraction = genesToSkew.get(genes[i].getName()).doubleValue();

				//for each exon fetch the intersecting alignments
				ArrayList<RegionScoreText> intersectingAlignments = new ArrayList<RegionScoreText>();
				ExonIntron[] exons = genes[i].getExons();
				for (int j=0; j< exons.length; j++){
					int start = exons[j].getStart();
					int stop = exons[j].getEnd();
					//scan the alignments
					for (int x=0; x< alignments.length; x++){
						if (alignments[x].intersects(start, stop) && alignments[x].getScore() != -1){
							intersectingAlignments.add(alignments[x]);
							alignments[x].setScore(-1);
						}
					}
				}

				//split the intersecting reads
				RegionScoreText[] ia = new RegionScoreText[intersectingAlignments.size()];
				intersectingAlignments.toArray(ia);
				Misc.randomize(ia, System.currentTimeMillis());
				int numA = (int)Math.round((((double)ia.length) * skewFraction));
				System.out.println(genes[i].getName()+"\t"+genes[i].coordinates()+"\t"+Num.formatNumber(skewFraction, 4)+ "\tTotalReads "+ia.length +" #T "+numA+" #C "+(ia.length-numA));
				for (int j=0; j< numA; j++) tOut.println(ia[j].getText());
				for (int j=numA; j< ia.length; j++) cOut.println(ia[j].getText());
			}

			//is it a gene with reads?
			else if (genesWithReads.contains(genes[i].getName())){
				//get skew fraction
				double skewFraction = fetchSkewFraction();

				//for each exon fetch the intersecting alignments
				ArrayList<RegionScoreText> intersectingAlignments = new ArrayList<RegionScoreText>();
				ExonIntron[] exons = genes[i].getExons();
				for (int j=0; j< exons.length; j++){
					int start = exons[j].getStart();
					int stop = exons[j].getEnd();
					//scan the alignments
					for (int x=0; x< alignments.length; x++){
						if (alignments[x].intersects(start, stop) && alignments[x].getScore() != -1){
							intersectingAlignments.add(alignments[x]);
							alignments[x].setScore(-1);
						}
					}
				}

				//split the intersecting reads
				RegionScoreText[] ia = new RegionScoreText[intersectingAlignments.size()];
				intersectingAlignments.toArray(ia);
				Misc.randomize(ia, System.currentTimeMillis());
				int numA = (int)Math.round((((double)ia.length) * skewFraction));
				System.out.println("\tODed\t"+genes[i].getName()+"\t"+genes[i].coordinates()+"\t"+Num.formatNumber(skewFraction, 4)+ "\tTotalReads "+ia.length +" #T "+numA+" #C "+(ia.length-numA));
				for (int j=0; j< numA; j++) tOut.println(ia[j].getText());
				for (int j=numA; j< ia.length; j++) cOut.println(ia[j].getText());
			}
		}

		//strip out negative scoring alignments (those that fell within a gene)
		ArrayList<RegionScoreText> nonIAs = new ArrayList<RegionScoreText>();
		for (int i=0; i< alignments.length; i++){
			if (alignments[i].getScore() != -1) nonIAs.add(alignments[i]);
		}
		//randomize
		alignments = new RegionScoreText[nonIAs.size()];
		nonIAs.toArray(alignments);
		Misc.randomize(alignments, System.currentTimeMillis());
		//print 1/2 to tOut and 1/2 to cOut
		int half = (int)Math.round(((double)alignments.length)/ 2.0);	
		for (int i=0; i< half; i++) tOut.println(alignments[i].getText());
		for (int i=half; i< alignments.length; i++) cOut.println(alignments[i].getText());	
	}

	private void scanGenes(){
		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){

			//does it intersect another gene?
			int st = genes[i].getTxStart();
			int sp = genes[i].getTxEnd();
			int numInts = 0;
			for (int x = 0; x< genes.length; x++){
				if (genes[x].intersects(st, sp)) numInts++;
			}
			if (numInts >1) continue;

			//calculate totals
			ExonIntron[] exons = genes[i].getExons();
			//fetch scores
			float tSumPlus = 0; 
			float tSumMinus = 0;
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				if (treatmentChromPlus != null) tSumPlus += treatmentChromPlus.sumScoreBP(start, stop); 
				if (treatmentChromMinus != null)tSumMinus += treatmentChromMinus.sumScoreBP(start, stop); 
			}
			float tSum = tSumPlus+ tSumMinus;
			if (tSum >= minimumNumberReads) genesWithMinimumNumberReads.add(genes[i].getName());
			if (tSum > 1) genesWithReads.add(genes[i].getName());
		}
	}

	public double fetchSkewFraction(){
		double random = Math.random();
		double range = maximumFraction - minimumFraction;
		double skew = (random * range) + minimumFraction;
		if (skew <= minimumExcludeFraction || skew >= maximumExcludeFraction) return skew;
		return fetchSkewFraction();
	}

	/**Fetchs the data for a particular chromosome.*/
	public boolean fetchData(){
		//merge treatment
		treatmentChromPlus = null;
		if (treatmentPlusPointData.containsKey(chromosome)) treatmentChromPlus = PointData.combinePointData(treatmentPlusPointData.get(chromosome), true);
		treatmentChromMinus = null;
		if (treatmentMinusPointData.containsKey(chromosome)) treatmentChromMinus = PointData.combinePointData(treatmentMinusPointData.get(chromosome), true);
		if (treatmentChromPlus == null && treatmentChromMinus== null) return false;
		return true;
	}

	/**Shifts the positions halfPeakShift (+ for sense, - for antisense) sets the positions into the data.
	 * May replace all scores with 1 if stripScores == true.*/
	private void shiftStripPointData(){
		//fetch data from treatments
		if (treatmentChromPlus !=null) {
			int[] p = treatmentChromPlus.getPositions();
			addShift(p,halfPeakShift);
			if (useWeightedReads == false) treatmentChromPlus.stripScores();
		}
		if (treatmentChromMinus!=null) {
			int[] p = treatmentChromMinus.getPositions();
			addShift(p, -1*halfPeakShift);
			if (useWeightedReads == false) treatmentChromMinus.stripScores();
		}
	}

	/**Adds the toAdd to each int.*/
	public static void addShift(int[] positions, int toAdd){
		for (int i=0; i< positions.length; i++){
			positions[i] += toAdd;
			if (positions[i]<0) positions[i] = 0;
		}
	}

	public static void main (String[] args){
		new OverdispersedRNASeqSimulatorCmdLine (args);
	}
}
