package edu.utah.seq.data;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;

import edu.utah.seq.useq.data.*;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;

/**Simulates over dispersed RNA-Seq datasets.
 * @author Nix
 * 
 * 
#Created pool of SAM RNA-Seq alignment files lacking headers in one directory, these can be from different experiments but should have the same read length

#Filter for high quality alignments
for x in $(ls *gz)
do
echo $x
gunzip -c $x | awk '$5 >= 30' >> filtered.sam
done

#Randomize the file, if too big, split the file in to 3, randomize, split in 3, merge splits, then randomize each merged split
java -jar -Xmx30G ~/AppsUSeq/RandomizeTextFile -f filtered.sam

#Split it into 3, use wc -l to get total then split
wc -l filtered.sam
java -jar -Xmx10G ~/AppsUSeq/FileSplitter -f filtered_Randomized.sam -g -n 5691699

#Put split files in a directory
mkdir SplitFiles
mv *_filtered* SplitFiles/
ls SplitFiles/
1_filtered_Randomized.sam.gz  2_filtered_Randomized.sam.gz  3_filtered_Randomized.sam.gz

#Run parser on all files to extract PointData 
java -jar -Xmx10G ~/AppsUSeq/SAMParser -v H_sapiens_Feb_2009 -f SplitFiles -r CombinePointData -p 0 -q 120

#Run RNASeqSimulator with overdispersion or without
mkdir OverdispersedSimulation
cp hg19EnsGenesUnique.ucsc.gz OverdispersedSimulation/
java -jar -Xmx10G ~/AppsUSeq/RNASeqSimulator -u OverdispersedSimulation/hg19EnsGenesUnique.ucsc.gz -p CombinePointData/PointData/ -n SplitFiles/ > OverdispersedSimulation/overdispersed.log
gzip OverdispersedSimulation/*txt

mkdir NonOverdispersedSimulation
cp hg19EnsGenesUnique.ucsc.gz NonOverdispersedSimulation/
java -jar -Xmx10G ~/AppsUSeq/RNASeqSimulator -o -u NonOverdispersedSimulation/hg19EnsGenesUnique.ucsc.gz -p CombinePointData/PointData/ -n SplitFiles/ > NonOverdispersedSimulation/nonOverdispersed.log
gzip NonOverdispersedSimulation/*txt

#What your left with is 3 or 4 relica SAM datasets (t or c).  Add back the SAM headers and run them through your favorite RNASeq differential gene detector and have it predict lists of differentially expressed genes at a given FDR.  Then intersect these picks against the key derived from the log files.

 */
public class RNASeqSimulator {

	//fields
	private HashMap<String,UCSCGeneLine[]> geneModels;
	private int numberGenesToSkew = 500;
	private File[] pointDataDirs;
	private HashMap<String,PointData[]> plusPointData;
	private HashMap<String,PointData[]> minusPointData;
	private String chromosome;
	private PointData chromPlus = null;
	private PointData chromMinus = null;
	private int minimumNumberReads = 50;
	private int minimumSubExonReads = 50;
	private ArrayList<String> genesWithMinimumNumberReads = new ArrayList<String>();
	private ArrayList<String> genesWithReads = new ArrayList<String>();
	private HashMap<String, Double>genesToSkew = null;
	private File ucscGeneTableFile = null;
	private File[] samFiles;
	private boolean overdisperse = true;
	private boolean skipIntersectingGenes = false;

	//initial skew values for differentially expressed genes
	private double minimumSkewFraction = 0.2;
	private double maximumSkewFraction = 0.8;
	private double minimumExcludeFraction = 0.45;
	private double maximumExcludeFraction = 0.55;

	//alignment file fields
	private File workingFile;
	private int chromosomeColumnIndex = 2;
	private int positionColumnIndex = 3;
	private int sequenceColumnIndex = 9;
	private Pattern tab = Pattern.compile("\\t");
	private HashMap <String, DataOutputStream> chromOut;
	private HashMap <String, File> chromFiles; 
	private File tempDirectory;
	private PrintWriter tOut;
	private PrintWriter cOut;

	//constructor
	public RNASeqSimulator (String[] args){
		processArgs(args);
		if (args.length == 0) Misc.printExit("UCSCGeneTable PointData DirectoryWithNovoAlignmentFiles");
		//load gene models
		System.out.println("Loading gene models... ");
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(ucscGeneTableFile, 0);
		reader.splitByChromosome();
		geneModels = reader.getChromSpecificGeneLines();

		//randomly order the gene names
		UCSCGeneLine[] lines = reader.getGeneLines();
		String[] uniqueGeneNames = new String[lines.length];
		for (int i=0; i< uniqueGeneNames.length; i++) uniqueGeneNames[i] = lines[i].getName();
		Misc.randomize(uniqueGeneNames, System.currentTimeMillis());


		//fetch point data, don't load
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (pointDataDirs);
		plusPointData = PointData.convertArrayList2Array(combo[0]);
		minusPointData = PointData.convertArrayList2Array(combo[1]);

		//scan gene models for those with minimum number of reads
		System.out.println("\nScanning genes for min reads...");
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

		//print skewed genes
		System.out.println("\nDifferentially expressed genes with skew factor:\n"+genesToSkew);


		//reset skew values for modifying non differentially expressed genes
		minimumSkewFraction = 0.45;
		maximumSkewFraction = 0.55;
		minimumExcludeFraction = 1;
		maximumExcludeFraction = 0;

		//for each randomized alignment file
		System.out.println("\nParsing each SAM file and splitting reads to t or c datasets based on class and skew factor:");
		for (int i=0; i<samFiles.length; i++){
			//make hashes
			chromOut = new HashMap <String, DataOutputStream>();
			chromFiles = new HashMap <String, File>();

			//make temp dir
			tempDirectory = new File (ucscGeneTableFile.getParentFile(), "TempDir");
			tempDirectory.mkdir();

			//parse file
			workingFile = samFiles[i];
			parseFile();

			//close the binary writers
			closeWriters();

			//rescan gene models, this time skewing select genes according to the percents, others by a random amount.
			try {
				it = geneModels.keySet().iterator();
				File tFile = new File (ucscGeneTableFile.getParentFile(), "tReads"+i+".sam");
				File cFile = new File (ucscGeneTableFile.getParentFile(), "cReads"+i+".sam");
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
		System.out.println("\t"+workingFile);
		//split file to chromosome strand specific temp files
		boolean parsed = parseWorkingFile(); 
		if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");

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
					if (line.startsWith("@")) continue;
					tokens = tab.split(line);

					//parse chromosome, start, and stop
					String chr = tokens[chromosomeColumnIndex];
					int start = Integer.parseInt(tokens[positionColumnIndex]);
					int stop = start + tokens[sequenceColumnIndex].length();

					//get PrintWriter
					if (currentChrom.equals(chr) == false){
						currentChrom = chr;
						if (chromOut.containsKey(currentChrom)) dos = chromOut.get(currentChrom);
						else {
							File f = new File(tempDirectory, currentChrom);
							chromFiles.put(currentChrom, f);
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

	public void scanChromosome(){
		//fetch data
		if (fetchData() == false) return;
		//fetch positions 
		loadPointData();
		//scan
		scanGenes();
	}

	public void scanChromosomeBedFiles(){
		//get UCSCGeneLine[], these are sorted by start position and length, shortest first
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//load split RegionScoreText data, sorted by position and length, shortest first
		if (chromFiles.get(chromosome) == null) return;
		RegionScoreText[] alignments = RegionScoreText.oldLoadBinary_DEPRECIATED(chromFiles.get(chromosome), true);

		//for each gene 
		System.out.println("\n#Class\tGeneName\tChr\tStrand\tStart\tStop\tSkewFactor\tTotalReads\t#T\t#C\t#Int Genes");
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

				//skew the intersecting reads 
				RegionScoreText[] ia = new RegionScoreText[intersectingAlignments.size()];
				intersectingAlignments.toArray(ia);
				Misc.randomize(ia, System.currentTimeMillis());
				int numA = (int)Math.round((((double)ia.length) * skewFraction));
				System.out.println("DiffExp\t"+genes[i].getName()+"\t"+genes[i].coordinates()+"\t"+Num.formatNumber(skewFraction, 4)+ "\t"+ia.length +"\t"+numA+"\t"+(ia.length-numA)+"\t"+genes[i].getTss());
				for (int j=0; j< numA; j++) tOut.println(ia[j].getText());
				for (int j=numA; j< ia.length; j++) cOut.println(ia[j].getText());
			}

			//is it a gene with reads?
			else if (genesWithReads.contains(genes[i].getName())){
				//get skew fraction
				double skewFraction; 
				if (overdisperse) skewFraction = fetchSkewFraction();
				else skewFraction = 0.5;

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
				//System.out.println("ODed\t"+genes[i].getName()+"\t"+genes[i].coordinates()+"\t"+Num.formatNumber(skewFraction, 4)+ "\t"+ia.length +"\t"+numA+"\t"+(ia.length-numA));
				for (int j=0; j< numA; j++) tOut.println(ia[j].getText());
				for (int j=numA; j< ia.length; j++) cOut.println(ia[j].getText());
			}
		}
		
		//print out nongenic reads equally
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
			genes[i].setTss(numInts);
			if (skipIntersectingGenes && numInts >1) continue;

			//calculate totals
			ExonIntron[] exons = genes[i].getExons();
			//fetch scores
			float total = 0; 
			float[] exonCounts = new float[exons.length];
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				if (chromPlus != null) {
					float num = chromPlus.sumScoreBP(start, stop);
					total += num;  
					exonCounts[j] = num;
				}
				if (chromMinus != null){
					float num = chromMinus.sumScoreBP(start, stop);
					total += num;  
					exonCounts[j] += num;
					
				}
			}
			if (total >= minimumNumberReads) genesWithMinimumNumberReads.add(genes[i].getName());
			if (total > 1) genesWithReads.add(genes[i].getName());
		}
	}
	
	public int[] findExonsToSkew(float[] countsPerExon){
		//try 100 times
		Random r = new Random();
		int halfNum = (int) Math.round(((double)countsPerExon.length)/2.0);
		for (int i=0; i< 100; i++){
			//pick an exon to start
			int index = r.nextInt(countsPerExon.length);
			//pick number of exons 1 to 1/2 total
			int num = r.nextInt(halfNum);
			if (num > countsPerExon.length) num = countsPerExon.length;
			//count counts in exon set
			float total = 0;
			for (int j=index; j< num; j++){
				total += countsPerExon[j];
			}
			if (total > minimumSubExonReads){
				return new int[]{index, num};
			}
		}
		return null;
	}

	public double fetchSkewFraction(){
		double random = Math.random();
		double range = maximumSkewFraction - minimumSkewFraction;
		double skew = (random * range) + minimumSkewFraction;
		if (skew <= minimumExcludeFraction || skew >= maximumExcludeFraction) return skew;
		return fetchSkewFraction();
	}

	/**Fetch the data for a particular chromosome.*/
	public boolean fetchData(){
		//merge 
		chromPlus = null;
		if (plusPointData.containsKey(chromosome)) chromPlus = PointData.combinePointData(plusPointData.get(chromosome), true);
		chromMinus = null;
		if (minusPointData.containsKey(chromosome)) chromMinus = PointData.combinePointData(minusPointData.get(chromosome), true);
		if (chromPlus == null && chromMinus== null) return false;
		return true;
	}

	/**Loads and strips the scores.*/
	private void loadPointData(){
		//fetch data 
		if (chromPlus !=null) {
			chromPlus.getPositions();
			chromPlus.stripScores();
		}
		if (chromMinus!=null) {
			chromMinus.getPositions();
			chromMinus.stripScores();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new RNASeqSimulator(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File pointData = null;
		File samDir = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'u': ucscGeneTableFile = new File(args[++i]); break;
					case 'p': pointData = new File(args[++i]); break;
					case 'n': samDir = new File(args[++i]); break;
					case 'g': numberGenesToSkew = Integer.parseInt(args[++i]); break;
					case 'r': minimumNumberReads = Integer.parseInt(args[++i]); break;
					case 'a': minimumSkewFraction = Double.parseDouble(args[++i]); break;
					case 'b': maximumSkewFraction = Double.parseDouble(args[++i]); break;
					case 'c': minimumExcludeFraction = Double.parseDouble(args[++i]); break;
					case 'd': maximumExcludeFraction = Double.parseDouble(args[++i]); break;
					case 'o': overdisperse = false; break;
					case 's': skipIntersectingGenes = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (pointData == null || pointData.exists() == false) Misc.printErrAndExit("\nError: cannot find your PointData directory -> "+pointData);
		if (samDir == null || samDir.exists() == false) Misc.printErrAndExit("\nError: cannot find your SAM directory -> "+pointData);

		pointDataDirs = IO.extractFiles(pointData);
		samFiles = IO.extractFiles(samDir);

		System.out.println("Params:");
		System.out.println("\t"+ucscGeneTableFile.getName()+"\tGene table");
		System.out.println("\t"+pointData.getName()+"\tPointData directory");
		System.out.println("\t"+samDir.getName()+"\tSAM alignment directory");
		System.out.println("\t"+numberGenesToSkew+"\tNumber genes to skew");
		System.out.println("\t"+minimumNumberReads+"\tMinimum number of reads in a skewed gene");
		System.out.println("\t"+minimumSkewFraction+"\tMinimum skew fraction");
		System.out.println("\t"+maximumSkewFraction+"\tMaximum skew fraction");
		System.out.println("\t"+minimumExcludeFraction+"\tMinimum exclude/ overdispersed fraction");
		System.out.println("\t"+maximumExcludeFraction+"\tMaximum exclude/ overdispersed fraction");
		System.out.println("\t"+overdisperse+"\tOverdisperse non diff expressed genes");
		System.out.println("\t"+skipIntersectingGenes+"\tSkip intersecting genes");
		System.out.println();

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            RNA Seq Simulator: Aug 2011                           **\n" +
				"**************************************************************************************\n" +
				"RSS takes SAM alignment files from RNA-Seq data and simulates over dispersed, multiple\n" +
				"replica, differential, non-stranded RNA-Seq datasets. \n\n" +

				"Options:\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (name1 name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds)\n"+
				"-p PointData directories, full path, comma delimited. These should contain parsed\n" +
				"       PointData (chromosome specific xxx_-/+_.bar.zip files) from running the\n" +
				"       NovoalignParser on all of your novoaligned RNA-Seq data. \n" +
				"-n A full path directory name containing 3 or 4 equally split, randomized alignment\n" +
				"       xxx.sam (.zip or .gz) files. One for each replica you wish to simulate. Use the\n" +
				"       RandomizeTextFile and FileSplitter apps to generate these.\n" +
				"\n"+
				"Default Options:\n"+
				"-g Number of genes to make differentially expressed, defaults to 500\n"+
				"-r Minimum number of mapped reads to include a gene in the differential expression\n" +
				"       defaults to 50.\n"+
				"-a Smallest skew factor for differential expression, defaults to 0.2\n"+
				"-b Largest skew factor for differential expression, defaults to 0.8\n"+
				"-c Smallest excluded skew factor for differential expression and for overdispersion,\n"+
				"       defaults to 0.45\n"+
				"-d Largest excluded skew factor for differential expression and for overdispersion,\n"+
				"       defaults to 0.55\n"+
				"-o Don't overdisperse datasets, defaults to overdispersing data using -c and -d params.\n"+
				"-s Skip intersecting genes.\n"+
				"\n"+

				"Example: java -Xmx12G -jar pathTo/USeq/Apps/RNASeqSimulator -u \n" +
				"       /anno/hg19RefFlatKnownGenes.ucsc.txt -p /Data/Heart/MergedPointData/ -n\n" +
				"       /Data/Heart/SplitSAM/ -s 46 -r 15 -g 1000 \n\n" +

		"**************************************************************************************\n");

	}
}
