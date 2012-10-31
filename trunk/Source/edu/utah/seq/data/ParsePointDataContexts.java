package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.BarParser;

/**For parsing PointData for that meeting particular contexts
 * @author Nix
 * */
public class ParsePointDataContexts {

	//user defined fields
	private File[] pointDataDirectories;
	private File saveDirectory;
	private HashMap<String, File> chromosomeFastaFiles = new HashMap<String, File>();
	private Pattern contextPattern = null;
	private boolean sumCGContexts = false;

	//internal fields
	private String[] chromosomes;
	private String genomicSequence = null;
	private int genomicSequenceLengthMinus3 = 0;
	private HashMap<String,PointData[]> plusPointData;
	private HashMap<String,PointData[]> minusPointData;
	private long matches = 0;
	private long nonMatches = 0;
	private HashSet<String> matchingContexts = new HashSet<String>();
	private HashSet<String> nonMatchingContexts = new HashSet<String>();

	//by chromosome
	private String chromosome;
	private PointData mergedChromPlus = null;
	private PointData mergedChromMinus = null;	

	//constructors
	/**Stand alone.*/
	public ParsePointDataContexts(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//scan each chromosome
		scan();
		
		//stats
		System.out.println("\nParsing stats:");
		System.out.println("\t"+matches+"\tMatching positions\n\t"+matchingContexts);
		System.out.println("\n\t"+nonMatches+"\tNon matching positions\n\t"+nonMatchingContexts);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods
	public void scan(){
		loadStrandHashes();
		
		fetchAllChromosomes();
		
		//scan for pvalues
		System.out.print("\nScanning each chromosome... ");
		String oldChrom = "";
		for (int i=0; i< chromosomes.length; i++){		
			chromosome = chromosomes[i];	
			//check if present
			if (chromosomeFastaFiles.get(chromosome) == null) {
				System.err.println("\n\tWarning, could not find a fasta file for "+chromosome+", skipping!");
				continue;
			}
			
			//fetch the chrom sequence?
			if (oldChrom.equals(chromosome) ==false){
				File seqFile = chromosomeFastaFiles.get(chromosome);
				System.out.print(chromosome+" ");
				MultiFastaParser mfp = new MultiFastaParser(seqFile);				
				genomicSequence = mfp.getSeqs()[0].toUpperCase();
				genomicSequenceLengthMinus3 = genomicSequence.length()-3;
				oldChrom = chromosome;
			}
			
			//fetch data
			fetchDataAndRemove();
			//calculate base level stats
			parsePointData(mergedChromPlus, false);
			parsePointData(mergedChromMinus, true);
		}
		System.out.println();

	}


	/**Fetches the names of all the chromosomes in the data excluding lambda and phiX if present.*/
	private void fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = plusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = minusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		chromosomes=  Misc.hashSetToStringArray(c);
	}
	
	/**Fetchs the data for a particular chromosome.*/
	private void fetchDataAndRemove(){
		ArrayList<PointData> al = null;
		PointData[] pd;
		//merge converted
		mergedChromPlus = null;
		if (plusPointData.containsKey(chromosome)) {
			pd = plusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			mergedChromPlus = PointData.mergePointData(al, false, true);
		}
		mergedChromMinus = null;
		if (minusPointData.containsKey(chromosome)) {
			pd = minusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			mergedChromMinus = PointData.mergePointData(al, false, true);
		}
		pd = null;
		al = null;
	}

	private void parsePointData(PointData pd, boolean negativeStrand){
		if (pd == null) return;
		//fetch arrays
		int[] positions = pd.getPositions();
		float[] observations = pd.getScores();
		
		//make containers for those to save
		ArrayList<Integer> indexesToSave = new ArrayList<Integer>();

		//for each position, parse sequence and check context.
		for (int i=0; i< positions.length; i++){
			//watch out for out of bounds sequence due to partial matches to sequence termini
			if (positions[i] < 2 || positions[i] > genomicSequenceLengthMinus3) continue;
			//fetch genomic sequence
			String genSeq = genomicSequence.substring(positions[i]-2,positions[i]+3);
			if (negativeStrand) genSeq = Seq.reverseComplementDNA(genSeq);
			//does it match?
			if (contextPattern.matcher(genSeq).matches()){
				indexesToSave.add(new Integer(i));
				matches++;
				matchingContexts.add(genSeq);
			}
			else {
				nonMatches++;
				nonMatchingContexts.add(genSeq);
			}
		}
		
		//save parsed data?
		int num = indexesToSave.size();
		if (num !=0){
			int[] pos = new int[num];
			float[] scores = new float[num];
			for (int i=0; i< num; i++){
				int index = indexesToSave.get(i).intValue();
				pos[i] = positions[index];
				scores[i] = observations[index];
			}
			pd.setPositions(pos);
			pd.setScores(scores);
			HashMap<String,String> notes = pd.getInfo().getNotes();
			if (notes == null) notes = new HashMap<String,String>();
			notes.put(BarParser.DESCRIPTION_TAG, "PointData parsed for particular genomic contexts using the following regular expression '"+contextPattern.toString()+"'.");
			pd.getInfo().setNotes(notes);
			pd.writePointData(saveDirectory);
		}
	}
	
	/*private void parsePointDataSumCGs(PointData pdPlus, PointData pdMinus){
		if (pdPlus == null) return;
		//fetch arrays
		int[] positionsPlus = pdPlus.getPositions();
		float[] observationsPlus = pdPlus.getScores();
		
		//make containers for those to save
		ArrayList<Point> pointAL = new ArrayList<Point>();

		//for each position, parse sequence and check context.
		for (int i=0; i< positionsPlus.length; i++){
			//watch out for out of bounds sequence due to partial matches to sequence termini
			if (positionsPlus[i] < 2 || positionsPlus[i] > genomicSequenceLengthMinus3) continue;
			//fetch genomic sequence
			here
			String genSeq = genomicSequence.substring(positionsPlus[i]-2,positionsPlus[i]+3);
			
			//does it match?
			if (genSeq.equals("CG")){
				//look for scores on G
				int pos = positionsPlus[i]+1;
				float nextScore = pdMinus.sumScoreBP(pos, pos +1);
				pointAL.add(new Point(pos, nextScore+ observationsPlus[i]));
				matches++;
				matchingContexts.add(genSeq);
			}
			else {
				nonMatches++;
				nonMatchingContexts.add(genSeq);
			}
		}
		
		//save parsed data?
		int num = indexesToSave.size();
		if (num !=0){
			int[] pos = new int[num];
			float[] scores = new float[num];
			for (int i=0; i< num; i++){
				int index = indexesToSave.get(i).intValue();
				pos[i] = positionsPlus[index];
				scores[i] = observationsPlus[index];
			}
			pdPlus.setPositions(pos);
			pdPlus.setScores(scores);
			HashMap<String,String> notes = pdPlus.getInfo().getNotes();
			if (notes == null) notes = new HashMap<String,String>();
			notes.put(BarParser.DESCRIPTION_TAG, "PointData parsed for particular genomic contexts using the following regular expression '"+contextPattern.toString()+"'.");
			pdPlus.getInfo().setNotes(notes);
			pdPlus.writePointData(saveDirectory);
		}
	}*/


	/**Collects and calculates a bunch of stats re the PointData.*/
	private void loadStrandHashes(){
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (pointDataDirectories);
		plusPointData = PointData.convertArrayList2Array(combo[0]);
		minusPointData = PointData.convertArrayList2Array(combo[1]);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ParsePointDataContexts(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		File[] fastas = null;
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fastas = IO.extractFiles(new File(args[++i])); break;
					case 'p': pointDataDirectories = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'c': contextPattern = Pattern.compile(args[++i].toUpperCase()); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//look for pattern
		if (contextPattern == null) Misc.printErrAndExit("\nError: please enter a regular expression that matches a 5bp sequence.\n");
		
		//look for fasta files
		if (fastas == null || fastas.length ==0) Misc.printErrAndExit("\nError: cannot find any fasta sequence files?\n");

		//look for point directories
		if (pointDataDirectories == null || pointDataDirectories[0].isDirectory() == false) Misc.printExit("\nError: cannot find your PointData directories(s)!\n");
		//only one directory look deeper
		if (pointDataDirectories.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(pointDataDirectories[0]);
			if (otherDirs != null && otherDirs.length > 0) pointDataDirectories = otherDirs;
		}

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		else if (saveDirectory.exists() == false) saveDirectory.mkdirs();

		//load fasta files into hash
		chromosomeFastaFiles = new HashMap<String,File>();
		Pattern chrom = Pattern.compile("(.+)\\.fa.*");
		for (int i=0; i< fastas.length; i++){
			Matcher mat = chrom.matcher(fastas[i].getName());
			if (mat.matches()) chromosomeFastaFiles.put(mat.group(1), fastas[i]);
		}
	}	
	
	public static void printStatLine(double numerator, double denomenator, String name){
		double fraction = numerator/ denomenator;
		System.out.println( "\t"+ Num.formatNumber(fraction, 3)+"\t("+(long)numerator+"/"+(long)denomenator+")\t"+name);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           ParsePointDataContexts: Feb 2011                       **\n" +
				"**************************************************************************************\n" +
				"Parses PointData for particular 5bp genomic sequence contexts.\n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-p PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. These will be merged before splitting by summing overlapping\n" +
				"       position scores.\n" +
				"-f Fasta files for each chromosome.\n"+
				"-c Context java regular expression, must be 5bp long, 5'->3', case insensitive, e.g.:\n"+
				"       '..CG.' for CG\n"+
				"       '..C[CAT]G' for CHG\n"+
				"       '..C[CAT][CAT]' for CHH\n"+
				"       '..C[CAT].' for nonCG\n"+
				"       '..C[^G].' for nonCG\n"+
				"\n"+



				"\n"+

				"Example: java -Xmx12G -jar pathTo/USeq/Apps/ParsePointDataContexts -c '..CG.' -s\n" +
				"      /Data/PointData/CG -f /Genomes/Hg18/Fastas -p /Data/PointData/All/  \n\n" +

		"**************************************************************************************\n");

	}
}
