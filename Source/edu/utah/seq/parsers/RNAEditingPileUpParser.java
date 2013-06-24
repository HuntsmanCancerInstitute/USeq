package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.ComparatorPointPosition;
import edu.utah.seq.data.Info;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import util.gen.*;


/**@author davidnix*/
public class RNAEditingPileUpParser {

	//user defined fields
	private File pileupFile;
	private float minimumReadCoverage = 5.0f;
	private String versionedGenome;
	private File saveDirectory;
	private File nonConvertedPointDataDirectory;
	private File convertedPointDataDirectory;
	private boolean nonStrandedParsing = true;
	
	//internal fields
	private Pattern space = Pattern.compile("\\t");
	private Pattern GBase = Pattern.compile("G");
	private Pattern cBase = Pattern.compile("c");
	private Pattern rBase = Pattern.compile("[\\.,]");
	private Pattern rBasePlus = Pattern.compile("\\.");
	private Pattern rBaseMinus = Pattern.compile(",");
	private int chromIndex = 0;
	private int positionIndex = 1;
	private int refseqIndex = 2;
	//private int readDepthIndex = 3;
	private int baseCallIndex = 4;
	//private int baseScoreIndex = 5;
	private File baseFractionEdited;
	private ArrayList<Integer> positionsPlus = new ArrayList<Integer>();
	private ArrayList<Float> conversionsPlus = new ArrayList<Float>();
	private ArrayList<Float> nonConversionsPlus = new ArrayList<Float>();
	private ArrayList<Integer> positionsMinus = new ArrayList<Integer>();
	private ArrayList<Float> conversionsMinus = new ArrayList<Float>();
	private ArrayList<Float> nonConversionsMinus = new ArrayList<Float>();
	public static Pattern INDEL = Pattern.compile("[-+](\\d+)");

	public RNAEditingPileUpParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		System.out.println("Processing pileup file");
		parseFile(pileupFile);
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	
	public void clearArrayLists(){
		positionsPlus.clear();
		conversionsPlus.clear();
		nonConversionsPlus.clear(); 
		positionsMinus.clear();
		conversionsMinus.clear();
		nonConversionsMinus.clear();
	}
	
	public void parseFile(File pileupFile){
		try {
			//stranded analysis?
			if (nonStrandedParsing){
				rBasePlus = rBase;
				rBaseMinus = rBase;
			}
			
			//input streams
			BufferedReader in = IO.fetchBufferedReader(pileupFile);
			
			//data output streams
			String chromosome = null;

			
			//for each line in the file
			String line;
			while ((line = in.readLine()) != null){
				if (line.length() == 0 || line.startsWith("#")) continue;
				String[] tokens = space.split(line);
				if (tokens.length < 5){
					System.err.println("\nMalformed pileup line, skipping -> "+line+"\n\t"+tokens.length+" Tokens");
					continue;
				}
				//set first chrom?
				if (chromosome == null) chromosome = tokens[chromIndex];
				
				//new chromosome?
				else if (tokens[chromIndex].equals(chromosome) == false){
					System.out.print(" "+chromosome);
					//process old
					processParsedData(positionsPlus, conversionsPlus, nonConversionsPlus, chromosome, "+");
					processParsedData(positionsMinus, conversionsMinus, nonConversionsMinus, chromosome, "-");
					//clear old
					clearArrayLists();
					chromosome = tokens[chromIndex];
				}
				
				//refseq base 
				String refSeqBase = tokens[refseqIndex].toUpperCase();
				
				//plus strand
				if (refSeqBase.equals("A")){
					Integer position = Integer.parseInt(tokens[positionIndex]) -1;
					//remove INDELS
					tokens[baseCallIndex] = removeINDELS(tokens[baseCallIndex]);
					//look for upper case G's, (conversion of A to G on plus strand)
					float numGs = countGs(tokens[baseCallIndex]);
					//look for periods ....
					float numNonGs = countRefs(tokens[baseCallIndex], rBasePlus);
					//save em
					positionsPlus.add(position);
					conversionsPlus.add(numGs);
					nonConversionsPlus.add(numNonGs); 
				}
				
				//minus strand,
				else if (refSeqBase.equals("T")){
					Integer position = Integer.parseInt(tokens[positionIndex]) -1;
					//remove INDELS
					tokens[baseCallIndex] = removeINDELS(tokens[baseCallIndex]);
					//look for lower case c's, (conversion of A to G on minus strand)
					float numCs = countCs(tokens[baseCallIndex]);
					//look for commas ,,,
					float numNonCs = countRefs(tokens[baseCallIndex], rBaseMinus);
					//save em
					positionsMinus.add(position);
					conversionsMinus.add(numCs);
					nonConversionsMinus.add(numNonCs); 
				}
				
			}
			//process final
			System.out.print(" "+chromosome);
			processParsedData(positionsPlus, conversionsPlus, nonConversionsPlus, chromosome, "+");
			processParsedData(positionsMinus, conversionsMinus, nonConversionsMinus, chromosome, "-");
			clearArrayLists();
			in.close();
			System.out.println();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static String removeINDELS(String bases){
		Matcher mat = INDEL.matcher(bases);
		int numDeleted = 0;
		while (mat.find()){
			//System.out.println("\nPre\t"+bases);
			int firstMatch = mat.start() - numDeleted;
			int lastMatch = mat.end() - numDeleted;
			int numToDelete = Integer.parseInt(mat.group(1));
			//System.out.println("\t"+numDeleted+"\t"+firstMatch+"\t"+lastMatch);
			bases = bases.substring(0, firstMatch) + bases.substring(lastMatch+numToDelete);
			//System.out.println("Post\t"+bases);
			numDeleted += numToDelete;
			numDeleted += (lastMatch-firstMatch);
		}
		return bases;
	}
	
	private void processParsedData(ArrayList<Integer> positions, ArrayList<Float> conversions, ArrayList<Float> nonConversions, String chromosome, String strand) {
		int[] pos = Num.arrayListOfIntegerToInts(positions);
		if (pos.length == 0) return;
		float[] con = Num.arrayListOfFloatToArray(conversions);
		float[] nonCon = Num.arrayListOfFloatToArray(nonConversions);
		
		//make non edited bases
		ArrayList<Point> pts = new ArrayList<Point>(10000);
		for (int i=0; i< pos.length; i++){
			if (nonCon[i] !=0) pts.add(new Point(pos[i], nonCon[i]));
		}
		if (pts.size() != 0){
			PointData pd = Point.extractPositionScores(pts);
			Info info = pd.getInfo();
			info.setName("NonEditedBases");
			info.setVersionedGenome(versionedGenome);
			info.setChromosome(chromosome);
			info.setStrand(strand);
			info.setReadLength(1);
			pd.writePointData(nonConvertedPointDataDirectory);
		}

		//make edited bases
		pts.clear();
		for (int i=0; i< pos.length; i++){
			if (con[i] !=0) pts.add(new Point(pos[i], con[i]));
		}
		if (pts.size() != 0){
			PointData pd = Point.extractPositionScores(pts);
			Info info = pd.getInfo();
			info.setName("EditedBases");
			info.setVersionedGenome(versionedGenome);
			info.setChromosome(chromosome);
			info.setStrand(strand);
			info.setReadLength(1);
			pd.writePointData(convertedPointDataDirectory);
		}
		
		//make base fraction bases
		pts.clear();
		for (int i=0; i< pos.length; i++){
			float total = con[i]+ nonCon[i];
			if (total < minimumReadCoverage) continue;
			float fraction = con[i]/ total;
			pts.add(new Point(pos[i], fraction));
		}
		
		if (pts.size() != 0){
			PointData pd = Point.extractPositionScores(pts);
			Info info = pd.getInfo();
			info.setName("BaseFractionEdited");
			info.setVersionedGenome(versionedGenome);
			info.setChromosome(chromosome);
			info.setStrand(strand);
			info.setReadLength(1);
			pd.writePointData(baseFractionEdited);
		}
		
	}

	public float countGs(String baseCalls){
		Matcher mat = GBase.matcher(baseCalls);
		float num = 0;
		while (mat.find()) num++;
		return num;
	}
	
	public float countCs(String baseCalls){
		Matcher mat = cBase.matcher(baseCalls);
		float num = 0;
		while (mat.find()) num++;
		return num;
	}
	
	public float countRefs(String baseCalls, Pattern rBase){
		Matcher mat = rBase.matcher(baseCalls);
		float num = 0;
		while (mat.find()) num++;
		return num;
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new RNAEditingPileUpParser(args);
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
					case 'p': pileupFile = new File(args[++i]); break;
					case 'r': minimumReadCoverage = Float.parseFloat(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); saveDirectory.mkdir(); break;
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 't': nonStrandedParsing = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for bam files
		if (pileupFile == null || pileupFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your pileup file?\n");
		if (versionedGenome == null) Misc.printExit("\nPlease enter a genome version recognized by UCSC, see http://genome.ucsc.edu/FAQ/FAQreleases.\n");
		
		if (saveDirectory == null ) {
			String baseName = Misc.capitalizeFirstLetter(Misc.removeExtension(pileupFile.getName()));
			saveDirectory = new File (pileupFile.getParentFile(), baseName);
		}
		String strand = "Stranded_";
		if (nonStrandedParsing) strand = "NonStranded_";
		File pointDataDirectory = new File (saveDirectory, strand+"PointData");
		pointDataDirectory.mkdirs();
		nonConvertedPointDataDirectory = new File (pointDataDirectory, "Reference");
		nonConvertedPointDataDirectory.mkdir();
		convertedPointDataDirectory = new File (pointDataDirectory, "Edited");
		convertedPointDataDirectory.mkdir();
		
		
		//where to save base fraction edited
		baseFractionEdited = new File (saveDirectory, strand+"BaseFractionEdited_" +(int)minimumReadCoverage +"RC");
		baseFractionEdited.mkdir();
		
		
		
		System.out.println(minimumReadCoverage+"\tMinimum Read Coverage");
		
		
	}	
	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          RNA Editing PileUp Parser: June 2013                    **\n" +
				"**************************************************************************************\n" +
				"Parses a SAMTools mpileup output file for refseq A bases that show evidence of\n" +
				"RNA editing via conversion to Gs, stranded. Base fraction editing is calculated for\n" +
				"bases passing the thresholds for viewing in IGB and subsequent clustering with\n" +
				"the RNAEditingScanSeqs app. The parsed PointData can be further processed using the\n" +
				"methylome analysis applications.\n\n"+

				"Options:\n"+
				"-p Path to a mpileup file (.gz or.zip OK, use 'samtools mpileup -Q 13 -A -B' params).\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-s Save directory, full path, defaults to pileup file directory.\n"+
				"-r Minimum read coverage, defaults to 5.\n"+
				"-t Generate stranded specific reference calls, defaults to non stranded. Required for\n" +
				"      stranded down stream analysis.\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/RNAEditingPileUpParser -t -p \n" +
				"      /Pileups/N2.mpileup.gz -v C_elegans_Oct_2010\n\n" +

		"**************************************************************************************\n");

	}

}
