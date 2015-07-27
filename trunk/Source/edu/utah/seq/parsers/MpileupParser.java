package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.Info;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import edu.utah.seq.useq.apps.Bar2USeq;
import trans.tpmap.WindowMaker;
import util.gen.*;


/**@author davidnix*/
public class MpileupParser {

	//user defined fields
	private File pileupFile;
	private float minimumReadCoverage = 20.0f;
	private String versionedGenome;
	private File saveDirectory;
	private File nonReferencePointDataDirectory;
	private File referencePointDataDirectory;
	private File baseFractionNonRef;
	
	//for window scanning
	private float maxError = 0.05f;
	private float maxFailing = 0.1f;
	private int windowSize = 50;
	private int minObs;
	
	//internal fields
	private Pattern space = Pattern.compile("\\t");
	private Pattern rBase = Pattern.compile("[\\.,]");
	private int chromIndex = 0;
	private int positionIndex = 1;
	private int refseqIndex = 2;
	private int readDepthIndex = 3;
	private int baseCallIndex = 4;
	private int baseScoreIndex = 5;
	private ArrayList<Integer> positions = new ArrayList<Integer>();
	private ArrayList<Float> reference = new ArrayList<Float>();
	private ArrayList<Float> nonReference = new ArrayList<Float>();
	private WindowMaker windowMaker = new WindowMaker(windowSize, minObs);
	private LinkedHashSet<String> uniRegions = new LinkedHashSet<String>();
	public static Pattern INDEL = Pattern.compile("[-+](\\d+)");
	private Histogram nonRefFreqHistogram = new Histogram(0,1,100);

	//constructor
	public MpileupParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		//print thresholds
		System.out.println("Thresholds:");
		System.out.println((int)minimumReadCoverage+"\tMinimum Read Coverage");
		System.out.println(windowSize + "\tWindow size");
		System.out.println(maxError + "\tMax error");
		System.out.println(maxFailing + "\tMax % bp in win failing max error");
		System.out.println(minObs+"\tMinimum obs in window");
		
		System.out.println("\nProcessing pileup file");
		parseFile(pileupFile);
		
		System.out.println("\nWriting out high non ref regions "+uniRegions.size());
		IO.writeHashSet(uniRegions, new File(saveDirectory, "highNonRefRegions.bed"));
		
		new Bar2USeq(saveDirectory, true);
		
		System.out.println("\nHistogram of non-zero non-ref base frequencies:");
		nonRefFreqHistogram.setSkipZeroBins(false);
		nonRefFreqHistogram.printScaledHistogram();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	
	public void clearArrayLists(){
		positions.clear();
		reference.clear();
		nonReference.clear(); 
	}
	
	public void parseFile(File pileupFile){
		try {
			//input streams
			BufferedReader in = IO.fetchBufferedReader(pileupFile);
			
			//data output streams
			String chromosome = null;

			//for each line in the file
			String line;
			int badLines = 0;
			while ((line = in.readLine()) != null){
				if (line.length() == 0 || line.startsWith("#")) continue;
				String[] tokens = space.split(line);
				if (tokens.length < 4){
					if (badLines++ < 20) System.err.println("\nMalformed pileup line, skipping -> "+line+"\n\t"+tokens.length+" Tokens");
					continue;
				}
				//fetch num bases
				float numBases = Float.parseFloat(tokens[readDepthIndex]);
				if (numBases == 0.0) continue; 
				
				//set first chrom?
				if (chromosome == null) chromosome = tokens[chromIndex];
				
				//new chromosome?
				else if (tokens[chromIndex].equals(chromosome) == false){
					System.out.print(" "+chromosome);
					//process old
					processParsedData(positions, reference, nonReference, chromosome);
					//clear old
					clearArrayLists();
					chromosome = tokens[chromIndex];
				}
				
				//position
				Integer position = Integer.parseInt(tokens[positionIndex]) -1;

				//count ref bases in baseCalls
				float numRef = countRefs(tokens[baseCallIndex]);
				
				//save em
				positions.add(position);
				reference.add(numRef);
				nonReference.add(numBases-numRef);
			}
			//process final
			System.out.print(" "+chromosome);
			processParsedData(positions, reference, nonReference, chromosome);
			clearArrayLists();
			in.close();
			System.out.println();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	


	
	private void processParsedData(ArrayList<Integer> positions, ArrayList<Float> refAL, ArrayList<Float> nonRefAL, String chromosome) {
		String chr = chromosome;
		if (chr.startsWith("chr") == false) chr = "chr"+chromosome;
		
		int[] pos = Num.arrayListOfIntegerToInts(positions);
		if (pos.length == 0) return;
		float[] refCount = Num.arrayListOfFloatToArray(refAL);
		float[] nonRefCount = Num.arrayListOfFloatToArray(nonRefAL);
		
		//make a Point for bases where non refs are observed
		ArrayList<Point> pts = new ArrayList<Point>(10000);
		for (int i=0; i< pos.length; i++){
			if (nonRefCount[i] !=0) pts.add(new Point(pos[i], nonRefCount[i]));
		}
		if (pts.size() != 0){
			PointData pd = Point.extractPositionScores(pts);
			Info info = pd.getInfo();
			info.setName("NonRefBases");
			info.setVersionedGenome(versionedGenome);
			info.setChromosome(chr);
			info.setStrand(".");
			info.setReadLength(1);
			pd.writePointData(nonReferencePointDataDirectory);
		}

		//make edited bases
		pts.clear();
		for (int i=0; i< pos.length; i++){
			if (refCount[i] !=0) pts.add(new Point(pos[i], refCount[i]));
		}
		if (pts.size() != 0){
			PointData pd = Point.extractPositionScores(pts);
			Info info = pd.getInfo();
			info.setName("RefBases");
			info.setVersionedGenome(versionedGenome);
			info.setChromosome(chr);
			info.setStrand(".");
			info.setReadLength(1);
			pd.writePointData(referencePointDataDirectory);
		}
		
		//make base fraction bases
		pts.clear();
		for (int i=0; i< pos.length; i++){
			float total = refCount[i]+ nonRefCount[i];
			if (total < minimumReadCoverage) continue;
			float fraction = nonRefCount[i]/ total;
			if (fraction !=0) nonRefFreqHistogram.count(fraction);
			pts.add(new Point(pos[i], fraction));
		}
		
		if (pts.size() != 0){
			PointData pd = Point.extractPositionScores(pts);
			Info info = pd.getInfo();
			info.setName("NonRefBaseFraction");
			info.setVersionedGenome(versionedGenome);
			info.setChromosome(chr);
			info.setStrand(".");
			info.setReadLength(1);
			pd.writePointData(baseFractionNonRef);
			windowScan(pd, chr);
		}	
	}
	
	
	
	private void windowScan(PointData pd, String chr) {
		
		float[] frac = pd.getScores();
		int[] bpPosition = pd.getPositions();
		int[][] startStops = windowMaker.makeWindows(pd.getPositions());
		
		//for each 
		for (int[] ss : startStops){
			int numFailing = 0;
			int bpPosStartFail = -1;
			int bpPosStopFail = -1;

			for (int i=ss[0]; i<= ss[1]; i++){
				if (frac[i] > maxError) {
					numFailing++;
					//set start?
					if (bpPosStartFail == -1) bpPosStartFail = bpPosition[i];
					//set end?
					else if (bpPosition[i] > bpPosStopFail) bpPosStopFail = bpPosition[i];
				}
			}
			if (numFailing >= minObs) {
				uniRegions.add(chr+"\t"+ bpPosStartFail+"\t"+ (bpPosStopFail+1));
			}
		}
	}


	public float countRefs(String baseCalls){
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
		new MpileupParser(args);
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
		saveDirectory.mkdirs();
		
		nonReferencePointDataDirectory = new File (saveDirectory, "NonReference");
		nonReferencePointDataDirectory.mkdir();
		referencePointDataDirectory = new File (saveDirectory, "Reference");
		referencePointDataDirectory.mkdir();
		
		//where to save base fraction edited
		baseFractionNonRef = new File (saveDirectory, "BaseFractionNonRef_" +(int)minimumReadCoverage +"RC");
		baseFractionNonRef.mkdir();
		
		//set min obs
		minObs = (int)Math.round((float)windowSize * maxFailing);
		
	}	
	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Mpileup Parser: July 2015                           **\n" +
				"**************************************************************************************\n" +
				"Parses a SAMTools mpileup output file for non reference bases generating PointData for\n"+
				"the reference, non reference, and fraction non reference for bases that pass the\n"+
				"minimum read coverage filter.\n\n"+

				"Options:\n"+
				"-p Path to a mpileup file (.gz or.zip OK, use 'samtools mpileup -Q 13 -A -B' params).\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-s Save directory, full path, defaults to pileup file directory.\n"+
				"-r Minimum read coverage, defaults to 20.\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MpileupParser -p /Pileups/N2.mpileup.gz -v\n"+
				"      C_elegans_Oct_2010\n\n" +

		"**************************************************************************************\n");

	}

}
