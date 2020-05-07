package edu.utah.seq.parsers.mpileup;

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
	private float minimumReadCoverage = 15.0f;
	private float maxError = 0.05f;
	private String versionedGenome;
	private File saveDirectory;
	private File baseFractionNonRef;
	private Gzipper failingBasesBed;
	
	//for window scanning
	private float maxFailing = 0.05f;
	private int windowSize = 50;
	private int minObs;
	
	//internal fields
	private Pattern rBase = Pattern.compile("[\\.,]");
	private int chromIndex = 0;
	private int positionIndex = 1;
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
		System.out.println(maxFailing + "\tMax fraction bp in win failing max error");
		System.out.println(minObs+"\tMinimum obs in window");
		
		System.out.println("\nProcessing chr");
		parseFile(pileupFile);
		try {
			failingBasesBed.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("\nWriting out high non ref regions "+uniRegions.size());
		IO.writeHashSet(uniRegions, new File(saveDirectory, "highNonRefRegions.bed"));
		
		new Bar2USeq(saveDirectory, true);
		
		System.out.println("\nHistogram of non-zero non-ref base frequencies:");
		nonRefFreqHistogram.setSkipZeroBins(false);
		nonRefFreqHistogram.printScaledHistogram();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" sec\n");
	}
	
	
	public void clearArrayLists(){
		positions.clear();
		reference.clear();
		nonReference.clear(); 
	}
	
	
	public void parseFile(File pileupFile){
		String line = null;
		try {
			BufferedReader in = IO.fetchBufferedReader(pileupFile);
			String chromosome = null;

			//for each line in the file
			while ((line = in.readLine()) != null){
				String[] tokens = Misc.TAB.split(line);
				if (tokens.length < 4) Misc.printErrAndExit("\nMalformed pileup line, skipping -> "+line+"\n\t"+tokens.length+" Tokens");
				
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
				
				Integer position = Integer.parseInt(tokens[positionIndex]) -1;
				
				//split by sample, first three indexes are chr, pos, refBase; then sets three for each bam: readDepth, baseString, baseQualString
				float numBases = 0f;
				float numRef = 0f;
				for (int i=3; i< tokens.length; i+=3){
					//check that there are three
					if ((i+2) < tokens.length){
						//fetch num bases
						numBases+= Float.parseFloat(tokens[i]);
						//count ref bases in baseCalls
						numRef+= countRefs(tokens[i+1]);
					}
				}
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
			
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing "+ line);
		}
	}
	
	private void processParsedData(ArrayList<Integer> positions, ArrayList<Float> refAL, ArrayList<Float> nonRefAL, String chromosome) throws IOException {
		String chr = chromosome;
		if (chr.startsWith("chr") == false) chr = "chr"+chromosome;
		
		int[] pos = Num.arrayListOfIntegerToInts(positions);
		if (pos.length == 0) return;
		float[] refCount = Num.arrayListOfFloatToArray(refAL);
		float[] nonRefCount = Num.arrayListOfFloatToArray(nonRefAL);
		
		//make a Point for bases where non refs are observed
		ArrayList<Point> pts = new ArrayList<Point>(10000);
		
		//make base fraction bases
		for (int i=0; i< pos.length; i++){
			float total = refCount[i]+ nonRefCount[i];
			if (total < minimumReadCoverage) continue;
			float fraction = nonRefCount[i]/ total;
			if (fraction !=0) nonRefFreqHistogram.count(fraction);
			if (fraction >= maxError) failingBasesBed.println(chr+"\t"+pos[i]+"\t"+(pos[i]+1)+"\t"+(int)nonRefCount[i]+"/"+(int)total+"\t"+fraction+"\t.");
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
					case 'e': maxError = Float.parseFloat(args[++i]); break;
					case 'f': maxFailing = Float.parseFloat(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
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
		
		//where to save base fraction edited
		baseFractionNonRef = new File (saveDirectory, "BaseFractionNonRef_" +(int)minimumReadCoverage +"RC");
		baseFractionNonRef.mkdir();
		
		//create per base bed with > maxError non ref
		File bad5 = new File (saveDirectory, "baseFracME"+maxError+"RC"+(int)minimumReadCoverage+".bed.gz");
		try {
			failingBasesBed = new Gzipper(bad5);
		} catch (Exception e) {
			e.printStackTrace();
		} 
		
		//set min obs
		minObs = (int)Math.round((float)windowSize * maxFailing);
		
	}	
	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              MpileUp Parser: April 2015                           **\n" +
				"**************************************************************************************\n" +
				"Parses a SAMTools mpileup output file for non reference bases generating bed files and\n" +
				"data tracks with information related to error prone bases. Multiple samples are merged.\n\n"+

				"Options:\n"+
				"-p Path to a mpileup file (.gz or.zip OK, use 'samtools mpileup -Q 20 -A -B *bam').\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-s Save directory, full path, defaults to pileup file directory.\n"+
				"-r Minimum read coverage, defaults to 15.\n"+
				"-e Max nonRef base fraction to score as a failing bp, defaults to 0.01\n"+
				"-w Window size, defaults to 50\n"+
				"-f Max fraction failing bp in window to score as a high error window, defaults to 0.05\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MpileupParser -p /Pileups/N2.mpileup.gz -v\n"+
				"      H_sapiens_Feb_2009 -e 0.1 -w 25\n\n" +

		"**************************************************************************************\n");

	}

}
