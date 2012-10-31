package trans.roc;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import trans.main.IntensityPrinter;
import util.gen.*;
/**
 * Application for scoring '.sgr' files for positives and negatives in relation to known (spike in) data.
 */
public class MethodTester {
	//fields
	private double start =0;
	private double stop =100;
	private double cutOff;
	private int numberBins = 102;
	private double targetFDR = 5;
	private File positiveRegionsFile = null;
	private Bin[] bins;
	private Positive[] positives;
	private File[] sgrFiles;
	private Sgr[] sgrLines;
	private int totalNumberSgrLinesInPositiveRegions;
	private int totalNumberSgrLinesNotInPositiveRegions; 
	private boolean zeroNaNScores = false;
	
	public MethodTester(String[] args){
		if (args.length==0) Misc.printExit("\nEnter a full path file/directory text containing 'xxx.sgr' files" +
				" and full path file text containing tab delimited positive regions (chrom, start, stop)\n");
		//parse params
		System.out.println("Thresholding at an FDR of "+targetFDR+"%");
		sgrFiles = IO.extractFiles(new File(args[0]), ".sgr");
		Arrays.sort(sgrFiles);
		zeroNaNScores = true;
		positiveRegionsFile = new File(args[1]);
		makePositives(positiveRegionsFile);
		System.out.println("Parsing "+positives.length+" positive regions -> "+positiveRegionsFile);
		
		//container for results
		String[] summaryLines = new String[sgrFiles.length];
		
		//run thru each .sgr file
		//make bins and positives
		for (int x =0; x<sgrFiles.length; x++){
			
			//load and check sgr file
			System.out.println("Parsing sgr file -> "+sgrFiles[x]);
			sgrLines = parseSgrFile(sgrFiles[x], zeroNaNScores);
			if (sgrLines == null){
				System.out.println("\nNot a number or infinity found, skipping file.\n");
				summaryLines[x] = sgrFiles[x].getName().substring(0, sgrFiles[x].getName().length()-4);
			}
			
			else{
				//threshold
				for (int i=0; i< 15; i++){
					System.out.println("\tBin and scoring from "+(float)start+ " to "+(float)stop);
					makeBins();
					//makePositives();
					blankPositives();
					binSgrs();
					calculateBinFDRs();
					if (resetStartStop() == false) break;
				}
				
				//print bins
				/*System.out.println("Printing Bins");
				 for (int j=bins.length-1; j>=0; j--){
				 System.out.println(bins[j]);
				 }*/
				
				cutOff = (start+stop)/2;
				
				//make summary line
				StringBuffer sb = new StringBuffer();
				sb.append(sgrFiles[x].getName().substring(0, sgrFiles[x].getName().length()-4));
				sb.append("\t");
				sb.append(cutOff);
				//for each region
				for (int i=0; i<positives.length; i++){
					sb.append("\t");
					sb.append(fractionPositive(positives[i].getScores(), cutOff, positives[i].getTotalNumberWindows()));
				}
				summaryLines[x] = sb.toString();
				
				//reset start stop
				start = 0;
				stop = 100;
			}
		}
		
		//print final report
		//print positive region names
		System.out.print("\nFile Name\tCutOff");
		for (int i=0; i<positives.length; i++){
			System.out.print("\t");
			System.out.print(positives[i].getChromosome()+":"+positives[i].getStart()+"-"+positives[i].getStop()+" ("+
					positives[i].getTotalNumberWindows()+")");
		}
		System.out.println();
		//print summary lines
		for (int i=0; i<summaryLines.length; i++){
			System.out.println(summaryLines[i]);
		}
		
		
	}
	
	/**Runs through all the bins calculating and setting fdrs for each.*/
	public void calculateBinFDRs(){
		//bins
		double cumPos =0;
		double cumNeg =0;
		double cumTot =0;
		for (int i=bins.length-1; i>=0; i--){
			cumPos += bins[i].getNumPositives();
			cumNeg += bins[i].getNumNegatives();
			cumTot = cumPos + cumNeg;
			double fdr = 100*cumNeg/cumTot;
			bins[i].setFdr(fdr);
		}
	}
	
	/**Runs through an array of Sgr[] noting whether it falls in the positive regions or not positive regions*/
	public void binSgrs(){
		totalNumberSgrLinesInPositiveRegions = 0;
		totalNumberSgrLinesNotInPositiveRegions = 0;
		double totalNeg =0;
		for (int i=0; i< sgrLines.length; i++){
			//find bin
			Bin bin = findBin(sgrLines[i]);
			//find whether it is in one of the positive regions
			int pos = positive(positives, sgrLines[i]);
			if (pos != -1) {
				bin.incrementPositives();
				positives[pos].getScores().add(new Double(sgrLines[i].getScore()));
				positives[pos].setTotalNumberWindows(positives[pos].getTotalNumberWindows()+1);
				totalNumberSgrLinesInPositiveRegions++;					
			}
			else {
				bin.incrementNegatives();
				totalNeg+= sgrLines[i].getScore();
				totalNumberSgrLinesNotInPositiveRegions++;
			}
		}
	}
	
	/**Checks to see if an sgr object is one of the positives, if not returns -1, otherwise returns index number.*/
	public static int positive(Positive[] pos, Sgr sgr){
		int numPos = pos.length;
		for (int i=0; i< numPos; i++){
			if (pos[i].matches(sgr)) return i;
		}
		return -1;
	}
	
	/**Finds an appropriate bin given an sgr object*/
	public Bin findBin(Sgr sgr){
		double score = sgr.getScore();
		//run thru bins
		for (int i=0; i<numberBins; i++){
			if (bins[i].contains(score))return bins[i];
		}
		System.out.println("\nERROR: no bin found for score "+score);
		System.out.println("Line: "+sgr.toString());
		System.exit(0);
		return null;
	}
	
	
	
	public static void main(String[] args){
		new MethodTester(args);
	}
	
	/**Given an ArrayList of Double, calculates the fraction that is >= scoreCutOff.*/
	public static double fractionPositive (ArrayList doubles, double scoreCutOff, double divider){
		double numValues = doubles.size();
		double counter = 0;
		double val;
		for (int i=0; i<numValues; i++){
			val = ((Double)doubles.get(i)).doubleValue();
			if ( val>= scoreCutOff ) counter++;
		}
		return counter/divider;
	}
	
	public void makeBins(){
		bins = new Bin[numberBins];
		double range = stop-start;
		double increment = range/(numberBins-2);
		int num = numberBins-1;
		//first bin and last bin are special
		bins[0]= new Bin(-1000000000,start);
		bins[num] = new Bin (stop, 1000000000);
		double oldInc = start;
		double newInc = start + increment;
		for (int i=1; i<num; i++) {
			bins[i]=new Bin(oldInc, newInc);
			oldInc += increment;
			newInc += increment;
		}
	}
	
	public void makePositives(){
		positives = new Positive[8];
		positives[0] = new Positive("chr3L", 14566109, 14745266 ); //1x
		positives[1] = new Positive("chr3R",12491763, 12640405);	//1x
		positives[2] = new Positive("chr2R",11316224, 11509485);	//4x
		positives[3] = new Positive("chr3R", 5128875, 5299312);		//4x
		positives[4] = new Positive("chr3L", 5474470, 5656680);		//10x
		positives[5] = new Positive("chr2R", 7672969, 7857589);		//10x
		positives[6] = new Positive("chr3L", 11803774, 11978404);	//20x
		positives[7] = new Positive("chr3R", 23311086, 23491584);	//20x
	}
	
	/**Resets the positive regions to their instantiation state.*/
	public void blankPositives(){
		for (int i=0; i<positives.length; i++){
			positives[i].setTotalNumberWindows(0);
			positives[i].getScores().clear();
		}
	}
	
	public void makePositives(File regions){
		try{
			ArrayList posAL = new ArrayList();
			String line;
			String[] tokens;
			BufferedReader in = new BufferedReader(new FileReader (regions));
			while ((line= in.readLine()) != null){
				line = line.trim();
				if (line.startsWith("#") || line.length()==0) continue;
				tokens = line.split("\\s+");
				posAL.add(new Positive(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2])));
			}
			positives = new Positive[posAL.size()];
			posAL.toArray(positives);
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Reads in an Sgr file. Chrom, position, score.
	 * Returns null if not a number or infinity is found and the zero flag is false.*/
	public static Sgr[] parseSgrFile(File sgrFile, boolean zeroNaNOrInfinityScores){
		ArrayList al = new ArrayList();
		String line;
		try{
			BufferedReader in = new BufferedReader(new FileReader(sgrFile));
			while ((line=in.readLine())!=null){
				if (line.trim().length()==0) continue;
				//make Sgr to hold line info
				Sgr sgr = new Sgr(line);
				if (Double.isNaN(sgr.getScore()) || Double.isInfinite(sgr.getScore())) {
					if (zeroNaNOrInfinityScores) sgr.setScore(0);
					else return null;
				}
				al.add(sgr);
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		Sgr[] s = new Sgr[al.size()];
		al.toArray(s);
		return s;
	}
	
	/**Scans the bins to determine the optimal start and stop for the next round.
	 * Returns true if start and stop were reset and another round should be run.
	 * Returns false if not reset.*/
	public boolean resetStartStop(){
		//find start score that gives an FDR just above target, set as start
		double newStart= -1;
		for (int i=bins.length-1; i>=0; i--){
			if (bins[i].getFdr() > targetFDR){
				newStart = bins[i].getStart();
				break;
			}
		}
		//find stop score that gives an FDR just below target, set as start
		double newStop=-1;
		for (int i=0; i<bins.length; i++){
			if (bins[i].getFdr() < targetFDR){
				newStop = bins[i].getStop();
				break;
			}
		}
		if (newStart == -1 || newStop == -1) {
			System.out.println("No reset start stop! Problem with initial bounds: newStart= "+newStart+ "  newStop= "+newStop );
			return false;
		}
		if ((float)newStop == (float)newStart ){
			System.out.println("No reset start stop! newStart == newStop "+(float)newStop);
			return false;
		}
		start = newStart;
		stop = newStop;
		return true;
	}
	
	/*
	public void processArgs(String[] args){
		File directory = null;
		File intensityDir = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': directory = new File(args[i+1]); i++; break;
					case 'i': intensityDir = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check to see if they entered required params
		if (directory==null || directory.isDirectory() == false){
			System.out.println("\nCannot find your bpmap directory!\n");
			System.exit(0);
		}
		if (intensityDir != null) intensityFiles = IO.extractFiles(intensityDir, "celp");
		else {
			System.out.println("\nCannot find your intensity files/ directory!\n");
			System.exit(0);
		}
		
		bpmapFile = new File(directory.toString(),"tpmap.fa"); 
		infoFile = new File(directory.toString(),"tpmap.faInfo");
		
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Intensity Printer: Jan 2006                              **\n" +
				"**************************************************************************************\n" +
				"IP prints to file  intensity scores in a .sgr text format for direct import into\n" +
				"Affy's IGB.\n" +
				"\n" +
				"Use the following options when running IP:\n\n" +
				"-b Full path directory text for the 'xxxTPMapFiles' generated by the TPMapProcessor\n" +
				"-i The full path director or file text for the serialized float[] intensity files.\n" +
				
				"\n" +
				"Example: java -Xmx256M trans/main/IntensityPrinter -b /affy/TPMapFiles/ -i \n" +
				"      /affy/t.celp \n" +
				"\n" +
				"You may need to increase the heap size of the java virtual machine using the -Xmx\n" +
				"flag. Adjust accordingly in response to out of memory errors and the available\n" +
				"resources. \n" +
				"\n" +
				"Questions, comments, suggestions? David_Nix@Affymetrix.com or SuperFly@lbl.gov\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length <4){
			printDocs();
			System.exit(0);
		}
		new IntensityPrinter(args);
	}
	*/
	
}
