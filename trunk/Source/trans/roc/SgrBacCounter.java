package trans.roc;
import java.io.*;
import java.util.*;

import util.gen.*;
/**
 * Application for scoring '.sgr' files for positives and negatives in relation to the BAC spike in data.
 */
public class SgrBacCounter {
	//fields
	private double start;
	private double stop;
	private double cutOff;
	private int numberBins = 102;
	private Bin[] bins;
	private Positive[] positives;
	private File[] sgrFiles;
	private int numberPositiveLines = 0;
	private int numberNegativeLines = 0;
	private Sgr[] sgrs;
	
	public SgrBacCounter(String[] args){
		if (args.length==0) {
			System.out.println("File/DirectoryOfSgrs, startScoreCutoff, stopScoreCutOff, cutOffToActuallyUseForPosNegCounting");
			System.exit(1);
		}
		//parse params
		start = Double.parseDouble(args[1]);
		stop = Double.parseDouble(args[2]);
		cutOff = Double.parseDouble(args[3]);
		sgrFiles = IO.extractFiles(new File(args[0]), ".sgr");
		
		//run thru each .sgr file
		//make bins and positives
		for (int x =0; x<sgrFiles.length; x++){
			System.out.println("\nProcessing: "+sgrFiles[x]);
			
			//read in sgr file
			sgrs = trans.misc.Util.loadSgrFile(sgrFiles[x]);
			
			//make bins and positives and scan
			makeBins();
			makePositives();
			scanSgrArray();
			
			
			//print report
			System.out.println("\n"+numberPositiveLines+" Num Positive Windows");
			System.out.println(numberNegativeLines+" Num Negative Windows\n");
			
			//bins
			double cumPos =0;
			double cumNeg =0;
			double cumTot =0;
			double fdr = 0;
			for (int i=bins.length-1; i>=0; i--){
				cumPos += bins[i].getNumPositives();
				cumNeg += bins[i].getNumNegatives();
				cumTot = cumPos + cumNeg;
				fdr = 100*cumNeg/cumTot;
				if (fdr > 0.5 && fdr< 3)System.out.println(bins[i]+"\t"+fdr+"\t*****");
				else System.out.println(bins[i]+"\t"+fdr);
			}
			
			//positives
			System.out.println("\nCutOff: "+cutOff+"\n");
			//print names
			for (int i=0; i<positives.length; i++){
				System.out.print(positives[i].getChromosome()+":"+positives[i].getStart()+"-"+positives[i].getStop()+"\t");
			}
			System.out.println();
			//print fraction above cut off
			for (int i=0; i<positives.length; i++){
				System.out.print(positives[i].getTotalNumberWindows()+"\t");
			}
			System.out.println();
			for (int i=0; i<positives.length; i++){
				System.out.print(fractionPositive(positives[i].getScores(), cutOff, positives[i].getTotalNumberWindows())+"\t");
			}
			
		}
	}
	
	/**Scans an sgr array, sorting into bins.*/
	public void scanSgrArray(){
		Bin bin = null;
		int numLines = sgrs.length;
		numberPositiveLines = 0;
		numberNegativeLines = 0;
		for (int i=0; i<numLines; i++){
			//find bin
			bin = findBin(sgrs[i]);
			//find whether it is a positive
			int pos = positive(positives, sgrs[i]);
			if (pos != -1) {
				bin.incrementPositives();
				positives[pos].getScores().add(new Double(sgrs[i].getScore()));
				positives[pos].setTotalNumberWindows(positives[pos].getTotalNumberWindows()+1);
				numberPositiveLines++;					
			}
			else {
				bin.incrementNegatives();
				numberNegativeLines++;
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
		new SgrBacCounter(args);
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
	
}
