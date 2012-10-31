package util.apps;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;

import util.gen.*;

/**Calculates the area of Ts and Ts under Cs to estimate bisulfite converted bases. */
public class ABITraceTCPeakCalculator {

	//fields
	private File abiTraceFile = null;
	private boolean printOnlyTCs = false;
	private int smoothingWindowSize = 16;
	private double minimumFractionTForSmoothing = 0.9;
	private int maxExpandArea = 10;
	private int maxExpandMode = 5;
	private Pattern tab = Pattern.compile("\\t");
	private ConsensusBase[] consensusBases;

	private double[] allTScores;

	public ABITraceTCPeakCalculator (String[] args) {
		//process args
		processArgs(args);

		//load peak indexes
		//loadSeq2PeakIndexes();

		//load consensus bases
		loadConsensusBases();

		//calculate area of Ts for consensus T and C Bases
		calculateTAreas();

		//smooth T's where the T is > 90% any other int
		double[] ts = fetchAreaTs();
		double[] smoothedTs = Num.windowTrimmedMeanScoresNoShrinkIgnoreZeros(ts, smoothingWindowSize);
		

		//assign smoothedTs
		setSmoothedAreaTs(smoothedTs);

		//print results
		printResults();

	}

	/**Sets the smoothed t in each consensus base.*/
	public void calculateTAreas(){
		for (int i=0; i< consensusBases.length; i++){
			char consensusBase = consensusBases[i].consensusBase;
			if (consensusBase == 'T' || consensusBase == 'C') {
				//System.out.println(consensusBases[i]);
				consensusBases[i].calculateAreaOfT();
			}
		}
	}

	/**Sets the smoothed t in each consensus base.*/
	public void printResults(){
		try {
			File outputFile = new File ( abiTraceFile.getParentFile(), Misc.removeExtension(abiTraceFile.getName())+"_TCCalc.xls"  );
			System.out.println("\nWriting output to "+outputFile+"\n");
			PrintWriter out = new PrintWriter ( new FileWriter( outputFile));

			out.println("ConsensusIndex\tTraceIndex\tG\tA\tT\tC\tConsBase\tObsAreaT\tAveAreaT\t(AveT-ObsT)/AveT");
			if (printOnlyTCs){
				for (int i=0; i< consensusBases.length; i++){
					if (consensusBases[i].consensusBase == 'C' || consensusBases[i].consensusBase == 'T') out.println((i+1)+"\t"+consensusBases[i]);
				}
			}
			else {
				for (int i=0; i< consensusBases.length; i++){
					out.println((i+1)+"\t"+consensusBases[i]);
				}
			}
			out.close();
		}
		catch (Exception e){
			e.printStackTrace();
		}
	}

	/**Sets the smoothed t in each consensus base.*/
	public void setSmoothedAreaTs(double[] smoothedAreaTs){
		for (int i=0; i< consensusBases.length; i++){
			consensusBases[i].smoothedAreaT = smoothedAreaTs[i];
		}
	}

	/**Fetches the ts that pass the minimumFractionTForSmoothing*/
	public double[] fetchAreaTs(){
		double[] t = new double[consensusBases.length];
		for (int i=0; i< consensusBases.length; i++){
			t[i] = consensusBases[i].fetchAreaT();
		}
		return t;
	}



	/**Reads in the trace file and build the array of consensus bases*/
	public void loadConsensusBases(){
		try {
			BufferedReader in = new BufferedReader(new FileReader(abiTraceFile));
			String line = null;
			String[] tokens = null;
			//advance to first data line
			boolean foundStart = false;
			while ((line = in.readLine() ) != null){
				if (line.trim().startsWith("Seq")){
					foundStart = true;
					in.readLine();
					break;
				}
			}
			if (foundStart == false) Misc.printErrAndExit("\nError: Could not find the start of your data?\n");
			ArrayList<ConsensusBase> cb = new ArrayList<ConsensusBase>();
			ArrayList<Double> allTs = new ArrayList<Double>();
			//for each data line
			int index = 0;
			while ((line = in.readLine() ) != null){
				tokens = tab.split(line);
				//is it empty
				if (tokens[0].equals("")) break;
				//parse T
				int t = Integer.parseInt(tokens[4]);
				allTs.add(new Double(t));
				//is it a consensus
				if (tokens[1].equals(".") == false){
					char consensusBase = tokens[1].charAt(0);
					int g = Integer.parseInt(tokens[2]);
					int a = Integer.parseInt(tokens[3]);
					int c = Integer.parseInt(tokens[5]);
					ConsensusBase x = new ConsensusBase(index,g,a,t,c,consensusBase);
					//System.out.println(x);
					cb.add(x);
				}
				index++;
			}
			consensusBases = new ConsensusBase[cb.size()];
			cb.toArray(consensusBases);
			allTScores = Num.arrayListOfDoubleToArray(allTs);			
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	private class ConsensusBase{
		double g;
		double a;
		double t;
		double c;
		char consensusBase;
		double smoothedAreaT;
		double areaT;
		int index;

		private ConsensusBase(int index, int g, int a, int t, int c, char consensusBase){
			this.index = index;
			this.g = g;
			this.a = a;
			this.t = t;
			this.c = c;
			this.consensusBase = consensusBase;
		}

		public String toString(){
			if (consensusBase == 'C' || consensusBase == 'T') {

				String strSmoothedT = Num.formatNumber(smoothedAreaT, 1);

				double fractionT = (smoothedAreaT - areaT)/ smoothedAreaT;

				String strFractionT = Num.formatNumber(fractionT, 3);

				return index+"\t"+g+"\t"+a+"\t"+t+"\t"+c+"\t"+consensusBase+"\t"+areaT+"\t"+strSmoothedT+"\t"+strFractionT;
			}
			return index+"\t"+g+"\t"+a+"\t"+t+"\t"+c+"\t"+consensusBase;
		}

		public void calculateAreaOfT(){
			//find mode
			int modeIndex = findMode(index);
			double modeValue = allTScores[modeIndex];
			//System.out.println("Index "+index+" Mode of T "+modeIndex); 
			/*
			//expand symetrically from modeIndex until up tick in value or hit maxExpandArea
			int expand = 0;
			
			double oldUpValue = modeValue;
			double oldDownValue = modeValue;
			while (expand++ < maxExpandArea){
				int modeUp = modeIndex - expand;
				if (modeUp < 0) modeIndex = 0;
				int modeDown = modeIndex + expand;
				if (modeDown >= allTScores.length) modeDown = allTScores.length-1;
				double upValue = allTScores[modeUp];
				double downValue = allTScores[modeDown];
				if (upValue ==0 || downValue ==0){
					expand--;
					break;
				}
				if (upValue>= oldUpValue || downValue >= oldDownValue){
					expand--;
					break;
				}
				oldUpValue = upValue;
				oldDownValue = downValue;
			}*/
			//try just setting expand = 2, adjacent T's are truncating the area
			int expand = 2;
			//sum values
			areaT = 0;
			int start = modeIndex-expand;
			if (start <0 ) start = 0;
			int end = modeIndex+expand+1;
			if (end > allTScores.length) end = allTScores.length;
			//System.out.println("\tStart "+start+"\tEnd "+(stop-1));
			for (int i= start; i< end; i++) areaT+= allTScores[i];
		}

		public double fetchAreaT(){
			if (consensusBase != 'T') return 0;
			else {
				if ((t/c) < minimumFractionTForSmoothing) return 0;
				if ((t/a) < minimumFractionTForSmoothing) return 0;
				if ((t/g) < minimumFractionTForSmoothing) return 0;
				return areaT;
			}
		}
	}

	public int findMode(int index){
		//find current value
		double indexValue = allTScores[index];
		//fetch up one value
		int upIndex = index -1;
		if (upIndex < 0) upIndex =0;
		double upValue = allTScores[upIndex];
		//fetch down one value
		int downIndex = index +1;
		if (downIndex == allTScores.length) downIndex = allTScores.length-1;
		double downValue = allTScores[downIndex];
		//inspect
		//index is highest or equal
		if (indexValue>=upValue && indexValue >= downValue) return index;
		//both are greater
		if (indexValue<upValue && indexValue<downValue) {
			System.out.println("\nWarning, cannot find mode of base call.\n");
			return index;
		}
		//up is greater
		if (upValue> indexValue){
			//look back
			int stop = index - maxExpandMode;
			if (stop < 0) stop = 0;
			for (int i = upIndex-1; i>= stop; i--){
				//is new up less than prior
				if (allTScores[i] <= upValue) return upIndex;
				else {
					upValue = allTScores[i];
					upIndex = i;
				}
			}
			//still going up
			return index;
		}
		//by default, down is greater
		else {
			//look forward
			int stop = index + maxExpandMode;
			if (stop >= allTScores.length) stop = allTScores.length;
			for (int i=downIndex+1; i<stop; i++){
				if (allTScores[i]<= downValue) return downIndex;
				else{
					downValue = allTScores[i];
					downIndex = i;
				}
			}
			//still going up
			return index;
		}
	}

	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': abiTraceFile = new File (args[++i]);  break;
					case 'c': printOnlyTCs = true; break;
					case 'w': smoothingWindowSize = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		System.out.println();
		System.out.println("Trace file:\t"+abiTraceFile);
		System.out.println("Window size:\t"+smoothingWindowSize);
		System.out.println("Print TCs:\t"+printOnlyTCs);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        ABI Trace TC Peak Calculator: July 2009                   **\n" +
				"**************************************************************************************\n" +
				"Uses a sliding window to estimate the mean T peak area, compares it to the observed\n" +
				"T area for a given T or C to estimate the fraction T.  Useful for calculating the \n" +
				"fraction of converted Cs from bisufite treated DNA in a methylation experiment.\n\n"+

				"Required Parameters:\n"+
				"-f Full path file text for tab delimited text ABI trace file.\n" +
				"-w Window size in bp for estimating mean T peak areas, defaults to 16.\n"+
				"-c Print only T and C bases, defaults to all.\n"+

				"\n" +
				"Example: java -jar pathTo/Apps/ABITraceTCPeakCalculator -f /MyBisSeqData/exp1.txt -c\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new ABITraceTCPeakCalculator(args);
	}

}


