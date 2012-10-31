package util.gen;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import edu.utah.seq.data.*;
import util.bio.parsers.*;

import edu.utah.seq.data.ComparatorSmoothingWindowScore;

/**For generating fuzzy non discrete binomial pvalues for trials with less than 21 observations.
 * @author Nix*/
public class FuzzyBinomialPValueGenerator {

	//fields
	private File tempDirectory;
	private File fullPathToR;
	private double expect;
	private double[][] pbinomTC;
	private double[][] dbinomTC;
	private int maxNumOb = 21;
	private Random random = new Random(0);
	private ArrayList<AdjPValQVal> alPValQVal = new ArrayList<AdjPValQVal>();
	private float[] adjPVals;
	private float[] qVals;
	private double maxPVal = 4.940666e-324;
	private double maxPValLogTrans = 3233.666;
	

	public FuzzyBinomialPValueGenerator(File tempDirectory, File fullPathToR, double expect){
		this.tempDirectory = tempDirectory;
		this.fullPathToR = fullPathToR;
		this.expect = expect;
		makeBinomialPValParts();
	}

	/**Returns the bionomial p-value with a bit of noise to make the distribution less discrite.
	 * Assumes the total obs < 21. Note there is a top stop to the pvals.*/
	public double fetchFuzzyPVal(int t, int c, boolean minusLog10Convert){
		//note, random returns 0 inclusive values so need to flip.
		double rand = random.nextDouble();
		double p = 1 - ( pbinomTC[t][c] + ( dbinomTC[t][c] * (1-rand)));
		//watch for pvals <= 0
		if (p<= 0){
			if (minusLog10Convert) return maxPValLogTrans;
			return maxPVal;
		} 
		if (minusLog10Convert) p = Num.minus10log10(p);
		return p;
	}

	/**Returns the binomial p-value, no noise.
	 * Assumes the total obs < 21.*/
	public double fetchPVal(int t, int c, boolean minusLog10Convert){
		double p = 1 - ( pbinomTC[t][c] + ( dbinomTC[t][c]));
		//watch for pvals <= 0
		if (p<= 0){
			if (minusLog10Convert) return maxPValLogTrans;
			return maxPVal;
		} 
		if (minusLog10Convert) p = Num.minus10log10(p);
		return p;
	}

	/**Adds a new AdjPVal to QVal point.*/
	public void addAdjPValQVal(float adjPVal, float qVal){
		alPValQVal.add(new AdjPValQVal(adjPVal, qVal));
	}
	
	/**Adds a new AdjPVal to QVal point.*/
	public void addAdjPValQVal(String adjPVal, String qVal) throws NumberFormatException{
		alPValQVal.add(new AdjPValQVal(Float.parseFloat(adjPVal), Float.parseFloat(qVal)));
	}


	/**Call after generatingSortedPAndQValArrays. Assumes scores are in -10Log10(value) space*/
	public void loadQValues(int index, SmoothingWindow[] toFixQVal){
		//sort windows on up pValue
		ComparatorSmoothingWindowScore comp = new ComparatorSmoothingWindowScore(index);
		Arrays.sort(toFixQVal, comp);
		//for every window pValue, look up closest adjPVal and it's associated qVal
		float currPVal = -10000000;
		float currQVal = 0; 		
		for (int i=0; i< toFixQVal.length; i++){
			//scores = pVal,qVal,empFDR,tSumPlus,tSumMinus,cSumPlus,cSumMinus
			float[] scores = toFixQVal[i].getScores();
			float testPVal = scores[index];			
			if (currPVal != testPVal) {
				currPVal = testPVal;
				//look up qVal, if exact match use that qVal otherwise take smaller one
				int mark = Arrays.binarySearch(adjPVals, currPVal);
				//not exact therefore interpolate
				if (mark < 0) {
					if (qVals.length ==1){
						currQVal = qVals[0];
					}
					else {
						mark = (-1 * mark) -2;
						if (mark < 0) mark = 0;
						int markPlusOne = mark +1;
						//currQVal = Num.interpolateY(adjPVals[mark], qVals[mark], adjPVals[markPlusOne], qVals[markPlusOne], currPVal);
						double trans = Num.interpolateY(Num.antiNeg10log10(adjPVals[mark]), Num.antiNeg10log10(qVals[mark]), Num.antiNeg10log10(adjPVals[markPlusOne]), Num.antiNeg10log10(qVals[markPlusOne]), Num.antiNeg10log10(currPVal));
						currQVal = Num.minus10log10Float(trans);
					}
				}
				else currQVal = qVals[mark];

			}
			scores[index] = currQVal;
		}
	}

	/**Call after generatingSortedPAndQValArrays. Assumes scores are in -10Log10(value) space*/
	public void loadQValues(int index, UCSCGeneLine[] toFixQVal){
		
		//sort windows on up pValue
		UCSCGeneLineComparatorScore comp = new UCSCGeneLineComparatorScore(index);
		Arrays.sort(toFixQVal, comp);
		//for every window pValue, look up closest adjPVal and it's associated qVal
		float currPVal = -10000000;
		float currQVal = 0; 
		for (int i=0; i< toFixQVal.length; i++){
			//scores = pVal,qVal,empFDR,tSumPlus,tSumMinus,cSumPlus,cSumMinus
			float[] scores = toFixQVal[i].getScores();
			float testPVal = scores[index];			
			if (currPVal != testPVal) {
				currPVal = testPVal;
				//look up qVal, if exact match use that qVal otherwise take smaller one
				int mark = Arrays.binarySearch(adjPVals, currPVal);
				//not exact therefore interpolate
				if (mark < 0) {
					mark = (-1 * mark) -2;
					if (mark < 0) mark = 0;
					int markPlusOne = mark +1;
					//currQVal = Num.interpolateY(adjPVals[mark], qVals[mark], adjPVals[markPlusOne], qVals[markPlusOne], currPVal);
					double trans = Num.interpolateY(Num.antiNeg10log10(adjPVals[mark]), Num.antiNeg10log10(qVals[mark]), Num.antiNeg10log10(adjPVals[markPlusOne]), Num.antiNeg10log10(qVals[markPlusOne]), Num.antiNeg10log10(currPVal));
					currQVal = Num.minus10log10Float(trans);
				}
				else currQVal = qVals[mark];
			}
			scores[index] = currQVal;
		}

	}

	/**Call after generatingSortedPAndQValArrays. Assumes scores are in -10Log10(value) space*/
	public void loadQValues(Point[] toFixQVal){
		//sort Point[] by value
		ComparatorPointDecendingScore comp = new ComparatorPointDecendingScore();
		Arrays.sort(toFixQVal, comp);
		//for every Point pValue, look up closest adjPVal and it's associated qVal
		float currPVal = -10000000;
		float currQVal = 0; 
		for (int i=0; i< toFixQVal.length; i++){
			float testPVal = toFixQVal[i].getScore();			
			if (currPVal != testPVal) {
				currPVal = testPVal;
				//look up qVal, if exact match use that qVal otherwise take smaller one
				int mark = Arrays.binarySearch(adjPVals, currPVal);
				//not exact therefore interpolate
				if (mark < 0) {
					mark = (-1 * mark) -2;
					if (mark < 0) mark = 0;
					int markPlusOne = mark +1;
					//currQVal = Num.interpolateY(adjPVals[mark], qVals[mark], adjPVals[markPlusOne], qVals[markPlusOne], currPVal);
					double trans = Num.interpolateY(Num.antiNeg10log10(adjPVals[mark]), Num.antiNeg10log10(qVals[mark]), Num.antiNeg10log10(adjPVals[markPlusOne]), Num.antiNeg10log10(qVals[markPlusOne]), Num.antiNeg10log10(currPVal));
					currQVal = Num.minus10log10Float(trans);
				}
				else currQVal = qVals[mark];
			}
			toFixQVal[i].setScore(currQVal);
		}
	}
	
	/**Call after adding all of the AdjPValQVals*/
	public void generateSortedPAndQValArrays(){
		//convert to array
		AdjPValQVal[] pq = new AdjPValQVal[alPValQVal.size()];
		alPValQVal.toArray(pq);
		//sort AdjPValQVal array
		Arrays.sort(pq);
		//eliminate dups
		ArrayList<AdjPValQVal> collapsedPQ = new ArrayList<AdjPValQVal>();
		float currQVal = pq[0].qVal;
		collapsedPQ.add(pq[0]);
		for (int i=1; i< pq.length; i++){			
			if (pq[i].qVal != currQVal) {
				collapsedPQ.add(pq[i]);
				currQVal = pq[i].qVal;
			}
		}
		//add last
		collapsedPQ.add(pq[pq.length-1]);
		
		pq = new AdjPValQVal[collapsedPQ.size()];
		collapsedPQ.toArray(pq);
		//split
		adjPVals= new float[pq.length];
		qVals = new float[pq.length];
		for (int i=0; i< pq.length; i++){
			adjPVals[i] = pq[i].adjPVal;
			qVals[i] = pq[i].qVal;
		}
	}

	/**Calculates p-values using a binomial probablility from built in R function. Set min10Log10Transform to true if you want
	 * to -10*Log10(p-val) transform the p-values.
	 * Note, don't set multiplyLastPValByRandomNumber = true unless 1) your max obs < 30 and 2) you want to generate a uniform distribution of binomial p-values using random error.*/
	private void makeBinomialPValParts(){
		//make array to feed to R
		int counter = 0;
		int[][] obs = new int[maxNumOb*maxNumOb][2];
		for (int x=0; x< maxNumOb; x++){
			for (int y=0; y< maxNumOb; y++){
				obs[counter++] = new int[]{x,y};
			}
		}
		//make random word
		String rndWrd = Passwords.createRandowWord(6);
		//write scores to file
		if (tempDirectory.exists() == false) tempDirectory.mkdir();
		File scores = new File(tempDirectory, rndWrd+"_Scores.txt");
		File rResults = new File (tempDirectory, rndWrd+"_RResults.txt");
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");
		try {
			PrintWriter out = new PrintWriter (new FileWriter (scores));
			//for each window
			for (int i=0; i< obs.length; i++) out.println(obs[i][0]+"\t"+obs[i][1]);
			out.close();
			//build R script
			StringBuilder script = new StringBuilder();
			script.append("e = "+expect+"\n");
			script.append("tc = read.table('"+scores.getCanonicalPath()+"')\n");
			script.append("numRow = nrow(tc); ");
			script.append("numObs = tc[,1] + tc[,2]\n");
			script.append("r = runif(numRow)\n");
			script.append("res = matrix (0,numRow,2)\n");
			//script.append("res = 1- (pbinom(tc[,1]-2, numObs, e) + (dbinom(tc[,1]-1, numObs, e) * r)); ");
			script.append("res[,1] = pbinom(tc[,1]-2, numObs, e)\n");
			script.append("res[,2] = dbinom(tc[,1]-1, numObs, e)\n");
			script.append("write.table(res, file='"+rResults.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE, sep = \"\t\")\n");

			//write script
			IO.writeString(script.toString(), scriptFile);
			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};			
			//execute
			IO.executeCommandLine(command);
			//read in results 
			double[] pbinom = new double[obs.length];
			double[] dbinom = new double[obs.length];
			if (rResults.exists() != false){
				String[] tokens;
				Pattern pat = Pattern.compile("\t");
				String line;
				BufferedReader in = new BufferedReader ( new FileReader(rResults));
				counter =0;
				while ((line=in.readLine()) != null){
					tokens = pat.split(line);
					pbinom[counter] = Double.parseDouble(tokens[0]);
					dbinom[counter++] = Double.parseDouble(tokens[1]);
				}
				if (counter != obs.length) throw new Exception("\nIncorrect length of R results file for fuzzy qval fix. See "+rResults);
			}
			else throw new Exception("\nR results file doesn't exist. Check tempFiles to debug.\n");

			//make and load matrix
			counter = 0;
			pbinomTC = new double[maxNumOb][maxNumOb];
			for (int x=0; x< maxNumOb; x++){
				for (int y=0; y< maxNumOb; y++){
					pbinomTC[x][y] = pbinom[counter++];
				}
			}
			counter = 0;
			dbinomTC = new double[maxNumOb][maxNumOb];
			for (int x=0; x< maxNumOb; x++){
				for (int y=0; y< maxNumOb; y++){
					dbinomTC[x][y] = dbinom[counter++];
				}
			}
			//cleanup
			scores.delete();
			rResults.delete();
			rOut.delete();
			scriptFile.delete();

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printExit("\nProblem with estimating binomial pvalues in R.\n");
		} 
	}	

	/*
	public static void main(String[] args){
		FuzzyBinomialPValueGenerator f = new FuzzyBinomialPValueGenerator( new File("/Users/nix/Desktop/TempFiles"), new File ("/usr/bin/R"), 0.3);
		double s = f.fetchPVal(12, 3);
		System.out.println("Straight "+s);
		for (int i=0; i< 5; i++){
			System.out.println(f.fetchFuzzyPVal(12, 3));
		}
	}*/
	private class AdjPValQVal implements Comparable{
		float adjPVal;
		float qVal;

		public AdjPValQVal(float adjPVal, float qVal){
			this.adjPVal = adjPVal;
			this.qVal = qVal;
		}

		public int compareTo(Object arg0) {
			AdjPValQVal other = (AdjPValQVal) arg0;
			if (other.adjPVal> adjPVal) return -1;
			if (other.adjPVal< adjPVal) return 1;
			return 0;
		}
	}

	public ArrayList<AdjPValQVal> getAlPValQVal() {
		return alPValQVal;
	}

	public float[] getAdjPVals() {
		return adjPVals;
	}
}
