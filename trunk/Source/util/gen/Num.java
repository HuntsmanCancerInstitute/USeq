package util.gen;
import java.util.*;
import java.util.regex.Pattern;
import java.text.*;
import java.io.*;

/**
 * Math and statistical static methods.
 */
public class Num {
	public static final double log10 = Math.log(10);
	public static final double log2 = Math.log(2);
	private static final double sRConstant = Math.log(Math.PI) / 2.0;

	/**Does a linear regression type interpolation, assumes 2's are > 1's.*/
	public static float interpolateY(float x1, float y1, float x2, float y2, float fixedX){
		float diffX = x2-x1;
		float diffFixed = fixedX- x1;
		float ratio = diffFixed/diffX;
		float diffY = y2-y1;
		float diff = diffY * ratio;
		return diff + y1;
	}
	
	/**Does a linear regression type interpolation, assumes 2's are > 1's.*/
	public static double interpolateY(double x1, double y1, double x2, double y2, double fixedX){
		double diffX = x2-x1;
		double diffFixed = fixedX- x1;
		double ratio = diffFixed/diffX;
		double diffY = y2-y1;
		double diff = diffY * ratio;
		return diff + y1;
	}
	
	/**Using interbase coordinates so length = stop - start.*/
	public static int calculateMiddleIntergenicCoordinates(int start, int end){
		if (start == end) return start;
		double length = end - start;
		double halfLength = length/2.0;
		return (int)Math.round(halfLength) + start;
	}


	/**@return number of trues in boolean[]*/
	public static int countNumberTrue (boolean[] b){
		int num =0;
		for (int i=0; i< b.length; i++) if (b[i]) num++;
		return num;
	}
	
	/**Takes t[reps][exonCounts] c[reps][exonCounts] and returns int[t/c][permutation][sumExonCounts] for permutated chi-square calculations*/
	public static int[][][] permutateReplicas(int[][] tReps, int[][] cReps){
		int numberTReplicas = tReps.length;
		int numberCReplicas = cReps.length;
		int numberExons = tReps[0].length;
		int numPermutations = numberTReplicas * numberCReplicas;
		int[][] tSum = new int[numPermutations][numberExons];
		int[][] cSum = new int[numPermutations][numberExons];
		int permutationIndex = 0;
		
		//for each t replica
		for (int t =0; t< numberTReplicas; t++){
			int[] tRepToSwap = tReps[t];
			//for each c replica
			for (int c=0; c< numberCReplicas; c++){
				int[] cRepToSwap = cReps[c];
				//make arrays to swap, don't want to modify original!
				int[][] tPerm = new int[numberTReplicas][];
				for (int i=0; i<numberTReplicas; i++) tPerm[i] = tReps[i];
				int[][] cPerm = new int[numberCReplicas][];
				for (int i=0; i<numberCReplicas; i++) cPerm[i] = cReps[i];
				//make swap
				tPerm[t] =cRepToSwap;
				cPerm[c] =tRepToSwap;
				//sum replica exon counts
				tSum[permutationIndex] = sumReplicaCounts (tPerm);
				cSum[permutationIndex++] = sumReplicaCounts (cPerm);
			}
		}
		return new int[][][]{tSum, cSum};
	}
	
	/**Takes int[replicas][exonCounts] and returns a sum of the replica counts at each exon*/
	public static int[] sumReplicaCounts (int[][] repsCounts){
		int[] summedCounts = new int[repsCounts[0].length];
		//for each replica
		for (int r=0; r<repsCounts.length; r++){
			//for each exon
			for (int e=0; e<summedCounts.length; e++){
				summedCounts[e] += repsCounts[r][e];
			}
		}
		return summedCounts;
	}

	/**Given a start bp (included) and stop bp (not included), returns start (included) and stop (not included) indexes.
	 * May return startIndex = endIndex, therefore nothing found. Be sure to sort the int[]!*/
	public static int[] findIndexes(int startBp, int stopBp, int[] sortedPositions){
		//find start index, included
		int startIndex = Arrays.binarySearch(sortedPositions, startBp);
		if (startIndex < 0) {
			startIndex = (startIndex*-1) -1;
		}
		else {
			//find first instance of startBp
			while (true){
				int minOne = startIndex - 1;
				if (minOne < 0) break;
				if (sortedPositions[minOne] != startBp) break;
				startIndex = minOne;
			}
		}
		//find stop index, not included
		int stopBPMinOne = stopBp-1;
		int endIndex = Arrays.binarySearch(sortedPositions, stopBPMinOne);		
		if (endIndex < 0) {
			endIndex = (endIndex*-1) -1;

		}
		else {
			//find last instance of endBp
			while (true){
				int addOne = endIndex +1;
				if (addOne >= sortedPositions.length) break;
				if (sortedPositions[addOne] != stopBPMinOne) break;
				endIndex = addOne;
			}
			//add one to stop index, it's not included
			endIndex++;
		}
		return new int[]{startIndex, endIndex};
	}
	
	/**Scales the array to 1:1000, converts all negative values to positive.*/
	public static void scale1To1000(float[] f){
		for (int i=0; i< f.length; i++) if (f[i]<0) f[i] = Math.abs(f[i]);
		float[] minMax = Num.findMinMaxFloatValues(f);
		for (int i=0; i< f.length; i++) f[i] -= minMax[0];
		float range = minMax[1]-minMax[0];
		float multi = 999.0f/range;
		for (int i=0; i< f.length; i++) f[i] = Math.round((f[i]) * multi)+1;
	}
	
	/**Scales the array to 100:1000*/
	public static int[] scale100To1000(double[] f){
		double[] minMax = Num.findMinMaxDoubleValues(f);
		double range = minMax[1]-minMax[0];
		double multi = 899.0/range;
		int[] scaled = new int[f.length];
		for (int i=0; i< f.length; i++) scaled[i] = (int)Math.round((f[i]-minMax[0]) * multi)+101;
		return scaled;
	}
	
	/**Provided a count matrix where the columns represent the counts from each gene in a given sample, countMatrix[gene][obs].
	 * Returns a scalar that a given genes counts should be multiplied by to normalize for differences in the total number of
	 * counts in a given sample.  normCount = scalar * observedCounts */
	public static double[] estimateSizeFactors(int[][] countMatrix){
		
		//get total obs on each column
		int numColumns = countMatrix[0].length;
		double[] columnTotals = new double[numColumns];
		//for each row
		for (int i=0; i<countMatrix.length; i++){
			//for each column
			for (int j=0; j< numColumns; j++){
				columnTotals[j] += countMatrix[i][j];
			}
		}
		
		//find smallest
		double smallest = columnTotals[0];
		for (int i=1; i< numColumns; i++) if (columnTotals[i] < smallest) smallest = columnTotals[i];
		
		//make ratios
		for (int i=0; i< numColumns; i++) columnTotals[i] = smallest/columnTotals[i];
	
		return columnTotals;
	}
	
	public static double[][] normalizeCountMatrix (int[][] countMatrix){
		//fetch scalars
		double[] scalars = estimateSizeFactors(countMatrix);
		
		int numColumns = countMatrix[0].length;
		double[][] normalizedData = new double[countMatrix.length][numColumns];
		
		//for each row
		for (int i=0; i<countMatrix.length; i++){
			//for each column
			for (int j=0; j< numColumns; j++){
				normalizedData[i][j] = scalars[j] * (double)countMatrix[i][j];
			}
		}
		return normalizedData;
	}
	
	/**Given two matrices of genes (rows) and their number of observations for each sample (column), normalizes based on total counts for each sample, and returns the log2(pseT+1/ pseC+1).
	 * Be sure to have at minimum 3 reps for each t and c. Assumes same number rows.
	 * Returns double[0][gene] = pseLog2Ratios, double[1][gene] = smallestPairwiseLog2Ratio*/
	public static double[][] calculatePseudoMedianLog2Ratios (int[][] treatmentCountMatrix, int[][] controlCountMatrix){
		//get numbers
		int numberGeneRows = treatmentCountMatrix.length;
		int numberTreatmentSamples = treatmentCountMatrix[0].length;
		int numberControlSamples = controlCountMatrix[0].length;
		int totalNumberSamples = numberTreatmentSamples + numberControlSamples;
		
		//merge samples	
		int[][] merged = new int[numberGeneRows][totalNumberSamples];
		for (int r=0; r<numberGeneRows; r++){
			//add treatments
			for (int t=0; t< numberTreatmentSamples; t++) merged[r][t] = treatmentCountMatrix[r][t];
			//add controls
			int index = 0;
			for (int m=numberTreatmentSamples; m< totalNumberSamples; m++) merged[r][m] = controlCountMatrix[r][index++];
		}
		
		//fetch normalized counts
		double[][] normMerged = normalizeCountMatrix(merged);
		//for each gene/ row, calculate the pseudoMedian for the t and c, then the ratio, then the log2Ratio
		double[] ratios = new double[numberGeneRows];
		double[] smallestAbsRatio = new double[numberGeneRows];
		double[] t = new double[numberTreatmentSamples];
		double[] c = new double[numberControlSamples];
		for (int r=0; r< numberGeneRows; r++){			
			//set t
			for (int x =0; x< numberTreatmentSamples; x++) t[x] = normMerged[r][x];
			//set c
			int index =0;
			for (int m=numberTreatmentSamples; m< totalNumberSamples; m++) c[index++] = normMerged[r][m];
			//calc
			double tPse = pseudoMedian (t);
			double cPse = pseudoMedian (c);
			double rto = (1.0 +tPse)/ (1.0 +cPse);
			ratios[r] = log2(rto);
			smallestAbsRatio[r] = findSmallestLog2RatioWithPseudoCounts(t, c);
		}

		return new double[][]{ratios, smallestAbsRatio};
	}
	
	/**Returns smallest abs all pair log2 ratio.*/
	public static double findSmallestLog2RatioWithPseudoCounts (double[] t, double[] c){
		double smallest = Double.MAX_VALUE;
		double smallestAbs = Double.MAX_VALUE;
		for (double a: t){
			for (double b: c){
				double ratio = log2((1.0 +a)/ (1.0 +b));
				double ratioAbs = Math.abs(ratio);
				if (ratioAbs < smallestAbs){
					smallestAbs = ratioAbs;
					smallest = ratio;
				}
			}
		}
		return smallest;
	}
	
	/**Calculates the permuted pval for all paired replica swaps with the chi-square test.
	 * @return -10Log10(numPerm>=RealChi / totalPerm) if numPerm == 0, numPerm set to 0.1
	 */
	public static double calculatePermutedChiSquarePValue (int[][] tReps, int[][] cReps){
		//calculate real chiSquare statistic
		int[] oriT = Num.sumReplicaCounts(tReps);
		int[] oriC = Num.sumReplicaCounts(cReps);
		int[][] combo = new int[][]{oriT,oriC};
		double oriChi = Num.chiSquareTestStatistic(Num.intArraysToLong(combo));
		
		int[][][] tc = Num.permutateReplicas(tReps, cReps);
		
		int[][] tSum = tc[0];
		int[][] cSum = tc[1];

		double numFail = 0;
		
		for (int i=0; i<tSum.length; i++){
			int[][] merge = new int[][]{tSum[i], cSum[i]};
			double permChi = Num.chiSquareTestStatistic(Num.intArraysToLong(merge));
			if (permChi >= oriChi) numFail++;
		}
		double totalPerm = (double)tSum.length;
		if (numFail == 0.0) numFail = 0.1;
		return Num.minus10log10(numFail/totalPerm);
	}
	
	 /**
     * @param counts array representation of 2-way table
     * @return chi-square test statistic
     * From apache Math package.
     */
    public static double chiSquareTestStatistic(long[][] counts) {

        int nRows = counts.length;
        int nCols = counts[0].length;

        // compute row, column and total sums
        double[] rowSum = new double[nRows];
        double[] colSum = new double[nCols];
        double total = 0.0d;
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                rowSum[row] += counts[row][col];
                colSum[col] += counts[row][col];
                total += counts[row][col];
            }
        }

        // compute expected counts and chi-square
        double sumSq = 0.0d;
        double expected = 0.0d;
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                expected = (rowSum[row] * colSum[col]) / total;
                //watch out for cases where expect is zero due to double zeros present in matrix (bad!) so just skip the addition to sumSq, a conservative approach
                if (expected == 0) continue;
                sumSq += ((counts[row][col] - expected) *
                        (counts[row][col] - expected)) / expected;
            }
        }
        return sumSq;
    }
    
    /**Casts ints to longs*/
    public static long[][] intArraysToLong (int[][] ints){
    	long[][] longs = new long[ints.length][];
    	for (int i=0; i< ints.length; i++){
    		longs[i] = new long[ints[i].length];
    		for (int j=0; j< ints[i].length; j++){
    			longs[i][j] = (long)ints[i][j];
    		}
    	}
    	return longs;
    }
    
    /**Rounds floats to ints*/
    public static int[][] floatArraysToInt (float[][] flt){
    	int[][] ints = new int[flt.length][];
    	for (int i=0; i< flt.length; i++){
    		ints[i] = new int[flt[i].length];
    		for (int j=0; j< flt[i].length; j++){
    			ints[i][j] = Math.round(flt[i][j]);
    		}
    	}
    	return ints;
    }
	
	/**@param treatmentsControls - int[window index][counts T1,T2,T..., counts C1,C2,C...]
	 * @param totalNumberObservationsInEachReplica - total number of observations in each treatment and control replica*/
	public static double[] quasipoisson (int[][] treatmentsControls, int numberTreatmentSamples, int[] totalNumberObservationsInEachReplica, File tempDirectory, File fullPathToR) throws Exception{
		//check number
		if (totalNumberObservationsInEachReplica.length != treatmentsControls[0].length) throw new Exception ("Number of totalNumberObservationsInEachReplica does not match number of treatment and control samples.");
		
		int[][] ids = new int[1][totalNumberObservationsInEachReplica.length];
		for (int i=0; i< numberTreatmentSamples; i++) ids[0][i] =1;
		
		//write arrays to file
		File tcFile = writeArrayToSingleLineFile(treatmentsControls, tempDirectory);
		int[][] obs = new int[1][];
		obs[0] = totalNumberObservationsInEachReplica;
		File observationsFile = writeArrayToSingleLineFile(obs, tempDirectory);
		File idsFile = writeArrayToSingleLineFile(ids, tempDirectory);
		
		File results = new File (tempDirectory, "tempFile_Results_"+Passwords.createRandowWord(6));
		
		//write script
		String script =
			"libraryPresent = library(aod, logical.return = T)\n"+
			"numRows = "+treatmentsControls.length+ "\n"+
			"numColumns = "+treatmentsControls[0].length+ "\n"+
			"y = matrix(scan('"+ tcFile+ "'), numRows, numColumns, byrow=T)\n"+
			"n = scan('"+observationsFile +"')\n"+
			"x = scan('"+idsFile +"')\n"+
			"res = 1:numRows\n"+
			"for (i in 1:numRows) {\n"+
			"	z = quasipois(y[i,] ~ x + offset(log(n)),data = data.frame(y[i,], x),tol=0.1)\n"+
			"	res[i] =  z@fm$null.deviance - z@fm$deviance \n"+
			"}\n"+
			"write.table(res, file='"+results+"',row.names = FALSE, col.names = FALSE, sep = \"\t\")\n";
		File rScriptFile = new File(tempDirectory, "tempFile_RScript_"+Passwords.createRandowWord(6));
		IO.writeString(script, rScriptFile);
		
		//make command
		File rOutFile = new File (tempDirectory, "tempFile_ROut_"+Passwords.createRandowWord(6));
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				rScriptFile.getCanonicalPath(),
				rOutFile.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		//load results
		double[] values = Num.loadDoubles(results);
		//check
		if (values == null || values.length != treatmentsControls.length) throw new Exception ("Number of results from R does not match the number of trials.");
		//clean up
		tcFile.deleteOnExit();
		observationsFile.deleteOnExit();
		idsFile.deleteOnExit();
		results.deleteOnExit();
		rScriptFile.deleteOnExit();
		rOutFile.deleteOnExit();
		return values;
	}
	
	/**Runs a modified chi square test in R looking for differences between two distributions of categories using a test for independence.
	 * Assumes that for any given rowNumber, the number of categories in each sample are the same. 
	 * R, given the screwed up "language" it is, cannot import ragged arrays, thus must add '-1' trailing place holders so the number of columns is the same, what the @#%$%@#$!!!!!
	 * @param sample int[rowNumber][counts of observations for each category]*/
	public static double[] chiSquareIndependenceTest (int[][] sampleA, int[][] sampleB, File tempDirectory, File fullPathToR, boolean useYatesCorrection) {
		try {
			//check sizes
			if (sampleA.length != sampleB.length) return null;
			int numberRows = sampleA.length;

			//make matrix
			int[][] data = new int[numberRows][];
			for (int i=0; i< numberRows; i++){
				data[i] = new int[sampleA[i].length *2];
				System.arraycopy(sampleA[i], 0, data[i], 0, sampleA[i].length);
				System.arraycopy(sampleB[i], 0, data[i], sampleA[i].length, sampleB[i].length);
			}

			//write matrix to file
			File matrixFile = writeArrayToFile(data, tempDirectory);

			//make holder for data
			File results = new File (tempDirectory, "tempFile_Results_"+Passwords.createRandowWord(6));

			String yates = "	res[i] = chisq.test(m,correct=FALSE)$p.value\n";
			if (useYatesCorrection) yates = "	res[i] = chisq.test(m,correct=TRUE)$p.value\n";
			String script = 
				"numberRows = "+numberRows+ "\n"+
				"bigM = read.table('"+matrixFile+"')\n"+
				"res = 1:numberRows\n"+
				"for (i in 1:numberRows){\n"+
				"	y = bigM[i,]\n"+
				"	y = y[y!=-1]\n"+
				"	m = matrix(y, nrow = 2, byrow = TRUE)\n"+
				yates+
				"	if (is.na(res[i])) res[i] = 1\n"+
				"}\n"+
				"res = -10 * log10(res)\n"+
				"write.table(res, file='"+results+"',row.names = FALSE, col.names = FALSE, sep = \"\t\")\n";
			

			File rScriptFile = new File(tempDirectory, "tempFile_RScript_"+Passwords.createRandowWord(6));
			IO.writeString(script, rScriptFile);

			//make command
			File rOutFile = new File (tempDirectory, "tempFile_ROut_"+Passwords.createRandowWord(6));
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					rScriptFile.getCanonicalPath(),
					rOutFile.getCanonicalPath()};			
			//execute
			IO.executeCommandLine(command);
			//load results
			double[] values = Num.loadDoubles(results, 3234.0);
			//check
			if (values == null || values.length != data.length) throw new Exception ("Number of results from R does not match the number of trials.");
			//clean up
			matrixFile.deleteOnExit();
			results.deleteOnExit();
			rOutFile.deleteOnExit();
			results.deleteOnExit();
			rScriptFile.deleteOnExit();
			rOutFile.deleteOnExit();
			
			return values;

		} catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}

	
	/**Runs a modified fisher exact test in R looking for differences between two distributions of GATC using a test for independence. 
	 * @param gatcSample int[rowNumber][4; counts of observations for each base: GATC]*/
	public static double[] fisherTest (int[][] gatcSampleA, int[][] gatcSampleB, File tempDirectory, File fullPathToR) {
		try {
		
		//make warped matrix
		int[][] data = new int[gatcSampleA.length][8];
		for (int i=0; i< data.length; i++){
			data[i] = new int[8];
			int counter = 0;
			for (int x=0; x< 4; x++){
				data[i][counter++] = gatcSampleA[i][x];
				data[i][counter++] = gatcSampleB[i][x];
			}
		}
		
		//write matrix to file
		File matrixFile = writeArrayToSingleLineFile(data, tempDirectory);
		
		//make holder for data
		File results = new File (tempDirectory, "tempFile_Results_"+Passwords.createRandowWord(6));
		
		//write script
		String script =
			"numberRows = "+data.length+ "\n"+
			"bigM = matrix ( scan (\""+ matrixFile +"\"),numberRows,8, byrow=T)\n"+
			"res = 1:numberRows\n"+
			"for (i in 1:numberRows){\n"+
			"	m = matrix (bigM[i,], nrow=2)\n"+
			"   if (max(m) ==0) res[i] = -1\n"+
			"	else res[i] = fisher.test(m)$p.value\n"+
			"}\n"+
			"write.table(res, file='"+results+"',row.names = FALSE, col.names = FALSE, sep = \"\t\")\n";
			
		File rScriptFile = new File(tempDirectory, "tempFile_RScript_"+Passwords.createRandowWord(6));
		IO.writeString(script, rScriptFile);
		
		//make command
		File rOutFile = new File (tempDirectory, "tempFile_ROut_"+Passwords.createRandowWord(6));
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				rScriptFile.getCanonicalPath(),
				rOutFile.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		//load results
		double[] values = Num.loadDoubles(results);
		//check
		if (values == null || values.length != data.length) throw new Exception ("Number of results from R does not match the number of trials.");
		//clean up
		matrixFile.deleteOnExit();
		results.deleteOnExit();
		rOutFile.deleteOnExit();
		results.deleteOnExit();
		rScriptFile.deleteOnExit();
		rOutFile.deleteOnExit();
		return values;
		
		} catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}

	
	
	public static File writeArrayToSingleLineFile (int[][] array, File tempDirectory) throws IOException {
		//make random word
		String rndWrd = "tempFile_Array_"+Passwords.createRandowWord(6);
		File file = new File (tempDirectory, rndWrd);
		//write to file
		PrintWriter out = new PrintWriter(new FileWriter(file));
		for (int i=0; i<array.length; i++){
			for (int j=0; j< array[i].length; j++){
				out.print(array[i][j]);
				out.print(" ");
			}
		}
		out.close();
		return file;
	}
	
	public static File writeArrayToFile (int[][] array, File tempDirectory) throws IOException {
		//make random word
		String rndWrd = "tempFile_Array_"+Passwords.createRandowWord(6);
		File file = new File (tempDirectory, rndWrd);
		//write to file
		PrintWriter out = new PrintWriter(new FileWriter(file));
		//for each row
		for (int i=0; i<array.length; i++){
			String line = intArrayToString(array[i], "\t");
			out.println(line);
		}
		out.close();
		return file;
	}

	/**Calculates the Poisson probablility for observed < 21.
	 * For observed > 20 call poissonProbabilityApproximation()
	 * @param rate - lambda, the average number of occurances
	 * @param observed - the actual number of occurances observed
	 * @return poissonProbability of observing the number of occurances by chance alone
	 * @author nix */
	public static double poissonProbablility (double rate, int observed){
		double a = Math.pow(rate, observed);
		double b = Math.pow(Math.E, -1*rate);
		double c = Num.factorial(observed);
		return (a*b)/c;
	}

	/**Calculates a natural log approximation of the Poisson probablility.
	 * @param mean - lambda, the average number of occurances
	 * @param observed - the actual number of occurances observed
	 * @return ln(poissonProbability) - the natural log of the poisson probability.*/
	public static double poissonProbabilityApproximation (double mean, int observed){
		double k = observed;
		double a = k * Math.log(mean);
		double b = -1.0 * mean * Math.log(Math.E);
		return a + b - Num.factorialApproximation(k);
	}
	/**Convert natural log to log 10.*/
	public static double ln2Log10 (double ln){
		return ln/ Math.log(10);
	}

	/**Convert natural log to log 10.*/
	public static double ln2Minus10Log10 (double ln){
		return -10* (ln/ Math.log(10));
	}
	/**Converter*/
	public static float[] relativeDifferenceToLog2Ratio (float[] relDiff){
		float[] log2Ratios = new float[relDiff.length];
		for (int i=0; i< relDiff.length; i++){
			double r = relDiff[i];
			Double log2Ratio = new Double(Math.log((r+2.0)/(2.0-r))/ log2);
			log2Ratios[i] = log2Ratio.floatValue();
		}
		return log2Ratios;
	}

	/**Converter*/
	public static double relativeDifferenceToLog2Ratio (double relDiff){
		double ratio = (relDiff+2)/(2-relDiff);
		return Math.log(ratio)/ log2;
	}


	/**Converts milliseconds to days.*/
	public static double millisecToDays (long ms){
		double current = (double)ms;
		current = current/1000;
		current = current/60;
		current = current/60;
		current = current/24;
		return current;
	}


	/**Writes a binary int[][2] array.
	 * @return true if sucessful, false if something bad happened.*/
	public static boolean writeBinaryInt2Array(int[][] c, File file){
		try {
			int num = c.length;
			FileOutputStream fos = new FileOutputStream(file);
			DataOutputStream dos = new DataOutputStream( new BufferedOutputStream (fos));
			//write length
			dos.writeInt(num);
			//write value
			for (int i=0; i<num; i++) { 
				dos.writeInt(c[i][0]);
				dos.writeInt(c[i][1]);
			}
			dos.close();
			fos.close();
			return true;
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return false; 
		}
	}

	/**Reads a binary int[][2] array file.
	 * @return null if something bad happened.*/
	public static int[][] readBinaryInt2Array(File file){
		try {
			FileInputStream fis = new FileInputStream(file);
			DataInputStream dis = new DataInputStream( new BufferedInputStream(fis ));
			//read lenth
			int num = dis.readInt();
			//make array
			int[][] c = new int[num][2];
			for (int i=0; i< num; i++) {
				c[i][0] = dis.readInt();
				c[i][1] = dis.readInt();
			}
			dis.close();
			fis.close();
			return c;
		}
		catch (Exception ioe){
			ioe.printStackTrace();
			return null;   
		} 
	}

	/**Inverts the float[].*/
	public static void invertArray(float[] x){
		int lenthMinusOne = x.length -1;
		int stop = x.length/2;
		for (int i=0; i<stop; i++){
			int switchIndex = lenthMinusOne - i;
			float first = x[i];
			float last = x[switchIndex];
			x[i] = last;
			x[switchIndex] = first;

		}
	}

	/**Converts a double[] to float[] */
	public static float[] convertToFloat(double[] d){
		float[] f = new float[d.length];
		for (int i=0; i< d.length; i++) f[i] = (float)d[i];
		return f;
	}

	/**Converts a double[][] to float[][] */
	public static float[][] convertToFloat(double[][] d){
		int sizeMain = d.length;
		float[][] f = new float[sizeMain][];
		for (int i=0; i< sizeMain; i++) {
			int sizeMinor = d[i].length;
			f[i] = new float[sizeMinor];
			for (int j=0; j< sizeMinor; j++){
				f[i][j] = (float)d[i][j];
			}
		}
		return f;
	}
	
	/**Converts a double[][] to Float[][] */
	public static Float[][] convertToFloatObjectArray(double[][] d){
		int sizeMain = d.length;
		Float[][] f = new Float[sizeMain][];
		for (int i=0; i< sizeMain; i++) {
			int sizeMinor = d[i].length;
			f[i] = new Float[sizeMinor];
			for (int j=0; j< sizeMinor; j++){
				f[i][j] = new Float(d[i][j]);
			}
		}
		return f;
	}

	/**Calculates N! modified from http://www.unix.org.ua/orelly/java-ent/jnut/ch01_03.htm
	 * Good for n < 21.*/
	public static double factorial(int n) {  
		if (n < 0) return 0;  
		if (n > 22 ) return -1;
		double fact = 1;                    
		while(n > 1) {                         
			fact = fact * n;                    
			n = n - 1;                          
		}                                      
		return fact;                           
	}  

	/**Srivivasa Ramanuja ln(n!) factorial estimation, better than Stirling's approximation.
	 * Good for larger values of n.
	 * @return ln(n!)*/
	public static double factorialApproximation(double n){
		if (n < 2) return 0;
		double a = (n * Math.log(n)) - n;
		double b = (Math.log(n * (1.0+ (4 * n * (1.0 + (2.0 * n)))))) / 6.0;
		return a + b + sRConstant;
	}


	/**Converts an ArrayList of int[]s to int[][].*/
	public static int[][] arrayList2IntArrayArray (ArrayList al){
		int length = al.size();
		int[][] intArray = new int[length][];
		for (int x= 0; x < length; x++){
			intArray[x] = (int[])al.get(x);
		}
		return intArray;
	}

	/**Calculates p-values using a Poisson probablility from built in R function.
	 * Set bonCorr to 1 if you don't want to Bonferroni correct the pvalues.
	 * You should check that the pvalues don't exceed 1 if you do correct.
	 * @return pvalues where the index is the observed number of hits. */
	public static double[] poissonPValues(int maxHit, double rate, File fullPathToR, File tempDir, double bonCorr){
		try {
			//make random word
			String rndWrd = Passwords.createRandowWord(6);
			//make temp files to hold data
			File scriptFile = new File (tempDir, rndWrd+"_TempFileRScript");
			File rOutFile = new File (tempDir, rndWrd+"_TempFileRScript.Rout");
			File resultsFile = new File (tempDir, rndWrd+"_TempFileRResults");
			//make script
			String script = "l="+rate+
			";x=0:"+maxHit+
			";y = ppois(x,l,lower.tail=FALSE) * "+bonCorr+";write.table(y,\""+resultsFile+
			"\",row.names=FALSE,col.names=FALSE)\n";
			//write script
			IO.writeString(script, scriptFile);
			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOutFile.getCanonicalPath()};			
			//execute
			IO.executeCommandLine(command);
			//load results
			double[] pvalues = Num.loadDoubles(resultsFile);
			//add a 1 to the zero index
			double[] pvaluesPlus = new double[pvalues.length+1];
			pvaluesPlus[0] = 1;
			System.arraycopy(pvalues, 0, pvaluesPlus, 1, pvalues.length);
			//delete files
			scriptFile.deleteOnExit();
			rOutFile.deleteOnExit();
			resultsFile.deleteOnExit();
			if (pvalues == null || pvalues.length != (maxHit+1)) return null;
			return pvaluesPlus;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	/**Run's Storey's qvalue FDR conversion of well behaved pvalues (uniformally distributed). 
	 * Set min10Log10Transformed to true if the pvalues are so transformed. Will return transformed qValueFDRs.
	 */
	public static float[] qValueFDR(File tempDirectory, float[] pValues, File fullPathToR, boolean min10Log10Transformed, boolean verbose){
		//make random word
		String rndWrd = Passwords.createRandowWord(8);
		//write scores to file
		File scores = new File(tempDirectory, rndWrd+"_Scores.txt");
		File rResults = new File (tempDirectory, rndWrd+"_RResults.txt");
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");
		float[] pvalues = null;
		try {
			PrintWriter out = new PrintWriter (new FileWriter (scores));
			//for each window
			for (int i=0; i< pValues.length; i++) out.println(pValues[i]);
			out.close();
			//build R script
			StringBuilder script = new StringBuilder("s = read.table('"+scores.getCanonicalPath()+"'); library(qvalue); ");
			if (min10Log10Transformed) script.append("s[,1] = 10^(-1*s[,1]/10) ; "	);
			script.append("qobj = qvalue (s[,1]); s= c(1); len = length (qobj$qvalues); start = 1; \n");
			script.append("while (TRUE){ \n");
			script.append("   stop = start + 1000000 \n");
			script.append("   if (stop > len) stop = len \n");
			script.append("   toPrint = qobj$qvalues[start:stop] \n");
			if (min10Log10Transformed) script.append("   toPrint = -10 * log10(toPrint) \n");
			script.append("   write.table(toPrint, file='"+rResults.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE, append = TRUE, sep = \"\t\") \n");
			script.append("   start = stop + 1 \n");
			script.append("   if (start > len) break \n");
			script.append("}");
			
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
			if (rResults.exists() != false){
				pvalues = Num.loadFloats(rResults);
				if (pvalues.length != pValues.length) throw new Exception("\nIncorrect length of R results file. See "+rResults);
			}
			else throw new Exception("\nR results file doesn't exist.\n");

			//cleanup
			scores.deleteOnExit();
			rResults.deleteOnExit();
			rOut.deleteOnExit();
			scriptFile.deleteOnExit();

		} catch (Exception e) {
			if (verbose) e.printStackTrace();
		} 
		return pvalues;
	}
	
	/**Run's Storey's qvalue FDR conversion of well behaved pvalues (uniformally distributed). 
	 * Set min10Log10Transformed to true if the pvalues are so transformed. Will return a file containing transformed qValueFDRs.
	 */
	public static File qValueFDR(File tempDirectory, File pValues, File fullPathToR, boolean min10Log10Transformed, boolean verbose){
		//make random word
		String rndWrd = Passwords.createRandowWord(8);
		//write scores to file
		File scores = pValues;
		File rResults = new File (tempDirectory, rndWrd+"_RResults.txt");
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");
		try {
			//build R script
			StringBuilder script = new StringBuilder("s = read.table('"+scores.getCanonicalPath()+"'); library(qvalue); ");
			if (min10Log10Transformed) script.append("s[,1] = 10^(-1*s[,1]/10) ; "	);
			script.append("qobj = qvalue (s[,1]); s= c(1); len = length (qobj$qvalues); start = 1; \n");
			script.append("while (TRUE){ \n");
			script.append("   stop = start + 1000000 \n");
			script.append("   if (stop > len) stop = len \n");
			script.append("   toPrint = qobj$qvalues[start:stop] \n");
			if (min10Log10Transformed) script.append("   toPrint = -10 * log10(toPrint) \n");
			script.append("   write.table(toPrint, file='"+rResults.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE, append = TRUE, sep = \"\t\") \n");
			script.append("   start = stop + 1 \n");
			script.append("   if (start > len) break \n");
			script.append("}");
			
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
			if (rResults.exists() == false)  throw new Exception("\nR results file doesn't exist.\n");

			//cleanup
			rOut.deleteOnExit();
			scriptFile.deleteOnExit();

		} catch (Exception e) {
			if (verbose) e.printStackTrace();
		} 
		return rResults;
	}
	
	/**Assumes sorted pvalues, descending. Alters the input array.*/
	public static void benjaminiHochbergCorrect(double[] sortedPValDecending){
		double num = sortedPValDecending.length;
		double prior = 1;
		for (int i=1; i< sortedPValDecending.length; i++){
			sortedPValDecending[i] = sortedPValDecending[i] * num / (num-i);
			if(sortedPValDecending[i] < prior) prior = sortedPValDecending[i]; 
			else sortedPValDecending[i] = prior;
		}
	}
	
	/**Assumes pvalues are -10Log10(pval) transformed and sorted in ascending order. Alters the input array.*/
	public static void benjaminiHochbergCorrect(float[] sortedMin10Log10PVals){
		double num = sortedMin10Log10PVals.length;
		double prior = 1;
		for (int i=1; i< sortedMin10Log10PVals.length; i++){
			double val = Num.antiNeg10log10(sortedMin10Log10PVals[i]);
			val = val * num / (num-i);
			if(val < prior) prior = val; 
			else val = prior;
			sortedMin10Log10PVals[i] = Num.minus10log10Float(val);
		}
	}
	
	/**Applies the Benjamini & Hochberg FDR correction to the incoming array. Use this method to correct threholded pvalues by setting the offset.
	 * @param sortedMin10Log10PVals -10Log10(pval) transformed and sorted in ascending order
	 * @param offset Number of pvalues with poor significance not included in the array
	 * @return nada, modifies the original array*/
	public static void benjaminiHochbergCorrect(float[] sortedMin10Log10PVals, long offset){
		double num = sortedMin10Log10PVals.length + offset;
		double offsetDouble = (double)offset;
		double prior = 1;
		for (int i=1; i< sortedMin10Log10PVals.length; i++){
			double val = Num.antiNeg10log10(sortedMin10Log10PVals[i]);
			val = val * num / (num-i-offsetDouble);
			if(val < prior) prior = val; 
			else val = prior;
			sortedMin10Log10PVals[i] = Num.minus10log10Float(val);
		}
		//check first value, if larger than second, replace with second
		if (sortedMin10Log10PVals[0] > sortedMin10Log10PVals[1]) sortedMin10Log10PVals[0] = sortedMin10Log10PVals[1];
	}


	/**Standard G-test, log-likelihood test, likelihood ratio test.
	 * Treatment and control values should be integers.
	 * The expectedProportion relates to the treatment.
	 * @return the G statistic.
	 * From http://udel.edu/~mcdonald/statgtestgof.html*/
	public static double gTest(double treatment, double control, double expectedProportionTreatment){
		double total = treatment + control;
		double expectedTreatment = expectedProportionTreatment * total;
		double expectedControl = total - expectedTreatment;
		double gTreatment = treatment * Math.log(treatment/expectedTreatment);
		double gControl = control * Math.log(control/expectedControl);
		return 2 * (gTreatment+gControl);
	}

	/**Calculates p-values using a binomial probablility from built in R function. Set min10Log10Transform to true if you want
	 * to -10*Log10(p-val) transform the p-values.
	 * Note, don't set multiplyLastPValByRandomNumber = true unless 1) your max obs < 30 and 2) you want to generate a uniform distribution of binomial p-values using random error.*/
	public static double[] binomialPValues(double expect, File tempDirectory, int[][] treatObsControlObs, File fullPathToR, boolean min10Log10Transform, boolean multiplyLastPValByRandomNumber, boolean leaveMaxValuesAsDouble_Min_Value){
		//make random word
		String rndWrd = Passwords.createRandowWord(6);
		//write scores to file
		File scores = new File(tempDirectory, rndWrd+"_Scores.txt");
		File rResults = new File (tempDirectory, rndWrd+"_RResults.txt");
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");
		double[] pvalues = null;
		try {
			PrintWriter out = new PrintWriter (new FileWriter (scores));
			//for each window
			for (int i=0; i< treatObsControlObs.length; i++) out.println(treatObsControlObs[i][0]+"\t"+treatObsControlObs[i][1]);
			out.close();
			//build R script
			StringBuilder script = new StringBuilder();
			if (multiplyLastPValByRandomNumber){		
				script.append("e = "+expect+"; ");
				script.append("tc = read.table('"+scores.getCanonicalPath()+"'); ");
				script.append("numRow = nrow(tc); ");
				script.append("numObs = tc[,1] + tc[,2] ;");
				script.append("r = runif(numRow); ");
				script.append("res = 1- (pbinom(tc[,1]-2, numObs, e) + (dbinom(tc[,1]-1, numObs, e) * r)); ");
				if (min10Log10Transform) script.append("res = -10 * log10(res); ");
				script.append("write.table(res, file='"+rResults.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE, sep = \"\t\") ");
			}
			else{
				script.append("s = read.table('"+scores.getCanonicalPath()+"'); ");
				if (min10Log10Transform) {
					//script.append("p = -10 * log10(pbinom(s[,1]-1, s[,1]+s[,2], "+expect+ ", lower.tail=FALSE)); ");
					script.append("logTen = log(10); ");
					script.append("p = -10 * ((pbinom(s[,1]-1, s[,1]+s[,2], "+expect+", lower.tail=FALSE, log=TRUE)/logTen)); ");
				}
				else script.append("p = pbinom(s[,1]-1, s[,1]+s[,2], "+expect+ ", lower.tail=FALSE); ");
				script.append("write.table(p, file='"+rResults.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE, sep = \"\t\") ");
			}
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
			//read in results watching for Inf values from R
			if (rResults.exists() != false){
				pvalues = Num.loadDoubles(rResults, treatObsControlObs, leaveMaxValuesAsDouble_Min_Value);
				if (pvalues.length != treatObsControlObs.length) throw new Exception("\nIncorrect length of R results file. See "+rResults);
			}
			else throw new Exception("\nR results file doesn't exist.\n");

			//cleanup
			scores.deleteOnExit();
			rResults.deleteOnExit();
			rOut.deleteOnExit();
			scriptFile.deleteOnExit();

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printExit("\nProblem with estimating binomial pvalues in R.\n");
		} 
		return pvalues;
	}
	
	/**Calculates p-values using a binomial probability from built in R function. Returns up and down pvals as -10log10(pval) for each.*/
	public static double[][] binomialPValues(float[][] treatmentControlExpect, File tempDirectory, File fullPathToR){
		//make random word
		String rndWrd = Passwords.createRandowWord(6);
		//write scores to file
		File scores = new File(tempDirectory, rndWrd+"_Scores.txt");
		File rResults = new File (tempDirectory, rndWrd+"_RResults.txt");
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");
		double[][] pvalues = null;
		try {
			PrintWriter out = new PrintWriter (new FileWriter (scores));
			//for each window
			for (int i=0; i< treatmentControlExpect.length; i++) out.println(treatmentControlExpect[i][0]+"\t"+treatmentControlExpect[i][1]+"\t"+treatmentControlExpect[i][2]);
			out.close();
			//build R script
			StringBuilder script = new StringBuilder();
			script.append("s = read.table('"+scores.getCanonicalPath()+"'); ");
			script.append("logTen = log(10); ");
			script.append("p = -10 * ((pbinom(s[,1]-1, s[,1]+s[,2], s[,3], lower.tail=FALSE, log=TRUE)/logTen)); ");
			script.append("q = -10 * ((pbinom(s[,2]-1, s[,1]+s[,2], 1-s[,3], lower.tail=FALSE, log=TRUE)/logTen)); ");
			script.append("z = matrix(nrow=nrow(s), ncol=2); z[,1] = p; z[,2] = q; ");
			script.append("write.table(z, file='"+rResults.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE, sep = \"\t\") ");
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
			//read in results watching for Inf values from R
			if (rResults.exists() != false){
				pvalues = Num.loadDoubleMatrix(rResults,2);
				if (pvalues.length != treatmentControlExpect.length) throw new Exception("\nIncorrect length of R results file. See "+rResults);
				Num.zeroNegativeValues(pvalues);
			}
			else throw new Exception("\nR results file doesn't exist.\n");
			//cleanup
			scores.deleteOnExit();
			rResults.deleteOnExit();
			rOut.deleteOnExit();
			scriptFile.deleteOnExit();

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printExit("\nProblem with estimating binomial pvalues in R.\n");
		} 
		return pvalues;
	}

	/**Returns a maxtrix of binomial pvalues for fast look up.*/
	public static double[][] binomialPValMatrix (int maxNumOb, double expectation, File tempDirectory, File R, boolean logTransPVals){
		int counter = 0;
		//make array to feed to R
		int[][] obs = new int[maxNumOb*maxNumOb][2];
		for (int x=0; x< maxNumOb; x++){
			for (int y=0; y< maxNumOb; y++){
				obs[counter++] = new int[]{x,y};
			}
		}
		//calc pvals
		double[] pvals = Num.binomialPValues(expectation, tempDirectory, obs, R, logTransPVals, false, false);
		//make and load matrix
		counter = 0;
		double[][] pvalMatrix = new double[maxNumOb][maxNumOb];
		for (int x=0; x< maxNumOb; x++){
			for (int y=0; y< maxNumOb; y++){
				pvalMatrix[x][y] = pvals[counter++];
			}
		}
		pvals = null;
		obs = null;
		return pvalMatrix;
	}

	/**Converts an array of double using -10 * (Math.log(double)/log10).
	 * Any doubles <=0 are assigned the minimal double from the array >0 then transformed.*/
	public static void convertToNeg10Log10(double[] d){
		double min = Double.MIN_VALUE;
		//find min
		for (int i=0; i< d.length; i++){
			if (d[i] < min && d[i] > 0.0) min = d[i];
		}
		min = -10 * (Math.log(min)/log10);
		//convert
		for (int i=0; i< d.length; i++){
			if (d[i] > 0.0) d[i] = -10 * (Math.log(d[i])/log10);
			else d[i] = min;
		}
	}

	/**Given a HashSet of Floats returns an array of float.*/
	public static float[] hashSetToFloat(HashSet<Float> hash){
		float[] f = new float[hash.size()];
		Iterator<Float> it = hash.iterator();
		int counter = 0;
		while (it.hasNext()) f[counter++] = it.next().floatValue();
		return f;
	}

	/**Given a HashSet of Integers returns an array of int.*/
	public static int[] hashSetToInt(HashSet<Integer> hash){
		int[] f = new int[hash.size()];
		Iterator<Integer> it = hash.iterator();
		int counter = 0;
		while (it.hasNext()) f[counter++] = it.next().intValue();
		return f;
	}

	/**Converts to unlogged value.*/
	public static double antiLog(double loggedValue, double base){
		return Math.pow(base, loggedValue);
	}

	/**Converts a -10Log10(val) back to the val.*/
	public static double antiNeg10log10(float fredScore){
		double s = -1.0 * fredScore/10.0;
		return Math.pow(10, s);
	}

	/**Converts an array float[] to unlogged values.*/
	public static float[] antiLog(float[] loggedValues, int base){
		float[] values = new float[loggedValues.length];
		for (int i=0; i< loggedValues.length; i++){
			values[i] = new Double(Math.pow(base, loggedValues[i])).floatValue();
		}
		return values;
	}

	/**Scans a float array for zeros, if found returns true, other wise returns false.*/
	public static boolean findZeros(float[] f){
		for (int i=0; i< f.length; i++){
			if (f[i] == 0) return true;
		}
		return false;
	}

	/**Rotates a square matrix 90 degrees CounterClockWise?*/
	public static float[][] rotateCounterClockwise(float [][] matrix){
		int numRows = matrix.length;
		float[][] rotated = new float[numRows][numRows];
		//for each row
		for (int i=0; i< numRows; i++){
			int rX = numRows - i -1;
			//for each column
			for (int j=0; j< numRows; j++){
				rotated[j][rX] = matrix[i][j];
			}
		}
		return rotated;
	}

	/**Rotates a square matrix 90 degrees Clock Wise?*/
	public static float[][] rotateClockwise(float [][] matrix){
		int numRows = matrix.length;
		float[][] rotated = new float[numRows][numRows];
		//for each row
		for (int i=0; i< numRows; i++){
			//for each column
			for (int j=0; j< numRows; j++){
				int rY = numRows - j -1;
				rotated[rY][i] = matrix[i][j];
			}
		}
		return rotated;
	}

	/**Returns the first index containing the maximum value. */
	public static int findMaxIntIndex(int[] ints) {
		int len = ints.length;
		int max = ints[0];
		int maxIndex = 0;
		for (int i = 1; i < len; i++) {
			if (ints[i]>max) {
				max=ints[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}

	/**Log2s every number in the array replacing the original values.*/
	public static void log2(float[][] f){
		double log2 = Math.log(2);
		int num = f.length;
		for (int i=0; i<num; i++){
			int num2 = f[i].length;
			for (int j=0; j<num2; j++){
				f[i][j] = new Double(Math.log(f[i][j])/log2).floatValue();
			}
		}
	}


	/**Log2s every number in the array returning the values, keep original intact.*/
	public static float[][] log2Return(float[][] f){
		double log2 = Math.log(2);
		int num = f.length;
		float[][] logged = new float[num][];
		for (int i=0; i<num; i++){
			int num2 = f[i].length;
			logged[i] = new float[num2];
			for (int j=0; j<num2; j++){
				logged[i][j] = new Double(Math.log(f[i][j])/log2).floatValue();
			}
		}
		return logged;
	}

	/**Log2s every number in the array replacing the original values.*/
	public static void log2(float[] f){
		double log2 = Math.log(2);
		int num = f.length;
		for (int i=0; i<num; i++){
			f[i] = new Double(Math.log(f[i])/log2).floatValue();
		}
	}

	/**Log2s every number in the array replacing the original values.*/
	public static void log2(double[] f){
		int num = f.length;
		for (int i=0; i<num; i++){
			f[i] = Math.log(f[i])/log2;
		}
	}

	/**Log2s every number in the array keeping original intact.*/
	public static float[] log2Return(float[] f){
		double log2 = Math.log(2);
		int num = f.length;
		float[] logged = new float[num];
		for (int i=0; i<num; i++){
			logged[i] = new Double(Math.log(f[i])/log2).floatValue();
		}
		return logged;
	}	

	/**Log2s every number in the array keeping original intact.*/
	public static double[] log2ReturnDouble(float[] f){
		double log2 = Math.log(2);
		int num = f.length;
		double[] logged = new double[num];
		for (int i=0; i<num; i++){
			logged[i] = Math.log(f[i])/log2;
		}
		return logged;
	}

	/**Loads columns of double into a HashMap, the first row is parsed as Strings and used as the key.
	 * The value is a double[].
	 * Can convert log10 values in the files to log2 (for Agilent data).
	 * Skips blank lines.*/
	public static LinkedHashMap loadDoubles(File dataFile1, boolean convertLog10ToLog2){
		LinkedHashMap idValues = new LinkedHashMap();
		String line = null;
		try {
			BufferedReader in = new BufferedReader ( new FileReader (dataFile1));
			//read first line
			line = in.readLine();
			String[] ids = line.split("\\t");
			ArrayList[] values = new ArrayList[ids.length];
			//instantiate ArrayLists
			for (int i=0; i< values.length; i++) values[i] = new ArrayList(20000);
			//read in values converting to log2
			while ((line = in.readLine())!=null){
				line = line.trim();
				if (line.length() == 0) continue;
				String[] cells = line.split("\\t");
				for (int i=0; i< ids.length; i++){
					double log = Double.parseDouble(cells[i]);
					if (convertLog10ToLog2){
						log = Num.antiLog(log, 10);
						log = Num.log2(log);
					}
					values[i].add(new Double(log));
				}
			}
			//add to hash
			for (int i=0; i<ids.length; i++){
				idValues.put(ids[i], Num.arrayListOfDoubleToArray(values[i]));
			}
			in.close();
		} catch (Exception e){
			System.out.println(dataFile1+ " Bad line "+line);
			e.printStackTrace();
		}
		return idValues;
	}

	/**Parses doubles from a String[] of doubles. Returns null if a parsing error is encountered.*/
	public static double[] parseDoubles(String[] doubles){
		int num = doubles.length;
		double[] d = new double[num];
		try {
			for (int i=0; i<num; i++){
				d[i] = Double.parseDouble(doubles[i]);
			}
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return d;
	}

	/**Parses ints from a String[] of ints. Returns null if a parsing error is encountered.*/
	public static int[] parseInts(String[] ints){
		int num = ints.length;
		int[] d = new int[num];
		try {
			for (int i=0; i<num; i++){
				d[i] = Integer.parseInt(ints[i]);
			}
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return d;
	}

	/**Converts an array of values to z scores given their mean and standard deviation.*/
	public static double[] convertToZScores(double[] values, double mean, double stndDev){
		int num = values.length;
		double[] z = new double[num];
		for (int i=0; i< num; i++){
			z[i] = (values[i]-mean)/stndDev;
		}
		return z;
	}
	/**Converts an array of values to z scores.*/
	public static double[] convertToZScores(double[] values){
		double mean = Num.mean(values);
		double stndDev = Num.standardDeviation(values,mean);
		return convertToZScores(values, mean, stndDev);
	}

	/**Converts an array of values to z scores using median as appose to mean.*/
	public static double[] convertToMedianZScores(double[] values){
		double median = Num.median(values);
		double stndDev = Num.standardDeviation(values,median);
		return convertToZScores(values, median, stndDev);
	}

	/**Returns a new array with the newValue appended on the stop.*/
	public static double[] appendDouble (double[] d, double newValue){
		double[] added = new double[d.length+1];
		System.arraycopy(d,0,added,0,d.length);
		added[d.length] = newValue;
		return added;
	}
	/**Returns a new array with the newValue appended on the end.*/
	public static float[] appendFloat (float[] d, float newValue){
		float[] added = new float[d.length+1];
		System.arraycopy(d,0,added,0,d.length);
		added[d.length] = newValue;
		return added;
	}

	/**Returns a new array with the newValue appended on the beginning.*/
	public static double[] prependDouble (double[] d, double newValue){
		double[] added = new double[d.length+1];
		System.arraycopy(d,0,added,1,d.length);
		added[0] = newValue;
		return added;
	}

	/**Calculates a one-step tukey biweight estimator using 
	 * 5 as the tuning constant and adding 0.0001 to the denominator
	 * to avoid dividing by zero, Affy's params.
	 * If all zeros are submitted then returns zero.*/
	public static double tukeyBiWeight(double[] x){
		int num = x.length;

		//find median
		Arrays.sort(x);
		double m = Num.median(x);		

		//find median diff
		double[] diffs = new double[num];
		for (int i=0; i< num; i++){
			diffs[i] = Math.abs(x[i]-m);
		}
		Arrays.sort(diffs);
		double mad = Num.median(diffs);	

		//calculate ts and bs
		double[] t = new double[num];
		double[] b = new double[num];
		double top = 0;
		double bottom = 0;
		for (int i=0; i<num; i++){
			//calc t
			t[i] = Math.abs((x[i] - m) / ((5 * mad)) + 0.0001);		
			//calc b
			if (t[i] <1){
				double w = t[i] * t[i];
				w = (1-w);
				b[i] = w * w;
			}
			else b[i] = 0;
			//increment sums
			top += (b[i] * x[i]);
			bottom += b[i];
		}
		if (bottom == 0) return 0;
		return top/bottom;
	}

	/**Geometric mean. Cannot have any zeros or negative numbers in the double[].*/
	public static double geometricMean(double[] x) {
		int num = x.length;
		//take logs
		double[] logs = new double[num];
		for (int i=0; i< num; i++){
			logs[i] = Math.log(x[i]);
		}
		//mean
		double mean = Num.mean(logs);
		//final
		return Math.pow(Math.E, mean);
	}

	/**Returns the value of a given percentile from a SORTED array.
	 * Percentile is from 0-1, ie 0.95 and is according to Lentner, 1982.*/
	public static float percentile(float[] sortedFloats, double percentile){
		//calculate index
		double index = ((percentile * (double)sortedFloats.length)) - 0.5;
		int trunkIndex = (int)index;
		//is it a whole number?
		double rnd = index%1;
		if (rnd<0.00001 || rnd> 0.99999){
			return sortedFloats[trunkIndex];
		}
		//otherwise average trunk and trunk +1
		else {
			return (sortedFloats[trunkIndex] + sortedFloats[trunkIndex+1])/2;
		}
	}


	/**Runs through float[][] setting any values that exceed the outlierThreshold
	 * to zero. Assumes float[i]s are equal length.*/
	public static void maskOutliers (float[][] f, float outlierThreshold){
		int numA = f.length;
		int numB = f[0].length;
		for (int i=0; i< numA; i++){
			for (int j=0; j< numB; j++){
				if (f[i][j]>outlierThreshold) f[i][j] = 0;
			}
		}
	}

	/**Returns the number of elements/ values in the array.*/
	public static int countValues(int[][] i){
		int num = i.length;
		int count = 0;
		for (int h=0; h<num; h++){
			count += i[h].length;
		}
		return count;
	}

	/**Returns the number of objects in the arrays.*/
	public static int countObjects(Object[][] i){
		int num = i.length;
		int count = 0;
		for (int h=0; h<num; h++){
			count += i[h].length;
		}
		return count;
	}

	/**Counts the number of values that are > or < threshold based on boolean.
	 * Returns - number of outliers > or < threshold
	 * Param - exceedsThreshold, boolean flag when true counts the number that
	 * exceed the threshold, when false the number less than the threshold.*/
	public static int countOutliers(float[] f, double threshold, boolean exceedsThreshold){
		int num = f.length;
		int count = 0;
		if (exceedsThreshold) {
			for (int i=0; i<num; i++) if (f[i]>threshold) count++;
		}
		else {
			for(int i=0; i<num; i++) if (f[i]<threshold) count++;
		}
		return count;
	}

	/**Loads a zeroed matrix with values that are > or < threshold based on boolean.
	 * Returns a loaded zeroed matrix.
	 * Param exceedsThreshold, boolean flag when true returns the values that
	 * exceed the threshold,  less than the threshold.*/
	public static float[][] identifyOutliers(float[][] f, int[][] controls, float threshold, boolean exceedsThreshold){
		//make a zeroed float[][] to load with outliers
		float[][] outliers = Num.zeroedFloatArray(f.length, f.length);
		//walk through control coordinates, fetch original matrix values and look for outliers, if found add to new
		int num = controls.length;
		if (exceedsThreshold) {
			for (int i=0; i<num; i++){
				int x = controls[i][0];
				int y = controls[i][1];
				if (f[x][y]> threshold)outliers[x][y] = f[x][y];
			}
		}
		else {
			for (int i=0; i<num; i++){
				int x = controls[i][0];
				int y = controls[i][1];
				if (f[x][y]< threshold)outliers[x][y] = f[x][y];
			}
		}
		return outliers;
	}

	/**Given a square table and a list of coordinates, will return the associated values.*/
	public static float[] fetchMatrixValues (float[][] matrix, int[][] xy){
		int numValues = xy.length;
		float[] floats = new float[numValues];
		for (int i=0; i<numValues; i++){
			floats[i]= matrix [xy[i][0]]  [xy[i][1]];
		}
		return floats;
	}

	/**Given a square table and a list of coordinates, save the associated values.*/
	public static void saveMatrixValues (float[][] matrix, int[][] xy, float threshold, File f){
		int numValues = xy.length;
		float[][] floats = new float[matrix.length][matrix[0].length];
		for (int i=0; i<numValues; i++){
			if (matrix [xy[i][0]]  [xy[i][1]]> threshold) floats [xy[i][0]]  [xy[i][1]]= matrix [xy[i][0]]  [xy[i][1]];
		}
		IO.saveObject(f, floats);
	}

	/**Given a square table and a list of coordinates, will return the associated values 
	 * in a new zeroed square table.*/
	public static float[][] loadMatrixValues (float[][] matrix, int[][] xy){
		int numValues = xy.length;
		float[][] floats = zeroedFloatArray(matrix.length,matrix.length);
		for (int i=0; i<numValues; i++){
			floats[xy[i][0]][xy[i][1]] = matrix[xy[i][0]][xy[i][1]];
		}
		return floats;
	}

	/**Calculates the Q1,Q2,and Q3 (25th, 50th, and  75th percentile).
	 * Ignores the middle value in calculating Q1 and Q3 if the array length is odd.
	 * Don't forget to sort the array!*/
	public static double[] quartiles(float[] sortedFloatArray){
		double p25th;
		double p50th;
		double p75th;
		int middle = sortedFloatArray.length/2;  // subscript of middle element
		int middleLeft = middle/2;
		int middleRight = middle+middleLeft+1;
		//Odd number of elements -- return the middle one.
		if (sortedFloatArray.length%2 == 1) {
			p50th = sortedFloatArray[middle];
			p25th = ((double)sortedFloatArray[middleLeft-1] + (double)sortedFloatArray[middleLeft]) / 2.0;
			p75th = ((double)sortedFloatArray[middleRight-1] + (double)sortedFloatArray[middleRight]) / 2.0;
		}
		// Even number -- return average of middle two
		else {
			//average
			p50th = ((double)sortedFloatArray[middle-1] + (double)sortedFloatArray[middle]) / 2.0;
			p25th = ((double)sortedFloatArray[middleLeft-1] + (double)sortedFloatArray[middleLeft]) / 2.0;
			p75th = ((double)sortedFloatArray[middleRight-1] + (double)sortedFloatArray[middleRight]) / 2.0;
		}
		return new double[] {p25th, p50th, p75th};
	}


	/**Makes a float[x][y] and fills it with zeros.*/
	public static float[][] zeroedFloatArray (int x, int y){	
		//make filtering array, fill with zeros
		float[][] intensities = new float[x][y];
		for (int i= x-1; i>=0; i--){
			for (int j= y-1; j>=0; j--){
				intensities[i][j]=0;
			}
		}
		return intensities;
	}
	/**Makes a float[x][y] and fills it with the value.*/
	public static float[][] loadFloatArray (int x, int y, float value){	
		//make filtering array, fill with zeros
		float[][] intensities = new float[x][y];
		for (int i= x-1; i>=0; i--){
			for (int j= y-1; j>=0; j--){
				intensities[i][j]=value;
			}
		}
		return intensities;
	}

	/**Returns an array of float[numToSample] populated with randomly picked values from x.*/
	public static float[] randomSample(float[] x, int numToSample){
		float[] z = new float[numToSample];
		Random generator = new Random();
		for (int i=0; i<z.length; i++){
			z[i] = x[ generator.nextInt(x.length)];
		}
		return z;
	}
	
	public static void randomize (float[] array, long seed){
	    Random rng = new Random(seed);       
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables
	        float tmp = array[k];
	        array[k] = array[n - 1];
	        array[n - 1] = tmp;
	    }
	}
	
	/**Randomized the paired arrays, keeping the indexed scores matched.*/
	public static void randomizePairedValues (float[] array, float[] array2, long seed){
	    Random rng = new Random(seed);       
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables, keeping array and array2 together
	        int nMinOne = n-1;
	        float tmp = array[k];
	        array[k] = array[nMinOne];
	        array[nMinOne] = tmp;
	        float tmp2 = array2[k];
	        array2[k] = array2[nMinOne];
	        array2[nMinOne] = tmp2;
	    }
	}
	
	public static void randomize (int[] array, long seed){
	    Random rng = new Random(seed);       
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables
	        int tmp = array[k];
	        array[k] = array[n - 1];
	        array[n - 1] = tmp;
	    }
	}


	/**Randomizes the first array by permutating indexes and swapping. Fast.
	 * Shuffles f[shuffle][].*/
	public static void randomizeFirstArray(float[][] f, long seed){
		int len = f.length;
		float[] current;
		float[] random;
		int index;
		Random generator = new Random(seed);
		for (int i=0; i<len; i++){
			index = generator.nextInt(len);
			current = f[i];
			random = f[index];
			f[i] = random;
			f[index] = current;
		}
	}

	/**Randomizes the intensities, not between replicas.
	 * @param float[replicas][intensities]*/
	public static void randomize (float[][] f){
		long seed = System.currentTimeMillis();
		for (int i=0; i<f.length; i++){
			seed +=i;
			randomize(f[i], seed);
		}
	}

	/**Given two float[replica][intensities] for treat and cont, for each oligo makes and returns all pairwise ratios 
	 * between treat and control replicas.
	 * @param float[replica][intensities] for treatment and control
	 * @return float[oligo index][all t/c ratios]*/
	public static float[][] layeredRatiosSeperate(float[][] treatment, float[][]control){
		int numOligos = treatment[0].length;
		int numTreat = treatment.length;
		int numCont = control.length;
		int totalRatiosPerOligo = numTreat*numCont;
		float[][] layers = new float[numOligos][];
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			float[] oligoRatios = new float[totalRatiosPerOligo];
			int counter = 0;
			//for each treatment
			for (int y=0; y< numTreat; y++){
				//for each control
				for (int z=0; z< numCont; z++){
					oligoRatios[counter++] = treatment[y][x]/control[z][x];
				}
			}
			layers[x] = oligoRatios;
		}
		return layers;
	}

	/**Given two float[replica][intensities] for treat and cont, for each oligo makes and returns all pairwise relative differences 
	 * between treat and control replicas.
	 * @param float[replica][intensities] for treatment and control
	 * @return float[oligo index][all relative differences]*/
	public static float[][] layeredRelativeDifferencesSeperate(float[][] treatment, float[][]control){
		int numOligos = treatment[0].length;
		int numTreat = treatment.length;
		int numCont = control.length;
		int totalScoresPerOligo = numTreat*numCont;
		float[][] layers = new float[numOligos][];
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			float[] oligoRelDiffs = new float[totalScoresPerOligo];
			int counter = 0;
			//for each treatment
			for (int y=0; y< numTreat; y++){
				//for each control
				for (int z=0; z< numCont; z++){
					oligoRelDiffs[counter++] = 2.0f *((treatment[y][x] - control[z][x]) / (treatment[y][x] + control[z][x]));
				}
			}
			layers[x] = oligoRelDiffs;
		}
		return layers;
	}


	/**Given two float[replica][values] for treat and cont, for each replica makes all pairwise relative differences 
	 * between treat and control. Returns as a big pool.*/
	public static double[] layeredRelativeDifferences(float[][] oligoTreat, float[][]oligoCont){
		int numOligos = oligoTreat[0].length;
		int numTreat = oligoTreat.length;
		int numCont = oligoCont.length;
		int totalRatios = numTreat*numCont*numOligos;
		double[] ratios = new double[totalRatios];
		int counter = 0;
		//for each oligo, calculate all relative differences between treatments and controls
		for (int x=0; x< numOligos; x++){
			for (int y=0; y< numTreat; y++){
				for (int z=0; z< numCont; z++){
					ratios[counter++] = Num.relativeDifference(oligoTreat[y][x], oligoCont[z][x]);
				}
			}
		}
		return ratios;
	}

	/**Given two float[replica][intensities] for treat and cont, for each oligo makes and returns all pairwise log2 ratios 
	 * between treat and control replicas.
	 * @param float[replica][intensities] for treatment and control
	 * @return float[oligo index][all log2(t/c) ratios]*/
	public static float[][] layeredLogRatiosSeperate(float[][] treatment, float[][]control){
		int numOligos = treatment[0].length;
		int numTreat = treatment.length;
		int numCont = control.length;
		int totalRatiosPerOligo = numTreat*numCont;
		double log2 = Math.log(2);
		float[][] layers = new float[numOligos][];
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			float[] oligoRatios = new float[totalRatiosPerOligo];
			int counter = 0;
			//for each treatment
			for (int y=0; y< numTreat; y++){
				//for each control
				for (int z=0; z< numCont; z++){
					oligoRatios[counter++] = new Double (Math.log(treatment[y][x]/control[z][x]) / log2).floatValue();
				}
			}
			layers[x] = oligoRatios;
		}
		return layers;
	}

	/**Given two float[replica][values] for treat and cont, for each replica makes all pairwise ratios 
	 * between treat and control.*/
	public static double[] layeredRatios(float[][] oligoTreat, float[][]oligoCont){
		int numOligos = oligoTreat[0].length;
		int numTreat = oligoTreat.length;
		int numCont = oligoCont.length;
		int totalRatios = numTreat*numCont*numOligos;
		double[] ratios = new double[totalRatios];
		int counter = 0;
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			for (int y=0; y< numTreat; y++){
				for (int z=0; z< numCont; z++){
					ratios[counter++] = oligoTreat[y][x]/oligoCont[z][x];
				}
			}
		}
		return ratios;
	}

	/**Given two float[replica][values] for treat and cont, for each replica makes matched pairwise ratios 
	 * between treat and control. Must have equal number of replicas for both T and C.*/
	public static double[] pairedRatios(float[][] treatments, float[][]controls){
		int numOligos = treatments[0].length;
		int numReplicas = treatments.length;
		int totalRatios = numReplicas*numOligos;
		double[] ratios = new double[totalRatios];
		int counter = 0;
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			//for each replica, make pair
			for (int y=0; y< numReplicas; y++){
				ratios[counter++] = treatments[y][x]/ controls[y][x];
			}
		}
		return ratios;
	}


	/**Given two float[replica][values] for treat and cont, for each replica makes matched pairwise relative differences 
	 * between treat and control. Must have equal number of replicas for both T and C.*/
	public static double[] pairedRelativeDifferences(float[][] treatments, float[][]controls){
		int numOligos = treatments[0].length;
		int numReplicas = treatments.length;
		int totalRatios = numReplicas*numOligos;
		double[] ratios = new double[totalRatios];
		int counter = 0;
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			//for each replica, make pair
			for (int y=0; y< numReplicas; y++){
				ratios[counter++] = Num.relativeDifference(treatments[y][x], controls[y][x]);
			}
		}
		return ratios;
	}


	/**Given two float[replica][values] for treat and cont, for each replica makes all pairwise 
	 * log2 ratios between treat and control.
	 * Zero ratios are assigned a log2 of 0.0000000001*/
	public static double[] layeredLogRatios(float[][] oligoTreat, float[][]oligoCont){
		int numOligos = oligoTreat[0].length;
		int numTreat = oligoTreat.length;
		int numCont = oligoCont.length;
		int totalRatios = numTreat*numCont*numOligos;
		double[] ratios = new double[totalRatios];
		int counter = 0;
		double log2 = Math.log(2);
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			for (int y=0; y< numTreat; y++){
				for (int z=0; z< numCont; z++){
					double r = oligoTreat[y][x]/oligoCont[z][x];
					if (r > 0 ) {
						ratios[counter] = Math.log(r) /log2;
					}
					else {
						ratios[counter] = 0.0000000001;
					}
					counter++;
				}
			}
		}
		return ratios;
	}

	/**Given two float[replica][values] for treat and cont, for each replica makes matched pairwise log ratios 
	 * between treat and control. Must have equal number of replicas for both T and C.*/
	public static double[] pairedLogRatios(float[][] treatments, float[][]controls){
		int numOligos = treatments[0].length;
		int numReplicas = treatments.length;
		int totalRatios = numReplicas*numOligos;
		double[] ratios = new double[totalRatios];
		int counter = 0;
		//for each oligo, calculate all ratios between treatments and controls
		for (int x=0; x< numOligos; x++){
			//for each replica, make pair
			for (int y=0; y< numReplicas; y++){
				double r = treatments[y][x]/ controls[y][x];
				if (r > 0 ) {
					ratios[counter] = Math.log(r) /log2;
				}
				else {
					ratios[counter] = 0.0000000001;
				}
				counter++;

			}
		}
		return ratios;
	}

	/**Takes unlogged MAT normalized t-values, sorts, logs and takes a variant of a number stabilized trimmed mean ratio.
	 * Returns 1 if trimmed mean will zero the number of values.*/
	public static double matScore(float[] t, float[] c){
		//proc T
		Arrays.sort(t);
		double[] logT = Num.log2ReturnDouble(t);
		int trimNumberT = (int)Math.round(((double)t.length) * 0.1);
		int numAfterTrimming = t.length-(2*trimNumberT);
		if (numAfterTrimming <= 0) return 1;
		double trMeanT = trimmedMean(logT, trimNumberT);

		//proc C
		Arrays.sort(c);
		double[] logC = Num.log2ReturnDouble(c);
		int trimNumberC = (int)Math.round(((double)c.length) * 0.1);
		numAfterTrimming = c.length-(2*trimNumberC);
		if (numAfterTrimming <= 0) return 1;

		double trMeanC = trimmedMean(logC, trimNumberC);
		return (Math.sqrt(trimNumberT) * trMeanT)  -  (Math.sqrt(trimNumberC) * trMeanC);
	}


	/**Returns a Hodges-Lehmann estimator, the median of all pairwise means.*/
	public static double pseudoMedian (double[] d){
		int len = d.length;
		int numPairs = (len*(len-1))/2;
		double[] means = new double[numPairs];
		int counter = 0;
		for (int i=0; i< len; i++){
			double one = d[i];
			for (int j=i+1; j< len; j++){
				double two = d[j];
				means[counter++] = (one+two)/2;
			}
		}
		Arrays.sort(means);
		return median(means);
	}

	/**Returns the integer square root from zero to num.*/
	public static float[] squareRoots(int num){
		float[] sqrs = new float[num];
		for (int i=1; i< num; i++){
			sqrs[i] = new Double(Math.sqrt(i)).floatValue();
		}
		return sqrs;
	}

	/**Returns all pairwise averages.*/
	public static float[] pairwiseMeans (float[] d){
		int len = d.length;
		int numPairs = (len*(len-1))/2;
		float[] means = new float[numPairs];
		int counter = 0;
		for (int i=0; i< len; i++){
			for (int j=i+1; j< len; j++){
				means[counter++] = (d[i]+d[j])/2.0f;
			}
		}
		return means;
	}

	/**Returns all pairwise averages.*/
	public static double[] pairwiseMeansDouble (float[] d){
		int len = d.length;
		int numPairs = (len*(len-1))/2;
		double[] means = new double[numPairs];
		int counter = 0;
		for (int i=0; i< len; i++){
			for (int j=i+1; j< len; j++){
				means[counter++] = (d[i]+d[j])/2.0f;
			}
		}
		return means;
	}

	/**Returns all pairwise averages between a and b but not self pairwise.*/
	public static double[] partialPairwiseMeans (float[] a, float[] b){
		int numPairs = a.length * b.length;
		double[] means = new double[numPairs];
		int counter = 0;
		for (int i=0; i< a.length; i++){
			for (int j=0; j< b.length; j++){
				means[counter++] = (a[i]+b[j])/2.0;
			}
		}
		return means;
	}

	/**Concatinates the float[]s, variable sizes OK.*/
	public static float[] concatinate(float[][] f){
		//calc size
		int size = 0;
		for (int i=0; i< f.length; i++){
			size+= f[i].length;
		}
		//make new
		float[] concat = new float[size];
		int index =0;
		for (int i=0; i< f.length; i++){
			System.arraycopy(f[i], 0, concat, index, f[i].length);
			index += f[i].length;
		}
		return concat;
	}

	/**Joins two float[]s using System.arraycopy().*/
	public static float[] concatinate(float[] one, float[] two){
		float[] merge = new float[one.length+two.length];
		System.arraycopy(one,0,merge,0,one.length);
		System.arraycopy(two,0,merge,one.length,two.length);
		return merge;
	}


	/** Builds a sub array from part of a float[replica][intensities].
	 * @param float[replicas][intensities]*/
	public static float[][] subArray(float[][] intensities, int startIndex, int stopIndex){
		int numIntensities = 1+stopIndex-startIndex;
		float[][] sub = new float[intensities.length][numIntensities];
		//for each intensity/ oligo
		int counter=0;
		for (int j=startIndex; j<=stopIndex; j++){
			//for each replica
			for (int i=0; i<intensities.length; i++) { 
				sub[i][counter] = intensities[i][j];
			}
			counter++;
		}
		return sub;
	}


	/**Calculates the pseudoMedian on each replica of intensities. If only one intensity,
	 * sets that as the pseudo median.
	 * @param repInt float[replica index][actual intensity values]*/
	public static double[] pseudoMedian(float[][] repInt){
		double[] sum = new double[repInt.length];
		//for each replica
		for (int i=0; i< repInt.length; i++){
			if (repInt[i].length == 1) sum[i] = repInt[i][0];
			else sum[i] = Num.pseudoMedian(repInt[i]);
		}
		return sum;
	}

	/**Returns a Hodges-Lehmann estimator, the median of all pairwise means.
	 * If only one pairwise mean, return the mean instead of the median on the means.*/
	public static double pseudoMedian (float[] d){
		int len = d.length;
		//only two values?
		if (len == 2) return mean(d);
		if (len == 1) return d[0];
		int numPairs = (len*(len-1))/2;
		double [] means = new double[numPairs];
		int counter = 0;
		for (int i=0; i< len; i++){
			double one = d[i];
			for (int j=i+1; j< len; j++){
				double two = d[j];
				means[counter++] = (one+two)/2;
			}
		}
		Arrays.sort(means);
		return median(means);
	}

	/**Removes zero values and their assoicated values from all arrays.
	 * Pass two float[numChips][oligoValues] representing treatment and control chips.
	 * The num oligoValues should be the same, the numChips can vary between t and c.
	 * Returns float[2][][] representing t[][] and c[][].*/
	public static float[][][] removeZeroValues (float[][] t, float[][] c){
		int windowSize = t[0].length;
		int numTreatments = t.length;
		int numControls = c.length;
		int numZeros = 0;
		float[][][] mod = new float[2][][];

		//first count the number of zeros found

		//for each oligo look for a zero
		for (int x=0; x< windowSize; x++){
			boolean zeroFound = false;
			//count number of zeros in treatment
			for (int y=0; y<numTreatments; y++){
				if (t[y][x] == 0){					
					numZeros++;
					zeroFound = true;
					break;
				}
			}
			//count number of zeros in control
			if (zeroFound == false){
				for (int y=0; y<numControls; y++){
					if (c[y][x] == 0){
						numZeros++;
						break;
					}
				}
			}	
		}

		//if any zeros present then make new arrays
		if (numZeros != 0){
			int newWindowSize = windowSize - numZeros;
			float[][] tMod = new float[numTreatments][newWindowSize];
			float[][] cMod = new float[numControls][newWindowSize];
			int counter = 0;
			for (int x=0; x< windowSize; x++){
				boolean zeroFound = false;
				//treatment
				for (int y=0; y<numTreatments; y++){
					if (t[y][x] == 0){
						zeroFound = true;
						break;
					}
					else{
						tMod[y][counter] = t[y][x];
					}
				}
				//count number of zeros in control
				if (zeroFound == false){
					for (int y=0; y<numControls; y++){
						if (c[y][x] == 0){
							zeroFound = true;
							break;
						}
						else{
							cMod[y][counter] = c[y][x];
						}
					}
				}	
				//only advance counter if no zeros found
				if (zeroFound == false) counter++;
				//watch out for zeros at the stop of the array
				if (counter == newWindowSize) break;
			}

			//set final
			mod[0] = tMod;
			mod[1] = cMod;
		}
		else {
			mod[0] = t;
			mod[1] = c;
		}
		return mod;
	}

	/**Sets values exceeding the cutoff to the default.*/
	public static float[] trimOutliers(float[] f, float cutoff, float replacement){
		int num = f.length;
		for (int i=0; i< num; i++){
			if (f[i] > cutoff) {
				f[i] = replacement;
			}
		}
		return f;
	}

	/**For each f[replica][]Sets values exceeding the cutoff to the default.*/
	public static float[][] trimOutliers(float[][] f, float cutoff, float replacement){
		int num = f.length;
		for (int i=0; i< num; i++){
			f[i] = trimOutliers(f[i], cutoff, replacement);
		}
		return f;
	}



	/**Median normalize unsorted array to a given number.*/
	public static float[] medianNormalize (float[] f, float targetMedian){
		int num = f.length;
		//find median
		float[] sorted = new float[num];
		System.arraycopy(f,0,sorted,0,num);
		Arrays.sort(sorted);
		double median;
		median = median(sorted);
		//make scalar
		double scalar = targetMedian/median;
		//correct new
		System.arraycopy(f,0,sorted,0,num);
		for (int i=0; i<num; i++){
			Double scaled = new Double((double)sorted[i] * scalar);
			sorted[i] = scaled.floatValue();
		}
		return sorted;
	}

	/**Median normalize unsorted array to a given number.*/
	public static double[] medianNormalize (double[] f, double targetMedian){
		int num = f.length;
		//find median
		double[] sorted = new double[num];
		System.arraycopy(f,0,sorted,0,num);
		Arrays.sort(sorted);
		double median;
		median = median(sorted);
		//make scalar
		double scalar = targetMedian/median;
		//correct new
		System.arraycopy(f,0,sorted,0,num);
		for (int i=0; i<num; i++){
			sorted[i] = sorted[i] * scalar;
		}
		return sorted;
	}

	/**Median normalizes an unsorted float[][] array to a given target.*/
	public static float[][] medianNormalize (float[][] f, double targetMedian){
		//find median
		float[] sorted = collapseFloatArray(f);
		Arrays.sort(sorted);
		double median;
		median = median(sorted);
		//make scalar
		double scalar = targetMedian/median;
		//correct new
		int num = f.length;
		for (int i=0; i<num; i++){
			int num2= f[i].length;
			for (int j=0; j<num2; j++){
				Double scaled = new Double(((double)f[i][j]) * scalar);
				f[i][j] = scaled.floatValue();
			}
		}
		return f;
	}

	/**Calculates the median value of a sorted float[] ignoring leading zero values!*/
	public static double medianIgnoreZeros(float[] m) {
		//find start of actual values not zeros
		int start = 0;
		for (int i=0; i< m.length; i++){
			if (m[i]!= 0){
				start = i;
				break;
			}
		}
		int length = m.length-start;
		int middle = length/2;  // subscript of middle element
		//Odd number of elements -- return the middle one.
		if (m.length%2 == 1) return m[start+ middle];
		// Even number -- return average of middle two
		return ((double)m[start + middle-1] + (double)m[start+ middle]) / 2.0;
	}	

	/**Convert to int[]*/
	public static int[] convertToInt(double[] d){
		int num = d.length;
		int[] ints = new int[num];
		for (int i=0; i<num; i++){
			ints[i] = (int)Math.round(d[i]);
		}
		return ints;
	}
	/**Convert to double[]*/
	public static double[] convertToDouble(int[] d){
		int num = d.length;
		double[] ints = new double[num];
		for (int i=0; i<num; i++){
			ints[i] = d[i];
		}
		return ints;
	}	
	/**Convert to double[]*/
	public static double[] convertToDouble(float[] d){
		int num = d.length;
		double[] ints = new double[num];
		for (int i=0; i<num; i++){
			ints[i] = d[i];
		}
		return ints;
	}
	/**Convert to int[]*/
	public static int[] convertToInt(float[] d){
		int num = d.length;
		int[] ints = new int[num];
		for (int i=0; i<num; i++){
			ints[i] = Math.round(d[i]);
		}
		return ints;
	}
	/**Convert to float[]*/
	public static float[] convertToFloat(int[] d){
		int num = d.length;
		float[] floats = new float[num];
		for (int i=0; i<num; i++){
			floats[i] = d[i];
		}
		return floats;
	}
	/**Convert to long[]*/
	public static long[] convertToLong(float[] d){
		int num = d.length;
		long[] ints = new long[num];
		for (int i=0; i<num; i++){
			ints[i] = Math.round(d[i]);
		}
		return ints;
	}

	/**Calculates relative differences. Uses 2 * ( (t-c)/(t+c) )  */
	public static double[] relativeDifferences(float[] t, float[] c){
		int num = t.length;
		double[] relDiffs = new double[num];
		for (int i=num-1; i>=0; i--){
			relDiffs[i] = 2 * ( (double)(t[i]-c[i]) / (double)(t[i]+c[i])   );
		}
		return relDiffs;
	}

	/**Calculates a relative difference. Uses 2 * ( (t-c)/(t+c) )  */
	public static double relativeDifference(double t, double c){
		return 2 * (  (t - c) / (t + c)   );
	}

	/**Calculates log base 2 ratios.*/
	public static double[] logRatios(double[] t, double[] c){
		int num = t.length;
		double[] ratios = new double[num];
		for (int i=num-1; i>=0; i--){
			ratios[i] = log2(t[i]/c[i]);
		}
		return ratios;
	}
	/**Calculates log base 2 ratios.*/
	public static double[] logRatios(float[] t, float[] c){
		int num = t.length;
		double[] ratios = new double[num];
		for (int i=num-1; i>=0; i--){
			ratios[i] = log2(t[i]/c[i]);
		}
		return ratios;
	}

	/**Rounds a double[] based on number of desired decimals, ie 1 = x.x, 2 = x.xx
	 * limited to 6*/
	public static double[] round(double[] d, int numDecimals ){
		if (numDecimals<0 || numDecimals>6) return null;
		int[] multipliers = {1,10,100,1000,10000,100000,1000000};
		double[] rounded = new double[d.length];
		for (int i=d.length-1; i>=0; i--) {
			double comp = d[i]*multipliers[numDecimals];
			comp = Math.round(comp);
			rounded[i] = comp/multipliers[numDecimals];
		}
		return rounded;
	}

	/**Rounds a double based on number of desired decimals, ie 1 = x.x, 2 = x.xx
	 * limited to 6*/
	public static double round(double d, int numDecimals ){
		if (numDecimals<0 || numDecimals>6) return d;
		int[] multipliers = {1,10,100,1000,10000,100000,1000000};
		double comp = d*multipliers[numDecimals];
		comp = Math.round(comp);
		return comp/multipliers[numDecimals];
	}

	/**Removes any NaN or Infinite*/
	public static double[] removeNaNInfinite(double[] d){
		int num = d.length;
		boolean ok = true;
		for (int i=0; i<num; i++){
			if (Double.isNaN(d[i]) || Double.isInfinite(d[i]) ){
				ok = false;
				break;
			}
		}
		if (ok) return d;
		//not ok clear bad ones
		ArrayList a = new ArrayList(num);
		for (int i=0; i<num; i++){
			if (Double.isNaN(d[i]) == false && Double.isInfinite(d[i]) == false ){
				a.add(new Double(d[i]));
			}
		}
		return arrayListOfDoubleToArray(a);
	}

	/**Returns the LOG base 10 of the number.*/
	public static double log10(double number){
		return Math.log(number)/log10;
	}

	/**Returns the LOG base 2 of the number.*/
	public static double log2(double number){
		return Math.log(number)/log2;
	}

	/**Returns the LOG base 2 of the number.*/
	public static float log2(float number){
		return (new Double (Math.log((double)number)/log2)).floatValue();
	}

	/**Returns the -10 * LOG base 10 of the number.*/
	public static double minus10log10(double pvalue){
		return -10 * (Math.log10(pvalue));
	}
	
	/**Returns the -10 * LOG base 10 of the number.*/
	public static float minus10log10Float(double pvalue){
		return new Double(-10 * (Math.log10(pvalue))).floatValue();
	}

	/**Loads matrix of double[row number][double columns], white space delimited.
	 * # comment lines and blank lines ignored.
	 * */
	public static double[][] loadDoubleMatrix(File f){
		double[][] matrix = null;
		try{
			ArrayList al = new ArrayList();
			BufferedReader in = new BufferedReader ( new FileReader (f) );
			String line;
			String[] tokens;
			while ((line=in.readLine()) != null){
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#"))continue;
				tokens = line.split("\\s+");
				double[] tds = new double[tokens.length];
				for (int i=0; i< tds.length; i++) tds[i] = Double.parseDouble(tokens[i]);
				al.add(tds);
			}
			//convert to array
			int size = al.size();
			matrix = new double[size][];
			for (int i=0; i< size; i++) matrix[i]= (double[]) al.get(i);

		} catch (Exception e){
			e.printStackTrace();
		}
		return matrix;
	}
	
	/**Loads columns of double from a file into an array.
	 * Will substitute defaultValue if value = "Inf"*/
	public static double[][] loadDoubleMatrix(File f, Double defaultValue, int numberOfColumns){
		BufferedReader in ;
		double[][] values = null;
		Pattern whiteSpace = Pattern.compile("\\s+");
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList<double[]> al = new ArrayList<double[]>();
			String line;
			String[] tokens;
			while((line= in.readLine()) != null){
				tokens = whiteSpace.split(line);
				double[] d = new double[numberOfColumns];
				for (int i=0; i< numberOfColumns; i++) {
					if (tokens[i].contains("Inf")) d[i] = defaultValue;
					else d[i] = Double.parseDouble(tokens[i]);
				}
				al.add(d);
			}
			values = new double[al.size()][numberOfColumns];
			al.toArray(values);
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}
	
	/**Loads columns of double from a file into an array.
	 * Will substitute max value for Infs */
	public static double[][] loadDoubleMatrix(File f, int numberOfColumns){
		BufferedReader in ;
		double[][] values = null;
		Pattern whiteSpace = Pattern.compile("\\s+");
		double maxValue = Double.MIN_VALUE;
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList<double[]> al = new ArrayList<double[]>();
			String line;
			String[] tokens;
			while((line= in.readLine()) != null){
				tokens = whiteSpace.split(line);
				double[] d = new double[numberOfColumns];
				for (int i=0; i< numberOfColumns; i++) {
					if (tokens[i].contains("Inf")) d[i] = Double.MIN_VALUE;
					else {
						d[i] = Double.parseDouble(tokens[i]);
						if (d[i] > maxValue) maxValue = d[i];
					}
				}
				al.add(d);
			}
			values = new double[al.size()][numberOfColumns];
			al.toArray(values);
			//swap out Double.MIN_VALUEs for maxValue
			for (int i=0; i< values.length; i++){
				for (int j=0; j< values[i].length; j++){
					if (values[i][j] == Double.MIN_VALUE) values[i][j] = maxValue;
				}
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}
	
	/**Converts and negative doubles to zero.*/
	public static void zeroNegativeValues( double[][] d){
		for (int i=0; i< d.length; i++){
			for (int j=0; j< d[i].length; j++){
				if (d[i][j] < 0) d[i][j] = 0;
			}
		}
	}

	/**Loads a column of ints from a file into an array.
	 * Returns null if nothing found.
	 * Skips blank lines.*/
	public static int[] loadInts(File f){
		BufferedReader in ;
		int[] values = null;
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList al = new ArrayList();
			String line;
			while((line= in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				al.add(Integer.valueOf(line));
			}
			values = Num.arrayListOfIntegerToInts(al);
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}

	/**Loads a column of double from a file into an array.
	 * Returns null if nothing found.
	 * Will set 'Inf' values to max found.
	 * Skips blank lines.*/
	public static double[] loadDoubles(File f){
		BufferedReader in ;
		double[] values = null;
		double maxValue = Double.MIN_VALUE;
		Double defaultValue = new Double (Double.MIN_VALUE);
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList<Double> al = new ArrayList<Double>();
			String line;
			while((line= in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				if (line.contains("Inf")) al.add(defaultValue);
				else {
					Double val = Double.valueOf(line);
					if (val.doubleValue() > maxValue) maxValue = val.doubleValue();
					al.add(val);
				}
			}
			values = Num.arrayListOfDoubleToArray(al);
			//convert Inf values?
			if (maxValue != Double.MAX_VALUE){
				for (int i=0; i< values.length; i++){
					if (values[i] == Double.MIN_VALUE) values[i] = maxValue;
				}
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}
	
	/**Loads a column of double from a file into an array.
	 * Returns null if nothing found.
	 * Will set 'Inf' values to max found.
	 * Will set 'NA' values to either max or 0 based on the tcObs[][], max if tc[x][0] > tc[x][1]
	 * Assumes the file contains the same number of lines as the tcObs
	 * Skips blank lines.*/
	public static double[] loadDoubles(File f, int[][]tcObs, boolean leaveMaxValuesAsDouble_Min_Value){
		BufferedReader in ;
		double[] values = null;
		double maxValue = Double.MIN_VALUE;
		Double defaultValue = new Double (Double.MIN_VALUE);
		Double zero = new Double(0);
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList<Double> al = new ArrayList<Double>();
			String line;
			int index = 0;
			while((line= in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				if (line.equals("NA") || line.equals("Inf")){
					//System.out.println("BadLine\t"+line+"\t"+tcObs[index][0]+"\t"+tcObs[index][1]);
					if (tcObs[index][0] < tcObs[index][1]) al.add(zero);
					else al.add(defaultValue);
				}
				else {
					Double val = Double.valueOf(line);
					if (val.doubleValue() > maxValue) maxValue = val.doubleValue();
					al.add(val);
				}
				index++;
			}
			values = Num.arrayListOfDoubleToArray(al);
			//convert Inf and NaN values?
			if (leaveMaxValuesAsDouble_Min_Value == false && maxValue != Double.MAX_VALUE){
				for (int i=0; i< values.length; i++){
					if (values[i] == Double.MIN_VALUE) {
						values[i] = maxValue;
						//System.out.println("BadLine coverting to max->"+maxValue+"\t"+tcObs[i][0]+"\t"+tcObs[i][1]);
					}
				}
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}

	/**Loads a column of double from a file into an array.
	 * Returns null if nothing found.
	 * Skips blank lines. Will substitute defaultValue if value = "Inf"*/
	public static double[] loadDoubles(File f, Double defaultValue){
		BufferedReader in ;
		double[] values = null;
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList<Double> al = new ArrayList<Double>();
			String line;
			while((line= in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				if (line.equals("Inf")) al.add(defaultValue);
				else al.add(Double.valueOf(line));
			}
			values = Num.arrayListToDoubles(al);
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}

	/**Loads a column of float from a file into an array.
	 * Returns null if nothing found.
	 * Skips blank lines.*/
	public static float[] loadFloats(File f, float defaultValue){
		BufferedReader in ;
		float[] values = null;
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList al = new ArrayList();
			String line;
			while((line= in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				if (line.equals("Inf")) al.add(defaultValue);
				else al.add(Float.valueOf(line));
			}
			values = Num.arrayListOfFloatToArray(al);
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}
	
	/**Loads a column of float from a file into an array.
	 * Returns null if nothing found.
	 * Will set 'INF' values to max found.
	 * Skips blank lines.*/
	public static float[] loadFloats(File f){
		BufferedReader in ;
		float[] values = null;
		float maxValue = Float.MIN_VALUE;
		Float defaultValue = new Float (Float.MIN_VALUE);
		try{
			in = new BufferedReader (new FileReader(f));
			ArrayList<Float> al = new ArrayList<Float>();
			String line;
			while((line= in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				if (line.contains("Inf")) al.add(defaultValue);
				else {
					Float val = Float.valueOf(line);
					if (val.floatValue() > maxValue) maxValue = val.floatValue();
					al.add(val);
				}
			}
			values = Num.arrayListOfFloatToArray(al);
			//convert Inf values?
			if (maxValue != Float.MAX_VALUE){
				for (int i=0; i< values.length; i++){
					if (values[i] == Float.MIN_VALUE) values[i] = maxValue;
				}
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return values;
	}

	/**Writes a float array to file, one per line.*/
	public static boolean writeToFile (float[] floats, File floatFile){
		try{
			PrintWriter out = new PrintWriter (new FileWriter (floatFile));
			for (int i=0; i< floats.length; i++) out.println(floats[i]);
			out.close();
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}

	/**Writes a float array to file, one per line.*/
	public static boolean writeToFile (double[] doubles, File doubleFile){
		try{
			PrintWriter out = new PrintWriter (new FileWriter (doubleFile));
			for (int i=0; i< doubles.length; i++) out.println(doubles[i]);
			out.close();
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}

	/**ArrayList of Integer to int[]*/
	public static int[] arrayListOfIntegerToInts(ArrayList integers){
		int size = integers.size();
		int[] values = new int[size];
		for (int i=0; i< size; i++){
			values[i] = ((Integer)integers.get(i)).intValue();
		}
		return values;
	}

	/**ArrayList of Double to double[]*/
	public static double[] arrayListToDoubles(ArrayList doubles){
		int size = doubles.size();
		double[] values = new double[size];
		for (int i=0; i< size; i++){
			values[i] = ((Double)doubles.get(i)).doubleValue();
		}
		return values;
	}

	/**Attemps to parse an int, returns failNumber if it fails.*/
	public static int parseInt(String stringInt, int failNumber){
		try {
			return Integer.parseInt(stringInt);
		}catch(Exception e){
		}//leave empty
		return failNumber;
	}
	
	/**Returns a tab delimited string of Mean Median StdDev Min Max 10th 90th for the Float[]*/
	public static String statFloatArray(float[] sortedFloat){
		//calc mean
		float mean = Num.mean(sortedFloat); 
		//calc min, max
		float[] minMax = Num.findMinMaxFloatValues(sortedFloat);
		//calc median, need array copy to not sort referenced float!
		float median = (float)Num.median(sortedFloat);
		//calc 10th
		float perc10th = Num.percentile(sortedFloat, 0.1);
		//calc 90th
		float perc90th = Num.percentile(sortedFloat,0.9);
		//standard deviation
		float stndDev = (float)Num.standardDeviation(sortedFloat,mean);
		return mean+"\t"+median+"\t"+stndDev+"\t"+minMax[0]+"\t"+minMax[1]+"\t"+perc10th+"\t"+perc90th;

		
	}

	/**Calculates Min, Max, Mean, Median, Mode, and Histogram/10 for a Float[]*/
	public static void statFloatArray(float[] f, boolean printHistogram){
		//calc mean
		float mean = Num.mean(f); 
		//calc min, max
		float[] minMax = Num.findMinMaxFloatValues(f);
		//calc median, need array copy to not sort referenced float!
		float[] fSort = new float[f.length];
		System.arraycopy(f,0,fSort,0,f.length);
		Arrays.sort(fSort);
		float median = (float)Num.median(fSort);
		//calc 10th
		float perc10th = Num.percentile(fSort, 0.1);
		//calc 90th
		float perc90th = Num.percentile(fSort,0.9);
		fSort = null;
		//standard deviation
		float stndDev = (float)Num.standardDeviation(f,mean);

		//make histogram int[], cap at 10,000
		int lengthHistogram = (int)(10+minMax[1]-minMax[0]);
		if (lengthHistogram>10000) lengthHistogram = 10001;
		int[] histo = new int[lengthHistogram];
		int numF = f.length;
		int intensity;
		for (int i=0; i<numF; i++){
			intensity = Math.round(f[i]-minMax[0]);
			if (intensity>10000) intensity = 10000;
			histo[intensity]++;
		}
		//find mode
		int mode = (int)minMax[0]+ Num.modeOfHistogram(histo)[0];
		//report
		System.out.println("Mean\tMedian\tMode\tStdDev\tMin\tMax\t10th\t90th");
		System.out.println(mean+"\t"+median+"\t"+mode+"\t"+stndDev+"\t"+minMax[0]+"\t"+minMax[1]+"\t"+perc10th+"\t"+perc90th);

		if (printHistogram){
			System.out.println("A Histogram of values (capped at 10,000) ");
			System.out.println("Value\tFrequency\tStars");
			Histogram.printHistogram(histo);
		}
		histo = null;
	}

	/**Remove zero values from float array.*/
	public static float[] removeZeroValues(float[] f){
		int num = f.length;
		ArrayList al = new ArrayList(num);
		for (int i=0; i<num; i++){
			if (f[i] != 0) al.add(new Float(f[i]));
		}
		return Num.arrayListOfFloatToArray(al);
	}

	/**Gets the average of the integers bracketed and including the start and stop.
	 * (ie 3,6  returns the average of 3+4+5+6/4= 4.5)*/
	public static float getAverageInts (int start, int end){
		int endOne = end+1;
		int len = endOne-start;
		int sum =0;
		for (int i= start; i<endOne; i++) sum+=i;
		return (float)sum/(float)len;
	}

	/**Any length OK of the int[]s*/
	public static int[] collapseIntArray (int[][] ints){
		int total = Num.countValues(ints);
		int[] combine = new int[total];
		int k =0;
		for (int i=0; i< ints.length; i++){
			for (int j=0; j< ints[i].length; j++){
				combine[k++] = ints[i][j];
			}
		}
		return combine;
	}
	
	/**Takes an ArrayList of int[]s and returns a sorted int[] of unique values.*/
	public static int[] returnUniques (ArrayList<int[]> ints){
		int[][] x = new int[ints.size()][];
		for (int i=0; i< x.length; i++) x[i] = ints.get(i);
		return returnUniques(x);
	}
	
	/**Takes arrays of int[] and returns a sorted int[] of unique values.*/
	public static int[] returnUniques (int[][] ints){
		int[] merged = ints[0];
		//must merge first array against itself, then others
		for (int j=0; j< ints.length; j++){
			int[] combine = new int[merged.length+ints[j].length];
			System.arraycopy(merged, 0, combine, 0, merged.length);
			System.arraycopy(ints[j], 0, combine, merged.length, ints[j].length);
			Arrays.sort(combine);
			merged = returnUniques(combine);
		}
		return merged;
	}
	
	/**Takes a sorted array of int[] and returns unique values.*/
	public static int[] returnUniques (int[] sortedArray){
		ArrayList<Integer> al = new ArrayList<Integer>();
		int oldInt = sortedArray[0];
		for (int i=1; i<sortedArray.length; i++){
			if (oldInt != sortedArray[i]){
				al.add(new Integer(oldInt));
				oldInt = sortedArray[i];
			}
		}
		//add last
		al.add(new Integer(oldInt));
		return Num.arrayListOfIntegerToInts(al);
	}

	/**@param replicasIntensities - float[replicaNumber][associatedIntensities]
	 * @return comma delimited list of medians, one for each replica.*/
	public static String medianFloatArray (float[][] replicasIntensities){
		StringBuffer sb = new StringBuffer();
		Arrays.sort(replicasIntensities[0]);
		double median = Num.median(replicasIntensities[0]);
		sb.append(median);
		for (int i=1; i< replicasIntensities.length; i++){
			sb.append(",");
			Arrays.sort(replicasIntensities[i]);
			median = Num.median(replicasIntensities[i]);
			sb.append(median);
		}
		return sb.toString();
	}

	/**Assumes equal lengths of the float[]s*/
	public static float[] collapseFloatArray (float[][] fs){
		int numI = fs.length;
		int numJ = fs[0].length;
		float[] combine = new float[numI * numJ];
		int k =0;
		for (int i=0; i< numI; i++){
			for (int j=0; j< numJ; j++){
				combine[k++] = fs[i][j];
			}
		}
		return combine;
	}
	
	public static float[] collapseFloatArray (ArrayList<float[]> al){
		//count total num
		int len = 0;
		for (float[] f: al) len += f.length;
		float[] combine = new float[len];
		int index = 0;
		
		//for each array
		for (float[] f: al){
			System.arraycopy(f, 0, combine, index, f.length);
			index += f.length;
		}
		return combine;
	}

	/**Assumes equal lengths of the float[]s*/
	public static float[] collapsePartOfAnArray (float[][] fs, int startIndex, int stopIndex){
		int numI = fs.length;
		int numJ = fs[0].length;
		float[] combine = new float[numI * numJ];
		int k =0;
		for (int i=0; i< numI; i++){
			for (int j=0; j< numJ; j++){
				combine[k++] = fs[i][j];
			}
		}
		return combine;
	}

	/**Assumes equal lengths of the double[]s*/
	public static double[] collapseDoubleArray (double[][] fs){
		int numI = fs.length;
		int numJ = fs[0].length;
		double[] combine = new double[numI * numJ];
		int k =0;
		for (int i=0; i< numI; i++){
			for (int j=0; j< numJ; j++){
				combine[k++] = fs[i][j];
			}
		}
		return combine;
	}

	/**Converts a float[] to float[][], must be certain about the number of rows and colmns!*/
	public static float[][] expandFloatArray(float[] f, int numRows, int numColumns){
		float[][] e = new float[numRows][numColumns];
		int counter = 0;
		for (int i=0; i<numRows; i++){
			for (int j=0; j<numColumns; j++){
				e[i][j]= f[counter++];
			}
		}
		return e;
	}

	/**Averages a int[]*/
	public static double averageIntArray(int[] f){
		double tTot = 0;
		int lenTreat = f.length;
		for (int i=0; i< lenTreat; i++) tTot+=f[i];
		return (double)tTot/(double)lenTreat;
	}	

	/**Averages a float[]*/
	public static double averageFloatArray(float[] f){
		double tTot = 0;
		int lenTreat = f.length;
		for (int i=0; i< lenTreat; i++) tTot+=f[i];
		return (double)tTot/(double)lenTreat;
	}

	/**Finds min and max values of a int array.*/
	public static int[] findMinMaxIntValues(int[] f){
		int min = f[0];
		int max = f[0];
		int num = f.length;
		for (int i=0; i<num; i++){
			if (f[i]< min) min=f[i];
			if (f[i]> max) max=f[i];
		}
		return new int[]{min,max};
	}


	/**Finds min and max values of an unsorted float array.*/
	public static float[] findMinMaxFloatValues(float[] f){
		float min = f[0];
		float max = f[0];
		int num = f.length;
		for (int i=0; i<num; i++){
			if (f[i]< min) min=f[i];
			if (f[i]> max) max=f[i];
		}
		return new float[]{min,max};
	}
	/**Finds min and max values of a float array.*/
	public static double[] findMinMaxDoubleValues(double[] f){
		double min = f[0];
		double max = f[0];
		int num = f.length;
		for (int i=0; i<num; i++){
			if (f[i]< min) min=f[i];
			if (f[i]> max) max=f[i];
		}
		return new double[]{min,max};
	}	

	/**Finds min and max values of a int array.*/
	public static int[] findMinMaxIntInt(int[][] f){
		int min = f[0][0];
		int max = f[0][0];
		int num = f.length;
		int[] minMax;
		for (int i=0; i<num; i++){
			minMax = findMinMaxIntValues(f[i]);
			if (minMax[0]<min) min = minMax[0];
			if (minMax[1]>max) max = minMax[1];
		}
		return new int[]{min,max};
	}

	/**Finds min and max values of a float array.*/
	public static float[] findMinMaxFloatArrays(float[][] f){
		float min = f[0][0];
		float max = f[0][0];
		int num = f.length;
		float[] minMax;
		for (int i=0; i<num; i++){
			minMax = findMinMaxFloatValues(f[i]);
			if (minMax[0]<min) min = minMax[0];
			if (minMax[1]>max) max = minMax[1];
		}
		return new float[]{min,max};
	}

	/**Finds the index of the key in the int[] but wont preceed the int[index]*/
	public static int findClosestStartIndex (int[] positions, int key){
		int index = Arrays.binarySearch(positions,key);
		if (index <0 ) {
			return (index+1) * -1;
		}
		return index;
	}

	/**Finds the index of the key in the int[] but wont exceed the int[index]*/
	public static int findClosestEndIndex (int[] positions, int key){
		int index = Arrays.binarySearch(positions,key);
		if (index <0 ) {
			return (index+2) * -1;
		}
		return index;
	}

	/**Returns the index of the closest values[index] to the key. Rounds up
	 * if the value is equi distant between two indexes.  Will not return an
	 * index <0 or > length-1 .
	 * Don't forget to SORT!.
	 * If identical values are present, returns the smallest index containing the key.*/
	public static int findClosestIndexToValue(int[] sortedValues, int key){
		int index = Arrays.binarySearch(sortedValues,key);
		//no exact match
		if (index<0){
			int indexAfter = (-1* index) -1;
			if (indexAfter >= sortedValues.length) return sortedValues.length-1;
			int indexBefore = indexAfter -1;
			if (indexBefore<0) return 0;
			int diffToAfter = Math.abs(sortedValues[indexAfter]-key);
			int diffToBefore = Math.abs(sortedValues[indexBefore]-key);
			if (diffToAfter >  diffToBefore) return indexBefore;
			else return indexAfter;
		}
		//exact match look for smaller indexes with same value
		else {
			float value = sortedValues[index];
			int testIndex = index;
			while (true){
				testIndex--;
				//look for negative index
				if (testIndex< 0) {
					index = 0;
					break;
				}
				//if different value break
				if (value != sortedValues[testIndex]) {
					index = ++testIndex;
					break;
				}
			}
			return index;
		}
	}

	/**Returns the index of the closest values[index] to the key. Rounds up
	 * if the value is equi distant between two indexes.  Will not return an
	 * index <0 or > length-1 .
	 * Don't forget to SORT!.
	 * If identical values are present, returns the smallest index containing the key.*/
	public static int findClosestIndexToValue(float[] sortedValues, float key){
		int index = Arrays.binarySearch(sortedValues,key);
		//no exact match
		if (index<0){
			int indexAfter = (-1* index) -1;
			if (indexAfter >= sortedValues.length) return sortedValues.length-1;
			int indexBefore = indexAfter -1;
			if (indexBefore<0) return 0;
			float diffToAfter = Math.abs(sortedValues[indexAfter]-key);
			float diffToBefore = Math.abs(sortedValues[indexBefore]-key);
			if (diffToAfter >  diffToBefore) return indexBefore;
			else return indexAfter;
		}
		//exact match look for smaller indexes with same value
		else {
			float value = sortedValues[index];
			int testIndex = index;
			while (true){
				testIndex--;
				//look for negative index
				if (testIndex< 0) {
					index = 0;
					break;
				}
				//if different value break
				if (value != sortedValues[testIndex]) {
					index = ++testIndex;
					break;
				}
			}
			return index;
		}
	}

	/**Finds the start and stop indexes given a sorted int[] of values and two values, one after the other, 
	 * close together, otherwise just use Arrays.binarySearch.  Good for huge int[] arrays. Set the notExact 
	 * boolean to true if the values are potentially not exact.*/
	public static int[] findStartStopIndexes(int[] positions, int startValue, int stopValue, boolean notExact){
		int startIndex;
		int stopIndex;
		if (notExact){
			startIndex = findClosestIndexToValue(positions, startValue);
			stopIndex = findClosestIndexToValue(positions, stopValue);
		}
		else {
			startIndex = Arrays.binarySearch(positions,startValue);
			stopIndex = startIndex;
			int length = positions.length-1;
			while (stopIndex!=length){
				if (stopValue == positions[stopIndex]) break;
				stopIndex++;
			}
		}
		return new int[]{startIndex, stopIndex};
	}

	/**Calculates the standard deviation of the difference between two int[]s 
	 * of the same length.*/
	public static double standardDeviationOfDifference(int[] t, int[] c){
		int num = t.length;
		int[] diff = new int[num];
		for (int i=0; i<num; i++){
			diff[i]= t[i]-c[i];
		}
		return standardDeviation(diff);
	}

	/**Average int[][] values to double[], averages repeats [numRepeats][values]*/ 
	public static double[] averageIntIntArray(int[][] ints){
		int num = ints.length;
		double[] ave = new double[num];
		double numSub;
		double sum;
		for (int i=0; i<num; i++){
			numSub = ints[i].length;
			sum =0;
			for (int j=0; j<numSub; j++){
				sum += ints[i][j];
			}
			ave[i] = sum/numSub;
		}
		return ave;
	}

	/**Average float[][] values to double[], averages repeats [values][numRepeats]*/ 
	public static double[] averageFloatArrays(float[][] ints){
		int num = ints.length;
		double[] ave = new double[num];
		double numSub;
		double sum;
		for (int i=0; i<num; i++){
			numSub = ints[i].length;
			sum =0;
			for (int j=0; j<numSub; j++){
				sum += ints[i][j];
			}
			ave[i] = sum/numSub;
		}
		return ave;
	}

	/**Average float[][] values to double[], averages repeats [numChips][OligoValues].
	 * Assumes equal number of oligos.*/ 
	public static double[] averageFloatArraysFlipped(float[][] ints){
		double numChips = ints.length;
		int numOligos = ints[0].length;
		double[] ave = new double[numOligos];
		double sum;
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//total the repeats
			for (int j=0; j<numChips; j++){
				sum += ints[j][i];
			}
			ave[i] = sum/numChips;
		}
		return ave;
	}

	/**Average float[][] values to double[], averages repeats [numChips][OligoValues].
	 * Assumes equal number of oligos.*/ 
	public static float[] averageFloatArraysFlippedToFloat(float[][] ints){
		float numChips = ints.length;
		int numOligos = ints[0].length;
		float[] ave = new float[numOligos];
		float sum;
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//total the repeats
			for (int j=0; j<numChips; j++){
				sum += ints[j][i];
			}
			ave[i] = sum/numChips;
		}
		return ave;
	}

	/**Takes the geometric mean of float[][] values to double[], [numChips][OligoValues], returns the values
	 * in log base 2.
	 * Assumes equal number of oligos.*/ 
	public static double[] geoMeanFloatArraysFlipped(float[][] ints){
		int numChips = ints.length;
		int numOligos = ints[0].length;
		double[] ave = new double[numOligos];
		double sum;
		//for each oligo
		double log2 = Math.log(2);
		for (int i=0; i<numOligos; i++){
			sum =0;
			//total the repeats
			for (int j=0; j<numChips; j++){
				sum += Math.log(ints[j][i])/log2;
			}
			ave[i] = sum/(double)numChips;
		}
		return ave;
	}

	/**Average float[][] values to float[], averages repeats [numChips][OligoValues].
	 * Assumes equal number of oligos.*/ 
	public static float[] averageFloatArraysFlippedToFloats(float[][] ints){
		float numChips = ints.length;
		int numOligos = ints[0].length;
		float[] ave = new float[numOligos];
		float sum;
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//total the repeats
			for (int j=0; j<numChips; j++){
				sum += ints[j][i];
			}
			ave[i] = sum/numChips;
		}
		return ave;
	}

	/**Takes a geometric average of  float[][] values to double[], averages repeats [numChips][OligoValues].
	 * Assumes equal number of oligos. Skips zero values. Return in log base 10.*/ 
	public static double[] geometricMeanSkipZeros(float[][] ints){
		int numChips = ints.length;
		int numOligos = ints[0].length;
		double[] ave = new double[numOligos];
		//for each oligo
		for (int i=0; i<numOligos; i++){
			double sum =0;
			double numLogs = 0;
			//sum the non zero intensities, all intensities >0 unless masked
			for (int j=0; j<numChips; j++){
				if (ints[j][i]>0) {
					sum += Math.log(ints[j][i]);
					numLogs++;
				}
			}
			//find mean of logs
			if (numLogs !=0){
				double meanLogs = sum/numLogs;
				ave[i]= Math.pow(Math.E, meanLogs);
			}
			else ave[i] = 0;
		}
		return ave;
	}

	/**Takes a geometric average of float[][] values [numChips][OligoValues]. 
	 * Assumes no intensities <= 0.
	 * Returns log base 10 values.*/ 
	public static double[] geometricMean(float[][] ints){
		int numChips = ints.length;
		int numOligos = ints[0].length;
		double[] ave = new double[numOligos];
		//double zeroLog = Math.log(0.0000000001);
		//for each oligo
		for (int i=0; i<numOligos; i++){
			double sum =0;
			//sum the logged intensities
			for (int j=0; j<numChips; j++){
				//if (ints[j][i]<=0 ) sum+= zeroLog;
				//else sum += Math.log(ints[j][i]);
				sum += Math.log(ints[j][i]);
			}
			//find mean of logs
			double meanLogs = sum/numOligos;
			ave[i]= Math.pow(Math.E, meanLogs);

		}
		return ave;
	}




	/**Average int[][] values to double[], averages values [numRepeats][values]*/ 
	public static double[] averageIntIntArray2(int[][] ints){
		int numOligos = ints[0].length;
		double[] ave = new double[numOligos];
		double numIntensities = ints.length;
		double sum;
		//for each oligo
		for (int i=0; i<numOligos; i++){
			//average the intensity values for the oligo
			sum =0;
			for (int j=0; j<numIntensities; j++){
				sum += ints[j][i];
			}
			ave[i] = sum/numIntensities;
		}
		return ave;
	}

	/**Uses a sliding window to smooth the scores applying a trimmed mean (drop 1 from ends) to the scores found within the window.
	 * The collection of values within the window is only an attempt to look windowSize/2 upstream and
	 * windowSize/2 downstream, not a requirement.  Thus at either stop of the scores possibly only 1/2 of
	 * the scores are used in making the calculation that is associated with the original score.*/
	public static double[] smoothScores(double[] scores, int[] positions, int windowSize){
		int num = scores.length;
		double[] smoothedScores = new double[num];
		int halfWindow = (int)Math.round((int)windowSize/2);
		for (int i=0; i< num; i++){
			ArrayList vals = new ArrayList();
			vals.add(new Double(scores[i]));
			//double total = scores[i];
			//double numPts =1;
			int start = positions[i]-halfWindow;
			int stop = positions[i]+ halfWindow;
			//take start, attempt to look up stream xbp and downstream x bp adding scores
			//look upstream
			boolean go = true;
			int advanceIndex = i;
			while (go){
				//attempt to advance one
				advanceIndex++;
				if (advanceIndex<num){
					//check position
					if (positions[advanceIndex]<=stop){
						//not too far thus add val and increment numPts
						//total += scores[advanceIndex];
						//numPts++;
						vals.add(new Double(scores[advanceIndex]));
					}
					else{
						//too far bp wise
						go = false;
					}
				}
				else{
					//at stop of interval cannot advance
					go = false;
				}
			}
			//look downstream
			go = true;
			advanceIndex = i;
			while (go){
				//attempt to back up one
				advanceIndex--;
				if (advanceIndex>-1){
					//check position
					if (positions[advanceIndex]>=start){
						//not too far thus add val and increment numPts
						//total += scores[advanceIndex];
						//numPts++;
						vals.add(new Double(scores[advanceIndex]));
					}
					else{
						//too far bp wise
						go = false;
					}
				}
				else{
					//at stop of interval cannot advance
					go = false;
				}
			}

			//make average,
			//smoothedScores[i]= total/numPts; //mean

			//median
			double[] fin = arrayListOfDoubleToArray(vals);
			Arrays.sort(fin);
			//smoothedScores[i] = Num.median(fin);

			//trimmed mean, watch out for too few in array, occurs at ends.
			if (fin.length>2) smoothedScores[i] = Num.trimmedMean(fin,1);
			else smoothedScores[i]= Num.mean(fin);

		}
		return smoothedScores;
	}
	/**Calculates a mean on the non NaN and Infinity floats, return zero if no non zero values were found.*/
	public static double meanIgnoreNaNInfinity(float[] f){
		double numVals = 0;
		double total =0;
		for (int i=0; i< f.length; i++){
			if (Float.isNaN(f[i]) == false && Float.isInfinite(f[i]) == false) {
				total += f[i];
				numVals++;
			}
		}
		if (numVals>0) return total/numVals;
		return 0;
	}
	/** Slides a window along an array, one index at a time, averaging the contents. */
	public static double[] windowAverageScores(double[] scores, int windowSize) {
		double window = windowSize;
		int num = 1+ scores.length - windowSize;
		if (num < 1) return null;
		double[] means = new double[num];
		for (int i=0; i<num; i++){
			//sum window
			double stop = i+window;
			double total = 0;
			for (int j=i; j< stop; j++){
				total += scores[j];
			}
			means[i] = total/window;
		}
		return means;
	}

	/** Slides a window along an array, one index at a time, averaging the contents, ignores zero values. */
	public static double[] windowAverageScoresIgnoreZeros(double[] scores, int windowSize) {
		double window = windowSize;
		int num = 1+ scores.length - windowSize;
		if (num < 1) return null;
		double[] means = new double[num];
		for (int i=0; i<num; i++){
			//sum window
			double stop = i+window;
			double total = 0;
			double count = 0;
			for (int j=i; j< stop; j++){
				if (scores[j] !=0){
					total += scores[j];
					count ++;
				}
			}
			means[i] = total/count;
		}
		return means;
	}

	/** Slides a window along an array, one index at a time, averaging the contents, ignores zero values. */
	public static float[] windowAverageScoresIgnoreScore(float[][] scores, float scoreToIgnore, int windowSize) {
		float window = windowSize;
		int num = 1+ scores[0].length - windowSize;
		float[] means = new float[num];
		for (int i=0; i<num; i++){
			//sum window
			float stop = i+window;
			float total = 0;
			float count = 0;
			//for each window
			for (int j=i; j< stop; j++){
				//for each score index look for non scoresToIgnore and add
				for (int k = 0; k< scores.length; k++){
					if (scores[k][j] != scoreToIgnore){
						total += scores[k][j];
						count ++;
					}
				}
			}
			if (count == 0) means[i] = scoreToIgnore;
			else means[i] = total/count;
		}
		return means;
	}

	/** Slides a window along an array, one index at a time, averaging the contents. */
	public static float[] windowAverageScoresIgnoreScore(float[] scores, float scoreToIgnore, int windowSize) {
		float window = windowSize;
		int num = 1+ scores.length - windowSize;
		if (num < 1) return null;
		float[] means = new float[num];
		for (int i=0; i<num; i++){
			//sum window
			float stop = i+window;
			float total = 0;
			float numScores = 0;
			for (int j=i; j< stop; j++){
				if (scores[j] != scoreToIgnore) {
					total += scores[j];
					numScores++;
				}
			}
			if (numScores !=0) means[i] = total/numScores;
			else means[i] = scoreToIgnore;
		}
		return means;
	}	

	/** Slides a window along an array, one index at a time, averaging the contents. */
	public static float[] windowAverageScores(float[] scores, int windowSize) {
		double window = windowSize;
		int num = 1+ scores.length - windowSize;
		if (num < 1) return null;
		float[] means = new float[num];
		for (int i=0; i<num; i++){
			//sum window
			double stop = i+window;
			double total = 0;
			for (int j=i; j< stop; j++){
				total += scores[j];
			}
			means[i] = new Double(total/window).floatValue();
		}
		return means;
	}

	/** Slides a window along an array, one index at a time, averaging the contents. 
	 * Scans all windows including start and ends. Note, if windowSize is odd then the scan size is actually windowSize-1*/
	public static float[] windowAverageScoresNoShrink(float[] scores, int windowSize) {
		if (windowSize == 1) return scores;
		double window = windowSize;
		int halfWindow = (int)(window/2);
		float[] means = new float[scores.length];
		for (int i=0; i<scores.length; i++){
			//set start and stop
			int start = i-halfWindow;
			if (start < 0) start = 0;
			int stop = i+halfWindow;
			if (stop > scores.length) stop = scores.length;
			//sum window
			double total = 0;
			for (int j=start; j< stop; j++){
				total += scores[j];
			}
			double num = stop - start;
			means[i] = new Double(total/num).floatValue();
			if (Float.isNaN(means[i])) means[i] = 0;
		}
		return means;
	}

	/** Slides a window along an array, one index at a time, averaging the contents. 
	 * Scans all windows including start and ends. Note, if windowSize is odd then the scan size is actually windowSize-1.
	 * Zero values are ignored in mean calculation.*/
	public static double[] windowAverageScoresNoShrinkIgnoreZeros(double[] scores, int windowSize) {
		double window = windowSize;
		int halfWindow = (int)(window/2);
		double[] means = new double[scores.length];
		for (int i=0; i<scores.length; i++){
			//set start and stop
			int start = i-halfWindow;
			if (start < 0) start = 0;
			int stop = i+halfWindow;
			if (stop > scores.length) stop = scores.length;
			//sum window
			double total = 0;
			double counter = 0;
			for (int j=start; j< stop; j++){
				if (scores[j] !=0){
					total += scores[j];
					counter++;
				}
			}
			if (counter !=0) means[i] = total/counter;
		}
		return means;
	}
	
	/** Slides a window along an array, one index at a time, calculating a trimmed mean on the contents. 
	 * Scans all windows including start and ends. Note, if windowSize is odd then the scan size is actually windowSize-1.
	 * Zero values are ignored in mean calculation.*/
	public static double[] windowTrimmedMeanScoresNoShrinkIgnoreZeros(double[] scores, int windowSize) {
		double window = windowSize;
		int halfWindow = (int)(window/2);
		double[] means = new double[scores.length];
		ArrayList<Double> subScoresAL = new ArrayList<Double>();
		double oldScore = 0;
		for (int i=0; i<scores.length; i++){
			//set start and stop
			int start = i-halfWindow;
			if (start < 0) start = 0;
			int stop = i+halfWindow;
			if (stop > scores.length) stop = scores.length;
			//collect scores
			 subScoresAL.clear();
			for (int j=start; j< stop; j++){
				if (scores[j] !=0){
					subScoresAL.add(new Double(scores[j]));
				}
			}
			//check number
			int numScores = subScoresAL.size();
			if (numScores == 0) means[i] = oldScore;
			else if (numScores == 1) means[i] = subScoresAL.get(0).doubleValue();
			else if (numScores == 2){
				double[] subScores = Num.arrayListOfDoubleToArray(subScoresAL);
				means[i] = Num.mean(subScores);
			}
			else {
				double[] subScores = Num.arrayListOfDoubleToArray(subScoresAL);
				Arrays.sort(subScores);
				means[i] = Num.trimmedMean(subScores, 1);
			}
			oldScore = means[i];
		}
		return means;
	}
	
	/** Slides a window along an array, one index at a time, calculating a median on the contents. 
	 * Scans all windows including start and ends. Note, if windowSize is odd then the scan size is actually windowSize-1.
	 * Zero values are ignored in median calculation.*/
	public static double[] windowMedianScoresNoShrinkIgnoreZeros(double[] scores, int windowSize) {
		double window = windowSize;
		int halfWindow = (int)(window/2);
		double[] means = new double[scores.length];
		ArrayList<Double> subScoresAL = new ArrayList<Double>();
		double oldScore = 0;
		for (int i=0; i<scores.length; i++){
			//set start and stop
			int start = i-halfWindow;
			if (start < 0) start = 0;
			int stop = i+halfWindow;
			if (stop > scores.length) stop = scores.length;
			//collect scores
			 subScoresAL.clear();
			for (int j=start; j< stop; j++){
				if (scores[j] !=0){
					subScoresAL.add(new Double(scores[j]));
				}
			}
			//check number
			int numScores = subScoresAL.size();
			if (numScores == 0) means[i] = oldScore;
			else if (numScores == 1) means[i] = subScoresAL.get(0).doubleValue();
			else if (numScores == 2){
				double[] subScores = Num.arrayListOfDoubleToArray(subScoresAL);
				means[i] = Num.mean(subScores);
			}
			else {
				double[] subScores = Num.arrayListOfDoubleToArray(subScoresAL);
				Arrays.sort(subScores);
				means[i] = Num.median(subScores);
			}
			oldScore = means[i];
		}
		return means;
	}

	/**Trims the fraction off the top and bottom, returns the mean of the remainder, rounds trimming # up.
	 * If trim number resolves to 0 will proceed returning mean.*/
	public static double trimmedMean(double[] sortedFloat, double fraction){
		int trimNumber = (int)Math.round(((double)sortedFloat.length) * fraction);
		return trimmedMean(sortedFloat, trimNumber);
	}

	/**Trims the fraction off the top and bottom, returns the mean of the remainder, rounds trimming # up.
	 * If trim number resolves to 0 will proceed returning mean.*/
	public static double trimmedMean(float[] sortedFloat, double fraction){
		int trimNumber = (int)Math.round(((double)sortedFloat.length) * fraction);
		if (trimNumber == 0) trimNumber = 1;
		return trimmedMean(sortedFloat, trimNumber);
	}

	/**Makes a trimmed mean, note the trim number is not a % but the number of values to drop
	 * from the beginning and stop of the ordered set. Set appropriately. 
	 * Be sure to submit an ordered array!*/
	public static double trimmedMean(double[] sortedDouble, int trimNumber){
		double num = sortedDouble.length - trimNumber;
		double total = 0;
		for (int i=trimNumber; i<num; i++){
			total += sortedDouble[i];
		}
		return total/(sortedDouble.length-(2*trimNumber));
	}

	/**Makes a trimmed mean, note the trim number is not a % but the number of values to drop
	 * from the beginning and stop of the ordered set. Set appropriately. 
	 * Be sure to submit an ordered array!*/
	public static double trimmedMean(float[] sortedFloat, int trimNumber){
		double num = sortedFloat.length - trimNumber;
		double total = 0;
		for (int i=trimNumber; i<num; i++){
			total += sortedFloat[i];
		}
		return total/(sortedFloat.length-(2*trimNumber));
	}

	/**ArrayList of Double to double[]*/
	public static double[] arrayListOfDoubleToArray(ArrayList dbl){
		int num = dbl.size();
		double[] d = new double[num];
		for (int i=0; i<num; i++)d[i]= ((Double)dbl.get(i)).doubleValue();
		return d;
	}
	/**ArrayList of Float to float[]*/
	public static float[] arrayListOfFloatToArray(ArrayList flt){
		int num = flt.size();
		float[] d = new float[num];
		for (int i=0; i<num; i++)d[i]= ((Float)flt.get(i)).floatValue();
		return d;
	}

	/**Calculates the average fold difference between two int[]s 
	 * of the same length.  Takes the mean of both arrays returns
	 * the ratio.*/
	public static double aveFoldDiffCombine(int[] t, int[] c){
		double tMean = mean(t);
		double cMean = mean(c);
		return tMean/cMean;
	}
	/**Calculates the average fold difference between two int[]s 
	 * of the same length. Calculates the ratio of each individually,
	 * then returns their average.*/
	public static double aveFoldDiffIndividual(int[] t, int[] c){
		int num = t.length;
		double[] ratios = new double[num];
		for (int i=0; i<num; i++){
			ratios[i] = (double)t[i]/(double)c[i];
		}
		return mean(ratios);
	}

	/**Returns the ratio of each pair.*/
	public static double[] ratio(double[] t, double[] c){
		int num = t.length;
		double[] ratios = new double[num];
		for (int i=0; i<num; i++){
			ratios[i] = t[i]/c[i];
		}
		return ratios;
	}

	/**Returns the ratio of each pair, averages replicas.
	 * @param float[replicas][intensities]*/
	public static float[] ratio(float[][] t, float[][] c){
		int numOligos = t[0].length;	
		float ratios[] = new float[numOligos];
		float numChipsT = t.length;
		float aveT;
		float numChipsC = c.length;
		float aveC;
		float sum;
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//average treatments
			for (int j=0; j<numChipsT; j++){
				sum += t[j][i];
			}
			aveT = sum/numChipsT;
			sum = 0;
			//average controls
			for (int j=0; j<numChipsC; j++){
				sum += c[j][i];
			}
			aveC = sum/numChipsC;
			//calc ratio
			ratios[i] = aveT/aveC;
		}
		return ratios;
	}

	/**Returns the difference of each pair of average, aveT - aveC.
	 * @param float[replicas][intensities]*/
	public static float[] difference(float[][] t, float[][] c){
		int numOligos = t[0].length;	
		float ratios[] = new float[numOligos];
		float numChipsT = t.length;
		float aveT;
		float numChipsC = c.length;
		float aveC;
		float sum;
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//average treatments
			for (int j=0; j<numChipsT; j++){
				sum += t[j][i];
			}
			aveT = sum/numChipsT;
			sum = 0;
			//average controls
			for (int j=0; j<numChipsC; j++){
				sum += c[j][i];
			}
			aveC = sum/numChipsC;
			//calc diff
			ratios[i] = aveT - aveC;
		}
		return ratios;
	}

	/**Returns the relative difference of each pair, averages replicas.
	 * @param float[replicas][intensities]*/
	public static float[] relativeDifferences(float[][] t, float[][] c){
		int numOligos = t[0].length;	
		float diffs[] = new float[numOligos];
		double numChipsT = t.length;
		double aveT;
		double numChipsC = c.length;
		double aveC;
		double sum;
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//average treatments
			for (int j=0; j<numChipsT; j++){
				sum += t[j][i];
			}
			aveT = sum/numChipsT;
			sum = 0;
			//average controls
			for (int j=0; j<numChipsC; j++){
				sum += c[j][i];
			}
			aveC = sum/numChipsC;
			//calc relative difference
			diffs[i] = new Double(2 * (  (aveT - aveC) / (aveT + aveC)   )).floatValue();
		}
		return diffs;
	}

	/**Returns the log2 ratio of each pair, takes a geometric mean of replicas.
	 * @param float[replicas][intensities]*/
	public static float[] geometricMeanRatio(float[][] t, float[][] c){
		int numOligos = t[0].length;	
		float ratios[] = new float[numOligos];
		double numChipsT = t.length;
		double aveT;
		double numChipsC = c.length;
		double aveC;
		double sum;
		double log2 = Math.log(2);
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//average treatments
			for (int j=0; j<numChipsT; j++){
				sum += Math.log(t[j][i])/log2;
			}
			aveT = sum/numChipsT;
			sum = 0;
			//average controls
			for (int j=0; j<numChipsC; j++){
				sum += Math.log(c[j][i])/log2;
			}
			aveC = sum/numChipsC;
			//calc ratio
			ratios[i] = new Double(aveT - aveC).floatValue();
		}
		return ratios;
	}

	/**Returns the log2 ratio of each pair, averaging replicas.
	 * @param float[replicas][intensities]*/
	public static float[] logRatios(float[][] t, float[][] c){
		int numOligos = t[0].length;	
		float ratios[] = new float[numOligos];
		float numChipsT = t.length;
		float aveT;
		float numChipsC = c.length;
		float aveC;
		float sum;
		double log2 = Math.log(2);
		//for each oligo
		for (int i=0; i<numOligos; i++){
			sum =0;
			//average treatments
			for (int j=0; j<numChipsT; j++){
				sum += t[j][i];
			}
			aveT = sum/numChipsT;
			sum = 0;
			//average controls
			for (int j=0; j<numChipsC; j++){
				sum += c[j][i];
			}
			aveC = sum/numChipsC;
			//calc ratio
			ratios[i] = new Double (Math.log(aveT/aveC) / log2).floatValue();
		}
		return ratios;
	}

	/**Returns the ratio of each pair.*/
	public static double[] ratio(float[] t, float[] c){
		int num = t.length;
		double[] ratios = new double[num];
		for (int i=0; i<num; i++){
			ratios[i] = (double)t[i]/ (double)c[i];
		}
		return ratios;
	}

	/**Returns the ratio of each pair.*/
	public static double[] ratio(int[] t, double[] c){
		int num = t.length;
		double[] ratios = new double[num];
		for (int i=0; i<num; i++){
			ratios[i] = ((double)t[i])/c[i];
		}
		return ratios;
	}
	/**Returns the ratio of each pair.*/
	public static double[] ratio(int[] t, int[] c){
		int num = t.length;
		double[] ratios = new double[num];
		for (int i=0; i<num; i++){
			ratios[i] = (double)t[i]/(double)c[i];
		}
		return ratios;
	}

	/**Returns the difference of each pair.*/
	public static double[] difference(int[] t, int[] c){
		int num = t.length;
		double[] diffs = new double[num];
		for (int i=0; i<num; i++){
			diffs[i] = t[i]-c[i];
		}
		return diffs;
	}
	/**Returns the difference of each pair.*/
	public static double[] difference(double[] t, double[] c){
		int num = t.length;
		double[] diffs = new double[num];
		for (int i=0; i<num; i++){
			diffs[i] = t[i]-c[i];
		}
		return diffs;
	}

	/** Builds a pooled from part of a float[replica][intensities] for treat and control.
	 * @param treatment and control float[replicas][intensities]
	 * @param start and stop indexes are inclusive
	 * @return float[treatment(0) or control(1)][pooled intensities]*/
	public static float[][] pooledSubArray(float[][] treatment, float[][] control, int startIndex, int stopIndex){
		int numIntensities = 1+stopIndex-startIndex;
		int numTreatments = treatment.length * numIntensities;
		float[] ts = new float[numTreatments];

		int numControls = control.length * numIntensities;
		float[] cs = new float[numControls];

		//for each intensity/ oligo
		int counter = 0;
		for (int j=startIndex; j<=stopIndex; j++){
			//make treatments
			for (int i=0; i<treatment.length; i++) { 
				ts[counter++] = treatment[i][j];
			}
			counter = 0;
			for (int i=0; i<control.length; i++) {
				cs[counter++] = control[i][j];
			}
		}
		float[][] tc = new float[2][];
		tc[0] = ts;
		tc[1] = cs;
		return tc;
	}

	/**Calculates the SAM d statistic, similar to a t-statistic with a bit of fudge added into the denominator to
	 * control for zero variance conditions. Set fudge to zero for standard unpaired T-Test.
	 * @return float[3] samScore, diff of means, variance (w/o fudge).*/
	public static float[] sam(float[] t, float[] c, double fudge){
		double meanT = mean(t);
		double meanC = mean(c);
		double top = meanT - meanC;

		double sT = Math.pow(standardDeviation(t, meanT), 2)/(double)t.length;
		double sC = Math.pow(standardDeviation(c, meanC), 2)/(double)c.length;
		double var = Math.sqrt(sT+sC);

		double sam = top/(var + fudge);
		float[] scores = new float[3];
		scores[0] = new Double(sam).floatValue();
		scores[1] = new Double(top).floatValue();
		scores[2] = new Double(var).floatValue();
		return scores;
	}

	/**Calculates standard deviation of an int[]*/
	public static double standardDeviation(int[] x){
		double mean = mean(x);
		return standardDeviation(x, mean);
	}
	/**Calculates standard deviation of an int[] and mean.*/
	public static double standardDeviation(int[] x, double mean){
		double N = x.length;
		double topTot =0;
		for (int i=0; i<N; i++){
			topTot += Math.pow(x[i]-mean,2);
		}
		return Math.sqrt( topTot/(N-1) );
	}

	/**Calculates standard deviation of an float[]*/
	public static double standardDeviation(float[] x){
		double mean = mean(x);
		return standardDeviation(x, mean);
	}
	/**Calculates standard deviation of an float[] given a mean.*/
	public static double standardDeviation(float[] x, double mean){
		double N = x.length;
		double topTot =0;
		for (int i=0; i<N; i++){
			topTot += Math.pow((double)x[i]-mean,2);
		}
		return Math.sqrt( topTot/(N-1) );
	}



	/**Calculates standard deviation of an double[]*/
	public static double standardDeviation(double[] x){
		double mean = mean(x);
		return standardDeviation(x, mean);
	}
	/**Calculates standard deviation of a double[] and mean.*/
	public static double standardDeviation(double[] x, double mean){
		double N = x.length;
		double topTot =0;
		for (int i=0; i<N; i++){
			topTot += Math.pow(x[i]-mean,2);
		}
		return Math.sqrt( topTot/(N-1) );
	}	

	/**Calculates standard error of a double[] given it's mean.
	 * Spead savings*/
	public static double standardError(double[] x, double mean){
		double N = x.length;
		double topTot =0;
		for (int i=0; i<N; i++){
			topTot += Math.pow(x[i]-mean,2);
		}
		return Math.sqrt( topTot/(N-1) )/ Math.sqrt(N);
	}

	/**Calculates the median absolute difference.*/
	public static double medianAbsoluteDifference(float[] x, float[] y){
		int num = x.length;
		float[] diff = new float[num];
		for (int i=0; i<num; i++){
			diff[i] = Math.abs(x[i]-y[i]);
		}
		Arrays.sort(diff);
		return median(diff);
	}

	/**Calculates the maximum median absolute difference between arrays of float where int[oligo index number][oligo intensity measurements]*/
	public static double calcMaxMedianAbsoluteDifference(float[][] oligosValues){
		//calculate medianAbsDiff for each pair
		double maxMAD = 0;
		double testMAD = 0;
		int numOligos = oligosValues.length;
		int numValues = oligosValues[0].length;
		float[] one = new float[numOligos];
		//get first array
		for (int k=0; k<numOligos; k++){
			one[k]= oligosValues[k][0];
		}
		//run thru remainder
		for (int j=1; j<numValues; j++){
			//get second array
			float[] two = new float[numOligos];
			for (int k=0; k<numOligos; k++){
				two[k]= oligosValues[k][j];
			}
			//calc medianAbsDiff
			testMAD = medianAbsoluteDifference(one, two);
			if (testMAD>maxMAD) maxMAD = testMAD;
			//reset one for next cycle
			one = two;
		}
		return maxMAD;
	}

	/**Calculates the mean of the median absolute differences between arrays of float where float[oligo index number][oligo intensity measurements]*/
	public static double calcMeanMedianAbsoluteDifference(float[][] oligosValues){
		//calculate medianAbsDiff for each pair
		double totMAD = 0;
		int numOligos = oligosValues.length;
		int numValues = oligosValues[0].length;
		float[] one = new float[numOligos];
		double mads = 0;
		//get first array
		for (int k=0; k<numOligos; k++){
			one[k]= oligosValues[k][0];
		}
		//run thru remainder
		for (int j=1; j<numValues; j++){
			//get second array
			float[] two = new float[numOligos];
			for (int k=0; k<numOligos; k++){
				two[k]= oligosValues[k][j];
			}
			//calc medianAbsDiff
			totMAD += medianAbsoluteDifference(one, two);
			mads++;
			//reset one for next cycle
			one = two;
		}
		return totMAD/mads;
	}

	/**Calculates Pearson correlation coefficient, r, from two int[]s. 
	 * Cannot have one int[] be uniform values, returns -2 if error.*/
	public static double correlationCoefficient (int[] x, int[] y){
		double N = x.length;
		double xTot=0;
		double yTot=0;
		double sqrXTot =0;
		double sqrYTot =0;
		double xYTot=0;
		for (int i=0; i<N; i++){
			xTot += x[i];
			yTot += y[i];
			sqrXTot += Math.pow(x[i],2);
			sqrYTot += Math.pow(y[i],2);
			xYTot += (x[i] * y[i]);
		}
		double top = (N * xYTot) - (xTot * yTot);
		double botLeft = Math.sqrt( (N * sqrXTot) - Math.pow(xTot,2) );
		double botRight = Math.sqrt( (N * sqrYTot) - Math.pow(yTot,2) );
		double bot = botLeft*botRight;
		if (bot == 0) {
			String message = "Warning: calc of corr coef error, Num.java";
			System.err.println(message);
			return -2;
		}
		return top/bot;
	}
	
	/**Calculates Pearson correlation coefficient, r, from two int[]s. 
	 * Cannot have one double[] be uniform values, returns -2 if error.*/
	public static double correlationCoefficient (double[] x, double[] y){
		double N = x.length;
		double xTot=0;
		double yTot=0;
		double sqrXTot =0;
		double sqrYTot =0;
		double xYTot=0;
		for (int i=0; i<N; i++){
			xTot += x[i];
			yTot += y[i];
			sqrXTot += Math.pow(x[i],2);
			sqrYTot += Math.pow(y[i],2);
			xYTot += (x[i] * y[i]);
		}
		double top = (N * xYTot) - (xTot * yTot);
		double botLeft = Math.sqrt( (N * sqrXTot) - Math.pow(xTot,2) );
		double botRight = Math.sqrt( (N * sqrYTot) - Math.pow(yTot,2) );
		double bot = botLeft*botRight;
		if (bot == 0) {
			String message = "Warning: calc of corr coef error, Num.java";
			System.err.println(message);
			return -2;
		}
		return top/bot;
	}

	/**Provide a float[2][0], will remove zeros from both and there matched
	 * pair.
	 * Does a pairwise zero removal, if either of the float arrays are zero
	 * both values are removed.*/
	public static float[][] removeZeroValues(float[][] f){
		int num = f[0].length;
		ArrayList one = new ArrayList(num);
		ArrayList two = new ArrayList(num);
		for (int i=0; i<num; i++){
			if (f[0][i] != 0 && f[1][i] != 0) {
				one.add(new Float(f[0][i]));
				two.add(new Float(f[1][i]));
			}
		}
		float[][] str = new float[2][];
		str[0] = Num.arrayListOfFloatToArray(one);
		str[1] = Num.arrayListOfFloatToArray(two);
		return str;
	}



	/**Converts 0.345543 to 34.6% */
	public static String formatPercentOneFraction(double num){
		NumberFormat f = NumberFormat.getPercentInstance();
		f.setMaximumFractionDigits(1);
		return f.format(num);
	}

	/**Converts a double ddd.dddddddd to sss.s */
	public static String formatNumberOneFraction(double num){
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(1);
		return f.format(num);
	}
	/**Converts a double ddd.dddddddd to a user determined number of decimal places right of the .  */
	public static String formatNumber(double num, int numberOfDecimalPlaces){
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(numberOfDecimalPlaces);
		f.setMinimumFractionDigits(numberOfDecimalPlaces);
		return f.format(num);
	}
	
	/**Converts a double[] to int[] using Math.round()*/
	public static int[] doubleArrayToIntArray(double[] d){
		int[] ints = new int[d.length];
		for (int i=0; i< d.length; i++) ints[i] = (int)Math.round(d[i]);
		return ints;
	}
	
	/**Converts a double[] to float[] */
	public static float[] doubleArrayToFloatArray(double[] d){
		float[] floats = new float[d.length];
		for (int i=0; i< d.length; i++) floats[i] = (float)d[i];
		return floats;
	}
	
	/**Converts a double[] to int[] using Math.round()*/
	public static float[] intArrayToFloat(int[] ints){
		float[] floats = new float[ints.length];
		for (int i=0; i< ints.length; i++) floats[i] = ints[i];
		return floats;
	}

	/**Converts an array of double to a String with a defined number of decimal places and a delimiter.  */
	public static String doubleArrayToString(double[] d, int numberOfDecimalPlaces, String delimiter){
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(numberOfDecimalPlaces);
		f.setMinimumFractionDigits(numberOfDecimalPlaces);
		//return f.format();
		StringBuffer sb = new StringBuffer();
		int num = d.length;
		sb.append(f.format(d[0]));
		for (int i=1; i<num; i++){
			sb.append(delimiter);
			sb.append(f.format(d[i]));
		}
		return sb.toString();
	}

	/**Converts an array of double to a String with a maximal number of defined decimal places and a delimiter. Min isn't set.  */
	public static String doubleArrayToStringOnlyMax(double[] d, int maxNumDec, String delimiter){
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(maxNumDec);
		StringBuffer sb = new StringBuffer();
		int num = d.length;
		sb.append(f.format(d[0]));
		for (int i=1; i<num; i++){
			sb.append(delimiter);
			sb.append(f.format(d[i]));
		}
		return sb.toString();
	}

	/**Converts to String using delimiter.*/
	public static String floatArrayToString(float[] f, String delimiter){
		StringBuilder sb = new StringBuilder();
		sb.append(f[0]);
		for (int i=1; i< f.length; i++){
			sb.append(delimiter);
			sb.append(f[i]);
		}
		return sb.toString();
	}
	
	/**Converts to String using delimiter.*/
	public static String doubleArrayToString(double[] f, String delimiter){
		StringBuilder sb = new StringBuilder();
		sb.append(f[0]);
		for (int i=1; i< f.length; i++){
			sb.append(delimiter);
			sb.append(f[i]);
		}
		return sb.toString();
	}

	/**Converts an array of int to a String seperated by the delimiter.  */
	public static String intArrayToString(int[] d, String delimiter){
		StringBuffer sb = new StringBuffer();
		int num = d.length;
		sb.append(d[0]);
		for (int i=1; i<num; i++){
			sb.append(delimiter);
			sb.append(d[i]);
		}
		return sb.toString();
	}

	/**Given a String of ints delimited by something, will parse or return null.*/
	public static int[] stringArrayToInts(String s, String delimiter){
		String[] tokens = s.split(delimiter);
		int[] num = new int[tokens.length];
		try {
			for (int i=0; i< tokens.length; i++){
				num[i] = Integer.parseInt(tokens[i]);
			}
			return num;
		} catch (Exception e){
			return null;
		}
	}

	/**Given a String of floats delimited by something, will parse or return null.*/
	public static float[] stringArrayToFloat(String s, String delimiter){
		String[] tokens = s.split(delimiter);
		float[] num = new float[tokens.length];
		try {
			for (int i=0; i< tokens.length; i++){
				num[i] = Float.parseFloat(tokens[i]);
			}
			return num;
		} catch (Exception e){
			return null;
		}
	}

	/**Given a String of double delimited by something, will parse or return null.*/
	public static double[] stringArrayToDouble(String s, String delimiter){
		String[] tokens = s.split(delimiter);
		double[] num = new double[tokens.length];
		try {
			for (int i=0; i< tokens.length; i++){
				num[i] = Double.parseDouble(tokens[i]);
			}
			return num;
		} catch (Exception e){
			return null;
		}
	}

	/**Calculates Median of a sorted double[]. Copied code from WWW.*/
	public static double median(double[] sorted) {
		int middle = sorted.length/2;  // subscript of middle element
		//Odd number of elements -- return the middle one.
		if (sorted.length%2 == 1) return sorted[middle];
		// Even number -- return average of middle two
		return (sorted[middle-1] + sorted[middle]) / 2.0;
	}

	/**Calculates Median of a sorted int[].*/
	public static double median (int[] sorted) {
		int middle = sorted.length/2;  // subscript of middle element
		//Odd number of elements -- return the middle one.
		if (sorted.length%2 == 1) return sorted[middle];
		// Even number -- return average of middle two
		return ((double)(sorted[middle-1] + sorted[middle])) / 2.0; 
	}	


	/**Calculates Median of a sorted short[]. Copied code from WWW.*/
	public static double median(short[] sorted) {
		int middle = sorted.length/2;  // subscript of middle element
		//Odd number of elements -- return the middle one.
		if (sorted.length%2 == 1) return sorted[middle];
		// Even number -- return average of middle two
		return ((double)sorted[middle-1] + (double)sorted[middle]) / 2.0;
	}	

	/**Calculates Median of a sorted float[]. Copied code from WWW.*/
	public static double median(float[] sorted) {
		int middle = sorted.length/2;  // subscript of middle element
		//Odd number of elements -- return the middle one.
		if (sorted.length%2 == 1) return sorted[middle];
		// Even number -- return average of middle two
		return ((double)sorted[middle-1] + (double)sorted[middle]) / 2.0;
	}	
	/**Averages a float array.*/
	public static float mean(float[] t){
		float sumT=0;
		for (int i=t.length-1; i>=0; i--) sumT+= t[i];
		return sumT/(float)t.length;
	}
	/**Averages a float array, including startIndex, not including stopIndex.*/
	public static float mean(int startIndex, int stopIndex, float[] t){
		float sumT=0;
		for (int i=startIndex; i< stopIndex; i++) sumT+= t[i];
		return sumT/(float)(stopIndex-startIndex);
	}
	/**Averages two float[]s .*/
	public static float[] mean(float[] one, float[] two){
		int len = one.length;
		if (len != two.length) return null;
		float[] ave = new float[len];
		for (int i=0; i<len; i++) ave[i] = (one[i]+two[i])/2.0F;
		return ave;
	}
	/**Averages two float[][]s, assumes a square.*/
	public static float[][] mean(float[][] one, float[][] two){
		int len = one.length;
		if (len != two.length) return null;
		float[][] ave = new float[len][len];
		for (int i=0; i<len; i++){
			for (int j=0; j<len; j++){
				ave[i][j] = (one[i][j]+two[i][j])/2.0F;
			}
		}
		return ave;
	}

	/**Averages a int array.*/
	public static double mean(int[] t){
		long sumT=0;
		for (int i=t.length-1; i>=0; i--) sumT+= t[i];
		return (double)sumT/(double)t.length;
	}
	/**Averages a double array.*/
	public static double mean(double[] t){
		double sumT=0;
		for (int i=t.length-1; i>=0; i--) sumT+= t[i];
		return sumT/(double)t.length;
	}
	/**Averages an ArrayList of Integer objects.*/
	public static int meanIntegers(ArrayList Integers){
		int len = Integers.size();
		if (len==0) return 0;
		long total = 0;
		for (int i=0; i<len; i++){
			int num = ((Integer)Integers.get(i)).intValue();
			total += num;
		}
		double ave = (double)total/(double)len;
		return (int)Math.round(ave);
	}

	/**Averages an ArrayList of Double objects.*/
	public static String meanDoubles(ArrayList Doubles){
		int len = Doubles.size();
		if (len==0) return "0";
		double total = 0;
		for (int i=0; i<len; i++){
			double num = ((Double)Doubles.get(i)).doubleValue();
			total += num;
		}
		double ave = total/(double)len;
		NumberFormat formatter = NumberFormat.getNumberInstance();
		formatter.setMinimumFractionDigits(2);
		return formatter.format(ave);
	}	
	/**Finds and returns the biggest int in an int[].*/
	public static int findHighestInt(int[] ints) {
		int len = ints.length;
		int max = ints[0];
		for (int i = 1; i < len; i++) {
			if (ints[i]>max) max=ints[i];
		}
		return max;
	}
	/**Finds and returns the smallest int in an int[].*/
	public static int findSmallestInt(int[] ints) {
		int len = ints.length;
		int min = ints[0];
		for (int i = 1; i < len; i++) {
			if (ints[i]<min) min=ints[i];
		}
		return min;
	}
	/**Finds and returns the biggest float in an float[].*/
	public static float findHighestFloat(float[] ints) {
		int len = ints.length;
		float max = ints[0];
		for (int i = 1; i < len; i++) {
			if (ints[i]>max) max=ints[i];
		}
		return max;
	}
	/**Returns the sum at each position, assumes equal lenght.*/
	public static float[] sum(float[] one, float[] two){
		float[] sum = new float[one.length];
		for (int i=0; i< sum.length; i++) sum[i] = one[i]+two[i];
		return sum;
	}
	
	/**Sums an int array*/
	public static int sumIntArray(int[] ints){
		int num = 0;
		for (int i= ints.length-1; i>=0; i--) num+= ints[i];
		return num;
	}
	/**Sums a float array*/
	public static float sumArray(float[] array){
		float num = 0;
		for (int i= array.length-1; i>=0; i--) num+= array[i];
		return num;
	}
	/**Sums a float[][] array*/
	public static float sumArray(float[][] array){
		float num = 0;
		for (int i= array.length-1; i>=0; i--) {
			for (int j= array[i].length-1; j>=0; j--) num+= array[i][j];
		}
		return num;
	}
	/**Sums a double array*/
	public static double sumArray(double[] array){
		double num = 0;
		for (int i= array.length-1; i>=0; i--) num+= array[i];
		return num;
	}
	/**Sums a float array*/
	public static double sumArrayReturnDouble(float[] array){
		double num = 0;
		for (int i= array.length-1; i>=0; i--) num+= array[i];
		return num;
	}

	/**Finds mode of a int[] histogram array, assumes one peak, returns the index position, value/frequency.*/
	public static int[] modeOfHistogram(int[] ints) {
		int len = ints.length;
		int max = ints[0];
		int index = 0;
		for (int i = 1; i < len; i++) {
			if (ints[i]>max){
				max=ints[i];
				index = i;
			}
		}
		return new int[]{index,max};
	}

	/**Returns a reversed array.*/
	public static short[] reverse(short[] array) {
		int len = array.length;
		short[] rev = new short[len];
		int last = len -1;
		for (int i=0; i< len; i++){
			rev[last--] = array[i];
			//System.out.println("i "+i+"  last "+last);
			//System.out.println("\t"+array[i]+" "+array[last]);
			
		}
		return rev;
	}

	public static short[] arrayListOfShortToArray(ArrayList<Short> s) {
		int size = s.size();
		short[] c = new short[size];
		for (int i=0; i< size; i++) c[i] = s.get(i);
		return c;
	}

	public static double mean(short[] c) {
		double total = 0;
		for (short s: c) total += s;
		return total/((double)c.length);
	}
}
