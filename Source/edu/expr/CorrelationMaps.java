package edu.expr;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;
import util.bio.annotation.*;
import trans.tpmap.*;
import trans.main.*;

public class CorrelationMaps {
	//fields
	private File geneFile;
	private File regionFilterFile;
	private ExpressedGene[] expressedGenes;
	private ExpressedGene[][] chromExprGenes;
	private SpearmanCorrelation spCorr;
	private int windowSize = 50000;
	private int minNumberGenes = 3;
	private int numberRandomTrials = 100;
	private Histogram histogram;
	private WindowMaker windowMaker;
	private CorrelationWindow[][] correlationWindows;
	public static final float[] SQUARE_ROOTS = Num.squareRoots(1000);
	private String genomeName = "C_elegans_May_2007";
	private boolean multiValuesPresent = false;
	private boolean multiplyScoreBySquareRoot = false;

	//constructors
	public CorrelationMaps(String[] args){
		processArgs(args);
		if (multiplyScoreBySquareRoot) System.out.println("Multiplying scores by square root of the number of genes...");
		//parse gene file
		expressedGenes = parseGeneFile(geneFile);
		System.out.println("\nNumber genes "+expressedGenes.length);
		//filter genes
		if (regionFilterFile != null) {
			filterGenes();
			System.out.println("Number filtered genes "+expressedGenes.length);
		}
		//more than one value?
		if (expressedGenes[0].getValues().length > 1 ) multiValuesPresent = true;
		else System.out.println("\tOnly one value found, taking mean...");
		//sort by chromosome, start, length
		Arrays.sort(expressedGenes);
		//split by chromosome
		chromExprGenes = ExpressedGene.splitByChromosome(expressedGenes);
		//correlate windows of genes
		correlate();
		//assign empirical p-values?
		if (numberRandomTrials > 0) assignPValues();
		//print spreadsheet
		printCorrelationWindows();
		//print bar files for p-values and actual correlations
		printBarFiles();
		System.out.println("\nDone!\n");
	}
	
	
	//methods
	
	public void filterGenes(){
		Coordinate[] operons = Coordinate.parseFile(regionFilterFile, 0, 0);
		if (operons == null) Misc.printExit("\nCannot parse region filter file. "+regionFilterFile +"\n");
		ArrayList cleanGenes = new ArrayList();
		for (int i=0; i< expressedGenes.length; i++){
			boolean noIntersection = true;
			Coordinate expressedGeneCoordinate = expressedGenes[i].getCoordinates();
			for (int j=0; j< operons.length; j++){
				if (operons[j].intersects(expressedGeneCoordinate)){
					//System.out.println(expressedGeneCoordinate);
					noIntersection = false;
					break;
				}
			}
			if (noIntersection) cleanGenes.add(expressedGenes[i]);
		}
		expressedGenes = new ExpressedGene[cleanGenes.size()];
		cleanGenes.toArray(expressedGenes);
	}
	
	public void printBarFiles(){
		//make Windows
		ArrayList al = new ArrayList();
		//for each CorrelationWindow, only save those with a positive score
		for (int i=0; i<correlationWindows.length; i++){
			for (int j=0; j< correlationWindows[i].length; j++){
				if (correlationWindows[i][j].getScore() > 0){
					Window win = new Window(correlationWindows[i][j].getChromosome(), correlationWindows[i][j].getStart(), 
							correlationWindows[i][j].getStop(), correlationWindows[i][j].getNames().length, 
							new double[]{correlationWindows[i][j].getScore(), correlationWindows[i][j].getPvalue()});
					al.add(win);
				}
			}
		}
		//convert to array
		Window[] windows = new Window[al.size()];
		al.toArray(windows);
		//make folder to hold bars
		String trunc = Misc.removeExtension(geneFile.getName());
		File scores = new File (geneFile.getParentFile(), trunc+"_MeanHM");
		if (scores.exists()==false) scores.mkdir();
		WindowBlockMaker blocks = new WindowBlockMaker(2, genomeName);
		blocks.makeHeatMapBarFiles(windows, scores);
		//write for pvalues?
		if (numberRandomTrials > 0) {
			blocks.setScoreIndex(1);
			File pvals = new File (geneFile.getParentFile(), trunc+"_PValHM");
			if (pvals.exists()==false) pvals.mkdir();
			blocks.makeHeatMapBarFiles(windows, pvals);
		}
	}
	/**Prints a spread sheet summary and a xxx.egr file for all of the windows.*/
	public void printCorrelationWindows(){
		try {
			File file = new File (geneFile.getParentFile(), Misc.removeExtension(geneFile.getName()) +"_Report.xls");
			File egrFile = new File (geneFile.getParentFile(), Misc.removeExtension(geneFile.getName()) +"_Win.egr");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			PrintWriter egrOut = new PrintWriter (new FileWriter (egrFile));
			out.println("#Names\tChromosome\tStart\tStop\tNumber\tMeanWindowScore\t-10Log10(pvalue)");
			//for each CorrelationWindow
			for (int i=0; i<correlationWindows.length; i++){
				for (int j=0; j< correlationWindows[i].length; j++){
					out.println(correlationWindows[i][j]);
					egrOut.println(correlationWindows[i][j].getChromosome()+"\t"+
							correlationWindows[i][j].getStart()+"\t"+
							correlationWindows[i][j].getStop()+"\t.\t"+
							correlationWindows[i][j].getScore());
				}
			}
			egrOut.close();
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Assigns -10Log10(pvalues) for each score > 0 in the CorrelationWindows*/
	public void assignPValues(){
		//check histogram for overflow
		//histogram.printScaledHistogram();
		//for each CorrelationWindow
		for (int i=0; i<correlationWindows.length; i++){
			for (int j=0; j< correlationWindows[i].length; j++){
				CorrelationWindow win = correlationWindows[i][j];
				if (win.getScore() > 0){
					double pval = histogram.pValueRightSide(win.getScore());
					//one tailed
					if (multiValuesPresent) pval = histogram.pValueRightSide(win.getScore());
					//two tailed
					else pval = histogram.pValue(win.getScore());
					pval = -10 * Num.log10(pval);
					win.setPvalue(new Double(pval).floatValue());
				}
			}
		}
	}
	
	public float[] minMax() {
		float min = 1000000000;
		float max = -1000000000;
		for (int i=0; i< expressedGenes.length; i++){
			float val = expressedGenes[i].getValues()[0];
			if (val > max) max = val;
			if (val < min) min = val;
		}
		
		/*Random randomGenerator = new Random();
		for (int i=0; i< expressedGenes.length; i++){
		int next = randomGenerator.nextInt(expressedGenes.length);
		float[] ori = expressedGenes[i].getValues();
		float rand[] = expressedGenes[next].getValues();
		expressedGenes[i].setValues(rand);
		expressedGenes[next].setValues(ori);
		}*/
		return new float[] {min, max};
	}
	
	/**Creates CorrelationWindows and performs random trials.
	 * Nulls various unused objects.*/
	public void correlate(){
		//instatiate objects
		spCorr = new SpearmanCorrelation();
		windowMaker = new WindowMaker(windowSize, minNumberGenes);
		if (numberRandomTrials > 0) {
			if (multiValuesPresent) histogram = new Histogram(-4,4,5000);
			else {
				//must find min and max
				System.out.print("\tFinding min and max... ");
				float[] minMax = minMax();
				System.out.println(minMax[0]+" "+minMax[1]);
				//multiply max by 2
				minMax[1] = minMax[1]*2;
				histogram = new Histogram(minMax[0],minMax[1],5000);
			}
		}
		correlationWindows = new CorrelationWindow[chromExprGenes.length][];

		//for each chromosome
		for (int i=0; i< chromExprGenes.length; i++){
			//make RankedFloatArrays? This call will null the values
			if (multiValuesPresent) chromExprGenes[i] = makeRankedFloatArrays(chromExprGenes[i]);
			String chrom = chromExprGenes[i][0].getCoordinates().getChromosome();
			System.out.println("Correlating "+chrom);
			//split ExpressedGene[] into windows
			ExpressedGene[][] exprGenWin = window(chromExprGenes[i]);
			correlationWindows[i] = new CorrelationWindow[exprGenWin.length];
			System.out.println("\t"+exprGenWin.length+" windows");
			//for each window, make CorrelationWindows
			for (int j=0; j< exprGenWin.length; j++){
				//calculate corr sum?
				float sum;
				if (multiValuesPresent) sum = pairwiseCorrelationScore (exprGenWin[j]);
				else sum = averageScores (exprGenWin[j]);
				//start
				int start = exprGenWin[j][0].getCoordinates().getStart();
				//stop
				int end = exprGenWin[j][exprGenWin[j].length-1].getCoordinates().getStop();
				//names
				String[] names = new String[exprGenWin[j].length];
				for (int k=0; k< exprGenWin[j].length; k++)	names[k] = exprGenWin[j][k].getName();
				//make window
				correlationWindows[i][j] = new CorrelationWindow(sum, chrom,start, end, names);
			}

			//calculate randomized scores and add to histogram
			if (numberRandomTrials > 0){
				System.out.print("\tRandom trials ");
				for (int r=0; r< numberRandomTrials; r++){
					System.out.print(".");
					//randomize chrom specific ExpressedGene.RankedFloatArray or float[] values, this, by reference, randomizes the ExpressedGene[][] exprGenWin too  
					randomize(i);
					//for each window
					for (int j=0; j< exprGenWin.length; j++){
						float sum;
						if (multiValuesPresent) sum = pairwiseCorrelationScore (exprGenWin[j]);
						else sum = averageScores (exprGenWin[j]);
						//add to histogram
						histogram.count(sum);
					}
				}
				System.out.println();
			}
			chromExprGenes[i] = null;
		}
		//null unused
		expressedGenes = null;
		chromExprGenes = null;
	}

	public void randomize(int index){
		Random randomGenerator = new Random();
		int length = chromExprGenes[index].length;
		if (multiValuesPresent){
			for (int i=0; i< length; i++){
				int next = randomGenerator.nextInt(length);
				RankedFloatArray ori = chromExprGenes[index][i].getRfa();
				RankedFloatArray rand = chromExprGenes[index][next].getRfa();
				chromExprGenes[index][i].setRfa(rand);
				chromExprGenes[index][next].setRfa(ori);
			}
		}
		else {
			for (int i=0; i< length; i++){
				int next = randomGenerator.nextInt(length);
				float[] ori = chromExprGenes[index][i].getValues();
				float[] rand = chromExprGenes[index][next].getValues();
				chromExprGenes[index][i].setValues(rand);
				chromExprGenes[index][next].setValues(ori);
			}
		}

	}

	/**Average coorelation coefficient * square root of number of genes.*/
	public float pairwiseCorrelationScore (ExpressedGene[] ex){
		float sum = 0; 
		int len = ex.length;
		float numPairs = (len*(len-1))/2;

		//calc all pair Spearman Rank Correlations
		for (int i=0; i< len; i++){
			for (int j=i+1; j< len; j++){
				//System.out.println("\t"+ex[i].getName()+" "+ ex[j].getName()+" : "+ spCorr.spearmanCorrelationCoefficient(ex[i].getRfa(), ex[j].getRfa()));
				sum += spCorr.spearmanCorrelationCoefficient(ex[i].getRfa(), ex[j].getRfa());
			}
		}
		if (multiplyScoreBySquareRoot) return SQUARE_ROOTS[len] * (sum/numPairs);
		else return sum/numPairs;
	}
	
	/**Average scores * square root of number of genes.*/
	public float averageScores (ExpressedGene[] ex){
		float sum = 0; 
		float len = ex.length;
		//average values across window
		for (int i=0; i< len; i++){
			sum += ex[i].getValues()[0];
		}
		if (multiplyScoreBySquareRoot) return SQUARE_ROOTS[(int)len] * (sum/len);
		else return (sum/len);
	}

	/**Splits a chromosome of ExpressedGenes into windows based on the windowSize and minNumberGenes*/
	public ExpressedGene[][] window(ExpressedGene[] ex){
		//find middles of each gene
		int[] middles = new int[ex.length];
		for (int i=0; i< ex.length; i++) {
			Coordinate coor = ex[i].getCoordinates();
			double len = (coor.getStop() - coor.getStart() +1);
			middles[i] = (int)Math.round(len/2)+ coor.getStart();
		}
		//make windows
		int[][] indexes = windowMaker.makeWindows(middles);
		//make final
		ExpressedGene[][] exSplit = new ExpressedGene[indexes.length][];
		for (int i=0; i< indexes.length; i++){
			int size = indexes[i][1] - indexes[i][0]+1;
			exSplit[i] = new ExpressedGene[size];
			int counter =0;
			for (int j=indexes[i][0]; j<= indexes[i][1]; j++){
				exSplit[i][counter++] = ex[j];
			}
		}
		return exSplit;
	}

	/**Makes and places RankedFloatArrays in each ExpressedGene,
	 * also nulls the ExpressedGene.values.*/
	public ExpressedGene[] makeRankedFloatArrays(ExpressedGene[] ex){
		for (int i=0; i< ex.length; i++){
			RankedFloatArray r = spCorr.rank(ex[i].getValues());
			ex[i].setRfa(r);
			ex[i].setValues(null);
		}
		return ex;
	}

	/**Converts a tab delimited file of (text chr start stop val1 val2 ...) into an ExpressedGene[]*/
	public static ExpressedGene[] parseGeneFile(File geneFile){
		try{
			ArrayList al = new ArrayList();
			BufferedReader in = new BufferedReader (new FileReader (geneFile));
			String line;
			String[] tokens;
			while ((line=in.readLine()) != null){
				if (line.startsWith("#")) continue;
				tokens = line.split("\\t");
				if (tokens.length < 4) continue;
				String name = tokens[0];
				String chr = tokens[1];
				int start = Integer.parseInt(tokens[2]);
				int stop = Integer.parseInt(tokens[3]);
				Coordinate coor = new Coordinate(chr, start, stop);
				//remainder are values
				int num = tokens.length - 4;
				float[] values = new float[num];
				int index = 0;
				for (int i=4; i< tokens.length; i++) values[index++] = Float.parseFloat(tokens[i]);
				//make obj
				al.add(new ExpressedGene(name, coor, values));
			}
			ExpressedGene[] expressedGenes = new ExpressedGene[al.size()];
			al.toArray(expressedGenes);
			in.close();
			return expressedGenes;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CorrelationMaps(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': geneFile = new File (args[i+1]); i++; break;
					case 'o': regionFilterFile = new File (args[i+1]); i++; break;
					case 'w': windowSize = Integer.parseInt(args[i+1]); i++; break;
					case 'n': minNumberGenes  = Integer.parseInt(args[i+1]); i++; break;
					case 'r': numberRandomTrials  = Integer.parseInt(args[i+1]); i++; break;
					case 'g': genomeName = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (geneFile == null || geneFile.canRead() == false) Misc.printExit("\nError: cannot find your gene file! "+geneFile);
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Correlation Maps:    Nov 2007                         **\n" +
				"**************************************************************************************\n" +
				"CM calculates a correlation score for each window of genes and using permutation, an\n" +
				"empirical p-value.  The correlation score is the mean of all pair Spearman ranks for\n" +
				"the gene expression profiles in each window. If a single value is given (unlogged!) for\n" +
				"each gene, a mean of the scores within each window is calculated.\n\n" +
				
				"To calculate p-values, X randomized datasets are created by shuffling the expression\n" +
				"profiles between genes, windows are scored and pooled.  P-values for each real\n" +
				"score are calculated based on the area under the right side of the randomized score\n" +
				"distribution. In addition to a spread sheet report summary, heat map xxx.bar files\n" +
				"for the p-values and mean correlation are created for visualization in IGB.\n" +
				"Note, this analysis is not stranded.  If so desired parse lists appropriately.\n\n" +
				
				"Parameters:\n" +
				"-f The full path file text for a tab delimited gene file (text,chr,start,stop,scores)\n" +
				"-o GenomicRegion filter file, full path file text for a tab delimited region file to use in\n"+
				"      removing genes from correlation analysis. (chrom, start, stop).\n"+
				"-g Genome version for IGB visualizations (e.g. C_elegans_May_2007).\n"+
				"-w Window size, default is 50000bp. Setting this too small may exclude some regions.\n" +
				"-n Minimum number of genes required in each window, defaults to 3. Setting this too\n" +
				"       high will exclude some regions.\n" +
				"-r Number random trials, defaults to 100\n\n" +
				
				"Example: java -Xmx256M -jar pathTo/T2/Apps/CorrelationMaps -f /Mango/geneFile.txt\n" +
				"       -w 30000 -n 2 -o /Mango/operons.txt\n\n" +
				
		"**************************************************************************************\n");		
	}	

}
