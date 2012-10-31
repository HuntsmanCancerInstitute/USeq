package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.*;


/**Takes PointData and generates sliding window scan statistics. Estimates FDRs using two different methods.
 * @author Nix
 * */
public class CorrelatePointData {

	//user defined fields
	private File firstDir;
	private File secondDir;
	private File outputFile;

	//internal fields
	private HashMap<String, ArrayList<PointData>> firstPointData;
	private HashMap<String, ArrayList<PointData>> secondPointData;
	private String chromosome;
	private PointData firstChrom = null;	
	private PointData secondChrom = null;
	private PearsonCorrelation correlation = new PearsonCorrelation();
	private long numberOfPairs = 0;
	private PrintWriter out = null;
	private boolean printPairs = false;

	//constructors

	/**Stand alone.*/
	public CorrelatePointData(String[] args){
		long startTime = System.currentTimeMillis();
		try {

			//set fields
			processArgs(args);

			//make print writer?
			if (outputFile != null) out = new PrintWriter (new FileWriter(outputFile));
				
			//for each chromosome
			System.out.println("Scanning Chromosomes:");
			String[] chromosomes = fetchAllChromosomes();
			for (int i=0; i< chromosomes.length; i++){
				chromosome = chromosomes[i];
				scanChromosome();
			}
			System.out.println();

			//print correlation
			System.out.println("Number of paired observations : "+numberOfPairs);
			if (numberOfPairs > 10) {
				double r = correlation.calculateAdditivePairCorrelation();
				System.out.println("\nPearson correlation (r) : "+r);
				System.out.println("            (r^2 * 100) : "+100*r*r);
			}
			else System.err.println("\nToo few paired observations, aborting!");

			//close writer
			if (printPairs) out.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods
	/**Fetches the names of all the chromosomes in the data.*/
	public String[] fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		c.addAll(firstPointData.keySet());
		c.addAll(secondPointData.keySet());
		return Misc.hashSetToStringArray(c);
	}

	/**Fetchs the data for a particular chromosome.*/
	public boolean fetchData(){
		//fetch treatment
		firstChrom = null;
		if (firstPointData.containsKey(chromosome)) firstChrom = firstPointData.get(chromosome).get(0);
		//fetch control
		secondChrom = null;
		if (secondPointData.containsKey(chromosome)) secondChrom = secondPointData.get(chromosome).get(0);
		if (firstChrom == null || secondChrom == null) return false;
		return true;
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void scanChromosome(){
		//fetch data
		if (fetchData() == false) {
			System.out.println("\t"+chromosome+" Skipping!! Failed to find PointData from both datasets.");
			return;
		}


		//fetch data
		int[] firstPositions = firstChrom.getPositions();
		float[] firstScores = firstChrom.getScores();
		int[] secondPositions = secondChrom.getPositions();
		float[] secondScores = secondChrom.getScores();
		int indexSecond = 0;

		System.out.println("\t"+chromosome+"\t"+firstPositions.length+"\t"+secondPositions.length);

		for (int firstIndex=0; firstIndex< firstPositions.length; firstIndex++){
			for (;indexSecond< secondPositions.length; indexSecond++){
				//is first < second, thus second has advanced past first
				if (firstPositions[firstIndex] < secondPositions[indexSecond]) break;
				//are they the same, thus a match
				if (firstPositions[firstIndex] == secondPositions[indexSecond]) {
					correlation.addPairsToCorrelate(firstScores[firstIndex], secondScores[indexSecond]);
					//print scores to file for scatter plot?
					if (printPairs) {
						out.print(firstScores[firstIndex]);
						out.print("\t");
						out.println(secondScores[indexSecond]);
					}
					indexSecond++;
					numberOfPairs++;
					break;
				}
				//must be that second is < first so do nothing and let it advance
			}
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CorrelatePointData(args);
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
					case 'f': firstDir = new File (args[++i]); break;
					case 's': secondDir = new File (args[++i]); break;
					case 'p': outputFile = new File (args[++i]); printPairs = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (firstDir == null || firstDir.isDirectory() == false || secondDir == null || secondDir.isDirectory() == false) Misc.printErrAndExit("\nPlease enter two PointData directories.\n");

		//fetch PointData
		firstPointData = PointData.fetchPointData (firstDir);
		secondPointData = PointData.fetchPointData (secondDir);
		if (firstPointData == null || firstPointData.size() ==0 || secondPointData == null || secondPointData.size() ==0 ) Misc.printErrAndExit("\nCannot find any PointData in one or both of your directories?\n");


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                CorrelatePointData: Aug 2011                      **\n" +
				"**************************************************************************************\n" +
				"Calculates a Pearson Correlation Coefficient on the values of PointData found with the\n" +
				"same positions in the two datasets. Do NOT use on stair-step/ heat-map graph data.\n" +
				"Only use on point representation data.\n\n" +

				"Options:\n"+
				"-f First PointData set. This directory should contain chromosome specific xxx.bar.zip\n" +
				"       files, stranded or unstranded.\n"+
				"-s Second PointData set, ditto. \n" +
				"-p Full path file name to use in saving paired scores, defaults to not printing.\n"+


				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/CorrelatePointData -f /BaseFracMethyl/X1\n" +
				"      -s /BaseFracMethyl/X2 \n\n" +

		"**************************************************************************************\n");

	}


}
