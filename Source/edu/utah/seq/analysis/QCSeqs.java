package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.bio.seq.*;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.main.Window;
import trans.tpmap.*;

/**Takes PointData and generates sliding window average scores.*/
public class QCSeqs {

	//fields
	private File[] pointDataDirectories;
	private String[] trunkatedNames;
	private File tempDirectory;
	//chromosomes found 
	private HashMap<String, Integer> chromosomes;
	//window scores for merged strand PointData
	private HashMap<String,ChromWindowScores>[] windowScores = null;
	//defaults
	private int windowSize = 500;
	private int stepSize = 250;
	private float minWindowScore = 5;
	private boolean sumNumberPositions = true;
	private File windowFile;



	//constructors
	public QCSeqs(String[] args){

		//set fields
		processArgs(args);

		//for each point data directory, merges stranded data and saves it to the temporary directory
		mergeAndSavePointData();

		//advance a window across each chromosome and save the window scores with any hits
		windowScorePointData();

		//load windowScores with missing chromosomes
		loadMissingChromosomes();

		//perform an all pair Pearson correlation
		correlate();

		//delete temp dir containing merged bar files
		IO.deleteDirectory(tempDirectory);

	}

	public void correlate(){
		System.out.println("\nAll pair correlations:\n\t\tr\t100*r^2\tFirst\tSecond");
		for (int i=0; i< windowScores.length; i++){
			//get text
			Iterator<String> x = windowScores[i].keySet().iterator();
			String firstName = windowScores[i].get(x.next()).name;
			for (int j= i+1; j< windowScores.length; j++){
				//get text
				x = windowScores[j].keySet().iterator();
				String secondName = windowScores[j].get(x.next()).name;
				//get r value
				double r = correlate (windowScores[i], windowScores[j]);
				double per = r * r ;
				String rS = Num.formatNumber(r, 3);
				String perS = Num.formatPercentOneFraction(per);
				System.out.println("\t\t"+rS+"\t"+perS+"\t"+firstName+"\t"+secondName);
			}
			System.out.println();
		}
	}

	public double correlate(HashMap<String,ChromWindowScores> first, HashMap<String,ChromWindowScores> second){
		try {
			PrintWriter out = null;
			if (windowFile != null){
				out = new PrintWriter (new FileWriter(windowFile));
			}
			PearsonCorrelation pc = new PearsonCorrelation();
			//for each chromosome
			Iterator<String> it = first.keySet().iterator();
			while (it.hasNext()){
				String chrom = it.next();
				float[] one = first.get(chrom).getWindowScores();
				float[] two = second.get(chrom).getWindowScores();
				pc.addPairsToCorrelate(one, two);
				if (out !=null){
					for (int i=0; i< one.length; i++) out.println(one[i]+"\t"+two[i]);
				}
			}
			if (out !=null) out.close();
			
			return pc.calculateAdditivePairCorrelation();
		}
		catch (Exception e){
			e.printStackTrace();
		}
		return 0;
	}

	public void loadMissingChromosomes(){
		System.out.println("\nLooking for and loading missing chromosomes...");
		//for each chromosome
		Iterator<String> it = chromosomes.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();

			//create an array to use as replacement, WARNING this will be shared by all ChromWindowScores
			int numWin = chromosomes.get(chrom).intValue();
			float[] f = new float[numWin];

			//for each windowScores
			for (int i=0; i< windowScores.length; i++){

				if (windowScores[i].containsKey(chrom) == false){
					//get text
					String name = null;
					Iterator<String> x = windowScores[i].keySet().iterator();
					name = windowScores[i].get(x.next()).name;

					//instantiate ChromWindowScores, add text and shared scores
					ChromWindowScores cw = new ChromWindowScores();
					cw.name = name;
					cw.windowScores = f;

					//add to 
					windowScores[i].put(chrom, cw);
					System.out.println("\tAdding blank "+chrom+" to "+name);
				}
			}
		}
	}

	/**Calculates window averages for merged PointData.*/
	public void windowScorePointData(){		
		System.out.println("Scoring windows...");

		//for each chromosome
		Iterator<String> it = chromosomes.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			System.out.println("\t"+ chrom);

			//calc max position
			int max = maxPosition (chrom);
			System.out.println("\t\t"+max+"\tMax bp");

			//fetch array of ChromWindowScores, chrom specific, possibly one per directory
			ChromWindowScores[] ws = fetchChromSpecificWindowScores(chrom);

			//make windows of set size and collect any hits underneath
			int numWindows = 0;
			for (int i=0; i<max; i+=stepSize){
				//define start and stop bp of window
				int start = i;
				int stop = i+windowSize;

				//define scores for each ChromWindowScores
				float[] scores = new float[ws.length];

				//for each ChromWindowScores, chrom specific, find if any hits
				boolean hitFound = false;
				for (int x=0; x< ws.length; x++){
					if (sumNumberPositions) scores[x] = ws[x].pointData.sumPositionBP(start, stop);
					else scores[x] = ws[x].pointData.sumScoreBP(start, stop);
					if (scores[x] >= minWindowScore) hitFound = true;
				}
				//saveScores if any hits
				if (hitFound){
					for (int x=0; x< ws.length; x++){
						ws[x].windowScoresAL.add(new Float(scores[x]));
					}
					numWindows++;					
				}
			}
			//convert ALs to float[] and null position[] and score[] by calling getWindowScores()
			for (int x=0; x< ws.length; x++) ws[x].getWindowScores();
			System.out.println("\t\t"+ numWindows +"\tTotal windows");
			//add num window
			chromosomes.put(chrom, new Integer(numWindows));
		}
	}

	private class ChromWindowScores{
		//fields
		PointData pointData;
		private float[] windowScores;
		ArrayList<Float> windowScoresAL = new ArrayList<Float>();
		String name;

		//constructor
		public ChromWindowScores (PointData pointData,String name){
			this.name= name;
			this.pointData = pointData;
		}
		public ChromWindowScores(){}

		public float[] getWindowScores(){
			if (windowScores == null){
				windowScores = Num.arrayListOfFloatToArray(windowScoresAL);
				windowScoresAL = null;
				pointData.nullPositionScoreArrays();
			}
			return windowScores;
		}
	}



	public ChromWindowScores[] fetchChromSpecificWindowScores(String chrom){
		ArrayList<ChromWindowScores> al = new ArrayList();
		for (int i=0; i<windowScores.length; i++){
			ChromWindowScores cws = windowScores[i].get(chrom);
			if (cws != null) al.add(cws);
		}
		ChromWindowScores[] toReturn = new ChromWindowScores[al.size()];
		al.toArray(toReturn);
		return toReturn;
	}

	public int maxPosition (String chrom){
		int max =0;
		for (int i=0; i<windowScores.length; i++){
			ChromWindowScores cws = windowScores[i].get(chrom);
			if (cws == null){
				System.out.println("\t\tMissing from "+pointDataDirectories[i].getName());
				continue;
			}
			PointData pd = cws.pointData;
			int[] positions = pd.getPositions();
			int val = positions[positions.length-1];
			if (val > max) max = val;
		}
		return max;
	}


	public void mergeAndSavePointData(){
		System.out.println("Merging and Saving Stranded Point Data...");

		//instantiate windowScores, one for every pointDataDir
		windowScores = new HashMap[pointDataDirectories.length];
		for (int i=0; i< pointDataDirectories.length; i++) windowScores[i] = new HashMap();

		//merged strand PointData stored in temp dir
		HashMap<String,PointData>[] pointData = new HashMap[pointDataDirectories.length];

		//instantiate Hash to hold all chromsomes
		chromosomes = new HashMap();

		//for each point data directory
		for (int i=0; i< pointDataDirectories.length; i++){
			//make random word to use in naming temp directories
			String randomWord = Passwords.createRandowWord(6);

			//make temp dir
			File td = new File (tempDirectory, pointDataDirectories[i].getName()+"_"+i+"_"+randomWord);
			td.mkdir();
			System.out.println("\t "+pointDataDirectories[i]);

			//merge point data and save in temp directory
			pointData[i] = PointData.fetchPointData(pointDataDirectories[i], td,true);

			//add chromosomes
			Iterator<String> it1 = pointData[i].keySet().iterator();
			while (it1.hasNext()) chromosomes.put(it1.next(), new Integer(0));

			//instantiate windowScores and load with merged point data
			Iterator it = pointData[i].keySet().iterator();
			while (it.hasNext()){
				String chrom = (String)it.next();
				windowScores[i].put(chrom, new ChromWindowScores(pointData[i].get(chrom), trunkatedNames[i]));
			}
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new QCSeqs(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		String dirs = null;
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': dirs = args[i+1]; i++; break;
					case 't': tempDirectory = new File(args[i+1]); i++; break;
					case 'e': windowFile = new File(args[i+1]); i++; break;
					case 'w': windowSize = Integer.parseInt(args[i+1]); i++; break;
					case 's': stepSize = Integer.parseInt(args[i+1]); i++; break;
					case 'm': minWindowScore = Float.parseFloat(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for point directories, this sorts the dirs
		pointDataDirectories = IO.extractFiles(dirs);
		if (pointDataDirectories == null || pointDataDirectories[0].isDirectory() == false) Misc.printExit("\nError: cannot find your Point Data directories(s)!\n");

		//make uniqueNames
		String[] fullPathNames = new String[pointDataDirectories.length];
		for (int y=0; y< pointDataDirectories.length; y++) fullPathNames[y] = pointDataDirectories[y].toString();
		trunkatedNames = Misc.trimCommon(fullPathNames);

		//look for and or create the save directory
		if (tempDirectory == null) Misc.printExit("\nError: enter a directory text for saving temporary results.\n");
		if (tempDirectory.exists() == false) tempDirectory.mkdir();
		else Misc.printExit("Error: the temp directory exists, select another.\n");

		//print params
		System.out.println("Window size: "+windowSize+"  window step size: "+stepSize+"  minimum # reads within window: "+minWindowScore+"\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 QCSeqs: Nov 2009                                 **\n" +
				"**************************************************************************************\n" +
				"QCSeqs takes directories of chromosome specific PointData xxx.bar.zip files that \n" +
				"represent replicas of signature sequencing data, merges the strands, uses a sliding\n" +
				"window to sum the hits, and calculate Pearson correlation coefficients for the window\n" +
				"sums between each pair of replicas.  Only windows with a sum score >= the minimum \n" +
				"are included in the correlation.\n\n" +

				"-d Split chromosome Point Data directories, full path, comma delimited. (These should\n" +
				"       contain chromosome specific xxx.bar.zip files). \n" +
				"-t Temp directory, full path. This will be created and then deleted.\n"+
				"-w Window size in bps, defaults to 500.\n"+
				"-s Window step size in bps, defaults to 250.\n"+
				"-m Minimum window sum score, defaults to 5.\n"+
				"-e (Optional) Provide a full path file name in which to write the window sums.\n"+

				"\nExample: java -Xmx1500M -jar pathTo/USeqs/Apps/QCSeqs -d /Solexa/PolII/Rep1PntData/,\n" +
				"      /Solexa/PolII/Rep2PntData/ -t /Solexa/PolII/TempDelMe -w 1000 -s 250 \n\n" +


		"**************************************************************************************\n");

	}
}
