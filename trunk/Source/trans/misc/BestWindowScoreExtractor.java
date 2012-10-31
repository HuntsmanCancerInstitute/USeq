package trans.misc;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import trans.main.Window;
import util.gen.*;

/**Finds the best overlapping Window for each region and prints that score with the region information.
 * Built to associate a qPCR test result with the most overlapping window.*/
public class BestWindowScoreExtractor {

	//fields
	private Window[] testWindows;
	private QPCRWindow[] qPCRWindows;
	private int scoreIndex = 1;
	private int sizeOfOligoMinusOne= 24;
	private boolean multiplyByNegativeOne = false;
	private boolean findBestScoringWindow = false;

	/**For each validation region, finds the single test window that 'best' overlaps.*/
	public BestWindowScoreExtractor(String[] args){
		processArgs(args);

		//add length of oligo on to ends of testWindows
		int numTestWin = testWindows.length;
		for (int i=0; i<numTestWin; i++){
			testWindows[i].setStartLastOligo( testWindows[i].getStartLastOligo() + sizeOfOligoMinusOne );
			if (multiplyByNegativeOne) {
				double[] scores = testWindows[i].getScores();
				scores[scoreIndex] = -1 * scores[scoreIndex];
			}
		}

		if (multiplyByNegativeOne) System.out.println("\tMultiplying scores by -1");

		System.out.println("Region_Chr\tRegion_Start\tRegion_Stop\tRegion_BestWindowScore\tWindow_Chr\tWindow_Start\tWindow_Stop\tWindow_Scores");
		
		if (findBestScoringWindow){
			//for each qPCR window find overlapping testWindows
			int numQPCRWin = qPCRWindows.length;			
			for (int i=0; i< numQPCRWin; i++){
				Window qPCR = qPCRWindows[i].window;

				//identify overlapping test window closest to center of qPCR region.
				Window bestWindow = null;
				double bestScore = -100000000;
				boolean overlapped = false;
				for (int j=0; j< numTestWin; j++){
					Window test = testWindows[j];
					if (qPCR.overlap(test)){
						double testScore = test.getScores()[scoreIndex];
						overlapped = true;
						if (testScore> bestScore) {
							bestWindow = test;
							bestScore = testScore;
						}
					}
					else if (overlapped) break;
				}
				//print out qPCR window info with maxScore
				String finalScore;
				if (bestWindow == null) finalScore = "No Overlapping Windows?";
				else finalScore = Num.formatNumber(bestScore, 4);
				//text rep qPCR region + it's region scores + best window text rep and its	
				StringBuffer sb = new StringBuffer();
				sb.append(qPCR.stringRep(0));
				if (qPCR.getScores() != null) sb.append("\t"+ Num.doubleArrayToStringOnlyMax(qPCR.getScores(), 4, "\t"));
				sb.append("\t");
				sb.append(finalScore);
				//best window found? print info
				if (bestWindow != null){
					sb.append("\t");
					sb.append(bestWindow.stringRep(0));
					sb.append("\t");
					sb.append(Num.doubleArrayToStringOnlyMax(bestWindow.getScores(), 4, "\t"));
				}
				System.out.println(sb);
			}
		}
		else {
			//this is sort of inefficient, lazy bastard!
			//for each qPCR window find overlapping testWindows
			int numQPCRWin = qPCRWindows.length;
			for (int i=0; i< numQPCRWin; i++){
				Window qPCR = qPCRWindows[i].window;
				int qPCRLength = qPCR.getStartLastOligo()- qPCR.getStart1stOligo()+1;
				double minMidPointDiff = 1000000;
				double maxScore = -100000;
				double qPCRMidPoint = windowMidPoint(qPCR);
				int minLengthDiff = 1000000;
				//identify overlapping test window closest to center of qPCR region.
				for (int j=0; j< numTestWin; j++){
					Window test = testWindows[j];
					if (qPCR.overlap(test)){
						double diff = Math.abs(qPCRMidPoint - windowMidPoint(test));
						//set better window?
						if (diff < minMidPointDiff ){
							//closer to center
							minMidPointDiff = diff;
							maxScore = test.getScores()[scoreIndex];
							qPCRWindows[i].bestOverlappingWindow = test;
							int lengthTest = test.getStartLastOligo()- test.getStart1stOligo()+1;
							minLengthDiff = Math.abs(lengthTest - qPCRLength);
						}
						else if (diff == minMidPointDiff){
							//same distance from center
							//check length diff
							int lengthTest = test.getStartLastOligo()- test.getStart1stOligo()+1;
							int lengthDiff = Math.abs(lengthTest - qPCRLength);
							if (lengthDiff< minLengthDiff) {
								qPCRWindows[i].bestOverlappingWindow = test;
								maxScore = test.getScores()[scoreIndex];
								minLengthDiff = lengthDiff;
							}
							else if (lengthDiff == minLengthDiff && maxScore < test.getScores()[scoreIndex]) {
								qPCRWindows[i].bestOverlappingWindow = test;
								maxScore = test.getScores()[scoreIndex];
							}
						}
						//otherwise just skip
					}
				}
				//print out qPCR window info with maxScore
				String finalScore;
				if (maxScore == -100000) finalScore = "No Overlapping Windows, possibly trimmed when scanning?";
				else finalScore = Num.formatNumber(maxScore, 4);
				//text rep qPCR region + it's region scores + best window text rep and its	
				StringBuffer sb = new StringBuffer();
				sb.append(qPCR.stringRep(0));
				if (qPCR.getScores() != null) sb.append("\t"+ Num.doubleArrayToStringOnlyMax(qPCR.getScores(), 4, "\t"));
				sb.append("\t");
				sb.append(finalScore);
				//best window found? print info
				Window best = qPCRWindows[i].bestOverlappingWindow;
				if (best != null){
					sb.append("\t");
					sb.append(best.stringRep(0));
					sb.append("\t");
					sb.append(Num.doubleArrayToStringOnlyMax(best.getScores(), 4, "\t"));
				}
				System.out.println(sb);
			}
		}

	}
	public static double windowMidPoint(Window w){
		return w.getStart1stOligo() + (  ((double)(w.getStartLastOligo()-w.getStart1stOligo())) / 2.0  );
	}

	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File windowFile = null;
		File regionFile = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'w': windowFile = new File(args[i+1]); i++; break;
					case 'r': regionFile = new File(args[i+1]); i++; break;
					case 'z': sizeOfOligoMinusOne =Integer.parseInt(args[i+1])-1;i++; break;
					case 'i': scoreIndex=Integer.parseInt(args[i+1]); i++; break;
					case 'm': multiplyByNegativeOne = true; break;
					case 's': findBestScoringWindow = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}

		//fetch windows and sort
		if (windowFile == null || regionFile == null){
			Misc.printExit("\nError: please enter a Window[] file and a text region file.\n");
		}
		System.out.println("\nLoading Window[] ("+windowFile.getName()+") and regions ("+regionFile.getName()+"), score index = "+scoreIndex+"\n");
		testWindows = (Window[])IO.fetchObject(windowFile);
		qPCRWindows = makeQPCRWindows( makeWindows(regionFile));

	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Best Window Score Extractor: September  2008                 **\n" +
				"**************************************************************************************\n" +
				"BWSE prints the best window score for every region. Preference is given to\n" +
				"the best centered and closest to the size to the region. Use the -s flag to bind the\n" +
				"best positive scoring window.\n\n"+

				"Parameters:\n" +
				"-w Full path file text for the serialized Window[] generated by ScanChip or\n" +
				"      ScanChromosome.\n" +
				"-r Full path file text for the tab delimited text regions file (chrom, start, stop, \n" +
				"      ... etc.\n" +
				"-z Size of oligo, defaults to 25.\n"+
				"-i Score index to use in extracting Window scores see ScanChip, defaults to 1.\n" +
				"-s Find best scoring window within each region, defaults to best centered, similar size.\n"+
				"-m Multiply scores by -1\n\n" +

				"Example: java -Xmx1500M -jar pathTo/Apps/BestWindowScoreExtractor -w\n" +
				"      /affy/wins -r /affy/qPCRRes.txt -i 0 -s\n\n" +

		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BestWindowScoreExtractor(args);
	}

	private class QPCRWindow{
		private Window window;
		private Window bestOverlappingWindow;

		public QPCRWindow(Window window){
			this.window = window;
		}
	}

	public QPCRWindow[] makeQPCRWindows(Window[] ws){
		QPCRWindow[] qs = new QPCRWindow[ws.length];
		for (int i=0; i<ws.length; i++){
			qs[i] = new QPCRWindow(ws[i]);
		}
		return qs;
	}

	/**Converts a tab delimited txt file of chr, start, stop, scores (multiple permitted or none)
	 * into an array of Windows.*/
	public static Window[] makeWindows(File txtFile){
		Window[] windows = null;
		try{
			BufferedReader in = new BufferedReader(new FileReader(txtFile));
			String line;
			String[] tokens;
			ArrayList al = new ArrayList(100);
			while ( (line=in.readLine()) !=null){
				if (line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 3) continue;
				String chr = tokens[0];
				int start = Integer.parseInt(tokens[1]);
				int end = Integer.parseInt(tokens[2]);
				//parse scores
				double[] scores = null;
				if (tokens.length > 3){
					String[] s = new String[tokens.length-3];
					System.arraycopy(tokens, 3, s, 0, s.length);
					scores = Num.parseDoubles(s);
				}
				Window win = new Window(chr, start, end, 0, scores);
				al.add(win);
			}
			in.close();
			windows = new Window[al.size()];
			al.toArray(windows);
		}catch (IOException e){
			e.printStackTrace();
		}
		return windows;
	}

}
