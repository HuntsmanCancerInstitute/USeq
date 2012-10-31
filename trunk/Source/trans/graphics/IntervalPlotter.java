package trans.graphics;

import javax.swing.*;
import trans.main.*;
import util.gen.*;
import java.awt.Color;
import java.awt.Dimension;
import java.io.*;
import java.util.*;
import java.util.regex.*;


/**
 * Application to draw graphical representations of loaded Intervals.
 *
 */
public class IntervalPlotter {
	//fields
	private File[] intervalFiles;
	private double scoreCutOff =0;
	private int rankCutOff = 0;
	private int rankCutOffBelow = 0;
	private boolean savePlots = false;
	private boolean showPSPM = true;
	private double scale=1;
	private boolean antiAlias = false;
	
	public IntervalPlotter (String[] args){
		processArgs(args);
		
		for (int i=0; i<intervalFiles.length; i++){
			//attempt to get intervals
			Interval[] intervals = (Interval[])IO.fetchObject(intervalFiles[i]);
			
			//did they score the intervals for hits to a pwm? if not hide em
			if (intervals[0].getBaseScores()==null) showPSPM = false;
			
			//sort intervals
			//set sortBy field
			int numIntervals = intervals.length;
			
			SubWindow sub;
			for (int j=0; j< numIntervals; j++){
				sub = intervals[j].getBestSubWindow();
				if (sub != null)intervals[j].setSortBy(sub.getMedianRatio());
				else intervals[j].setSortBy(0);
			}	
			Arrays.sort(intervals);
			
			//set number cut offs
			if (rankCutOff != 0 && rankCutOff<numIntervals) numIntervals = rankCutOff;
			if (rankCutOffBelow !=0) rankCutOffBelow--;
			
			//make directory to hold images named after Interval[] + Plots
			String intervalName = intervalFiles[i].getName();
			
			//make window for each interval
			if (savePlots){
				//kill the need to display, don't instantiate a Frame
				System.setProperty("java.awt.headless","true");
				
				File saveDirectory = new File(intervalFiles[i].getParent(), intervalName+"Plots");
				saveDirectory.mkdir();
				IntervalDrawPanel panel = null;
				System.out.print("Rendering and saving plots ");
				for (int j=rankCutOffBelow; j< numIntervals; j++){
					if (intervals[j].getSortBy()>= scoreCutOff) {
						System.out.print(".");
						panel = new IntervalDrawPanel(intervals[j], j+1, null, showPSPM);
						panel.setBackground(Color.BLACK);
						panel.setPreferredSize(new Dimension((int)panel.getRunningX()+10, (int)panel.getRunningY()+10));
						panel.saveBufferedImage(scale, new File(saveDirectory.getPath(), (j+1)+"_"+intervals[j].getChromosome()+
								"_"+intervals[j].getStart1stOligo()+"_"+ (intervals[j].getStartLastOligo()+
										intervals[j].getSizeOfOligoMinusOne())+".png"), antiAlias);
					}
					else break;
				}
				System.out.println();
			}
			else {
				TextFrame textFrame = new TextFrame(400, 0, 250, 600, "Coordinate");
				for (int j=numIntervals-1; j>= rankCutOffBelow; j--){
					if (intervals[j].getSortBy()>= scoreCutOff) {
						IntervalDrawFrame frame = new IntervalDrawFrame(intervals[j], j+1, textFrame, showPSPM);
						frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						frame.show();
					}
					else break;
				}
				textFrame.show();
			}
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Interval Plotter:     Nov 2005                          **\n" +
				"**************************************************************************************\n" +
				"IP plots and saves graphs for serialized Interval[] arrays.  Load the Intervals with\n" +
				"data from the FindSubBindingRegion, LoadIntervalOligoInformation, and optionally\n" +
				"ScoreIntervals programs prior to running the IntervalPlotter. Plotted are graphs for\n" +
				"averaged treatment intensity, the averaged control intensity, the average difference,\n" +
				"the average fold difference, a smoothed trimmed mean ratio, peak picks, the best\n" +
				"window, the best sub window, centered PSPM hit scores, the number of 1bp mis/matches\n" +
				"for each oligo in the genome, and the treatment and control intensities for each\n" +
				"individual processed cel file. Click and or drag to fetch the coordinates and sequence\n" +
				"for a selected region. Columns in console: interval rank, median ratio best sub window,\n" +
				"trimmed mean score for the closest oligo if just one click or the max trimmed mean\n" +
				"score within the dragged box, chromosome, start, stop, and sequence. Intervals are\n" +
				"sorted by the best median ratio sub window.\n\n" +
				
				"-f Full path file or directory text for the data loaded Interval[](s).\n" +
				"-s Score cut off, plot everything above this score, defaults to all.\n" +
				"-r Rank cut off, plot everything above this rank, ie the top 200, defaults to all.\n" +
				"-q Rank cut off, plot everything below this rank, defaults to all.\n"+
				"-p Save plots to disk, default is no.\n"+
				"-m Magnify saved plots (2 twice as big, 3 three times as big...), default is none.\n"+
				"-a Anti alias saved plots, good for printed figures, bad for computer display,\n" +
				"      default is no.\n"+
				"-b Hide the PSPM graphs.\n"+
				"\n" +
				"Example: java -jar pathTo/T2/Apps/IntervalPlotter -f /affy/res/Z.resAll -t -s 1.5 -r\n" +
				"      200 -p -m 1.5 -a\n" +
				"\n" +
				"Questions? Contact David_Nix@Affymetrix.com or SuperFly@lbl.gov\n"+
		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File directory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': directory = new File(args[i+1]); i++; break;
					case 's': scoreCutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'r': rankCutOff =Integer.parseInt(args[i+1]); i++; break;
					case 'q': rankCutOffBelow =Integer.parseInt(args[i+1]); i++; break;
					case 'p': savePlots=true; break;
					case 'a': antiAlias=true; break;
					case 'b': showPSPM=false; break;
					case 'm': scale = Double.parseDouble(args[i+1]); i++; break;
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
		if (directory == null || directory.exists() == false){
			System.out.println("\nCould not find your interval file or directory?!\n");
			System.exit(1);
		}
		intervalFiles = IO.extractFiles(directory);
	}
	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}
		new IntervalPlotter (args);
	}
}



