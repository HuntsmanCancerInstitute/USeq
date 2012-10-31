package trans.graphics;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.xy.*;

import util.gen.*;
import java.io.*;
import java.text.NumberFormat;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;


/**Performs set analysis on genomic regions with a visual display.*/
public class RankedSetAnalysis {
	//fields
	private File[] files;
	private File setOneFile;
	private File setTwoFile;
	private GenomicRegion[] regionsOne;
	private GenomicRegion[] regionsTwo;
	private int maxGap = -100;
	private TextFrame textFrame;
	private boolean savePngs = false;
	private int[] rankedHitsOne;
	private int[] rankedHitsTwo;
	private ChartPanel chartPanel;
	private JFreeChart chart;
	private StringBuffer summary = new StringBuffer();
	
	//constructor
	public RankedSetAnalysis(String[] args){
		//process args
		processArgs(args);
		
		//start summary buffer
		summary.append("Fraction Total Intersection:\nComparison\tIntersection Relative to First\t# Regions First\tIntersection Relative to Second\t# Regions Second\n");
		
		//make text frame for notes
		textFrame = new TextFrame(400, 0, 250, 250, "Info");
		textFrame.setVisible(true);
		
		//single pair? or multiple comparisons
		if (files == null) intersect();
		else {
			//set files then call intersect
			for (int i=0; i<files.length; i++){
				setOneFile = files[i];
				for (int j=i+1; j<files.length; j++){
					setTwoFile = files[j];
					intersect();
				}
			}
		}
		
		//print summary
		System.out.println("\n"+summary);
	}
	
	public void intersect(){
		System.out.println("\nComparing "+setOneFile.getName() + " against "+setTwoFile.getName()+ " (MaxGap= "+maxGap+"bp)");
		
		//parse files
		regionsOne = parseGenomicRegionsFile(setOneFile);
		regionsTwo = parseGenomicRegionsFile(setTwoFile);
		
		//intersect regions, loading regionsOne with references to intersecting regionsTwo
		if (intersectRegions() == false){
			System.out.println("\nSorry no regions were found to intersect! Try increasing the max gap.\n");
			return;
		}
		
		//print ranked intersection percents
		//System.out.println("Ranked Intersection for "+setOneFile.getName());
		//System.out.println("Rank\t# With Int\t%");
		rankedHitsOne = rankedIntersection(regionsOne);
		//printRankedIntersection(rankedHitsOne);

System.out.println("\nRanked Intersection for "+setTwoFile.getName());
System.out.println("Rank\t# With Int\t%");
		rankedHitsTwo = rankedIntersection(regionsTwo);
printRankedIntersection(rankedHitsTwo);

System.exit(0);
		
		//print tables of intersecting regions from both perspectives
		System.out.println("\nHits to "+ setOneFile.getName());
		for (int i=0; i< regionsOne.length; i++){
			ArrayList hits = regionsOne[i].getIntersectingRegions();
			if (hits.size() !=0){
				System.out.print(regionsOne[i]);
				for (int j=0; j<hits.size(); j++){
					System.out.print("\t"+(GenomicRegion)hits.get(j));
				}
				System.out.println();
			}
		}
		
		System.out.println("\nHits to "+ setTwoFile.getName());
		for (int i=0; i< regionsTwo.length; i++){
			ArrayList hits = regionsTwo[i].getIntersectingRegions();
			if (hits.size() !=0){
				System.out.print(regionsTwo[i]);
				for (int j=0; j<hits.size(); j++){
					System.out.print("\t"+(GenomicRegion)hits.get(j));
				}
				System.out.println();
			}
		}
		
		//make plot
		makeRankedHitsGraph();
		
		//make panel and display
		RankedSetDrawFrame frame = new RankedSetDrawFrame(this);
		
		//save summary line
		double fractIntOne = fractionTotalIntersection(regionsOne);
		double fractIntTwo = fractionTotalIntersection(regionsTwo);
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(3);
		f.setMinimumFractionDigits(3);
		summary.append(setOneFile.getName()+" Vs "+setTwoFile.getName()+"\t"+f.format(fractIntOne)+"\t"+regionsOne.length+"\t"+f.format(fractIntTwo)+"\t"+regionsTwo.length+"\n");
		
		
	}
	
	public void makeRankedHitsGraph(){
		//data objects
		XYSeries series = new XYSeries(setOneFile.getName());
		XYSeries series2 = new XYSeries(setTwoFile.getName());
		for (int i=0; i<rankedHitsOne.length; i++){
			double rank = i+1;
			series.add(rank, ((double)rankedHitsOne[i])/rank);
		};
		for (int i=0; i<rankedHitsTwo.length; i++){
			double rank = i+1;
			series2.add(rank, ((double)rankedHitsTwo[i])/rank);
		};
		//Add the series to your data set
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		dataset.addSeries(series2);
		//Generate the graph
		chart = ChartFactory.createXYLineChart(
				"Ranked Intersection AnalysisBean", // Title
				"Rank", // x-axis Label
				"Fraction of Intersection", // y-axis Label
				dataset, // Dataset
				PlotOrientation.VERTICAL, // Plot Orientation
				true, // Show Legend
				true, // Use tooltips
				false // Configure chart to generate URLs?
		);
		chartPanel = new ChartPanel(chart);
	}
	
	/**Calculates the number of hits in the ranked matched list.  As if one submitted a standard intersection 
	 * analysis for every incremental rank.*/
	public static int[] rankedIntersection(GenomicRegion[] r){
		int[] hits = new int[r.length];
		for (int i=0; i< r.length; i++){
			int numIntersect =0;
			for (int j = 0; j < i+1; j++){
				if (containsMatchedRankIntersect(r[j], i)) {
					numIntersect++;
				}
			}
			hits[i] = numIntersect;
		}
		return hits;
	}
	
	/**Prints to screen the results of rankedIntersection(GenomicRegion[] r).*/
	public void printRankedIntersection(int[] hits){
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(2);
		f.setMinimumFractionDigits(2);
		for (int i=0; i< hits.length; i++){
			double total = i+1;
			double numIntersect = hits[i];
			if (total == 0) System.out.println((int)total+"\t"+ (int)numIntersect+"\t0");
			else System.out.println((int)total+"\t"+ (int)numIntersect+"\t"+f.format(100*numIntersect/total));
		}
	}
	
	/**Checks to see if among the hits to a particular GenomicRegion they have an equal or smaller rank.*/
	public static boolean containsMatchedRankIntersect(GenomicRegion r, int testRank){
		//any hits?
		ArrayList hits = r.getIntersectingRegions();
		if (hits.size() == 0) return false;
		//scan hits for equal to or better ranked match
		int numHits = hits.size();
		for (int i=0; i<numHits; i++){
			GenomicRegion hit = (GenomicRegion)hits.get(i);	
			if (hit.getRank() <= testRank) return true;
		}
		return false;
	}
	
	/**Returns the total fraction of intersection.*/
	public double fractionTotalIntersection(GenomicRegion[] r){
		double total = 0;
		double numIntersect =0;
		for (int i=0; i< r.length; i++){
			total++;
			if (r[i].getIntersectingRegions().size() !=0) numIntersect++;
		}
		return numIntersect/total;
	}
	
	/**Rank analysis where each regions is check against the entire other list.*/
	public void printRankedIntersectionAll(GenomicRegion[] r){
		double total = 0;
		double numIntersect =0;
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(2);
		f.setMinimumFractionDigits(2);
		for (int i=0; i< r.length; i++){
			total++;
			if (r[i].getIntersectingRegions().size() !=0) numIntersect++;
			System.out.println((int)total+"\t"+ (int)numIntersect+"\t"+f.format(100*numIntersect/total));
		}
	}
	
	/**Loads regions with references to intersection partners.*/
	public boolean intersectRegions(){
		boolean someIntersections = false;
		for (int i=0; i< regionsOne.length; i++){
			for (int j=0; j< regionsTwo.length; j++){
				if (regionsOne[i].overlap(regionsTwo[j], maxGap)){
					regionsOne[i].getIntersectingRegions().add(regionsTwo[j]);
					regionsTwo[j].getIntersectingRegions().add(regionsOne[i]);
					someIntersections = true;
				}
			}
		}
		return someIntersections;
	}
	

	
	/**Parses a regions file.*/
	public static GenomicRegion[] parseGenomicRegionsFile(File f){
		GenomicRegion[] regions =null;
		try{
			BufferedReader in = new BufferedReader(new FileReader(f));
			String line;
			ArrayList regionsAL = new ArrayList();
			int rank = 0;
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.equals("") || line.startsWith("#") || line.startsWith("//")) continue;
				regionsAL.add(new GenomicRegion(line, rank));
				rank++;
			}
			regions = new GenomicRegion[regionsAL.size()];
			regionsAL.toArray(regions);
		}catch (IOException e){
			e.printStackTrace();
		}
		return regions;
	}
	

	
	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		new RankedSetAnalysis(args);
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Ranked Set Analysis: Jan 2006                           **\n" +
				"**************************************************************************************\n" +
				"RSA performs set analysis (intersection, union, difference) on lists of\n" +
				"genomic regions (tab delimited: chrom, start, stop, score, (optional notes)).\n\n" +
				
				"-a Full path file text for the first list of genomic regions.\n" +
				"-b Full path file text for the second list of genomic regions.\n" +
				"-d (Optional) Full path directory containing region files for all pair analysis.\n"+
				"-m Max gap, bps, set negative to force an overlap, defaults to -100\n"+
				"-s Save comparison as a PNG, default is no.\n"+
				"\n" +
				
				"Example: java -jar pathTo/T2/Apps/RankedSetAnalysis -a /affy/nonAmpA.txt -b\n" +
				"      /affy/nonAmpB.txt -s\n"+
				"\n" +

		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File directory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': setOneFile = new File(args[i+1]); i++; break;
					case 'b': setTwoFile = new File(args[i+1]); i++; break;
					case 'd': directory = new File(args[i+1]); i++; break;
					case 'm': maxGap = Integer.parseInt(args[i+1]); i++; break;
					case 's': savePngs = true; break;
					case 'h': printDocs(); System.exit(0);
					default: {
						Misc.printExit("\nProblem, unknown option! " + mat.group()+"\n");
					}
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test+"\n");
				}
			}
		}
		//extract cel files?
		if (directory == null){
			if (setOneFile == null || setTwoFile == null || setOneFile.canRead()==false || setTwoFile.canRead()== false){
				Misc.printExit("\nCannot find or read one of your files!\n");
			}
		}
		else {
			//check directory
			if (directory.canRead() == false || directory.isDirectory()== false){
				Misc.printExit("\nCannot find or read your directory!\n");
			}
			files = IO.extractFiles(directory);
			if (files.length<2) Misc.printExit("\nToo few files in your directory!\n");
		}
	}

	public GenomicRegion[] getRegionsOne() {
		return regionsOne;
	}

	public GenomicRegion[] getRegionsTwo() {
		return regionsTwo;
	}

	public File getSetOneFile() {
		return setOneFile;
	}

	public File getSetTwoFile() {
		return setTwoFile;
	}

	public TextFrame getTextFrame() {
		return textFrame;
	}

	public boolean isSavePngs() {
		return savePngs;
	}

	public ChartPanel getChartPanel() {
		return chartPanel;
	}

	public void setChartPanel(ChartPanel chartPanel) {
		this.chartPanel = chartPanel;
	}

	public JFreeChart getChart() {
		return chart;
	}
}