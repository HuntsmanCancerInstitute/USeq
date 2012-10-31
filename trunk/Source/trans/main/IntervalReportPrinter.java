package trans.main;
import java.io.*;
import java.text.NumberFormat;
import java.util.regex.*;
import java.util.*;

import util.gen.*;

/**
 * Prints all things contained within an interval to a spreadsheet or a report card.
 *
 */
public class IntervalReportPrinter {
	//fields
	private File[] files;
	private double scoreCutOff = -10000000;
	private int rankCutOff = 0;
	private boolean printSequence = false;
	private boolean printReports = true;
	private boolean printBestWindowSequence = false;
	private double coefVarTreatementIntensities =0;	
	private double coefVarControlIntensities =0;	
	private double meanTreatment =0;			
	private double meanControl =0;				
	private double meanSubWindowTreatment = 0;
	private double meanSubWindowControl = 0;
	private double meanNumberOligoMatches = 0;
	private int scoreIndexForDefaultSort = 1;
	private boolean setDefaultSort = true;
	private boolean printCoordinates = false;
	private int numberOfScores;
	
	public void printReport(){
		//print header
		System.out.println("\nInterval Report: "+new Date());
		System.out.println("\nLikely best window scores : "+Misc.stringArrayToString(ScanChromosomesCNV.getSCORE_DESCRIPTION(),", "));
		System.out.println();
	
		//for each file fetch intervals
		int numFiles = files.length;
		Interval[] intervals;
		int numIntervals;
		for (int i=0; i< numFiles; i++){
			try{
				//make print writer
				PrintWriter out = new PrintWriter(new FileWriter(files[i].getCanonicalPath()+".xls"));
				//attempt to get intervals
				intervals = (Interval[])IO.fetchObject(files[i]);
				System.out.println("Processing Interval[] file: "+files[i].getCanonicalPath());
				if (intervals.length == 0) {
					System.out.println("\tNo Intervals, skipping!");
					continue;
				}
				//check default score index
				int numScores = intervals[0].getBestWindow().getScores().length;
				if (scoreIndexForDefaultSort>= numScores) scoreIndexForDefaultSort = numScores-1;
				//sort intervals
				numIntervals = intervals.length;
				//use best sub window or zero
				SubWindow sub;
				boolean warningIssued = false;
					for (int j=0; j< numIntervals; j++){
						sub = intervals[j].getBestSubWindow();
						if (sub != null)intervals[j].setSortBy(sub.getMedianRatio());
						else {
							if (setDefaultSort){
								Window best = intervals[j].getBestWindow();
								intervals[j].setSortBy(best.getScores()[scoreIndexForDefaultSort]);
								if (warningIssued == false){
									System.out.println("Warning: No sub windows, setting sort by to best window default score index "+scoreIndexForDefaultSort);
									warningIssued = true;
								}
							}
							else {
								if (warningIssued == false){
									System.out.println("Warning: No sub window, setting sort by to zero.");
									warningIssued = true;
								}
								intervals[j].setSortBy(0);
							}
						}
					}	
				
				Arrays.sort(intervals);
				//set number cut offs
				if (rankCutOff != 0 && rankCutOff<numIntervals) numIntervals = rankCutOff;
				
				//print header for spread sheet?
				if (printReports==false && printCoordinates == false) {
					StringBuffer sb = new StringBuffer();
					sb.append("Interval #\t");
					sb.append("Chromosome\t");
					sb.append("Start interval\t");
					sb.append("Stop interval\t");
					sb.append("Length interval\t");
					sb.append("# Oligos in interval\t");
					sb.append("# Windows in interval\t");
					
					sb.append("Size best window (BW)\t");
					sb.append("Start BW\t");
					sb.append("End BW\t");
					sb.append("# Oligos in BW\t");
					sb.append("Mean # Oligo Matches in BW\t");
					sb.append("Mean BW treatment ints\t");
					sb.append("Mean BW control ints\t");
					sb.append("Coef of var treat ints, BW\t");
					sb.append("Coef of var control ints, BW\t");
				
					//scores
					double[] scores = intervals[0].getBestWindow().getScores();
					numberOfScores=scores.length;  //printing just first two
					for (int x=0; x<numberOfScores; x++){
						sb.append("BW Score[");
						sb.append(x);
						sb.append("]\t");
					}
					
					sb.append("Start best sub window (BSW)\t");
					sb.append("Stop BSW\t");
					sb.append("# Oligos in BSW\t");
					sb.append("Mean BSW treatment ints\t");
					sb.append("Mean BSW control ints\t");
					sb.append("Median ratio, BSW\t");
					
					sb.append("BindingPeaks, position (smoothed ratio)\t");
					
					if (intervals[0].getBaseScores() != null){
					sb.append("Motif hits per KB\t");
					sb.append("Motif hits in max cluster\t");
					sb.append("Sum motif hit scores/ Length\t");
					}
					
					//region intersections?
					if (intervals[0].getFractionIntersections() != null){
						String[] names = intervals[0].getRegionNames();
						sb.append("Fraction ");
						sb.append(names[0]);
						for (int z=1; z< names.length; z++){
							sb.append("\tFraction ");
							sb.append(names[z]);
						}
						sb.append("\t");
					}
					
					if (Misc.isNotEmpty(intervals[0].getSequence())){
						sb.append("Sequence (Interval or Best BAIW)");
					}
					out.println(sb);
				}
				for (int j=0; j< numIntervals; j++){
					if (intervals[j].getSortBy()>= scoreCutOff) {
						if (printCoordinates){
							out.println(intervals[j].getChromosome()+"\t"+intervals[j].getStart1stOligo()+"\t"+
									(intervals[j].getStartLastOligo()+intervals[j].getSizeOfOligoMinusOne()));
						}
						else {
							generateIntervalStats(intervals[j]);
							if (printReports) out.println(fetchIntervalReport(intervals[j],j+1));
							else out.println(fetchIntervalLine(intervals[j],j+1));
						}
					}
					else break;
				}
				out.close();
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
		
	}

	
	
	/**Calculate stats on oligos from best window.*/
	public void generateIntervalStats (Interval interval){
		Oligo[] oligosBestWindow = interval.extractSubRegionOligos(interval.getBestWindow().getStart1stOligo(),interval.getBestWindow().getStartLastOligo());	
		if (oligosBestWindow==null) {
			meanTreatment =0;		
			meanControl =0;				
			meanSubWindowTreatment = 0;
			meanSubWindowControl = 0;
			meanNumberOligoMatches = 0;
			return;
		}
		meanNumberOligoMatches = calculateMeanNumberOligoMatches (oligosBestWindow);
		int numOligos = oligosBestWindow.length;
		float[][] treatments = new float[numOligos][];
		float[][] controls = new float[numOligos][];	
		for (int i=0; i<numOligos; i++){
			treatments[i] = oligosBestWindow[i].getTreatmentIntensities(interval.getNumberTreatmentIntensities());
			controls[i] = oligosBestWindow[i].getControlIntensities(interval.getNumberControlIntensities());
		}
		//calc standard deviation/ mean, fold diff
		float[] tCombine = Num.collapseFloatArray(treatments);
		float[] cCombine = Num.collapseFloatArray(controls);

		meanTreatment = Num.averageFloatArray(tCombine);
		meanControl = Num.averageFloatArray(cCombine);
		coefVarTreatementIntensities = Num.standardDeviation(tCombine, meanTreatment)/meanTreatment;
		coefVarControlIntensities = Num.standardDeviation(cCombine, meanControl)/meanControl;
		
		//calc stats for sub window
		//resets treatment and control [][]s!
		SubWindow sub = interval.getBestSubWindow();
		if (sub != null){
			oligosBestWindow = sub.getOligos();
			numOligos = oligosBestWindow.length;
			treatments = new float[numOligos][];
			controls = new float[numOligos][];
			for (int i=0; i<numOligos; i++){
				treatments[i] = oligosBestWindow[i].getTreatmentIntensities(interval.getNumberTreatmentIntensities());
				controls[i] = oligosBestWindow[i].getControlIntensities(interval.getNumberControlIntensities());
			}
			tCombine = Num.collapseFloatArray(treatments);
			cCombine = Num.collapseFloatArray(controls);
			meanSubWindowTreatment = Num.averageFloatArray(tCombine);
			meanSubWindowControl = Num.averageFloatArray(cCombine);
		}
		else {
			meanSubWindowTreatment = 0;
			meanSubWindowControl = 0;
		}
	}
	
	/**Calculates the average number of oligo matches to the genome.*/
	public static double calculateMeanNumberOligoMatches(Oligo[] oligos){
		double total = 0;
		for (int i=0; i< oligos.length; i++) total += oligos[i].getMatches();
		return total/((double)oligos.length);
	}
	
	/**Used when making mock Intervals from which you want to calculate various oligo based measurements.*/
	public void printRatios (Interval interval){
		Oligo[] oligos = interval.getOligos();
		int numOligos = oligos.length;
		double[] aveRatios = new double[numOligos]; 
		float[][] treatments = new float[numOligos][];
		float[][] controls = new float[numOligos][];	
		for (int i=0; i<numOligos; i++){
			treatments[i] = oligos[i].getTreatmentIntensities(interval.getNumberTreatmentIntensities());
			controls[i] = oligos[i].getControlIntensities(interval.getNumberControlIntensities());
			double aveT = Num.mean(treatments[i]);
			double aveC = Num.mean(controls[i]);
			aveRatios[i] = aveT/aveC;		
		}
		double trMean = Num.trimmedMean(aveRatios, 1);
		System.out.println(interval.getChromosome()+"\t"+interval.getStart1stOligo()+"\t"+interval.getStartLastOligo()+"\t"+interval.getBestWindow().getScores()[0]+"\t"+trMean);
		
	}
	
	public String fetchIntervalLine(Interval i, int intervalNumber){
		Window window = i.getBestWindow();
		SubWindow subWindow = i.getBestSubWindow();
		Oligo[] oligos = i.getOligos();
		Oligo[] subOligos = null;
		if (subWindow != null) subOligos = subWindow.getOligos();
		StringBuffer sb = new StringBuffer();
		String tab = "\t";
		//Interval #
		sb.append(intervalNumber); sb.append(tab);
		//Chromosome
		sb.append(i.getChromosome()); sb.append(tab);
		//Start Interval
		sb.append(i.getStart1stOligo()); sb.append(tab);
		//Stop Interval
		sb.append((i.getStartLastOligo()+ i.getSizeOfOligoMinusOne())); sb.append(tab);
		//Length Interval
		sb.append((1+ i.getStartLastOligo()+ i.getSizeOfOligoMinusOne() - i.getStart1stOligo())); sb.append(tab);
		//# Oligos in Interval
		if (oligos != null) sb.append(oligos.length); 
		sb.append(tab);//always add to keep spacing
		//# Windows in Interval
		sb.append(i.getNumberOfWindows()); sb.append(tab);
		
		//Size BW
		sb.append((1+ window.getStartLastOligo()+ i.getSizeOfOligoMinusOne() - window.getStart1stOligo())); sb.append(tab);
		//Start BW
		sb.append((window.getStart1stOligo())); sb.append(tab);
		//Stop BW
		sb.append((window.getStartLastOligo()+ i.getSizeOfOligoMinusOne())); sb.append(tab);
		//# Oligos
		sb.append(window.getNumberOligos()); sb.append(tab);
		//mean number of oligo matches to genome, should be 1, much more indicates a repeat
		sb.append(Num.formatNumber(meanNumberOligoMatches,1)); sb.append(tab);
		//averageMeanTreatment, Mean BAIW Treatment Ints
		sb.append(Num.formatNumber(meanTreatment,1)); sb.append(tab);
		//averageMeanControl, Mean BAIW Control Ints
		sb.append(Num.formatNumber(meanControl,1)); sb.append(tab);
		
		//stats on BW window
		if (oligos != null){
			//coefficient of Variation for Treatment Intensities
			sb.append(Num.formatNumber(coefVarTreatementIntensities,4)); sb.append(tab);
			//coefficient of Variation for Control Intensities
			sb.append(Num.formatNumber(coefVarControlIntensities,4)); sb.append(tab);
		}
		else {
			sb.append(tab); sb.append(tab);  
		}
		
		//scores print scores
		double[] scores = window.getScores();
		//sb.append(Num.doubleArrayToStringOnlyMax(scores,3,"\t")); sb.append(tab);
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(3);
		
		sb.append(f.format(scores[0]));
		for (int k=1; k<numberOfScores; k++){
			sb.append("\t");
			if (k<scores.length)sb.append(f.format(scores[k]));
		}
		sb.append(tab);
			
		
		//sub window
		if (subWindow != null){
			//start
			sb.append(subOligos[0].getStart()); sb.append(tab);
			//stop
			sb.append(subOligos[subOligos.length-1].getStart()+ i.getSizeOfOligoMinusOne()); sb.append(tab);
			//number of oligos
			sb.append(subOligos.length);sb.append(tab);
			//Mean treatment values sub window
			sb.append(Num.formatNumber(meanSubWindowTreatment, 1));sb.append(tab);
			//Mean control values sub window
			sb.append(Num.formatNumber(meanSubWindowControl, 1));sb.append(tab);
			//Median Ratio Score, Sub Window
			sb.append(Num.formatNumber(subWindow.getMedianRatio(),2)); sb.append(tab);
		}
		else {
			sb.append(tab); sb.append(tab); sb.append(tab); sb.append(tab); sb.append(tab); sb.append(tab); 
		}
		
		//binding peaks
		BindingPeak[] bp = i.getBindingPeaks();
		if (bp != null){
			int num = bp.length;
			sb.append(bp[0].getPeakBP());
			sb.append("(");
			sb.append(Num.formatNumberOneFraction(bp[0].getScore()));
			sb.append(")");
			for (int x= 1; x<num; x++){
				sb.append(", ");
				sb.append(bp[x].getPeakBP());
				sb.append("(");
				sb.append(Num.formatNumberOneFraction(bp[x].getScore()));
				sb.append(")");
			}
		}
		sb.append(tab);
		
		//motif info
		if (i.getBaseScores() != null){
			int sizeSeq;
			if (i.isBestWindowScored()) sizeSeq = window.getStartLastOligo()+ i.getSizeOfOligoMinusOne() - window.getStart1stOligo();
			else sizeSeq = i.getSequence().length();
			//Motif Hits Per KB
			sb.append(Num.formatNumber((1000* (double)i.getNumberMotifHits() / (double)sizeSeq),3)); sb.append(tab);
			//Motif Hits in Max Cluster
			sb.append(i.getMaxCluster()); sb.append(tab);
			//Sum motif hits
			sb.append(Num.formatNumber(i.getMotifHitSum()/sizeSeq,5));
			sb.append(tab);
		}
		
		//region intersections?
		if (i.getFractionIntersections() != null){
			double[] intersections = i.getFractionIntersections();
			if (intersections[0] == 0)sb.append ("0");
			else sb.append (Num.formatNumber(intersections[0], 2));
			for (int z=1; z< intersections.length; z++){
				sb.append("\t");
				if (intersections[z] == 0)sb.append ("0");
				else sb.append (Num.formatNumber(intersections[z], 2));
			}
			sb.append(tab);
		}
		
		//print sequence?
		if (Misc.isNotEmpty(i.getSequence()) && printSequence){
				if (printBestWindowSequence) sb.append(i.getBestWindowSequence());
				else sb.append(i.getSequence());

		}
		
		return sb.toString();
	}	
	
	public String fetchIntervalReport(Interval i, int intervalNumber){
		Window window = i.getBestWindow();
		File[] files = i.getCelFiles();
		int numFiles = 0;
		if (files != null) numFiles = files.length;
		SubWindow subWindow = i.getBestSubWindow();
		Oligo[] oligos = i.getOligos();
		StringBuffer sb = new StringBuffer();
		
		//coordinates
		sb.append("\nInterval "); sb.append(intervalNumber); sb.append(" "); sb.append(i.getChromosome()); sb.append(":"); 
		sb.append(i.getStart1stOligo()); sb.append("-"); sb.append((i.getStartLastOligo()+ i.getSizeOfOligoMinusOne()));
		int lengthInterval = i.getStartLastOligo()+ i.getSizeOfOligoMinusOne() +1 - i.getStart1stOligo();
		sb.append(", "); sb.append(lengthInterval); sb.append("bp");
		if (oligos != null){
			sb.append(", "); sb.append(oligos.length); sb.append(" oligos, ");
		}
		//window info
		sb.append(i.getNumberOfWindows());sb.append(" windows\n");
		
		sb.append("\nBest Window: ");   
		sb.append(window.getStart1stOligo()); sb.append("-");
		sb.append((window.getStartLastOligo()+ i.getSizeOfOligoMinusOne()));
		lengthInterval = 1 + window.getStartLastOligo()+ i.getSizeOfOligoMinusOne() - window.getStart1stOligo();
		sb.append(", "); sb.append(lengthInterval); sb.append("bp, "); 
		sb.append(window.getNumberOligos()); sb.append(" # oligos, ");
		sb.append(" Scores: "); sb.append(Num.doubleArrayToStringOnlyMax(window.getScores(), 3, ", "));
		
		
		//sub window information
		if (subWindow != null){
			Oligo[] subOligos = subWindow.getOligos();
			sb.append("\n\nBest Sub Window: "); sb.append(subOligos[0].getStart()); sb.append("-"); 
			sb.append((subOligos[subOligos.length-1].getStart() + i.getSizeOfOligoMinusOne()));
			lengthInterval = 1+ subOligos[subOligos.length-1].getStart() + i.getSizeOfOligoMinusOne() - subOligos[0].getStart();
			sb.append(", "); sb.append(lengthInterval); sb.append("bp, "); sb.append(subOligos.length);sb.append(" oligos, ");
			sb.append("Median Ratio: "); sb.append(Num.formatNumber(subWindow.getMedianRatio(), 3));
			//Mean treatment values sub window
			sb.append(", Mean Treatment: "); sb.append(Num.formatNumber(meanSubWindowTreatment, 1));
			//Mean control values sub window
			sb.append(", Mean Control: "); sb.append(Num.formatNumber(meanSubWindowControl, 1));
		}
		else sb.append("\n\nNo Sub Window.");
		//oligo info
		if (oligos!=null){
			//files
			sb.append("\n\nOligo intensity information:");
			sb.append("\n\tOrder of Extracted Cel Files: ");
			sb.append(files[0].getName());
			for (int x=1; x<numFiles; x++){
				sb.append(", ");
				sb.append(files[x].getName());
			}
			//for each oligo print start, mismatches, treatmentIntensities, controlIntensities, sequence
			sb.append("\n\nStart\tTreatment Intensities\tControl Intensities\tSequence");
			int numOligos = oligos.length;
			for (int x=0; x< numOligos; x++){
				sb.append("\n");
				sb.append(oligos[x].getStart()); sb.append("\t");
				sb.append(Misc.floatArrayToString(oligos[x].getTreatmentIntensities(i.getNumberTreatmentIntensities()),","));
				sb.append("\t");
				sb.append(Misc.floatArrayToString(oligos[x].getControlIntensities(i.getNumberControlIntensities()),","));
				sb.append("\t"); sb.append(oligos[x].getSequence());
			}
			// oligo stats
			sb.append("\n\nStats on best Window:\n"); 
			//stndDivTreatementIntDivMean
			sb.append(Num.formatNumber(coefVarTreatementIntensities,2)); sb.append("\tCoefficient of the variation for all treatment intensities (Standard Deviation Treatment Intensities/ Mean)\n");
			//stndDivControlIntDivMean
			sb.append(Num.formatNumber(coefVarControlIntensities,2)); sb.append("\tCoefficient of the variation for all control intensities (Standard Deviation Control Intensities/ Mean)\n");
			//averageMeanTreatment
			sb.append(Num.formatNumber(meanTreatment,1)); sb.append("\tMean of Treatment Intensities\n");
			//averageMeanControl		
			sb.append(Num.formatNumber(meanControl,1)); sb.append("\tMean of Control Intensities\n");
		}
		else sb.append("\n\nNo oligo information.  Run LoadIntervalOligoInfo if so desired and reprint the report.");
		
		//binding peaks
		BindingPeak[] bps = i.getBindingPeaks();
		if (bps!=null){
			sb.append("\nBinding Peaks (Smoothed Ratio Score, Left Boundary Estimate, Peak, Right Boundary Estimate):\n");
			for (int x=0; x<bps.length;x++){
				sb.append("\t");
				sb.append(bps[x].toString(oligos, i.getSizeOfOligoMinusOne()));
				sb.append("\n");
			}
		}
		
		//motif hit info
		if (i.getBaseScores() != null){
			String word = "Interval";
			if (i.isBestWindowScored()) word = "Best Window";
			sb.append("\nMotif information for "+word+":\n");
			sb.append(i.getNumberMotifHits());
			sb.append(" # Motif hits above cut off\n");
			sb.append(i.getMaxCluster());
			sb.append(" Max cluster containing motif hits above cut off\n");
			int sizeSeq;
			if (i.isBestWindowScored()) sizeSeq = window.getStartLastOligo()+ i.getSizeOfOligoMinusOne() - window.getStart1stOligo();
			else sizeSeq = i.getSequence().length();
			sb.append(Num.formatNumber((1000* (double)i.getNumberMotifHits() / (double)sizeSeq),3));
			sb.append(" Hits above cut off per KB\n");
			sb.append(Num.formatNumber(i.getMotifHitSum()/sizeSeq,5));
			sb.append(" Sum of the logged motif hit scores over the "+word+" / length.\n");
		}
		else sb.append("\n\nNo motif or sequence information.  Run ScoreIntervals if so desired and reprint the report.");
		
		//print seq?
		if (printSequence){
			sb.append("\n\nInterval Sequence:\n");
			sb.append(i.getSequence());
			sb.append("\n\nBest Window Sequence:\n");
			sb.append(i.getBestWindowSequence());
			if (i.getBestSubWindow() != null){
				sb.append("\n\nBest Sub Window Sequence:\n");
				sb.append(i.getBestSubWindowSequence());
			}
		}

		sb.append("\n\n");	
		return sb.toString();
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Interval Report Printer: April 2007                       **\n" +
				"**************************************************************************************\n" +
				"IRP prints reports in spread sheet format or as detailed pages for Interval[] arrays.\n" +
				"Intervals are sorted by the median ratio of the best sub window.\n" +
				"Use the following options:\n\n" +
				
				"-f Full path file text for the Interval[], if a directory is specified, all files\n" +
				"      within will be processed. (Required)\n" +
				"-s Score cut off, print everything above this score, defaults to all\n" +
				"-r Rank cut off, prints everything above this rank, ie the top 200, defaults to all\n" +
				"-i Score index to use in sorting when there is no sub window, defaults to 1.\n"+
				"-p Print sequences, default is no\n"+
				"-a Print best window sequence in summary line instead of interval sequence.\n"+
				"-b Print tab delimited summary line, default is to print a detailed report.\n"+
				"-c Print coordinates (chrom, start, stop).\n"+
				
				"\n" +
				"Example: java -jar pathTo/T2/Apps/IntervalReportPrinter -f /my/affy/res/ -s 1.5\n" +
				"      -r 200 -p -b\n" +
				"\n" +
				
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
					case 'i': scoreIndexForDefaultSort = Integer.parseInt(args[i+1]); break;
					case 'p': printSequence=true; break;
					case 'a': printBestWindowSequence=true; printSequence=true; break;
					case 'b': printReports=false; break;
					case 'c': printCoordinates=true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		//check to see if they entered required params
		if (directory==null || directory.exists() == false){
			System.out.println("\nEnter a serialized Interval[] array file or directory!\n");
			System.exit(0);
		}
		
		//get files to process
		files = IO.extractFiles(directory);
	}
	
	public IntervalReportPrinter (String[] args){
		processArgs(args);
		printReport();
		
	}
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new IntervalReportPrinter(args);
	}
	
	/**Calculate Median ratio of SubWindow where the treatments are averaged, and the controls are averaged*/
	public static double medianRatioSubWindow(Interval  i){
		Oligo[] oligos = i.getBestSubWindow().getOligos();
		int numTreatments = i.getNumberTreatmentIntensities();
		int numControls = i.getNumberControlIntensities();
		//calculate ratios	
		double[] ratios = new double[oligos.length];
		for (int j=oligos.length-1; j>=0; j--){
			ratios[j] = Num.mean(oligos[j].getTreatmentIntensities(numTreatments))/ Num.mean(oligos[j].getControlIntensities(numControls));
		}
		
		Arrays.sort(ratios);
		return Num.median(ratios);
	}
	/**Calculate maximum median absolute diff of SubWindow treatment values*/
	public static double maxMedianAbsDiffSubWinTreatment(Interval  i){
		int numTreatments = i.getNumberTreatmentIntensities();
		if (numTreatments <2) return 0;
		Oligo[] oligos = i.getBestSubWindow().getOligos();
		int numOligos = oligos.length;
		//build arrays int[oligo][values]
		float[][] ovT = new float[numOligos][];
		for (int j=0; j<numOligos; j++){
			ovT[j] = oligos[j].getTreatmentIntensities(numTreatments);
		}
		return Num.calcMaxMedianAbsoluteDifference(ovT);
	}

	/**Calculate maximum median absolute diff of SubWindow control values*/
	public static double aveMedianAbsDiffSubWinControl(Interval  i){
		int numControls = i.getNumberControlIntensities();
		if (numControls <2) return 0;
		Oligo[] oligos = i.getBestSubWindow().getOligos();
		int numOligos = oligos.length;
		//build arrays int[oligo][values]
		float[][] ov = new float[numOligos][];
		for (int j=0; j<numOligos; j++){
			ov[j] = oligos[j].getControlIntensities(numControls);
		}
		return Num.calcMeanMedianAbsoluteDifference(ov);
	}
}
