package trans.main;
import java.io.*;
import java.util.regex.*;
import trans.misc.*;
import util.gen.*;


/**
 * Prints intervals as a GFF3 File.
 *
 */
public class IntervalGFFPrinter {
	//fields
	private File[] files;
	private String name;
	private int scoreIndex = 1;
	private boolean printNonIntervalLines = true;
	private boolean makeSimpleGFF = false;

	public void printGFF(){
		//for each file fetch intervals
		int numFiles = files.length;
		Interval[] intervals;
		int numIntervals;
		try{
			//for each file
			for (int i=0; i< numFiles; i++){
				//make print writer
				String gffExt = ".gff3";
				if (makeSimpleGFF) gffExt = ".gff";
				PrintWriter out = new PrintWriter(new FileWriter(files[i].getAbsolutePath()+gffExt));
				
				//print track info
				//String text = Misc.removeExtension(files[i].getName());
				//out.println("track text="+text+" description=\"GFF export of interval file\" useScore=1");
				
				//print header?
				if (makeSimpleGFF == false){
					out.println("## GFF version 3");
					out.println("## "+Misc.getDate());
					out.println("## Summary scores: intervals use median ratio of the best sub window when present or the user defined best window score, " +
							"best windows use a user defined score index, best sub window uses the median ratio, oligos use log2(aveT/ aveC), " +
					"binding peak use a trimmed mean of window ratios.\n\n");
				}
				
				//attempt to get intervals
				intervals = (Interval[])IO.fetchObject(files[i]);
				System.out.println("Processing Interval[] file: "+files[i].getCanonicalPath());
				name = files[i].getName();
				//sort intervals
				numIntervals = intervals.length;
				Util.sortIntervalsBySubWindowMedianRatio(intervals);

				//for each interval
				for (int j=0; j< numIntervals; j++){
					StringBuffer s = new StringBuffer();
					//print Interval line
					gff(intervals[j], j+1, name, s);
					s.append("\n");

					if (printNonIntervalLines){
						//print Window line
						gff(intervals[j], name+"_"+(j+1), s);
						s.append("\n");

						//print others if oligos were loaded
						if (intervals[j].getOligos()!=null){

							//print SubWindow line
							SubWindow sub = intervals[j].getBestSubWindow();
							if (sub!=null){
								gff(sub, intervals[j], name+(j+1), s);
								s.append("\n");
							}
							//print binding peak lines
							BindingPeak[] peaks = intervals[j].getBindingPeaks();
							if (peaks != null){
								for (int k=0; k<peaks.length; k++){
									gff(peaks[k], intervals[j], name+(j+1), k+1, s);
									s.append("\n");
								}
							}
							//print oligo lines
							Oligo[] oligos = intervals[j].getOligos();
							for (int k=0; k<oligos.length; k++){
								gff(oligos[k], intervals[j], name+(j+1), k+1, s);
								s.append("\n");
							}
						}
					}
					out.print(s.toString());
				}
				out.close();
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}

	/**Appends a GFF3 line for an Oligo. Summary Score is a log2Ratio*/
	public static void gff(Oligo o, Interval i, String name, int oligoNumber, StringBuffer s){
		s.append(i.getChromosome());
		s.append("\t");
		s.append("TiMAT2");
		s.append("\t");
		s.append("Oligo");
		s.append("\t");
		s.append(o.getStart());
		s.append("\t");
		s.append(o.getStart()+i.getSizeOfOligoMinusOne());
		s.append("\t");
		float[] ts = o.getTreatmentIntensities(i.getNumberTreatmentIntensities());
		float[] cs = o.getControlIntensities(i.getNumberControlIntensities());
		double t = Num.mean(ts);
		double c = Num.mean(cs);
		s.append(Num.log2(t/c));
		s.append("\t");
		s.append(".");
		s.append("\t");
		s.append(".");
		s.append("\t");
		//attributes
		s.append("ID=");
		s.append(name);
		s.append("O");
		s.append(oligoNumber);
		s.append("; Name=O");
		s.append(oligoNumber);
		s.append("; Parent=");
		s.append(name);
		s.append("; type=oligo; treatmentScores=");
		s.append(Misc.floatArrayToString(ts,","));
		s.append("; controlScores=");
		s.append(Misc.floatArrayToString(cs,","));
		s.append("; seq=");
		s.append(o.getSequence());
		s.append(";");
	}

	/**Appends a GFF3 line for a BindingPeak*/
	public static void gff(BindingPeak peak, Interval i, String name, int peakNumber, StringBuffer s){
		s.append(i.getChromosome());
		s.append("\t");
		s.append("TiMAT2");
		s.append("\t");
		s.append("BindingPeak");
		s.append("\t");
		s.append(peak.getPeakBP());
		s.append("\t");
		s.append(peak.getPeakBP());
		s.append("\t");
		s.append(peak.getScore());
		s.append("\t");
		s.append(".");
		s.append("\t");
		s.append(".");
		s.append("\t");
		//attributes
		s.append("ID=");
		s.append(name);
		s.append("P");
		s.append(peakNumber);
		s.append("; Name=P");
		s.append(peakNumber);
		s.append("; Parent=");
		s.append(name);
		s.append("; type=peak; leftEdge=");
		Oligo[] oligos = i.getOligos();
		s.append(oligos[peak.getLeftFlankingIndex()].getStart());
		s.append("; rightEdge=");
		s.append(oligos[peak.getRightFlankingIndex()].getStart() + i.getSizeOfOligoMinusOne());
		s.append(";");
	}


	/**Appends a GFF3 line for a SubWindow. Summary score is median ratio*/
	public static void gff(SubWindow sub, Interval i, String name, StringBuffer s){
		s.append(i.getChromosome());
		s.append("\t");
		s.append("TiMAT2");
		s.append("\t");
		s.append("SubWindow");
		s.append("\t");
		Oligo[] oligos = sub.getOligos();
		s.append(oligos[0].getStart());
		s.append("\t");
		s.append(oligos[oligos.length-1].getStart()+i.getSizeOfOligoMinusOne());
		s.append("\t");
		s.append(sub.getMedianRatio());
		s.append("\t");
		s.append(".");
		s.append("\t");
		s.append(".");
		s.append("\t");
		//attributes
		s.append("ID=");
		s.append(name);
		s.append("SW; Name=SW; Parent=");
		s.append(name);
		s.append("; type=bestSubWindow;");	
	}

	/**Appends a GFF3 line for an Interval.*/
	public void gff(Interval i, int rank, String name, StringBuffer s){
		//columns
		s.append(i.getChromosome());
		s.append("\t");
		s.append(name);
		s.append("\t");
		s.append("Interval");
		s.append("\t");
		s.append(i.getStart1stOligo());
		s.append("\t");
		s.append(i.getStartLastOligo()+i.getSizeOfOligoMinusOne());
		s.append("\t");
		SubWindow sub = i.getBestSubWindow();
		if (sub !=null) s.append(sub.getMedianRatio());
		else s.append(i.getBestWindow().getScores()[scoreIndex]);
		s.append("\t");
		s.append(".");
		s.append("\t");
		s.append(".");
		s.append("\t");
		String nameRank = name+"_"+rank;
		if (makeSimpleGFF) s.append(nameRank);
		else{
			//attributes
			s.append("ID="); 
			s.append(nameRank);
			s.append("; Name="); 
			s.append(nameRank);
			s.append("; type=interval; scores="); 
			s.append(Num.doubleArrayToStringOnlyMax(i.getBestWindow().getScores(), 3, ",") );
			s.append(";");
		}
	}

	/**Appends a GFF3 line for a Window.*/
	public void gff(Interval i, String name, StringBuffer s){
		//columns
		s.append(i.getChromosome());
		s.append("\t");
		s.append("TiMAT2");
		s.append("\t");
		s.append("Window");
		s.append("\t");
		Window bw = i.getBestWindow();
		s.append(bw.getStart1stOligo());
		s.append("\t");
		s.append(bw.getStartLastOligo()+i.getSizeOfOligoMinusOne());
		s.append("\t");
		s.append(bw.getScores()[scoreIndex]);
		s.append("\t");
		s.append(".");
		s.append("\t");
		s.append(".");
		s.append("\t");
		//attributes
		s.append("ID=");
		s.append(name);
		s.append("W; Name=W; Parent=");
		s.append(name);
		s.append("; type=bestWindow; scores="); 
		s.append(Num.doubleArrayToStringOnlyMax(bw.getScores(), 3, ",") );
		s.append(";");
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Interval GFF Printer: Feb 2007                            **\n" +
				"**************************************************************************************\n" +
				"IGP prints a GFF3 txt file given a serialized Interval file, sorts by the median ratio\n" +
				"of the best sub window. Use the following options:\n\n" +

				"-f Full path file text for the Interval[], if a directory is specified, all files\n" +
				"      within will be processed and saved to the same GFF3 file.\n" +
				"-d Don't print non interval lines (ie window, oligo).\n"+
				"-s Print simple GFF not GFF3.\n"+
				"-i Score index to use in assigning the best window summary score. See ScanChip.\n\n" +

				"Example: java -jar pathTo/T2/Apps/IntervalGFFPrinter -f /my/affy/res/ -i 2\n" +
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
					case 'i': scoreIndex = Integer.parseInt(args[i+1]); i++; break;
					case 'd': printNonIntervalLines = false; break;
					case 's': makeSimpleGFF = true; break;
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
		//check to see if they entered required params
		if (directory==null || directory.exists() == false){
			System.out.println("\nEnter a serialized Interval[] array file or directory!\n");
			System.exit(0);
		}

		//get files to process
		if (directory.isDirectory()){
			files = IO.extractFiles(directory);
			if 	(files == null || files.length==0){
				System.out.println("Cannot find the directory or files in the directory?!\n\n");
				System.exit(0);	
			}
		}
		else files = new File[]{directory};	

	}

	public IntervalGFFPrinter (String[] args){
		processArgs(args);
		printGFF();
	}
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new IntervalGFFPrinter(args);
	}


}
