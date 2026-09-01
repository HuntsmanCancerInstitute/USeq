package edu.utah.seq.cnv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.hci.misc.Gzipper;
import util.apps.MergeRegions;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class CopyAnalysisCalledAnnoSegFilter {

	//User fields
	private File[] segFilesToParse = null;
	private File resultsDir = null;
	private double minTNLogRatio = Double.MAX_VALUE;
	private double minTLogRatio = Double.MAX_VALUE;
	private double maxNLogRatio = Double.MAX_VALUE;
	private double minGatkLgGeoMean = Double.MAX_VALUE;
	private boolean passGATK = false;
	private String headerLine = "Sample\tChr\tStart\tEnd\tCR_NumPoints\tCR_LgTumorMean\tCR_LgNormalMean\tCR_LgTNRatio\tAF_TumorSize\tAF_TumorMean\tAF_NormalSize\tAF_NormalMean\tAF_TNRatio\tGATK_Call\tGATK_LgTumorGeoMean";
	private boolean outputBed = false;
	
	public CopyAnalysisCalledAnnoSegFilter(String[] args) {

		try {

			processArgs(args);
			
			printSettings();

			parseSegFiles();

			
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: running the CopyAnalysisCalledAnnoSegFilter\n");
		}
		IO.pl("\nDone");

	}
	private void parseSegFiles() throws Exception {
		IO.pl("Seg File\t# Passing\t# Failing");
		String line = null;

		//parse seg files
		for (File segFile: segFilesToParse) {
			IO.p(segFile.getName()+"\t");
			//make IO
			BufferedReader in = IO.fetchBufferedReader(segFile);
			String name = Misc.removeExtension(segFile.getName());
			if (outputBed) name = name+".pass.bed.gz";
			else name = name+".pass.seg.gz";
			Gzipper out = new Gzipper( new File(resultsDir, name));
			int numPass = 0;
			int numFail = 0;

			while ((line = in.readLine())!=null) {
				//header? sometimes it's not present like after running AnnotateBedWithGenes
				if (line.startsWith("Sample") || line.startsWith("#")) {
					if (line.equals(headerLine)==false)	{
						throw new IOException("Failed to find the proper header line? Are these xxx.called.anno.seg files? See:\n"+
								headerLine+"\nVerse:\n"+line);
					}
					if (outputBed == false) out.println(line);
					continue;
				}

				line = line.trim();
				if (line.length()!=0) {
					String[] t = Misc.TAB.split(line);
					boolean pass = check(t);
					if (pass) {
						numPass++;
						if (outputBed) out.println(toBed(t));
						else out.println(line);
					}
					else numFail++;
				}
			}
			in.close();
			out.close();
			IO.pl(numPass+"\t"+numFail);
		}
	}
	
	private String toBed(String[] t) {
		//chr16	  29751	3169750	numOb=305;lg2Tum=-0.1976;lg2Norm=-0.0064;genes=	  -0.1909	   -	
		StringBuilder sb = new StringBuilder();
		sb.append(t[1]); sb.append("\t");
		sb.append(t[2]); sb.append("\t");
		sb.append(t[3]); sb.append("\t");
		//name
		sb.append("numOb="); sb.append(t[4]);
		sb.append(";lg2Tum="); sb.append(t[5]);
		sb.append(";lg2Norm="); sb.append(t[6]);
		if (t.length == 16) {
			sb.append(";genes="); 
			if (t[15].equals(".")==false) sb.append(t[15]);
		}
		sb.append("\t");
		//score
		sb.append(t[5]); sb.append("\t");
		//strand
		if (t[5].startsWith("-")) sb.append("-");
		else sb.append("+");
		return sb.toString();
	}
	
	public boolean check(String[] t) {
		//Sample Chr Start End CR_NumPoints CR_LgTumorMean CR_LgNormalMean CR_LgTNRatio AF_TumorSize AF_TumorMean
		//   0    1    2    3        4            5               6             7              8           9             
		//AF_NormalSize AF_NormalMean AF_TNRatio GATK_Call GATK_LgTumorGeoMean	'Genes'
		//      10             11         12         13             14			   15
		if (checkGT(Double.parseDouble(t[5]), minTLogRatio) == false) return false;
		
		//check that GATK made a call?
		if (passGATK && t[13].equals("0")) return false;
		
		if (checkGT(Double.parseDouble(t[14]), minGatkLgGeoMean) == false) return false;
		
		//normal values present?
		if (t[6].equals("0")==false && t[7].equals("0")==false) {
			if (checkGT(Double.parseDouble(t[7]), minTNLogRatio) == false) return false;
			if (checkLT(Double.parseDouble(t[6]), maxNLogRatio) == false) return false;
		}
		return true;
	}

	private boolean checkGT(double value, double threshold) {
		//was the threshold set by user?
		if (threshold == Double.MAX_VALUE) return true;
		double abs = Math.abs(value);
		return (abs >= threshold); 
	}
	private boolean checkLT(double value, double threshold) {
		//was the threshold set by user?
		if (threshold == Double.MAX_VALUE) return true;
		double abs = Math.abs(value);
		return (abs <= threshold); 
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CopyAnalysisCalledAnnoSegFilter(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File segDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': segDir = new File (args[++i]); break;
					case 'r': resultsDir = new File (args[++i]); break;
					case 'm': minTNLogRatio = Double.parseDouble(args[++i]); break;
					case 't': minTLogRatio = Double.parseDouble(args[++i]); break;
					case 'n': maxNLogRatio = Double.parseDouble(args[++i]); break;
					case 'c': passGATK = true; break;
					case 'b': outputBed = true; break;
					case 'g': minGatkLgGeoMean = Double.parseDouble(args[++i]); break;
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull bed files
		if (segDir == null || segDir.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a bed file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(segDir, ".seg");
		tot[1] = IO.extractFiles(segDir,".seg.gz");
		tot[2] = IO.extractFiles(segDir,".seg.zip");
		segFilesToParse = IO.collapseFileArray(tot);
		if (segFilesToParse == null || segFilesToParse.length ==0 || segFilesToParse[0].canRead() == false) {
			Misc.printExit("\nError: cannot find your xxx.bed(.zip/.gz OK) file(s)!\n");
		}
		if (resultsDir == null) Misc.printExit("\nError: please provide a directory path for saving the results.\n");
		resultsDir.mkdirs();
	}	
	
	private void printSettings() throws IOException {
		IO.pl("Settings:");
		IO.pl("\tSeg Dir:\t"+segFilesToParse[0].getCanonicalPath());
		IO.pl("\tRes Dir:\t"+resultsDir.getCanonicalPath());
		if (minTNLogRatio != Double.MAX_VALUE) IO.pl("\tminTNLogRatio:\t"+minTNLogRatio);
		if (minTLogRatio != Double.MAX_VALUE) IO.pl("\tminTLogRatio:\t"+minTLogRatio);
		if (maxNLogRatio != Double.MAX_VALUE) IO.pl("\tmaxNLogRatio:\t"+maxNLogRatio);
		if (minGatkLgGeoMean != Double.MAX_VALUE) IO.pl("\tminGatkLgGeoMean:\t"+minGatkLgGeoMean);
		if (passGATK) IO.pl("\tGATK Called:\ttrue");
		IO.pl();
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                    CopyAnalysisCalledAnnoSegFilter: June 2026                    **\n" +
				"**************************************************************************************\n" +
				"Filters the xxx.called.anno.seg files from the TNRunner CopyAnalysis workflow. Only\n"+
				"thresholds set below will be used in filtering. For bed output, run the USeq\n"+
				"AnnotateBedWithGenes to annotate the input seg files first if you want to include\n"+
				"intersecting genes info.\n"+

				"\nOptions:\n"+
				"-s Directory containing the xxx.called.anno.seg files.\n"+
				"-r Directory to save the filtered files.\n"+
				"-m Minimum abs TN Log2Rto, only checked if a paired normal was used in the analysis.\n"+
				"-n Maximum abs N Log2Rto, ditto.\n"+
				"-t Minimum abs T Log2Rto.\n"+
				"-c Check that GATK issued a + or - call.\n"+
				"-g Minimum abs GATK geometric mean T Log Rto.\n"+
				"-b Output bed formal.\n"+
				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/CopyAnalysisCalledAnnoSegFilter -s Segs/ \n" +
				"     -r PassingSegs/ -c -g 0.1 -n 0.5 -m 0.1 -t 0.1 \n\n"+ 


				"************************************************************************************\n");

	}
}
