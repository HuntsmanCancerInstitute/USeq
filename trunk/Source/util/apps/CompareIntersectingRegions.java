package util.apps;
import java.io.*;

import util.bio.annotation.Bed;
import util.gen.*;
import edu.utah.seq.useq.data.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**Takes a master bed file and intersects it with one or more noting shared regions with the master.*/
public class CompareIntersectingRegions {

	//fields
	private File masterBedFile;
	private Bed[] masterRegions;
	private File[] testRegionsFiles;
	private HashMap<String,RegionScoreText[]> testRegions;
	private int maximumGap = 0;
	private StringBuilder[] results;
	private int[] hits;
	private StringBuilder header = new StringBuilder();
	

	public CompareIntersectingRegions(String[] args) {
		//process args
		processArgs(args);

		//load regions
		System.out.println("Loading master regions...");
		masterRegions = Bed.parseFile(masterBedFile, 0,0);
		Arrays.sort(masterRegions);
		
		//build StringBuilder[] to hold results
		results = new StringBuilder[masterRegions.length];
		hits = new int[masterRegions.length];
		for (int i=0; i< masterRegions.length; i++) {
			results[i] = new StringBuilder();
			results[i].append(masterRegions[i].toStringNoStrand());
			//increase size if max gap specified
			if (maximumGap > 0) {
				int start = masterRegions[i].getStart();
				start = start - maximumGap;
				if (start < 0) start = 0;
				masterRegions[i].setStart(start);
				masterRegions[i].setStop(masterRegions[i].getStop()+ maximumGap);
			}
			//shrink size
			else if (maximumGap <0){
				int start = masterRegions[i].getStart();
				start = start - maximumGap;
				int middle = masterRegions[i].getMiddle();
				// too big a subtraction? if so make it one base in length around center position
				if (start > middle){
					masterRegions[i].setStart(middle);
					masterRegions[i].setStop(middle+1);
				}
				else {
					masterRegions[i].setStart(start);
					masterRegions[i].setStop(masterRegions[i].getStop() + maximumGap);
				}
			}
		}
		header.append("Chr\tStart\tStop\tName\tScore");
		
		System.out.println("\nIntersecting with:");
		//for each file to split
		for (File test: testRegionsFiles){
			
			System.out.println("\t"+test.getName());
			//parse regions and split by chromosome
			testRegions = Bed.parseBedFile(test, true);
			//add name to header
			header.append("\t"+Misc.removeExtension(test.getName()));

			//for each master regions
			String oldChrom = "";
			RegionScoreText[] chromTestRegions = null;
			for (int i=0; i< masterRegions.length; i++){
				
				//load chrom regions?
				if (masterRegions[i].getChromosome().equals(oldChrom) == false){
					//master has switched
					chromTestRegions = testRegions.get(masterRegions[i].getChromosome());
					oldChrom = masterRegions[i].getChromosome();
				}
				//see if it intersects
				RegionScoreText hit = null;
				if (chromTestRegions != null) {
					hit = intersect(masterRegions[i], chromTestRegions);
				}
				results[i].append("\t");
				if (hit != null) {
					results[i].append(hit.toStringUnderscore());
					hits[i]++;
				}
			}
		}	
		
		//append hits
		header.append("\tHits");
		for (int i=0; i< results.length; i++){
			results[i].append("\t");
			results[i].append(hits[i]);
		}
		
		//print out results
		printResults();
		

			
	}

	private void printResults() {
		File outFile = new File (masterBedFile.getParentFile(), Misc.removeExtension(masterBedFile.getName())+"_CIR.xls");
		System.out.println("\nSaving results to "+outFile+"\n");
		try {
			PrintWriter out = new PrintWriter( new FileWriter(outFile));
			out.println(header);
			for (StringBuilder sb: results) out.println(sb);
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(0);
		}
		
	}

	private RegionScoreText intersect(Bed bed, RegionScoreText[] chromTestRegions) {
		for (RegionScoreText test: chromTestRegions){
			if (bed.intersects(test.getStart(), test.getStop())) return test;
		}
		return null;
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CompareIntersectingRegions(args);
	}		
	

	/**This method will process each argument and assign new varibles*/
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
					case 'm': masterBedFile = new File (args[++i]); break;
					case 't': testRegionsFiles = IO.extractFiles(new File(args[++i])); break;
					case 'g': maximumGap = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (masterBedFile == null || masterBedFile.canRead() == false) Misc.printErrAndExit("\nCannot find your master region file! "+masterBedFile+"\n");
		if (testRegionsFiles == null || testRegionsFiles[0].canRead() == false) Misc.printErrAndExit("\nCannot find your test region file(s) to compare! "+testRegionsFiles[0]+"\n");
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Compare Intersecting Regions: Nov 2012                    **\n" +
				"**************************************************************************************\n" +
				"Compares test region file(s) against a master set of regions for intersection.\n" +
				"Reports the results as columns relative to the master. Assumes interbase coordinates.\n" +

				"\nOptions:\n"+
				"-m Full path for the master bed file (tab delim: chr start stop ...).\n"+
				"-t Full path to the test bed file to intersect or directory of files.\n"+
				"-g Maximum bp gap allowed for scoring an intersection, defaults to 0 bp. Negative gaps\n" +
				"     force overlaps, positive gaps allow non intersecting bases between regions.\n"+

				"\nExample: java -Xmx4G -jar pathTo/Apps/CompareIntersectingRegions -g 1000\n" +
				"        -m /All/mergedRegions.bed.gz -t /IndividualERs/\n\n" +

		"************************************************************************************\n");
	}

}
