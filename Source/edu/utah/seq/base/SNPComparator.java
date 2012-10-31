package edu.utah.seq.base;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.parsers.VCFLookUp;
import edu.utah.seq.parsers.VCFParser;
import edu.utah.seq.parsers.VCFRecord;
import edu.utah.seq.useq.data.RegionScoreText;

import util.bio.annotation.Bed;
import util.gen.*;

/**Compares variant lists
 * @author Nix
 * */
public class SNPComparator {

	//user fields
	private File sequencingVCFFile;
	private File arraySnpBedFile;
	private float minimumArrayScore = 0.9f;
	private double minimumVCFScore = 0.01;
	private int vcfScoreIndex = 7;
	
	private HashMap<String, VCFLookUp> chromVCFRecords;
	private HashMap<String,RegionScoreText[]> chromRegions;
	private int numberArraySnpsFailingMinimumScore = 0;
	private int numberVCFsFailingMinimumScore = 0;
	private int totalNumberArraySnps = 0;
	private boolean excludeFailingCalls = true;



	//constructor
	public SNPComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}

	//methods
	
	
	private void doWork(){
		//parse vcf file
		VCFParser vcfParser = new VCFParser(sequencingVCFFile, excludeFailingCalls);
		chromVCFRecords = vcfParser.getChromVCFRecords();
		System.out.println(vcfParser.getNumberVCFRecords() + "\t# Parsed VCF records");
		
		//parse bed file of array snp calls
		chromRegions = Bed.parseBedFile(arraySnpBedFile, true);
		
		//for each chromosome of snp calls (the truth), compare to the seq calls
		for (String chr: chromRegions.keySet()){
			
			if (chromVCFRecords.containsKey(chr)){
				System.out.println(chr);
				compare (chromRegions.get(chr), chromVCFRecords.get(chr));
			}
		}
		
	}

	private void compare(RegionScoreText[] regionScoreTexts, VCFLookUp vcfLookUp) {
		totalNumberArraySnps += regionScoreTexts.length;
		for (RegionScoreText snp: regionScoreTexts){
			//check score
			if (snp.getScore()< minimumArrayScore){
				numberArraySnpsFailingMinimumScore++;
				continue;
			}
			
			//an sequence calls? should be just one or null
			VCFRecord[] vcf = vcfLookUp.fetchVCFRecords(snp.getStart(), snp.getStop());
			if (vcf == null) continue;
			
			//check score
			double score = Double.parseDouble(vcf[0].getSampleField(vcfScoreIndex));
			if (score < minimumVCFScore){
				numberVCFsFailingMinimumScore++;
				continue;
			}
			
			//
			//compare bases
		}
		
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SNPComparator(args);
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
					case 'v': sequencingVCFFile = new File(args[++i]); break;
					case 'b': arraySnpBedFile = new File(args[++i]); break;
					case 's': minimumArrayScore = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//checkfiles
		if (sequencingVCFFile == null || sequencingVCFFile.canRead() == false) Misc.printExit("\nError: cannot find your sequencing variant file!\n");
		if (arraySnpBedFile == null || arraySnpBedFile.canRead() == false) Misc.printExit("\nError: cannot find your array variant bed file!\n");
		

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Base Classifier : Oct 2012                            **\n" +
				"**************************************************************************************\n" +
				"Beta.\n\n" +

				"Options:\n"+
				"-v Full path to a sorted bam aligment file or directory containing such. Multiple files\n" +
				"       are merged. xxx.bai indexes required.\n"+
				"-v Full path to a bed file containing variants to classify.\n"+
				"-n Length of the N mer, defaults to 5.\n" +
				"-q Minimum base quality score, defaults to 20. Only N-mers where all bases pass the\n"+
				"       threshold are scored.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/BamIntensityParser -f /Data/BamFiles/\n" +
				"       -n 7 -q 30\n\n"+

		"**************************************************************************************\n");

	}
}
