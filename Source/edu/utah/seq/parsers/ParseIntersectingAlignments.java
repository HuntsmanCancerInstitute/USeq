
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.data.*;
import edu.utah.seq.useq.data.RegionScoreText;
import util.gen.*;
import java.util.*;


import util.bio.annotation.*;
import util.bio.seq.Seq;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class ParseIntersectingAlignments {
	
	//fields
	private File alleleBedFile;
	private File alignmentBedFile;
	private int plusIndex;
	private int minusIndex;
	private ArrayList<AlleleWithAlignments> hitsAL = new ArrayList<AlleleWithAlignments>();
	private HashMap<String,RegionScoreText[]> chromAlleles;
	private HashMap<String,File> chromAlignments;
	private double minimumBaseQuality = 13;
	private boolean saveAlleleWithAlignments = true;
	
	//constructors
	public ParseIntersectingAlignments(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		//split alignment file by chromosome, not sorted
		System.out.println("Splitting alignment file by chromosome...");
		File tempDirectory = new File (alignmentBedFile.getParentFile(), alignmentBedFile.getName()+"_tempDirDeleteMe");
		tempDirectory.mkdir();
		chromAlignments = Bed.splitBedFileByChromStrandToTempDir(alignmentBedFile, tempDirectory);
		
		//parse allele file, assumes interbase coordinates, returns sorted
		System.out.println("Loading alleles...");
		chromAlleles = Bed.parseBedFile(alleleBedFile, false);
		
		//find alignments and record results in AlleleWithAlignments array
		System.out.println("Fetching intersections...");
		findIntersectingAlignments();
		
		//print results to screen
		System.out.println("Saving spreadsheet results...");
		printResults();
		
		//save array?
		if (saveAlleleWithAlignments){
			System.out.println("Saving serialized alleles...");
			File f = new File(alleleBedFile.getParentFile(), Misc.removeExtension(alleleBedFile.getName())+"_PIA.alleles");
			AlleleWithAlignments[] alleles = new AlleleWithAlignments[hitsAL.size()];
			hitsAL.toArray(alleles);
			IO.saveObject(f, alleles);
		}
		
		//cleanup
		IO.deleteDirectory(tempDirectory);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public void printResults(){
		try {
			File spreadSheetFile = new File(alleleBedFile.getParentFile(), Misc.removeExtension(alleleBedFile.getName())+"_PIA.xls");
			PrintWriter outSS = new PrintWriter ( new FileWriter(spreadSheetFile));
			File bedFile = new File(alleleBedFile.getParentFile(), Misc.removeExtension(alleleBedFile.getName())+"_PIA.bed");
			PrintWriter outBed = new PrintWriter ( new FileWriter(bedFile));
			outSS.println("#Chromosome\tStart\tStop\tName\tConScore\tStrand\t#ParsedAlignments\t#UniqueAlignments\tFracG\tFracA\tFracT\tFracC\tFracN\tParsedBPs\tParsedBPScores");
			for (int i=0; i< hitsAL.size(); i++){
				AlleleWithAlignments al = hitsAL.get(i);
				outSS.println(al);
				Bed[] b = al.alignments;
				outBed.println("#\t"+al.toStringBed());
				if (b != null && b.length !=0){
					for (int j=0; j< b.length; j++) outBed.println(b[j]);
				}
			}
			outBed.close();
			outSS.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void findIntersectingAlignments(){
		//find intersecting alignments
		System.out.println("Finding overlapping alignments...");
		Iterator<String> it = chromAlleles.keySet().iterator();
		while (it.hasNext()){
			//reset indexes
			plusIndex = 0;
			minusIndex = 0;
			
			//get chrom alleles
			String chromStrand = it.next();
			String chrom = chromStrand.substring(0, chromStrand.length()-1);
			String strand = chromStrand.substring(chromStrand.length()-1);
			RegionScoreText[] alleles = chromAlleles.get(chromStrand);
			
			//get sorted chrom alignments
			File plusFile = chromAlignments.get(chrom+"+");
			RegionScoreText[] plus = null;;
			if (plusFile != null) plus = RegionScoreText.oldLoadBinary_DEPRECIATED(plusFile, false);
			File minusFile = chromAlignments.get(chrom+"-");
			RegionScoreText[] minus = null;
			if (minusFile != null) minus = RegionScoreText.oldLoadBinary_DEPRECIATED(minusFile, false);

			
			//for each allele
			for (int i=0; i< alleles.length; i++){
				ArrayList<Double> baseScores = new ArrayList<Double>();
				ArrayList<String> parsedBPs = new ArrayList<String>();
				ArrayList<Bed> hits = new ArrayList<Bed>();
				
				ArrayList<RegionScoreText> plusIntAligns = fetchIntersectingAlignments(true, plus, alleles[i]);
				if (plusIntAligns != null && plusIntAligns.size()!=0) {
					for (int j=0; j<plusIntAligns.size(); j++){
						RegionScoreText alignment = plusIntAligns.get(j);
						String[] seqQual = fetchBPs(alignment, alleles[i], true);
						//check quality scores
						if (seqQual[1].equals("") == false){
							int[] vals = Seq.convertScores(seqQual[1]);
							double ave = Num.averageIntArray(vals);
							if (ave < minimumBaseQuality) continue;
							baseScores.add(new Double (ave));
						}

						//add parsed bps
						parsedBPs.add(seqQual[0]);

						//add bed
						hits.add(new Bed(chrom, '+', alignment));
					}
				}
				ArrayList<RegionScoreText> minusIntAligns = fetchIntersectingAlignments(false, minus, alleles[i]);
				if (minusIntAligns != null && minusIntAligns.size()!=0) {
					for (int j=0; j<minusIntAligns.size(); j++){
						RegionScoreText alignment = minusIntAligns.get(j);
						String[] seqQual = fetchBPs(alignment, alleles[i], false);
						//check quality scores
						if (seqQual[1].equals("") == false){
							int[] vals = Seq.convertScores(seqQual[1]);
							double ave = Num.averageIntArray(vals);
							if (ave < minimumBaseQuality) continue;
							baseScores.add(new Double (ave));
						}

						//add parsed bps
						parsedBPs.add(seqQual[0]);
						//add bed
						hits.add(new Bed(chrom, '-', alignment));
					}
				}
				//create hit public AlleleWithAlignments (RegionScoreText allele, double[] baseScores, String[] parsedBPs, String[] alignments){
				
				Bed[] bed = new Bed[hits.size()];
				hits.toArray(bed);
				hitsAL.add(new AlleleWithAlignments(chrom, strand, alleles[i], Num.arrayListOfDoubleToArray(baseScores), Misc.stringArrayListToStringArray(parsedBPs), bed));
			}
		}

	}
	
	/**Returns a String[sequence,qualityScores], assumes they intersect. QualityScores may be "" if none present.*/
	public static String[] fetchBPs (RegionScoreText read, RegionScoreText allele, boolean readPlusStrand){
		int readLength = read.getLength();
		String readName = read.getText();
		int diffStart = allele.getStart()- read.getStart();
		int diffEnd = diffStart+allele.getLength();
		if (diffStart < 0) {
			diffEnd = allele.getLength()+diffStart;
			diffStart = 0;
		}
		else if (diffEnd > readLength) diffEnd = readLength;
		String qual ="";
		int firstSemiIndex = readName.indexOf(";");
		String[] seqQual;
		if (firstSemiIndex != -1){
			seqQual = new String[]{readName.substring(0, firstSemiIndex), readName.substring(firstSemiIndex+1)};
		}
		else seqQual = new String[]{readName};
		if (readPlusStrand == false){
			seqQual[0] = Seq.reverseComplementDNA(seqQual[0]);
			if (seqQual.length>1) {
				seqQual[1] = new StringBuilder(seqQual[1]).reverse().toString();
			}
		}
		int lenSeq = seqQual[0].length();
		if (diffStart < 0 || diffStart> lenSeq) diffStart =0;
		if (diffEnd> lenSeq) diffEnd = lenSeq;
		String seq = seqQual[0].substring(diffStart, diffEnd);
		if (seqQual.length > 1){
			qual = seqQual[1].substring(diffStart, diffEnd);
		}
		return new String[]{seq,qual};
	}

	public ArrayList<RegionScoreText> fetchIntersectingAlignments(boolean plusStrand, RegionScoreText[] alignments, RegionScoreText allele){
		//any alignments?
		if (alignments == null) return null;
		
		//set starting index
		int startingIndex = plusIndex;
		if (plusStrand == false) startingIndex = minusIndex;
		
		//System.out.println("PlusIndex "+plusIndex);
		//System.out.println("MinusIndex "+minusIndex);
		
		//initialize objects
		ArrayList<RegionScoreText> ints = new ArrayList<RegionScoreText>();
		int start = allele.getStart();
		int end = allele.getStop();
		boolean foundIntersectingAlignments = false;
		
		//for each alignment look for intersection, these are sorted so use indexing
		int firstIntersectionIndex = startingIndex;
		for (int i=startingIndex; i< alignments.length; i++){
			//does it intersect?
			if (alignments[i].intersects(start, end)) {
				ints.add(alignments[i]);
				if (foundIntersectingAlignments == false){
					foundIntersectingAlignments = true;
					firstIntersectionIndex = i;
				}
			}
			//is it past
			else if (end<= alignments[i].getStart()){
				//System.out.println("Past "+alignments[i]); 
				
				break;
			}
		}
		
		//set index
		if (foundIntersectingAlignments){
			if (plusStrand) plusIndex = firstIntersectionIndex;
			else minusIndex = firstIntersectionIndex;
		}
		return ints;
	}
	

	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ParseIntersectingAlignments(args);
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
					case 's': alleleBedFile = new File(args[++i]); break;
					case 'a': alignmentBedFile = new File(args[++i]); break;
					case 'm': minimumBaseQuality = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (alleleBedFile == null || alignmentBedFile == null) Misc.printErrAndExit("\nCouldn't find one or both of your xxx.bed files?\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        ParseIntersectingAlignments: June 2010                    **\n" +
				"**************************************************************************************\n" +
				"Parses bed alignment files for intersecting reads provided another bed file of alleles.\n" +

				"\nOptions:\n"+
				"-s Full path file text for your SNP allele five column bed file (tab delimited chr,\n" +
				"      start,stop,text,score,strand)\n" +
				"-a Full path file text for your alignment bed file from the NovoalignParser.\n" +
				"-m Minimum base quality, defaults to 13\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/ParseIntersectingAlignments \n" +
				"     -s /LympAlleles/ex1.bed -a /SeqData/lymphAlignments.bed -m 13\n\n" +

		"**************************************************************************************\n");

	}	

}
