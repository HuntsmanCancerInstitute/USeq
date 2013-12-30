package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.mes.*;
import util.bio.annotation.ExonIntron;
import util.bio.annotation.ExportIntergenicRegions;
import util.bio.parsers.MultiFastaParser;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.bio.seq.Seq;
import util.gen.*;

/**Annotates vcf records with Max Ent Scan predictions for the gain or loss of a 5' or 3' splice junction.
 * @author Nix
 * */
public class VCFSpliceAnnotator {

	//user fields
	private File[] vcfFiles;
	private File chromDirectory;
	private File spliceModelDirectory;
	private File refSeqFile;
	private double min5Threshold = 5;
	private double min5DeltaThreshold = 5;
	private double min3Threshold = 5;
	private double min3DeltaThreshold = 5;
	private boolean scoreNovelIntronJunctions = true;
	private boolean scoreNovelExonJunctions = true;
	private boolean scoreNovelSpliceJunctionsInSplice = true;
	
	//internal fields
	private HashMap<String, File> chromFile;
	private MaxEntScanScore5 score5;
	private MaxEntScanScore3 score3;
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private String workingChromosomeName = "";
	private String workingSequence = null;
	private UCSCGeneLine[] workingTranscripts;
	private Histogram exonicHistogram5;
	private Histogram exonicHistogram3;
	private Histogram intronicHistogram5;
	private Histogram intronicHistogram3;
	private Histogram spliceHistogram5;
	private Histogram spliceHistogram3;
	private boolean printHistograms = true;
	
	private int numTranscripts = 0; //
	private HashSet<String> intersectingTranscriptNames = new HashSet<String>();
	private int numVariantsScanned = 0;
	private int numVariantsIntersectingTranscripts = 0;
	private int numVariantsIntersectingExons = 0;
	private int numVariantsIntersectingIntrons = 0;
	private int numVariantsIntersectingSJs = 0;
	private int numExonSJsGained = 0;
	private int numIntronSJsGained = 0;
	private int numSpliceJunctionSJsGained = 0;
	private int numSJsLost = 0;
	
	//constructor
	public VCFSpliceAnnotator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//fetch seqs by chrom
		chromFile = Seq.fetchChromosomeFastaFileHashMap(chromDirectory);

		//start up mes
		score5 = new MaxEntScanScore5(spliceModelDirectory);
		score3 = new MaxEntScanScore3(spliceModelDirectory);

		//load transcripts
		System.out.println("Loading transcripts...");
		loadTranscripts();
		
		//histograms of novel exonic splices
		System.out.println("Loading null score histograms...");
		loadNullScoreHistograms();

		//for each vcf file
		System.out.println("Processing...");
		for (int i=0; i< vcfFiles.length; i++){
			System.out.println("\t"+vcfFiles[i]);
			annotateVCFWithSplices(vcfFiles[i]);
		}
		
		//summaries
		System.out.println();
		printThresholds();
		printSummary();
		
		
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone "+Math.round(diffTime)+" seconds!\n");
	}
	
	//methods
	
	private void printThresholds(){
		StringBuilder sb = new StringBuilder();
		sb.append("Threholds:\n");
		sb.append(min5Threshold +"\tMinimum 5' splice junction threshold for scoring the presence of a junction\n");
		sb.append(min3Threshold +"\tMinimum 3' splice junction threshold for scoring the presence of a junction\n");
		sb.append(min5DeltaThreshold +"\tMinimum score difference for loss or gain of a 5' splice junction\n");
		sb.append(min3DeltaThreshold +"\tMinimum score difference for loss or gain of a 3' splice junction\n");
		System.out.println(sb);
	}
	
	private void printSummary(){
		StringBuilder sb = new StringBuilder();
		sb.append("Summary stats:\n");
		sb.append(numTranscripts+"\tTranscripts\n"); 
		sb.append(intersectingTranscriptNames.size()+"\tTranscripts intersecting variants\n"); 
		sb.append(numVariantsScanned+"\tVariants (including alternates)\n");
		sb.append(numVariantsIntersectingTranscripts+"\tVariants intersecting transcripts\n");
		sb.append(numVariantsIntersectingExons+"\tVariants intersecting exons (non splice junction)\n"); 
		sb.append(numVariantsIntersectingIntrons+"\tVariants intersecting introns (non splice junction)\n"); 
		sb.append(numVariantsIntersectingSJs+"\tVariants intersecting splice junctions\n"); 
		sb.append(numExonSJsGained+"\tExonic splice junctions gained\n"); 
		sb.append(numIntronSJsGained+"\tIntronic splice junctions gained\n"); 
		sb.append(numSpliceJunctionSJsGained+"\tSplice junctions with novel splice junction gained\n"); 
		sb.append(numSJsLost+"\tSplice junctions damaged\n");
		System.out.println(sb);
	}

	public void loadTranscripts(){
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(refSeqFile, 0);
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your transcript's coordinates are reversed. Check that each start is less than the stop.\n");
		chromGenes = reader.getChromSpecificGeneLines();
		numTranscripts = reader.getGeneLines().length;
	}

	/**Adds splice junction information to each record. */
	public void annotateVCFWithSplices(File vcfFile) {
		File modFile = new File (vcfFile.getParentFile(), Misc.removeExtension(vcfFile.getName())+".splice.txt");
		BufferedReader in = null;
		PrintWriter out = null;
		VCFParser parser = new VCFParser();
		try {
			in  = IO.fetchBufferedReader(vcfFile);
			out = new PrintWriter ( new FileWriter (modFile));

			//add ##INFO line and find "#CHROM" line 
			loadAndModifyHeader(in, out);

			//For each record
			String line;
			while ((line=in.readLine()) != null){
				//could delete for speed
				if (line.length()== 0) continue;
				if (line.startsWith("#")) {
					System.out.println(line);
					continue;
				}

				//parse record adding chr to chrom name if not present
				VCFRecord vcf = new VCFRecord(line, parser, true, true);
				vcf.appendChr();
				vcf.correctChrMTs();
				
				//load new sequence and transcripts?
				if (vcf.getChromosome().equals(workingChromosomeName) == false) {
					//load new working seq
					loadChromosomeData(vcf.getChromosome());
					System.out.print(workingChromosomeName+" ");
				}
				
				//check if any transcripts were found, otherwise skip, might want to write out this vcf record somewhere?
				if (workingTranscripts == null) continue;

				//for each alternate allele, some vcf records have several
				String[] alts = vcf.getAlternate();
				for (int i=0; i< alts.length; i++){	
					numVariantsScanned++;
					vcf.setAlternate(new String[]{alts[i]});
					
					//fetch intersecting transcripts
					UCSCGeneLine[] trans = fetchIntersectingTranscripts(vcf);
					if (trans == null) continue;
					numVariantsIntersectingTranscripts++;
					
					//for each transcript score effect
					for (UCSCGeneLine l: trans){
						String effect = scoreVariantTranscript(vcf, l);
						if (effect != null) {
							out.println("Hit to gene -> "+l.toString());
							out.println("By variant  -> "+vcf.toString());
							out.println(effect);
						}
					}
				}
				vcf.setAlternate(alts);
			}
			//last
			System.out.println();
			out.close();
			in.close();
		}catch (Exception e) {
			System.err.println("Error -> "+e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

	/**Be sure to first mask the known splices in the workingSequence*/
	private void scoreAllExonsOrIntronsForSplices(boolean plusStrand, boolean scoreExons) {
		//flatten exons
		boolean[] bases = new boolean[workingSequence.length()];
		Arrays.fill(bases, true);
		//for each transcript
		for (UCSCGeneLine g: workingTranscripts){
			//correct strand?
			boolean isPlus = g.getStrand().equals("+");
			if ( isPlus != plusStrand) continue;
			//get exons
			ExonIntron[] regions;
			if (scoreExons) regions = g.getExons();
			else regions = g.getIntrons();
			if (regions == null) continue;
			//for each region
			for (ExonIntron e: regions){
				int start = e.getStart();
				int end = e.getEnd();
				//mask bases
				for (int i=start; i< end; i++){
					bases[i] = false;
				}
			}
		}
		//fetch blocks of false, 
		int[][] startStop = ExportIntergenicRegions.fetchFalseBlocks(bases,0,0);
		//for each exon block
		for (int[] ss: startStop){
			//fetch sequence with masked splices
			String seq = workingSequence.substring(ss[0], ss[1]+1);
			if (plusStrand == false) seq = Seq.reverseComplementDNA(seq);
			double[] scores5 = score5.scanSequence(seq);
			double[] scores3 = score3.scanSequence(seq);
			//add scores to histograms
			if (scoreExons){
				exonicHistogram5.countAll(scores5);
				exonicHistogram3.countAll(scores3);
			}
			else {
				intronicHistogram5.countAll(scores5);
				intronicHistogram3.countAll(scores3);
			}
		}
	}


	

	/**Deletion 20     2 .         TC      T      .   PASS  DP=100"
	 * Insertion 20     2 .         TC      TCA    .   PASS  DP=100  
	 */
	private String scanIntronsForNovelSpliceJunctions (VCFRecord vcf, UCSCGeneLine transcript, boolean isPlusStrand){
		ExonIntron[] introns = transcript.getIntrons();
		String results = "noIntersection";
		StringBuilder sb = new StringBuilder();
		if (introns != null) {
			int leftAdd = 6;
			int rightSub = 20;
			if (isPlusStrand == false){
				leftAdd = 20;
				rightSub = 6;
			}
			int varStart = vcf.getPosition();
			int varEnd = varStart+ vcf.getReference().length();
			for (int i=0; i< introns.length; i++){
				//mod coordinates of intron
				int intronStart = introns[i].getStart()+leftAdd;
				int intronEnd = introns[i].getEnd()-rightSub;
				if (intronEnd <= intronStart) continue;
				//if a hit then score 
				if (intersects(varStart, varEnd, intronStart, intronEnd) == false) continue;
				//ok its a hit, set to noChange
				numVariantsIntersectingIntrons++;
				if (results.equals("noIntersection")) results = "noChange";
				//scan for novel, reference and alternate			
				String[] fiveSeqs = fetchNMerSeqs(vcf, 9, isPlusStrand == false);
				String fiveRes = scoreNewJunction(fiveSeqs, true, false);
				String[] threeSeqs = fetchNMerSeqs(vcf, 23, isPlusStrand == false);
				String threeRes = scoreNewJunction(threeSeqs, false, false);
				if (fiveRes != null || threeRes != null){
					StringBuilder alerts = new StringBuilder();
					alerts.append ("Novel junction in intron "+introns[i].getStartStopString()+"\n");
					if (fiveRes != null) {
						alerts.append (fiveRes+"\n");
						numIntronSJsGained++;
					}
					if (threeRes != null) {
						alerts.append (threeRes+"\n");
						numIntronSJsGained++;
					}
					//save results 
					sb.append(alerts.toString()+"\n");
				}
				//can only have one intersection if snp
				if (vcf.isSNP()) break;
			}
			if (sb.length() != 0) results = sb.toString();
		}
		return results;
	}



	private String scanExonsForNovelSpliceJunctions (VCFRecord vcf, UCSCGeneLine transcript, boolean isPlusStrand){
		ExonIntron[] exons = transcript.getExons();
		String results = "noIntersection";
		StringBuilder sb = new StringBuilder();
		if (exons.length != 0) {
			int varStart = vcf.getPosition();
			int varEnd = varStart+ vcf.getReference().length();
			int lastExonIndex = exons.length-1;
			for (int i=0; i< exons.length; i++){
				int start = exons[i].getStart();
				int end = exons[i].getEnd();
				//if more than one exon, need to modify ends to avoid splice junctions
				if (exons.length !=1){
					if (i!=0) start+=3; 
					if (i!= lastExonIndex) end-=3; 
					if (end<= start) continue;
				}
				//if a hit then score 
				if (intersects(varStart, varEnd, start, end) == false) continue;
				numVariantsIntersectingExons++;
				if (results.equals("noIntersection")) results = "noChange";
				//scan for novel, reference and alternate			
				String[] fiveSeqs = fetchNMerSeqs(vcf, 9, isPlusStrand == false);
				String fiveRes = scoreNewJunction(fiveSeqs, true, true);
				String[] threeSeqs = fetchNMerSeqs(vcf, 23, isPlusStrand == false);
				String threeRes = scoreNewJunction(threeSeqs, false, true);
				if (fiveRes != null || threeRes != null){
					StringBuilder alerts = new StringBuilder();
					alerts.append ("Novel junction in exon "+exons[i].getStartStopString()+"\n");
					if (fiveRes != null) {
						alerts.append (fiveRes+"\n");
						numExonSJsGained++;
					}
					if (threeRes != null) {
						alerts.append (threeRes+"\n");
						numExonSJsGained++;
					}
					//save results 
					sb.append(alerts.toString()+"\n");
				}
				//can only have one intersection if snp
				if (vcf.isSNP()) break;
			}
			if (sb.length() != 0) results = sb.toString();
			
		}
		return results;
	}
	
	private boolean intersects(int varStart, int varEnd, int annoStart, int annoEnd){
		if (varStart >= annoEnd) return false;
		if (varEnd <= annoStart) return false;
		return true;
	}
	
	private String scanKnownSplices (VCFRecord vcf, UCSCGeneLine transcript, boolean isPlusStrand){
		ExonIntron[] introns = transcript.getIntrons();
		int varStart = vcf.getPosition();
		int varEnd = varStart+ vcf.getReference().length();
		if (introns!= null) {
			//find junction
			int junc = 0;
			int startJunc = 0;
			int endJunc = 0;
			boolean found = false;
			StringBuilder sb = new StringBuilder();
			for (int i=0; i< introns.length; i++){
				//plus strand 
				if (isPlusStrand){
					//5'?
					startJunc = introns[i].getStart()-3;     
					endJunc = introns[i].getStart()+6;
					if (intersects(varStart, varEnd, startJunc, endJunc)) {
						junc = introns[i].getStart();
						String spliceResult = scoreLossOfKnownJunction(junc, vcf, true, isPlusStrand);
						if (spliceResult != null) {
							sb.append(spliceResult);
							numSJsLost++;
						}
						found = true;
					}
					//3'?
					else {
						startJunc = introns[i].getEnd()-20;
						endJunc = introns[i].getEnd()+3;
						if (intersects(varStart, varEnd, startJunc, endJunc)) {
							junc = introns[i].getEnd();
							String spliceResult = scoreLossOfKnownJunction(junc, vcf, false, isPlusStrand);
							if (spliceResult != null) {
								sb.append(spliceResult);
								numSJsLost++;
							}
							found = true;
						}
					}
				}
				//minus strand
				else {
					//3'?
					startJunc = introns[i].getStart()-3;
					endJunc = introns[i].getStart()+20;
					if (intersects(varStart, varEnd, startJunc, endJunc)) {
						junc = introns[i].getStart();
						String spliceResult = scoreLossOfKnownJunction(junc, vcf, false, isPlusStrand);
						if (spliceResult != null) {
							sb.append(spliceResult);
							numSJsLost++;
						}
						found = true;
					}
					//5'?
					else {
						startJunc = introns[i].getEnd()-6;
						endJunc = introns[i].getEnd()+3;
						if (intersects(varStart, varEnd, startJunc, endJunc)) {
							junc = introns[i].getEnd();
							String spliceResult = scoreLossOfKnownJunction(junc, vcf, true, isPlusStrand);
							if (spliceResult != null) {
								sb.append(spliceResult);
								numSJsLost++;
							}
							found = true;
						}
					}
				}
			}

			//did it hit a splice junction?
			if (found) numVariantsIntersectingSJs++;
			else  return "noIntersection";

			//now scan for novel new splice in splice site
			if (scoreNovelSpliceJunctionsInSplice){
				String[] fiveSeqs = fetchNMerSeqs(vcf, 9, isPlusStrand == false);
				String fiveRes = scoreNewJunction(fiveSeqs, true, false);
				String[] threeSeqs = fetchNMerSeqs(vcf, 23, isPlusStrand == false);
				String threeRes = scoreNewJunction(threeSeqs, false, false);
				if (fiveRes != null || threeRes != null){
					sb.append ("Novel junction found near splice junction ("+junc+")\n");
					if (fiveRes != null) {
						sb.append (fiveRes+"\n");
						numSpliceJunctionSJsGained++;
					}
					if (threeRes != null) {
						sb.append (threeRes+"\n");
						numSpliceJunctionSJsGained++;
					}
				}
			}

			//any results?
			if (sb.length() == 0) return "noChange";
			return sb.toString();
		}

		return "noIntersection";

	}


	/**Returns null if no actionable info.  Otherwise with message.*/
	private String scoreVariantTranscript(VCFRecord vcf, UCSCGeneLine transcript) {
		boolean isPlusStrand = transcript.getStrand().equals("+");
		boolean isSNP = vcf.isSNP();

		//can hit multiple annotation classes so must scan all
		StringBuilder sb = new StringBuilder();
		//intronic and away from splice junctions
		if (scoreNovelIntronJunctions){
//System.out.println("Scanning Introns!");
			String intronResults = scanIntronsForNovelSpliceJunctions(vcf, transcript, isPlusStrand);
//System.out.println(intronResults);
			if (isSNP){
				if (intronResults.equals("noChange")) return null;
				if (intronResults.equals("noIntersection") == false) return intronResults;
			}
			//if (intronResults.length()>20 && vcf.isDeletion()) System.out.println("intronINDELhit \n"+vcf.getOriginalRecord()+"\n"+intronResults);
			if (intronResults.equals("noChange") == false && intronResults.equals("noIntersection") == false) sb.append(intronResults);
		}
		//exonic and away from splice junction
		if (scoreNovelExonJunctions){
//System.out.println("Scanning Exons!");			
			String exonResults = scanExonsForNovelSpliceJunctions(vcf, transcript, isPlusStrand);
//System.out.println(exonResults);
			if (isSNP){
				if (exonResults.equals("noChange")) return null;
				if (exonResults.equals("noIntersection") == false) return exonResults;
			}
			//if (exonResults.length()>20 && vcf.isDeletion()) System.out.println("exonINDELhit \n"+vcf.getOriginalRecord()+"\n"+exonResults);
			if (exonResults.equals("noChange") == false && exonResults.equals("noIntersection") == false) sb.append(exonResults);
		}
		//scan hits to splices, gain and loss
//System.out.println("Scanning Splices!");		
		String spliceResults = scanKnownSplices(vcf, transcript, isPlusStrand);
//System.out.println(spliceResults);
		//if (sb.length() == 0 && spliceResults.equals("noIntersection")) Misc.printErrAndExit("Hmm variant didn't land on an annotation, something is wrong!\n"+transcript.toStringAll()+"\n"+ vcf.toString());
		if (spliceResults.equals("noChange") == false && spliceResults.equals("noIntersection") == false) sb.append(spliceResults);		
		
		if (sb.length() == 0) return null;
		return sb.toString();
	}
	
	/**Returns three sequences, the reference splice junction, the modified splice junction where the position has been fixed or follows the original.
	 * {refSeq, fixedSeq, shiftedSeq} */
	private String[] fetchModifiedSpliceSequences(VCFRecord vcf, int junctionPosition, boolean isPlusStrand, boolean fetch5Prime){
		//set cut points
		int leftSub;
		int rightAdd;
		if (isPlusStrand){
			if (fetch5Prime){
				leftSub = 3;
				rightAdd = 6;
			}
			//3'
			else {
				leftSub = 20;
				rightAdd = 3;
			}
		}
		//minus strand
		else {
			if (fetch5Prime){
				leftSub = 6;
				rightAdd = 3;
			}
			//3'
			else {
				leftSub = 3;
				rightAdd = 20;
			}
		}
		
		int position = vcf.getPosition();
		String ref = vcf.getReference();
		String alt = vcf.getAlternate()[0];
		String modChromSeq = workingSequence.substring(0,position) + alt + workingSequence.substring(position + ref.length() );
		
		//if indel precedes or is at the junction, must shift junction location
		int modJunction = junctionPosition;
		if (position <= junctionPosition){
			int diff = alt.length() - ref.length();
			modJunction += diff;
		}
		
		String refSeq = workingSequence.substring(junctionPosition - leftSub, junctionPosition + rightAdd);
		String fixedSeq = modChromSeq.substring(junctionPosition - leftSub, junctionPosition + rightAdd);
		String shiftedSeq = modChromSeq.substring(modJunction - leftSub, modJunction+ rightAdd);
		
		if (isPlusStrand == false){
//System.out.println("NoRCRS  \t"+refSeq);
			refSeq = Seq.reverseComplementDNA(refSeq);
			fixedSeq = Seq.reverseComplementDNA(fixedSeq);
			shiftedSeq = Seq.reverseComplementDNA(shiftedSeq);
		}
		
		/*
System.out.println("RefJun  \t"+refSeq );
System.out.println("FixedJun\t"+fixedSeq+" "+junctionPosition);
System.out.println("ShiftedJun\t"+shiftedSeq+" "+modJunction);
*/
		
		return new String[]{refSeq, fixedSeq, shiftedSeq};
	}

	private String scoreLossOfKnownJunction(int junctionPosition, VCFRecord vcf, boolean score5Junction, boolean isPlusStrand){
		//fetch sequences and check for non GATC
		String[] rfs = fetchModifiedSpliceSequences(vcf, junctionPosition, isPlusStrand, score5Junction);
		//watch out for non GATC bases
		Matcher mat;
		for (String t: rfs){
			mat = MaxEntScanScore5.NonGATC.matcher(t);
			if (mat.find() ) return null;
		}
		boolean lostKnown = false;
		double refScore = 0;
		double altScoreFixed = 0;
		double altScoreShift = 0;
		double altScoreMax = 0;
		double d;
		String name;
		if (score5Junction){
			refScore = score5.scoreSequenceNoChecks(rfs[0]);
			//does ref meet threshold?
			if (refScore < min5Threshold) return null;
			name = "5'";
			//calc alts
			altScoreFixed = score5.scoreSequenceNoChecks(rfs[1]);
			altScoreMax = altScoreFixed;
			if (rfs[1].equals(rfs[2]) == false) {
				altScoreShift = score5.scoreSequenceNoChecks(rfs[2]);
				if (altScoreShift > altScoreMax) altScoreMax = altScoreShift;
			}
			else altScoreShift = altScoreFixed;
			
			//calc delta
			d= refScore- altScoreMax;
			if (d >= min5DeltaThreshold) lostKnown = true;
			
		}
		else {
			refScore = score3.scoreSequenceNoChecks(rfs[0]);
			//does ref meet threshold?
			if (refScore < min3Threshold) return null;
			name = "3'";
			//calc alts
			altScoreFixed = score3.scoreSequenceNoChecks(rfs[1]);
			altScoreMax = altScoreFixed;
			if (rfs[1].equals(rfs[2]) == false) {
				altScoreShift = score3.scoreSequenceNoChecks(rfs[2]);
				if (altScoreShift > altScoreMax) altScoreMax = altScoreShift;
			}
			else altScoreShift = altScoreFixed;
			//calc delta
			d= refScore- altScoreMax;
			if (d >= min3DeltaThreshold) lostKnown = true;
			
		}

		if (lostKnown){
			//insert spaces into seqs to represent junctions
			insertJunctionSpace(rfs, score5Junction);
			StringBuilder sb = new StringBuilder();
			sb.append (name+" splice junction ("+junctionPosition+") possibly damaged\n");
			sb.append ("Ref\t"+rfs[0]+"\t"+refScore+"\n");
			if (vcf.isSNP()) sb.append ("Alt\t"+rfs[1]+"\t"+altScoreFixed+"\n");
			else {
				sb.append ("AltF\t"+rfs[1]+"\t"+altScoreFixed+"\n");
				sb.append ("AltS\t"+rfs[2]+"\t"+altScoreShift+"\n");
			}
			sb.append ("Delta\t"+d+"\n");
			return sb.toString();
		}

		return null;
	}
	
	/**inserts a space in the seqs to represent the junction. Assumes correct strand*/
	private void insertJunctionSpace(String[] seqs, boolean score5Junction){
		if (score5Junction){
			//xxx xxxxxx
			for (int i=0; i< seqs.length; i++){
				seqs[i] = seqs[i].substring(0, 3)+" "+seqs[i].substring(3);
			}
		}
		else {
			//xxxxxxxxxxxxxxxxxxxx xxx
			for (int i=0; i< seqs.length; i++){
				seqs[i] = seqs[i].substring(0, 20)+" "+seqs[i].substring(20);
			}
		}
	}
	
	private String scoreNewJunction(String[] refAltSeqsToScan, boolean score5Junction, boolean exonic) {
		double[] r;
		double[] a;
		double minThres;
		double minDelta;
		String name;
		Histogram histogram;
		if (score5Junction){
			r = score5.scanSequence(refAltSeqsToScan[0]);
			a = score5.scanSequence(refAltSeqsToScan[1]);
			minThres = min5Threshold;
			minDelta = min5DeltaThreshold;
			name = "5'";
			if (exonic) histogram = exonicHistogram5;
			else histogram = intronicHistogram5;
		}
		else {
			r = score3.scanSequence(refAltSeqsToScan[0]);
			a = score3.scanSequence(refAltSeqsToScan[1]);
			minThres = min3Threshold;
			minDelta = min3DeltaThreshold;
			name = "3'";
			if (exonic) histogram = exonicHistogram3;
			else histogram = intronicHistogram3;
		}
		
		//compare scores
		StringBuilder sb = new StringBuilder();
		//snp, direct 1:1 comparison
		if (r.length == a.length) {
			for (int i=0; i< r.length; i++){
				//does alt exceed min threshold
				if (a[i] < minThres) continue;
				double testR = r[i];
				if (testR < 0) testR = 0;
				double delta = a[i] - testR;
				if (delta >= minDelta) {
					//calc pval
					double pval = Num.minus10log10(histogram.pValue(a[i]));
					sb.append(name+" "+Num.formatNumber(r[i],1)+" -> "+Num.formatNumber(a[i],1)+" index: "+i+" delta: "+Num.formatNumber(delta, 1)+" pval: "+pval+"\n");
//System.out.println(name+" "+Num.formatNumber(r[i],1)+" -> "+Num.formatNumber(a[i],1)+" index: "+i+" delta: "+Num.formatNumber(delta, 1)+" pval: "+pval);					
				}
			}
		}
		//indel, round scores and trim ends of identical scores then look at remainder for novel
		else{
			//trim 5' and 3' ends of identical scores
			int[][] ra = Num.trimIdenticalEnds(r, a);
			//find max remaining ref score or 0
			int maxRef = 0;
			if (ra[0].length !=0) maxRef = Num.maxValue(ra[0]);
			if (maxRef < 0) maxRef = 0;
			//find max remaining alt score or 0
			int maxAlt = 0;
			if (ra[1].length !=0) maxAlt = Num.maxValue(ra[1]);
			if (maxAlt < 0) maxAlt = 0;
			//if alt meets minThres and delta meets minDelta
			if (maxAlt >= minThres){
				int delta = maxAlt - maxRef;
				if (delta >= minDelta) {
					//calc pval
					double pval = Num.minus10log10(histogram.pValue(maxAlt));
					sb.append(name+" "+maxRef+" -> "+maxAlt+" delta: "+delta +" pval: "+pval+"\n");
//System.out.println(name+" "+maxRef+" -> "+maxAlt+" delta: "+delta +" pval: "+pval);
				}
			}
		}
		//add on score arrays if hits were found
		if (sb.length() !=0){
			sb.append("Ref "+refAltSeqsToScan[0]+"\t"+Num.doubleArrayToString(r, 1, "\t")+"\n");
			sb.append("Alt "+refAltSeqsToScan[1]+"\t"+Num.doubleArrayToString(a, 1, "\t")+"\n");
			return sb.toString();
		}
		return null;
	}
	
	
	/*This is going to be really slow....hmm.*/
	private UCSCGeneLine[] fetchIntersectingTranscripts(VCFRecord vcf){
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		int start = vcf.getPosition();
		int stop = start+ vcf.getAlternate()[0].length();
		for (int i=0; i< workingTranscripts.length; i++){
			if (workingTranscripts[i].intersects(start, stop)) {
				al.add(workingTranscripts[i]);
				intersectingTranscriptNames.add(workingTranscripts[i].getDisplayNameThenName()); 
			}
		}
		if (al.size() == 0) return null;
		UCSCGeneLine[] l = new UCSCGeneLine[al.size()];
		al.toArray(l);
		return l;
	}

	private void loadChromosomeData(String chrom) throws IOException {
		//set new name
		workingChromosomeName = chrom;
		//find and load sequence
		File seqFile = chromFile.get(workingChromosomeName);
		if (seqFile == null || seqFile.canRead() == false) throw new IOException ("\n\nFailed to find or load a fasta sequence for '"+workingChromosomeName+"', aborting.\n");
		workingSequence = new MultiFastaParser(seqFile).getSeqs()[0];
		workingSequence = workingSequence.toUpperCase();
		//fetch transcripts
		workingTranscripts = chromGenes.get(workingChromosomeName);
		if (workingTranscripts == null) System.out.println("\n\nWARNING: no transcripts found for this chromosome, skipping all associated vcf records.");
	}
	
	private void loadNullScoreHistograms(){
		spliceHistogram5 = new Histogram(-12, 12, 1000);
		spliceHistogram3 = new Histogram(-16, 16, 1000);
		
		if (scoreNovelExonJunctions){
			exonicHistogram5 = new Histogram(min5Threshold, 12, 1000);
			exonicHistogram3 = new Histogram(min3Threshold, 16, 1000);
		}
		if (scoreNovelIntronJunctions || scoreNovelSpliceJunctionsInSplice){
			intronicHistogram5 = new Histogram(min5Threshold, 12, 1000);
			intronicHistogram3 = new Histogram(min3Threshold, 16, 1000);
		}
		
		//for each chrom of transcripts
		System.out.println("\tChr\tSeqLength\tTranscripts");
		for (String chr: chromGenes.keySet()){
			workingChromosomeName = chr;
			workingTranscripts = chromGenes.get(workingChromosomeName);
			//fetch seq
			File seqFile = chromFile.get(workingChromosomeName);
			if (seqFile == null || seqFile.canRead() == false) {
				System.out.println ("\t"+workingChromosomeName+"\tNo fasta, skipping!");
				continue;
			}
			workingSequence = new MultiFastaParser(seqFile).getSeqs()[0];
			workingSequence = workingSequence.toUpperCase();
			System.out.println("\t"+workingChromosomeName+"\t"+workingSequence.length()+"\t"+workingTranscripts.length);
			//mask known splices
			maskKnownSplices();
			//score plus and minus strands exons
			if (scoreNovelExonJunctions){
				scoreAllExonsOrIntronsForSplices(true, true);
				scoreAllExonsOrIntronsForSplices(false, true);
			}
			if (scoreNovelIntronJunctions || scoreNovelSpliceJunctionsInSplice){
				scoreAllExonsOrIntronsForSplices(true, false);
				scoreAllExonsOrIntronsForSplices(false, false);
			}
		}
		//reset so these load
		workingSequence = "";
		workingChromosomeName = "";
		
		//print histograms
		if (printHistograms){
			System.out.println("\nHistograms of novel splices");
			if (scoreNovelExonJunctions){
				System.out.println("5' Exon Splice:");
				exonicHistogram5.printScaledHistogram();
				System.out.println("\n3' Exon Splice:");
				exonicHistogram3.printScaledHistogram();
			}
			if (scoreNovelIntronJunctions || scoreNovelSpliceJunctionsInSplice){
				System.out.println("5' Intron Splice:");
				intronicHistogram5.printScaledHistogram();
				System.out.println("\n3' Intron Splice:");
				intronicHistogram3.printScaledHistogram();
			}
			System.out.println("\nHistograms of known splices");
			System.out.println("5' Splice Junctions:");
			spliceHistogram5.printScaledHistogram();
			System.out.println("3' Splice Junctions:");
			spliceHistogram3.printScaledHistogram();
		}
		System.out.println();
	}

	private void maskKnownSplices() {
		//make array to hold info about presence of a splice base pair, default is false
		boolean[] spliceBase = new boolean[workingSequence.length()];
		int lenMinOne = spliceBase.length-1;
		HashSet<String> scored5Plus = new HashSet<String>();
		HashSet<String> scored3Plus = new HashSet<String>();
		HashSet<String> scored5Minus = new HashSet<String>();
		HashSet<String> scored3Minus = new HashSet<String>();
		//for each transcript
		for (UCSCGeneLine g: workingTranscripts){
			ExonIntron[] introns = g.getIntrons();
			if (introns == null) continue;
			boolean isPlusStrand = g.getStrand().equals("+");
			//mask both 5' and 3' junctions
			int startJunc = 0;
			int endJunc = 0;
			for (int i=0; i< introns.length; i++){
				//plus strand 
				if (isPlusStrand){
					//5'
					startJunc = introns[i].getStart()-3;     
					endJunc = introns[i].getStart()+6;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Plus.contains(coor) == false){
						scored5Plus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						double score = score5.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) spliceHistogram5.count(score);
						//System.out.println("+ 5\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);					
					}
					//3'
					startJunc = introns[i].getEnd()-20;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
					coor = startJunc+"-"+endJunc;
					if (scored3Plus.contains(coor) == false){
						scored3Plus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						double score = score3.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) spliceHistogram3.count(score);
						//System.out.println("+ 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);	
					}
				}
				//minus strand
				else {
					//3'
					startJunc = introns[i].getStart()-3;
					endJunc = introns[i].getStart()+20;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
					
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored3Minus.contains(coor) == false){
						scored3Minus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						seq = Seq.reverseComplementDNA(seq);
						double score = score3.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) spliceHistogram5.count(score);
						//System.out.println("- 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);	
					}
					
					//5'
					startJunc = introns[i].getEnd()-6;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
					
					coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Minus.contains(coor) == false){
						scored5Minus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						seq = Seq.reverseComplementDNA(seq);
						double score = score5.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) spliceHistogram5.count(score);
						//System.out.println("- 5\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);	
					}
				}
			}
		}
		char[] bases = workingSequence.toCharArray();
		for (int i=0; i< bases.length; i++){
			if (spliceBase[i]) bases[i] = 'N';
		}
		workingSequence = new String(bases);
		bases = null;
	}

	private String[] fetchNMerSeqs(VCFRecord vcf, int nMer, boolean reverseComplement) {
		String refSeq = null;
		String altSeq = null;
		//set coordinates
		int pos = vcf.getPosition();
		int start = pos - nMer +1;
		int stop = 0;
//System.out.println ("\nVar\t"+vcf.getReference()+" -> "+vcf.getAlternate()[0]);
			stop = vcf.getPosition()+ nMer + vcf.getReference().length() -1;
			refSeq = workingSequence.substring(start, stop);
//System.out.println ("RefSeq\t"+refSeq);
//System.out.println ("AltSeq\t"+workingSequence.substring(start, pos) +" "+ vcf.getAlternate()[0] +" "+ workingSequence.substring(pos+vcf.getReference().length(), stop));
			altSeq = workingSequence.substring(start, pos) + vcf.getAlternate()[0] + workingSequence.substring(pos+vcf.getReference().length(), stop);

		//reverse comp it?
		if (reverseComplement){
			refSeq = Seq.reverseComplementDNA(refSeq);
			altSeq = Seq.reverseComplementDNA(altSeq);
		}
		return new String[]{refSeq, altSeq};
	}

	public void loadAndModifyHeader(BufferedReader in, PrintWriter out) throws IOException{
		boolean addedInfo = false;
		boolean foundChrom = false;
		String line;
		while ((line=in.readLine()) != null){
			//comments
			if (line.startsWith("#")){
				//add info lines?
				if (addedInfo == false && line.startsWith("##INFO=")){
					addedInfo = true;
					//out.println("##INFO=<ID=MES3,Number=1,Type=Float,Description=\"MaxEntScan 3' splice score.\">");
					//out.println("##INFO=<ID=MES5,Number=1,Type=Float,Description=\"MaxEntScan 5' splice score.\">");
				}
				//out.println(line);
				if (line.startsWith("#CHROM")){
					foundChrom = true;
					break;
				}
			}
		}
		if (foundChrom == false) throw new IOException("\tError: Failed to find the #CHROM header line? Aborting.\n");
		if (addedInfo == false) throw new IOException("\tError: Failed to find any ##INFO header lines? Aborting.\n");
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFSpliceAnnotator(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'f': chromDirectory = new File(args[++i]); break;
					case 'm': spliceModelDirectory = new File(args[++i]); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'e': scoreNovelExonJunctions = false; break;
					case 'i': scoreNovelIntronJunctions = false; break;
					case 's': scoreNovelSpliceJunctionsInSplice = false; break;
					case 'a': min5Threshold = Double.parseDouble(args[++i]); break;
					case 'b': min3Threshold = Double.parseDouble(args[++i]); break;
					case 'c': min5DeltaThreshold = Double.parseDouble(args[++i]); break;
					case 'd': min3DeltaThreshold = Double.parseDouble(args[++i]); break;
					
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//checkfiles
		if (chromDirectory == null || chromDirectory.isDirectory() == false ) {
			Misc.printErrAndExit("\nError: please provide a path to a directory containing chromosome specific fasta files (e.g. chr1.fasta.gz, chr2.fasta.zip, chr3.fa, ...)\n");
		}
		if (spliceModelDirectory == null || spliceModelDirectory.isDirectory() == false ) {
			Misc.printErrAndExit("\nError: please provide a path to a directory containing the splice models.\n");
		}
		if (refSeqFile == null || refSeqFile.canRead()== false){
			Misc.printErrAndExit("\nPlease enter a transcript table file in ucsc refflat format.\n");
		}
		
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a vcf file or directory containing such to annotate splice junctions.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            VCF Splice Annotator : Dec 2013                       **\n" +
				"**************************************************************************************\n" +
				"WARNING: beta!\n"+
				"\n"+
				"Scores variants for changes in splicing using the MaxEntScan algorithms. See Yeo and\n"+
				"Burge 2004, http://www.ncbi.nlm.nih.gov/pubmed/15285897 for details. Known splice\n"+
				"acceptors and donors are scored for loss of a junction.  Exonic, intronic, and splice\n"+
				"bases are scanned for novel junctions.  Note, indels in exons causing\n"+
				"frameshifts are not annotated. This app only looks for changes in splicing.\n\n" +

				"Required Options:\n"+
				"-v VCF file or directory containing such (xxx.vcf(.gz/.zip OK)).\n"+
				"-f Fasta file directory, chromosome specific xxx.fa/.fasta(.zip/.gz OK) files.\n" +
				"-u UCSC RefFlat or RefSeq transcript (not merged genes) file, full path. See RefSeq \n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 .\n"+
				"-m Full path directory name containing the me2x3acc1-9, splice5sequences and me2x5\n"+
				"       splice model files. See USeq/Documentation/ or \n"+
				"       http://genes.mit.edu/burgelab/maxent/download/ \n"+
				
				"\n"+
				"Optional options:\n"+
				"-e Don't scan exonic bases for novel splice junctions.\n"+
				"-i Don't scan intronic bases for novel splice junctions.\n"+
				"-s Don't scan known splice junctions for novel splice junctions.\n"+
				"-a Minimum 5' threshold for scoring the presence of a splice junction, defaults to 5.\n"+
				"-b Minimum 3' threshold for scoring the presence of a splice junction, defaults to 5.\n"+
				"-c Minimum difference for loss or gain of a 5' splice junction, defaults to 5.\n"+
				"-d Minimum difference for loss or gain of a 3' splice junction, defaults to 5.\n"+
				"\n"+

				"Example: java -Xmx10G -jar ~/USeq/Apps/VCFSpliceAnnotator -f ~/Hg19/Fa/ -v ~/exm2.vcf\n"+
				"       -m ~/USeq/Documentation/splicemodels -i -u ~/Hg19/hg19EnsTran.ucsc.zip \n\n"+

				"**************************************************************************************\n");

	}
	
}
