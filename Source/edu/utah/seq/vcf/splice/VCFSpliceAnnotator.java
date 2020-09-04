package edu.utah.seq.vcf.splice;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.mes.*;
import edu.utah.seq.vcf.VCFLookUp;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;
import util.bio.annotation.ExonIntron;
import util.bio.annotation.ExportIntergenicRegions;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.bio.seq.Seq;
import util.gen.*;

/**Depreciated, don't use!  Use VCFSpliceScanner
 * @author Nix
 * */
public class VCFSpliceAnnotator {
	
	/*@TODO:
	 * don't reprocess sjs that have already been scored
	 * memory issues?
	 * correct the pvalues
	 * thread*/

	//user fields
	private File[] vcfFiles;
	private File spliceModelDirectory;
	private IndexedFastaSequenceFile fasta; 
	private File histogramFile;
	private File transcriptSeqFile;
	private File ccdsTranscriptFile;
	private File samSpliceFile;
	private File saveDirectory;
	private int minimumSpliceJunctionCoverage = 10;

	private double min5Threshold = 2;
	private double min5DeltaThreshold = 2;
	private double min3Threshold = 2;
	private double min3DeltaThreshold = 2;
	private double minPValue = 10;
	
	private boolean scoreNovelIntronJunctions = true;
	private boolean scoreNovelExonJunctions = true;
	private boolean scoreNovelSpliceJunctionsInSplice = true;
	private short vcfExportCategory = 0;
	
	//internal fields
	private double minimumFractionCorrectlySpliced = 0.9;
	private MaxEntScanScore5 score5;
	private MaxEntScanScore3 score3;
	private LinkedHashMap<String, ArrayList<SpliceHit>> geneNameSpliceHits = null;
	private ArrayList<SpliceJunction> allSpliceJunctions = new ArrayList<SpliceJunction>();
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private HashMap<String,UCSCGeneLine[]> ccdsChromGenes;
	private String workingChromosomeName = "";
	private String workingSequence = null;
	private boolean workingTranscriptIsPlusStrand;
	private IntervalST<ArrayList<UCSCGeneLine>> workingGeneTree;
	private UCSCGeneLine[] workingTranscripts;
	private Histogram exonicHistogram5;
	private Histogram exonicHistogram3;
	private Histogram intronicHistogram5;
	private Histogram intronicHistogram3;
	private Histogram spliceHistogram5;
	private Histogram spliceHistogram3;
	private HashMap<String, Histogram> typeHist = new HashMap<String, Histogram>();
	private boolean printHistograms = true;
	private SamReader  samSpliceReader;
	private Gzipper vcfOut;
	
	private int numTranscripts = 0; 
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
	
	//for multiple testing correction, to do! not yet implemented for vcf file, just spreadsheet
	private double num5GainedTested = 0;
	private double num3GainedTested = 0;
	private double num5LostTested = 0;
	private double num3LostTested = 0; 
	
	//constructor
	public VCFSpliceAnnotator(String[] args){
		try {
			//start clock
			long startTime = System.currentTimeMillis();

			//process args
			processArgs(args);
			printThresholds();

			//start up mes
			score5 = new MaxEntScanScore5(spliceModelDirectory);
			score3 = new MaxEntScanScore3(spliceModelDirectory);

			//load transcripts
			System.out.println("Loading transcripts...");
			loadTranscripts();

			//histograms of novel exonic splices
			System.out.println("Loading splice score histograms...");
			loadNullScoreHistograms();

			//for each vcf file
			System.out.println("Processing...");
			for (int i=0; i< vcfFiles.length; i++){
				System.out.println("\t"+vcfFiles[i]);
				geneNameSpliceHits = new LinkedHashMap<String, ArrayList<SpliceHit>>();
				annotateVCFWithSplices(vcfFiles[i]);
				correctPValues();
				printSpreadSheetResults(vcfFiles[i]);
				System.out.println();
				printSummary();
				zeroSummaryStats();
			}

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
			System.out.println("\nDone "+Math.round(diffTime)+" min!\n");
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	private void correctPValues() {
		double logG5 = Num.log10(num5GainedTested);
		double logG3 = Num.log10(num3GainedTested);
		double logL5 = Num.log10(num5LostTested);
		double logL3 = Num.log10(num3LostTested);
		System.out.println("\nApplying a Bonferroni correction to the spreadsheet pvalues...");
		System.out.println("\tCorr log vals:\tG3 "+logG3+"\tG5 "+logG5+"\tD3 "+logL3+"\tD5 "+logL5);
		for (SpliceJunction sj: allSpliceJunctions) {
			//G3 G5 D5 D3
			double sub;
			String type = sj.getType();
			if (type.startsWith("G3")) sub = logG3;
			else if (type.startsWith("G5")) sub = logG5;
			else if (type.startsWith("D3")) sub = logL3;
			else sub = logL5;
			//watchout for cases when the log val is infinite due to taking log of zero. 
			if (Double.isInfinite(sub) || Double.isNaN(sub)) sub = 0;
			double pval = sj.getTransPValue() - sub;
			if (pval < 0) pval = 0;
			sj.setTransPValue(pval);
		}
		
	}


	//methods
	private void printSpreadSheetResults(File vcfFile) {
		File modFile = new File (saveDirectory, Misc.removeExtension(vcfFile.getName())+".xls");
		PrintWriter out;
		try {
			out = new PrintWriter( new FileWriter (modFile));
			//save header
			out.println("GeneName\tTranscriptNames\tChrom\tVCF RefSeq\tVCF AltSeq\tVCF Pos\tSJ Type\tSJ -10Log10(adjPVal)\tSJ Pos\tSJ RefSeq\tSJ AltSeq\tSJ RefScore\tSJ AltScore\tSJ RefZScore\tSJ AltZScore\tVCF Record Fields...");
			LinkedHashMap<String, ArrayList<String>> dataTranscripts = new LinkedHashMap <String, ArrayList<String>>();
			//walk through hash collapsing transcripts with same hit
			for (String geneName : geneNameSpliceHits.keySet()) {
				dataTranscripts.clear();
				//for each hit
				ArrayList<SpliceHit> spliceHits = geneNameSpliceHits.get(geneName);
				//for each splice hit				
				for (SpliceHit sh: spliceHits){
					String transcriptName = sh.getTranscript().getName();
					StringBuilder sb = new StringBuilder();
					//chrom
					sb.append(sh.getTranscript().getChrom()); sb.append("\t");
					//vcf refseq
					sb.append(sh.getVcf().getReference()); sb.append("\t");
					//vcf altseq
					sb.append(sh.getVcf().getAlternate()[sh.getVcfAltIndex()]); sb.append("\t");
					//vcf pos
					sb.append(sh.getVcf().getPosition()); sb.append("\t");
					String head = sb.toString();
					
					//for each splice junction affect, typically only one
					ArrayList<SpliceJunction> sjAL = sh.getAffectedSpliceJunctions();					
					for (SpliceJunction sj: sjAL){
						//pass thresholds after multiple testing correction?
						if (sj.getTransPValue() < minPValue) continue;
						sb = new StringBuilder();
						sb.append(head); 
						//sj type
						sb.append(sj.getType()); sb.append("\t");
						//sj pval
						sb.append(sj.getTransPValue()); sb.append("\t");
						//sj pos
						sb.append(sj.getPosition()); sb.append("\t");
						//sj ref seq
						sb.append(sj.getReferenceSequence()); sb.append("\t");
						//sj alt seq
						sb.append(sj.getAlternateSequence()); sb.append("\t");
						//sj ref score
						sb.append(sj.getReferenceScore()); sb.append("\t");
						//sj alt score
						sb.append(sj.getAlternateScore()); sb.append("\t");
						
						//sj ref zscore
						StandardDeviation sd = typeHist.get(sj.getType()).getStandardDeviation();
						sb.append(sd.getZScore(sj.getReferenceScore())); sb.append("\t");
						//sj alt zscore
						sb.append(sd.getZScore(sj.getAlternateScore())); sb.append("\t");
						
						//vcf record
						sb.append(sh.getVcf().toString());
						String data = sb.toString();
						//add to hash
						ArrayList<String> trans = dataTranscripts.get(data);
						if (trans == null){
							trans = new ArrayList<String>();
							dataTranscripts.put(data, trans);
						}
						trans.add(transcriptName);
					}
				}
				//print out each data line with concat of transcriptNames
				for (String data: dataTranscripts.keySet()){
					ArrayList<String> trans = dataTranscripts.get(data);
					String transConcat = Misc.stringArrayListToString(trans, ",");
					out.println(geneName+"\t"+transConcat+"\t"+data);
				}
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void zeroSummaryStats(){
		intersectingTranscriptNames.clear();
		numVariantsScanned = 0;
		numVariantsIntersectingTranscripts = 0;
		numVariantsIntersectingExons = 0;
		numVariantsIntersectingIntrons = 0;
		numVariantsIntersectingSJs = 0;
		numExonSJsGained = 0;
		numIntronSJsGained = 0;
		numSpliceJunctionSJsGained = 0;
		numSJsLost = 0;
		num5GainedTested = 0;
		num3GainedTested = 0;
		num5LostTested = 0;
		num3LostTested = 0; 
		allSpliceJunctions = new ArrayList<SpliceJunction>();
	}

	private void printThresholds(){
		StringBuilder sb = new StringBuilder();
		sb.append("Threholds:\n");
		sb.append(min5Threshold +"\tMinimum 5' splice junction threshold for scoring the presence of a junction\n");
		sb.append(min3Threshold +"\tMinimum 3' splice junction threshold for scoring the presence of a junction\n");
		sb.append(min5DeltaThreshold +"\tMinimum score difference for loss or gain of a 5' splice junction\n");
		sb.append(min3DeltaThreshold +"\tMinimum score difference for loss or gain of a 3' splice junction\n");
		sb.append(minPValue +"\tMinimum -10Log10(pval) for reporting a loss or gain of a splice junction\n");
		sb.append(scoreNovelIntronJunctions +"\tLook for novel intron junctions, outside of known splice junctions.\n");
		sb.append(scoreNovelExonJunctions +"\tLook for novel exon junctions, outside of known splice junctions.\n");
		sb.append(scoreNovelSpliceJunctionsInSplice +"\tLook for novel splice junctions inside known splice junctions.\n");
		sb.append(vcfExportCategory +"\tExport catagory for adding info to vcf records.\n");
		System.out.println(sb);
	}
	
	private void printSummary(){
		int total = numVariantsIntersectingIntrons+ numVariantsIntersectingExons + numVariantsIntersectingSJs;
		StringBuilder sb = new StringBuilder();
		sb.append("Summary stats:\n");
		sb.append(numTranscripts+"\tTranscripts\n"); 
		sb.append(intersectingTranscriptNames.size()+"\tTranscripts intersecting variants\n"); 
		sb.append(numVariantsScanned+"\tVariants (including alternates)\n");
		sb.append(numVariantsIntersectingTranscripts+"\tVariants intersecting transcripts, no repeat counting.\n");
		sb.append("\nCounts including repeats with overlapping transcripts:\n");
		sb.append(total + "\tVariants intersecting annotations\n");
		sb.append(numVariantsIntersectingExons+"\tVariants intersecting exons (non splice junction)\n"); 
		sb.append(numVariantsIntersectingIntrons+"\tVariants intersecting introns (non splice junction)\n"); 
		sb.append(numVariantsIntersectingSJs+"\tVariants intersecting splice junctions\n"); 
		sb.append(numExonSJsGained+"\tExonic splice junctions gained\n"); 
		sb.append(numIntronSJsGained+"\tIntronic splice junctions gained\n"); 
		sb.append(numSpliceJunctionSJsGained+"\tKnown splice junctions with novel splice junction gained\n"); 
		sb.append(numSJsLost+"\tKnown splice junctions damaged\n");
		System.out.println(sb);
	}

	/**Adds splice junction information to each record. */
	public void annotateVCFWithSplices(File vcfFile) {
		BufferedReader in = null;
		VCFParser parser = new VCFParser();
		try {
			in  = IO.fetchBufferedReader(vcfFile);
			vcfOut = new Gzipper(new File(saveDirectory, Misc.removeExtension(vcfFile.getName())+"_VCFSA.vcf.gz"));
			//add ##INFO line and find "#CHROM" line 
			loadAndModifyHeader(in);

			//For each record
			String line;
			HashSet<String> effects = new HashSet<String>();
			while ((line=in.readLine()) != null){
				//parse record
				VCFRecord vcf = new VCFRecord(line, parser, true, true);
				
				//load new sequence and transcripts?
				if (vcf.getChromosome().equals(workingChromosomeName) == false) {
					//load new working seq
					loadChromosomeData(vcf.getChromosome());
				}
				
				//check if any transcripts were found
				if (workingTranscripts == null) {
					vcfOut.println(line);
					continue;
				}

				//for each alternate allele, some vcf records have several
				String[] alts = vcf.getAlternate();
				effects.clear();
				for (int i=0; i< alts.length; i++){	
					numVariantsScanned++;
					vcf.setAlternate(new String[]{alts[i]});
					
					//fetch intersecting transcripts
					UCSCGeneLine[] trans = fetchIntersectingTranscripts(vcf);
					if (trans == null) continue;
					numVariantsIntersectingTranscripts++;
					
					//for each transcript score effect
					for (UCSCGeneLine l: trans){
						SpliceHit sh = new SpliceHit (vcf, i, l);
						scoreVariantTranscript(sh);
						
						//any changes? add to hash for spreadsheet output
						if (sh.getAffectedSpliceJunctions() != null) {
							String geneName = sh.getTranscript().getDisplayName();
							ArrayList<SpliceHit> al = geneNameSpliceHits.get(geneName);
							if (al == null) {
								al = new ArrayList<SpliceHit>();
								geneNameSpliceHits.put(geneName, al);
							}
							al.add(sh);
							//add vcf entry
							effects.addAll(sh.getVcfEntries(vcfExportCategory));
						}
					}
				}
				vcf.setAlternate(alts);
				//modify INFO field?
				if (effects.size() == 0) vcfOut.println(vcf.getOriginalRecord());
				else {
					String e = Misc.hashSetToString(effects, ":");
					String[] fields = VCFParser.TAB.split(vcf.getOriginalRecord());
					fields[7] = fields[7]+";VCFSA="+e;
					String f = Misc.stringArrayToString(fields, "\t");
					vcfOut.println(f);
				}
			}
			//last
			System.out.println();
			in.close();
			vcfOut.close();
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

	private void scanIntronsForNovelSpliceJunctions (SpliceHit hit){
		VCFRecord vcf = hit.getVcf();
		ExonIntron[] introns = hit.getTranscript().getIntrons();		
		if (introns != null) {
			int leftAdd = 6;
			int rightSub = 20;
			if (workingTranscriptIsPlusStrand == false){
				leftAdd = 20;
				rightSub = 6;
			}
			int varStart = vcf.getPosition();
			int varEnd = varStart+ vcf.getReference().length();
			boolean intersection = false;
			for (int i=0; i< introns.length; i++){
				//mod coordinates of intron
				int intronStart = introns[i].getStart()+leftAdd;
				int intronEnd = introns[i].getEnd()-rightSub;
				if (intronEnd <= intronStart) continue;
				//if a hit then score 
				if (intersects(varStart, varEnd, intronStart, intronEnd) == false) continue;
				intersection = true;
				//scan for novel, reference and alternate
				String[] fiveSeqs = fetchNMerSeqs(vcf, 9);
				scoreNewJunction(fiveSeqs, true, false, false, hit);
				String[] threeSeqs = fetchNMerSeqs(vcf, 23);
				scoreNewJunction(threeSeqs, false, false, false, hit);
				//can only have one intersection if snp or insertion
				if (vcf.isSNP() || vcf.isInsertion()) break;
			}
			if (intersection) numVariantsIntersectingIntrons++;
			
		}
	}

	private void scanExonsForNovelSpliceJunctions (SpliceHit hit){
		VCFRecord vcf = hit.getVcf();
		ExonIntron[] exons = hit.getTranscript().getExons();
		if (exons.length != 0) {
			int varStart = vcf.getPosition();
			int varEnd = varStart+ vcf.getReference().length();
			int lastExonIndex = exons.length-1;
			boolean intersection = false;
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
				intersection = true;
				//scan for novel, reference and alternate			
				String[] fiveSeqs = fetchNMerSeqs(vcf, 9);
				scoreNewJunction(fiveSeqs, true, true, false, hit);
				String[] threeSeqs = fetchNMerSeqs(vcf, 23);
				scoreNewJunction(threeSeqs, false, true, false, hit);
				//can only have one intersection if snp or insertion
				if (vcf.isSNP() || vcf.isInsertion()) break;
			}
			if (intersection) numVariantsIntersectingExons++;
		}
	}
	
	private boolean intersects(int varStart, int varEnd, int annoStart, int annoEnd){
		if (varStart >= annoEnd) return false;
		if (varEnd <= annoStart) return false;
		return true;
	}
	
	private void scanKnownSplices (SpliceHit hit){
		VCFRecord vcf = hit.getVcf();
		ExonIntron[] introns = hit.getTranscript().getIntrons();
		
		int varStart = vcf.getPosition();
		int varEnd = varStart+ vcf.getReference().length();
		boolean intersectsSplice = false;
		if (introns!= null) {
			//find junction
			int junc = 0;
			int startJunc = 0;
			int endJunc = 0;
			for (int i=0; i< introns.length; i++){
				SpliceJunction sj3 = null;
				SpliceJunction sj5 = null;
				//plus strand 
				if (workingTranscriptIsPlusStrand){
					//5'?
					startJunc = introns[i].getStart()-3;     
					endJunc = introns[i].getStart()+6;
					if (intersects(varStart, varEnd, startJunc, endJunc)) {
						junc = introns[i].getStart();
						sj5 = scoreLossOfKnownJunction(junc, vcf, true);
						intersectsSplice = true;
					}
					//3'?
					startJunc = introns[i].getEnd()-20;
					endJunc = introns[i].getEnd()+3;
					if (intersects(varStart, varEnd, startJunc, endJunc)) {
						junc = introns[i].getEnd();
						sj3 = scoreLossOfKnownJunction(junc, vcf, false);
						intersectsSplice = true;
					}
				}
				//minus strand
				else {
					//3'?
					startJunc = introns[i].getStart()-3;
					endJunc = introns[i].getStart()+20;
					if (intersects(varStart, varEnd, startJunc, endJunc)) {
						junc = introns[i].getStart();
						sj3 = scoreLossOfKnownJunction(junc, vcf, false);
						intersectsSplice = true;
					}
					//5'?
					startJunc = introns[i].getEnd()-6;
					endJunc = introns[i].getEnd()+3;
					if (intersects(varStart, varEnd, startJunc, endJunc)) {
						junc = introns[i].getEnd();
						sj5 = scoreLossOfKnownJunction(junc, vcf, true);
						intersectsSplice = true;
					}
				}
				//damage?
				if (sj3 != null) hit.saveSpliceJunction(sj3);
				if (sj5 != null) hit.saveSpliceJunction(sj5);
			}
			if (intersectsSplice){
				numVariantsIntersectingSJs++;
				//now scan for novel new splice in splice site
				if (scoreNovelSpliceJunctionsInSplice){
					String[] fiveSeqs = fetchNMerSeqs(vcf, 9);
					scoreNewJunction(fiveSeqs, true, false, true, hit);
					String[] threeSeqs = fetchNMerSeqs(vcf, 23);
					scoreNewJunction(threeSeqs, false, false, true, hit);
				}
			}
			
		}
	}

	/**Returns null if no actionable info.  Otherwise with message.*/
	private void scoreVariantTranscript(SpliceHit hit) {
		//set working vals
		workingTranscriptIsPlusStrand = hit.getTranscript().getStrand().equals("+");
		//check to see if variant is a snp or insertion, if so then a hit is mutually exclusive to that annotation class so can skip. For deletions, must scan all.
		boolean isSNPInsertion = true;
		if (hit.getVcf().isDeletion()) isSNPInsertion = false;
		//intronic and away from splice junctions
		if (scoreNovelIntronJunctions) scanIntronsForNovelSpliceJunctions(hit);
		if (isSNPInsertion && hit.getAffectedSpliceJunctions() != null) return;
		//exonic and away from splice junction
		if (scoreNovelExonJunctions) scanExonsForNovelSpliceJunctions(hit);
		if (isSNPInsertion && hit.getAffectedSpliceJunctions() != null) return;
		//scan hits to splices, gain and loss
		scanKnownSplices(hit);
	}
	
	/**Returns three sequences, the reference splice junction, the modified splice junction where the position has been fixed or follows the original.
	 * {refSeq, fixedSeq, shiftedSeq} */
	private String[] fetchModifiedSpliceSequences(VCFRecord vcf, int junctionPosition, boolean fetch5Prime){
		//set cut points
		int leftSub;
		int rightAdd;
		if (workingTranscriptIsPlusStrand){
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
		
		if (workingTranscriptIsPlusStrand == false){
			refSeq = Seq.reverseComplementDNA(refSeq);
			fixedSeq = Seq.reverseComplementDNA(fixedSeq);
			shiftedSeq = Seq.reverseComplementDNA(shiftedSeq);
		}
		return new String[]{refSeq, fixedSeq, shiftedSeq};
	}

	/**Returns null if no damage passing thresholds or SpliceJunction.*/
	private SpliceJunction scoreLossOfKnownJunction(int junctionPosition, VCFRecord vcf, boolean score5Junction){
		//fetch sequences and check for non GATC
		String[] rfs = fetchModifiedSpliceSequences(vcf, junctionPosition, score5Junction);
		//watch out for non GATC bases
		Matcher mat;
		for (String t: rfs){
			mat = MaxEntScanScore5.NonGATC.matcher(t);
			if (mat.find() ) return null;
		}
		SpliceJunction sj = null;
		double refScore = 0;
		double altScoreMax = 0;
		String altScoreMaxSeq = null;
		double pval = 0;
		if (score5Junction){
			refScore = score5.scoreSequenceNoChecks(rfs[0]);
			//does ref meet threshold?
			if (refScore < min5Threshold) return null;
			//calc alts
			double altScoreFixed = score5.scoreSequenceNoChecks(rfs[1]);
			altScoreMax = altScoreFixed;
			altScoreMaxSeq = rfs[1];
			if (rfs[1].equals(rfs[2]) == false) {
				double altScoreShift = score5.scoreSequenceNoChecks(rfs[2]);
				if (altScoreShift > altScoreMax) {
					altScoreMax = altScoreShift;
					altScoreMaxSeq = rfs[2];
				}
			}
			//calc delta
			double d= refScore- altScoreMax;
			if (d >= min5DeltaThreshold) {
				sj =  new SpliceJunction('D', '5', 'S', junctionPosition);
				pval = Num.minus10log10(spliceHistogram5.pValue(altScoreMax, false));
				num5LostTested++;
			}
		}
		else {
			refScore = score3.scoreSequenceNoChecks(rfs[0]);
			//does ref meet threshold?
			if (refScore < min3Threshold) return null;
			//calc alts
			double altScoreFixed = score3.scoreSequenceNoChecks(rfs[1]);
			altScoreMax = altScoreFixed;
			altScoreMaxSeq = rfs[1];
			if (rfs[1].equals(rfs[2]) == false) {
				double altScoreShift = score3.scoreSequenceNoChecks(rfs[2]);
				if (altScoreShift > altScoreMax) {
					altScoreMax = altScoreShift;
					altScoreMaxSeq = rfs[2];
				}
			}
			//calc delta
			double d= refScore- altScoreMax;
			if (d >= min3DeltaThreshold) {
				sj =  new SpliceJunction('D', '3', 'S', junctionPosition);
				pval = Num.minus10log10(spliceHistogram3.pValue(altScoreMax, false));
				num3LostTested++;
			}
		}
		//damaged?
		if (pval >= minPValue && sj != null){
			numSJsLost++;
			allSpliceJunctions.add(sj);
			//add info to sj
			sj.setReferenceSequence(rfs[0]);
			sj.setReferenceScore(refScore);
			sj.setAlternateSequence(altScoreMaxSeq);
			sj.setAlternateScore(altScoreMax);
			sj.setTransPValue(pval);
			return sj;
		}
		//nope
		return null;
	}

	
	private void scoreNewJunction(String[] refAltSeqsToScan, boolean score5Junction, boolean exonic, boolean splice, SpliceHit spliceHit) {
		double[] r;
		double[] a;
		double minThres;
		double minDelta;
		Histogram histogram;
		SpliceJunction spliceJunction;
		if (score5Junction){
			r = score5.scanSequence(refAltSeqsToScan[0], -1000);
			a = score5.scanSequence(refAltSeqsToScan[1], -1000);
			minThres = min5Threshold;
			minDelta = min5DeltaThreshold;
			if (exonic) {
				histogram = exonicHistogram5;
				spliceJunction = new SpliceJunction('G', '5', 'E');
			}
			else {
				histogram = intronicHistogram5;
				spliceJunction = new SpliceJunction('G', '5', 'I');
			}
		}
		else {
			r = score3.scanSequence(refAltSeqsToScan[0], -1000);
			a = score3.scanSequence(refAltSeqsToScan[1], -1000);
			minThres = min3Threshold;
			minDelta = min3DeltaThreshold;
			if (exonic) {
				histogram = exonicHistogram3;
				spliceJunction = new SpliceJunction('G', '3', 'E');
			}
			else {
				histogram = intronicHistogram3;
				spliceJunction = new SpliceJunction('G', '3', 'I');
			}
		}
		
		//compare scores
		//snp, direct 1:1 comparison
		if (r.length == a.length) {
			double maxPval = 0;
			//scan all for maxPval, must first meet minimum delta
			for (int i=0; i< r.length; i++){
				//does alt exceed min threshold, if it's -1000 it'll be skipped
				if (a[i] < minThres) continue;
				//skip due to bad base in reference?
				if (r[i] == -1000) continue;
				double testR = r[i];
				if (testR < 0) testR = 0;
				double delta = a[i] - testR;
				if (delta >= minDelta) {
					//calc pval
					double pval = Num.minus10log10(histogram.pValue(a[i], true));
					if (score5Junction) num5GainedTested++;
					else num3GainedTested++;
					if (pval < maxPval) continue;
					maxPval = pval;
					//good so save info in spliceJunction
					spliceJunction.setTransPValue(pval);
					spliceJunction.setReferenceScore(r[i]);
					spliceJunction.setAlternateScore(a[i]);
					int posOfZeroSeq;
					int relPosSpliceInSeq;
					if (score5Junction) {
						spliceJunction.setReferenceSequence(refAltSeqsToScan[0].substring(i, i+9));
						spliceJunction.setAlternateSequence(refAltSeqsToScan[1].substring(i,i+9));
						posOfZeroSeq = spliceHit.getVcf().getPosition() - 9 +1;
						relPosSpliceInSeq = i+3;
					}
					else {
						spliceJunction.setReferenceSequence(refAltSeqsToScan[0].substring(i, i+23));
						spliceJunction.setAlternateSequence(refAltSeqsToScan[1].substring(i,i+23));
						posOfZeroSeq = spliceHit.getVcf().getPosition() -23 +1;
						relPosSpliceInSeq = i+20;
					}
					//if neg strand then flip
					if (workingTranscriptIsPlusStrand == false){
						relPosSpliceInSeq = refAltSeqsToScan[0].length() - relPosSpliceInSeq;
					}
					spliceJunction.setPosition(posOfZeroSeq + relPosSpliceInSeq);
					
					//System.out.println("\nHereeeeSNP "+workingTranscriptIsPlusStrand+"\n"+spliceJunction);
					//System.out.println(Num.formatNumber(r[i],1)+" -> "+Num.formatNumber(a[i],1)+" index: "+i+" delta: "+Num.formatNumber(delta, 1)+" pval: "+pval);					
					
 
				}
			}
		}
		
		//indel, round scores and trim ends of identical scores then look at remainder for novel
		else{
			//trim 5' and 3' ends of identical scores after rounding to ints
			double[][] ra = Num.trimIdenticalEnds(r, a);
			//find max remaining ref score or 0
			double maxRef = 0;
			double actualMaxRef = 0;
			if (ra[0].length !=0) {
				maxRef = Num.maxValue(ra[0]);
				actualMaxRef = maxRef;
			}
			if (maxRef < 0) maxRef = 0;
			//find max remaining alt score or 0
			double maxAlt = 0;
			if (ra[1].length !=0) maxAlt = Num.maxValue(ra[1]);
			if (maxAlt < 0) maxAlt = 0;
			//if alt meets minThres and delta meets minDelta
			if (maxAlt >= minThres){
				double delta = maxAlt - maxRef;
				if (delta >= minDelta) {
					//calc pval
					double pval = Num.minus10log10(histogram.pValue(maxAlt, true));
					if (score5Junction) num5GainedTested++;
					else num3GainedTested++;
					spliceJunction.setTransPValue(pval);
					spliceJunction.setReferenceScore(actualMaxRef);
					spliceJunction.setReferenceSequence(refAltSeqsToScan[0]);
					spliceJunction.setAlternateScore(maxAlt);
					spliceJunction.setAlternateSequence(refAltSeqsToScan[1]);
					spliceJunction.setPosition(spliceHit.getVcf().getPosition());
//System.out.println("\nHereeeeINDEL"+workingTranscriptIsPlusStrand+"\n"+spliceJunction);
//System.out.println(maxRef+" -> "+maxAlt+" delta: "+delta +" pval: "+pval);
				}
			}
		}
		
		//pass pval threshold? save splice junction in splice hit
		if (spliceJunction.getTransPValue() >= minPValue){
			spliceHit.saveSpliceJunction(spliceJunction);
			allSpliceJunctions.add(spliceJunction);
			//increment counters
			if (exonic) numExonSJsGained++;
			else if (splice) numSpliceJunctionSJsGained++;
			else numIntronSJsGained++;
//System.out.println("Ref "+refAltSeqsToScan[0]+"\t"+Num.doubleArrayToString(r, 1, "\t")+"\n");
//System.out.println("Alt "+refAltSeqsToScan[1]+"\t"+Num.doubleArrayToString(a, 1, "\t")+"\n");
		}
	}

	
	private UCSCGeneLine[] fetchIntersectingTranscripts(VCFRecord vcf){
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		int start = vcf.getPosition();
		int stop = start+ vcf.getAlternate()[0].length();
		
		//end is included in search so subtract 1
		Iterable<Interval1D> it = workingGeneTree.searchAll(new Interval1D(start, stop-1));
		for (Interval1D x : it) {
			ArrayList<UCSCGeneLine> genes = workingGeneTree.get(x);
			al.addAll(genes);
			for (UCSCGeneLine l: genes)	intersectingTranscriptNames.add(l.getDisplayNameThenName()); 
		}
		if (al.size() == 0) return null;
		UCSCGeneLine[] l = new UCSCGeneLine[al.size()];
		al.toArray(l);
		return l;
	}
	
	/*TODO: replace with IntervalST*/
	public void loadTranscripts(){
		//main transcripts
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(transcriptSeqFile, 0);
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your transcript's coordinates are reversed. Check that each start is less than the stop.\n");
		chromGenes = reader.getChromSpecificGeneLines();
		numTranscripts = reader.getGeneLines().length;
		//load ccsg
		if (ccdsTranscriptFile != null){
			UCSCGeneModelTableReader ccsg = new UCSCGeneModelTableReader(ccdsTranscriptFile, 0);
			//check ordering
			if (ccsg.checkStartStopOrder() == false) Misc.printExit("\nOne of your ccds transcript's coordinates are reversed. Check that each start is less than the stop.\n");
			ccdsChromGenes = ccsg.getChromSpecificGeneLines();
		}
	}

	private void loadChromosomeData(String chrom) {
		try {
			//set new name
			workingChromosomeName = chrom;

			//find and load sequence
			ReferenceSequence p = fasta.getSequence(workingChromosomeName);
			if (p == null ) throw new IOException ("\n\nFailed to find or load a fasta sequence for '"+workingChromosomeName+"', aborting.\n");
			workingSequence = new String(p.getBases());
			workingSequence = workingSequence.toUpperCase();

			//fetch transcripts
			workingTranscripts = chromGenes.get(workingChromosomeName);
			if (workingTranscripts == null) {
				System.out.println("\tWARNING: no transcripts found for "+workingChromosomeName+" chromosome, skipping all associated vcf records.");
				workingGeneTree = null;
				return;
			}
			System.out.println("\tLoading: "+workingChromosomeName+"\tLen: "+workingSequence.length()+"\t Trans: "+workingTranscripts.length);

			//create interval tree, watch out for duplicates
			workingGeneTree = new IntervalST<ArrayList<UCSCGeneLine>>();
			for (UCSCGeneLine line : workingTranscripts){
				//the end is included in IntervalST so sub 1 from end
				int start = line.getTxStart();
				int stop = line.getTxEnd() -1;
				Interval1D it = new Interval1D(start, stop);
				ArrayList<UCSCGeneLine> al;
				if (workingGeneTree.contains(it)) al = workingGeneTree.get(it);
				else {
					al = new ArrayList<UCSCGeneLine>();
					workingGeneTree.put(it, al);
				}
				al.add(line);
			}
			if (workingGeneTree.size() ==0 ) workingGeneTree = null;
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nDo your chromosome names in the gene table and vcf file match? Watch for 'chr'.\n");
		}
	}




	

	private void loadNullScoreHistograms() throws IOException{
		
		if (histogramFile != null && histogramFile.exists()){
			System.out.println("\tFrom serialized object file");
			//Histogram[]{exonicHistogram5, exonicHistogram3, intronicHistogram5, intronicHistogram3, spliceHistogram5, spliceHistogram3};
			Histogram[] hist = (Histogram[]) IO.fetchObject(histogramFile);
			exonicHistogram5 = hist[0];
			exonicHistogram3 = hist[1];
			intronicHistogram5 = hist[2];
			intronicHistogram3 = hist[3];
			spliceHistogram5 = hist[4];
			spliceHistogram3 = hist[5];
			loadTypeHisto();
			return;
		}
		
		//for known splice
		spliceHistogram5 = new Histogram(-14, 0, 2500);
		spliceHistogram3 = new Histogram(-14, 0, 2500);
		exonicHistogram5 = new Histogram(0, 13, 2500);
		exonicHistogram3 = new Histogram(0, 17, 2500);
		intronicHistogram5 = new Histogram(0, 13, 2500);
		intronicHistogram3 = new Histogram(0, 17, 2500);
		
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		samSpliceReader = factory.open(samSpliceFile);
		
		//for each chrom of transcripts
		System.out.println("\tChr\tSeqLength\tTranscripts");
		for (String chr: chromGenes.keySet()){
			workingChromosomeName = chr;
			workingTranscripts = chromGenes.get(workingChromosomeName);
			//fetch seq
			ReferenceSequence p = fasta.getSequence(workingChromosomeName);
			if (p == null ) {
				System.out.println ("\t"+workingChromosomeName+"\tNo fasta, skipping!");
				continue;
			}
			workingSequence = new String(p.getBases());
			workingSequence = workingSequence.toUpperCase();
			System.out.println("\t"+workingChromosomeName+"\t"+workingSequence.length()+"\t"+workingTranscripts.length);
			
			//score known sj for background histogram
			if (ccdsChromGenes.containsKey(workingChromosomeName)) scoreKnownSplices();
			
			//mask known splices
			maskKnownSplices();
			
			//score plus and minus strands for novel splices in exons
			scoreAllExonsOrIntronsForSplices(true, true);
			scoreAllExonsOrIntronsForSplices(false, true);
			scoreAllExonsOrIntronsForSplices(true, false);
			scoreAllExonsOrIntronsForSplices(false, false);
		
		}
		samSpliceReader.close();
		
		//reset so these load
		workingSequence = "";
		workingChromosomeName = "";
		
		//print histograms
		if (printHistograms){
			System.out.println("\nHistograms of novel splice scores");
			System.out.println("5' Exon Splice:");
			exonicHistogram5.printScaledHistogram();
			System.out.println("\n3' Exon Splice:");
			exonicHistogram3.printScaledHistogram();
			System.out.println("5' Intron Splice:");
			intronicHistogram5.printScaledHistogram();
			System.out.println("\n3' Intron Splice:");
			intronicHistogram3.printScaledHistogram();
			System.out.println("\nHistograms of known splice scores");
			System.out.println("5' Splice Junctions:");
			spliceHistogram5.printScaledHistogram();
			System.out.println("\n3' Splice Junctions:");
			spliceHistogram3.printScaledHistogram();
			System.out.println("\nMean\tStandard deviation");
			System.out.println(exonicHistogram5.getStandardDeviation()+"\tExonic 5' Scores");
			System.out.println(exonicHistogram3.getStandardDeviation()+"\tExonic 3' Scores");
			System.out.println(intronicHistogram5.getStandardDeviation()+"\tIntronic 5' Scores");
			System.out.println(intronicHistogram3.getStandardDeviation()+"\tIntronic 3' Scores");
			System.out.println(spliceHistogram5.getStandardDeviation()+"\tKnown Splice Junction  5' Scores");
			System.out.println(spliceHistogram3.getStandardDeviation()+"\tKnown Splice Histogram 3' Scores");
		}
		System.out.println();
		//save them
		Histogram[] hist = new Histogram[]{exonicHistogram5, exonicHistogram3, intronicHistogram5, intronicHistogram3, spliceHistogram5, spliceHistogram3};
		histogramFile = new File (saveDirectory, "spliceHistograms.sjo");
		System.out.println("Saving histogram serialized object file. Use this to speed up subsequent matched species analysis: "+ histogramFile);
		IO.saveObject(histogramFile, hist);
		
		loadTypeHisto();
	}
	
	private void loadTypeHisto() {
		/*Three letter code: Gain or Damaged; 5' or 3' relative to gene not genomic; Exonic or Intronic or Splice
		 * e.g. G5E - gain 5' splice in exonic; D3S - damaged 3' splice in splice 
		 * G5E, G3E, G5I, G3I, D5S, D3S*/
		typeHist.clear();
		typeHist.put("G5E", exonicHistogram5);
		typeHist.put("G3E", exonicHistogram3);
		typeHist.put("G5I", intronicHistogram5);
		typeHist.put("G3I", intronicHistogram3);
		typeHist.put("D5S", spliceHistogram5);
		typeHist.put("D3S", spliceHistogram3);
	}

	private boolean checkCoverage(int position, String toLookFor, boolean leftRight){
		SAMRecordIterator i = samSpliceReader.queryOverlapping(workingChromosomeName, position, position);
		HashSet<String> names = new HashSet<String>();
		double numberOverlaps = 0;
		double numberCorrect = 0;

		ArrayList<String> allNames = new ArrayList<String>();
		while (i.hasNext()) {
			SAMRecord sam = i.next();
			allNames.add(sam.getReadName());
			//look for match to splice junction
			String sj = (String)sam.getAttribute("SJ");
			if (sj.contains(toLookFor)) {
				//does it end/ start at the junction
				if (leftRight){
					int alignEnd = sam.getAlignmentEnd()-1;
					if (position == alignEnd) continue;
				}
				else {
					// hmmm doesn't appear to be hit?
					int alignStart = sam.getAlignmentStart()-1;
					if (alignStart == position){
						Misc.printErrAndExit("SJ startIssue "+position+" "+toLookFor+"  as "+alignStart);
						continue;
					}
				}
				//check for ending read block (1based, end excluded; so add 2 to position)
				List<AlignmentBlock> blocks = sam.getAlignmentBlocks();
				//System.out.println(sam.getReadName());
				int numBlocks = blocks.size();
				if (numBlocks == 1) continue;
				numberOverlaps++;
				boolean add = false;
				for (int x=0; x< numBlocks; x++){
					AlignmentBlock ab = blocks.get(x);
					int pos;
					if (leftRight) pos = ab.getReferenceStart()+ab.getLength()-2;
					else pos = ab.getReferenceStart() -1;
					//correct position?
					if (pos == position){
						//is there another block downstream?
						if (leftRight){
							if (++x < numBlocks) add = true;
						}
						//is there a block upstream?
						else if (x != 0) add = true;
						break;
					}
				}
				if (add){
					numberCorrect++;
					names.add(sam.getReadName());
				}
			}
		}
		i.close();

		double fractionCorrect = numberCorrect/ numberOverlaps;
		if (names.size() < minimumSpliceJunctionCoverage || fractionCorrect < minimumFractionCorrectlySpliced) {
			//if (leftRight == false && names.size() < 10 && fractionCorrect < minimumFractionCorrectlySpliced) {
				//System.out.println("SJ "+position+" "+toLookFor+" "+fractionCorrect+" "+numRecords+" "+names);
			//}
			return false;
		}
		return true;
	}
	
	private void scoreKnownSplices() {
		//make array to hold info about presence of a splice base pair, default is false
		int lenMinOne = workingSequence.length() - 1;
		HashSet<String> scored5Plus = new HashSet<String>();
		HashSet<String> scored3Plus = new HashSet<String>();
		HashSet<String> scored5Minus = new HashSet<String>();
		HashSet<String> scored3Minus = new HashSet<String>();
		//for each transcript
		for (UCSCGeneLine g: ccdsChromGenes.get(workingChromosomeName)){
			ExonIntron[] introns = g.getIntrons();
			if (introns == null) continue;
			boolean plusStrand = g.getStrand().equals("+");
			//mask both 5' and 3' junctions
			int startJunc = 0;
			int endJunc = 0;
			for (int i=0; i< introns.length; i++){
				//plus strand 
				if (plusStrand){
					//5'
					startJunc = introns[i].getStart()-3;     
					endJunc = introns[i].getStart()+6;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Plus.contains(coor) == false){
						scored5Plus.add(coor);
						if (checkCoverage(introns[i].getStart()-1, "-"+introns[i].getStart()+"_", true)){
							String seq = workingSequence.substring(startJunc, endJunc);
							double score = score5.scoreSequenceWithChecks(seq);
							if (score != Double.MIN_VALUE) {
								spliceHistogram5.count(score);
								//System.out.println("+ 5\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
							}
							
						}
					}
					//3'
					startJunc = introns[i].getEnd()-20;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					coor = startJunc+"-"+endJunc;
					if (scored3Plus.contains(coor) == false){
						scored3Plus.add(coor);
						if (checkCoverage(introns[i].getEnd(), "_"+introns[i].getEnd(), false)){
							String seq = workingSequence.substring(startJunc, endJunc);
							double score = score3.scoreSequenceWithChecks(seq);
							if (score != Double.MIN_VALUE) spliceHistogram3.count(score);
							//System.out.println("+ 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);	
						}
					}
				}
				//minus strand
				else {
					//3'
					startJunc = introns[i].getStart()-3;
					endJunc = introns[i].getStart()+20;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored3Minus.contains(coor) == false){
						scored3Minus.add(coor);
						if (checkCoverage(introns[i].getStart()-1, "-"+introns[i].getStart()+"_", true)){
							String seq = workingSequence.substring(startJunc, endJunc);
							seq = Seq.reverseComplementDNA(seq);
							double score = score3.scoreSequenceWithChecks(seq);
							if (score != Double.MIN_VALUE) spliceHistogram3.count(score);
							//System.out.println("- 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
						}
					}
					
					//5'
					startJunc = introns[i].getEnd()-6;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Minus.contains(coor) == false){
						scored5Minus.add(coor);
						if (checkCoverage(introns[i].getEnd(), "_"+introns[i].getEnd(), false)){
							String seq = workingSequence.substring(startJunc, endJunc);
							seq = Seq.reverseComplementDNA(seq);
							double score = score5.scoreSequenceWithChecks(seq);
							if (score != Double.MIN_VALUE) {
								spliceHistogram5.count(score);
								//System.out.println("- 5\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
							}
								
						}
					}
				}
			}
		}
	}

	private void maskKnownSplices() {
		//make array to hold info about presence of a splice base pair, default is false
		boolean[] spliceBase = new boolean[workingSequence.length()];
		int lenMinOne = spliceBase.length-1;
		//for each transcript
		for (UCSCGeneLine g: workingTranscripts){
			ExonIntron[] introns = g.getIntrons();
			if (introns == null) continue;
			boolean plusStrand = g.getStrand().equals("+");
			//mask both 5' and 3' junctions
			int startJunc = 0;
			int endJunc = 0;
			for (int i=0; i< introns.length; i++){
				//plus strand 
				if (plusStrand){
					//5'
					startJunc = introns[i].getStart()-3;     
					endJunc = introns[i].getStart()+6;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
					//3'
					startJunc = introns[i].getEnd()-20;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
				}
				//minus strand
				else {
					//3'
					startJunc = introns[i].getStart()-3;
					endJunc = introns[i].getStart()+20;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
					//5'
					startJunc = introns[i].getEnd()-6;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					for (int j=startJunc; j< endJunc; j++) spliceBase[j] = true;
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

	private String[] fetchNMerSeqs(VCFRecord vcf, int nMer) {
		boolean reverseComplement = (workingTranscriptIsPlusStrand == false);
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

	public void loadAndModifyHeader(BufferedReader in) throws IOException{
		boolean addedInfo = false;
		boolean foundChrom = false;
		String line;
		while ((line=in.readLine()) != null){
			//comments
			if (line.startsWith("#")){
				//add info lines?
				if (addedInfo == false && line.startsWith("##INFO=")){
					addedInfo = true;
					vcfOut.println("##INFO=<ID=VCFSA,Number=.,Type=String,Description=\"USeq VCFSpliceAnnotator output. "
							+ "One or more splice junction (SJ) annotations delimited by a : each containing comma "
							+ "delimited SJType,SJPositon,GeneName,VCFAltSeq,SJ-10Log10(pval),SJRefScore,SJAltScore. "
							+ "SJTypes are a three char string with Gain or Damage; 3 or 5 prime; Intron, Exon or Splice "
							+ "(e.g. D3S, G5E, G3S). The SJPos is the interbase coordinate position of the damaged or "
							+ "gained SJ for SNPs.  For INDELS, it is the variant position.\">");
					
				}
				vcfOut.println(line);
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


	/**This method will process each argument and assign new variables
	 * @throws FileNotFoundException */
	public void processArgs(String[] args) throws FileNotFoundException{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		File indexedFasta = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'f': indexedFasta = new File(args[++i]); break;
					case 'm': spliceModelDirectory = new File(args[++i]); break;
					case 'u': transcriptSeqFile = new File(args[++i]); break;
					case 't': ccdsTranscriptFile = new File(args[++i]); break;
					case 'r': saveDirectory = new File(args[++i]); break;
					case 'j': samSpliceFile = new File(args[++i]); break;
					case 'k': minimumSpliceJunctionCoverage = Integer.parseInt(args[++i]); break;
					case 'e': scoreNovelExonJunctions = false; break;
					case 'i': scoreNovelIntronJunctions = false; break;
					case 's': scoreNovelSpliceJunctionsInSplice = false; break;
					case 'x': vcfExportCategory = Short.parseShort(args[++i]); break;
					case 'p': minPValue = Double.parseDouble(args[++i]); break;
					case 'a': min5Threshold = Double.parseDouble(args[++i]); break;
					case 'b': min3Threshold = Double.parseDouble(args[++i]); break;
					case 'c': min5DeltaThreshold = Double.parseDouble(args[++i]); break;
					case 'd': min3DeltaThreshold = Double.parseDouble(args[++i]); break;
					case 'h': histogramFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//checkfiles
		//Create fasta fetcher
		if (indexedFasta == null || indexedFasta.exists() == false)  Misc.printErrAndExit("\nError: cannot find indexed fasta file? -> "+ indexedFasta);
		fasta = new IndexedFastaSequenceFile(indexedFasta);
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);

		if (spliceModelDirectory == null || spliceModelDirectory.isDirectory() == false ) {
			Misc.printErrAndExit("\nError: please provide a path to a directory containing the splice models.\n");
		}
		if (transcriptSeqFile == null || transcriptSeqFile.canRead()== false){
			Misc.printErrAndExit("\nPlease enter a transcript table file in ucsc refflat format.\n");
		}
		//histogram
		if (histogramFile == null || histogramFile.exists() == false ){
			//look for bam and ccds
			if (ccdsTranscriptFile == null || samSpliceFile == null) Misc.printErrAndExit("\nPlease enter a histogram object file or an "
					+ "indexed bam file containing splice junction reads and a CCDS transcript table to generate background score histograms.\n");
		}
		//save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nPlease enter a directory to save the results.\n");
		saveDirectory.mkdirs();
		
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
				"**                            VCF Splice Annotator : Sept 2020                      **\n" +
				"**************************************************************************************\n" +
				"DON't USE, see VCFSpliceScanner!\n"+

				"**************************************************************************************\n");

	}
	
}
