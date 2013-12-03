package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.mes.*;
import util.bio.annotation.ExonIntron;
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
		loadTranscripts();

		//for each vcf file
		System.out.println("Processing:");
		for (int i=0; i< vcfFiles.length; i++){
			System.out.println("\t"+vcfFiles[i]);
			annotateVCFWithSplices(vcfFiles[i]);
		}
		
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
		//check gene name is unique
		if (reader.uniqueGeneNames() == false) Misc.printExit("\nDuplicate transcript names were found in your gene file, these must be unique.\n");
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
				//parse record
				VCFRecord vcf = new VCFRecord(line, parser, true, true);
				vcf.appendChr();
				//load new sequence and transcripts?
				if (vcf.getChromosome().equals(workingChromosomeName) == false) {
					loadChromosomeData(vcf);
					System.out.print(workingChromosomeName+" ");
				}
				//check if any transcripts were found, otherwise skip, might want to write out this vcf record somewhere?
				if (workingTranscripts == null) continue;

				//for each alternate allele
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
			System.out.println();
			out.close();
			in.close();
		}catch (Exception e) {
			System.err.println("Error -> "+e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

	private String scanIntronsWithSNP (VCFRecord vcf, UCSCGeneLine transcript, boolean isPlusStrand){
		ExonIntron[] introns = transcript.getIntrons();
		int varPos = vcf.getPosition();
		String results = "noIntersection";
		if (introns != null) {
			int leftAdd = 6;
			int rightSub = 20;
			if (isPlusStrand == false){
				leftAdd = 20;
				rightSub = 6;
			}
			for (int i=0; i< introns.length; i++){
				int start = introns[i].getStart()+leftAdd;
				int end = introns[i].getEnd()-rightSub;
				if (end <= start) continue;
				//if a hit then score 
				if (varPos >= start && varPos < end){
					numVariantsIntersectingIntrons++;
					//scan for novel, reference and alternate
					String[] fiveSeqs = fetchNMerSeqs(vcf, 9, isPlusStrand == false);
					String fiveRes = scoreNewJunction_SNP(fiveSeqs, true);
					String[] threeSeqs = fetchNMerSeqs(vcf, 23, isPlusStrand == false);
					String threeRes = scoreNewJunction_SNP(threeSeqs, false);
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
						//can only be one hit to an intron with a snp.
						results = alerts.toString();
					}
					else results = "noChange";
				}
			}
		}
		return results;
	}

	private String scanExonsWithSNP (VCFRecord vcf, UCSCGeneLine transcript, boolean isPlusStrand){
		ExonIntron[] exons = transcript.getExons();
		int varPos = vcf.getPosition();
		String results = "noIntersection";
		if (exons.length != 0) {
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
				if (varPos>= start && varPos < end){
					numVariantsIntersectingExons++;
					//scan for novel, reference and alternate
					String[] fiveSeqs = fetchNMerSeqs(vcf, 9, isPlusStrand == false);
					String fiveRes = scoreNewJunction_SNP(fiveSeqs, true);
					String[] threeSeqs = fetchNMerSeqs(vcf, 23, isPlusStrand == false);
					String threeRes = scoreNewJunction_SNP(threeSeqs, false);
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
						//can only be one exon hit with a snp
						results = alerts.toString();
					}
					else results = "noChange";
				}
			}
		}
		return results;
	}

	private String scanSplicesWithSNP (VCFRecord vcf, UCSCGeneLine transcript, boolean isPlusStrand){
		ExonIntron[] introns = transcript.getIntrons();
		int varPos = vcf.getPosition();

		if (introns!= null) {
			//find junction
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
					if (varPos >= startJunc && varPos < endJunc) {
						String spliceResult = scoreLossOfKnownJunction_SNP(startJunc, endJunc, vcf, true, isPlusStrand);
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
						if (varPos >= startJunc && varPos < endJunc) {
							String spliceResult = scoreLossOfKnownJunction_SNP(startJunc, endJunc, vcf, false, isPlusStrand);
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
					if (varPos >= startJunc && varPos < endJunc) {
						String spliceResult = scoreLossOfKnownJunction_SNP(startJunc, endJunc, vcf, false, isPlusStrand);
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
						if (varPos >= startJunc && varPos < endJunc) {
							String spliceResult = scoreLossOfKnownJunction_SNP(startJunc, endJunc, vcf, true, isPlusStrand);
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
				String fiveRes = scoreNewJunction_SNP(fiveSeqs, true);
				String[] threeSeqs = fetchNMerSeqs(vcf, 23, isPlusStrand == false);
				String threeRes = scoreNewJunction_SNP(threeSeqs, false);
				if (fiveRes != null || threeRes != null){
					sb.append ("Novel juction found in splice junction "+startJunc+"-"+endJunc+"\n");
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

		//is it a snp?
		if (vcf.isSNP()){
			//intronic and away from splice junctions
			if (scoreNovelIntronJunctions){
				String intronResults = scanIntronsWithSNP(vcf, transcript, isPlusStrand);
				if (intronResults.equals("noChange")) return null;
				if (intronResults.equals("noIntersection") == false) return intronResults;
			}

			//exonic and away from splice junctions
			if (scoreNovelExonJunctions){
				String exonResults = scanExonsWithSNP(vcf, transcript, isPlusStrand);
				if (exonResults.equals("noChange")) return null;
				if (exonResults.equals("noIntersection") == false) return exonResults;
			}

			//must be a splice hit, scan for loss or gain of novel
			String spliceResults = scanSplicesWithSNP(vcf, transcript, isPlusStrand);
			if (spliceResults.equals("noChange")) return null;
			if (spliceResults.equals("noIntersection") ){
				if (scoreNovelIntronJunctions && scoreNovelExonJunctions)  Misc.printErrAndExit("Hmm variant didn't land on an annotation, something is wrong!\n"+transcript.toStringAll()+"\n"+ vcf.toString());
				return null;
			}
			return spliceResults;

		}

		//indel handling
		return null;
	}

	/**
	 * 1) Need to test if it is a junction.  
	 * So get a distribution of non junction scores (all introns from merged gene table) and calc pvalue.
	 * 
	 * 2) Need to test it is not a junction. 
	 * So get a distribution of junction scores (all from transcript table 
	 * */


	private String scoreLossOfKnownJunction_SNP(int startJunc, int endJunc, VCFRecord vcf, boolean score5Junction, boolean isPlusStrand){
		String refSeq = workingSequence.substring(startJunc, endJunc);
		//watch out for non GATC bases
		Matcher mat = MaxEntScanScore5.NonGATC.matcher(refSeq);
		if (mat.find() ) return null;
		int varPos = vcf.getPosition();
		String altSeq = workingSequence.substring(startJunc, varPos) + vcf.getAlternate()[0] + workingSequence.substring(varPos+1, endJunc);

		//reverse comp?
		if (isPlusStrand == false) {
			refSeq = Seq.reverseComplementDNA(refSeq);
			altSeq = Seq.reverseComplementDNA(altSeq);
		}

		boolean lostKnown = false;
		double refScore;
		double altScore;
		double d;
		String name;
		if (score5Junction){
			refScore = score5.scoreSequence(refSeq);
			altScore = score5.scoreSequence(altSeq);			
			//does ref meet threshold?
			if (refScore < min5Threshold) return null;
			//does the delta
			d= refScore- altScore;
			if (d >= min5DeltaThreshold) lostKnown = true;
			name = "5'";
		}
		else {
			refScore = score3.scoreSequence(refSeq);
			altScore = score3.scoreSequence(altSeq);
			//does ref meet threshold?
			if (refScore < min3Threshold) return null;
			//does the delta
			d= refScore- altScore;
			if (d >= min3DeltaThreshold) lostKnown = true;
			name = "3'";
		}

		if (lostKnown){
			StringBuilder sb = new StringBuilder();
			sb.append (name+" splice junction ("+startJunc+"-"+endJunc+") possibly damaged\n");
			sb.append ("Ref\t"+refSeq+"\t"+refScore+"\n");
			sb.append ("Alt\t"+altSeq+"\t"+altScore+"\n");
			sb.append ("Delta\t"+d+"\n");
			return sb.toString();
		}

		return null;
	}

	private String scoreNewJunction_SNP(String[] refAltSeqsToScan, boolean score5Junction) {
		double[] r;
		double[] a;
		double minThres;
		double minDelta;
		String name;
		if (score5Junction){
			r = score5.scanSequence(refAltSeqsToScan[0]);
			a = score5.scanSequence(refAltSeqsToScan[1]);
			minThres = min5Threshold;
			minDelta = min5DeltaThreshold;
			name = "5'";
		}
		else {
			r = score3.scanSequence(refAltSeqsToScan[0]);
			a = score3.scanSequence(refAltSeqsToScan[1]);
			minThres = min3Threshold;
			minDelta = min3DeltaThreshold;
			name = "3'";
		}

		StringBuilder sb = new StringBuilder();
		for (int i=0; i< r.length; i++){
			//does alt exceed min threshold
			if (a[i] < minThres) continue;
			double testR = r[i];
			if (testR < 0) testR = 0;
			double delta = a[i] - testR;
			if (delta >= minDelta) {
				sb.append(name+" "+Num.formatNumber(r[i],1)+" -> "+Num.formatNumber(a[i],1)+" index: "+i+" delta: "+Num.formatNumber(delta, 1)+"\n");
			}
		}
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

	private void loadChromosomeData(VCFRecord vcf) throws IOException {
		//set new name
		workingChromosomeName = vcf.getChromosome();
		//find and load sequence
		File seqFile = chromFile.get(workingChromosomeName);
		if (seqFile == null || seqFile.canRead() == false) throw new IOException ("\n\nFailed to find or load a fasta sequence for '"+workingChromosomeName+"', aborting.\n");
		workingSequence = new MultiFastaParser(seqFile).getSeqs()[0];
		workingSequence.toUpperCase();
		//fetch transcripts
		workingTranscripts = chromGenes.get(workingChromosomeName);
		if (workingTranscripts == null) System.out.println("\n\nWARNING: no transcripts found for this chromosome, skipping all associated vcf records.");
	}

	private String[] fetchNMerSeqs(VCFRecord vcf, int nMer, boolean reverseComplement) {
		String refSeq = null;
		String altSeq = null;
		//set coordinates
		int pos = vcf.getPosition();
		int start = pos - nMer +1;
		int stop = 0;
		//if snp
		if (vcf.isSNP()) {
			//System.out.println ("\nSNP\t"+vcf.getReference()+" -> "+vcf.getAlternate()[0]);
			stop = vcf.getPosition() + nMer;
			refSeq = workingSequence.substring(start, stop);
			//System.out.println ("RefSeq\t"+refSeq);
			//System.out.println ("AltSeq\t"+sequence.substring(start, pos) +" "+ vcf.getAlternate()[0] +" "+ sequence.substring(pos+1, stop));
			altSeq = workingSequence.substring(start, pos) + vcf.getAlternate()[0] + workingSequence.substring(pos+1, stop);
		}
		//insertion, alt > ref (1base)
		else if (vcf.isInsertion()) {
			//System.out.println ("\nIns\t"+vcf.getReference()+" -> "+vcf.getAlternate()[0]);
			stop = vcf.getPosition()+ nMer;
			refSeq = workingSequence.substring(start, stop);
			//System.out.println ("RefSeq\t"+refSeq);
			//System.out.println ("AltSeq\t"+sequence.substring(start, pos) +" "+ vcf.getAlternate()[0] +" "+ sequence.substring(pos+1, stop));
			altSeq = workingSequence.substring(start, pos) + vcf.getAlternate()[0] + workingSequence.substring(pos+1, stop);
		}
		//deletion, ref > alt, CCT	C
		else if (vcf.isDeletion()) {
			//System.out.println ("\nDel\t"+vcf.getReference()+" -> "+vcf.getAlternate()[0]);
			stop = vcf.getPosition()+ nMer + vcf.getReference().length() -1;
			refSeq = workingSequence.substring(start, stop);
			//System.out.println ("RefSeq\t"+refSeq);
			//System.out.println ("AltSeq\t"+sequence.substring(start, pos) +" "+ vcf.getAlternate()[0] +" "+ sequence.substring(pos+vcf.getReference().length(), stop));
			altSeq = workingSequence.substring(start, pos) + vcf.getAlternate()[0] + workingSequence.substring(pos+vcf.getReference().length(), stop);
		}
		//?
		else {
			Misc.printErrAndExit("\nHmm, shouldn't hit this "+vcf.getOriginalRecord());
			return null;
		}
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
				"WARNING: beta, only supports SNPs.\n"+
				"\n"+
				"Scores variants for changes in splicing using the MaxEntScan algorithms. Known splice\n"+
				"acceptors and donors are scored for loss of a junction.  Exonic, intronic, and splice\n"+
				"bases are scanned for novel junctions. See Yeo and Burge 2004,\n"+
				"http://www.ncbi.nlm.nih.gov/pubmed/15285897 for details.\n\n" +

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
				"-c Minimum difference for loss or gain of a 3' splice junction, defaults to 5.\n"+
				"\n"+

				"Example: java -Xmx10G -jar ~/USeq/Apps/VCFSpliceAnnotator -f ~/Hg19/Fa/ -v ~/exm2.vcf\n"+
				"       -m ~/USeq/Documentation/splicemodels -u ~/Hg19/hg19EnsTran.ucsc.zip \n\n"+

				"**************************************************************************************\n");

	}
	
}
