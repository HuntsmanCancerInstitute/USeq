package edu.utah.seq.vcf.splice;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.mes.MaxEntScanScore3;
import edu.utah.seq.mes.MaxEntScanScore5;
import edu.utah.seq.vcf.xml.foundation.SimpleVcf;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.UCSCGeneLine;
import util.bio.seq.Seq;
import util.gen.Gzipper;
import util.gen.Misc;
import util.gen.Num;


/**This is a loader for splice junction info for use with the VCFSpliceAnnotator app.*/
public class SpliceAnnotationLoader implements Runnable {

	//fields
	private VCFSpliceScanner vsa;
	private MaxEntScanScore3 score3;
	private MaxEntScanScore5 score5;
	private boolean failed = false;
	private ArrayList<SimpleVcf> vcfRecords = new ArrayList<SimpleVcf>();
	private ArrayList<UCSCGeneLine> transcripts = new ArrayList<UCSCGeneLine>();
	private HashSet<String> effects = new HashSet<String>();
	private HashMap<String, SpliceJunction[]> scoredExons = new HashMap<String, SpliceJunction[]>();
	private HashMap<String, SpliceJunction[]> scoredIntrons = new HashMap<String, SpliceJunction[]>();
	private HashMap<String, SpliceJunction[]> scoredNovelSplices = new HashMap<String, SpliceJunction[]>();
	private SpliceJunction dummySpliceJunction = new SpliceJunction();
	private SpliceJunction[] knownSplice = null;
	private boolean removeInfoDropNonAffected = false;
	private IntervalST<ArrayList<UCSCGeneLine>> workingGeneTree;
	private boolean workingTranscriptIsPlusStrand;
	private String workingSequence;
	private int vcfExportCategory;
	private Gzipper vcfOut;
	
	//thresholds
	private double minNewAltScore;
	private double minNewScoreDelta;
	private double maxDamagedAltScore;
	private double minDamagedScoreDelta;
	
	//counters to upload to main thread
	private int numVariantsScanned = 0;
	private int numVariantsIntersectingTranscripts = 0;
	private int numVariantsEffectingSJs = 0;
	private int numVariantsIntersectingIntrons = 0;
	private int numVariantsIntersectingSJs = 0;
	private int numVariantsIntersectingExons = 0;
	private int numSJsLost = 0;
	private int numExonSJsGained = 0;
	private int numSpliceJunctionSJsGained = 0;
	private int numIntronSJsGained = 0;
	private HashSet<String> intersectingTranscriptNames = new HashSet<String>();
	

	//constructor
	public SpliceAnnotationLoader(VCFSpliceScanner vsa, File tempVcf) {
		try {
			this.vsa = vsa;

			//start up mes
			score5 = new MaxEntScanScore5(vsa.getSpliceModelDirectory());
			score3 = new MaxEntScanScore3(vsa.getSpliceModelDirectory());

			//pull thresholds and params
			this.minNewAltScore = vsa.getMinNewAltScore();
			this.minNewScoreDelta = vsa.getMinNewScoreDelta();
			this.maxDamagedAltScore = vsa.getMaxDamagedAltScore();
			this.minDamagedScoreDelta = vsa.getMinDamagedScoreDelta();
			this.vcfExportCategory = vsa.getVcfExportCategory();
			this.removeInfoDropNonAffected = vsa.isRemoveInfoDropNonAffected();

			//results writer, MUST close at end of loader life!!!!!!!!!
			vcfOut = new Gzipper (tempVcf);
			
		} catch (IOException e) {
			try { vcfOut.close(); } catch (IOException e1) {}
			e.printStackTrace();
			failed = true;
			tempVcf.delete();
		}
	}
	
	private void loadChromData() {
		//reference sequence
		workingSequence = vsa.getWorkingSequence();
		
		//create interval tree, watch out for duplicates
		workingGeneTree = new IntervalST<ArrayList<UCSCGeneLine>>();
		for (UCSCGeneLine line : vsa.getWorkingTranscripts()){
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
		
	}
	
	private void fetchIntersectingTranscripts(SimpleVcf vcf, ArrayList<UCSCGeneLine> al){
		int start = vcf.getPos();
		int stop = start+ vcf.getAlt().length();
		
		//end is included in search so subtract 1
		Iterable<Interval1D> it = workingGeneTree.searchAll(new Interval1D(start, stop-1));
		for (Interval1D x : it) {
			ArrayList<UCSCGeneLine> genes = workingGeneTree.get(x);
			al.addAll(genes);
			for (UCSCGeneLine l: genes)	intersectingTranscriptNames.add(l.getDisplayNameThenName()); 
		}
	}
	
	public void run() {	
		try {
			//another chromosome call so reset and pull fields
			loadChromData();
			//get next chunk of work
			while (vsa.loadRecords(vcfRecords)){ 
				
				//for each record
				for (SimpleVcf vcf: vcfRecords){
					String oriAlts = vcf.getAlt();
					String[] alts = Misc.COMMA.split(oriAlts);
					
					//for each alt
					for (int i=0; i< alts.length; i++){	
						numVariantsScanned++;
						vcf.setAlt(alts[i]);
						
						//load intersecting transcripts
						transcripts.clear();
						fetchIntersectingTranscripts(vcf, transcripts);
						if (transcripts.size() == 0) continue;
						numVariantsIntersectingTranscripts++;
						
						//for each transcript score effect
						for (UCSCGeneLine l: transcripts){
							SpliceHitThreaded sh = new SpliceHitThreaded (vcf, i, l);
							scoreVariantTranscript(sh); 
							if (sh.getAffectedSpliceJunctions() != null) effects.addAll(sh.getVcfEntries(vcfExportCategory));
						}
						
						//vcf is changing so clear scored seqs
						scoredExons.clear();
						scoredIntrons.clear();
						scoredNovelSplices.clear();
						knownSplice = null;
					}
					vcf.setAlt(oriAlts);
					
					//modify INFO field?
					if (effects.size() == 0 ) {
						if (removeInfoDropNonAffected == false) vcfOut.println(vcf.getOriginalRecord());
					}
					else {
						numVariantsEffectingSJs++;
						String e = consolidateEffects(effects);
						String[] fields = Misc.TAB.split(vcf.getOriginalRecord());
						if (removeInfoDropNonAffected) {
							//remove info from ID, FILTER, QUAL columns
							fields[2] = ".";
							fields[5] = ".";
							fields[6] = ".";
							//set INFO
							fields[7] = "VCFSS="+e;
							//wipe out any FORMAT or sample columns
							for (int i=8; i< fields.length; i++) fields[i] = null;
						}
						else fields[7] = fields[7]+";VCFSS="+e; 
						String f = Misc.stringArrayToStringSkipNulls(fields, "\t");
						vcfOut.println(f);
						//System.out.println("\n"+f);
						effects.clear();
					}
				}
				
				//clear and then load more
				vcfRecords.clear();
				
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem annotating splices" );
			e.printStackTrace();
		}
	}
	
	private static String consolidateEffects(HashSet<String> combineEff){
		//geneName+":"+sjType+","+sjPos+","+refSeq+","+altSeq+","+sjRefScore+","+sjAltScore+","+scoreDelta
		Iterator<String> it = combineEff.iterator();
		HashMap<String, ArrayList<String>> effName = new HashMap<String, ArrayList<String>>();
		while (it.hasNext()){
			String[] geneNamePlusEffs = Misc.COLON.split(it.next());
			//has this effect been seen?
			ArrayList<String> geneNames = effName.get(geneNamePlusEffs[1]);
			if (geneNames == null){
				geneNames = new ArrayList<String>();
				effName.put(geneNamePlusEffs[1], geneNames);
			}
			geneNames.add(geneNamePlusEffs[0]);
		}
		
		//build combine genes with each effect
		StringBuilder sb = new StringBuilder();
		
		int numEff = effName.keySet().size()-1;
		int counter = 0;
		for (String singleEffect: effName.keySet()){
			ArrayList<String> gn = effName.get(singleEffect);
			String geneNames = Misc.stringArrayListToString(gn, ",");
			sb.append(geneNames);
			sb.append(":");
			sb.append(singleEffect);
			if (counter++ != numEff) sb.append("&");
		}
		
		return sb.toString();
	}

	/**Returns null if no actionable info.  Otherwise with message.*/
	private void scoreVariantTranscript(SpliceHitThreaded hit) throws Exception{
		//set working vals
		workingTranscriptIsPlusStrand = hit.getTranscript().getStrand().equals("+");
		//check to see if variant is a snp or insertion, if so then a hit is mutually exclusive to that annotation class so can skip. For deletions, must scan all.
		boolean isSNPInsertion = true;
		if (hit.getVcf().isDeletion()) isSNPInsertion = false;
		//intronic and away from splice junctions
		if (vsa.isScoreNovelIntronJunctions()) scanIntronsForNovelSpliceJunctions(hit);
		if (isSNPInsertion && hit.getAffectedSpliceJunctions() != null) return;
		//exonic and away from splice junction
		if (vsa.isScoreNovelExonJunctions()) scanExonsForNovelSpliceJunctions(hit);
		if (isSNPInsertion && hit.getAffectedSpliceJunctions() != null) return;
		//scan hits to splices, gain and loss
		scanKnownSplices(hit);
	}
	
	private void scanIntronsForNovelSpliceJunctions (SpliceHitThreaded hit){
		SimpleVcf vcf = hit.getVcf();
		ExonIntron[] introns = hit.getTranscript().getIntrons();		
		if (introns != null) {
			int leftAdd = 6;
			int rightSub = 20;
			if (workingTranscriptIsPlusStrand == false){
				leftAdd = 20;
				rightSub = 6;
			}
			int varStart = vcf.getPos();
			int varEnd = varStart+ vcf.getRef().length();
			boolean intersection = false;
			for (int i=0; i< introns.length; i++){
				//mod coordinates of intron
				int intronStart = introns[i].getStart()+leftAdd;
				int intronEnd = introns[i].getEnd()-rightSub;
				if (intronEnd <= intronStart) continue;
				//if a hit then score 
				if (intersects(varStart, varEnd, intronStart, intronEnd) == false) continue;
				intersection = true;
				
				//already scored?
				String scoreKey = intronStart+"-"+intronEnd;  //the vcf info is the same
				SpliceJunction[] sjs = scoredIntrons.get(scoreKey);

				if (sjs == null){
					//scan for novel, reference and alternate			
					String[] fiveSeqs = fetchNMerSeqs(vcf, 9);
					SpliceJunction sj5 = scoreNewJunction(fiveSeqs, true, false, false, hit);
					String[] threeSeqs = fetchNMerSeqs(vcf, 23);
					SpliceJunction sj3 = scoreNewJunction(threeSeqs, false, false, false, hit);

					//save it
					sjs = new SpliceJunction[]{sj5, sj3};
					scoredIntrons.put(scoreKey, sjs);
				}

				//add em
				for (SpliceJunction sj: sjs) addSpliceJunction(sj, hit, false, false);
				
				//can only have one intersection if snp or insertion
				if (vcf.isSnv() || vcf.isInsertion()) break;
			}
			if (intersection) numVariantsIntersectingIntrons++;
			
		}
	}

	private void scanExonsForNovelSpliceJunctions (SpliceHitThreaded hit){
		SimpleVcf vcf = hit.getVcf();
		ExonIntron[] exons = hit.getTranscript().getExons();
		if (exons.length != 0) {
			int varStart = vcf.getPos();
			int varEnd = varStart+ vcf.getRef().length();
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
				
				//already scored?
				//String scoreKey = vcf.getChr()+":"+start+"-"+end+" "+vcf.getPos()+" "+vcf.getRef()+" "+vcf.getAlt();
				String scoreKey = start+"-"+end;  //the vcf info is the same
				SpliceJunction[] sjs = scoredExons.get(scoreKey);
				if (sjs == null){
					//OK need to scan
					
					//scan for novel, reference and alternate			
					String[] fiveSeqs = fetchNMerSeqs(vcf, 9);
					SpliceJunction sj5 = scoreNewJunction(fiveSeqs, true, true, false, hit);
					String[] threeSeqs = fetchNMerSeqs(vcf, 23);
					SpliceJunction sj3 = scoreNewJunction(threeSeqs, false, true, false, hit);

					//save it
					sjs = new SpliceJunction[]{sj5, sj3};
					scoredExons.put(scoreKey, sjs);
				}

				//add em
				for (SpliceJunction sj: sjs) addSpliceJunction(sj, hit, true, false);
				
				//can only have one intersection if snp or insertion so skip remaining exons
				if (vcf.isSnv() || vcf.isInsertion()) break;
			}
			if (intersection) numVariantsIntersectingExons++;
		}
	}
	
	private void addSpliceJunction(SpliceJunction spliceJunction, SpliceHitThreaded spliceHit, boolean exonic, boolean splice){
		if (spliceJunction == null) return;
		
		spliceHit.saveSpliceJunction(spliceJunction);
		//increment counters
		if (exonic) numExonSJsGained++;
		else if (splice) numSpliceJunctionSJsGained++;
		else numIntronSJsGained++;
	}
	
	private String[] fetchNMerSeqs(SimpleVcf vcf, int nMer) {
		boolean reverseComplement = (workingTranscriptIsPlusStrand == false);
		String refSeq = null;
		String altSeq = null;
		//set coordinates
		int pos = vcf.getPos();
		int start = pos - nMer +1;
		int stop = 0;
			stop = vcf.getPos()+ nMer + vcf.getRef().length() -1;
			refSeq = workingSequence.substring(start, stop);
			altSeq = workingSequence.substring(start, pos) + vcf.getAlt() + workingSequence.substring(pos+vcf.getRef().length(), stop);
		//reverse comp it?
		if (reverseComplement){
			refSeq = Seq.reverseComplementDNA(refSeq);
			altSeq = Seq.reverseComplementDNA(altSeq);
		}
		return new String[]{refSeq, altSeq};
	}

	
	private boolean intersects(int varStart, int varEnd, int annoStart, int annoEnd){
		if (varStart >= annoEnd) return false;
		if (varEnd <= annoStart) return false;
		return true;
	}
	
	
	private void scanKnownSplices (SpliceHitThreaded hit) throws Exception{
		SimpleVcf vcf = hit.getVcf();
		ExonIntron[] introns = hit.getTranscript().getIntrons();

		int varStart = vcf.getPos();
		int varEnd = varStart+ vcf.getRef().length();
		boolean intersectsSplice = false;

		if (introns!= null) {
			//find junction
			int junc = 0;
			int startJunc = 0;
			int endJunc = 0;
			for (int i=0; i< introns.length; i++){

				int intronStart = introns[i].getStart();
				int intronEnd = introns[i].getEnd();

				//already scored?
				String scoreKey = intronStart+"-"+intronEnd+"-"+workingTranscriptIsPlusStrand;  //the vcf info is the same
				SpliceJunction[] sjs = scoredNovelSplices.get(scoreKey);

				if (sjs == null){
					SpliceJunction sj3 = null;
					SpliceJunction sj5 = null;

					//plus strand 
					if (workingTranscriptIsPlusStrand){
						//5'?
						startJunc = intronStart-3;     
						endJunc = intronStart+6;
						if (intersects(varStart, varEnd, startJunc, endJunc)) {
							junc = intronStart;
							sj5 = scoreLossOfKnownJunction(junc, vcf, true);
							intersectsSplice = true;
						}
						//3'?
						startJunc = intronEnd-20;
						endJunc = intronEnd+3;
						if (intersects(varStart, varEnd, startJunc, endJunc)) {
							junc = intronEnd;
							sj3 = scoreLossOfKnownJunction(junc, vcf, false);
							intersectsSplice = true;
						}
					}
					//minus strand
					else {
						//3'?
						startJunc = intronStart-3;
						endJunc = intronStart+20;
						if (intersects(varStart, varEnd, startJunc, endJunc)) {
							junc = intronStart;
							sj3 = scoreLossOfKnownJunction(junc, vcf, false);
							intersectsSplice = true;
						}
						//5'?
						startJunc = intronEnd-6;
						endJunc = intronEnd+3;
						if (intersects(varStart, varEnd, startJunc, endJunc)) {
							junc = intronEnd;
							sj5 = scoreLossOfKnownJunction(junc, vcf, true);
							intersectsSplice = true;
						}
					}
					//need to save intersect splice boolean so add a dummy sj if it is true
					SpliceJunction dummy = null;
					if (intersectsSplice) dummy = dummySpliceJunction;
					sjs = new SpliceJunction[] {sj3, sj5, dummy};
					scoredNovelSplices.put(scoreKey, sjs);
				}
				//read in if it intersects the splice boolean
				else if (sjs[2] != null) intersectsSplice = true;

				//damaged?
				if (sjs[0] != null) {
					hit.saveSpliceJunction(sjs[0]);
					numSJsLost++;
				}
				if (sjs[1] != null) {
					hit.saveSpliceJunction(sjs[1]);
					numSJsLost++;
				}
			}
			
			//did any vars hit the splice junction?
			if (intersectsSplice){
				numVariantsIntersectingSJs++;
				//now scan for novel new splice in splice site
				if (vsa.isScoreNovelSpliceJunctionsInSplice()){
				
					if (knownSplice == null){
						String[] fiveSeqs = fetchNMerSeqs(vcf, 9);
						SpliceJunction sj5 = scoreNewJunction(fiveSeqs, true, false, true, hit);
						String[] threeSeqs = fetchNMerSeqs(vcf, 23);	
						SpliceJunction sj3 = scoreNewJunction(threeSeqs, false, false, true, hit);
						knownSplice = new SpliceJunction[]{sj5, sj3};
					}

					//add em, ok if null
					for (SpliceJunction sj: knownSplice) addSpliceJunction(sj, hit, false, true);
				}
			}
		}
	}	
	
	
	

	

	/**Returns three sequences, the reference splice junction, the modified splice junction where the position has been fixed or follows the original.
	 * {refSeq, fixedSeq, shiftedSeq} */
	private String[] fetchModifiedSpliceSequences(SimpleVcf vcf, int junctionPosition, boolean fetch5Prime) throws Exception{
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
		
		int position = vcf.getPos();
		String ref = vcf.getRef();
		String alt = vcf.getAlt();
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
	private SpliceJunction scoreLossOfKnownJunction(int junctionPosition, SimpleVcf vcf, boolean score5Junction) throws Exception{
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
		double delta = 0;
		String altScoreMaxSeq = null;

		if (score5Junction){
			//calculate ref
			refScore = score5.scoreSequenceNoChecks(rfs[0]);
			
			//calculate max alt
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
			//check alt
			if (altScoreMax <= maxDamagedAltScore){
				//calc delta
				if (altScoreMax < 0) delta = refScore;
				else delta = refScore- altScoreMax;
				if (delta >= minDamagedScoreDelta) sj =  new SpliceJunction('D', '5', 'S', junctionPosition);
			}
		}
		else {
			//calculate ref
			refScore = score3.scoreSequenceNoChecks(rfs[0]);

			//calc max alts
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
			//check alt
			if (altScoreMax <= maxDamagedAltScore){
				//calc delta
				if (altScoreMax < 0) delta = refScore;
				else delta = refScore- altScoreMax;
				if (delta >= minDamagedScoreDelta) sj =  new SpliceJunction('D', '3', 'S', junctionPosition);
			}
		}
		//damaged?
		if (sj != null){
			//add info to sj
			sj.setReferenceSequence(rfs[0]);
			sj.setReferenceScore(refScore);
			sj.setAlternateSequence(altScoreMaxSeq);
			sj.setAlternateScore(altScoreMax);
			sj.setScoreDelta(delta);
			return sj;
		}
		//nope
		return null;
	}

	
	private SpliceJunction scoreNewJunction(String[] refAltSeqsToScan, boolean score5Junction, boolean exonic, boolean splice, SpliceHitThreaded spliceHit) {
		double[] r;
		double[] a;

		SpliceJunction spliceJunction;
		if (score5Junction){
			r = score5.scanSequence(refAltSeqsToScan[0], -1000);
			a = score5.scanSequence(refAltSeqsToScan[1], -1000);
			if (exonic) spliceJunction = new SpliceJunction('G', '5', 'E');
			else if (splice) spliceJunction = new SpliceJunction('G', '5', 'S');
			else spliceJunction = new SpliceJunction('G', '5', 'I');
		}
		else {
			r = score3.scanSequence(refAltSeqsToScan[0], -1000);
			a = score3.scanSequence(refAltSeqsToScan[1], -1000);
			if (exonic) spliceJunction = new SpliceJunction('G', '3', 'E');
			else if (splice) spliceJunction = new SpliceJunction('G', '3', 'S');
			else spliceJunction = new SpliceJunction('G', '3', 'I');
		}
		
		//compare scores
		//snp, direct 1:1 comparison
		if (r.length == a.length) {

			//scan all for maxDelta
			double maxDelta = 0;
			int indexMaxDelta = -1;
			for (int i=0; i< r.length; i++){
				//does alt exceed min threshold, if it's -1000 it'll be skipped
				if (a[i] < minNewAltScore) continue;
				//skip due to bad base in reference?
				if (r[i] == -1000) continue;
				//is the ref less than min threshold
				if (r[i] > minNewAltScore) continue;
				double testR = r[i];
				if (testR < 0) testR = 0;
				double delta = a[i] - testR;
				if (delta > maxDelta) {
					maxDelta = delta;
					indexMaxDelta = i;
				}
			}
			//anything pass min delta?
			if (maxDelta >= minNewScoreDelta){

				//good so save info in spliceJunction
				spliceJunction.setReferenceScore(r[indexMaxDelta]);
				spliceJunction.setAlternateScore(a[indexMaxDelta]);
				spliceJunction.setScoreDelta(maxDelta);
				int posOfZeroSeq;
				int relPosSpliceInSeq;
				if (score5Junction) {
					spliceJunction.setReferenceSequence(refAltSeqsToScan[0].substring(indexMaxDelta, indexMaxDelta+9));
					spliceJunction.setAlternateSequence(refAltSeqsToScan[1].substring(indexMaxDelta,indexMaxDelta+9));
					posOfZeroSeq = spliceHit.getVcf().getPos() - 8;
					relPosSpliceInSeq = indexMaxDelta+3;
				}
				else {
					spliceJunction.setReferenceSequence(refAltSeqsToScan[0].substring(indexMaxDelta, indexMaxDelta+23));
					spliceJunction.setAlternateSequence(refAltSeqsToScan[1].substring(indexMaxDelta,indexMaxDelta+23));
					posOfZeroSeq = spliceHit.getVcf().getPos() -22;
					relPosSpliceInSeq = indexMaxDelta+20;
				}
				//if neg strand then flip
				if (workingTranscriptIsPlusStrand == false){
					relPosSpliceInSeq = refAltSeqsToScan[0].length() - relPosSpliceInSeq;
				}
				spliceJunction.setPosition(posOfZeroSeq + relPosSpliceInSeq);
			}
			//nope, nada passes so return
			else return null;
		}

		//indel, this is a partial direct comparison, not necessarily 1:1.
		else{
			//insert zeros at the center position of the smaller array
			if (a.length > r.length) r = Num.expand(r, a.length);
			else a = Num.expand(a, r.length);
			
			//for each, find max delta
			int bestIndex = -1;
			double maxDelta = -1;
			for (int i=0; i< a.length; i++){
				//skip due to bad base in reference?
				if (r[i] == -1000) continue;
				//is the alt score big enough?
				if (a[i] < minNewAltScore) continue;
				//is the ref small enought?
				if (r[i] > minNewAltScore) continue;
				//calc diff
				double ref = r[i];
				if (ref < 0) ref = 0;
				double delta = a[i] - ref;
				//delta big enough?
				if (delta > maxDelta){
					bestIndex = i;
					maxDelta = delta;
				}
			}

			
			if (maxDelta < minNewScoreDelta) return null;
			
			//OK it's good, save info
			spliceJunction.setReferenceScore(r[bestIndex]);
			spliceJunction.setReferenceSequence(refAltSeqsToScan[0]);
			spliceJunction.setAlternateScore(a[bestIndex]);
			spliceJunction.setAlternateSequence(refAltSeqsToScan[1]);
			spliceJunction.setPosition(spliceHit.getVcf().getPos());
			spliceJunction.setScoreDelta(maxDelta);
		}
		return spliceJunction;
	}
	
	public boolean isFailed() {
		return failed;
	}

	public Gzipper getVcfOut() {
		return vcfOut;
	}

	public int getNumVariantsScanned() {
		return numVariantsScanned;
	}

	public int getNumVariantsIntersectingTranscripts() {
		return numVariantsIntersectingTranscripts;
	}

	public int getNumVariantsIntersectingIntrons() {
		return numVariantsIntersectingIntrons;
	}

	public int getNumVariantsIntersectingSJs() {
		return numVariantsIntersectingSJs;
	}

	public int getNumVariantsIntersectingExons() {
		return numVariantsIntersectingExons;
	}

	public int getNumSJsLost() {
		return numSJsLost;
	}

	public int getNumExonSJsGained() {
		return numExonSJsGained;
	}

	public int getNumSpliceJunctionSJsGained() {
		return numSpliceJunctionSJsGained;
	}

	public int getNumIntronSJsGained() {
		return numIntronSJsGained;
	}

	public HashSet<String> getIntersectingTranscriptNames() {
		return intersectingTranscriptNames;
	}

	public int getNumVariantsEffectingSJs() {
		return numVariantsEffectingSJs;
	}

}
