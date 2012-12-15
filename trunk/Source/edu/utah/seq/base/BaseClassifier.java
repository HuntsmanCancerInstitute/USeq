package edu.utah.seq.base;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.util.FastMath;

import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import net.sf.samtools.*;
import net.sf.samtools.SAMRecord.SAMTagAndValue;

import util.bio.annotation.Bed;
import util.gen.*;

/**
 * Copyright 2012 Huntsman Cancer Institute.  Use of this software requires a license from HCI.
 * @author Nix
 * */
public class BaseClassifier {

	//user fields
	private File[] bamFiles;
	private File bedFile;
	private int nMerSize = 5;
	private byte minimumBaseQuality = 20;
	private byte numberFlankingBaseQualitiesToCheck = 3;
	private float minimumMappingQuality = 150;
	private float maximumAlignmentScore = 0;
	private String cigar = "101M";
	private float minimumCorrelation= 0.85f;
	private float minimumCorrelationDifference = 0.01f;
	private double minimumLikelihood= 0.95;
	private double minimumLikelihoodDifference = 0.5;
	private int maxAdjacentDistance = 3;
	private boolean printOnlyChanged = true;
	private boolean useLikelihood = true;
	private boolean printLikelihoods = false;
	private boolean setRecalledBasesToN = false;

	//internal

	private int nMerAdder;
	private HashMap<String, HashMap<String, NMer>> firstReadMatrixLookup = new HashMap<String, HashMap<String, NMer>>(1000000);
	private HashMap<String, HashMap<String, NMer>> secondReadMatrixLookup = new HashMap<String, HashMap<String, NMer>>(1000000);
	private HashMap<String, SAMRecord> changedAlignments = new HashMap<String, SAMRecord>();

	private Bed[] variants;
	private VariantAlignment[][] variantAlignments;
	private SAMFileReader[] bamReaders;

	private long numberModels = 0;
	private long numReads = 0;
	private long numFlaggedReads = 0;
	private long numFailedCigar = 0;
	private long numCtrlChromosomes = 0;
	private long numPassingReads = 0;
	private long numFailedScore = 0;
	private long numVarAlignFailingMinCor = 0;
	private long numVarAlignFailingMinCorDiff = 0;
	private long numVarAlign = 0;
	private long numVarAlignChanged2N = 0;
	private long numVarAlignChanged2Novel = 0;

	private WeightedPearsonCorrelation wpc = new WeightedPearsonCorrelation();

	//constructor
	public BaseClassifier(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//make +/-
		nMerAdder = (int)((float)nMerSize/2.0f);

		//load variants
		System.out.println("Fetching variant associated alignments...");
		variants = Bed.parseFile(bedFile, 0, 0);
		Arrays.sort(variants);
		variantAlignments = new VariantAlignment[variants.length][];

		//make alignment readers
		bamReaders = new SAMFileReader[bamFiles.length];
		for (int i=0; i< bamReaders.length; i++) {
			bamReaders[i] = new SAMFileReader(bamFiles[i]);
			bamReaders[i].enableIndexMemoryMapping(false);
		}

		//load alignments for variants
		loadVariantAlignments();

		System.out.print("Loading intensity models");
		//scan alignments for building lookup matricies
		scanVariantAlignments();

		//load matricies with high quality values
		loadModels();
		
		//print good alignment stats
		printAlignmentStats();
		
		//print observations per model
		//printModels();
		
		System.out.println("Calling variant alignments with models...");
		scoreVariants();


		System.out.println("\nPrinting variants with alignments...");
		printVariants();

		//print stats
		printBaseClassifierStats();
		
		//write out new bam file?
		writeBAM();

		//close alignment readers
		for (int i=0; i< bamReaders.length; i++) bamReaders[i].close();


		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}

	//methods
	
	public void writeBAM(){
		if (numVarAlignChanged2N ==0 && numVarAlignChanged2Novel == 0) {
			System.out.println("\nNo changes made to bam files so not writing out new file.\n");
			return;
		}
		
		System.out.print("\nWriting modified BAM file(s)");
		int counter = 0;

		//for each bam reader
		for (int i=0; i< bamFiles.length; i++){
			//make writer
			String name = Misc.removeExtension(bamFiles[i].getName()) +"_Ml"+Num.formatNumber(minimumLikelihood, 1)+"Mld"+Num.formatNumber(minimumLikelihoodDifference, 1);
			 
			File bam = new File(bamFiles[i].getParentFile(), name +".bam");
			SAMFileWriterFactory f = new SAMFileWriterFactory();
			f.setCreateIndex(true);
			SAMFileWriter out = f.makeBAMWriter(bamReaders[i].getFileHeader(), true, bam);
			
			//for each sam record
			SAMRecordIterator it = bamReaders[i].iterator();
			
			while (it.hasNext()) {
				//counter
				if (counter++ > 1000000) {
					counter = 0;
					System.out.print(".");
				}
				//load record 
				SAMRecord sam = it.next();
				
				
				//has it been changed?
				String samName = "1";
				if (sam.getSecondOfPairFlag()) samName = "2";
				samName = sam.getReadName()+samName;
				if (changedAlignments.containsKey(samName)) sam = changedAlignments.get(samName);
				
				//strip intensities
				sam.setAttribute("IA", null);
				sam.setAttribute("IC", null);
				sam.setAttribute("IG", null);
				sam.setAttribute("IT", null);

				out.addAlignment(sam);
				

			}
			it.close();
			out.close();
		}
		System.out.println();
		
	}

	public void printAlignmentStats(){
		System.out.println("\nAlignment statistics for building models:");
		System.out.println(numReads+"\tNumber alignments");

		double frac = (double) numFlaggedReads/(double)numReads;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numFlaggedReads+")\tAlignments failing vender QC, duplicate, unmapped mate, unmapped, proper pair");

		frac = (double) numCtrlChromosomes/(double)numReads;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numCtrlChromosomes+")\tAlignments to chrPhiX or chrAdapters");

		frac = (double) numFailedScore/(double)numReads;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numFailedScore+")\tAlignments failing either the alignment ("+maximumAlignmentScore+") or mapping ("+minimumMappingQuality+") scores");

		frac = (double) numFailedCigar/(double)numReads;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numFailedCigar+")\tAlignments failing to match the CIGAR ("+cigar+") field");

		frac = (double) numPassingReads/(double)numReads;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numPassingReads+")\tAlignments passing all filters and parsed for models.");
		System.out.println(numberModels +"\tNumber models matched for strand, position, and first/second sequence.");
		System.out.println();
	}

	public void printBaseClassifierStats(){
		System.out.println("\nBase classifier statistics:");
		System.out.println(numVarAlign+"\tNumber variant alignments");
		double frac = (double) numVarAlignFailingMinCor/(double)numVarAlign;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numVarAlignFailingMinCor+")\tNumber variant alignments failing minimum likelihood ("+minimumLikelihood+")");
		frac = (double)numVarAlignFailingMinCorDiff/(double)numVarAlign;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numVarAlignFailingMinCorDiff+")\tNumber variant alignments failing minimum likelihood difference ("+minimumLikelihoodDifference+")");
		frac = (double)numVarAlignChanged2N/(double)numVarAlign;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numVarAlignChanged2N+")\tNumber variant center base changed to N");
		frac = (double)numVarAlignChanged2Novel/(double)numVarAlign;
		System.out.println(Num.formatNumber(frac, 3)+" ("+numVarAlignChanged2Novel+")\tNumber variant center base changed to different base");
	}

	public void scoreVariants(){

		//for each variant
		for (int i=0; i< variants.length; i++){
			VariantAlignment[] alignments = variantAlignments[i];

			if (alignments == null) continue;

			//for each VariantAlignment
			if (useLikelihood) for (VariantAlignment va: alignments)  calculateLikelihoods(va);
			else for (VariantAlignment va: alignments)  correlate(va);

		}

	}

	public void printVariants(){
		
		System.out.println("Chr\tStart\tStop\tName\tScore\t#Align\t#A\t#C\t#G\t#T\t#N\t#cA\t#cC\t#cG\t#cT\t#cN\tChanged");
		//for each variant
		for (int i=0; i< variants.length; i++){
			Bed variant = variants[i];	
			
			//any alignments?
			VariantAlignment[] alignments = variantAlignments[i];
			if (alignments != null){
				int[] countsACGTN = new int[5];
				int[] cCountsACGTN = new int[5];
				//for each VariantAlignment
				boolean changedVariant = false;
				for (VariantAlignment va: alignments) {
					//increment base counters
					countBase(va.getCenterSeq(), countsACGTN);
					countBase(va.getClassifiedCenterSeq(), cCountsACGTN);
					if (va.isChanged() && changedVariant == false) changedVariant = true;
				}
				if ((printOnlyChanged == false) || (printOnlyChanged == true && changedVariant == true)){
					StringBuilder sb = new StringBuilder();
					//chr start stop name score strand
					sb.append(variant.toStringNoStrand());
					sb.append("\t");
					//num overlapping alignments
					sb.append(alignments.length);
					//Original calls
					for (int x=0; x<countsACGTN.length; x++){
						sb.append("\t");
						sb.append(countsACGTN[x]);
					}
					//re calls
					for (int x=0; x<cCountsACGTN.length; x++){
						sb.append("\t");
						sb.append(cCountsACGTN[x]);
					}
					//any changes
					sb.append("\t");
					sb.append(changedVariant);
					sb.append("\n");
					
					//print individual calls?
					if (changedVariant == true && printLikelihoods == true){
						
						for (VariantAlignment va: alignments) {
							sb.append("\t"+va.getIndexStrand()+"\t"+va.getVariantSeq()+"\t"+va.isChangedToN()+"\t"+va.isChangedToDiff()+"\t"+va.getAlignment().getReadName());
							sb.append("\n\t");
							for (String s: va.getVariantSeqs()) {
								sb.append(s);
								sb.append("\t");
							}
							sb.append("\n\t");
							for (double s: va.getLikelihoodACGTRatios()) {
								sb.append(s);
								sb.append("\t");
							}
							sb.append("\n");
						}
					}
					System.out.println(sb.toString());
					
				}

			}


		}
	}

	public void countBase(String base, int[] countsACGTN){
		if (base.equals("A")) countsACGTN[0]++;
		else if (base.equals("C")) countsACGTN[1]++;
		else if (base.equals("G")) countsACGTN[2]++;
		else if (base.equals("T")) countsACGTN[3]++;
		else countsACGTN[4]++;
	}

	public void correlate(VariantAlignment va){
		//get first or second matrix set
		HashMap<String, HashMap<String, NMer>> matri;
		if (va.getAlignment().getFirstOfPairFlag() == false) matri = secondReadMatrixLookup;
		else matri = firstReadMatrixLookup;

		double[] alignmentInts = va.fetchACGTIntensities();
		String[] varSeqs = va.getVariantSeqs();
		float[] corr = new float[varSeqs.length];

		//for each different ACTG seq at center base
		for (int i=0; i< varSeqs.length; i++){
			//get NMer
			NMer nmer = matri.get(varSeqs[i]).get(va.getIndexStrand());
			//fetch intensity arrays
			double[] matrixInts = nmer.fetchACTGMeans();
			double[] matrixS2N = nmer.fetchACTGSig2Noise();
			corr[i] = (float)fetchWeightedCorrelation(alignmentInts, matrixInts, matrixS2N);
		}		
		va.setWeightedCorrelationsACGT(corr);
		va.callCenterBaseWithCorrelation();

	}
	
	public void calculateLikelihoods(VariantAlignment va){
		//get first or second matrix set
		HashMap<String, HashMap<String, NMer>> matri;
		if (va.getAlignment().getFirstOfPairFlag() == false) matri = secondReadMatrixLookup;
		else matri = firstReadMatrixLookup;

		double[] alignmentInts = va.fetchACGTIntensities();
		String[] varSeqs = va.getVariantSeqs();

		float[][] meansACGT = new float[varSeqs.length][];
		float[][] stdsACGT = new float[varSeqs.length][];
		
		//for each different ACTG seq at center base, fetch intensities and stndev
		for (int i=0; i< varSeqs.length; i++){
			//get NMer
			NMer nmer = matri.get(varSeqs[i]).get(va.getIndexStrand());
			//fetch intensity arrays, note this will null the StandardDeviation objects in the NMer so other NMer method calls may return null pointer errors.
			meansACGT[i] = nmer.getMeansACGT();
			stdsACGT[i] = nmer.getStandardDeviationsACGT();
		}	
		
		Classifier c = new Classifier(meansACGT, stdsACGT);
		double[] l = c.calculateLikelihoods(alignmentInts);
		
		va.setLikelihoodACGTRatios(l);
		va.callCenterBaseWithLikelihood();

	}
	

	public double fetchWeightedCorrelation(double[] a, double [] b, double[] w){
		wpc.reset();
		for (int i=0; i< a.length; i++) wpc.put(a[i], b[i], w[i]);
		return wpc.getCorrelation();
	}

	public void printModels(){
		System.out.println("First read models");
		printNMers(firstReadMatrixLookup);

		System.out.println("\nSecond read models");
		printNMers(secondReadMatrixLookup);
	}

	public void printNMers(HashMap<String, HashMap<String, NMer>> matrix){
		//for each sequence
		for (String seq:matrix.keySet()){
			System.out.println(seq);
			//for each nmer
			HashMap<String, NMer> nmers = matrix.get(seq);
			for (String indexStrand: nmers.keySet()){
				NMer nmer = nmers.get(indexStrand);
				System.out.println("\t"+indexStrand+" "+nmer.getNumberObservations());

				System.out.println(nmer.getMeanMatrix());
				System.out.println(nmer.getCVMatrix());
			}

		}
	}

	public boolean checkAlignment(SAMRecord sam){

		//any bad flags?
		if (sam.getReadFailsVendorQualityCheckFlag() || sam.getDuplicateReadFlag() || sam.getMateUnmappedFlag() || sam.getReadUnmappedFlag()) {
			numFlaggedReads++;
			return false;
		}

		//proper pair?
		if (sam.getReadPairedFlag()){
			if (sam.getProperPairFlag() == false){
				numFlaggedReads++;
				return false;
			}
		}

		//skip phiX and adapter
		if (sam.getReferenceName().startsWith("chrPhiX") || sam.getReferenceName().startsWith("chrAdapt")){
			numCtrlChromosomes++;
			return false;
		}

		//does it pass the score thresholds?
		List<SAMTagAndValue> attributes = sam.getAttributes();
		int alignmentScore = Integer.MIN_VALUE;
		for (SAMTagAndValue tagVal : attributes){
			String tag = tagVal.tag;
			if (tag.equals("AS")){
				alignmentScore = (Integer)tagVal.value;
				break;
			}
		}
		if (alignmentScore != Integer.MIN_VALUE){
			if (alignmentScore > maximumAlignmentScore){
				numFailedScore++;
				return false;
			}
		}
		int mappingQuality = sam.getMappingQuality();
		if (mappingQuality < minimumMappingQuality){
			numFailedScore++;
			return false;
		}

		//check CIGAR
		if (sam.getCigarString().equals(cigar) == false){
			numFailedCigar++;
			return false;
		}

		return true;
	}

	public void loadModels(){
		int counter = 0;

		//for each bam reader
		for (SAMFileReader bamReader: bamReaders){
			SAMRecordIterator it = bamReader.iterator();

			while (it.hasNext()) {
				//load record and check
				SAMRecord sam = it.next();
				numReads++;
				if (counter++ > 1000000) {
					counter = 0;
					System.out.print(".");
				}
				if (checkAlignment(sam) == false) continue;
				numPassingReads++;

				//scan alignment
				scanQualityAlignment(sam);

			}
			it.close();
		}
		System.out.println();
	}

	public void scanQualityAlignment(SAMRecord sam){
		try {
			//first or second read matri?
			HashMap<String, HashMap<String, NMer>> matri;
			if (sam.getProperPairFlag() && sam.getFirstOfPairFlag() == false)matri = secondReadMatrixLookup;
			else matri = firstReadMatrixLookup;

			//strand
			String strand = "+";
			if (sam.getReadNegativeStrandFlag()) strand = "-";

			//fields
			String seq = sam.getReadString();
			byte[] qualities = null;
			short[] aIntensities = null;
			short[] cIntensities = null;
			short[] gIntensities = null;
			short[] tIntensities = null;

			int seqLength = seq.length();
			int startIndex;
			int endIndex;

			//for each index position
			for (int i=0; i< seqLength; i++){

				//define a substring
				startIndex = i - nMerAdder;
				if (startIndex < 0) startIndex = 0;
				endIndex = i + nMerAdder +1;
				if (endIndex > seqLength) endIndex = seqLength;
				String subSeq = seq.substring(startIndex, endIndex);

				//look for it in the hash
				HashMap<String, NMer> nmers = matri.get(subSeq);
				if (nmers != null){
					String indexStrand = i+strand;
					//look for the proper index and strand
					NMer nmer = nmers.get(indexStrand);
					if (nmer != null){
						//System.out.println("\t\tFOUND"+indexStrand+" "+subSeq);

						//check qualities
						//look for bad qualities
						if (minimumBaseQuality !=0){
							if (qualities == null) qualities = sam.getBaseQualities();
							if (qualitiesOK (startIndex, endIndex, qualities)){
								//load intensity arrays?
								if (aIntensities == null){
									//System.out.println(sam.getSAMString());
									aIntensities = sam.getSignedShortArrayAttribute("IA");
									cIntensities = sam.getSignedShortArrayAttribute("IC");
									gIntensities = sam.getSignedShortArrayAttribute("IG");
									tIntensities = sam.getSignedShortArrayAttribute("IT"); 
									//if negative strand must reverse the intensity arrays, these then represent the complement sequence 
									if (sam.getReadNegativeStrandFlag()){
										aIntensities = Num.reverse(aIntensities);
										cIntensities = Num.reverse(cIntensities);
										gIntensities = Num.reverse(gIntensities);
										tIntensities = Num.reverse(tIntensities); 
										//Misc.printArray(aIntensities);
										//Misc.printArray(cIntensities);
										//Misc.printArray(gIntensities);
										//Misc.printArray(tIntensities);
									}
								}
								//add intensities to nmer!
								//System.out.println("\t\t\t adding "+startIndex + "-" + endIndex);
								nmer.count(startIndex, endIndex, aIntensities, cIntensities, gIntensities, tIntensities);
							}
						}

					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR counting intensities for nmer.\n");
		}



	}

	public boolean qualitiesOK (int startIndex, int endIndex, byte[] qualities){
		int start = startIndex - numberFlankingBaseQualitiesToCheck;
		if (start < 0) start = 0;
		int end = endIndex + numberFlankingBaseQualitiesToCheck;
		if (end > qualities.length) end = qualities.length;
		for (int x=start; x< end; x++){
			if (qualities[x] < minimumBaseQuality) {
				return false;
			}
		}
		return true;
	}
	
	

	public void scanVariantAlignments(){
		//for each variant
		for (int i=0; i< variants.length; i++){

			VariantAlignment[] alignments = variantAlignments[i];

			//any alignments
			if (alignments == null) continue;

			//for each alignment
			for (VariantAlignment va: alignments){

				//get first or second matrix set
				HashMap<String, HashMap<String, NMer>> matri;
				if (va.getAlignment().getFirstOfPairFlag() == false) matri = secondReadMatrixLookup;
				else matri = firstReadMatrixLookup;

				//for each ACGT sequence, TT[ACGT]CC
				String[] varSeqs = va.getVariantSeqs();
				for (int x=0; x< varSeqs.length; x++){
					//fetch nmers indexStrand, 23+ with nmer object
					HashMap<String, NMer> nmers = matri.get(varSeqs[x]);
					//not seen before so make it new
					if (nmers == null){
						nmers = new HashMap<String, NMer>();
						matri.put(varSeqs[x], nmers);
						numberModels++;
					}
					//look to see if the IndexStand NMer exists
					if (nmers.containsKey(va.getIndexStrand()) == false){
						nmers.put(va.getIndexStrand(), new NMer((byte)varSeqs[x].length()));
					}


				}
			}
		}
	}
	
	/*
	public void scoreForStrandBias(){
		System.out.println("Scoring variants for strand and position bias...");
		double numVars = variants.length;
		
		//for each variant
		for (int i=0; i< variants.length; i++){
			Bed variant = variants[i];			
			VariantAlignment[] vas = variantAlignments[i];
			
			//strand bias?
			int numPlusStrand = 0;
			int numMinusStrand = 0;
			int maxReadLengthInt = 0;
			for (VariantAlignment va : vas){
				SAMRecord sam = va.getAlignment();
				if (sam.getReadNegativeStrandFlag()) numMinusStrand++;
				else numPlusStrand++;
				if (maxReadLengthInt < sam.getReadLength()) maxReadLengthInt = sam.getReadLength();
			}
			double maxReadLength = maxReadLengthInt;
			BinomialDistribution bd = new BinomialDistribution(numPlusStrand+numMinusStrand, 0.5);
			double pvalStrand = bd.getProbabilityOfSuccess() * numVars;
			
			//position bias? 
			int numBasesPerBlock = maxReadLength/5.0;
			double expect = ((double)(numPlusStrand+numMinusStrand)) / 5 ;
			long[] hits = new long[5];
			double[] expects = new double[5];
			Arrays.fill(expects, expect);
			for (VariantAlignment va : vas){
				int centerIndex = va.getCenterIndex();
				//scale it?
				int readLength = va.getAlignment().getReadLength();
				if (readLength != maxReadLength) {
					float scalar = (float)maxReadLength/ (float)readLength;
					centerIndex = Math.round(scalar * (float)centerIndex);
				}
				hits[centerIndex]++;
			}
			//probably not a good test since there are alot of cells with zero, better to bin by quartile and score
			ChiSquareTest chi = new ChiSquareTest();
			double pvalPosition = chi.chiSquareTest(expects, hits) * numVars;
			
			//adjacent variants?
			boolean adjacentFound = false;
			//look left
			if (i!=0){
				Bed pre = variants[i-1];
				adjacentFound = checkForAdjacentVariant(pre, variant);
			}
			//look right
			if (adjacentFound == false){
				int index = i+1;
				if (index != variants.length) {
					Bed post = variants[index];
					adjacentFound = checkForAdjacentVariant(variant, post);
				}
			}
			
		}
	}*/
	
	public boolean checkForAdjacentVariant(Bed left, Bed right){
		if (left.getChromosome().equals(right.getChromosome())){
			if (left.intersects(right.getStart()- maxAdjacentDistance, right.getStop()+ maxAdjacentDistance)) return true;
		}
		return false;
	}

	public void loadVariantAlignments(){
		System.out.println("Loading variants...");

		//for each variant
		for (int i=0; i< variants.length; i++){
			Bed variant = variants[i];			

			//fetch alignments, might be null if none found
			ArrayList<SAMRecord> al = fetchAllAlignments(variant);

			//create VariantAlignments
			if (al != null) {
				
				//position to fetch,this works for snps but how about INDELs?
				int centerBasePosition = variant.getMiddle();
				
				//for each alignment
				ArrayList<VariantAlignment> alVA = new ArrayList<VariantAlignment>();
				int num = al.size();
				for (int x=0; x< num; x++){
					SAMRecord sam = al.get(x);
					VariantAlignment va = new VariantAlignment(this, sam, centerBasePosition, nMerAdder);
					//watch out for deletion and N alignments
					if (va.getIndexStrand() != null) alVA.add(va);
				}
				if (alVA.size() !=0){
					variantAlignments[i] = new VariantAlignment[alVA.size()];
					alVA.toArray(variantAlignments[i]);
					numVarAlign += variantAlignments[i].length;
				}
			}
		}
	}

	/**Checks that the alignment actually touches down on at least one base of the region to avoid spanners.*/
	public ArrayList<SAMRecord> fetchAlignments (Bed ei, SAMFileReader reader){
		ArrayList<SAMRecord> al = new ArrayList<SAMRecord>();
		SAMRecordIterator i = reader.queryOverlapping(ei.getChromosome(), (ei.getStart()+1), ei.getStop());
		while (i.hasNext()) {
			SAMRecord sam = i.next();
			//fetch blocks of actual alignment
			ArrayList<int[]> blocks = SamAlignmentExtractor.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
			//check to see if any actual touch down and intersect the region
			for (int[] b : blocks){
				if (ei.intersects(b[0], b[1])){
					al.add(sam);
					break;
				}
			}
		}
		i.close();
		i = null;
		if (al.size() == 0) return null;
		return al;
	}

	public ArrayList<SAMRecord> fetchAllAlignments(Bed variant){
		ArrayList<SAMRecord> al = new ArrayList<SAMRecord>();
		//for each bam reader
		for (SAMFileReader bamReader: bamReaders){
			ArrayList<SAMRecord> sub = fetchAlignments(variant, bamReader);
			if (sub != null) al.addAll(sub);
		}
		if (al.size() == 0) return null;
		return al;
	}

	public void incrementNumVarAlignFailingMinCor(){
		numVarAlignFailingMinCor++;
	}
	public void incrementNumVarAlignFailingMinCorDiff(){
		numVarAlignFailingMinCorDiff++;
	}
	public void incrementNumVarAlignChanged2N(){
		numVarAlignChanged2N++;
	}
	public void incrementNumVarAlignChanged2Novel(){
		numVarAlignChanged2Novel++;
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BaseClassifier(args);
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
					case 'a': bamFiles = IO.extractFiles(new File(args[++i]), ".bam"); break;
					case 'v': bedFile = new File(args[++i]); break;
					case 'n': setRecalledBasesToN = true; break;
					case 'l': minimumLikelihood = Double.parseDouble(args[++i]); break;
					case 'd': minimumLikelihoodDifference = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check bam files
		if (bamFiles == null || bamFiles.length ==0 || bamFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam file(s)!\n");
		OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, false);

		//check variant file
		if (bedFile == null || bedFile.canRead() == false) Misc.printExit("\nError: cannot find or read your bed formated variant file!\n");


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Base Classifier : Oct 2012                            **\n" +
				"**************************************************************************************\n" +
				"Beta.\n\n" +

				"Options:\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/BamIntensityParser -f /Data/BamFiles/\n" +
				"       -n 7 -q 30\n\n"+

		"**************************************************************************************\n");

	}

	public float getMinimumCorrelation() {
		return minimumCorrelation;
	}

	public float getMinimumCorrelationDifference() {
		return minimumCorrelationDifference;
	}

	public HashMap<String, SAMRecord> getChangedAlignments() {
		return changedAlignments;
	}

	public double getMinimumLikelihood() {
		return minimumLikelihood;
	}

	public double getMinimumLikelihoodDifference() {
		return minimumLikelihoodDifference;
	}

	public boolean isSetRecalledBasesToN() {
		return setRecalledBasesToN;
	}


}
