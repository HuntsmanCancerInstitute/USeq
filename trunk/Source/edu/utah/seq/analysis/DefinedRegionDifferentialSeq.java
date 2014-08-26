package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import org.apache.poi.ss.usermodel.*;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import htsjdk.samtools.*;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.analysis.multi.Condition;
import edu.utah.seq.analysis.multi.GeneCount;
import edu.utah.seq.analysis.multi.GeneResult;
import edu.utah.seq.analysis.multi.PairedCondition;
import edu.utah.seq.analysis.multi.Replica;
import edu.utah.seq.data.HeatMapMakerPosNeg;
import edu.utah.seq.data.Info;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.apps.Bar2USeq;


/** Compares datasets for differential counts over user defined regions (e.g. gene models).  
 * Wraps Simon Anders/ Michael Love's et al. DESeq2 package with lots of modifications. Estimates alternative splicing.
 * @author Nix
 * */
public class DefinedRegionDifferentialSeq {

	//user defined fields
	private File[] conditionDirectories;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private File bedFile;
	private File refSeqFile;
	private String genomeVersion;
	private float minAdjP = 20;
	private float minLog2Ratio = 1f;
	private boolean scoreIntrons = false;
	private boolean removeOverlappingRegions = false;
	private int minimumCounts = 10;
	private boolean deleteTempFiles = true;
	private int minimumSpliceCounts = 10;
	private boolean performStrandedAnalysis = false;
	private boolean performReverseStrandedAnalysis = false;
	private int maxNumAlignments = 1;
	private boolean verbose = true;
	private int maxAlignmentsDepth = 50000;
	private boolean secondStrandFlipped = false;
	private boolean useSamSeq = false;
	
	//hidden options
	private boolean trimUTRBPsFromExons = false;
	private boolean printFirstLastCountTable = false;
	private int maxFirstLast = 150;

	//internal fields
	private UCSCGeneLine[] genes;
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private HashMap<String, UCSCGeneLine> name2Gene;
	private Condition[] conditions;
	private HashSet<String> flaggedGeneNames = new HashSet<String>();
	private HashSet<String> geneNamesPassingThresholds = new HashSet<String>();
	private LinkedHashSet<String> geneNamesWithMinimumCounts = new LinkedHashSet<String>();
	private String[] geneNamesToAnalyze = null;
	private File serializedConditons = null;
	private String url;
	private static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MSDHN])");
	public static final Pattern BAD_NAME = Pattern.compile("(.+)/[12]$");
	public boolean saveCounts = true;
	public File rLogValues;

	//for loading data
	private int workingFragmentNameIndexPlus = 1;
	private int workingFragmentNameIndexMinus = -1;
	private HashMap<String, Integer> workingFragNameIndex = new HashMap<String, Integer>(10000);

	//for paired diff expression
	private File geneCountTable;
	private String[] conditionNamesPerReplica;
	private String[] replicaNames;
	private PairedCondition[] pairedConditions;
	private File spliceGraphDirectory;
	private File workingSpliceGraphDir;

	//from RNASeq app integration
	private File treatmentBamDirectory;
	private File controlBamDirectory;

	//spreadsheet
	private Workbook workbook = null;

	//constructors
	/**Stand alone.*/
	public DefinedRegionDifferentialSeq(String[] args){	
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);
		//launch
		run();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	/**For integration with RNASeq app.*/
	public DefinedRegionDifferentialSeq(File treatmentBamDirectory, File controlBamDirectory, String genomeVersion, File saveDirectory,  File fullPathToR, File processedRefSeqFile, boolean scoreIntrons, boolean performStrandedAnalysis, 
			boolean verbose, boolean reverse, boolean flipped, int maxAlignDepth, boolean useSamSeq){
		this.treatmentBamDirectory = treatmentBamDirectory;
		this.controlBamDirectory = controlBamDirectory;
		this.saveDirectory = saveDirectory;
		this.genomeVersion = genomeVersion;
		url = "http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
		this.fullPathToR = fullPathToR;
		this.refSeqFile = processedRefSeqFile;
		this.scoreIntrons = scoreIntrons;
		this.performStrandedAnalysis = performStrandedAnalysis;
		this.verbose = verbose;
		removeOverlappingRegions = false;
		this.secondStrandFlipped = flipped;
		this.performReverseStrandedAnalysis = reverse;
		this.maxAlignmentsDepth = maxAlignDepth;
		this.useSamSeq = useSamSeq;
		run();
	}

	public void run(){
		//run in non X11 mode
		System.setProperty("java.awt.headless", "true");

		//load gene models
		if (verbose) System.out.println("Loading regions/ gene models...");
		loadGeneModels();

		//load count data by replica 
		loadConditions();

		//write out count table of genes passing minimum counts
		writeCountTable();

		//launch DESeq2/SAMSeq for differential expression
		String thresholdName = "-10Log10(AdjP) "+Num.formatNumber(minAdjP,1)+", Log2Ratio "+Num.formatNumber(minLog2Ratio, 1);

		if (useSamSeq) {
			if (verbose)  System.out.println("Running SAMseq analysis...");
			analyzeUsingSamSeq();
		} else {
			if (verbose) System.out.println("Running DESeq2 analysis...");
			runDESeq2();
		}

		if (verbose) System.out.println("\n\t"+geneNamesPassingThresholds.size()+" / "+geneNamesToAnalyze.length+"\tgenes differentially expressed ("+thresholdName+")\n");

		//launch Diff splice
		if (verbose) System.out.println("Running chi-square tests for differential splicing...");
		spliceGraphDirectory = new File (saveDirectory, "DiffSpliceLog2RtoGraphs");
		spliceGraphDirectory.mkdirs();
		analyzeForDifferentialSplicing();
		new Bar2USeq(spliceGraphDirectory,true);

		//set max deseqAdjP and Log2Ratio in geneNamesToAnalyze
		setMaxScores();

		//print spreadsheet 
		printStatSpreadSheet();


	}

	private void setMaxScores() {
		float[] maxAdjPs = new float[geneNamesToAnalyze.length];
		float[] maxLog2 = new float[geneNamesToAnalyze.length];
		//for each pair
		for (int i=0; i< pairedConditions.length; i++){
			//gene : float[]{AdjP, log2Rto}
			float[][] scores = pairedConditions[i].getParsedDiffExpResults();
			for (int j=0; j< geneNamesToAnalyze.length; j++){
				if (scores[j][0] > maxAdjPs[j]) maxAdjPs[j] = scores[j][0];
				float abs = Math.abs(scores[j][1]);
				if (abs > maxLog2[j]) maxLog2[j] = abs;
			}
		}
		//set in genes
		for (int i=0; i< geneNamesToAnalyze.length; i++){
			UCSCGeneLine gene = name2Gene.get(geneNamesToAnalyze[i]);
			gene.setMaxAbsLog2Ratio(maxLog2[i]);
			gene.setFdr(maxAdjPs[i]);  //note misnaming!
		}
	}
	
	private void analyzeForDifferentialSplicing() {
		//for each PairedCondition
		for (PairedCondition pc : pairedConditions){
			//make dir to hold splice graph
			workingSpliceGraphDir = new File (spliceGraphDirectory, pc.getName());
			workingSpliceGraphDir.mkdirs();
			workingSpliceGraphDir.deleteOnExit();
			estimateDifferencesInReadDistributions(pc);
		}

	}

	private void loadConditions(){

		//Any saved conditions
		serializedConditons = new File (saveDirectory, "conditions.ser");
		if (serializedConditons.exists()){
			conditions = (Condition[])IO.fetchObject(serializedConditons);
			if (verbose) System.out.println("\nWARNING: Loading "+conditions.length+" conditions from cached data file, delete "+serializedConditons+" if you'd like to recount gene exon/ regions....\n");
			//load min count hash, this is normally done in the scanGene() method
			loadMinimumCountsHash();
		}

		else {
			String strand = "";
			if (performStrandedAnalysis) strand = " stranded ";
			System.out.println("\nCollecting "+ strand +"counts for each gene exon/ region...");
			//from RNASeq app?
			if (treatmentBamDirectory != null){
				conditions = new Condition[2];
				conditions[0] = new Condition(treatmentBamDirectory);
				conditions[1] = new Condition(controlBamDirectory);
			}
			else {
				conditions = new Condition[conditionDirectories.length];				
				for (int i=0; i< conditionDirectories.length; i++) conditions[i] = new Condition(conditionDirectories[i]);
			}

			//If SAMseq, make sure there are least two replicas
			if (useSamSeq) {
				for (Condition c: conditions) {
					if (c.getReplicas().length < 2) {
						System.out.println("You must have at least two replicates per condition in order to run SAMseq");
						System.exit(1);
					}
				}
			}

			//load em with data
			for (int i=0; i< conditions.length; i++) {				
				for (Replica r: conditions[i].getReplicas()){			
					if (verbose) System.out.print("\t"+r.getNameNumber());
					loadReplica(r);
				}
			}
			//any genes with too many reads that were excluded?
			if (flaggedGeneNames.size() !=0) {
				String[] badGeneNames = Misc.hashSetToStringArray(flaggedGeneNames);			

				//remove flagged genes from all replicas
				//for each condition
				for (Condition c: conditions) {
					//for each replica
					for (Replica replica: c.getReplicas()) replica.removeFlaggedGenes(badGeneNames);
				}
				//remove them from geneNamesWithMinimumCounts
				for (String baddie: badGeneNames){
					geneNamesWithMinimumCounts.remove(baddie);
				}
			}
			if (verbose) {
				System.out.println();
				System.out.println(geneNamesWithMinimumCounts.size()+" genes, with >= "+minimumCounts+" and < "+maxAlignmentsDepth+" counts, will be examined for differential expression.\n ");
			}
			
			//save conditions?
			if (saveCounts) IO.saveObject(serializedConditons, conditions);

		}
		geneNamesToAnalyze = Misc.hashSetToStringArray(geneNamesWithMinimumCounts);

		//print conditions
		if (verbose) {
			System.out.println("Conditions, replicas, and total gene region mapped counts:");
			for (Condition c: conditions) System.out.println(c);
		}
	}

	public void loadBlocks(SAMRecord sam, int chrStartBp, ArrayList<Integer>[] bpNames, boolean[] badBases){

		NameInteger nameIndex = null;
		boolean addIt;
		ArrayList<int[]> blocks = fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
		//add name to each bp
		for (int[] b : blocks){
			int start = b[0] - chrStartBp;
			int stop = b[1] - chrStartBp;
			for (int i=start; i < stop; i++){
				//bad base?
				if (badBases[i]) continue;
				addIt = true;

				//never seen before?
				if (bpNames[i] == null) bpNames[i] = new ArrayList<Integer>();
				//old so check size
				else if (bpNames[i].size() == maxAlignmentsDepth) {
					badBases[i] = true;
					addIt = false;
					bpNames[i] = null;	
					if (nameIndex !=null) workingFragNameIndex.remove(nameIndex.name);
				}

				// add it?
				if (addIt) {
					if (nameIndex == null) nameIndex = fetchFragmentNameIndex(sam);
					bpNames[i].add(nameIndex.index);
				}
			}
		}
	}

	/**Assumes interbase coordinates for start and returned blocks.*/
	public static ArrayList<int[]> fetchAlignmentBlocks(String cigar, int start){
		//for each cigar block
		Matcher mat = CIGAR_SUB.matcher(cigar);
		ArrayList<int[]> blocks = new ArrayList<int[]>();
		while (mat.find()){
			String call = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (call.equals("M")) {
				blocks.add(new int[]{start, start+numberBases});
			}
			//just advance for all but insertions which should be skipped via the failure to match
			start += numberBases;
		}
		return blocks;
	}

	/**Fetches an old or makes a new Integer to represent the sam read name (e.g. fragment name)*/
	public NameInteger fetchFragmentNameIndex(SAMRecord sam){
		String samReadName = sam.getReadName();
		Integer index;
		Matcher mat = BAD_NAME.matcher(samReadName);
		if (mat.matches()) samReadName = mat.group(1);

		if (workingFragNameIndex.containsKey(samReadName)) {
			index = workingFragNameIndex.get(samReadName);
		}

		else {
			if (!secondStrandFlipped && sam.getReadPairedFlag()) {
				if ((sam.getFirstOfPairFlag() && sam.getReadNegativeStrandFlag()) || (sam.getSecondOfPairFlag() && !(sam.getReadNegativeStrandFlag()))) {
					index = new Integer(workingFragmentNameIndexMinus--);
				}
				else {
					index = new Integer(workingFragmentNameIndexPlus++);
				}
			} else {
				if (sam.getReadNegativeStrandFlag()) {
					index = new Integer(workingFragmentNameIndexMinus--);
				} else {
					index = new Integer(workingFragmentNameIndexPlus++);
				}
			}

			workingFragNameIndex.put(samReadName, index);
		}
		return new NameInteger (samReadName, index);
	}

	private class NameInteger {
		String name;
		Integer index;

		public NameInteger(String name, Integer index){
			this.name = name;
			this.index = index;
		}
	}

	public void loadReplica(Replica replica){
		//the idea here is to take a sorted bam file and add the reads to an ArrayList, basically every base should have a ArrayList of reads that overlap that base
		//one can then count the number of overlapping fragments for an exon or a collection of exons by hashing the ArrayList
		//assumes each read (first or second) has the same name

		//make reader
		SAMFileReader reader = new SAMFileReader(replica.getBamFile());	
		reader.setValidationStringency(ValidationStringency.SILENT);


		//fetch chromName: length for all chroms
		HashMap<String, Integer> chromLength = new HashMap<String, Integer>();
		List<SAMSequenceRecord> seqs =  reader.getFileHeader().getSequenceDictionary().getSequences();
		for (SAMSequenceRecord sr: seqs) chromLength.put(sr.getSequenceName(), sr.getSequenceLength());

		SAMRecordIterator iterator = reader.iterator();

		HashSet<String> priorChroms = new HashSet<String>();
		String chrom = null;

		//make first chromosome
		SAMRecord sam = null;
		int chrStartBp = -1;
		ArrayList<Integer>[] bpNames = null;
		boolean[] badBases = null;
		while (iterator.hasNext()){
			sam = iterator.next();

			//unaligned? too many hits?
			if (alignmentFails(sam)) continue;
			chrom = sam.getReferenceName();

			if (chromGenes.containsKey(chrom) == false) {
				continue;
			}
			if (verbose) System.out.print(".");			
			priorChroms.add(chrom);
			chrStartBp = sam.getUnclippedStart()-1;
			badBases = new boolean[chromLength.get(chrom) - chrStartBp];
			bpNames = new ArrayList[badBases.length];
			break;
		}
		//any data found?
		if (bpNames == null) {
			reader.close();
			return;
		} 

		//reset working fields
		workingFragNameIndex.clear();
		workingFragmentNameIndexPlus = 1;
		workingFragmentNameIndexMinus = -1;

		//load first alignment into bpNames
		loadBlocks (sam, chrStartBp, bpNames, badBases);

		//for each record
		while (iterator.hasNext()){
			sam = iterator.next();

			//unaligned? too many hits?
			if (alignmentFails(sam)) continue;

			//same chrom?
			if (sam.getReferenceName().equals(chrom)){
				loadBlocks (sam, chrStartBp, bpNames, badBases);
			}
			else {
				//different chrom so time to scan
				if (printFirstLastCountTable) loadFirstLastGeneCounts(replica, bpNames, chrStartBp, chrom);
				else loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);

				//check that new chrom from SAM is something interrogated by their gene list
				boolean reset = false;
				if (chromGenes.containsKey(sam.getReferenceName()) == false) {
					//advance until it does
					while (iterator.hasNext()){
						sam = iterator.next();
						if (chromGenes.containsKey(sam.getReferenceName())) {
							reset = true;
							break;
						}
					}
				}
				else reset = true;

				//reset
				if (reset){
					//reset working fields
					workingFragNameIndex.clear();
					workingFragmentNameIndexPlus = 1;
					workingFragmentNameIndexMinus = -1;

					chrom = sam.getReferenceName();
					if (verbose) System.out.print(".");
					if (priorChroms.contains(chrom)) Misc.printErrAndExit("\nError: your sam file isn't sorted by chromosome! Aborting.\n");
					priorChroms.add(chrom);
					chrStartBp = sam.getUnclippedStart()-1;
					badBases = new boolean[chromLength.get(chrom) - chrStartBp];
					bpNames = new ArrayList[badBases.length];

					//load
					loadBlocks (sam, chrStartBp, bpNames, badBases);
				}

			}
		}
		//add last
		if (verbose) System.out.println();
		if (chromGenes.containsKey(sam.getReferenceName())) {
			if (printFirstLastCountTable) loadFirstLastGeneCounts(replica, bpNames, chrStartBp, chrom);
			else loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);
		}

		reader.close();
		bpNames = null;
		iterator = null;
		seqs = null;
		chromLength = null;
	}

	private boolean alignmentFails(SAMRecord sam){
		//aligned?
		if (sam.getReadUnmappedFlag()) return true;
		//limit to max matches?
		if (maxNumAlignments !=0){
			Object o = sam.getAttribute("NH");
			if (o != null)  {
				int num = (Integer)o;
				if (num > maxNumAlignments) return true;
			}
		}

		return false;
	}

	private void loadGeneCounts(Replica replica, ArrayList<Integer>[] bpNames, int chrStartBp, String chromosome, boolean[] badBases){
		HashMap<String, GeneCount> geneCounts = replica.getGeneCounts();
		HashSet<Integer> allReads = new HashSet<Integer>();
		HashSet<Integer> exonReads = new HashSet<Integer>();
		int lengthBpNames = bpNames.length -1;

		//for each gene in the chromosome 
		UCSCGeneLine[] chrGenes = chromGenes.get(chromosome);

		for (int i=0; i< chrGenes.length; i++){
			//flagged gene
			if (chrGenes[i].isFlagged()) continue;

			String geneName = chrGenes[i].getDisplayNameThenName();
			boolean plusStrand = chrGenes[i].getStrand().equals("+");					
			allReads.clear();

			//get exons
			ExonIntron[] exons = chrGenes[i].getExons();
			int[] exonCounts = new int[exons.length];

			//for each exon
			exonLoop:
				for (int x=0; x< exons.length; x++){
					exonReads.clear();
					int start = exons[x].getStart() - chrStartBp;
					if (start < 0) start = 0;

					int end = exons[x].getEnd() - chrStartBp;
					if (end > lengthBpNames) end = lengthBpNames;
					//for each base in the exon, see if there is a read
					for (int y=start; y< end; y++){
						//bad base? if so then flag entire gene
						if (badBases[y]) {
							chrGenes[i].setFlagged(true);
							flaggedGeneNames.add(geneName);
							allReads.clear();
							break exonLoop;
						}

						if (bpNames[y] != null) {
							if (performStrandedAnalysis){
								for (Integer id : bpNames[y]){
									if ((id > 0 && plusStrand == true) || (id < 0 && plusStrand == false )) exonReads.add(id);
								}
							} else if (performReverseStrandedAnalysis) {
								for (Integer id: bpNames[y]) {
									if ((id > 0 && plusStrand == false) || (id < 0 && plusStrand == true)) exonReads.add(id);
								}
							}
							else exonReads.addAll(bpNames[y]);
						}
					}

					exonCounts[x] = exonReads.size();
					allReads.addAll(exonReads);
				}

			int numCounts = allReads.size();

			if (numCounts !=0){
				GeneCount tcg = new GeneCount(numCounts, exonCounts);
				geneCounts.put(geneName, tcg);
				replica.setTotalCounts(replica.getTotalCounts() + numCounts);
				if (numCounts >= minimumCounts) geneNamesWithMinimumCounts.add(geneName);
			}
		}	

		//clean up
		allReads = null;
		exonReads = null;
	}

	private void loadFirstLastGeneCounts(Replica replica, ArrayList<Integer>[] bpNames, int chrStartBp, String chromosome){
		HashMap<String, GeneCount> geneCounts = replica.getGeneCounts();
		HashSet<Integer> allReads = new HashSet<Integer>();
		
		int lengthBpNames = bpNames.length -1;

		//for each gene in the chromosome 
		UCSCGeneLine[] chrGenes = chromGenes.get(chromosome);

		for (int i=0; i< chrGenes.length; i++){
			//flagged gene
			if (chrGenes[i].isFlagged()) continue;

			String geneName = chrGenes[i].getDisplayNameThenName();
			boolean plusStrand = chrGenes[i].getStrand().equals("+");					

			//get exons
			ExonIntron[] geneExons = chrGenes[i].getExons();
			
			//get exon set describing first third and last third of gene
			ExonIntron[] firstThird = ExonIntron.fetchFirstExonSet(geneExons, maxFirstLast, 3);
			ExonIntron[] lastThird = ExonIntron.fetchLastExonSet(geneExons, maxFirstLast, 3);
			
			//for each first third exons
			for (int x=0; x< firstThird.length; x++){
				//set first and last
				int start = firstThird[x].getStart() - chrStartBp;
				if (start < 0) start = 0;
				int end = firstThird[x].getEnd() - chrStartBp;
				if (end > lengthBpNames) end = lengthBpNames;
				//for each base in the exon, see if there is a read
				for (int y=start; y< end; y++){
					//any reads?
					if (bpNames[y] != null) {
						if (performStrandedAnalysis){
							for (Integer id : bpNames[y]){
								if ((id > 0 && plusStrand == true) || (id < 0 && plusStrand == false )) allReads.add(id);
							}
						} else if (performReverseStrandedAnalysis) {
							for (Integer id: bpNames[y]) {
								if ((id > 0 && plusStrand == false) || (id < 0 && plusStrand == true)) allReads.add(id);
							}
						}
						else allReads.addAll(bpNames[y]);
					}
				}
			}
			int numCountsFirst = allReads.size();
			allReads.clear();
			
			//for each last third exons
			for (int x=0; x< lastThird.length; x++){
				//set first and last
				int start = lastThird[x].getStart() - chrStartBp;
				if (start < 0) start = 0;
				int end = lastThird[x].getEnd() - chrStartBp;
				if (end > lengthBpNames) end = lengthBpNames;
				//for each base in the exon, see if there is a read
				for (int y=start; y< end; y++){
					//any reads?
					if (bpNames[y] != null) {
						if (performStrandedAnalysis){
							for (Integer id : bpNames[y]){
								if ((id > 0 && plusStrand == true) || (id < 0 && plusStrand == false )) allReads.add(id);
							}
						} else if (performReverseStrandedAnalysis) {
							for (Integer id: bpNames[y]) {
								if ((id > 0 && plusStrand == false) || (id < 0 && plusStrand == true)) allReads.add(id);
							}
						}
						else allReads.addAll(bpNames[y]);
					}
				}
			}
			int numCountsSecond = allReads.size();
			allReads.clear();
			
			int totalCounts = numCountsFirst + numCountsSecond;
			//order counts by gene strand
			int[] counts;
			if (chrGenes[i].getStrand().equals("-")) counts = new int[]{numCountsSecond, numCountsFirst};
			else counts = new int[]{numCountsFirst, numCountsSecond};
			
			if (totalCounts !=0){
				GeneCount tcg = new GeneCount(totalCounts, counts);
				geneCounts.put(geneName, tcg);
				replica.setTotalCounts(replica.getTotalCounts() + totalCounts);
				if (totalCounts >= minimumCounts) geneNamesWithMinimumCounts.add(geneName);
			}
		}	
		//clean up
		allReads = null;
	}

	public void loadGeneLineWithExonCounts(UCSCGeneLine gl, Condition treatment, Condition control){
		int numTReps = treatment.getReplicas().length;
		int numCReps = control.getReplicas().length;
		int numExons = gl.getExons().length;
		float[][] tCounts = new float[numTReps][numExons];
		float[][] cCounts = new float[numCReps][numExons];
		float totalTCounts = 0;
		float totalCCounts = 0;

		//for each rep
		Replica[] r = treatment.getReplicas();
		for (int i=0; i< numTReps; i++){
			GeneCount gc = r[i].getGeneCounts().get(gl.getDisplayNameThenName());
			if (gc != null) {
				tCounts[i] = Num.intArrayToFloat(gc.getExonCounts());
				totalTCounts+= gc.getCount();
			}
			else {
				tCounts[i] = new float[numExons];
				Arrays.fill(tCounts[i], 0);
			}
		}

		//for each rep
		r = control.getReplicas();
		for (int i=0; i< numCReps; i++){
			GeneCount gc = r[i].getGeneCounts().get(gl.getDisplayNameThenName());
			if (gc != null) {
				cCounts[i] = Num.intArrayToFloat(gc.getExonCounts());
				totalCCounts+= gc.getCount();
			}
			else {
				cCounts[i] = new float[numExons];
				Arrays.fill(cCounts[i], 0);
			}
		}
		gl.setTreatmentExonCounts(tCounts);
		gl.setControlExonCounts(cCounts);
		gl.setScores(new float[]{totalTCounts, totalCCounts});
	}

	public void estimateDifferencesInReadDistributions(PairedCondition pair){
		int maxNumberExons = -1;
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();

		//for each gene with min counts
		for (String geneName: geneNamesToAnalyze){
			UCSCGeneLine gl = name2Gene.get(geneName);
			//any exons?
			if (gl.getExons().length < 1) continue;

			//load it with counts
			loadGeneLineWithExonCounts(gl, pair.getFirstCondition(), pair.getSecondCondition());

			//null splice info
			gl.zeroNullSpliceScores();

			//check exon counts and modify if too few, sets a bunch of splice info too 
			if (checkForMinimums(gl)){
				al.add(gl);
				int numEx = gl.getExonCounts()[0].length;
				if (numEx > maxNumberExons) maxNumberExons = numEx;
			}
		}
		//any genes?
		if (al.size() == 0) return;

		UCSCGeneLine[] genesWithExonsAndReads = new UCSCGeneLine[al.size()];
		al.toArray(genesWithExonsAndReads);

		//collect counts
		int[][] treatment = new int[genesWithExonsAndReads.length][maxNumberExons];
		int[][] control = new int[genesWithExonsAndReads.length][maxNumberExons];

		for (int i=0; i< genesWithExonsAndReads.length; i++){
			float[][] tc = genesWithExonsAndReads[i].getExonCounts();
			Arrays.fill(treatment[i], -1);
			Arrays.fill(control[i], -1);
			int[] t = Num.convertToInt(tc[0]);
			int[] c = Num.convertToInt(tc[1]);
			System.arraycopy(t, 0, treatment[i], 0, t.length);
			System.arraycopy(c, 0, control[i], 0, c.length);
		}

		//estimate chi-square pvalues using R for resolution of extreemly small p-values, radiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, saveDirectory, fullPathToR, true);

		//bonferroni correction
		float bc = (float)Num.minus10log10(genesWithExonsAndReads.length);

		//add back
		for (int i=0; i< genesWithExonsAndReads.length; i++){
			//set corrected p-value 
			float pAdj = (float) pVals[i] + bc;
			if (pAdj > 0) genesWithExonsAndReads[i].setSplicingPValue(pAdj);
		}

		//collect all scored exons for diff splicing
		createSpliceGraphs(genesWithExonsAndReads);

		//extract info and put back in PairedCondition
		HashMap<String, float[]> geneNameSpliceStats = new HashMap<String, float[]>();
		for (UCSCGeneLine gene: genesWithExonsAndReads){
			ExonIntron maxLog2RtoSplicedExonIntron = gene.getMaxLog2RtoSplicedExon();
			float[] scores = new float[]{gene.getSplicingPValue(), gene.getSplicingLog2Ratio(), maxLog2RtoSplicedExonIntron.getStart(), maxLog2RtoSplicedExonIntron.getEnd()};
			geneNameSpliceStats.put(gene.getDisplayNameThenName(), scores);
		}
		pair.setGeneNameSpliceStats(geneNameSpliceStats);

	}

	/**Creates a square wave graph of exon log2Rtos for diff splicing.*/
	private void createSpliceGraphs(UCSCGeneLine[] genesWithExonsAndReads) {
		//load chromName : ExonIntron hashmap
		HashMap<String, ArrayList<ExonIntron>> chromExons = new HashMap<String, ArrayList<ExonIntron>>();
		String oldChr = "";
		ArrayList<ExonIntron> al = null;
		//walk through each gene
		for (int i=0; i< genesWithExonsAndReads.length; i++){
			//check thresholds, already thresholded at minLog2, check pval
			if (genesWithExonsAndReads[i].getSplicingPValue() < minAdjP) continue;
			//make ArrayList?
			if (genesWithExonsAndReads[i].getChrom().equals(oldChr) == false){
				oldChr = genesWithExonsAndReads[i].getChrom();
				al = chromExons.get(oldChr);
				if (al == null) {
					al = new ArrayList<ExonIntron>();
					chromExons.put(oldChr, al);
				}
			}
			//add Exons 
			ExonIntron[] exons = genesWithExonsAndReads[i].getScoredExons();
			for (ExonIntron e: exons) al.add(e);
		}
		//for each chromosome of exons 
		for (String chrom : chromExons.keySet()){
			al = chromExons.get(chrom);
			ExonIntron[] exons = new ExonIntron[al.size()];
			al.toArray(exons);
			Arrays.sort(exons);
			saveSmoothedHeatMapData(exons, chrom);
		}

	}


	/**Saves bar heatmap/ stairstep graph files*/
	public void saveSmoothedHeatMapData (ExonIntron[] exons, String chrom){
		HashMap<String,String> map = new HashMap<String,String>();	
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color 
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#CC9900"); //mustard
		//description
		map.put(BarParser.DESCRIPTION_TAG, "Normalized diff splice count log2Rto");

		//(String name, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes)
		Info info = new Info(workingSpliceGraphDir.getName(), genomeVersion, chrom, ".", 0, map);
		
		//convert to SmoothingWindow
		SmoothingWindow[] sm = new SmoothingWindow[exons.length];
		for (int i=0; i< exons.length; i++){
			float[] scores = {exons[i].getScore()};
			sm[i] = new SmoothingWindow(exons[i].getStart(), exons[i].getEnd(), scores);
		}

		//get heatmap positions and values
		HeatMapMakerPosNeg hm = new HeatMapMakerPosNeg(0, 0, 0);
		PointData pd = hm.makeHeatMapPositionValues(sm);
		

		pd.setInfo(info);
		pd.writePointData(workingSpliceGraphDir);
		pd.nullPositionScoreArrays();
	}

	/**Looks for minimum number of reads and minimum log2Ratio difference between exon counts.*/
	private boolean checkForMinimums(UCSCGeneLine gene){
		//get total treatment and total control
		float[][] tExonCounts = gene.getTreatmentExonCounts();
		float[][] cExonCounts = gene.getControlExonCounts();
		int numExons = tExonCounts[0].length;
		ExonIntron[] exons = gene.getExons();
		
		//collapse
		float[] tCounts = new float[numExons];
		float[] cCounts = new float[numExons];
		//for each exon
		for (int i=0; i< numExons; i++){
			//for each replica
			for (int j=0; j< tExonCounts.length; j++){
				tCounts[i] += tExonCounts[j][i];
			}
			//sum from c
			for (int j=0; j< cExonCounts.length; j++){
				cCounts[i] += cExonCounts[j][i];
			}
		}

		//need to estimate scalars
		float[] totalTC = gene.getScores();
		double scalarTC = totalTC[0]/ totalTC[1];
		double scalarCT = totalTC[1]/totalTC[0];
		
		//for each exon
		ArrayList<Integer> goodIndexes = new ArrayList<Integer>();
		ArrayList<ExonIntron> goodExonsAL = new ArrayList<ExonIntron>();
		float maxLogRatio = 0;
		int maxLogRatioIndex = 0;
		float[] log2Ratios = new float[numExons];
		for (int i=0; i< numExons; i++){
			//require minimum counts in each exon from each sample
			if (tCounts[i] < minimumSpliceCounts || cCounts[i] < minimumSpliceCounts) continue;
			//check ratio
			log2Ratios[i] = calculateLog2Ratio(tCounts[i], cCounts[i], scalarTC, scalarCT);
			float logRatio = Math.abs(log2Ratios[i]);
			if (logRatio > maxLogRatio) {
				maxLogRatio = logRatio;
				maxLogRatioIndex = i;
			}
			goodIndexes.add(new Integer(i));
			exons[i].setScore(log2Ratios[i]);
			goodExonsAL.add(exons[i]);
		}

		//check ratio and good exons, need at least two for alt splicing
		int numGoodExons = goodIndexes.size();
		if (numGoodExons < 2 || maxLogRatio < minLog2Ratio) {
			gene.setExonCounts(null);
			return false;
		}
		//set maxLog2Ratio and maxExon
		gene.setSplicingLog2Ratio(log2Ratios[maxLogRatioIndex]);
		gene.setMaxLog2RtoSplicedExon(exons[maxLogRatioIndex]);

		//set good exons, for graphing splicing
		ExonIntron[] goodExons = new ExonIntron[goodExonsAL.size()];
		goodExonsAL.toArray(goodExons);
		gene.setScoredExons(goodExons);

		//set exon counts
		if (numGoodExons == numExons) {
			gene.setExonCounts(new float[][]{tCounts, cCounts});
			return true;
		}
		//nope need to reset the counts
		float[] tExonCountsSub = new float[numGoodExons];
		float[] cExonCountsSub = new float[numGoodExons];

		for (int i=0; i< numGoodExons; i++){
			int goodIndex = goodIndexes.get(i);
			tExonCountsSub[i] = tCounts[goodIndex];
			cExonCountsSub[i] = cCounts[goodIndex];
		}
		gene.setExonCounts(new float[][]{tExonCountsSub, cExonCountsSub});
		return true;

	}

	private void analyzeUsingSamSeq(){
		ArrayList<PairedCondition> pairedConditionsAL = new ArrayList<PairedCondition>();
		try {
			StringBuilder sb = new StringBuilder();
			sb.append("library(samr)\n");

			//diff express
			pairedConditionsAL = new ArrayList<PairedCondition>();
			for (int i=0; i< conditions.length; i++){
				for (int j=i+1; j< conditions.length; j++){
					//make condition
					PairedCondition pc = new PairedCondition(conditions[i], conditions[j]);
					pairedConditionsAL.add(pc);
					File results = new File (saveDirectory, pc.getName()+"DESeq2.txt");
					if (deleteTempFiles) results.deleteOnExit();
					pc.setDiffExpResults(results);

					//make R script


					//Load up data
					sb.append("dataTable <- read.table('" + geneCountTable.getCanonicalPath()+"', header=TRUE,row.names=1)\n");
					sb.append("geneNames <- rownames(dataTable)\n");
					sb.append("conds <- as.factor(c('"+ Misc.stringArrayToString(conditionNamesPerReplica, "','") + "'))\n");

					//Create data model
					sb.append("dataModel <- list(x=dataTable,y=conds,geneid=geneNames)\n");

					//Run SAMSeq
					sb.append("results <- SAMseq(dataTable,conds,resp.type='Two class unpaired',geneid=geneNames,fdr.output=1)\n");

					//Create output table
					sb.append("outputTables <- samr.compute.siggenes.table(results$samr.obj,results$del,dataModel,results$delta.table,all.genes=TRUE)\n");

					//Merge output tables
					sb.append("mergeTable <- rbind(outputTables$genes.up,outputTables$genes.lo)\n");

					//Clean up tables
					sb.append("finalTable <- mergeTable[,c(3,4,7,8)]\n");
					sb.append("finalTable[,4] <- as.numeric(finalTable[,4]) / 100\n");
					sb.append("finalTable[,3] <- log(as.numeric(finalTable[,3]),2)\n");
					sb.append("finalTable <- finalTable[match(geneNames,finalTable[,1]),]\n");

					//Write Table
					sb.append("write.table(finalTable,file='"+results.getCanonicalPath()+"',sep='\t',quote=FALSE,row.names=FALSE)");

				}
			}

			//write script to file
			File scriptFile = new File (saveDirectory,"samseq_RScript.txt");
			File rOut = new File(saveDirectory, "samseq_RScript.txt.Rout");
			IO.writeString(sb.toString(), scriptFile);

			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};

			//execute command
			IO.executeCommandLine(command);

			//check for warnings (Not sure what to check at this point)
			//String[] res = IO.loadFile(rOut);

			//Make conditions
			pairedConditions = new PairedCondition[pairedConditionsAL.size()];
			pairedConditionsAL.toArray(pairedConditions);

			//look for results files and parse results
			for (PairedCondition pc: pairedConditions){
				if (pc.getDiffExpResults().exists() == false ) throw new IOException("\n\nR results file doesn't exist. Check temp files in save directory for error.\n");
				pc.parseSamSeqStatResults(geneNamesToAnalyze, minAdjP, minLog2Ratio);
				if (verbose) System.out.println("\t"+pc.getVsName()+"\t"+pc.getDiffExpGeneNames().size());
				geneNamesPassingThresholds.addAll(pc.getDiffExpGeneNames());
			}

			//cleanup
			if (deleteTempFiles) {
				rOut.deleteOnExit();
				scriptFile.deleteOnExit();
			}

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error: failed to execute SAMSeq.\n");
		}

	}


	public void runDESeq2(){
		ArrayList<PairedCondition> pairedConditionsAL = new ArrayList<PairedCondition>();
		
		//make File objects
		rLogValues = new File(saveDirectory, "sampleRlogValues.txt");
		if (deleteTempFiles) rLogValues.deleteOnExit();
		File sampleClusterPlot = new File (saveDirectory, "sampleClusterPlot.pdf");
		
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("#load libraries\n");
			sb.append("library(DESeq2)\n");
			sb.append("library(gplots)\n");
			sb.append("library(RColorBrewer)\n");
			
			sb.append("\n#load count table, replace rownames with first column\n");
			sb.append("countTable = read.delim('"+geneCountTable.getCanonicalPath()+"', header=TRUE)\n");
			sb.append("rownames(countTable) = countTable[,1]\n");
			sb.append("countTable = countTable[,-1]\n");
			
			sb.append("\n#make info object describing different conditions, swap out row names with actual names from table\n");
			sb.append("sampleInfo = data.frame(condition=as.factor(c('"+ Misc.stringArrayToString(conditionNamesPerReplica, "','") + "')))\n");
			sb.append("rownames(sampleInfo) = colnames(countTable)\n");
			
			sb.append("\n#run DESeq2, this might throw warnings but should restart itself\n");
			sb.append("cds = DESeqDataSetFromMatrix(countData=countTable, colData=sampleInfo, design = ~condition)\n");
			sb.append("cds = DESeq(cds)\n");
			
			sb.append("\n#get and save rlog transformed data\n");
			sb.append("rld = rlog(cds)\n");
			sb.append("write.table(assay(rld), file = '"+rLogValues.getCanonicalPath()+"', quote=FALSE, sep ='\t')\n");
			
			sb.append("\n#make sample heat map\n");
			sb.append("sampleDists = dist(t(assay(rld)))\n");
			sb.append("sampleDistMatrix <- as.matrix( sampleDists )\n");
			sb.append("colours = colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)\n");
			sb.append("pdf('"+sampleClusterPlot.getCanonicalPath()+"')\n");
			sb.append("heatmap.2( sampleDistMatrix, trace='none', col=colours)\n");
			sb.append("dev.off()\n");
			
			sb.append("\n#for each pairing, order and number of rows same in output as in starting count table\n");
			for (int i=0; i< conditions.length; i++){
				for (int j=i+1; j< conditions.length; j++){
					//make condition
					PairedCondition pc = new PairedCondition(conditions[i], conditions[j]);
					pairedConditionsAL.add(pc);
					File results = new File (saveDirectory, pc.getName()+"DESeq2.txt");
					if (deleteTempFiles) results.deleteOnExit();
					pc.setDiffExpResults(results);
					sb.append("res = results(cds, contrast = c('condition', '"+conditions[i].getName()+"', '"+conditions[j].getName()+"'))\n");
					sb.append("res[,6] = -10 * log10(res[,6])\n"); //transform pvals since java can't handle some of these big numbers
					sb.append("write.table(res, file = '"+results.getCanonicalPath()+"', quote=FALSE, sep ='\t')\n");
				}
			}

			//write script to file
			File scriptFile = new File (saveDirectory,"deseq2_RScript.txt");
			File rOut = new File(saveDirectory, "deseq2_RScript.txt.Rout");
			IO.writeString(sb.toString(), scriptFile);

			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};

			//execute command
			IO.executeCommandLine(command);

			//check for warnings?
			String[] res = IO.loadFile(rOut);

			/*
			if (res.length !=0) System.out.println("\nMessages from R...");
			for (String s: res) {
				System.out.println(s);
				
				if (s.contains("did not converge") || s.contains("fit failed")) {
					System.err.println("\n\t\tWarning, DESeq2's GLM dispersion fit failed. Relaunching using fitType='local'");
					analyzeForDifferentialExpression(true);
					return;
				
			}
			}*/

			//Make conditions
			pairedConditions = new PairedCondition[pairedConditionsAL.size()];
			pairedConditionsAL.toArray(pairedConditions);

			//look for results files and parse results
			for (PairedCondition pc: pairedConditions){
				if (pc.getDiffExpResults().exists() == false ) throw new IOException("\n\nR results file doesn't exist. Check temp files in save directory for error.\n");
				pc.parseDESeq2Results(geneNamesToAnalyze, minAdjP, minLog2Ratio);
				if (verbose) System.out.println("\t"+pc.getVsName()+"\t"+pc.getDiffExpGeneNames().size());
				geneNamesPassingThresholds.addAll(pc.getDiffExpGeneNames());
			}

			//cleanup
			if (deleteTempFiles) {
				rOut.deleteOnExit();
				scriptFile.deleteOnExit();
			}

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error: failed to execute DESeq2.\n");
		}

	}

	/**Calculates a log2( (tSum+1)/(cSum+1) ) on linearly scaled tSum and cSum based on the total observations.*/
	public float calculateLog2Ratio( double tSum, double cSum, double scalarTC, double scalarCT){
		double t;
		double c;
		if (tSum !=0 ) {
			t = tSum * scalarCT;
			c = cSum;
		}
		else {
			c = cSum * scalarTC;
			t = tSum;
		}
		double ratio = (t+1)/(c+1);
		return (float)Num.log2(ratio);
	}

	/**Prints a variety of spreadsheets.*/
	public void printStatSpreadSheet(){
		try {
			//create a workbook
			workbook = new XSSFWorkbook();
			FileOutputStream fileOut = new FileOutputStream(new File(saveDirectory, "geneStats.xlsx"));

			//create a worksheet
			Sheet sheet0 = workbook.createSheet("Analyzed Genes");
			sheet0.createFreezePane( 0, 1, 0, 1 );

			//create rows, one for each gene with counts plus header
			Row[] rows = new Row[geneNamesToAnalyze.length+1];
			for (int i=0; i< rows.length; i++) rows[i] = sheet0.createRow(i);

			//make header line
			makeHeaderLine(rows[0], true);

			//add gene info
			int cellIndex = addGeneInfo(rows);

			//add PairedCondition info
			cellIndex = addPairedConditionInfo(rows, cellIndex);

			//add count values
			cellIndex = addCountInfo(rows, cellIndex);

			if (!this.useSamSeq) {
				//add rLog values
				cellIndex = allRLogValues(rows, cellIndex);
			}

			//add FPKM
			cellIndex = addFPKMInfo(rows, cellIndex);

			//new sheet for diff expressed genes at relaxed thresholds
			Sheet sheet1 = workbook.createSheet("AdjP 10 Lg2Rt 0.585");
			sheet1.createFreezePane( 0, 1, 0, 1 );
			addDiffExpressedGenes(sheet1, 10f, 0.5849625f);

			//new sheet for diff expressed genes at standard thresholds
			Sheet sheet2 = workbook.createSheet("AdjP 20 Lg2Rt 1");
			sheet2.createFreezePane( 0, 1, 0, 1 );
			addDiffExpressedGenes(sheet2, 20f, 1f);

			//new sheet for diff expressed genes at standard thresholds
			Sheet sheet3 = workbook.createSheet("AdjP 30 Lg2Rt 1.585");
			sheet3.createFreezePane( 0, 1, 0, 1 );
			addDiffExpressedGenes(sheet3, 30f, 1.584963f);

			//new sheet for flagged genes
			Sheet sheet4 = workbook.createSheet("Flagged Genes");
			addFlaggedGenes(sheet4);

			//new sheet for descriptions of headings
			Sheet sheet5= workbook.createSheet("Column Descriptions");
			addColumnDescriptors(sheet5);

			//write out and close
			workbook.write(fileOut);
			fileOut.close();


		} catch (Exception e){
			System.err.println("\nProblem printing spreadsheet!");
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void addDiffExpressedGenes(Sheet sheet, float minAdjP, float minLog2Ratio) {
		//add header row
		makeHeaderLine(sheet.createRow(0), false);

		//for each paired condition
		int columnIndex = 0;
		for (int i=0; i< pairedConditions.length; i++){
			PairedCondition pc = pairedConditions[i];
			float[][] geneScores = pc.getParsedDiffExpResults();

			//for each gene find those that pass
			ArrayList<GeneResult> passing = new ArrayList<GeneResult>();
			for (int j=0; j< geneNamesToAnalyze.length; j++){
				//fdr and log2 diff exp
				float[] fdrLog2 = geneScores[j];
				//pass?
				if (Math.abs(fdrLog2[1]) >= minLog2Ratio && fdrLog2[0] >= minAdjP){
					passing.add(new GeneResult(geneNamesToAnalyze[j], fdrLog2[0], fdrLog2[1]));
				}
			}

			//convert to array and sort by abs of log2
			GeneResult[] gr = new GeneResult[passing.size()];
			passing.toArray(gr);
			Arrays.sort(gr);

			//for each gene result
			int rowIndex = 1;
			for (int j=0; j< gr.length; j++){
				Row r = sheet.getRow(rowIndex);
				if (r == null) r = sheet.createRow(rowIndex);
				gr[j].addInfo(r, columnIndex);
				//advance row
				rowIndex++;
			}

			//advance column index
			columnIndex+=3;
		}

	}

	private void addColumnDescriptors(Sheet sheet) {
		//make some styles
		CellStyle bold = workbook.createCellStyle();
		Font font = workbook.createFont();
		font.setBoldweight(Font.BOLDWEIGHT_BOLD);
		bold.setFont(font);

		int rowIndex = 0;

		String[][] usedDescriptors;

		if (useSamSeq) {
			usedDescriptors = this.samSeqDescriptors;
		} else {
			usedDescriptors = this.deSeqDescriptors;
		}

		for (int i=0; i< usedDescriptors.length; i++){
			Row row = sheet.createRow(rowIndex++);
			Cell cell = row.createCell(0);
			cell.setCellStyle(bold);
			cell.setCellValue(usedDescriptors[i][0]);
			if (usedDescriptors[i][1] != null) row.createCell(1).setCellValue(usedDescriptors[i][1]);
		}

		sheet.autoSizeColumn((short)0);
	}

	private String[][] deSeqDescriptors = {
			{"Gene/Region info:", null},
			{"IGB HyperLink", "Click each to reposition a running instance of IGB, http://bioviz.org/igb/"},
			{"Alt Name", "Alternative name for this gene/region"},
			{"Coordinates", "Chromosome : Start of the gene/region - end of the gene/region"},
			{"Strand", "Strand of the gene/region"},
			{"Total BPs", "Total base pairs of gene/region after masking for overlapping annotations if so indicated"},
			{"Max Abs Lg2Rto", "Maximum absolute log2 ratio observed"},
			{"Max DESeq2 AdjP", "Maximum observed -10Log10( AdjP ) between any pairwise DESeq2 comparisons.  Note the MaxAbsLg2Rto and MaxDESeqAdjP may be from different comparisons.  "
					+ "Moreover, all p-values and AdjPs in USeq output, both printed and in graph form have been phred transformed where the value is actually -10*Log10(AdjP or pval). Thus -10*Log10(0.001) = 30, -10*Log10(0.01) = 20, -10*Log10(0.05) = 13, and -10*Log10(0.1) = 10. "},
					{"", ""},
					{"For each Paired Comparison:", null},
					{"Lg2Rto", "DESeq2's log2 ratio for differential expression."},
					{"AdjP", "DESeq2's negative binomial -10Log10(adjusted p-value) following the Benjamini and Hochberg multiple testing correction."},
					{"Spli Lg2Rto", "The maximum log2 ratio difference between two conditions and a gene's exon after correcting for differences in counts.  Use this to gauge the degree of potential differential splicing."},
					{"Spli AdjPVal", "Chi-square test of independence between two conditions and a gene's exons that contain 5 or more counts.  A Bonferroni correction is applied to the p-values."},
					{"Spli Coor", "Coordinates of the exon displaying the SpliLg2Rto."},
					{"", ""},
					{"For each Replica:", null},
					{"Counts", "Number of alignments that overlap a given gene/region by at least one base pair. Paired alignments are only counted once.  Alignments that span a gene/region but don't actually touch down are ignored. "},
					{"RLog", "Per gene/region regularized log counts from DESeq2.  Use these for subsequent pca and clustering analysis."},
					{"FPKM", "Transformed count data normalized to library size and gene/region length (counts/ exonicBasesPerKB/ millionTotalMappedReadsToGeneTable)"}

	};

	private String[][] samSeqDescriptors = {
			{"Gene/Region info:", null},
			{"IGB HyperLink", "Click each to reposition a running instance of IGB, http://bioviz.org/igb/"},
			{"Alt Name", "Alternative name for this gene/region"},
			{"Coordinates", "Chromosome : Start of the gene/region - end of the gene/region"},
			{"Strand", "Strand of the gene/region"},
			{"Total BPs", "Total base pairs of gene/region after masking for overlapping annotations if so indicated"},
			{"Max Abs Lg2Rto", "Maximum absolute log2 ratio observed"},
			{"Max SAMseq AdjP", "Maximum observed -10Log10( AdjP ) between any pairwise DESeq2 comparisons.  Note the MaxAbsLg2Rto and MaxDESeqAdjP may be from different comparisons.  "
					+ "Moreover, all p-values and AdjPs in USeq output, both printed and in graph form have been phred transformed where the value is actually -10*Log10(AdjP or pval). Thus -10*Log10(0.001) = 30, -10*Log10(0.01) = 20, -10*Log10(0.05) = 13, and -10*Log10(0.1) = 10. "},
					{"", ""},
					{"For each Paired Comparison:", null},
					{"Lg2Rto", "SAMseq's log2 ratio for differential expression"},
					{"AdjP", "SAMseq's -10Log10(AdjP)"},
					{"Spli Lg2Rto", "The maximum log2 ratio difference between two conditions and a gene's exon after correcting for differences in counts.  Use this to gauge the degree of potential differential splicing."},
					{"Spli AdjPVal", "Chi-square test of independence between two conditions and a gene's exons that contain 5 or more counts.  A Bonferroni correction is applied to the p-values."},
					{"Spli Coor", "Coordinates of the exon displaying the SpliLg2Rto."},
					{"", ""},
					{"For each Replica:", null},
					{"Counts", "Number of alignments that overlap a given gene/region by at least one base pair. Paired alignments are only counted once.  Alignments that span a gene/region but don't actually touch down are ignored. "},
					{"FPKM", "Transformed count data normalized to library size and gene/region length (counts/ exonicBasesPerKB/ millionTotalMappedReadsToGeneTable)"}

	};

	private void addFlaggedGenes(Sheet sheet) {
		//generate gene lists
		String[] tooMany = Misc.hashSetToStringArray(flaggedGeneNames);
		Arrays.sort(tooMany);
		ArrayList<String> tooFewAL = new ArrayList<String>();
		for (UCSCGeneLine gene: genes){
			String name = gene.getDisplayNameThenName();
			if (flaggedGeneNames.contains(name) || geneNamesWithMinimumCounts.contains(name)) continue;
			else tooFewAL.add(name);
		}
		String[] tooFew = Misc.stringArrayListToStringArray(tooFewAL);
		Arrays.sort(tooFew);
		int maxNum = tooFew.length;
		if (tooMany.length > maxNum) maxNum = tooMany.length;	

		//bold styles
		CellStyle bold = workbook.createCellStyle();
		Font font = workbook.createFont();
		font.setBoldweight(Font.BOLDWEIGHT_BOLD);
		bold.setFont(font);

		int rowIndex = 0;
		Row header = sheet.createRow(rowIndex++);
		header.createCell(0).setCellValue("< "+minimumCounts+ " counts");
		header.createCell(1).setCellValue("> "+maxAlignmentsDepth+ " counts");
		header.getCell(0).setCellStyle(bold);
		header.getCell(1).setCellStyle(bold);

		for (int i=0; i< maxNum; i++){
			String tooManyCell = "";
			String tooFewCell = "";
			if (tooMany.length > i) tooManyCell = tooMany[i];
			if (tooFew.length > i) tooFewCell = tooFew[i];

			Row r = sheet.createRow(rowIndex++);
			r.createCell(0).setCellValue(tooFewCell);
			r.createCell(1).setCellValue(tooManyCell);
		}
		//autosize
		sheet.autoSizeColumn(0);
		sheet.autoSizeColumn(1);
	}

	private int addCountInfo(Row[] rows, int cellIndex) {
		int ci = cellIndex;

		//for each gene / row
		for (int j=0; j< geneNamesToAnalyze.length; j++){
			int rowIndex = j+1;
			ci = cellIndex;
			for (Condition c: conditions){
				//for each replica
				for (Replica r: c.getReplicas()){
					int num = 0;
					if (r.getGeneCounts().containsKey(geneNamesToAnalyze[j])){
						num = r.getGeneCounts().get(geneNamesToAnalyze[j]).getCount();
					}
					rows[rowIndex].createCell(ci++).setCellValue(num);
				}
			}

		}
		return ci;
	}

	private int allRLogValues(Row[] rows, int cellIndex) {
		//First parse the table of numbers
		float[][] geneRLog = parseRLogData();

		//write them to the excel sheet
		int ci = cellIndex;
		for (int j=0; j< geneNamesToAnalyze.length; j++){
			int rowIndex = j+1;
			ci = cellIndex;
			for (float c: geneRLog[j]){
				rows[rowIndex].createCell(ci++).setCellValue(c);
			}
		}

		return ci;
	}

	private int addFPKMInfo(Row[] rows, int cellIndex) {
		int ci = cellIndex;
		//for each gene
		for (int j=0; j< geneNamesToAnalyze.length; j++){
			int rowIndex = j+1;
			ci = cellIndex;
			double totalBps = name2Gene.get(geneNamesToAnalyze[j]).getTotalExonicBasePairs();
			//calculate FPKM values for each condition			
			for (Condition c: conditions){
				//for each replica
				for (Replica r: c.getReplicas()){
					double num = 0;
					if (r.getGeneCounts().containsKey(geneNamesToAnalyze[j])){
						num = r.getGeneCounts().get(geneNamesToAnalyze[j]).calculateFPKM(r.getTotalCounts(), totalBps);
					}
					rows[rowIndex].createCell(ci++).setCellValue(num);
				}
			}
		}
		return ci;
	}
	
	/*Rips the deseq rlog file checking and fixing negative values.*/
	private float[][] parseRLogData() {
		Pattern tab = Pattern.compile("\\t");
		float[][] geneRLog = new float[geneNamesToAnalyze.length][replicaNames.length];
		try {
			BufferedReader in = IO.fetchBufferedReader(rLogValues);
			//skip header 
			in.readLine();
			int numFields = replicaNames.length +1;
			for (int i=0; i< geneNamesToAnalyze.length; i++){
				String line = in.readLine();
				String[] scoresStringWithName = tab.split(line);
				//check number of scores
				if (scoresStringWithName.length != numFields) throw new Error("\nError: one of the rLog data rows doesn't contain the appropriate number of values! See row "+line);
				//check name
				if (scoresStringWithName[0].equals(geneNamesToAnalyze[i]) == false) throw new Error("\nError: one of the rLog data rows name doesn't match expected! See row "+line);
				//parse scores
				String[] scoresString = new String[replicaNames.length];
				System.arraycopy(scoresStringWithName, 1, scoresString, 0, scoresString.length);
				float[] scores = Num.parseFloats(scoresString);
				if (scores == null) throw new Error("\nError: couldn't parse floats from one of the rLog data rows! See row "+line);
				geneRLog[i] = scores;
			}
			in.close();
		} catch (Exception e) {
			Misc.printErrAndExit("\nError: problem parsing rLog data\n"+e.getMessage());
		}
		return geneRLog;
	}

	/*Adds DiffExpLg2Rto DiffExpAdjP DiffSpliceLg2Rto DiffSpliceAdjP DiffSpliceCoordinates*/
	private int addPairedConditionInfo(Row[] rows, int cellIndex) {
		float[] noScores = {0f,0f,0f,0f,0f};
		int ci = 0;

		//for each paired condition
		for (int i=0; i< pairedConditions.length; i++){
			PairedCondition pc = pairedConditions[i];
			float[][] geneScores = pc.getParsedDiffExpResults();
			HashMap<String, float[]> spliceScores = pc.getGeneNameSpliceStats();

			//for each gene / row
			for (int j=0; j< geneNamesToAnalyze.length; j++){
				int rowIndex = j+1;
				ci = cellIndex;
				//fdr and log2 diff exp
				float[] fdrLog2 = geneScores[j];
				//fdr, lg2, index
				float[] splice = noScores;
				if (spliceScores != null && spliceScores.containsKey(geneNamesToAnalyze[j]))  splice = spliceScores.get(geneNamesToAnalyze[j]); 
				//lg2Rto deseq
				rows[rowIndex].createCell(ci++).setCellValue(fdrLog2[1]);
				//fdr deseq
				rows[rowIndex].createCell(ci++).setCellValue(fdrLog2[0]);
				//lg2 splice
				rows[rowIndex].createCell(ci++).setCellValue(splice[1]);
				//fdr splice
				rows[rowIndex].createCell(ci++).setCellValue(splice[0]);
				//coordinates of max log2 rto spliced exon
				String coor = "-";
				if (splice[2] !=0 || splice[3] !=0){
					UCSCGeneLine gene = name2Gene.get(geneNamesToAnalyze[j]);
					coor = gene.getChrom()+":"+(int)splice[2]+"-"+(int)splice[3];
				}
				rows[rowIndex].createCell(ci++).setCellValue(coor);
			}
			cellIndex = ci;
		}
		return ci;
	}

	/*Adds DisplayName Name Chr Strand Start Stop TotalBPs MaxAbsLg2Rto MaxDESeqAdjP*/
	private int addGeneInfo(Row[] rows) {
		CreationHelper createHelper = workbook.getCreationHelper();

		//make style for hyperlinks
		CellStyle hlStyle = workbook.createCellStyle();
		Font hlFont = workbook.createFont();
		hlFont.setUnderline(Font.U_SINGLE);
		hlFont.setColor(IndexedColors.BLUE.getIndex());
		hlStyle.setFont(hlFont);

		int cellIndex = 0;
		for (int i=0; i < geneNamesToAnalyze.length; i++){
			cellIndex = 0;
			int rowIndex = i+1;
			UCSCGeneLine gene = name2Gene.get(geneNamesToAnalyze[i]);

			//hyperlink with displayName
			Hyperlink link = createHelper.createHyperlink(Hyperlink.LINK_URL);
			link.setAddress(fetchIGBHyperLink(gene));
			link.setLabel(gene.getDisplayNameThenName());

			Cell hlCell = rows[rowIndex].createCell(cellIndex++);
			hlCell.setCellStyle(hlStyle);
			hlCell.setHyperlink(link);
			hlCell.setCellValue(gene.getDisplayNameThenName());

			//alt name
			if (gene.getName() != null) rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getName());
			else rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getDisplayName());
			//chr:start-stop
			rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getChrStartStop());
			//Strand
			rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getStrand());
			//TotalBPs
			rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getTotalExonicBasePairs());
			//MaxAbsLg2Rto
			rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getMaxAbsLog2Ratio());
			//MaxDESeqAdjP
			rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getFdr());
		}
		return cellIndex;
	}

	private String fetchIGBHyperLink(UCSCGeneLine gene){
		StringBuilder text = new StringBuilder();
		//url
		int start = gene.getTxStart() - 50000;
		if (start < 0) start = 0;
		int end = gene.getTxEnd() + 50000;

		text.append(url);
		text.append(gene.getChrom());
		text.append("&start=");
		text.append(start);
		text.append("&end=");
		text.append(end);
		return text.toString();
	}

	private int createCell (Row row, String value, int index, CellStyle style){
		Cell cell = row.createCell(index);
		cell.setCellType(Cell.CELL_TYPE_STRING);
		if (style != null) cell.setCellStyle(style);
		cell.setCellValue(value);
		return index +1;
	}

	private CellStyle createHeaderStyle(short colorIndex){
		CellStyle style = workbook.createCellStyle();
		style.setWrapText(true);
		style.setAlignment(CellStyle.ALIGN_CENTER);
		Font font = workbook.createFont();
		font.setBoldweight(Font.BOLDWEIGHT_BOLD);
		style.setFont(font);
		if (colorIndex != -1){
			style.setFillForegroundColor(colorIndex);
			style.setFillPattern(CellStyle.SOLID_FOREGROUND);
		}
		return style;
	}
	private void makeHeaderLine(Row row, boolean all) {
		//create header
		CellStyle headerCellStyle = createHeaderStyle((short)-1);

		int index = 0;
		if (all){
			index = createCell(row, "IGB HyperLink", index, headerCellStyle);
			index = createCell(row, "Alt Name", index, headerCellStyle);
			index = createCell(row, "Coor "+genomeVersion, index, headerCellStyle);
			index = createCell(row, "Strand", index, headerCellStyle);
			index = createCell(row, "Total BPs", index, headerCellStyle);
			index = createCell(row, "Max Abs Lg2Rto", index, headerCellStyle);
			index = createCell(row, "Max DESeq2 AdjP", index, headerCellStyle);
		}

		//for each PairedCondition
		int colorIndex = 0;
		short[] colorIndexes = new short[]{IndexedColors.LIGHT_BLUE.getIndex(), IndexedColors.LIGHT_GREEN.getIndex(), IndexedColors.LIGHT_ORANGE.getIndex(), IndexedColors.LIGHT_TURQUOISE.getIndex(), IndexedColors.LIGHT_YELLOW.getIndex()};
		for (PairedCondition pc : pairedConditions){
			CellStyle pcCellStyle = createHeaderStyle(colorIndexes[colorIndex++]);
			if (colorIndex == colorIndexes.length) colorIndex = 0;
			String name = pc.getName();
			if (all == false) index = createCell(row, name, index, pcCellStyle);
			index = createCell(row, "Lg2Rto "+name, index, pcCellStyle);
			index = createCell(row, "AdjP "+name, index, pcCellStyle);
			if (all){
				index = createCell(row, "Spli Lg2Rto "+name, index, pcCellStyle);
				index = createCell(row, "Spli AdjPVal "+name, index, pcCellStyle);
				index = createCell(row, "Spli Coor "+name, index, pcCellStyle);
			}
		}

		if (all){
			//count values for each replica
			CellStyle countCellStyle = createHeaderStyle(IndexedColors.LIGHT_CORNFLOWER_BLUE.getIndex());
			for (String repName: replicaNames){
				index = createCell(row, "Counts "+repName, index, countCellStyle);
			}

			if (!this.useSamSeq) {
				//rLog values for each replica
				CellStyle rLogCellStyle = createHeaderStyle(IndexedColors.WHITE.getIndex());
				for (String repName: replicaNames){
					index = createCell(row, "RLog "+repName, index, rLogCellStyle);
				}
			}

			//FPKM values for each replica	    
			CellStyle fpkmCellStyle = createHeaderStyle(IndexedColors.SKY_BLUE.getIndex());
			for (String repName: replicaNames){
				index = createCell(row, "FPKM "+repName, index, fpkmCellStyle);
			}
		}
	}

	public void writeCountTable(){
		ArrayList<String> conditionNamesPerReplicaAL = new ArrayList<String>();
		ArrayList<String> replicaNamesAL = new ArrayList<String>();
		geneCountTable = new File(saveDirectory, "geneCountTableMin"+minimumCounts+".txt");
		if (deleteTempFiles) geneCountTable.deleteOnExit();

		try {
			//write matrix of name, t,t,t...c,c,c,c...d,d,d, to file for genes with observations
			PrintWriter out = new PrintWriter( new FileWriter(geneCountTable));

			//print header
			out.print("GeneName");
			//for each condition
			for (Condition c: conditions){
				//for each replica
				String conditionName = c.getName();
				for (Replica r: c.getReplicas()){
					out.print("\t");
					out.print(r.getNameNumber());
					conditionNamesPerReplicaAL.add(conditionName);
					replicaNamesAL.add(r.getNameNumber());
				}
			}
			out.println();

			//print counts for genes with minimal counts in any replica
			for (String geneName: geneNamesToAnalyze){
				out.print(geneName);
				for (Condition c: conditions){
					//for each replica
					for (Replica r: c.getReplicas()){
						out.print("\t");
						//printing 5' and 3' counts or standard?
						if (printFirstLastCountTable){
							String num = "0:0";
							if (r.getGeneCounts().get(geneName) != null) {
								int[] fiveThree = r.getGeneCounts().get(geneName).getExonCounts();
								num = fiveThree[0]+":"+ fiveThree[1];
							}
							out.print(num);
						}
						else {
							int num = 0;
							//remember that genes with zero counts are not stored in a replica
							if (r.getGeneCounts().get(geneName) != null) num = r.getGeneCounts().get(geneName).getCount();
							out.print(num);
						}
					}
				}
				out.println();
			}
			out.close();
		}
		catch (Exception e){
			System.err.println("Problem writing out count table.");
			e.printStackTrace();
		}

		conditionNamesPerReplica = Misc.stringArrayListToStringArray(conditionNamesPerReplicaAL);
		replicaNames = Misc.stringArrayListToStringArray(replicaNamesAL);
		
		//exit?
		if (printFirstLastCountTable) System.exit(0);
	}

	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = null;
		if (refSeqFile != null){
			reader = new UCSCGeneModelTableReader(refSeqFile, 0);
			//just coding?
			if (trimUTRBPsFromExons) reader.trimExonsOfUTRBPs();
			//introns?
			if (scoreIntrons){
				if (verbose) System.out.println("\tParsing/scoring introns instead of exons from gene models.");
				reader.swapIntronsForExons();
			}
			//remove overlap
			if (removeOverlappingRegions) {
				if (verbose) System.out.print("\tRemoving overlapping regions from gene models");
				String deletedGenes = reader.removeOverlappingExons();
				int numDelGenes = deletedGenes.split(",").length;
				if (deletedGenes.length() !=0) {
					File deleted = new File (saveDirectory, Misc.removeExtension(refSeqFile.getName())+"_TrimmedDeletedGenes.txt") ;
					if (verbose) System.out.println("\tWARNING: "+ numDelGenes +" genes had more than 1/2 of their exonic bps removed. See "+ deleted);
					IO.writeString(deletedGenes, deleted);
				}
				//write out gene table sans overlaps
				File f = new File (saveDirectory, Misc.removeExtension(refSeqFile.getName())+"_NoOverlappingExons.txt") ;
				if (verbose) System.out.println("\tWrote the trimmed gene table to the save directory. Use this table to speed up subsequent processing.");
				reader.writeGeneTableToFile(f);
			}
			genes = reader.getGeneLines();
		}
		//or from bed file
		else if (bedFile != null) {
			Bed[] bed = Bed.parseFile(bedFile, 0, 0);
			genes = new UCSCGeneLine[bed.length];
			boolean addName = bed[0].getName().trim().equals("");
			for (int i=0; i< bed.length; i++){
				if (addName) bed[i].setName((i+1)+"");
				genes[i] = new UCSCGeneLine(bed[i]);
			}
			reader = new UCSCGeneModelTableReader();
			reader.setGeneLines(genes);
		}
		if (genes == null || genes.length == 0) Misc.printExit("\nProblem loading your USCS gene model table or bed file? No genes/ regions?\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check gene name is unique
		if (reader.uniqueGeneNames() == false) Misc.printExit("\nDuplicate gene names were found in your gene / bed file, these must be unique.\n");
		//check that genes are stranded
		if (performStrandedAnalysis && reader.checkStrand() == false) Misc.printExit("\nError: you have indicated to perform a stranded analysis yet one or more of your genes/ regions is unstranded, aborting.\n");
		chromGenes = reader.getChromSpecificGeneLines();
		name2Gene = new HashMap<String, UCSCGeneLine>();
		for (UCSCGeneLine gl: genes) {
			String name = gl.getDisplayNameThenName();
			name2Gene.put(name, gl);
		}
	}


	/**Loads the geneNamesWithMinimumCounts hash with gene names that pass the minimumCounts array.
	 * No need to call if scanGenes called.*/
	private void loadMinimumCountsHash(){
		geneNamesWithMinimumCounts.clear();
		//for each condition
		for (Condition c: conditions){
			//for each replica
			for (Replica r: c.getReplicas()){
				HashMap<String, GeneCount> counts = r.getGeneCounts();
				for (String geneName: counts.keySet()){
					if (counts.get(geneName).getCount() >= minimumCounts) geneNamesWithMinimumCounts.add(geneName);
				}
			}
		}
	}

	public static void main(String[] args) {

		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new DefinedRegionDifferentialSeq(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': conditionDirectories = IO.extractOnlyDirectories(new File (args[++i])); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'e': minimumCounts = Integer.parseInt(args[++i]); break;
					case 'm': removeOverlappingRegions = true; break;
					case 'n': maxNumAlignments = Integer.parseInt(args[++i]); break;
					case 'x': maxAlignmentsDepth = Integer.parseInt(args[++i]); break;
					case 'p': performStrandedAnalysis = true; break;
					case 'i': scoreIntrons = true; break;
					case 't': deleteTempFiles = false; break;
					case 'g': genomeVersion = args[++i]; break;
					case 'j': performReverseStrandedAnalysis = true; break;
					case 'k': secondStrandFlipped = true; break;
					case 'a': useSamSeq = true; break;
					//hidden options
					case 'o': trimUTRBPsFromExons = true; break;
					case 'z': printFirstLastCountTable = true; break;
					
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}

		//fetch genome version and make url
		if (genomeVersion == null) Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");
		url = "http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";

		//look for bam files
		if (conditionDirectories == null || conditionDirectories.length == 0) Misc.printErrAndExit("\nError: cannot find any condition directories?\n");
		if (conditionDirectories.length < 2) Misc.printErrAndExit("\nError: must provide at least two Conditions for analysis.\n");
		for (File dir: conditionDirectories){
			File[] bamFiles = IO.extractFiles(dir, ".bam");
			OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, true);
		}

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdir();

		//check for R and required libraries, don't need it if they just want the first and last 1/3 count table
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else if (printFirstLastCountTable == false) {
			String errors = IO.runRCommandLookForError("library(DESeq2); library(gplots)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printErrAndExit("\nError: Cannot find the required R library.  Did you install DESeq2 " +
						"(http://www-huber.embl.de/users/anders/DESeq2/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(DESeq2)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
			}
		}

		//Make sure stranded and reverse aren't both set
		if (performStrandedAnalysis == true && performReverseStrandedAnalysis == true) {
			Misc.printErrAndExit("\nError: Stranded and reverse-stranded analyses both set to true, pick one and restart\n\n");
		}

		//look for bed file
		if (refSeqFile == null && bedFile == null){
			Misc.printErrAndExit("\nPlease enter a regions or gene file to use in scoring regions.\n");
		}
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Defined Region Differential Seq:   Aug 2014                  **\n" +
				"**************************************************************************************\n" +
				"DRDS takes sorted bam files, one per replica, minimum one per condition, minimum two\n" +
				"conditions (e.g. treatment and control or a time course/ multiple conditions) and\n" +
				"identifies differentially expressed genes using DESeq2 or SAMTools. DESeq2's rLog\n" +
				"normalized count data is used to heirachically cluster the samples. Differential\n"+
				"splicing is estimated using a chi-square test of independence. When testing only a\n"+
				"few genes or regions, append these onto a full gene table so that DESeq2 can\n"+
				"appropriately estimate the library size and replica variance.\n"+

				"\nOptions:\n"+
				"-s Save directory.\n"+
				"-c Conditions directory containing one directory for each condition with one xxx.bam\n" +
				"       file per biological replica and their xxx.bai indexs. 3-4 reps recommended per\n" +
				"       condition. The BAM files should be sorted by coordinate using Picard's SortSam.\n" +
				"       All spice junction coordinates should be converted to genomic coordinates, see\n" +
				"       USeq's SamTranscriptomeParser.\n" +
				"-r Full path to R (version 3+) loaded with DESeq2, samr, and gplots defaults to\n"+
				"       '/usr/bin/R' file, see http://www.bioconductor.org . Type 'library(DESeq2);\n"+
				"       library(samr); library(gplots)' in R to see if they are installed. \n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. Tab delimited, see RefSeq Genes\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 . NOTE:\n" +
				"       this table should contain only ONE composite transcript per gene (e.g. use\n" +
				"       Ensembl genes NOT transcripts). Use the MergeUCSCGeneTable app to collapse\n" +
				"       transcripts. See http://useq.sourceforge.net/usageRNASeq.html for details.\n"+
				"-b (Or) a bed file (chr, start, stop,...), full path, See,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+
				"-g Genome Version  (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +

				"\nAdvanced Options:\n"+
				"-m Mask overlapping gene annotations, recommended for well annotated genomes.\n"+
				"-x Max per base alignment depth, defaults to 50000. Genes containing such high\n"+
				"       density coverage are ignored.\n"+
				"-n Max number alignments per read. Defaults to 1, unique.  Assumes 'NH' tags have\n"+
				"      been set by processing raw alignments with the SamTranscriptomeProcessor.\n"+
				"-e Minimum number alignments per gene-region per replica, defaults to 10.\n"+
				"-i Score introns instead of exons.\n"+
				"-p Perform a stranded analysis. Only collect reads from the same strand as the\n" +
				"      annotation.\n" +
				"-j Reverse stranded analysis.  Only collect reads from the opposite strand of the\n" +
				"      annotation.  This setting should be used for the Illumina's strand-specific\n" +
				"      dUTP protocol.\n" +
				"-k Second read's strand is flipped. Otherwise, assumes this was not done in the \n" +
				"      SamTranscriptomeParser.\n" +
				"-t Don't delete temp files (R script, R results, Rout, etc..).\n"+
				"-a Run SAMseq in place of DESeq2.  This is only recommended with five or more\n" +
				"      replicates per condition.\n" +
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/DefinedRegionDifferentialSeq -c\n" +
				"      /Data/TimeCourse/ESCells/ -s /Data/TimeCourse/DRDS -g H_sapiens_Feb_2009\n" +
				"     -u /Anno/mergedHg19EnsemblGenes.ucsc.gz\n\n" +

				"**************************************************************************************\n");

	}
}
