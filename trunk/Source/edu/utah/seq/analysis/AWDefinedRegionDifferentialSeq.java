package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.analysis.multi.Condition;
import edu.utah.seq.analysis.multi.GeneCount;
import edu.utah.seq.analysis.multi.Replica;
import edu.utah.seq.useq.apps.Text2USeq;

/** Compares datasets for differential counts over user defined regions (e.g. gene models).  Wraps Simon Anders' DESeq package with lots of modifications. Estimates alternative splicing.
 * @author Nix
 * */
public class AWDefinedRegionDifferentialSeq {

	private boolean debugVerbose = false;

	private boolean shouldReportPValsInLogSpace = true; // By default, we report p-values and FDR in -log10 scale (where, for example, 13.01 is equivalent to a P value of 0.05). This can be changed with the "--plain" command-line flag. Added by Alex Williams at Gladstone in Apr. 2013.
	private boolean showFragAssignments = false; // By default, just do the calculations. If this is TRUE, then we actually report which gene each individual fragment was assigned to. This allows the fragment assignment in USeq to be compared to fragment assignment used by Cufflinks and HTseq-counts (and potentially other programs). Added by Alex Williams at Gladstone in Apr. 2013. The equivalent command in the program "htseq-counts" is the "--samout" command line argument.

	private BufferedWriter fragAssignBuffer;

	private File[] conditionDirectories;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private File bedFile;
	private File refSeqFile;
	private String genomeVersion;
	private float minFDR = 10; // Note the default value is 10 here (an FDR of 0.1) here! This is in the -10log10 space
	private float minLog2Ratio = 1f; // Note the default value is 1 (a 2-fold change) here! This is in LOG 2 space, not log10.
	private boolean scoreIntrons = false;
	private boolean removeOverlappingRegions = false;
	private int minimumCounts = 20; // Default value is TWENTY reads are required for a gene to be considered expressed/interesting.
	private boolean deleteTempFiles = true;
	private boolean filterOutliers = false;
	private float minimumSpliceLog2Ratio = 1f; // This CANNOT be set anywhere by the user
	private boolean performStrandedAnalysis = false;
	private boolean performReverseStrandedAnalysis = false;
	private boolean usePermutationSpliceTest = false;
	private int maxRepeats = 0;
	private boolean verbose = true; // This can't be set directly on the command line, but it can be set in the constructor when called from another file.
	private int maxAlignmentDepth = 50000;
	private boolean secondStrandFlipped = false;
	private boolean isPaired = false;

	//internal fields
	private UCSCGeneLine[] genes;
	private HashMap<String, UCSCGeneLine[]> chromGenes;
	private HashMap<String, UCSCGeneLine> genesWithExons;
	private Condition[] conditions;
	private HashSet<String> flaggedGeneNames = new HashSet<String>();
	private HashSet<String> geneNamesPassingThresholds = new HashSet<String>();
	private HashSet<String> geneNamesWithMinimumCounts = new HashSet<String>();
	private File serializedConditions = null;
	private String url;
	private static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MSDHN])");
	public static final Pattern BAD_NAME = Pattern.compile("(.+)/[12]$");
	public boolean saveCounts = false;

	//paired diff expression
	private UCSCGeneLine[] workingGenesToTest = null;
	private ArrayList<String> workingTNs = null;
	private ArrayList<String> workingConditionNames = null;
	private String workingName = "";
	private File workingCountTable = null;
	private File[] workingDESeqFiles = null;
	private ArrayList<UCSCGeneLine> workingGenesPassingFilters = new ArrayList<UCSCGeneLine>();
	private int numberWorkingGenesPassingFilters =0;
	private StringBuilder spreadSheetHeader = null;
	private int workingNumberFirstReplicas;
	private int workingNumberSecondReplicas;
	private int workingFragmentNameIndexPlus = 1;
	private int workingFragmentNameIndexMinus = -1;
	private HashMap<String, Integer> workingFragNameIndex = new HashMap<String, Integer>(10000); // Maps from the FULL fragment name to the integer index that is used internally. Apparently there is no way to map back?

	//from RNASeq app integration
	private File treatmentBamDirectory;
	private File controlBamDirectory;

	//constructors
	/**Stand alone.*/
	public AWDefinedRegionDifferentialSeq(String[] args){	
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
	public AWDefinedRegionDifferentialSeq(File treatmentBamDirectory, File controlBamDirectory, String genomeVersion, File saveDirectory,  File fullPathToR, File processedRefSeqFile, boolean scoreIntrons, boolean performStrandedAnalysis, 
			boolean verbose, boolean reverse, boolean flipped, int maxAlignDepth){
		this.treatmentBamDirectory = treatmentBamDirectory;
		this.controlBamDirectory = controlBamDirectory;
		this.saveDirectory = saveDirectory;
		this.genomeVersion = genomeVersion;
		url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
		this.fullPathToR = fullPathToR;
		this.refSeqFile = processedRefSeqFile;
		this.scoreIntrons = scoreIntrons;
		this.performStrandedAnalysis = performStrandedAnalysis;
		this.verbose = verbose;
		removeOverlappingRegions = false;
		saveCounts = false;
		this.secondStrandFlipped = flipped;
		this.performReverseStrandedAnalysis = reverse;
		this.maxAlignmentDepth = maxAlignDepth;
		run();
	}

	public void run(){
		//load gene models

		if (this.showFragAssignments) {
			if (verbose) System.out.println("Creating a new fragment assignment file (fragment_assignments.txt) in the output directory...");
			try {
				fragAssignBuffer = new BufferedWriter(new FileWriter(new File(saveDirectory, "fragment_assignments.txt")));
				fragAssignBuffer.write("## Fragment assignments are being output to this tab-delimited file with UNIX-style line-endings.\n");
				fragAssignBuffer.write("## (This is because --showfragments was included in the command line call to DefinedRegionDifferentialSeq)\n");
				fragAssignBuffer.write("## Note that although a fragment may appear in MORE THAN ONE exon (example: a fragmen spanning a splice junction),\n");
				fragAssignBuffer.write("## each unique fragment only counts once toward a gene's total count (even for fragments appearing in 2+ exons).\n");
				fragAssignBuffer.write("## Note that per-gene counts that do not match the table format are also printed below (these lines also start with '##'). You may want to filter those lines out before analyzing this data.\n");
				fragAssignBuffer.write("FRAGMENT_NAME" 
						+ "\t" + "USEQ_FRAGMENT_INTERNAL_ID"
						+ "\t" + "GENE_NAME"
						+ "\t" + "EXON_NUM"
						+ "\t" + "TOTAL_EXONS"
						+ "\t" + "FROM_CONDITION"
						+ "\t" + "FROM_FILENAME"
						+ "\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		if (verbose) System.out.println("Loading regions / gene models...");
		loadGeneModels();

		//make TimeCourseConditions and print
		loadConditions();

		//launch pairwise analysis between conditions
		String varianceFiltered = "without";
		if (filterOutliers) varianceFiltered = "with";
		System.out.println("Running pairwise DESeq analysis "+varianceFiltered+" variance outlier filtering to identify differentially expressed and spliced genes...");
		analyzeForDifferentialExpression();
		System.out.println("\n\t"+geneNamesPassingThresholds.size()+"\tgenes differentially expressed (FDR "+Num.formatNumber(minFDR, 2)+", log2Ratio "+Num.formatNumber(minLog2Ratio, 2)+")\n");

		//too many?
		if (geneNamesPassingThresholds.size() > 5000 && minFDR <=13) if (verbose) System.out.println("WARNING: too many differentially expressed genes?! This may freeze the R clustering. Consider more stringent thresholds (e.g. FDR of 20 or 30)\n");

		//sort by max pvalue
		Arrays.sort(genes, new UCSCGeneLineComparatorMaxPValue());
		//sort by max abs log2 ratio, this is only set for diff expressed genes so non diffs are at the bottom
		Arrays.sort(genes, new UCSCGeneLineComparatorMaxAbsLog2Ratio());

		//launch cluster analysis
		if (verbose) System.out.println("Clustering genes and samples, BETA, ...");
		clusterDifferentiallyExpressedGenes();

		//print final spreadsheet for all genes
		printStatSpreadSheet();

		try {
			if (null != fragAssignBuffer) {
				fragAssignBuffer.close(); // As long as it isn't null, it doesn't matter if it was ever opened or not. It will of course throw a fit if it IS null!
				if (verbose) { System.out.println("Finished writing the fragment assignements."); }
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private void loadConditions(){
		//Any saved conditions
		this.serializedConditions = new File(saveDirectory, "conditions.ser");
		if (this.serializedConditions.exists()){
			conditions = (Condition[])IO.fetchObject(this.serializedConditions);
			if (verbose) System.out.println("\nWARNING: Loading "+conditions.length+" conditions from cached data file, delete "+this.serializedConditions+" if you'd like to recount gene exon/ regions....\n");

			//load min count hash, this is normally done in the scanGene() method
			loadMinimumCountsHash();
		} else {
			if (performStrandedAnalysis) {
				System.out.println("\nCollecting STRANDED counts for each gene exon / region...");
			} else {
				System.out.println("\nCollecting counts for each gene exon / region (not strand-specific)...");
			}

			//from RNASeq app?
			if (treatmentBamDirectory != null) {
				// Looks like this is hard-coed for only TWO conditions (treatment/control)
				conditions = new Condition[2]; // In this case, looks like we hard-code conditions to TREATMENT and CONTROL
				conditions[0] = new Condition(treatmentBamDirectory);
				conditions[1] = new Condition(controlBamDirectory);
			} else {
				conditions = new Condition[conditionDirectories.length];
				for (int i=0; i < conditionDirectories.length; i++) {
					// Looks like we can have as many conditions as we want here... one per condition directory
					conditions[i] = new Condition(conditionDirectories[i]);
				}
			}
			//load em with data
			for (int i=0; i < conditions.length; i++) {				
				for (Replica r: conditions[i].getReplicas()){			
					if (verbose) { System.out.print("\t"+r.getNameNumber()); }
					loadReplica(r);
				}
			}

			//any genes with too many reads that were excluded?
			if (flaggedGeneNames.size() > 0) {
				System.err.println("\nWARNING: The following genes/ regions were excluded from the analysis due to one or more bps exceeding the maximum read coverage of "+maxAlignmentDepth+" . If these are genes/ regions you wish to interrogate, increase the maximum read coverage threshold. " +
						"Realize this will likely require additional memory. Often, these are contaminants (e.g. rRNA) or super abundant transcripts that should be dropped from the analysis.\n"+flaggedGeneNames.toString());
				String[] badGeneNames = Misc.hashSetToStringArray(flaggedGeneNames);

				//remove flagged genes from all replicas
				//for each condition
				for (Condition c: conditions) {
					//for each replica
					for (Replica replica: c.getReplicas()) {
						replica.removeFlaggedGenes(badGeneNames);
					}
				}
			}

			if (verbose) { System.out.println(); }
			//save conditions?
			if (saveCounts) { IO.saveObject(this.serializedConditions, conditions); }

		}

		//print conditions
		if (verbose) {
			System.out.println("Conditions, replicas, and mapped counts:");
			for (Condition c: conditions) { System.out.println(c); }
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
				if (bpNames[i] == null) {
					bpNames[i] = new ArrayList<Integer>();
					//old so check size
				} else if (bpNames[i].size() >= maxAlignmentDepth && (maxAlignmentDepth >= 0)) {
					// The check for 0 here is to allow a 'maxAlignmentDepth' of -1 to suppress checking entirely. If maxAlignmentDepth is -1, then we DO NOT disquality anything ever.
					// uhh... too many I guess?
					badBases[i] = true;
					addIt = false;
					bpNames[i] = null;
					if (nameIndex != null) { workingFragNameIndex.remove(nameIndex.name); }
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
				if ((sam.getFirstOfPairFlag() && sam.getReadNegativeStrandFlag()) || (!sam.getSecondOfPairFlag() && !(sam.getReadNegativeStrandFlag()))) {
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
		return new NameInteger(samReadName, index);
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

		HashSet<String> priorChroms = new HashSet<String>();
		String chrom = null;

		//make first chromosome
		SAMRecord sam = null;
		int chrStartBp = -1;
		ArrayList<Integer>[] bpNames = null;
		boolean[] badBases = null;
		SAMRecordIterator iterator = reader.iterator();
		while (iterator.hasNext()){
			sam = iterator.next();
			if (alignmentFails(sam)) { continue; } //unaligned? too many hits?
			chrom = sam.getReferenceName();
			//if (chrom.equals("chr7") == false) continue;
			if (!chromGenes.containsKey(chrom)) {
				continue;
			} else {
				if (verbose) { System.out.print("."); }
				priorChroms.add(chrom);
				chrStartBp = sam.getUnclippedStart()-1;
				badBases = new boolean[chromLength.get(chrom) - chrStartBp];
				bpNames = new ArrayList[badBases.length];
				break; // EXIT the loop entirely!!
			}
		}
		//any data found?
		if (bpNames == null) {
			return;
		} 

		//reset working fields
		workingFragNameIndex.clear();
		workingFragmentNameIndexPlus = 1;
		workingFragmentNameIndexMinus = -1;

		//load first alignment into bpNames
		loadBlocks(sam, chrStartBp, bpNames, badBases); // note that sam is still used from up above! It's whatever it was on its last iteration!

		//for each record
		while (iterator.hasNext()){
			sam = iterator.next();
			if (alignmentFails(sam)) { continue; } //unaligned? too many hits?

			//same chrom?
			if (sam.getReferenceName().equals(chrom)){
				loadBlocks(sam, chrStartBp, bpNames, badBases);
			} else {
				//different chrom so time to scan
				loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);

				//check that new chrom from SAM is something interrogated by their gene list
				boolean reset = false;
				if (!chromGenes.containsKey(sam.getReferenceName())) {
					while (iterator.hasNext()) { //advance until it does
						sam = iterator.next();
						if (chromGenes.containsKey(sam.getReferenceName())) {
							reset = true;
							break;
						}
					}
				} else {
					reset = true;
				}

				//reset
				if (reset){
					//reset working fields
					workingFragNameIndex.clear();
					workingFragmentNameIndexPlus = 1;
					workingFragmentNameIndexMinus = -1;
					chrom = sam.getReferenceName();
					if (verbose) { System.out.print("."); }
					if (priorChroms.contains(chrom)) { Misc.printErrAndExit("\nError: your sam file isn't sorted by chromosome! Aborting.\n"); }
					priorChroms.add(chrom);
					chrStartBp = sam.getUnclippedStart()-1;
					badBases = new boolean[chromLength.get(chrom) - chrStartBp];
					bpNames = new ArrayList[badBases.length];

					//load
					loadBlocks(sam, chrStartBp, bpNames, badBases);
				}

			}
		}
		//add last
		if (verbose) { System.out.println(); }
		if (chromGenes.containsKey(sam.getReferenceName())) {
			loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases); // Run it AGAIN -- but on different genes this time!
		}
		reader.close();
	}

	private boolean alignmentFails(final SAMRecord sam) { // Returns true if the alignment FAILED
		//aligned?
		if (sam.getReadUnmappedFlag()) { return true; }
		//limit to max matches?
		if (maxRepeats != 0) {
			final Object o = sam.getAttribute("IH");
			if (o != null) {
				final int numRepeats = (Integer)o;
				if (numRepeats > maxRepeats) { return true; } // failure!
			}
		}
		return false;
	}


	private void loadGeneCounts(Replica replica, ArrayList<Integer>[] bpNames, int chrStartBp, String chromosome, boolean[] badBases){
		// Useful comment describing this function, from an email from Tim:
		// "Fragments are assigned to a gene in the method loadGeneCounts.  Earlier in
		// the app, the list of reads that 'cover' each base is stored in a array
		// (bpNames in loadGeneCounts). LoadGeneCounts then queries each base of each
		// exon in the gene and adds the list of read  names to a gene-wide HashSet.
		// HashSets only store unique values, so reads that cover more than one base
		// in the gene are only counted once.
		// David uses an integer index instead of the full fragment name throughout
		// the app.  The map between the two is the HashMap workingFragNameIndex."

		HashMap<Integer, String> indexToNameMap = null; // Translates the internal Integer index for each read back to a string (the string is the read's long / ugly name that you'd see in a FASTA/SAM file)
		if (this.showFragAssignments) { // <-- ONLY run this code when the user wants fragment assignments!
			indexToNameMap = new HashMap<Integer, String>(workingFragNameIndex.size()); // Maps from the internal ID back to the original read name. Only useful if we are going to laboriously print out all the fragment assignments; otherwise it should be ignored!
			for (Map.Entry<String, Integer> entry : workingFragNameIndex.entrySet()) {
				// This might take a little while, so we only run this when the user REALLY wants to know everything about the specific read assignments. Don't just run it normally!
				final String  theReadName = entry.getKey();
				final Integer theReadId   = entry.getValue();
				indexToNameMap.put((Integer)theReadId, (String)theReadName);
			}
		}

		HashMap<String, GeneCount> geneCounts = replica.getGeneCounts();
		HashSet<Integer> allReads = new HashSet<Integer>();
		HashSet<Integer> exonReads = new HashSet<Integer>();
		final int lengthBpNames = (bpNames.length - 1);

		//for each gene in the chromosome 
		UCSCGeneLine[] chrGenes = chromGenes.get(chromosome);

		for (int geneNum = 0; geneNum < chrGenes.length; geneNum++){ // For each gene...
			if (chrGenes[geneNum].isFlagged()) continue; // SKIP if... the gene is flagged

			final String geneName = chrGenes[geneNum].getDisplayNameThenName();
			final boolean plusStrand = chrGenes[geneNum].getStrand().equals("+");					
			allReads.clear();

			//get exons
			final ExonIntron[] exons = chrGenes[geneNum].getExons();
			int[] exonCounts = new int[exons.length];

			//for each exon
			exonLoop:
				for (int exonNum=0; exonNum < exons.length; exonNum++){ // For each exon within a gene...
					exonReads.clear();
					int start = exons[exonNum].getStart() - chrStartBp;
					if (start < 0) { start = 0; }

					int end = exons[exonNum].getEnd() - chrStartBp;
					if (end > lengthBpNames) { end = lengthBpNames; }
					//for each base in the exon, see if there is a read

					for (int y = start; y < end; y++){
						//bad base? if so then flag entire gene
						if (badBases[y]) {
							chrGenes[geneNum].setFlagged(true);
							flaggedGeneNames.add(geneName); // Apparently this is now a BAD GENE
							allReads.clear();
							break exonLoop; // exit the exon loop ENTIRELY! (not just this inner loop)
						}

						if (bpNames[y] != null) {
							if (performStrandedAnalysis){
								for (Integer id : bpNames[y]){
									if (((id > 0) && plusStrand) || ((id < 0) && !plusStrand)) { exonReads.add(id); }
								}
							} else if (performReverseStrandedAnalysis) {
								for (Integer id: bpNames[y]) {
									if (((id > 0) && !plusStrand) || ((id < 0) && plusStrand)) { exonReads.add(id); }
									//System.out.println("Added id " + id + " to gene " + geneName + "... total length is " + exonReads.length());
								}
							} else {
								exonReads.addAll(bpNames[y]); // This is the normal thing to do, unless there is stranded analysis
							}
						}
					}

					if (this.showFragAssignments) {
						// Ok, now let's print out all the fragment assignment stuff, assuming we want that info in the first place.
						// Note that this will typically be a TON of data! Probably several gigabytes worth!
						// Note that a read can actually be part of MULTIPLE exons!
						// But it only counts once toward the final gene count.
						for (Integer id: exonReads) {
							final String theName = indexToNameMap.get(id);
							if (null == theName) { Misc.printErrAndExit("\nProgramming error in the part of DefinedRegionDifferentialSeq.java that implements --showfragments. This key should NEVER be null!\n"); }
							try {
								fragAssignBuffer.write(theName 
										+ "\t" + id // USEQ_FRAGMENT_INTERNAL_ID
										+ "\t" + geneName // GENE_NAME
										+ "\t" + (1+exonNum) // EXON_NUM (note the +1!)
										+ "\t" + exons.length // TOTAL_EXONS
										+ "\t" + replica.getNameNumber() // FROM_CONDITION (e.g., treatment, control)
										+ "\t" + replica.getBamFile() // FROM_FILENAME (e.g., my/filepath/file.bam)
										+ "\n"); // Don't forget the newline!
								//System.out.println("      Added new id " + id + " (\"" + theName + "\") --> \"" + geneName + "\" (Exon " + (1+exonNum) + "/" + exons.length + "), which now has " + exonReads.size() + " fragments.");
							} catch (IOException e) {
								e.printStackTrace();
							}
						}
					}

					exonCounts[exonNum] = exonReads.size();
					allReads.addAll(exonReads);
					//if (ALEX_DEBUG) { System.out.println("   --- end of exon " + (exonNum+1) + " ---"); }
					// <-- When we get to here, we finished an EXON
				}

			final int numCounts = allReads.size();
			if (numCounts > 0) {
				GeneCount tcg = new GeneCount(numCounts, exonCounts);
				geneCounts.put(geneName, tcg);
				replica.setTotalCounts(replica.getTotalCounts() + numCounts);
				if (numCounts >= minimumCounts) { geneNamesWithMinimumCounts.add(geneName); }
			}

			if (this.showFragAssignments) {
				// Ok, now let's print out all the fragment assignment stuff, assuming we want that info in the first place.
				// Note that this will typically be a TON of data! Probably several gigabytes worth!
				try {
					fragAssignBuffer.write("## GENE_DATA_FOR" // Starts with a '##' to make it easier to filter these out. Note that these lines do NOT follow the same format as the others
							+ "\t" + geneName
							+ "\t" + "TOTAL_COUNTS_IS"
							+ "\t" + numCounts
							+ "\t" + "FROM_CONDITION"
							+ "\t" + replica.getNameNumber()
							+ "\t" + "FROM_FILENAME"
							+ "\t" + replica.getBamFile()
							//+ "\t" + "TOTAL_COUNTS_THIS_REPLICATE_IS" // <-- this is not accurate since it hasn't actually finished tallying at this point
							//+ "\t" + replica.getTotalCounts()
							+ "\n");
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

			//if (ALEX_DEBUG) { System.out.println("=-=-=-=-=-=-= End of GENE " + geneName + " =-=-=-=-=-=-=-=-="); }

			// <-- When we get to here, we finished a GENE
		} // end of the "for each gene" loop	
	}

	public void printGeneModelsBed(){
		try {
			if (numberWorkingGenesPassingFilters == 0) return;

			File bed = new File(saveDirectory, workingName+"FDR"+(int)minFDR+"Lg2Rto"+ Num.formatNumber(minLog2Ratio, 1)+".bed.gz");
			Gzipper out = new Gzipper(bed);
			//genome version
			out.println("# genome_version = "+genomeVersion);
			//score names
			out.println("# score = VarCorrLog2Ratio");
			//print
			for (int i=0; i< numberWorkingGenesPassingFilters; i++) {
				UCSCGeneLine gene = workingGenesPassingFilters.get(i);
				String name = gene.getDisplayNameThenName()+"_FDR"+(int)gene.getFdr()+"Lg2Rto"+Num.formatNumber(gene.getLog2Ratio(), 2);
				out.println(gene.toStringBed12Float(gene.getLog2Ratio(), name));
			}
			out.close();

			//convert to useq
			new Text2USeq(bed, genomeVersion, "#FF0000");


		} catch (Exception e){
			e.printStackTrace();
		}
	}

	private void clusterDifferentiallyExpressedGenes(){
		//write counts for all genes to file
		File countTable = new File(saveDirectory, "geneCountTable.txt");
		String[] allGenes = new String[genes.length];
		for (int i=0; i < genes.length; i++){
			allGenes[i] = genes[i].getDisplayNameThenName();
		}
		ArrayList<String> conditionNames = writeCountTable(allGenes, countTable);

		//execute clustering
		File varCorrData = executeDESeqCluster(countTable, conditionNames, false);

		//append count and varCorr data onto genes
		appendCountVarCorrData(countTable, varCorrData, conditionNames);

		//this doesn't appear to work that well?  need test data!
		//write counts for genes with > 20 reads
		//File selectCountTable = new File (saveDirectory, "gene20CountTable.txt");
		//String[] selectGenes = new String[geneNamesWithMinimumCounts.size()];
		//int index = 0;
		//for (String geneName: geneNamesWithMinimumCounts) selectGenes[index++] = geneName;

		//id diff express genes when using whole table
		//executeDESeqDiffExpAll(selectCountTable, writeCountTable(selectGenes, selectCountTable));

	}

	public void appendCountVarCorrData(File countTable, File varCorrData, ArrayList<String> conditionNames){
		//add onto header
		spreadSheetHeader.append("\t"); // <-- required spacer
		spreadSheetHeader.append("Counts_"+Misc.stringArrayListToString(conditionNames, "\tCounts_"));
		spreadSheetHeader.append("\t");
		spreadSheetHeader.append("VarCorCounts_"+Misc.stringArrayListToString(conditionNames, "\tVarCorCounts_"));
		spreadSheetHeader.append("\t");
		spreadSheetHeader.append("FPKM_"+Misc.stringArrayListToString(conditionNames, "\tFPKM_"));

		try { 		//append onto genes
			BufferedReader counts = new BufferedReader( new FileReader(countTable));
			//skip header
			counts.readLine();
			BufferedReader varCounts = new BufferedReader (new FileReader(varCorrData));
			//append, these are tab delimited
			for (int i=0; i< genes.length; i++){
				StringBuilder sb = genes[i].getText();
				//sb.append("\t");
				String line = counts.readLine();
				int index = line.indexOf("\t");
				sb.append(line.substring(index));
				sb.append("\t");
				sb.append(varCounts.readLine());
				sb.append("\t");
				//calculate FPKM values
				for (final Condition c: conditions){
					//for each replica
					for (final Replica r: c.getReplicas()){
						double num = 0;
						if (null != r.getGeneCounts().get(genes[i].getDisplayNameThenName())) {
							num = r.getGeneCounts().get(genes[i].getDisplayNameThenName()).calculateFPKM(r.getTotalCounts(), genes[i].getTotalExonicBasePairs());
						}
						sb.append(num);
						sb.append("\t");
					}
				}
			}

			//close readers
			counts.close();
			varCounts.close();

		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing counts and varCorr counts.\n");
		}
	}

	public File executeDESeqCluster(File countTable, ArrayList<String> conditionNames, boolean useLocalFitType){
		int numDiffExp = geneNamesPassingThresholds.size();
		File varCorrData = new File (saveDirectory, "allGene_DESeqVarCorrData.txt");
		File clusteredDiffExpressGenes = new File (saveDirectory, "clusterPlot"+numDiffExp+"DiffExpGenes.pdf");
		File sampleClustering = new File (saveDirectory, "clusterPlotSamples.pdf");
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("numDiffExpGenes = "+numDiffExp+"\n");
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+countTable.getCanonicalPath()+"', header=TRUE)\n");
			sb.append("rownames(countsTable) = countsTable[,1]\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c('"+ Misc.stringArrayListToString(conditionNames, "','") + "')\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("cds = estimateSizeFactors( cds )\n");
			if (useLocalFitType) sb.append("cds = estimateDispersions( cds, method='blind', sharingMode='fit-only', fitType='local' )\n");
			else sb.append("cds = estimateDispersions( cds, method='blind', sharingMode='fit-only' )\n");
			sb.append("vsd = getVarianceStabilizedData( cds )\n");
			//write out the vsd data
			sb.append("write.table(vsd, file = '"+varCorrData.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			//cluster diff expressed genes?
			if (numDiffExp >= 10){
				sb.append("colors = colorRampPalette(c('white','darkblue'))(100)\n");
				sb.append("pdf('"+clusteredDiffExpressGenes.getCanonicalPath()+"', height=10, width=10)\n");
				sb.append("heatmap( vsd[1:numDiffExpGenes,], col=colors, scale='none')\n");
				sb.append("dev.off()\n");
			}
			//cluster by sample
			sb.append("dists = dist( t( vsd ) )\n");
			sb.append("pdf('"+ sampleClustering.getCanonicalPath() +"', height=10, width=10)\n");
			sb.append("heatmap( as.matrix( dists ), symm=TRUE, scale='none', margins=c(10,10), col = colorRampPalette(c('darkblue','white'))(100), labRow = paste( pData(cds)$condition, pData(cds)$type ) )\n");
			sb.append("dev.off()\n");


			//write script to file
			File scriptFile = new File (saveDirectory, "allGene_RScript.txt");
			File rOut = new File(saveDirectory, "allGene_RScript.txt.Rout");
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

			//check for warnings
			String[] res = IO.loadFile(rOut);

			for (String s: res) {
				if (s.contains("Dispersion fit did not converge") || s.contains("Parametric dispersion fit failed")) {
					System.err.println("\tWarning, DESeq's GLM dispersion fit failed. Relaunching using fitType='local'");
					if (deleteTempFiles == false) {
						// If we are NOT deleting temp files... uhh, then apparently delete these files?
						// Seems backwards, but maybe it's correct!
						rOut.delete();
						scriptFile.delete();
					}

					if (useLocalFitType == false) {
						executeDESeqCluster(countTable, conditionNames, true);
					} else {
						throw new IOException();
					}
				}
			}

			//any problems?
			if (varCorrData.exists() == false) throw new IOException();

			//cleanup
			if (deleteTempFiles) {
				rOut.deleteOnExit();
				countTable.deleteOnExit();
				scriptFile.deleteOnExit();
				varCorrData.deleteOnExit();
			}

		} catch (IOException e) {
			Misc.printErrAndExit("\nError failed to cluster data. Check temp files in save directory for error.\n"+ e.getMessage());
		}
		return varCorrData;
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



	public void estimateDifferencesInReadDistributionsWithPermutationTest(Condition one, Condition two){
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();

		//for each gene with min counts
		//note all of the scores have been zeroed or nulled
		for (UCSCGeneLine gene: workingGenesToTest){
			//does it have 2 or more exons?
			UCSCGeneLine gl = genesWithExons.get(gene.getDisplayNameThenName());
			if (gl == null) {
				continue;
			} else {
				//load it with counts
				loadGeneLineWithExonCounts(gl, one, two);
				//calc permutated pvalue
				double pval = Num.calculatePermutedChiSquarePValue (Num.floatArraysToInt(gl.getTreatmentExonCounts()), Num.floatArraysToInt(gl.getControlExonCounts()));
				gl.setSplicingPValue((float)pval);
				setMaxLog2RatioSplice(gl); //set ratio
			}

		}
	}

	public void estimateDifferencesInReadDistributions(Condition one, Condition two){
		int maxNumberExons = -1;
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();

		//for each gene with min counts
		//note all of the scores have been zeroed or nulled
		for (UCSCGeneLine gene: workingGenesToTest){

			//does it have 2 or more exons?
			UCSCGeneLine gl = genesWithExons.get(gene.getDisplayNameThenName());
			if (gl == null) continue;

			//load it with counts
			loadGeneLineWithExonCounts(gl, one, two);

			//check exon counts and modify if too few 
			if (checkForMinimums(gl)){
				al.add(gl);
				int numEx = gl.getExonCounts()[0].length;
				if (numEx > maxNumberExons) maxNumberExons = numEx;
			}
		}
		//any genes?
		if (al.size() == 0) { return; } // nothing, so uh, return early apparently

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
	}

	/**Looks for minimum number of reads and minimum log2Ratio difference between exon counts.*/
	private boolean checkForMinimums(UCSCGeneLine gene){
		//get total treatment and total control
		float[][] tExonCounts = gene.getTreatmentExonCounts();
		float[][] cExonCounts = gene.getControlExonCounts();
		int numExons = tExonCounts[0].length;

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
		float maxLogRatio = 0;
		int maxLogRatioIndex = 0;
		float[] log2Ratios = new float[numExons];
		for (int i=0; i< numExons; i++){
			//enough reads?
			if ((tCounts[i] + cCounts[i]) < 10) continue;
			//check ratio
			log2Ratios[i] = calculateLog2Ratio(tCounts[i], cCounts[i], scalarTC, scalarCT);
			//if (gene.getName().equals("ENSG00000211896_IGHG1")) if (verbose) System.out.println("XXX "+log2Ratios[i]+" "+tCounts[i]+ " "+cCounts[i]+ " "+scalarTC+" "+scalarCT);
			float logRatio = Math.abs(log2Ratios[i]);
			if (logRatio > maxLogRatio) {
				maxLogRatio = logRatio;
				maxLogRatioIndex = i;
			}
			goodIndexes.add(new Integer(i));
		}

		//check ratio and good exons, need at least two for alt splicing
		int numGoodExons = goodIndexes.size();
		if (numGoodExons < 2 || maxLogRatio < minimumSpliceLog2Ratio) {
			gene.setExonCounts(null);
			return false;
		}
		//set log2Ratio 
		gene.setSplicingLog2Ratio(log2Ratios[maxLogRatioIndex]);
		gene.setSplicingExon(maxLogRatioIndex);
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

	/**Sets max log2Ratio difference between exon counts.*/
	private void setMaxLog2RatioSplice(UCSCGeneLine gene){
		//get total treatment and total control
		float[][] tExonCounts = gene.getTreatmentExonCounts();
		float[][] cExonCounts = gene.getControlExonCounts();
		int numExons = tExonCounts[0].length;

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
		float maxLogRatio = 0;
		int maxLogRatioIndex = 0;
		float[] log2Ratios = new float[numExons];
		for (int i=0; i< numExons; i++){
			//enough reads?
			if ((tCounts[i] + cCounts[i]) < 10) continue;
			//check ratio
			log2Ratios[i] = calculateLog2Ratio(tCounts[i], cCounts[i], scalarTC, scalarCT);
			float logRatio = Math.abs(log2Ratios[i]);
			if (logRatio > maxLogRatio) {
				maxLogRatio = logRatio;
				maxLogRatioIndex = i;
			}
		}
		//set log2Ratio 
		gene.setSplicingLog2Ratio(log2Ratios[maxLogRatioIndex]);
		gene.setSplicingExon(maxLogRatioIndex);

	}

	public void analyzeForDifferentialExpression(){
		for (int i=0; i< this.conditions.length; i++){
			Condition first = this.conditions[i];
			for (int j=i+1; j< this.conditions.length; j++){
				Condition second = this.conditions[j];
				differentialExpress(first, second);
			}
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

	public void differentialExpress(Condition first, Condition second){
		if (verbose) System.out.print("\tComparing "+first.getName()+" vs "+ second.getName());

		workingName = first.getName() +"_"+ second.getName();
		workingNumberFirstReplicas = first.getReplicas().length;
		workingNumberSecondReplicas = second.getReplicas().length;

		//write count table
		if (writePairedCountTable(first, second) == false) return;

		//execute DESeq
		executeDESeqDiffExp(false);
		if (workingDESeqFiles == null) return;

		//parse results and populate gene scores
		if (parseDESeqStatResults() == false) return;
		if (verbose) System.out.println(" : " + numberWorkingGenesPassingFilters);

		//differential splice
		if (usePermutationSpliceTest) estimateDifferencesInReadDistributionsWithPermutationTest(first, second);
		else estimateDifferencesInReadDistributions(first, second);

		//print results spread sheet
		if (spreadSheetHeader == null) {
			saveFirstWorkingGeneModels();
		} else {
			saveWorkingGeneModelsAfterTheFirstOne(); // Save the SECOND and onward genes
		}

		//print bed file
		printGeneModelsBed();

	}

	/**Prints a spread sheet of counts from all of the genes*/
	public void printCountSpreadSheet(){
		File f = new File (saveDirectory, "geneCounts.xls");
		String[] allGenes = new String[genes.length];
		for (int i=0; i< genes.length; i++){
			allGenes[i] = genes[i].getDisplayNameThenName();
		}
		writeCountTable(allGenes, f);
	}

	/**Prints a spread sheet from all of the genes sorted by max log2 ratio for those found differentially expressed.*/
	public void printStatSpreadSheet(){
		try {
			File f = new File (saveDirectory, "geneStats.xls.gz");
			Gzipper spreadSheetOut = new Gzipper (f);
			spreadSheetOut.println(spreadSheetHeader.toString() +"\tGenomeVersion=" +genomeVersion);
			for (UCSCGeneLine gene : genes) spreadSheetOut.println(gene.getText().toString());
			spreadSheetOut.close();
		} catch (Exception e){
			System.err.println("\nProblem printing gene models.");
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void saveFirstWorkingGeneModels(){
		// Note: this only affects the FIRST comparison. If there are multiple comparisons (say, DrugX vs DrugY vs Wildtype), then you will also need to change "saveWorkingGeneModelsAfterTheFirstOne" below.

		//build header line
		spreadSheetHeader = new StringBuilder();
		boolean secondNamePresent = false;
		if (workingGenesToTest[0].getDisplayName() != null && workingGenesToTest[0].getName() != null) {
			spreadSheetHeader.append("#DisplayName\tName\t");
			secondNamePresent = true;
		} else {
			spreadSheetHeader.append("#Name\t");
		}
		final String splicePValString = (usePermutationSpliceTest) ? "\tSplicePermChiPVal_" : "\tSpliceChiPVal_";
		spreadSheetHeader.append("Chr\tStrand\tStart\tStop\tTotalBPs\tPVal_"+workingName+"\tFDR_"+workingName+"\tVarCorLg2Rto_"+workingName+ splicePValString +workingName+"\tSpliceMaxLg2Rto_"+workingName+"\tSpliceMaxExon_"+workingName);

		//for each gene
		for (int i=0; i< genes.length; i++){
			//instantiate new SB
			StringBuilder text = new StringBuilder();

			//name
			String name = genes[i].getDisplayNameThenName();

			//url
			int start = genes[i].getTxStart() - 10000;
			if (start < 0) start = 0;
			int end = genes[i].getTxEnd() + 10000;

			text.append(url);
			text.append(genes[i].getChrom());
			text.append("&start=");
			text.append(start);
			text.append("&end=");
			text.append(end);
			text.append("\",\"");
			text.append(name);
			text.append("\")\t");

			//print second name?
			if (secondNamePresent) {
				text.append(genes[i].getName());
				text.append("\t");
			}

			//status? OK or sectioned or read depth
			if (genes[i].isFlagged()){
				text.append("Too many reads");
			}

			//coordinates
			text.append(genes[i].coordinates());
			text.append("\t");

			text.append(genes[i].getTotalExonicBasePairs()); // The "Total BPs" field
			text.append("\t");

			appendResultStatisticsForOneGene(text, genes[i]); // See below. Adds the statistics for a single gene.

			//set text
			genes[i].setText(text);
		}

	}

	public void saveWorkingGeneModelsAfterTheFirstOne(){
		//add to header line
		spreadSheetHeader.append("\tPVal_"+workingName
				+"\tFDR_"+workingName
				+"\tVarCorLg2Rto_"+workingName
				+"\tSplicePVal_"+workingName
				+"\tSpliceMaxLg2Rto_"+workingName
				+"\tSpliceMaxExon_"+workingName);
		//for each gene
		for (int i=0; i< genes.length; i++){
			StringBuilder text = genes[i].getText();
			text.append("\t");
			appendResultStatisticsForOneGene(text, genes[i]); // See below. Adds the statistics for a single gene.
		}
	}

	private void appendResultStatisticsForOneGene(StringBuilder theText, UCSCGeneLine theGene) {
		//scores pval, fdr, log2rto, splicePVal, spliceLog2
		float regularP = 1, regularFdr = 1, regularSplicingP = 1; // This section on reporting non-log P-values was added by Alex Williams at Gladstone in Apr. 2013.
		if (this.shouldReportPValsInLogSpace) {
			// Nothing to do, the P-values are ALREADY in log space! So we can just leave this alone.
		} else {
			// The conversion formula from USeq's log-converted P-value (the "x" variable here) to a standard one is:
			//   * EXCEL FORMULA:       =POWER(10, x/-10)
			//   * R (or Perl) FORMULA: 10**(x/-10)
			//   * Java formula:        Math.pow(10.0, x/-10.0)
			// This gives you a REGULAR p-value from USeq's converted P-value, where the input x is Useq's converted P-value. If you try that with "x" as 13.01, you should get approximately 0.05 as a result.
			regularP         = (float)Math.pow(10.0, theGene.getpValue() / -10.0);
			regularFdr       = (float)Math.pow(10.0, theGene.getFdr() / -10.0);
			regularSplicingP = (float)Math.pow(10.0, theGene.getSplicingPValue() / -10.0);
		}
		theText.append((this.shouldReportPValsInLogSpace ? theGene.getpValue() : regularP));  // <-- by default, in log10 space, but we can ALSO report it in "regular" space
		theText.append("\t");
		theText.append((this.shouldReportPValsInLogSpace ? theGene.getFdr() : regularFdr));  // <-- by default, in log10 space, but we can ALSO report it in "regular" space
		theText.append("\t");
		theText.append(theGene.getLog2Ratio());  // log 2
		theText.append("\t");
		theText.append((this.shouldReportPValsInLogSpace ? theGene.getSplicingPValue() : regularSplicingP));   // <-- by default, in log10 space, but we can ALSO report it in "regular" space
		theText.append("\t");
		theText.append(theGene.getSplicingLog2Ratio());  // log2
		theText.append("\t");
		theText.append(theGene.getSplicingExon());	
	}

	private boolean parseDESeqStatResults(){
		try {
			BufferedReader inStats = new BufferedReader (new FileReader (workingDESeqFiles[0]));
			BufferedReader inVarCorr = new BufferedReader (new FileReader (workingDESeqFiles[1]));
			numberWorkingGenesPassingFilters = 0;
			workingGenesPassingFilters.clear();
			String line;
			String[] stats;
			String[] varCorr;
			Pattern tab = Pattern.compile("\t");
			float maxPVal = 0;
			float maxAdjPVal = 0;
			int totalReps = workingNumberFirstReplicas + workingNumberSecondReplicas;

			for (int x=0; x< workingGenesToTest.length; x++){

				//parse stats line: pval, padj, (meanVarCorT- meanVarCorC)
				line= inStats.readLine();
				stats = tab.split(line);
				if (stats.length!=3) Misc.printErrAndExit("One of the DESeq stats R results rows is malformed -> "+line);
				float[] scores = new float[stats.length];
				for (int i=0; i< stats.length; i++) {
					if (stats[i].equals("Inf") || stats[i].equals("NA") || stats[i].equals("-Inf")) scores[i] = Float.MIN_VALUE;
					else scores[i] = Float.parseFloat(stats[i]);
				}
				//pval
				if (scores[0]> maxPVal) maxPVal = scores[0];
				//adjPval
				if (scores[1]> maxAdjPVal) maxAdjPVal = scores[1];

				//parse variance corrected counts line and recalculate varCorr diff using pseudo median, otherwise use meanT - meanC from R
				if (workingNumberFirstReplicas > 2 || workingNumberSecondReplicas > 2){
					line= inVarCorr.readLine();
					varCorr = tab.split(line);
					if (varCorr.length!=totalReps) Misc.printErrAndExit("One of the DESeq varCorr R results rows is malformed -> "+line);
					scores[2] = calculateDifference(varCorr);
				}
				//add
				workingGenesToTest[x].setpValue(scores[0]);
				workingGenesToTest[x].setFdr(scores[1]);
				workingGenesToTest[x].setLog2Ratio(scores[2]);
				workingGenesToTest[x].setScores(null);
			}
			inStats.close();
			inVarCorr.close();

			//convert Inf to max values * 1%
			maxPVal = maxPVal *1.01f;
			maxAdjPVal = maxAdjPVal * 1.01f;
			for (int i=0; i< workingGenesToTest.length; i++){
				//check pval
				if (workingGenesToTest[i].getpValue() == Float.MIN_VALUE) workingGenesToTest[i].setpValue(maxPVal);
				if (workingGenesToTest[i].getpValue() > workingGenesToTest[i].getMaxPValue()) workingGenesToTest[i].setMaxPValue(workingGenesToTest[i].getpValue());
				//check adjPVal
				if (workingGenesToTest[i].getFdr() == Float.MIN_VALUE) workingGenesToTest[i].setFdr(maxAdjPVal);
				//does it log2ration difference  thresholds?
				float absLog2Ratio = Math.abs(workingGenesToTest[i].getLog2Ratio());
				if (workingGenesToTest[i].getFdr() >= minFDR && absLog2Ratio >= minLog2Ratio){
					geneNamesPassingThresholds.add(workingGenesToTest[i].getDisplayNameThenName());
					numberWorkingGenesPassingFilters++;
					workingGenesPassingFilters.add(workingGenesToTest[i]);
					if (workingGenesToTest[i].getMaxAbsLog2Ratio() < absLog2Ratio) workingGenesToTest[i].setMaxAbsLog2Ratio(absLog2Ratio);
				}
			}
			//clean up
			if (deleteTempFiles) {
				workingDESeqFiles[0].delete();
				workingDESeqFiles[1].delete();
			}
			return true;
		} catch (Exception e){
			System.err.println("\nProblem parsing DESeq stats results from R.");
			e.printStackTrace();
			System.exit(1);
		}
		return false;
	}

	public float calculateDifference (String[] varCorr){
		//parse values
		float[] t = new float[workingNumberFirstReplicas];
		float[] c = new float[workingNumberSecondReplicas];
		for (int i=0; i< workingNumberFirstReplicas; i++){
			t[i] = Float.parseFloat(varCorr[i]);
		}
		int index =0;
		for (int i=workingNumberFirstReplicas; i < varCorr.length; i++){
			c[index++] = Float.parseFloat(varCorr[i]);
		}
		return (float)(Num.pseudoMedian(t) - Num.pseudoMedian(c));
	}

	public float parseFloat(String f){
		if (f.equals("Inf") || f.equals("NA")) return Float.MIN_VALUE;
		else return Float.parseFloat(f);
	}

	public void executeDESeqDiffExp(boolean useLocalFitType){

		File rResultsStats = new File (saveDirectory, workingName+"_DESeqResults.txt");
		File rResultsData = new File (saveDirectory, workingName+"_DESeqResultsData.txt");
		workingDESeqFiles = null;
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+workingCountTable.getCanonicalPath()+"', header=TRUE)\n");
			sb.append("rownames(countsTable) = countsTable[,1]\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c('"+ Misc.stringArrayListToString(workingTNs, "','") + "')\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("cds = estimateSizeFactors( cds )\n");
			if (workingTNs.size() == 2) {
				if (useLocalFitType) sb.append("cds = estimateDispersions( cds, method='blind', sharingMode='fit-only', fitType='local' )\n");
				else sb.append("cds = estimateDispersions( cds, method='blind', sharingMode='fit-only' )\n");
			}

			else {
				String outlier = "";
				String local = "";
				if (filterOutliers == false) outlier= ", sharingMode='fit-only'";
				if (useLocalFitType) local = ", fitType='local'";
				sb.append("cds = estimateDispersions( cds" +outlier + local +")\n");
			}
			sb.append("res = nbinomTest( cds, 'N', 'T', pvals_only = FALSE)\n");
			//Recalculate log2 ratio using moderated values
			sb.append("vsd = getVarianceStabilizedData( cds )\n");
			sb.append("res[,6] = (rowMeans( vsd[, conditions(cds)=='T', drop=FALSE] ) - rowMeans( vsd[, conditions(cds)=='N', drop=FALSE] ))\n");
			//Fred  pvalues
			sb.append("res[,7] = -10 * log10(res[,7])\n");
			sb.append("res[,8] = -10 * log10(res[,8])\n");
			//Parse pval, padj, log2ratio; note flip of A and B back to T and C
			sb.append("res = res[,c(7,8,6)]\n");
			//note, the order of the rows is the same as the input
			sb.append("write.table(res, file = '"+rResultsStats.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			sb.append("write.table(vsd, file = '"+rResultsData.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");

			//write script to file
			File scriptFile = new File (saveDirectory, workingName +"_RScript.txt");
			File rOut = new File(saveDirectory, workingName +"_RScript.txt.Rout");
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


			//check for warnings
			String[] res = IO.loadFile(rOut);

			for (String s: res) {
				if (s.contains("Dispersion fit did not converge") || s.contains("Parametric dispersion fit failed")) {
					System.err.println("\n\t\tWarning, DESeq's GLM dispersion fit failed. Relaunching using fitType='local'");
					if (deleteTempFiles == false) {
						rOut.delete();
						scriptFile.delete();
					}
					executeDESeqDiffExp(true);
					return;
				}
			}

			//look for results file
			if (rResultsStats.exists() == false ) throw new IOException("\n\nR results file doesn't exist. Check temp files in save directory for error.\n");


			//cleanup
			if (deleteTempFiles) {
				rOut.deleteOnExit();
				scriptFile.deleteOnExit();
			}

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error: failed to execute DESeq.\n");
		}
		workingDESeqFiles = new File[]{rResultsStats,rResultsData};
	}

	public void setTNs(final Condition first, final Condition second){
		this.workingTNs = new ArrayList<String>(); // Clear it out!
		this.workingConditionNames = new ArrayList<String>(); // Clear it out!
		for (final Replica r: first.getReplicas()) {
			this.workingTNs.add("T");
			this.workingConditionNames.add(r.getNameNumber());
		}
		for (final Replica r: second.getReplicas()) {
			this.workingTNs.add("N");
			this.workingConditionNames.add(r.getNameNumber());
		}
	}

	public boolean writePairedCountTable(Condition first, Condition second){
		try {
			HashSet<String> genesToTest = new HashSet<String>();
			Condition[] twoConditions = new Condition[]{first, second};

			//write matrix of name, t1,t2,t3...c1,c2,c3 to file for genes with observations
			workingCountTable = new File(saveDirectory, "countTable_"+workingName+".txt");
			if (deleteTempFiles) { workingCountTable.deleteOnExit(); }
			PrintWriter out = new PrintWriter( new FileWriter(workingCountTable));

			//print header also look for genes with sufficient reads
			out.print("GeneName");
			for (Condition c: twoConditions){
				for (Replica r: c.getReplicas()){
					out.print("\t");
					out.print(r.getNameNumber());

					//scan for genes that pass minimum read coverage
					HashMap<String, GeneCount> geneCounts = r.getGeneCounts();
					Iterator<String> it = geneCounts.keySet().iterator();
					while (it.hasNext()){
						final String geneName = it.next();
						if (flaggedGeneNames.contains(geneName)) {
							continue;
						} else {
							final GeneCount gc = geneCounts.get(geneName);
							if (gc != null && gc.getCount() >= minimumCounts) {
								genesToTest.add(geneName);
								if (debugVerbose) { System.err.println("DEBUG: The gene " + geneName + " had enough counts (" + gc.getCount() + ")"); }
							} else {
								if (debugVerbose) { System.err.println("DEBUG: The gene " + geneName + " did not have the minimum number of required counts (" + minimumCounts + ")"); }
							}
						}
					}
				}

			}
			out.println();

			//any genes to test?
			if (genesToTest.size() == 0) {
				out.close();
				System.err.println("\nWARNING: no genes were found with minimum counts.  Skipping "+workingName);
				return false;
			}

			setTNs(first, second);

			//fetch ucsc gene lines
			workingGenesToTest = new UCSCGeneLine[genesToTest.size()];
			int index = 0;
			for (UCSCGeneLine gene: genes){
				//zero scores
				gene.zeroNullScores();
				if (genesToTest.contains(gene.getDisplayNameThenName())) {
					workingGenesToTest[index] = gene;
					index++;
				}
			}

			//fetch counts, print and set in UCSCGeneLine
			for (UCSCGeneLine gene: workingGenesToTest){
				String geneName = gene.getDisplayNameThenName();
				out.print(geneName);
				for (Condition c: twoConditions){
					//for each replica
					for (Replica r: c.getReplicas()){
						out.print("\t");
						int count = 0;
						if (r.getGeneCounts().get(geneName) != null) count = r.getGeneCounts().get(geneName).getCount();
						out.print(count);
					}
				}
				out.println();
			}
			out.close();

			return true;
		}
		catch (Exception e){
			System.err.println("Problem writing out count table.");
			e.printStackTrace();
		}
		return false;
	}

	public ArrayList<String> writeCountTable(String[] genesNamesToWrite, File countTable){
		ArrayList<String> conditionNames = new ArrayList<String>();
		try {
			//write matrix of name, t1,t2,t3...c1,c2,c3 to file for genes with observations
			PrintWriter out = new PrintWriter( new FileWriter(countTable));

			//print header
			out.print("GeneName");
			//for each condition
			for (Condition c: conditions){
				//for each replica
				for (Replica r: c.getReplicas()){
					out.print("\t");
					out.print(r.getNameNumber());
					conditionNames.add(r.getNameNumber());
				}
			}
			out.println();

			//print counts
			for (String geneName: genesNamesToWrite){
				out.print(geneName);
				for (Condition c: conditions){
					//for each replica
					for (Replica r: c.getReplicas()){
						out.print("\t");
						int num = 0;
						if (r.getGeneCounts().get(geneName) != null) num = r.getGeneCounts().get(geneName).getCount();
						out.print(num);
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
		return conditionNames;
	}

	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = null;
		if (refSeqFile != null){
			reader = new UCSCGeneModelTableReader(refSeqFile, 0);
			if (scoreIntrons){
				if (verbose) System.out.println("\tParsing/scoring introns instead of exons from gene models.");
				reader.swapIntronsForExons();
			}
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
				if (verbose) System.out.println("\tWrote the trimmed gene table to the save directory. Use this table and the -o option to speed up subsequent processing.");
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

		genesWithExons = new HashMap<String, UCSCGeneLine>();
		for (UCSCGeneLine gl: genes) {
			if (gl.getExons().length > 1) genesWithExons.put(gl.getDisplayNameThenName(), gl);
		}
	}


	/**Loades the geneNamesWithMinimumCounts hash with gene names that pass the minimumCounts array.
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
		new AWDefinedRegionDifferentialSeq(args);
	}		

	/* This method will process each command line argument */
	public void processArgs(String[] args){
		// Test for running this program: java -jar ./Releases/USeq_vcf/Apps/DefinedRegionDifferentialSeq --plainp

		// This is some kind of hand-crafted version of GetOpts / GetOptions.
		// There are a lot of unusual programming specifics that make it work, so beware if you change any of it!
		// In particular, note the "args[++i]" code---this is our way of getting the next argument! Be careful!
		// Java has no built-in "getopts", the closest thing is probably: http://commons.apache.org/proper/commons-cli/

		// BEWARE: appears to treat -a and -A the same way. I guess this limits you to 26 options.

		final Pattern pat = Pattern.compile("-[a-z]");
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			final String lcArg = args[i].toLowerCase(); // lower-case-ify the argument
			final Matcher mat = pat.matcher(lcArg);

			// Before we check any of the single-character options, we check to see if any of these
			// LONG command-line options were specified.

			if (lcArg.equals("--help") || lcArg.equals("-help") || lcArg.equals("--man") || lcArg.equals("-man")) {
				printDocs(); System.exit(0);
			} else if (lcArg.equals("--showfragments") || lcArg.equals("-showfragments")) {
				this.showFragAssignments = true;
				if (verbose) { System.out.println("Note that due to the --showfragments command line option, all fragment assignments (i.e., which fragments are assigned to which genes/exons) are being printed out. This is VERY verbose and is usually only useful for debugging! This typically generates several gigabytes of data in the output file \"fragment_assignments.txt\".\n"); }
				//Misc.printErrAndExit("\nProblem, Alex added an option for printing SAM output (like htseq has) but did NOT add any support for it!!!!");
			} else if (lcArg.equals("--plainp") || lcArg.equals("-plainp")) {
				this.shouldReportPValsInLogSpace = false; // We actually want to report p-values in REGULAR (not log-transformed) space.
				if (verbose) System.out.println("Note that due to the --plainp command line option, all P-values are being reported as REGULAR numbers (on a 0-to-1 scale) and not the log-transformed space that is the USeq default.\n");
			} else if (lcArg.equals("--debug")) {
				this.debugVerbose = true;
				this.verbose = true; // ALWAYS verbose when debugging
				System.out.println("DEBUGGING VERBOSE MODE ON: Due to the --debug option being specified, we will print a huge quantity of data to the terminal, for debugging purposes.");
			} else if (mat.matches()){
				final char test = args[i].charAt(1);

				if (((String)args[i]).length() != 2) {
					// If we get here, then everything SHOULD be a two-character option like '-a' or '-f' or whatever.
					// If it is NOT two characters exactly, then this is a huge problem and we should report that the command line arguments are invalid!
					Misc.printErrAndExit("ERROR: Exiting. We encountered an unknown command line option, specifically: " + args[i] + " . You should probably check to make sure everything is in order in your command line call---maybe some option that requires a value has a blank value?");
				}

				try{
					switch (test){
					case 'a': usePermutationSpliceTest = true; break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'c': conditionDirectories = IO.extractOnlyDirectories(new File (args[++i])); break;
					case 'd': saveCounts = true; break; // <-- D is a HIDDEN OPTION!!!!!
					case 'e': minimumCounts = Integer.parseInt(args[++i]); break;
					case 'f': minFDR = Float.parseFloat(args[++i]); break;
					case 'g': genomeVersion = args[++i]; break;
					case 'h': printDocs(); System.exit(0); break;
					case 'i': scoreIntrons = true; break;
					case 'j': performReverseStrandedAnalysis = true; break;
					case 'k': secondStrandFlipped = false; break;
					case 'l': minLog2Ratio = Float.parseFloat(args[++i]); break;
					case 'm': removeOverlappingRegions = true; break;
					case 'n': maxRepeats = Integer.parseInt(args[++i]); break;
					// o is unused
					case 'p': performStrandedAnalysis = true; break;
					// q is unused
					case 'r': fullPathToR = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 't': deleteTempFiles = false; break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'v': filterOutliers = true; break;
					// w is unused
					case 'x': {
						if (args[i+1].toLowerCase().equals("all")) { // <-- do NOT increment i here!
							maxAlignmentDepth = -1; // Any negative value is the same as "allow any depth"
						} else {
							maxAlignmentDepth = Integer.parseInt(args[i+1]); // do NOT increment i here!
							if (maxAlignmentDepth < 1) {
								Misc.printErrAndExit("\nProblem! You specified a number LESS THAN 1 for maximum alignment depth (the -x option)! That will disqualify EVERY SINGLE READ, so you will not get any results. You need to change your '-x' option to a valid value. See the documentation.");
							}
						}
						i++; // <-- Skip the next argument, since it was a value for -x. This is CRITICAL!
						break; }
					// y is unused
					// z is unused
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			} else {
				// Apparently we don't do anything if we get here / should never get here
			}
		}

		//fetch genome version and make url
		if (genomeVersion == null) { Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n"); }
		url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";

		//look for bam files
		if (conditionDirectories == null || conditionDirectories.length == 0) Misc.printErrAndExit("\nError: cannot find any condition directories?\n");
		if (conditionDirectories.length < 2) Misc.printErrAndExit("\nError: must provide at least two Conditions for analysis.\n");
		for (File dir: conditionDirectories){
			File[] bamFiles = IO.extractFiles(dir, ".bam"); // Note from Alex Williams: interesting: this appears to actually also accept SAM files, as long as they are named '.bam' for some reason.
			OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, true);
		}

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdir();

		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else {
			String errors = IO.runRCommandLookForError("library(DESeq)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printErrAndExit("\nError: Cannot find the required R library.  Did you install DESeq " +
						"(http://www-huber.embl.de/users/anders/DESeq/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(DESeq)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
			}
		}

		//Make sure stranded and reverse aren't both set
		if (performStrandedAnalysis == true && performReverseStrandedAnalysis == true) {
			Misc.printErrAndExit("\nError: Stranded and reverse-stranded analyses both set to true, pick one and restart\n\n");
		}

		//look for estimateDispersions() function
		if (estimateDispersions(fullPathToR, saveDirectory) == false){
			Misc.printErrAndExit("\nError: Please upgrade DESeq to the latest version, see http://www-huber.embl.de/users/anders/DESeq/ \n");
		}

		//look for bed file
		if (refSeqFile == null && bedFile == null){
			Misc.printErrAndExit("\nPlease enter a regions or gene file to use in scoring regions.\n");
		}
	}	

	/**Returns true if DESeq uses the estimateDispersions function otherwise false.*/
	public static boolean estimateDispersions(File fullPathToR, File tempDir){
		//this will return an error
		String error = IO.runRCommandLookForError("library(DESeq); estimateDispersions(5);", fullPathToR, tempDir);
		if (error.contains("could not find function")) {
			System.err.println("\nWARNING: You have installed an obsolete version of DESeq. Update R, Bioconductor, and DESeq to latest versions. Key changes have been implemented to control for variance outliers. Don't use this old version!\n");
			return false;
		}
		return true;
	}

	public static void printDocs(){
		System.out.println("\n"
				+ "**************************************************************************************\n"
				+ "**                        Defined Region Differential Seq: Apr 2013                 **\n"
				+ "**************************************************************************************\n"
				+ "DRDS takes bam files, one per replica, minimum one per condition, minimum two\n"
				+ "conditions (e.g. treatment and control or a time course/ multiple conditions) and\n"
				+ "identifies differentially expressed genes under any pairwise comparison using DESeq.\n"
				+ "DESeq's variance corrected count data is used to heirachically cluster the \n"
				+ "differentially expressed genes as well as the samples. See the DESeq manual for\n"
				+ "details. Alternative splicing is estimated using a chi-square test of independence.\n"
				+ "In addition to the cluster plots, a spread sheet is created with the pValue,\n"
				+ "FDR, and variance corrected log2Ratios for each of the pairwise comparisons as well as\n"
				+ "the raw, FPKM, and log2 variance corrected alignment counts.  Use the later for\n"
				+ "subsequent clustering and distance estimations.\n"
				+ "\n"
				+ "Options:\n"
				+ "-s PATH : Output save directory (mandatory). USeq will write to this directory.\n"
				+ "-c PATH : Input conditions directory (mandatory).\n"
				+ "       * This directory must contain a single sub-directory for each condition\n"
				+ "       * And one BAM file per biological replicate (plus .BAI index files)\n"
				+ "       * 3 or more replicates are recommended per condition.\n"
				+ "       * The input BAM files MUST be sorted by COORDINATE using Picard's SortSam.\n"
				+ "       * Splice junction coordinates must be converted to genomic coordinates (see\n"
				+ "         USeq's SamTranscriptomeParser).\n"
				+ "       * Example directory structure for the sample command '-c MY_EXPERIMENT':\n"
				+ "           MY_EXPERIMENT/\n"
				+ "                         CONTROL/\n"
				+ "                                 control1.bam\n"
				+ "                                 control2.bam\n"
				+ "                                 control3.bam\n"
				+ "                          DRUG_X/\n"
				+ "                                 drugX_replicate1.bam\n"
				+ "                                 drugX_replicate2.bam\n"
				+ "                                 drugX_replicate3.bam\n"
				+ "                          DRUG_Y/\n"
				+ "                                 drugY_replicate1.bam\n"
				+ "                                 drugY_replicate2.bam\n"
				+ "                                 drugY_replicate3.bam\n"
				+ "                                 drugY_replicate4.bam\n"
				+ "       * Note that you can also use symbolic links to refer to your BAM files.\n"
				+ "\n"
				+ "-r FILEPATH : Full path to R loaded with DESeq library. Default: '/usr/bin/R'\n"
				+ "       See http://www-huber.embl.de/users/anders/DESeq/ . Type 'library(DESeq)' in\n"
				+ "       an R terminal to see if DESeq is already installed. (If not, you will need to install\n"
				+ "       DESeq through Bioconductor. Search online for 'bioconductor' and 'DESeq' for details.)\n"
				+ "\n"
				+ "ANNOTATION FILE (there are 2 options, you must include exactly one):\n"
				+ "-u FILEPATH : UCSC RefFlat or RefSeq gene table file, full path. Tab delimited, see RefSeq Genes\n"
				+ "       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n"
				+ "       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n"
				+ "       (commaDelimited)exonEnds).\n"
				+ "       Example: ENSG00000183888 C1orf64 chr1 16203317 16207889 16203385 16205428 2 16203317,16205000 16203467,16207889\n"
				+ "       NOTE: this table should contain only ONE composite transcript per gene (e.g.\n"
				+ "       use Ensembl genes, NOT transcripts).\n"
				+ "       Use the MergeUCSCGeneTable app to collapse transcripts.\n"
				+ "       See http://useq.sourceforge.net/usageRNASeq.html for details.\n"
				+ " OR instead of a UCSC-formatted file, you can use a BED file:\n"
				+ "-b FILEPATH : A bed file (chr, start, stop,...), full path, See\n"
				+ "       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"
				+ "\n"
				+ "-g VERSION : Genome Version (ie H_sapiens_Mar_2006). Important for generating clickable\n"
				+ "             links in the USeq output. See http://genome.ucsc.edu/FAQ/FAQreleases.\n"
				+ "             If you are using USeq on a species that is not supported by the UCSC browser,\n"
				+ "             then you can use any text. Example: -g my_unmapped_species_build_19.\n"
				+ "\n"
				+ "--plainp : Instead of reporting P-values in log space, report them as standard P-values.\n"
				+ "           (By default, USeq reports -10Log10(P-value) instead of the P-value.)\n"
				+ "           Example: without this option, the value 13.01 might be reported.\n"
				+ "                    whereas WITH this option, the value 0.05 would be reported.\n"
				+ "           Note: with --plainp, LOWER numbers are more significant (1.0 is worst).\n"
				+ "           Without it, HIGHER numbers are more significant (0.0 is worst).\n"
				+ "\n"
				+ "Advanced Options:\n"
				+ "-v : Filter for variance outliers in DESeq. Default: no filtering.\n"
				+ "-m : Mask overlapping gene annotations. Recommended for well annotated genomes.\n"
				+ "-x NUMBER : Max per base alignment depth, defaults to 50000. Genes containing such high\n"
				+ "            density coverage are ignored. Warnings are thrown.\n"
				+ "            If you really want to allow ANY number of reads, set '-x ALL'.\n"
				+ "            Negative numbers and 0 are not allowed, as that would disqualify all genes.\n"
				+ "-n INTEGER : Max number repeat alignments. Defaults to all.  Assumes 'IH' tags have been set by\n"
				+ "             processing raw alignments with the SamTranscriptomeProcessor.\n"
				+ "-f NUMBER : Minimum FDR threshold for sorting. Default is 10 (-10Log10(FDR=0.1)).\n"
				+ "            Example values: 10 means FDR<=0.1, 13.01 means FDR<=0.05, 20 means FDR<=0.01.\n"
				+ "            Set it to 0 to not filter ANY values.\n"
				+ "-l NUMBER : Minimum absolute varCorLog2Rto threshold for sorting. Defaults is 1.0 (2x fold change).\n"
				+ "            Example values: 0 means 'no filtering at all'. 1.0 means '2x change',\n"
				+ "            2.0 means '4x change', 3.0 means '8x change', etc.\n"
				+ "-e INTEGER : (Default: 20 reads). Minimum number alignments required for a gene/region to be displayed.\n"
				+ "-i : Score introns instead of exons.\n"
				+ "-p : Perform a stranded analysis. Only collect reads from the same strand as the\n"
				+ "      annotation. By default, the analysis is UNstranded.\n"
				+ "-j : Reverse stranded analysis.  Only collect reads from the opposite strand of the\n"
				+ "      annotation.  This setting should be used for the Illumina's strand-specific dUTP protocol.\n"
				+ "-k : Second read flipped. This setting can be used to flip the strand of the second read in a pair.\n"
				+ "      This setting makes it easier to view in IGB, but can break other downstream applications.\n"
				+ "-a : Perform a permutation based chi-square test for differential exon usage.\n"
				+ "      Needs 4 or more replicas per condition.\n"
				+ "\n"
				+ "DEBUGGING OPTIONS\n"
				+ "-t : Don't delete temp files (R script, R results, Rout, etc..). By default, they are deleted..\n"
				+ "--showfragments : Show the assignments of each fragment -> gene. May be useful for debugging.\n"
				+ "                  Assignments are written to the file \"fragment_assignments.txt\"\n"
				+ "                  (which goes into the save directory specified with the '-s' option).\n"
				+ "                  Note that this file may be many gigabytes in size!\n"
				+ "--debug : Turns on an extra-verbose mode that prints out diagnostic information.\n"
				+ "\n"
				+ "USAGE EXAMPLES\n"
				+ "Example command to run with 4 GB of memory in java (-Xmx4G):\n"
				+ "   java -Xmx4G -jar pathTo/USeq/Apps/DefinedRegionDifferentialSeq\n"
				+ "               -c /Data/TimeCourse/ESCells  -s /Data/TimeCourse/My_ES_Save_Directory\n"
				+ "               -g H_sapiens_Feb_2009  -u /Anno/mergedHg19EnsemblGenes.ucsc.gz\n"
				+ "   (That command should be written all on one line, with no line breaks.)\n"
				+ "\n"
				+ "Example command supposing you want even genes with no coverage or differential expression:\n"
				+ "   java -Xmx4G -jar useq/path/DefinedRegionDifferentialSeq\n"
				+ "               -c /my/condition/input  -s /my/save/directory\n"
				+ "               -g New_species_version_12  -u new_species_annotation.ucsc\n"
				+ "               -l 0  -f 0  -e 0  -x ALL\n"
				+ "   (Warning: the command above will take a long time to run AND generate a lot of superfluous data.)\n"
				+ "   (That command should also be written all on one line, with no line breaks.)\n"
				+ "\n"
				+ "**************************************************************************************\n");

	}
}
