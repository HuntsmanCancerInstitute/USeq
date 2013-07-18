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
public class DefinedRegionDifferentialSeq {

	//user defined fields
	private File[] conditionDirectories;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private File bedFile;
	private File refSeqFile;
	private String genomeVersion;
	private float minFDR = 13;
	private float minLog2Ratio = 1f;
	private boolean scoreIntrons = false;
	private boolean removeOverlappingRegions = false;
	private int minimumCounts = 20;
	private boolean deleteTempFiles = true;
	private boolean filterOutliers = false;
	private float minimumSpliceLog2Ratio = 1f;
	private boolean performStrandedAnalysis = false;
	private boolean performReverseStrandedAnalysis = false;
	private boolean usePermutationSpliceTest = false;
	private int maxRepeats = 0;
	private boolean verbose = true;
	private int maxAlignmentsDepth = 50000;
	private boolean secondStrandFlipped = false;

	//internal fields
	private UCSCGeneLine[] genes;
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private HashMap<String, UCSCGeneLine> genesWithExons;
	private HashMap<String, UCSCGeneLine> name2Gene;
	private Condition[] conditions;
	private HashSet<String> flaggedGeneNames = new HashSet<String>();
	private HashSet<String> geneNamesPassingThresholds = new HashSet<String>();
	private HashSet<String> geneNamesWithMinimumCounts = new HashSet<String>();
	private ArrayList<String> geneNamesPassingThresholdsEdgeR;
	private File serializedConditons = null;
	private String url;
	private static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MSDHN])");
	public static final Pattern BAD_NAME = Pattern.compile("(.+)/[12]$");
	public boolean saveCounts = true;

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
	private HashMap<String, Integer> workingFragNameIndex = new HashMap<String, Integer>(10000);

	//from RNASeq app integration
	private File treatmentBamDirectory;
	private File controlBamDirectory;


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
		this.maxAlignmentsDepth = maxAlignDepth;
		run();
	}

	public void run(){
		//load gene models
		if (verbose) System.out.println("Loading regions/ gene models...");
		loadGeneModels();

		//make TimeCourseConditions and print
		loadConditions();

		//launch pairwise analysis between conditions
		String varianceFiltered = "without";
		if (filterOutliers) varianceFiltered = "with";
		System.out.println("Running pairwise DESeq analysis "+varianceFiltered+" variance outlier filtering to identify differentially expressed and spliced genes...");
		analyzeForDifferentialExpression();
		System.out.println("\n\t"+geneNamesPassingThresholds.size()+"\tgenes differentially expressed (FDR "+Num.formatNumber(minFDR, 2)+", log2Ratio "+Num.formatNumber(minLog2Ratio, 2)+")\n");
		String thresholdName = Num.formatNumber(minFDR,1)+"FDR"+Num.formatNumber(minLog2Ratio, 1)+"Lg2Rto";
		File f = new File (saveDirectory, "diffExprGenes_DESeqAllPair_"+thresholdName+".txt");
		IO.writeHashSet(geneNamesPassingThresholds, f);


		//launch all pair ANOVA-like edgeR analysis
		if (conditions.length > 2){
			System.out.println("Launching ANOVA-like edgeR any condition analysis...");
			if (runEdgeRAllConditionAnalysis()){
				System.out.println("\t"+geneNamesPassingThresholdsEdgeR.size()+"\tGenes/ regions identified as differentially expressed.\n");
				File e = new File (saveDirectory, "diffExprGenes_EdgeRANOVALike_"+thresholdName+".txt");
				IO.writeArrayList(geneNamesPassingThresholdsEdgeR, e);
			}
			else System.err.println("Skipped edgeR analysis!\n");
		}

		//sort by max pvalue
		Arrays.sort(genes, new UCSCGeneLineComparatorMaxPValue());

		//launch cluster analysis
		if (verbose) System.out.println("Clustering genes and samples...");
		clusterDifferentiallyExpressedGenes();

		//print final spreadsheet for all genes
		printStatSpreadSheet();


	}


	private boolean runEdgeRAllConditionAnalysis() {
		//find genes with at least 20 reads across all conditions and haven't been flagged
		ArrayList<String> genesToTest = new ArrayList<String>();
		Iterator<String> it = geneNamesWithMinimumCounts.iterator();
		while (it.hasNext()){
			String geneName = it.next();
			if (flaggedGeneNames.contains(geneName) == false) genesToTest.add(geneName);
		}

		//write it out
		workingCountTable = new File(saveDirectory, "countTable_EdgeRANOVA.txt");
		writeCountTable(Misc.stringArrayListToStringArray(genesToTest), workingCountTable);

		//create and execute script
		String[] results = launchEdgeRANOVA(genesToTest.size());

		if (results != null){
			//parse results and append FDRs to genes
			parseEdgeRANOVAResults(results);

			//append result to gene txt
			appendGeneTxtWithEdgeRANOVA();
		}
		else return false;

		return true;
	}

	private void appendGeneTxtWithEdgeRANOVA() {
		geneNamesPassingThresholdsEdgeR = new ArrayList<String>();
		//append header
		spreadSheetHeader.append("\tEdgeR_ANOVALike_FDR\tMaxAbsVarCorLog2Rto");
		//for each gene
		for (int i=0; i< genes.length; i++){
			StringBuilder txt = genes[i].getText();
			txt.append("\t");
			txt.append(genes[i].getFdrEdgeR());
			txt.append("\t");
			txt.append(genes[i].getMaxAbsLog2Ratio());
			if (genes[i].getFdrEdgeR() >= minFDR && genes[i].getMaxAbsLog2Ratio() >= minLog2Ratio) geneNamesPassingThresholdsEdgeR.add(genes[i].getDisplayNameThenName());
		}
	}

	/*Parses edgeR ANOVA like results table
	 * genes	logFC.GV	logFC.ICM	logFC.MI	logFC.MII	logFC.MOR	logFC.PN	logFC.TROPH	logCPM	LR	PValue	FDR
	13350	ENSG00000179046_TRIML2	-1.328269e+01	 9.853460e-01	-1.328269e+01	-1.328269e+01	 7.211705e-01	-6.308912e+00	 5.338231e-01	 5.0623416637	4.832502e+03	 0.000000e+00	 0.000000e+00
	5808	ENSG00000129824_RPS4Y1	-1.258739e+01	 2.077325e+00	-1.258739e+01	-1.258739e+01	 1.975965e+00	-5.544523e+00	 2.389599e+00	 5.6126126090	5.412190e+03	 0.000000e+00	 0.000000e+00
	4908	ENSG00000121570_DPPA4	-1.252600e+01	 2.049318e+00	-1.056383e+01	-1.252600e+01	 1.466659e+00	-5.974433e+00	 1.463207e+00	 5.1144759329	7.345167e+03	 0.000000e+00	 0.000000e+00
	 */
	private void parseEdgeRANOVAResults(String[] results) {
		Pattern tab = Pattern.compile("\\t");
		//find max fdr, not transformed so smallest
		String[] geneNames = new String[results.length];
		double[] fdrs = new double[results.length];

		//parse first one skipping header line
		String[] tokens = tab.split(results[1]);
		int numCol = tokens.length;
		if (numCol != (conditions.length + 5)) Misc.printErrAndExit("\nError: too few columns found in edgeR ANOVA- like output?! Aborting, see -> "+results[1]+"\n");
		int fdrColIndex = numCol - 1; 
		geneNames[1] = tokens[1];
		fdrs[1] = Double.parseDouble(tokens[fdrColIndex]);
		double maxFDR = fdrs[1];

		//parse remainder
		for (int i=2; i< results.length; i++){
			tokens = tab.split(results[i]);
			if (tokens.length != numCol) Misc.printErrAndExit("\nError: different number of columns in edgeR ANOVA-like results line?! Aborting, see -> "+results[i]+"\n");
			geneNames[i] = tokens[1];
			fdrs[i] = Double.parseDouble(tokens[fdrColIndex]);
			if (maxFDR == 0) maxFDR = fdrs[i];
			else if (fdrs[i] !=0 && fdrs[i] < maxFDR) maxFDR = fdrs[i];
		}
		//transform maxFDR
		float transMaxFDR = (float)(Num.minus10log10(maxFDR) * 1.01);

		//modify genes edgeR fdr, watch for NaN or infinity?
		for (int i=1; i< geneNames.length; i++){
			UCSCGeneLine gene = name2Gene.get(geneNames[i]);
			if (fdrs[i] == 0) gene.setFdrEdgeR(transMaxFDR);
			else gene.setFdrEdgeR((float)(Num.minus10log10(fdrs[i])));
		}
	}

	private String[] launchEdgeRANOVA(int numGenesToTest) {
		//create & write out edgeR script
		File edgeRResults = new File(saveDirectory, "results_EdgeRANOVA.txt");
		String script = createEdgeRANOVAScript(edgeRResults);
		File edgeRScript = new File(saveDirectory, "script_EdgeRANOVA.txt");
		IO.writeString(script, edgeRScript);
		File rOut = new File(saveDirectory, "script_EdgeRANOVA.txt.Rout");

		//make command
		String[] command = new String[] {
				fullPathToR.toString(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				edgeRScript.toString(),
				rOut.toString()};

		//execute command
		IO.executeCommandLine(command);

		//check for errors and warnings in the rOut file
		String[] res = IO.loadFile(rOut);
		for (String s: res) {
			if (s.contains("Error")){
				System.err.println("\nError thrown by edgeR application:\n"+Misc.stringArrayToString(res, "\n")+"\n");
				return null;
			}
		}
		for (String s: res) {
			if (s.contains("Warning")){
				System.err.println("\nWarning thrown by edgeR application:\n"+Misc.stringArrayToString(res, "\n")+"\n");
			}
		}

		//any results? problems?
		if (edgeRResults.exists() == false) {
			System.err.println("\nError: no edgeR results file found? Check temp files for problems -> "+rOut.toString()+"\n");
			return null;
		}

		res = IO.loadFile(edgeRResults);
		if ((res.length -1) != numGenesToTest) {
			System.err.println("\nError: edgeR results file contains a different number of rows than tested genes?! Aborting, see temp files for issues -> "+rOut.toString()+"\n");
			return null;
		}

		//clean up
		if (deleteTempFiles) {
			workingCountTable.deleteOnExit();
			edgeRResults.deleteOnExit();
			edgeRScript.deleteOnExit();
			rOut.deleteOnExit();
		}

		return res;

	}

	private String createEdgeRANOVAScript(File results) {
		StringBuilder sb = new StringBuilder();
		sb.append("library(edgeR)\n");
		sb.append("rawdata = read.delim('"+workingCountTable.toString()+"', check.names=FALSE, stringsAsFactors=FALSE)\n");
		//created groups
		ArrayList<String> conNames = new ArrayList<String>();
		//for each condition
		for (Condition c: conditions){
			String name = c.getName();
			//for each replica
			for (Replica r: c.getReplicas()) conNames.add(name);		
		}
		sb.append("group = factor(c('"+Misc.stringArrayListToString(conNames, "','")+"'))\n");
		sb.append("y = DGEList(counts=rawdata[,2:length(colnames(rawdata))], genes=rawdata[,1], group=group)\n");
		sb.append("y = calcNormFactors(y)\n");
		sb.append("design = model.matrix(~group,y$samples)\n");
		sb.append("colnames(design) = levels(y$samples$group)\n");
		sb.append("y <- estimateGLMCommonDisp(y,design)\n");
		sb.append("y <- estimateGLMTrendedDisp(y,design)\n");
		sb.append("y <- estimateGLMTagwiseDisp(y,design)\n");
		sb.append("fit <- glmFit(y, design)\n");
		sb.append("lrt = glmLRT(fit, coef=2:length(levels(y$samples$group)))\n");
		sb.append("write.table(as.matrix(topTags(lrt,n=100000000)$table),file='"+results.toString()+"',sep='\\t',quote=FALSE)\n");
		return sb.toString();
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
			//load em with data
			for (int i=0; i< conditions.length; i++) {				
				for (Replica r: conditions[i].getReplicas()){			
					if (verbose) System.out.print("\t"+r.getNameNumber());
					loadReplica(r);
				}
			}

			//any genes with too many reads that were excluded?
			if (flaggedGeneNames.size() !=0) {
				System.err.println("\nWARNING: The following genes/ regions were excluded from the analysis due to one or more bps exceeding the maximum read coverage of "+maxAlignmentsDepth+" . If these are genes/ regions you wish to interrogate, increase the maximum read coverage threshold. " +
						"Realize this will likely require additional memory. Often, these are contaminants (e.g. rRNA) or super abundant transcripts that should be dropped from the analysis.\n"+flaggedGeneNames.toString());
				String[] badGeneNames = Misc.hashSetToStringArray(flaggedGeneNames);

				//remove flagged genes from all replicas
				//for each condition
				for (Condition c: conditions) {
					//for each replica
					for (Replica replica: c.getReplicas()) replica.removeFlaggedGenes(badGeneNames);
				}
			}

			if (verbose) System.out.println();
			//save conditions?
			if (saveCounts) IO.saveObject(serializedConditons, conditions);

		}

		//print conditions
		if (verbose) {
			System.out.println("Conditions, replicas, and mapped counts:");
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
			//if (chrom.equals("chr7") == false) continue;
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
				loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);

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
			loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);
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
		if (maxRepeats !=0){
			Object o = sam.getAttribute("IH");
			if (o != null)  {
				int numRepeats = (Integer)o;
				if (numRepeats > maxRepeats) return true;
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
		File countTable = new File (saveDirectory, "geneCountTable.txt");
		String[] allGenes = new String[genes.length];
		for (int i=0; i< genes.length; i++){
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
		//spacer
		spreadSheetHeader.append("\t");
		//counts
		spreadSheetHeader.append("Counts_"+Misc.stringArrayListToString(conditionNames, "\tCounts_"));
		//spacer
		spreadSheetHeader.append("\t");
		//varCoor
		spreadSheetHeader.append("VarCorCounts_"+Misc.stringArrayListToString(conditionNames, "\tVarCorCounts_"));
		//spacer
		spreadSheetHeader.append("\t");
		//FPKM
		spreadSheetHeader.append("FPKM_"+Misc.stringArrayListToString(conditionNames, "\tFPKM_"));

		//append onto genes
		try {
			BufferedReader counts = new BufferedReader( new FileReader (countTable));
			//skip header
			counts.readLine();
			BufferedReader varCounts = new BufferedReader (new FileReader (varCorrData));
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
				for (Condition c: conditions){
					//for each replica
					for (Replica r: c.getReplicas()){
						double num = 0;
						if (r.getGeneCounts().get(genes[i].getDisplayNameThenName()) != null) num = r.getGeneCounts().get(genes[i].getDisplayNameThenName()).calculateFPKM(r.getTotalCounts(), genes[i].getTotalExonicBasePairs());
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
					if (deleteTempFiles == false){
						rOut.delete();
						scriptFile.delete();
					}
					if (useLocalFitType == false) executeDESeqCluster(countTable, conditionNames, true);
					else throw new IOException();
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
			if (gl == null) continue;

			//load it with counts
			loadGeneLineWithExonCounts(gl, one, two);

			//calc permutated pvalue
			double pval = Num.calculatePermutedChiSquarePValue (Num.floatArraysToInt(gl.getTreatmentExonCounts()), Num.floatArraysToInt(gl.getControlExonCounts()));
			gl.setSplicingPValue((float)pval);

			//set ratio
			setMaxLog2RatioSplice(gl);

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
		for (int i=0; i< conditions.length; i++){
			Condition first = conditions[i];
			for (int j=i+1; j< conditions.length; j++){
				Condition second = conditions[j];
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
		if (verbose) System.out.println(" : "+numberWorkingGenesPassingFilters);

		//differential splice
		if (usePermutationSpliceTest) estimateDifferencesInReadDistributionsWithPermutationTest(first, second);
		else estimateDifferencesInReadDistributions(first, second);

		//print results spread sheet
		if (spreadSheetHeader == null) saveFirstWorkingGeneModels();
		else saveWorkingGeneModels();

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
		//build header line
		spreadSheetHeader = new StringBuilder();
		boolean secondNamePresent = false;
		if (workingGenesToTest[0].getDisplayName() != null && workingGenesToTest[0].getName() != null) {
			spreadSheetHeader.append("#DisplayName\tName\t");
			secondNamePresent = true;
		}
		else spreadSheetHeader.append("#Name\t");
		String splicePVal = "\tSpliceChiPVal_";
		if (usePermutationSpliceTest) splicePVal = "\tSplicePermChiPVal_";
		spreadSheetHeader.append("Chr\tStrand\tStart\tStop\tTotalBPs\tPVal_"+workingName+"\tFDR_"+workingName+"\tVarCorLg2Rto_"+workingName+ splicePVal +workingName+"\tSpliceMaxLg2Rto_"+workingName+"\tSpliceMaxExon_"+workingName);

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
				text.append("\"");
				text.append(genes[i].getName());
				text.append("\"\t");
			}

			//status? OK or sectioned or read depth
			if (genes[i].isFlagged()){
				text.append("Too many reads");
			}

			//coordinates
			text.append(genes[i].coordinates());
			text.append("\t");

			//total bases
			text.append(genes[i].getTotalExonicBasePairs());
			text.append("\t");

			//scores pval, fdr, log2rto, splicePVal, spliceLog2
			text.append(genes[i].getpValue()); 
			text.append("\t");
			text.append(genes[i].getFdr()); 
			text.append("\t");
			text.append(genes[i].getLog2Ratio()); 
			text.append("\t");
			text.append(genes[i].getSplicingPValue()); 
			text.append("\t");
			text.append(genes[i].getSplicingLog2Ratio()); 
			text.append("\t");
			text.append(genes[i].getSplicingExon());

			//set text
			genes[i].setText(text);
		}

	}

	public void saveWorkingGeneModels(){
		//add to header line
		spreadSheetHeader.append("\tPVal_"+workingName+"\tFDR_"+workingName+"\tVarCorLg2Rto_"+workingName+"\tSplicePVal_"+workingName+"\tSpliceMaxLg2Rto_"+workingName+"\tSpliceMaxExon_"+workingName);

		//for each gene
		for (int i=0; i< genes.length; i++){
			StringBuilder text = genes[i].getText();
			text.append("\t");

			//scores pval, fdr, log2rto, splicePVal, spliceLog2
			text.append(genes[i].getpValue()); 
			text.append("\t");
			text.append(genes[i].getFdr()); 
			text.append("\t");
			text.append(genes[i].getLog2Ratio()); 
			text.append("\t");
			text.append(genes[i].getSplicingPValue()); 
			text.append("\t");
			text.append(genes[i].getSplicingLog2Ratio()); 
			text.append("\t");
			text.append(genes[i].getSplicingExon());
		}

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
				//does it log2ratio difference  thresholds?
				float absLog2Ratio = Math.abs(workingGenesToTest[i].getLog2Ratio());
				if (workingGenesToTest[i].getFdr() >= minFDR && absLog2Ratio >= minLog2Ratio){
					geneNamesPassingThresholds.add(workingGenesToTest[i].getDisplayNameThenName());
					numberWorkingGenesPassingFilters++;
					workingGenesPassingFilters.add(workingGenesToTest[i]);
				}
				//set max log2ratio?
				if (workingGenesToTest[i].getMaxAbsLog2Ratio() < absLog2Ratio) workingGenesToTest[i].setMaxAbsLog2Ratio(absLog2Ratio);
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
					if (deleteTempFiles == false){
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

	public void setTNs(Condition first, Condition second){
		workingTNs = new ArrayList<String>();
		workingConditionNames = new ArrayList<String>();
		for (Replica r: first.getReplicas()) {
			workingTNs.add("T");
			workingConditionNames.add(r.getNameNumber());
		}
		for (Replica r: second.getReplicas()) {
			workingTNs.add("N");
			workingConditionNames.add(r.getNameNumber());
		}
	}

	public boolean writePairedCountTable(Condition first, Condition second){
		try {

			HashSet<String> genesToTest = new HashSet<String>();
			Condition[] twoConditions = new Condition[]{first, second};

			//write matrix of name, t1,t2,t3...c1,c2,c3 to file for genes with observations
			workingCountTable = new File(saveDirectory, "countTable_"+workingName+".txt");
			if (deleteTempFiles) workingCountTable.deleteOnExit();
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
						String geneName = it.next();
						//has this been flagged? Hmmm. I thought these were already removed from all of the conditions?
						if (flaggedGeneNames.contains(geneName)) {
							Misc.printErrAndExit("\nFlagged gene found! Exiting\n");
							continue;
						}
						GeneCount gc = geneCounts.get(geneName);
						if (gc != null && gc.getCount()> minimumCounts) {
							genesToTest.add(geneName);
						}
					}
				}

			}
			out.println();

			//any genes to test?
			if (genesToTest.size() ==0) {
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
		name2Gene = new HashMap<String, UCSCGeneLine>();
		for (UCSCGeneLine gl: genes) {
			String name = gl.getDisplayNameThenName();
			if (gl.getExons().length > 1) genesWithExons.put(name, gl);
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
					case 'f': minFDR = Float.parseFloat(args[++i]); break;
					case 'l': minLog2Ratio = Float.parseFloat(args[++i]); break;
					case 'e': minimumCounts = Integer.parseInt(args[++i]); break;
					case 'm': removeOverlappingRegions = true; break;
					case 'n': maxRepeats = Integer.parseInt(args[++i]); break;
					case 'x': maxAlignmentsDepth = Integer.parseInt(args[++i]); break;
					case 'p': performStrandedAnalysis = true; break;
					case 'i': scoreIntrons = true; break;
					case 't': deleteTempFiles = false; break;
					case 'v': filterOutliers = true; break;
					case 'g': genomeVersion = args[++i]; break;
					case 'a': usePermutationSpliceTest = true; break;
					case 'j': performReverseStrandedAnalysis = true; break;
					case 'k': secondStrandFlipped = true; break;
					case 'h': printDocs(); System.exit(0);
					//hidden options!
					case 'd': saveCounts = true; break;
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
		url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";

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
			Misc.printErrAndExit("\nError: Please upgrade R, Bioconductor, DESeq and edgeR to the latest versions, see http://www.bioconductor.org \n");
		}

		//look for bed file
		if (refSeqFile == null && bedFile == null){
			Misc.printErrAndExit("\nPlease enter a regions or gene file to use in scoring regions.\n");
		}
	}	

	/**Returns true if DESeq uses the estimateDispersions function otherwise false.*/
	public static boolean estimateDispersions(File fullPathToR, File tempDir){
		//this will return an error
		String error = IO.runRCommandLookForError("library(DESeq); estimateDispersions(5); library(edgeR)", fullPathToR, tempDir);
		if (error.contains("could not find function")) {
			System.err.println("\nWARNING: Please install and update R, Bioconductor, DESeq, and edgeR to latest versions.\n");
			return false;
		}
		return true;
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Defined Region Differential Seq: July 2013                **\n" +
				"**************************************************************************************\n" +
				"DRDS takes bam files, one per replica, minimum one per condition, minimum two\n" +
				"conditions (e.g. treatment and control or a time course/ multiple conditions) and\n" +
				"identifies differentially expressed genes under any pairwise comparison using DESeq.\n" +
				"DESeq's variance corrected count data is used to heirachically cluster the \n" +
				"differentially expressed genes as well as the samples. See the DESeq manual for\n" +
				"details. An ANOVA-like any condition differential expression analysis is also run using\n" +
				"edgeR. Alternative splicing is estimated using a chi-square test of independence.\n" +
				"In addition to the cluster plots, a spread sheet is created with the pValue,\n" +
				"FDR, and variance corrected log2Ratios for each of the pairwise comparisons as well as\n" +
				"the raw, FPKM, and log2 variance corrected alignment counts.  Use the later for\n" +
				"subsequent clustering and distance estimations.\n"+

				"\nOptions:\n"+
				"-s Save directory.\n"+
				"-c Conditions directory containing one directory for each condition with one xxx.bam\n" +
				"       file per biological replica and their xxx.bai indexs. 3-4 reps recommended per\n" +
				"       condition. The BAM files should be sorted by coordinate using Picard's SortSam.\n" +
				"       All spice junction coordinates should be converted to genomic coordinates, see\n" +
				"       USeq's SamTranscriptomeParser.\n" +
				"-r Full path to R (version 3+) loaded with DESeq and edgeR, defaults to '/usr/bin/R'\n" +
				"       file, see http://www.bioconductor.org . Type 'library(DESeq) and library(edgeR)'\n" +
				"       in an R terminal to see if they are installed. \n"+
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
				"-v Filter for variance outliers in DESeq, defaults to not filtering.\n"+
				"-m Mask overlapping gene annotations, recommended for well annotated genomes.\n"+
				"-x Max per base alignment depth, defaults to 50000. Genes containing such high\n"+
				"       density coverage are ignored. Warnings are thrown.\n"+
				"-n Max number repeat alignments. Defaults to all.  Assumes 'IH' tags have been set by\n" +
				"       processing raw alignments with the SamTranscriptomeProcessor.\n"+
				"-f Minimum FDR threshold for sorting, defaults to 13 (-10Log10(FDR=0.05)).\n"+
				"-l Minimum absolute varCorLog2Rto threshold for sorting, defaults to 1 (2x).\n"+
				"-e Minimum number alignments per gene/ region, defaults to 20.\n"+
				"-i Score introns instead of exons.\n"+
				"-p Perform a stranded analysis. Only collect reads from the same strand as the\n" +
				"      annotation.\n" +
				"-j Reverse stranded analysis.  Only collect reads from the opposite strand of the\n" +
				"      annotation.  This setting should be used for the Illumina's strand-specific\n" +
				"      dUTP protocol.\n" +
				"-k Second read's strand is flipped. Otherwise, assumes this was not done in the \n" +
				"      SamTranscriptomeParser.\n" +
				"-a Perform a permutation based chi-square test for differential exon usage.\n"+
				"      Needs 4 or more replicas per condition.\n"+
				"-t Don't delete temp files (R script, R results, Rout, etc..).\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/DefinedRegionDifferentialSeq -c\n" +
				"      /Data/TimeCourse/ESCells/ -s /Data/TimeCourse/DRDS -g H_sapiens_Feb_2009\n" +
				"     -u /Anno/mergedHg19EnsemblGenes.ucsc.gz\n\n" +

		"**************************************************************************************\n");

	}
}
