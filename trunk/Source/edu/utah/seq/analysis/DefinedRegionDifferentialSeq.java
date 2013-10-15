package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import org.apache.poi.ss.usermodel.*;
import org.apache.poi.ss.*;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.analysis.multi.Condition;
import edu.utah.seq.analysis.multi.GeneCount;
import edu.utah.seq.analysis.multi.GeneResult;
import edu.utah.seq.analysis.multi.PairedCondition;
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
	private int minimumCounts = 10;
	private boolean deleteTempFiles = true;
	private boolean filterOutliers = false;
	private float minimumSpliceLog2Ratio = 1f;
	private boolean performStrandedAnalysis = false;
	private boolean performReverseStrandedAnalysis = false;
	private int maxRepeats = 0;
	private boolean verbose = true;
	private int maxAlignmentsDepth = 50000;
	private boolean secondStrandFlipped = false;

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
	public boolean saveCounts = false;
	
	//for loading data
	private int workingFragmentNameIndexPlus = 1;
	private int workingFragmentNameIndexMinus = -1;
	private HashMap<String, Integer> workingFragNameIndex = new HashMap<String, Integer>(10000);
	
	//for paired diff expression
	private File geneCountTable;
	private String[] conditionNamesPerReplica;
	private String[] replicaNames;
	private PairedCondition[] pairedConditions;
	private File blindedVarCorData;

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
			boolean verbose, boolean reverse, boolean flipped, int maxAlignDepth){
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

		//load count data by replica 
		loadConditions();

		//write out count table of genes passing minimum counts
		writeCountTable();
		
		//cluster samples based on blinded normalization
		executeDESeqCluster(false);
		
		//launch DESeq for differential expression
		String thresholdName = "-10Log10(FDR) "+Num.formatNumber(minFDR,1)+", Log2Ratio "+Num.formatNumber(minLog2Ratio, 1);
		String varianceFiltered = "without";
		if (filterOutliers) varianceFiltered = "with";
		if (verbose) System.out.println("Running pairwise DESeq analysis "+varianceFiltered+" variance outlier filtering to identify differentially expressed genes...");
		analyzeForDifferentialExpression(false);
		if (verbose) System.out.println("\n\t"+geneNamesPassingThresholds.size()+" / "+geneNamesToAnalyze.length+"\tgenes differentially expressed ("+thresholdName+")\n");
			
		//launch all pair ANOVA-like edgeR analysis, appends edgeR FDR onto gene line
		if (conditions.length > 2){
			if (verbose) System.out.println("Running ANOVA-like edgeR all condition analysis...");
			if (runEdgeRAllConditionAnalysis() == false) System.err.println("Skipped edgeR analysis!\n");
		}
		
		//launch Diff splice
		if (verbose) System.out.println("Running pairwise chi-square tests for differential splicing...");
		analyzeForDifferentialSplicing();
		
		//set max deseqFDR and varcorLog2Ratio in geneNamesToAnalyze
		setMaxScores();

		//print spreadsheet 
		printStatSpreadSheet();


	}

	private void setMaxScores() {
		float[] maxFDRs = new float[geneNamesToAnalyze.length];
		float[] maxLog2 = new float[geneNamesToAnalyze.length];
		//for each pair
		for (int i=0; i< pairedConditions.length; i++){
			//gene : float[]{FDR, log2VarCorRto}
			float[][] scores = pairedConditions[i].getParsedDiffExpResults();
			for (int j=0; j< geneNamesToAnalyze.length; j++){
				if (scores[j][0] > maxFDRs[j]) maxFDRs[j] = scores[j][0];
				float abs = Math.abs(scores[j][1]);
				if (abs > maxLog2[j]) maxLog2[j] = abs;
			}
		}
		//set in genes
		for (int i=0; i< geneNamesToAnalyze.length; i++){
			UCSCGeneLine gene = name2Gene.get(geneNamesToAnalyze[i]);
			gene.setMaxAbsLog2Ratio(maxLog2[i]);
			gene.setFdr(maxFDRs[i]);  //note misnaming!
		}
	}

	private boolean runEdgeRAllConditionAnalysis() {

		//create and execute script
		String[] results = launchEdgeRANOVA();

		if (results != null){
			//parse results and append FDRs to UCSCGeneLines
			parseEdgeRANOVAResults(results);
		}
		else return false;

		return true;
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

	private String[] launchEdgeRANOVA() {
		//create & write out edgeR script
		File edgeRResults = new File(saveDirectory, "edgeRANOVA_Results.txt");
		String script = createEdgeRANOVAScript(edgeRResults);
		File edgeRScript = new File(saveDirectory, "edgeRANOVA_RScript.txt");
		IO.writeString(script, edgeRScript);
		File rOut = new File(saveDirectory, "edgeRANOVA_RScript.txt.Rout");

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
		if ((res.length -1) != geneNamesToAnalyze.length) {
			System.err.println("\nError: edgeR results file contains a different number of rows than tested genes?! Aborting, see temp files for issues -> "+rOut.toString()+"\n");
			return null;
		}

		//clean up
		if (deleteTempFiles) {
			edgeRResults.deleteOnExit();
			edgeRScript.deleteOnExit();
			rOut.deleteOnExit();
		}

		return res;

	}

	private String createEdgeRANOVAScript(File results) {
		StringBuilder sb = new StringBuilder();
		sb.append("library(edgeR)\n");
		sb.append("rawdata = read.delim('"+geneCountTable+"', check.names=FALSE, stringsAsFactors=FALSE)\n");
		sb.append("group = factor(c('"+Misc.stringArrayToString(conditionNamesPerReplica, "','")+"'))\n");
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

	private void analyzeForDifferentialSplicing() {
		//for each PairedCondition
		for (PairedCondition pc : pairedConditions){
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



	public void executeDESeqCluster(boolean useLocalFitType){
		File sampleClustering = new File (saveDirectory, "sampleClusterPlot.pdf");
		blindedVarCorData = new File (saveDirectory, "blindedDESeqVarCorData.txt");
		if (deleteTempFiles) blindedVarCorData.deleteOnExit();
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+geneCountTable.getCanonicalPath()+"', header=TRUE)\n");
			sb.append("rownames(countsTable) = countsTable[,1]\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c('"+ Misc.stringArrayToString(replicaNames, "','") + "')\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("cds = estimateSizeFactors( cds )\n");
			if (useLocalFitType) sb.append("cds = estimateDispersions( cds, method='blind', sharingMode='fit-only', fitType='local' )\n");
			else sb.append("cds = estimateDispersions( cds, method='blind', sharingMode='fit-only' )\n");
			sb.append("vsd = getVarianceStabilizedData( cds )\n");
			//write out the vsd data
			sb.append("write.table(vsd, file = '"+blindedVarCorData.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			//cluster by sample
			sb.append("dists = dist( t( vsd ) )\n");
			sb.append("pdf('"+ sampleClustering.getCanonicalPath() +"', height=10, width=10)\n");
			sb.append("heatmap( as.matrix( dists ), symm=TRUE, scale='none', margins=c(10,10), col = colorRampPalette(c('darkblue','white'))(100), labRow = paste( pData(cds)$condition, pData(cds)$type ) )\n");
			sb.append("dev.off()\n");

			//write script to file
			File scriptFile = new File (saveDirectory, "cluster_RScript.txt");
			File rOut = new File(saveDirectory, "cluster_RScript.txt.Rout");
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
					if (useLocalFitType == false) executeDESeqCluster(true);
					else throw new IOException();
				}
			}
			//any problems?
			if (sampleClustering.exists() == false || blindedVarCorData.exists() == false) throw new IOException();

			//cleanup
			if (deleteTempFiles) {
				rOut.deleteOnExit();
				scriptFile.deleteOnExit();
			}

		} catch (IOException e) {
			Misc.printErrAndExit("\nError failed to cluster data. Check temp files in save directory for error.\n"+ e.getMessage());
		}
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
		
		//extract info and put back in PairedCondition
		HashMap<String, float[]> geneNameSpliceStats = new HashMap<String, float[]>();
		for (UCSCGeneLine gene: genesWithExonsAndReads){
			float[] scores = new float[]{gene.getSplicingPValue(), gene.getSplicingLog2Ratio(), gene.getSplicingExon()};
			geneNameSpliceStats.put(gene.getDisplayNameThenName(), scores);
		}
		pair.setGeneNameSpliceStats(geneNameSpliceStats);
		
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


	public void analyzeForDifferentialExpression(boolean useLocalFitType){
		ArrayList<PairedCondition> pairedConditionsAL = new ArrayList<PairedCondition>();
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+geneCountTable.getCanonicalPath()+"', header=TRUE)\n");
			sb.append("rownames(countsTable) = countsTable[,1]\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c('"+ Misc.stringArrayToString(conditionNamesPerReplica, "','") + "')\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("cds = estimateSizeFactors( cds )\n");
			//estimate dispersions
			//no replicas?
			boolean noReplicas = true;
			for (int i=0; i< conditions.length; i++){
				if (conditions[i].getReplicas().length > 1){
					noReplicas = false;
					break;
				}
			}
			if (noReplicas) {
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
			sb.append("vsd = getVarianceStabilizedData( cds )\n");
			
			//write out variance stabilized data, actually we don't want this, need it to be blinded!
			//deseqVarCorData = new File(saveDirectory,"deseqVarCorData.txt");
			//if (deleteTempFiles) deseqVarCorData.deleteOnExit();
			//sb.append("write.table(vsd, file = '"+deseqVarCorData.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			
			//diff express
			pairedConditionsAL = new ArrayList<PairedCondition>();
			for (int i=0; i< conditions.length; i++){
				for (int j=i+1; j< conditions.length; j++){
					//make condition
					PairedCondition pc = new PairedCondition(conditions[i], conditions[j]);
					pairedConditionsAL.add(pc);
					File results = new File (saveDirectory, pc.getName()+"DESeq.txt");
					if (deleteTempFiles) results.deleteOnExit();
					pc.setDiffExpResults(results);
					
					sb.append("res = nbinomTest( cds, '"+conditions[i].getName()+"', '"+conditions[j].getName()+"', pvals_only = FALSE)\n");
					sb.append("res[,6] = (rowMeans( vsd[, conditions(cds)=='"+conditions[i].getName()+"', drop=FALSE] ) - rowMeans( vsd[, conditions(cds)=='"+conditions[j].getName()+"', drop=FALSE] ))\n");
					//Fred adjP
					sb.append("res[,8] = -10 * log10(res[,8])\n");
					//Parse padj, log2ratio; note flip of A and B back to T and C
					sb.append("res = res[,c(8,6)]\n");
					//note, the order of the rows is the same as the input
					sb.append("write.table(res, file = '"+results.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
				}
			}
			
			//write script to file
			File scriptFile = new File (saveDirectory,"deseq_RScript.txt");
			File rOut = new File(saveDirectory, "deseq_RScript.txt.Rout");
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
					analyzeForDifferentialExpression(true);
					return;
				}
			}

			//Make conditions
			pairedConditions = new PairedCondition[pairedConditionsAL.size()];
			pairedConditionsAL.toArray(pairedConditions);
			
			//look for results files and parse results
			for (PairedCondition pc: pairedConditions){
				if (pc.getDiffExpResults().exists() == false ) throw new IOException("\n\nR results file doesn't exist. Check temp files in save directory for error.\n");
				pc.parseDESeqStatResults(geneNamesToAnalyze, minFDR, this.minLog2Ratio);
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
			Misc.printErrAndExit("Error: failed to execute DESeq.\n");
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
		    
		    //add varcor values
		    cellIndex = addVarCorInfo(rows, cellIndex);
		    
		    //add FPKM
		    cellIndex = addFPKMInfo(rows, cellIndex);
		    
		    //new sheet for diff expressed genes at relaxed thresholds
		    Sheet sheet1 = workbook.createSheet("FDR 10 Lg2Rt 0.585");
		    sheet1.createFreezePane( 0, 1, 0, 1 );
		    addDiffExpressedGenes(sheet1, 10f, 0.5849625f);
		    
		    //new sheet for diff expressed genes at standard thresholds
		    Sheet sheet2 = workbook.createSheet("FDR 13 Lg2Rt 1");
		    sheet2.createFreezePane( 0, 1, 0, 1 );
		    addDiffExpressedGenes(sheet2, 13f, 1f);
		    
		    //new sheet for diff expressed genes at standard thresholds
		    Sheet sheet3 = workbook.createSheet("FDR 30 Lg2Rt 1.585");
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

	private void addDiffExpressedGenes(Sheet sheet, float minFDR, float minLog2Ratio) {
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
				if (Math.abs(fdrLog2[1]) >= minLog2Ratio && fdrLog2[0] >= minFDR){
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
		for (int i=0; i< descriptors.length; i++){
			Row row = sheet.createRow(rowIndex++);
			Cell cell = row.createCell(0);
			cell.setCellStyle(bold);
			cell.setCellValue(descriptors[i][0]);
			if (descriptors[i][1] != null) row.createCell(1).setCellValue(descriptors[i][1]);
		}
		
		sheet.autoSizeColumn((short)0);
	}
	
	private String[][] descriptors = {
			{"Gene/Region info:", null},
			{"IGB HyperLink", "Click each to reposition a running instance of IGB, http://bioviz.org/igb/"},
			{"Alt Name", "Alternative name for this gene/region"},
			{"Coordinates", "Chromosome : Start of the gene/region - end of the gene/region"},
			{"Strand", "Strand of the gene/region"},
			{"Total BPs", "Total base pairs of gene/region after masking for overlapping annotations if so indicated"},
			{"Max Abs Lg2Rto", "Maximum absolute log2 ratio observed using DESeq's varCor values"},
			{"Max DESeq FDR", "Maximum observed -10Log10( FDR ) between any pairwise DESeq comparisons.  Note the MaxAbsLg2Rto and MaxDESeqFDR may be from different comparisons.  "
					+ "Moreover, all p-values and FDRs in USeq output, both printed and in graph form have been phred transformed where the value is actually -10*Log10(FDR or pval). Thus -10*Log10(0.001) = 30, -10*Log10(0.01) = 20, -10*Log10(0.05) = 13, and -10*Log10(0.1) = 10. "},
			{"EdgeR FDR", "Maximum -10Log10( FDR ) from running edgeR's ANOVA-like analysis designed to score genes for differential expression under any conditions."},
			{"", ""},
			{"For each Paired Comparison:", null},
			{"Lg2Rto", "DESeq's log2 ratio for differential expression, log2(mean(varCorT1,varCorT2, ...)/mean(varCorC1,varCorC2,...)). Note these varCor values will differ somewhat from the blinded varCor values exported to the spreadsheet since they use a different variance estimation."},
			{"FDR", "DESeq's negative binomial p-value following the Benjamini and Hochberg multiple testing correction."},
			{"Spli Lg2Rto", "The maximum log2 ratio difference between two conditions and a gene's exon after correcting for differences in counts.  Use this to gauge the degree of potential differential splicing."},
			{"Spli AdjPVal", "Chi-square test of independence between two conditions and a gene's exons that contain 5 or more counts.  A Bonferroni correction is applied to the p-values."},
			{"Spli Index", "The zero based index relative to the plus genomic strand for the exon displaying the SpliLg2Rto."},
			{"", ""},
			{"For each Replica:", null},
			{"Counts", "Number of alignments that overlap a given gene/region by at least one base pair. Paired alignments are only counted once.  Alignments that span a gene/region but don't actually touch down are ignored. "},
			{"VarCor", "Per gene/region variance corrected counts in log space from DESeq estimated using its blind method.  Use these for subsequent clustering.  These values are a bit different from the varCor values used in scoring paired differential expression. "},
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

	private int addVarCorInfo(Row[] rows, int cellIndex) {
		//First parse the table of numbers
		float[][] geneVarCor = parseVarCorData(blindedVarCorData);
				
		//write them to the excel sheet
		int ci = cellIndex;
		for (int j=0; j< geneNamesToAnalyze.length; j++){
			int rowIndex = j+1;
			ci = cellIndex;
			for (float c: geneVarCor[j]){
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

	/*Rips the deseqVarCorData file checking and fixing negative values.*/
	private float[][] parseVarCorData(File deseqVarCorData) {
		Pattern tab = Pattern.compile("\\t");
		float[][] geneVarCor = new float[geneNamesToAnalyze.length][replicaNames.length];
		float minimum = 0;
		try {
			BufferedReader in = IO.fetchBufferedReader(deseqVarCorData);
			for (int i=0; i< geneNamesToAnalyze.length; i++){
				String[] scoresString = tab.split(in.readLine());
				if (scoresString.length != replicaNames.length) Misc.printErrAndExit("\nError: one of the varCor data rows doesn't contain the appropriate number of values! See row "+i+" "+deseqVarCorData);
				float[] scores = Num.parseFloats(scoresString);
				if (scores == null) Misc.printErrAndExit("\nError: couldn't parse floats from one of the varCor data rows! See row "+i+" "+deseqVarCorData);
				//look for negative values
				for (float f : scores) if (f< minimum) minimum = f;
				geneVarCor[i] = scores;
			}
			in.close();
		} catch (Exception e) {
			Misc.printErrAndExit("\nError: problem parsing varCor data, see "+deseqVarCorData+"\n"+e.getMessage());
		}
		
		//fix negatives? remember these are in log space so it is OK to add a constant
		if (minimum < 0){
			minimum = -1 * minimum;
			for (int i=0; i< geneVarCor.length; i++){
				for (int j=0; j< geneVarCor[i].length; j++) geneVarCor[i][j] += minimum;
			}
		}
		return geneVarCor;
	}

	/*Adds DiffExpLg2Rto DiffExpFDR DiffSpliceLg2Rto DiffSpliceFDR DiffSpliceIndex*/
	private int addPairedConditionInfo(Row[] rows, int cellIndex) {
		float[] noScores = {0f,0f,0f};
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
				//index splice
				rows[rowIndex].createCell(ci++).setCellValue(splice[2]);
			}
			cellIndex = ci;
		}
		return ci;
	}

	/*Adds DisplayName Name Chr Strand Start Stop TotalBPs MaxAbsLg2Rto MaxDESeqFDR MaxEdgeRFDR*/
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
			//MaxDESeqFDR
			rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getFdr());
			//MaxEdgeRFDR
			rows[rowIndex].createCell(cellIndex++).setCellValue(gene.getFdrEdgeR());
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
	    	index = createCell(row, "Max DESeq FDR", index, headerCellStyle);
	    	index = createCell(row, "EdgeR FDR", index, headerCellStyle);
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
	    	index = createCell(row, "FDR "+name, index, pcCellStyle);
	    	if (all){
	    		index = createCell(row, "Spli Lg2Rto "+name, index, pcCellStyle);
	    		index = createCell(row, "Spli AdjPVal "+name, index, pcCellStyle);
	    		index = createCell(row, "Spli Index "+name, index, pcCellStyle);
	    	}
	    }
	    
	    if (all){
	    	//count values for each replica
	    	CellStyle countCellStyle = createHeaderStyle(IndexedColors.LIGHT_CORNFLOWER_BLUE.getIndex());
	    	for (String repName: replicaNames){
	    		index = createCell(row, "Counts "+repName, index, countCellStyle);
	    	}

	    	//VarCor values for each replica
	    	CellStyle varCorCellStyle = createHeaderStyle(IndexedColors.WHITE.getIndex());
	    	for (String repName: replicaNames){
	    		index = createCell(row, "VarCor "+repName, index, varCorCellStyle);
	    	}

	    	//FPKM values for each replica	    
	    	CellStyle fpkmCellStyle = createHeaderStyle(IndexedColors.SKY_BLUE.getIndex());
	    	for (String repName: replicaNames){
	    		index = createCell(row, "FPKM "+repName, index, fpkmCellStyle);
	    	}
	    }
	}

	public float parseFloat(String f){
		if (f.equals("Inf") || f.equals("NA")) return Float.MIN_VALUE;
		else return Float.parseFloat(f);
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
					out.print(conditionName);
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
						int num = 0;
						//remember that genes with zero counts are not stored in a replica
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
		
		conditionNamesPerReplica = Misc.stringArrayListToStringArray(conditionNamesPerReplicaAL);
		replicaNames = Misc.stringArrayListToStringArray(replicaNamesAL);
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
				"**                       Defined Region Differential Seq: Sept 2013                 **\n" +
				"**************************************************************************************\n" +
				"DRDS takes sorted bam files, one per replica, minimum one per condition, minimum two\n" +
				"conditions (e.g. treatment and control or a time course/ multiple conditions) and\n" +
				"identifies differentially expressed genes under any condition using edgeR's ANOVA-like\n"+
				"analysis and an all pairwise comparison with DESeq. DESeq's blinded variance corrected\n"+
				"count data is used to heirachically cluster the samples. Alternative splicing\n"+
				"is estimated using a chi-square test of independence. Note, when interested in only a\n"+
				"few genes or regions, append these onto a full gene table so that DESeq can\n"+
				"appropriately estimate the library size and replica variance.\n"+

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
				"       density coverage are ignored.\n"+
				"-n Max number repeat alignments. Defaults to all.  Assumes 'IH' tags have been set by\n" +
				"       processing raw alignments with the SamTranscriptomeProcessor.\n"+
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
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/DefinedRegionDifferentialSeq -c\n" +
				"      /Data/TimeCourse/ESCells/ -s /Data/TimeCourse/DRDS -g H_sapiens_Feb_2009\n" +
				"     -u /Anno/mergedHg19EnsemblGenes.ucsc.gz\n\n" +

		"**************************************************************************************\n");

	}
}
