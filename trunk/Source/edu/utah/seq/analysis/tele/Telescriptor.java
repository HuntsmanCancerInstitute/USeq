package edu.utah.seq.analysis.tele;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import htsjdk.samtools.*;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.gen.*;
import edu.utah.seq.analysis.DefinedRegionDifferentialSeq;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;

/** Application for identifying genes with poor splicing and telescripting.
 * @author Nix
 * 
 * */
public class Telescriptor {

	/* Issues!
	 * Genes with an un annotated exon produce a misSplice read for every alignment to the missing exon.*/

	//user defined fields
	private File[] treatmentBamFiles;
	private File[] controlBamFiles;
	private File refSeqFile;
	private File resultsDirectory;
	private File geneResultsDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private int minimumGeneReadCoverage = 50;
	private int minimumBaseReadCoverage = 10;
	private int minimumTranscriptLength = 250;
	private int length5PrimeWindowScan = 125;
	private double fractionLengthBackground = 0.5;
	private int minimumWindowAndBackgroundCount = 25;
	private double minimumSkew = 2;
	private boolean isStranded = true;

	//internal fields
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private ArrayList<String> skippedGenes = new ArrayList<String>();
	private SamReader[] treatmentSamReaders;
	private SamReader[] controlSamReaders;
	private SamReader[] samReaders;
	private Gzipper treatmentMisSplicedSamOut;
	private Gzipper controlMisSplicedSamOut;
	private StringBuilder rScript = new StringBuilder("library(ggplot2)\n");
	public static final Pattern CIGAR = Pattern.compile("(\\d+)([MND])");
	public ArrayList<TeleGene> scoredGenes = new ArrayList<TeleGene>();
	private long startTime;

	//working fields
	private String chromosome;
	private UCSCGeneLine[] genes;
	private HashSet<String> unSplicedTreatmentAlignments = new HashSet<String>();
	private HashSet<String> unSplicedControlAlignments = new HashSet<String>();
	private HashSet<String> unSplicedAlignments;
	private ArrayList<String>[] treatmentReadCoverage; 
	private ArrayList<String>[] controlReadCoverage; 
	private ArrayList<String>[] readCoverage;
	private ArrayList<String>[] toCheck = new ArrayList[6];
	private UCSCGeneLine gene;
	private int firstGeneBase;
	private int lastGeneBase;
	private TeleStats treatmentGeneStats;
	private TeleStats controlGeneStats;
	

	//constructor
	public Telescriptor(String[] args){
		startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load genes
		System.out.println("Loading gene models...");
		loadGeneModels();

		//create sam readers and gzippers
		System.out.println("\nOpening SAM readers and writers...");
		openSams();
		writeOutSamHeaders();

		//walk through each chromosome of genes building ReadCoverageScoredTranscript
		System.out.println("\nWalking chromosomes and genes...");
		walkChromosomes();

		//close readers and writers
		closeSams();

		//analyze genes with chiSquare tests
		System.out.println("\n# Genes to Analyze, "+scoredGenes.size());
		if (scoredGenes.size() == 0) Misc.printExit("\nNo Genes to analyze? Aborting!\n");
		analyseGenes();
		
		//print spreadsheet results
		printSpreadSheet();

		//execute R graphing
		System.out.println("\nGenerating GGPlots...");
		runRGraphics();

		//sort and index misspliced reads
		System.out.println("\nSorting and writing mis-spliced sam alignments to file...");
		new PicardSortSam(treatmentMisSplicedSamOut.getGzipFile());
		new PicardSortSam(controlMisSplicedSamOut.getGzipFile());
		
		//write out skipped genes
		IO.writeArrayList(skippedGenes, new File (resultsDirectory, "skippedGenes.txt"));

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes\n");
	}

	private void runRGraphics() {
		String errors = IO.runRCommandLookForError(rScript.toString(), fullPathToR, resultsDirectory);
		if (errors == null || errors.length() !=0){
			Misc.printErrAndExit("\nError: Found when generating relative read coverage graphs in R. See R error message:\n\t\t"+errors+"\n\n");
		}

	}

	private void writeOutSamHeaders() {
		try {
			String head = treatmentSamReaders[0].getFileHeader().getTextHeader().trim();
			treatmentMisSplicedSamOut.println(head);
			head = controlSamReaders[0].getFileHeader().getTextHeader().trim();
			controlMisSplicedSamOut.println(head);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void analyseGenes() {

			//run chisquare test on spiced vs misSplice for t and c
			System.out.println("\tRunning spliced vs mis-spliced chiSquare test...");
			calcDiffMisSplicePVal();

			//run chisquare test on spiced vs misSplice for t and c
			System.out.println("\tRunning read count skew chiSquare test...");
			calcDiffSkewSplicePVal();
			
			//run chisquare test on 3'UTR diff for t and c
			System.out.println("\tRunning read count 3'UTR chiSquare test...");
			calcDiffUTRPVal();
			
	} 

	
	private void printSpreadSheet() {
		try {
			//make writer for spreadsheet and write header
			PrintWriter spreadSheetOut = new PrintWriter( new FileWriter( new File (resultsDirectory, "teleStatsSummary.xls")));
			spreadSheetOut.println("GeneName\tAltName: Description\tcDNA Length\t"
					+ "# A align\t# A Unspliced\t# B align\t# B Unspliced\tlog2RtoMisSplice\tpAdjMisSplice\t"
					+ "# A 3'UTR\t# B 3'UTR\tlog2Rto((AUtrFPKM/NonUtrFPKM)/ (BUtrFPKM/NonUtrFPKM))\tpAdj UTR\t"
					+ "5' Index\tA 5' Median\t"
					+ "A 3' Median\tA Median Log2Rto\tB 5' Median\tB 3' Median\tB Median Log2Rto"
					+ "\tMedian Log2 (aSkew/ bSkew)\tA 5' Count\tA 3' Count\tA Count Log2Rto\tA Bkgrnd Coeff Var\tB 5' Count"
					+ "\tB 3' Count\tB Count Log2Rto\tB Bkgrnd Coeff Var\tCount Log2 (aSkew/ bSkew)\tpAdj Count Skew");

			for (TeleGene g : scoredGenes){
				spreadSheetOut.println(g);
			}
			
			spreadSheetOut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void calcDiffMisSplicePVal(){
		//collect counts
		int numGenes = scoredGenes.size();
		int[][] treatment = new int[numGenes][2];
		int[][] control = new int[numGenes][2];
		
		for (int i=0; i< numGenes; i++){
			TeleGene tg = scoredGenes.get(i);
			TeleStats tts = tg.getTreatment();
			TeleStats cts = tg.getControl();
			int numMisT = tts.getNumberUnsplicedAlignments();
			int numNonMisT = Math.abs(tts.getNumberExonicAlignments() - numMisT);
			int numMisC = cts.getNumberUnsplicedAlignments();
			int numNonMisC = Math.abs(cts.getNumberExonicAlignments() - numMisC);
			
			treatment[i] = new int[]{numMisT, numNonMisT};
			control[i] = new int[]{numMisC, numNonMisC};

			//calc log2Rto
			double t = (double)(numMisT+1)/(double)(tts.getNumberExonicAlignments()+1);
			double c = (double)(numMisC+1)/(double)(cts.getNumberExonicAlignments()+1);
			double lrto = Num.log2( t / c );
			tg.setMisSpliceLog2Rto(lrto);
		}

		//estimate chi-square pvalues using R for resolution of extreemly small p-values, radiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, resultsDirectory, fullPathToR, true);

		//bonferroni correction
		double bc = Num.minus10log10(numGenes);

		//add back
		for (int i=0; i< numGenes; i++){
			//set corrected p-value 
			double pAdj = pVals[i] + bc;
			if (pAdj > 0) scoredGenes.get(i).setTransMisSplicePVal(pAdj);
		}
	} 
	
	private void calcDiffUTRPVal(){
		//collect counts
		int numGenes = scoredGenes.size();
		int[][] treatment = new int[numGenes][2];
		int[][] control = new int[numGenes][2];
		
		for (int i=0; i< numGenes; i++){
			TeleGene tg = scoredGenes.get(i);
			TeleStats tts = tg.getTreatment();
			TeleStats cts = tg.getControl();
			int numUtrT = tts.getNumber3UTRAlignments();
			int numNonUtrT = Math.abs(tts.getNumberExonicAlignments()-numUtrT);
			int numUtrC = cts.getNumber3UTRAlignments();
			int numNonUtrC = Math.abs(cts.getNumberExonicAlignments()-numUtrT);
			
			treatment[i] = new int[]{numUtrT, numNonUtrT};
			control[i] = new int[]{numUtrC, numNonUtrC};

			//calc log2Rto
			//double t = (double)(numUtrT+1)/(double)(tts.getNumberExonicAlignments()+1);
			double t = tts.getUtrFpkm() / tts.getNonUtrFpkm();
			double c = cts.getUtrFpkm() / cts.getNonUtrFpkm();
			//double c = (double)(numUtrC+1)/(double)(cts.getNumberExonicAlignments()+1);
			double lrto = Num.log2( t / c );
			tg.setUtrLog2Rto(lrto);
		}

		//estimate chi-square pvalues using R for resolution of extreemly small p-values, radiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, resultsDirectory, fullPathToR, true);

		//bonferroni correction
		double bc = Num.minus10log10(numGenes);

		//add back
		for (int i=0; i< numGenes; i++){
			//set corrected p-value 
			double pAdj = pVals[i] + bc;
			if (pAdj > 0) scoredGenes.get(i).setpAdjUTR(pAdj);
		}
	} 


	
	private void calcDiffSkewSplicePVal(){
		//count number of trans
		int numGenes = scoredGenes.size();

		int[][] treatment = new int[numGenes][2];
		int[][] control = new int[numGenes][2];
		for (int i=0; i< numGenes; i++){
			TeleGene tg = scoredGenes.get(i);
			TeleStats tt = tg.getTreatment();
			treatment[i] = new int[]{tt.getCountWindow(), tt.getCountBackground()};
			tt = tg.getControl();
			control[i] = new int[]{tt.getCountWindow(), tt.getCountBackground()};	
		}

		//estimate chi-square pvalues using R for resolution of extremely small p-values, ridiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, resultsDirectory, fullPathToR, true);

		//bonferroni correction
		double bc = Num.minus10log10(numGenes);

		//add back
		for (int i=0; i< numGenes; i++){
			TeleGene tg = scoredGenes.get(i);
			double pAdj = pVals[i] + bc;
			if (pAdj > 0) tg.setpAdjSkewedReadCount(pAdj);
		}
	}

	private void walkChromosomes() {
		for (String chromName: chromGenes.keySet()){
			chromosome = chromName;
			genes = chromGenes.get(chromosome);
			System.out.print("\t"+chromosome+"\t"+genes.length + " ");
			int counter = 0;
			//for each gene
			for (UCSCGeneLine g: genes){
				if (counter++ > 500) {
					System.out.print(".");
					counter = 0;
				}
				gene = g;
				firstGeneBase = gene.getTxStart();
				lastGeneBase = gene.getTxEnd();
				
				//long enought?
				if (gene.getTotalExonicBasePairs() < minimumTranscriptLength){
					skippedGenes.add(gene.getNames("_"));
					continue;
				}

				//load sam alignments for the gene, skip genes with too low counts
				if (loadTCSams() == false){
					skippedGenes.add(gene.getNames("_"));
					continue;
				}
			
				//score genes, fails if too few counts in windows.
				if (scoreGene() == false){
					skippedGenes.add(gene.getNames("_"));
					continue;
				};
				
				TeleGene ftg = new TeleGene(gene, treatmentGeneStats, controlGeneStats);
				
				//check scores?
				if (minimumSkew != 0) {
					//check treatment skew
					double skewMed = ftg.getTreatment().getMedianSkewLog2Rto();
					double skewCount = ftg.getTreatment().getCountSkewLog2Rto();
					if (skewMed < minimumSkew && skewCount < minimumSkew) {
						skippedGenes.add(gene.getNames("_"));
						continue;
					}
					//check treatment/ control skew
					skewMed = ftg.getMedianSkewLog2Rto();
					skewCount = ftg.getCountSkewLog2Rto();
					if (skewMed < minimumSkew && skewCount < minimumSkew) {
						skippedGenes.add(gene.getNames("_"));
						continue;
					}
				}
				
				//print out graphs
				printGraphs();
				
				//save gene and results
				scoredGenes.add(ftg);
				
				//cleanup
				treatmentGeneStats.nullBigArrays();
				controlGeneStats.nullBigArrays();
				
			}
			//writeOut unspliced sams
			int[] minMax = UCSCGeneLine.findMinMax(genes);
			writeOutUnspliceSams(treatmentSamReaders, treatmentMisSplicedSamOut, unSplicedTreatmentAlignments, minMax);
			writeOutUnspliceSams(controlSamReaders, controlMisSplicedSamOut, unSplicedControlAlignments, minMax);
			System.out.println();
		}
	}

	private void printGraphs() {
		File geneDir = new File (geneResultsDirectory, gene.getDisplayName());
		geneDir.mkdir();
		printSgrGraphs(geneDir);
	}
	
	public void printSgrGraphs(File saveDir){
		File f;
		String geneName = gene.getDisplayName();
		//treatment exonic bc
		float[] tbc = Num.intArrayToFloat(treatmentGeneStats.getBaseCoverage());
		f = new File (saveDir, geneName+"_TreatmentExonBC.sgr.gz");
		writeSgr(f, tbc, true);
		//control exonic bc
		float[] cbc = Num.intArrayToFloat(controlGeneStats.getBaseCoverage());
		f = new File (saveDir, geneName+"_ControlExonBC.sgr.gz");
		writeSgr(f, cbc, true);
		//log2Rto with introns
		f = new File (saveDir, geneName+"_ExonNormLog2Rto.sgr.gz");
		writeSgr(f, getScaledRatios(tbc, cbc), false);
		
		//print ggplot
		printExonicGGPlot(saveDir, tbc, cbc);
	}
	
	//calculates log2(numT+1/numC+1) for bases passing the minimumBaseReadCoverage threshold
	public float[] getScaledRatios(float[] tbc, float[] cbc) {
		//scalar to multiply C or divide T
		float scalar = (float)(treatmentGeneStats.getNumberExonicAlignments()+1) / (float)(controlGeneStats.getNumberExonicAlignments()+1);

		float[] ratios = new float[tbc.length];
		for (int i=0; i< ratios.length; i++){
			//fails base read coverage? if so the defaults to zero
			if ((tbc[i] + cbc[i]) < minimumBaseReadCoverage) continue;
			//divide T?
			if (tbc[i] >= cbc[i]) {
				float t = tbc[i]/ scalar;
				ratios[i] = Num.log2( (t+1.0f)/ (cbc[i]+1.0f));
			}
			else {
				float c = cbc[i] * scalar;
				ratios[i] = Num.log2( (tbc[i]+1.0f)/ (c+1.0f));
			}
		}
		return ratios;
	}

	private void writeSgr(File f, float[] vals, boolean makeRelative) {
		try {
			Gzipper out = new Gzipper(f);
			String chr = gene.getChrom() +"\t";
			float divider = 1;
			if (makeRelative) divider = Num.maxValue(vals);
			for (int i=0; i< vals.length; i++){
				if (vals[i] != 0) {
					out.println(chr+ (firstGeneBase+i)+ "\t"+ (vals[i]/divider));
				}
			}
			out.close();
		} catch (Exception e) {
			System.err.println("Problem writing "+f);
			e.printStackTrace();
		} 
	}
	
	public void printExonicGGPlot(File saveDir, float[] t, float[] c) {
		try {
			//generate table
			File table = new File (saveDir, gene.getDisplayName()+"_ggplot2Data.txt");
			table.deleteOnExit();
			PrintWriter out = new PrintWriter( new FileWriter (table));
			out.println("Group\tPosition\tRelativeCoverage");
			
			//find max
			float maxT = Num.findHighestFloat(t);
			float maxC = Num.findHighestFloat(c);
			//print
			for (int i=0; i< t.length; i++){
				if (t[i] == 0.0f) continue;
				float val = t[i]/maxT;
				out.println("T\t"+i+"\t"+val);
			}
			for (int i=0; i< c.length; i++){
				if (c[i] == 0.0f) continue;
				float val = c[i]/maxC;
				out.println("C\t"+i+"\t"+val);
			}
			//append R script
			File png = new File (saveDir, gene.getDisplayName()+"_Exonic.png");
			rScript.append("png('"+png+"', height=400, width=1200)\n");
			rScript.append("dfn = read.table(header=T, file='"+table+"')\n");
			rScript.append("ggplot(data=dfn, aes(x=Position, y=RelativeCoverage, group=Group, colour=Group)) + geom_line() + geom_point()\n");
			rScript.append("dev.off()\n");
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	private void writeOutUnspliceSams(SamReader[] samReaders, Gzipper out, HashSet<String> unsplicedNames, int[] startStop) {
		//fetch all chrom alignments from in start and stop
		try {
			//for each samReader
			for (SamReader sr: samReaders){
				SAMRecordIterator it = sr.queryOverlapping(chromosome, startStop[0], startStop[1]);
				while (it.hasNext()){
					SAMRecord sam = it.next();
					//make new String so can null sa
					String name = sam.getReadName();
					if (sam.getReadPairedFlag() == true && sam.getFirstOfPairFlag() == false) {
						name = name+"s";
					}					
					if (unsplicedNames.contains(name))  out.println(sam.getSAMString().trim());
				}
				it.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		} 
		//clear em!
		unsplicedNames.clear();
	}

	private boolean scoreGene() {

		//for treatment
		readCoverage = treatmentReadCoverage;
		unSplicedAlignments = unSplicedTreatmentAlignments;
		//count stat the gene
		treatmentGeneStats = countStatGene();		
		//calculate background
		treatmentGeneStats.calculateMedianBackground(fractionLengthBackground);
		//window scan 5' end
		treatmentGeneStats.windowScan(length5PrimeWindowScan);
		//count reads in background
		treatmentGeneStats.countBackgroundReads();
		//count reads in window
		treatmentGeneStats.countWindowReads(length5PrimeWindowScan);
		//check counts
		int total = treatmentGeneStats.getCountWindow() + treatmentGeneStats.getCountBackground();
		if (total < minimumWindowAndBackgroundCount) return false;
		//scan background
		treatmentGeneStats.windowScanBackground(length5PrimeWindowScan);

		//for control
		readCoverage = controlReadCoverage;
		unSplicedAlignments = unSplicedControlAlignments;
		//count stat the gene
		controlGeneStats = countStatGene();
		//calculate background
		controlGeneStats.calculateMedianBackground(fractionLengthBackground);
		//calculate median using best window coordinates from treatment 
		controlGeneStats.calculateMedianWindow(treatmentGeneStats.getMedianWindowIndex(), length5PrimeWindowScan);
		//count reads in background
		controlGeneStats.countBackgroundReads();
		//count reads in window
		controlGeneStats.countWindowReads(length5PrimeWindowScan);
		//check counts
		total = controlGeneStats.getCountWindow() + controlGeneStats.getCountBackground();
		if (total < minimumWindowAndBackgroundCount) return false;
		//scan background
		controlGeneStats.windowScanBackground(length5PrimeWindowScan);
		
		return true;
	}

		
	
	
	

	




	private TeleStats countStatGene() {
		ExonIntron[] exons = gene.getExons();
		int cDNALength = gene.getTotalExonicBasePairs();
		
		//containers for data
		HashSet<String> uniqueNames = new HashSet<String>();
		HashSet<String> unSpliced = new HashSet<String>();
		int[] baseCoverage = new int[cDNALength];
		ArrayList<String>[] baseCoverageNames = new ArrayList[cDNALength];
		
		//3'UTRs
		int utrStart = 0; 
		int utrEnd = 0;
		int lengthUtr = 0;
		int lengthNonUtr = 0;
		if (gene.isPlusStrand()){
			utrStart = gene.getCdsEnd() - firstGeneBase;
			utrEnd = gene.getTxEnd() - firstGeneBase;
		}
		else utrEnd = gene.getCdsStart()- firstGeneBase;	
		
		HashSet<String> unique3UTRNames = new HashSet<String>();
		
		int indexLast = exons.length-1;
		int index = 0;
		boolean examineStart;
		boolean examineEnd;	
		
		//for each exon
		for (int i=0; i< exons.length; i++){
			int start = exons[i].getStart() - firstGeneBase;
			int end = exons[i].getEnd() - firstGeneBase;
			
			//for each exonic base
			for (int j=start; j< end; j++) {
				if (readCoverage[j]!=null) {
					uniqueNames.addAll(readCoverage[j]);
					baseCoverage[index] = readCoverage[j].size();
					baseCoverageNames[index] = readCoverage[j];
					//check if in utr
					if (j >= utrStart && j < utrEnd) {
						unique3UTRNames.addAll(readCoverage[j]);
						lengthUtr++;
					}
					else lengthNonUtr++;
				} 
				index++;
			}
			
			//skip exons with less than 6 bp for diff splice
			if (exons[i].getLength() < 6) continue;
			
			//scan for unspliced
			examineStart = true;
			examineEnd = true;
			if (i==0) examineStart = false;
			else if (i==indexLast) examineEnd = false;
			if (examineStart){
				//copies
				if ((start-3) >= 0){ 
					System.arraycopy(readCoverage, start-3, toCheck, 0, 6);
					ArrayList<String> common = Misc.commonToAll(toCheck);
					unSpliced.addAll(common);	
				}
			}
			if (examineEnd){
				//copies
				//check end
				if ((end+3) < readCoverage.length) {
					System.arraycopy(readCoverage, end-3, toCheck, 0, 6);
					ArrayList<String> common = Misc.commonToAll(toCheck);
					unSpliced.addAll(common);
				}
			}
			
		}
		unSplicedAlignments.addAll(unSpliced);
		int numUniqueAlignments = uniqueNames.size();
		int numUnspliced = unSpliced.size();
		if (gene.isMinusStrand()) {
			Num.invertArray(baseCoverage);
			Misc.invertArray(baseCoverageNames);
		}

		return new TeleStats(lengthNonUtr, lengthUtr, numUniqueAlignments, unique3UTRNames.size(), numUnspliced, baseCoverage, baseCoverageNames);
		
	}

	private boolean loadTCSams() {
		//for treatment
		treatmentReadCoverage = new ArrayList[gene.getLength()];
		samReaders = treatmentSamReaders;
		readCoverage = treatmentReadCoverage;
		int numAlignmens = loadSams();
		if (numAlignmens < this.minimumGeneReadCoverage) return false;
		//for control
		controlReadCoverage = new ArrayList[treatmentReadCoverage.length];
		samReaders = controlSamReaders;
		readCoverage = controlReadCoverage;
		numAlignmens = loadSams();
		if (numAlignmens < this.minimumGeneReadCoverage) return false;
		return true;
	}


	private int loadSams() {
		int numReads = 0;
		boolean isNegStrand = gene.isMinusStrand();
		try {
			//for each samReader
			for (SamReader sr: samReaders){
				SAMRecordIterator it = sr.queryOverlapping(chromosome, firstGeneBase, lastGeneBase);
				while (it.hasNext()){
					SAMRecord sam = it.next();
					//check strand?
					if (isStranded){
						if (sam.getReadNegativeStrandFlag() != isNegStrand) continue;
					}
					SamAlignment sa = new SamAlignment(sam.getSAMString().trim(), true);
					//make new String so can null sa
					String name = new String(sa.getName());
					if (sa.isSecondPair()) name = name+"s";
					layoutSam(sa, name);
					sa = null;
					sam = null;
					numReads++;
				}
				it.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		} 
		return numReads;
	}

	private void layoutSam(SamAlignment sam, String name) {
		//for each cigar block in first, looking for MDNSH but not I
		Matcher mat = CIGAR.matcher(sam.getCigar());
		int index = sam.getPosition() - firstGeneBase;		
		int lastBase = lastGeneBase- firstGeneBase -1;
		while (mat.find()){
			String cCall = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (cCall.equals("M")) {
				//layout Ms
				for (int i = 0; i< numberBases; i++){
					//past end?
					if (index > lastBase) return;
					//past beginning?
					if (index >= 0) {
						//addit!
						ArrayList<String> names = readCoverage[index];
						if (names == null) {
							names = new ArrayList<String>();
							readCoverage[index] = names;
						}
						names.add(name);
					}
					index++;
				}
			}
			//N D (H,S, and I's won't match)
			else {
				for (int i = 0; i< numberBases; i++){
					//past end?
					if (index > lastBase) return;
					//advance 
					index++;
				}
			}
		}

	}


	/**Looks at each gene in a particular chromosome and sets the appropriate base type
	 * 0 integenic, 1 exonic, 2 intronic, maximizes exonic, minimizes intergenic.
	 *
	private void loadAnnotationBpArray() {
		int length = lastGeneBase - firstGeneBase;
		baseAnnotation = new byte[length]; //defaults to all 0 so intergenic
		for (UCSCGeneLine gene: genes){
			ExonIntron[] exons = gene.getExons();
			//load intronic if currently intergenic
			int start = exons[0].getEnd() - firstLastGeneBase[0];
			int end = exons[exons.length-1].getStart() - firstLastGeneBase[0];
			for (int i= start; i< end; i++) if (baseAnnotation[i] == 0) baseAnnotation[i] = 2;
			//load exonic
			for (ExonIntron e: exons){
				start = e.getStart()- firstLastGeneBase[0];
				end = e.getEnd()- firstLastGeneBase[0];
				for (int i= start; i< end; i++) baseAnnotation[i] = 1;
			}
		}
	}*/

	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader= new UCSCGeneModelTableReader(refSeqFile, 0);
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check gene name is unique
		if (reader.uniqueGeneNames() == false)  Misc.printExit("\nOne or more of your gene names is not unique.\n");
		chromGenes = reader.getChromSpecificGeneLines();
	}

	private void closeSams() {
		try {
			//readers
			for (SamReader sr: treatmentSamReaders) sr.close();
			for (SamReader sr: controlSamReaders) sr.close();
			//writers
			treatmentMisSplicedSamOut.close();
			controlMisSplicedSamOut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void openSams() {
		try {
			//Readers
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			treatmentSamReaders = new SamReader[treatmentBamFiles.length];
			controlSamReaders = new SamReader[controlBamFiles.length];
			//treatment
			System.out.println("Treatment:");
			for (int i=0; i< treatmentBamFiles.length; i++){
				System.out.println("\t"+treatmentBamFiles[i]);
				treatmentBamFiles[i].setLastModified(0);
				treatmentSamReaders[i] = factory.open(treatmentBamFiles[i]);
			}
			System.out.println("Control:");
			for (int i=0; i< controlBamFiles.length; i++){
				System.out.println("\t"+controlBamFiles[i]);
				controlBamFiles[i].setLastModified(0);
				controlSamReaders[i] = factory.open(controlBamFiles[i]);
			}

			//mis splice writers
			File t = new File (resultsDirectory, "treatmentMisSplice.sam.gz");
			t.deleteOnExit();
			treatmentMisSplicedSamOut = new Gzipper(t);
			File c = new File (resultsDirectory, "controlMisSplice.sam.gz");
			c.deleteOnExit();
			controlMisSplicedSamOut = new Gzipper(c);

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to open sam readers or writers.\n");
		} 
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Telescriptor(args);
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
					case 't': treatmentBamFiles = IO.extractFiles(new File(args[++i]), ".bam"); break;
					case 'c': controlBamFiles = IO.extractFiles(new File(args[++i]), ".bam"); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 's': resultsDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'i': isStranded = false; break;
					case 'g': minimumGeneReadCoverage = Integer.parseInt(args[++i]); break;
					case 'a': minimumWindowAndBackgroundCount = Integer.parseInt(args[++i]); break;
					case 'b': minimumBaseReadCoverage = Integer.parseInt(args[++i]); break;
					case 'l': minimumTranscriptLength = Integer.parseInt(args[++i]); break;
					case 'w': length5PrimeWindowScan = Integer.parseInt(args[++i]); break;
					case 'f': fractionLengthBackground = Float.parseFloat(args[++i]); break;
					case 'k': minimumSkew = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (treatmentBamFiles == null || treatmentBamFiles[0].canRead() == false) Misc.printErrAndExit("\nError: can't find your treatment bam alignment files?\n");
		if (controlBamFiles == null || controlBamFiles[0].canRead() == false) Misc.printErrAndExit("\nError: can't find your control bam alignment files?\n");
		if (refSeqFile == null || refSeqFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your refseq ucsc gene table files?\n");
		if (resultsDirectory == null) Misc.printErrAndExit("\nError: can't find your results directory?\n");
		resultsDirectory.mkdirs();
		geneResultsDirectory = new File(resultsDirectory, "Genes");
		geneResultsDirectory.mkdir();
		if (fullPathToR == null || fullPathToR.canExecute() == false) Misc.printErrAndExit("\nError: can't find or execute R? "+fullPathToR+"\n");
		
		//check for R and required libraries, don't need it if they just want the first and last 1/3 count table
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}

		String errors = IO.runRCommandLookForError("library(ggplot2);", fullPathToR, resultsDirectory);
		if (errors == null || errors.length() !=0){
			Misc.printErrAndExit("\nError: Cannot find the required R library.  Did you install ggplot2? R error message:\n\t\t"+errors+"\n\n");
		}



	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Telescriptor:  Sept 2014                            **\n" +
				"**************************************************************************************\n" +
				"Compares two RNASeq datasets for possible telescripting. Generates a spreadsheet of\n"+
				"statistics for each gene as well as a variety of graphs in exonic bp space. The\n"+
				"ordering of A and B is important since A is window scanned to identify the maximal 5'\n"+
				"region. Thus A should be where you suspect telescripting, B where you do not.\n"+

				"\nOptions:\n"+
				"-t Directory of bam files representing the first condition A.\n"+
				"-c Directory of bam files representing the second condition B.\n"+
				"-u UCSC refflat formatted Gene table. Run MergeUCSCGeneTable on a transcript table.\n"+
				"-s Director in which to save the results.\n"+
				"-r Full path to R, defaults to '/usr/bin/R', with installed ggplot2 package.\n"+

				"\nDefault Options:\n"+
				"-g Minimum gene alignment count, defaults to 50\n"+
				"-a Minimum window + background alignment count, defaults to 25\n"+
				"-k Minimum Log2(ASkew/BSkew), defaults to 2. Set to 0 to print all.\n"+
				"-b Minimum base read coverage for log2Ratio graph output, defaults to 10\n"+
				"-l Minimum transcript exonic length, defaults to 250\n"+
				"-w Size of 5' window for scanning, defaults to 125\n"+
				"-f Fraction of exonic gene length to calculate background, defaults to 0.5\n"+
				"-i Data is not stranded, assumes both first and second reads follow annotation.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/Telescriptor -u hg19EnsTrans.ucsc -t Bam/T\n"+
				"       -c Bam/C -s GV_MOR \n\n" +

				"**************************************************************************************\n");
	}


}
