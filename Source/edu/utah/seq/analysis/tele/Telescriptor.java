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
	private File fullPathToR = new File ("/usr/bin/R");
	private int minimumGeneReadCoverage = 50;
	private int minimumBaseReadCoverage = 10;
	private int minimumTranscriptLength = 250;
	private int length5PrimeWindowScan = 100;
	private float fractionLengthBackground = 0.5f;
	private float minIntronicToScoreMisSplice = 3;
	private double minimumTreatmentSkew = 2;

	//internal fields
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private HashMap<String,Integer> readNameAlignmentTreatment = new HashMap<String,Integer>();
	private HashMap<String,Integer> readNameAlignmentControl = new HashMap<String,Integer>();
	private HashMap<String,SamAlignment> treatmentSams = new HashMap<String,SamAlignment>();
	private HashMap<String,SamAlignment> controlSams = new HashMap<String,SamAlignment>();
	private HashSet<String> treatmentPrintedSamNames = new HashSet<String>();
	private HashSet<String> controlPrintedSamNames = new HashSet<String>();
	private ArrayList<String> skippedTranscripts = new ArrayList<String>();
	private HashMap<String,ArrayList<UCSCGeneLine>> genes = new HashMap<String,ArrayList<UCSCGeneLine>>();
	private SamReader[] treatmentSamReaders;
	private SamReader[] controlSamReaders;
	private Gzipper treatmentMisSplicedSamOut;
	private Gzipper controlMisSplicedSamsOut;
	private StringBuilder rScript = new StringBuilder("library(ggplot2)\n");
	public static final Pattern CIGAR = Pattern.compile("(\\d+)([MND])");
	public ArrayList<TeleGene> scoredGenes = new ArrayList<TeleGene>();
	private boolean debug = false;

	//working fields
	private String chromosome;
	private TeleTranscript rcst;
	private ExonIntron[] exons;
	private int firstBase;


	//constructor
	/**Stand alone.*/
	public Telescriptor(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load genes
		System.out.println("Loading transcript models...");
		loadGeneModels();

		//create sam readers and gzippers
		System.out.println("\nOpening SAM readers and writers...");
		openSams();
		writeOutSamHeaders();

		//walk through each chromosome of genes building ReadCoverageScoredTranscript
		System.out.println("\nExtracting alignments over each gene...");
		walkChromosomes();

		//close readers and writers
		closeSams();

		//print out skipped transcripts

		//analyze ReadCoverageScoredTranscripts
		System.out.println("\n# Genes to Analyze, "+scoredGenes.size());
		analyseGenes();
		
		//execute R graphing
		System.out.println("\nGenerating GGPlots...");
		runRGraphics();

		//sort and index misspliced reads
		System.out.println("\nSorting and writing mis-spliced sam alignments to file...");
		new PicardSortSam(treatmentMisSplicedSamOut.getGzipFile());
		new PicardSortSam(controlMisSplicedSamsOut.getGzipFile());

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
			controlMisSplicedSamsOut.println(head);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void analyseGenes() {
		try {
			//calculate skew stats and coverage
			for (TeleGene g : scoredGenes){
				ArrayList<TeleTranscript> ttAL = g.getScoredTranscripts();
				//for each transcript
				for (TeleTranscript tt: ttAL) tt.calculateSkewStats(fractionLengthBackground, length5PrimeWindowScan);
			}
			
			//run chisquare test on spiced vs misSplice for t and c
			System.out.println("\tRunning spliced vs mis-spliced chiSquare test...");
			calcDiffMisSplicePVal();
			
			//run chisquare test on spiced vs misSplice for t and c
			System.out.println("\tRunning read cound skew chiSquare test...");
			calcDiffSkewSplicePVal();
			
			//make writer for spreadsheet and write header
			PrintWriter spreadSheetOut = new PrintWriter( new FileWriter( new File (resultsDirectory, "teleStatsSummary.xls")));
			spreadSheetOut.println("GeneName\tTreatmentTypes\tControlTypes\tpAdjMisSplice\tlog2RtoMisSplice\tTranscriptName\tLength\t5' Index\tTreatment 5' Median\tTreatment 3' Median\tTreatment Median Log2Rto\tControl 5' Median\tControl 3' Median\tControl Median Log2Rto\tMedian Log2 (tSkew/ cSkew)\tTreatment 5' Count\tTreatment 3' Count\tTreatment Count Log2Rto\tControl 5' Count\tControl 3' Count\tControl Count Log2Rto\tCount Log2 (tSkew/ cSkew)\tpAdj Count Skew");
			
			File genesDir = new File (resultsDirectory, "Genes");
			genesDir.mkdirs();
			
			//print to file....
			for (TeleGene g : scoredGenes){
				ArrayList<TeleTranscript> ttAL = g.getScoredTranscripts();
				String geneInfo = g.toStringTabLine();
				//for each transcript
				for (TeleTranscript tt: ttAL){
					//print sgr graphs?
					double tSkew = tt.getTreatment().getLog2MedianSkew();
					if (tSkew > minimumTreatmentSkew  && tt.getControl().getLog2MedianSkew()< tSkew) {
						tt.printSgrGraphs(genesDir, minimumBaseReadCoverage);
					}
					tt.printExonicGGPlot(genesDir, rScript);
					//print summary line to spreadsheet
					String sumLine = tt.toStringTabLine(length5PrimeWindowScan, fractionLengthBackground);
					spreadSheetOut.println(geneInfo +"\t" +sumLine);
				}	
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
		//0 All Exon, 1 Splice Exon, 2 All Intron, 3 MisSpliced
		for (int i=0; i< numGenes; i++){
			TeleGene tg = scoredGenes.get(i);
			int[] t = tg.getMasterTypesTreatment();
			int[] c = tg.getMasterTypesControl();
			treatment[i] = new int[]{t[1], t[3]};
			control[i] = new int[]{c[1], c[3]};
			//calc log2Rto
			double lrto = Num.log2( ((double)(t[3]+1)/(double)(t[1]+2+t[3])) / ((double)(c[3]+1)/((double)c[1]+2+c[3])) );
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
	
	private void calcDiffSkewSplicePVal(){
		//count number of trans
		int numTrans = 0;
		for (int i=0; i< scoredGenes.size(); i++){
			numTrans += scoredGenes.get(i).getScoredTranscripts().size();
		}
				
		int[][] treatment = new int[numTrans][2];
		int[][] control = new int[numTrans][2];
		//0 All Exon, 1 Splice Exon, 2 All Intron, 3 MisSpliced
		int index = 0;
		for (int i=0; i< scoredGenes.size(); i++){
			ArrayList<TeleTranscript> al = scoredGenes.get(i).getScoredTranscripts();
			for (TeleTranscript tt: al){
				TeleStats ts = tt.getTreatment();
				treatment[index] = new int[]{ts.getFivePrimeAlignmentFragmentCount(), ts.getThreePrimeAlignmentFragmentCount()};
				ts = tt.getControl();
				control[index] = new int[]{ts.getFivePrimeAlignmentFragmentCount(), ts.getThreePrimeAlignmentFragmentCount()};
				index++;
			}
		}

		//estimate chi-square pvalues using R for resolution of extreemly small p-values, radiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, resultsDirectory, fullPathToR, true);

		//bonferroni correction
		double bc = Num.minus10log10(numTrans);

		//add back
		index = 0;
		for (int i=0; i< scoredGenes.size(); i++){
			ArrayList<TeleTranscript> al = scoredGenes.get(i).getScoredTranscripts();
			for (TeleTranscript tt: al){
				double pAdj = pVals[index++] + bc;
				if (pAdj > 0) {
					TeleStats ts = tt.getTreatment();
					ts.setpAdjSkewedReadCount((float)pAdj);
				}
			}
		}
	}

	private void walkChromosomes() {
		for (String chromName: chromGenes.keySet()){
			chromosome = chromName;
			UCSCGeneLine[] lines = chromGenes.get(chromosome);
			//group transcripts by gene
			groupTranscriptsByGene(lines);
			System.out.print("\t"+chromosome+"\t"+lines.length+"\t"+ genes.size()+" ");
			//clear sam names
			treatmentPrintedSamNames.clear();
			controlPrintedSamNames.clear();
			walkGenes();
			System.out.println();
		}
	}

	private void groupTranscriptsByGene(UCSCGeneLine[] lines) {
		genes.clear();
		ArrayList<UCSCGeneLine> al = null;
		for (UCSCGeneLine trans : lines){
			String geneName = trans.getDisplayName();
			if (genes.containsKey(geneName)){
				al = genes.get(geneName);
			}
			else {
				al= new ArrayList<UCSCGeneLine>();
				genes.put(geneName, al);
			}
			al.add(trans);
		}
	}

	private void walkGenes() {
		//for each gene group
		int counter = 0;
		for (ArrayList<UCSCGeneLine> transAL : genes.values()){
			if (counter++ > 100){
				counter = 0;
				System.out.print(".");
			}
			//System.out.println(transAL.get(0).getDisplayName());

			//check for at least one transcript >= min length and > 1 exon
			transAL = filterTranscripts(transAL);
			if (transAL.size() == 0) continue;

			//clear master types and alignments;
			readNameAlignmentTreatment.clear();
			readNameAlignmentControl.clear();
			treatmentSams.clear();
			controlSams.clear();

			ArrayList<TeleTranscript> scoredTranscripts = new ArrayList<TeleTranscript>();

			//for each transcript
			for (UCSCGeneLine transcript: transAL){
				exons = transcript.getExons();

				//fetch alignments over whole gene for treatment
				HashMap<String, SamAlignment> treatmentSams = fetchSams(transcript.getTxStart(), transcript.getTxEnd(), treatmentSamReaders);
				if (treatmentSams.size() < minimumGeneReadCoverage) {
					skippedTranscripts.add(transcript.getNames("_"));
					continue;
				}
				int numExonicT = countExonicAlignments(treatmentSams.values());
				if (numExonicT < minimumGeneReadCoverage) {
					skippedTranscripts.add(transcript.getNames("_"));
					continue;
				}
				//for control
				HashMap<String, SamAlignment> controlSams = fetchSams(transcript.getTxStart(), transcript.getTxEnd(), controlSamReaders);
				if (controlSams.size() < minimumGeneReadCoverage) {
					skippedTranscripts.add(transcript.getNames("_"));
					continue;
				}
				int numExonicC = countExonicAlignments(controlSams.values());
				if (numExonicC < minimumGeneReadCoverage) {
					skippedTranscripts.add(transcript.getNames("_"));
					continue;
				}
				//make new results container
				rcst = new TeleTranscript(transcript, numExonicT, numExonicC);

				//lay alignments out for read coverage
				layout(transcript.getTxStart(), transcript.getTxEnd(), treatmentSams, controlSams);

				//save
				scoredTranscripts.add(rcst);
			}
			//any transcripts? 
			if (scoredTranscripts.size() == 0) continue;

			//count master types
			int[] masterTypesTreatment = countMasterTypes(readNameAlignmentTreatment);
			int[] masterTypesControl = countMasterTypes(readNameAlignmentControl);

			TeleGene gene = new TeleGene(scoredTranscripts, masterTypesTreatment, masterTypesControl);
			scoredGenes.add(gene);

			//write out mis spliced sam alignments
			writeMisSplicedSams(readNameAlignmentTreatment, treatmentSams, treatmentPrintedSamNames, treatmentMisSplicedSamOut);
			writeMisSplicedSams(readNameAlignmentControl, controlSams, controlPrintedSamNames, controlMisSplicedSamsOut);
		}
	}

	private ArrayList<UCSCGeneLine> filterTranscripts(ArrayList<UCSCGeneLine> transAL) {
		ArrayList<UCSCGeneLine> good = new ArrayList<UCSCGeneLine>();
		for (UCSCGeneLine transcript: transAL){
			exons = transcript.getExons();
			//at least two exons?
			if (exons.length > 1 && transcript.getTotalExonicBasePairs() >= minimumTranscriptLength) good.add(transcript);
			else skippedTranscripts.add(transcript.getNames("_"));
		}
		return good;
	}

	private void writeMisSplicedSams( HashMap<String, Integer> nameType, HashMap<String, SamAlignment> nameSam, HashSet<String> printedNames, Gzipper out) {
		try {
			//for each nameType
			for (String name: nameType.keySet()){
				int type = nameType.get(name);
				//misSpliced?
				if (type == 3){
					//printed before?
					if (printedNames.contains(name) == false){
						out.println(nameSam.get(name));
					}
					printedNames.add(name);
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	private int[] countMasterTypes(HashMap<String, Integer> nameType) {
		//0 All Exon, 1 Splice Exon, 2 All Intron, 3 MisSpliced, 4 ambiguous
		int[] types = new int[5];
		for (Integer i : nameType.values()){
			types[i]++;
		}
		return types;
	}

	/**Returns the Alignment type:
	 * 0 All Exon
	 * 1 Splice Exon
	 * 2 All Intron
	 * 3 MisSpliced
	 * 4 Ambiguous (often exonic then spliced to a novel exon in intron)
	 * And -1 or the bpPosition of the misSpliced alignment.  Sometimes this is ambiguous and left at -1.*/
	private int[] mapMisSplices(float[] baseCoverage){
		float total = Num.sumArray(baseCoverage);
		if (debug) System.out.println("Total base coverage "+total);	
		float[] exonSum = new float[exons.length];
		float exonTotal = 0;
		int exonsWithCoverage = 0;
		int positionMisSplice = -1;
		int type = 4; //default to ambiguous

		//for each exon
		for (int i=0; i< exons.length; i++){
			if (debug) System.out.println("Exon "+exons[i]);
			int start = exons[i].getStart() - firstBase;
			if (start < 0) start = 0;
			int stop = exons[i].getEnd() - firstBase;
			if (stop > baseCoverage.length) stop = baseCoverage.length;
			//for each base in exon
			for (int j= start; j< stop; j++){
				exonSum[i] += baseCoverage[j];
			}
			exonTotal+= exonSum[i];
			if (exonSum[i] !=0) exonsWithCoverage++;
			if (debug) System.out.println(start+"\t"+stop+"\t"+exonSum[i]);
		}
		
		//score it: AllExon 0, SplicedExon 1, AllIntron 2, MisSpliced 3
		//must have >=3bp intronic for 3

		float numBpIntronic = total-exonTotal;
		float totalMin = total-minIntronicToScoreMisSplice;
		//entirely intronic?
		if (numBpIntronic >= totalMin){
			if (debug) System.out.println("Intronic");
			type = 2;
		}
		//entirely exonic
		else if (exonTotal >= totalMin){
			//spliced?
			if (exonsWithCoverage > 1){
				if (debug) System.out.println("Spliced");
				type = 1;
			}
			else {
				if (debug) System.out.println("Exonic");
				type = 0;
			}
		}
		//may be misSpliced
		else {
			if (debug) System.out.println("MisSpliced");
			
			//find it! doesn't always work, some are ambiguous
			int exonLengthMinOne = exons.length -1;
			for (int i=0; i< exons.length; i++){
				int start = exons[i].getStart() - firstBase;
				//if (start < 0) start = 0;
				int stop = exons[i].getEnd() - firstBase;
				//if (stop > baseCoverage.length) stop = baseCoverage.length;
				boolean look3 = true;
				boolean look5 = true;	
				//first exon?
				if (i == 0) look5 = false;
				//last exon?
				else if (i == exonLengthMinOne) look3 = false;
				//look 3' side
				if (look3 && baseCoverage[stop] !=0) {
					positionMisSplice = stop;
					type = 3;
					break;
				}
				//look 5' side
				if (look5 && baseCoverage[start-1] !=0) {
					positionMisSplice = start;
					type = 3;
					break;
				}
			}
			//defaults to type = 4
		}

		if (debug) {
			System.out.println(type +"\tType");
			if (type == 3){
				System.out.println("PosMisSplice "+positionMisSplice);
				System.exit(0);
			}
			//Misc.printArray(baseCoverage);
		}

		return new int[]{type, positionMisSplice};
	}


	private int countExonicAlignments(Collection<SamAlignment> values) {
		int hits = 0;
		for (SamAlignment sam: values){
			ArrayList<int[]> blocks = DefinedRegionDifferentialSeq.fetchAlignmentBlocks(sam.getCigar(), sam.getUnclippedStart());
			//for each block
			boolean hit = false;
			scan : {
				for (int[] ss : blocks){
					//for each exon
					for (ExonIntron e: exons){
						if (e.intersects(ss[0], ss[1])) {
							hit = true;
							break scan;
						}
					}
				}
			}
			if (hit) hits++;
		}
		return hits;
	}

	private void layout(int start, int stop, HashMap<String, SamAlignment> treatmentSams, HashMap<String, SamAlignment> controlSams) {
		//define the firstBase for the transcript
		firstBase = start;
		int len = stop-firstBase;
		//layout each set of alignments
		layout (true, treatmentSams, len);
		layout (false, controlSams, len);
	}

	private void layout(boolean isTreatment, HashMap<String, SamAlignment> sams, int lengthToLayout) {
		float[] baseCoverage = new float[lengthToLayout];
		ArrayList<String>[] baseCoverageNames = new ArrayList[lengthToLayout];
		int[] types = new int[5];
		ArrayList<Integer> positionsMisSplices = new ArrayList<Integer>();
		//for each alignment
		for (SamAlignment sa: sams.values()){
			//if (sa.getName().equals("HWI-ST179R:412:D1YHYACXX:2:2114:11181:70480")) debug = true;
			//else debug = false;

			float[] readBases = layout (sa, lengthToLayout);
			//add to total
			for (int i=0; i<lengthToLayout; i++) {
				baseCoverage[i]+= readBases[i];
				//add name, keep it to the fragment so will collapse with hashing
				if (readBases[i] !=0){
					if (baseCoverageNames[i] == null) baseCoverageNames[i] = new ArrayList<String>();
					baseCoverageNames[i].add(sa.getName());
				}
			}

			if (debug) {
				System.out.println(sa);
				System.out.println(rcst.getTranscript().getDisplayName()+"\t"+rcst.getTranscript().getName());
			} 
			int[] typePosMisMatch = mapMisSplices(readBases);
			types[typePosMisMatch[0]]++;
			if (typePosMisMatch[1] != -1) positionsMisSplices.add(new Integer(typePosMisMatch[1]));
			sa.setMisc(typePosMisMatch[0]);

		} 

		//for each alignment check type against master
		HashMap<String,Integer> nameType;
		HashMap<String,SamAlignment> nameAlignment;
		if (isTreatment) {
			nameType = readNameAlignmentTreatment;
			nameAlignment = treatmentSams;
		}
		else {
			nameType = readNameAlignmentControl;
			nameAlignment = controlSams;
		}
		for (SamAlignment sa: sams.values()){
			int currentType = sa.getMisc();
			//fetch name
			String alignmentName = sa.getName()+sa.isFirstPair();
			//misSpliced?
			if (currentType == 3) nameAlignment.put(alignmentName, sa);
			//new?
			if (nameType.containsKey(alignmentName) == false){
				nameType.put(alignmentName, new Integer(currentType));
			}
			else {
				//check if types differ, 0 All Exon, 1 Splice Exon, 2 All Intron, 3 MisSpliced
				int masterType = nameType.get(alignmentName);
				if (currentType != masterType){
					//skip if currentType is ambiguous
					if (currentType != 4){
						//reset if master is intronic or current is spliced exonic
						if (masterType == 2 || currentType == 1) nameType.put(alignmentName, new Integer(currentType));
						//reset if master is mis spliced and current is exonic or spliced exonic (not intronic or ambiguous)
						else if (masterType == 3 && currentType != 2) nameType.put(alignmentName, new Integer(currentType));
					}
					

				}
			}
		}
		//set results for this transcript 
		TeleStats rcs;
		if (isTreatment) rcs = rcst.getTreatment();
		else rcs = rcst.getControl();
		rcs.setBaseCoverage(baseCoverage);	
		rcs.setBaseCoverageNames(baseCoverageNames);
		rcs.setTypes(types);
		if (positionsMisSplices.size() !=0) rcs.setMisSplicePositions(Num.arrayListOfIntegerToInts(positionsMisSplices));
	}

	private float[] layout(SamAlignment sam, int lengthToLayout) {
		float[] bases = new float[lengthToLayout];
		//for each cigar block in first, looking for MDNSH but not I
		Matcher mat = CIGAR.matcher(sam.getCigar());
		int index = sam.getPosition() - firstBase;
		if (debug) System.out.println("Start "+sam.getPosition());		
		int lastBase = bases.length -1;
		while (mat.find()){
			String cCall = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (cCall.equals("M")) {
				//layout Ms
				for (int i = 0; i< numberBases; i++){
					//past end?
					if (index > lastBase) return bases;
					//past beginning?
					if (index >= 0) {
						//addit!
						bases[index]++;
						if (debug) System.out.println(sam.getReferenceSequence()+"\t"+(index+firstBase));
					}
					index++;
				}
			}
			//N D (H,S, and I's won't match)
			else {
				for (int i = 0; i< numberBases; i++){
					//past end?
					if (index > lastBase) return bases;
					//advance 
					index++;
				}
			}
		}
		return bases;
	}

	private HashMap<String, SamAlignment> fetchSams(int start, int stop, SamReader[] samReaders) {
		HashMap<String, SamAlignment> sams = new HashMap<String, SamAlignment>();
		try {
			//for each samReader
			for (SamReader sr: samReaders){
				SAMRecordIterator it = sr.queryOverlapping(chromosome, start, stop);
				while (it.hasNext()){
					SAMRecord sam = it.next();
					SamAlignment sa = new SamAlignment(sam.getSAMString().trim(), true);
					String name = sa.getName();
					if (sa.isSecondPair()) name = name+"s";
					sams.put(name, sa);
				}
				it.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		} 
		return sams;
	}

	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader= new UCSCGeneModelTableReader(refSeqFile, 0);
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check gene name is unique
		chromGenes = reader.getChromSpecificGeneLines();
	}

	private void closeSams() {
		try {
			//readers
			for (SamReader sr: treatmentSamReaders) sr.close();
			for (SamReader sr: controlSamReaders) sr.close();
			//writers
			treatmentMisSplicedSamOut.close();
			controlMisSplicedSamsOut.close();
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
			controlMisSplicedSamsOut = new Gzipper(c);

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
					case 'g': minimumGeneReadCoverage = Integer.parseInt(args[++i]); break;
					case 'b': minimumBaseReadCoverage = Integer.parseInt(args[++i]); break;
					case 'l': minimumTranscriptLength = Integer.parseInt(args[++i]); break;
					case 'w': length5PrimeWindowScan = Integer.parseInt(args[++i]); break;
					case 'f': fractionLengthBackground = Float.parseFloat(args[++i]); break;
					case 'i': minIntronicToScoreMisSplice = Float.parseFloat(args[++i]); break;
					case 'j': minimumTreatmentSkew = Double.parseDouble(args[++i]); break;
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
		if (fullPathToR == null || fullPathToR.canExecute() == false) Misc.printErrAndExit("\nError: can't find or execute R? "+fullPathToR+"\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Telescriptor:  August 2014                          **\n" +
				"**************************************************************************************\n" +
				"Compares two RNASeq datasets for possible telescripting. Generates a spreadsheet of\n"+
				"statistics for each transcript as well as a variety of graphs in genomic and exonic\n"+
				"bp space. The ordering of T and C is important since T is window scanned to identify\n"+
				"the maximal 5' region. Thus T should be where you suspect telescripting, C where you\n"+
				"do not.\n"+

				"\nOptions:\n"+
				"-t Directory of bam files representing the first condition.\n"+
				"-c Directory of bam files representing the second condition.\n"+
				"-u UCSC refflat formatted transcript table. First column is gene name, second column\n"+
				"       transcript name.\n"+
				"-s Director in which to save the results.\n"+
				"-r Full path to R, defaults to '/usr/bin/R'\n"+
				
				"\nDefault Options:\n"+
				"-g Minimum gene alignment count, defaults to 50\n"+
				"-b Minimum base read coverage for log2Ratio graph output, defaults to 10\n"+
				"-l Minimum transcript exonic length, defaults to 250\n"+
				"-w Size of 5' window for scanning, defaults to 100\n"+
				"-f Fraction of exonic gene length to calculate background, defaults to 0.5\n"+
				"-i Minimum # intronic bases to the score alignment as mis spliced, defaults to 3\n"+
				"-j Minimum log2(5'TBaseCov/3'TBaseCov) treatment skew to write out sgr graphs,\n"+
				"        defaults to 2.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/Telescriptor -u hg19EnsTrans.ucsc -t Bam/T\n"+
				"       -c Bam/C -s GV_MOR \n\n" +

				"**************************************************************************************\n");
	}
	

}
