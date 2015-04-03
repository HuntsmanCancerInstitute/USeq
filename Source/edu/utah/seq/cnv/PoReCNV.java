package edu.utah.seq.cnv;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.analysis.DefinedRegionDifferentialSeq;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.analysis.multi.GeneCount;
import edu.utah.seq.analysis.multi.Replica;
import edu.utah.seq.data.HeatMapMaker;
import edu.utah.seq.data.HeatMapMakerPosNeg;
import edu.utah.seq.data.Info;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.apps.Bar2USeq;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;


/** App to detect CNV's in normal datasets.  Wraps Alun Thomas's algorithm.
 * @author Nix
 * */
public class PoReCNV {

	//user defined fields
	private File saveDirectory;
	private File bamDirectory;
	private File vcfDirectory;
	private String[] caseSampleFileNames;
	private String[] controlSampleFileNames;
	private File fullPathToR = new File ("/usr/bin/R");
	private File refSeqFile;
	private int minimumCounts = 20;
	private int maxAlignmentsDepth = 50000;
	private File alunRScript = null;
	private boolean deleteTempFiles = true;
	private String genomeVersion;
	private float minimumLog2Ratio = 0.585f;
	private float minimumAdjPVal = 0.01f;
	private int minNumInEachChunk = 1500;
	private int numberConcurrentThreads = 0;
	private boolean mergeExonCounts = false;
	private String annoType = "Exon";
	private boolean excludeSexChromosomes = true;
	private int minimumGenotypeQuality = 20;

	//internal fields
	private File graphDirectory;
	private UCSCGeneLine[] genes;
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private HashMap<String, UCSCGeneLine> name2Gene;
	public static final Pattern BAD_NAME = Pattern.compile("(.+)/[12]$");
	private HashSet<String> flaggedGeneNames = new HashSet<String>();
	private LinkedHashSet<String> geneNamesWithMinimumCounts = new LinkedHashSet<String>();
	private Replica[] cases = null;
	private Replica[] controls = null;
	private PoReDataChunk[] chunkThreads;

	//container for cases and their exons that pass the minimumCounts threshold
	private GeneExonSample[] ges = null;
	private File countTableFile;
	private int numExonsProcessed = 0;
	private int numberPassingExons = 0;
	private int[] passingExonsByCase = null;

	//from R
	private float significanceThreshold;

	//for loading data
	private int workingFragmentNameIndexPlus = 1;
	private int workingFragmentNameIndexMinus = -1;
	private HashMap<String, Integer> workingFragNameIndex = new HashMap<String, Integer>(10000);


	//constructors
	/**Stand alone.*/
	public PoReCNV(String[] args){	
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);
		//launch
		run();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	//methods
	public void run(){

		//load gene models
		System.out.println("Loading regions/ gene models...");
		loadGeneModels();		

		//load samples 
		loadSamples();

		//write out count table of genes passing minimum counts
		if (controls != null) writeCountTableMatrix();
		else {
			makeGeneExonSamples();

			chunkGeneExonSamples();

			System.out.println("\nProcessing data chunks...");
			boolean ok = runThreads();
			if (ok == false) {
				System.err.println("\nError processing chunks, check logs and temp files. Aborting.\n");
				return;
			}
			//set sig threshold
			setSignificanceThreshold();

			//write out graph files for each sample
			System.out.println("\nSaving graphs and spreadsheets...");
			exportGeneExonSpreadSheets();
			exportCaseSampleGraphs();
			
			System.out.println("\n"+numExonsProcessed+" "+annoType+"s processed, "+numberPassingExons+" passed thresholds ("+minimumLog2Ratio+" Lg2Rto, "+minimumAdjPVal+" AdjPVal)\n");
		}

		//convert graphs to useq format
		new Bar2USeq(graphDirectory, true);
	}

	private void setSignificanceThreshold() {
		float[] sigs = new float[chunkThreads.length];
		for (int i=0; i< sigs.length; i++) sigs[i] = chunkThreads[i].getSignificanceThreshold();
		//System.out.println("\nSignificance thresholds for each data chunk:\n"+ Misc.floatArrayToString(sigs, ", "));
		this.significanceThreshold = Num.mean(sigs);

	}


	/**Runs each in its own thread to max threads defined by user. Returns whether all completed without issue.*/
	private boolean runThreads() {
		ArrayList<PoReDataChunk> waitingJobs = new ArrayList<PoReDataChunk>();
		String old = "";
		while (true){ 
			try {
				//check status
				int waiting = 0;
				int running = 0;
				int complete = 0;
				int error = 0;
				waitingJobs.clear();
				for (PoReDataChunk c: chunkThreads) {
					int status = c.getStatus();
					if (status == 0) waitingJobs.add(c);
					else if (status == 1) running++;
					else if (status == 2) complete++;
					else error++;
				}
				waiting = waitingJobs.size();
				String currentStatus = "\tW:"+waiting+" R:"+running+" C:"+complete+" E:"+error;
				if (currentStatus.equals(old) == false) {
					System.out.println(currentStatus);
					old = currentStatus;
				}

				//all complete?
				if (waiting == 0 && running == 0) return (error == 0);
				//start jobs?
				if (waiting !=0){
					int toStart = numberConcurrentThreads - running;
					if (toStart > waiting) toStart = waiting;
					if (toStart >0){
						for (int i=0; i< toStart; i++) {
							PoReDataChunk c = waitingJobs.get(i);
							if (c.isAlive() == false && c.isInterrupted() == false) c.start();
							else System.out.println("\nWARNING, already started ("+c.isAlive()+
									") or interrupted ("+c.isInterrupted()+"). Skipping -> "+c.getRealName()+"\n");
						}
					}
				}
				
				//sleep 2 seconds
				Thread.sleep(4000);
			} catch (InterruptedException ie) {
				ie.printStackTrace();
				return false;
			}
		}
	}

	private void exportGeneExonSpreadSheets() {
		try {
			String currGeneName = ges[0].getGene().getDisplayName();
			ArrayList<GeneExonSample> al = new ArrayList<GeneExonSample>();
			al.add(ges[0]);

			Gzipper outAll = new Gzipper(new File(saveDirectory, "results"+annoType+"All.xls.gz"));
			Gzipper outPass = new Gzipper(new File(saveDirectory, "results"+annoType+"Pass.xls.gz"));

			//print headers
			String gn = "Gene Name";
			if (ges[0].getGene().getName() != null) gn = gn+"\tDescription";
			StringBuilder sb = new StringBuilder(gn+"\tCoordinates\tExon Index\tIGB Link\tIGV Link\tPassing Samples");
			for (int i=0; i< cases.length; i++){
				String name = cases[i].getNameNumber();
				sb.append("\t"+name+" Lg2(Ob/Ex)\t"+name+" Res\t"+name+" Counts");
			}
			sb.append("\tLog2Rto "+minimumLog2Ratio);
			sb.append("\tResidual "+minimumAdjPVal+" sig threshold +/- "+significanceThreshold);

			outAll.println(sb);
			outPass.println(sb);

			for (int i=1; i< ges.length; i++){
				String nextName = ges[i].getGene().getDisplayNameThenName();
				if (nextName.equals(currGeneName) == false){
					printGeneSamples(outAll, outPass, al);
					currGeneName = nextName;
					al.clear();
				}
				al.add(ges[i]);
			}
			//print out last!
			printGeneSamples(outAll, outPass, al);
			
			outAll.close();
			outPass.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nError printing spreadsheets.\n");
		}
	}

	private void printGeneSamples(Gzipper outAll, Gzipper outPass, ArrayList<GeneExonSample> al) throws IOException {
		//print geneName, coordinates, exon index, igb, igv, passsingSampleNames, sample1_Lg2Rt, sample1_Res, sample1_counts, sample2_.....
		GeneExonSample g = al.get(0);
		int exonIndex = g.getExonIndex();

		StringBuilder sb = new StringBuilder (fetchGeneInfo(g));
		//any samples that pass?
		String passingSamples = exonPassThresholds(exonIndex, al);
		if (passingSamples != null) sb.append (passingSamples+"\t");
		else sb.append("\t");

		//add first sample
		sb.append(g.getDataString());

		//for each subsequent sample
		for (int i=1; i< al.size(); i++){
			GeneExonSample t = al.get(i);

			//different exon index?
			if (t.getExonIndex() != exonIndex){
				//print out old
				outAll.println(sb);
				if (passingSamples != null) outPass.println(sb);

				//reset for new
				exonIndex = t.getExonIndex();
				sb = new StringBuilder (fetchGeneInfo(t));
				passingSamples = exonPassThresholds(exonIndex, al);
				if (passingSamples != null) sb.append (passingSamples);
			}

			sb.append("\t");
			sb.append(t.getDataString());
		}

		//print last
		outAll.println(sb);
		if (passingSamples != null) outPass.println(sb);
	}

	/**Returns null if no sample exon passes filters, otherwise comma delimited sample name list of those that pass.*/
	public String exonPassThresholds(int exonIndex, ArrayList<GeneExonSample> al){
		ArrayList<String> sampPass = new ArrayList<String>();
		for (int i=0; i<al.size(); i++){
			GeneExonSample g = al.get(i);
			
			//correct index?
			if (g.getExonIndex() == exonIndex){
				if (passesThresholds(g)) {
					sampPass.add(cases[g.getSampleIndex()].getNameNumber());
					passingExonsByCase[g.getSampleIndex()]++;
				}
			}
		}
		if (sampPass.size() !=0) {
			numberPassingExons++;
			return Misc.stringArrayListToString(sampPass, ",");
		}
		return null;
	}

	public String fetchGeneInfo(GeneExonSample g){
		UCSCGeneLine gene = g.getGene();
		ExonIntron[] exons = gene.getExons();
		String chr = gene.getChrom();
		int exonIndex = g.getExonIndex();
		//name
		StringBuilder sb = new StringBuilder (g.getGene().getDisplayName()); sb.append("\t");
		//description
		if (g.getGene().getName() != null) sb.append(g.getGene().getName()); sb.append("\t");
		//exon coordinates
		sb.append(chr);
		sb.append(":");
		int start;
		int stop;
		if (mergeExonCounts) {
			start = gene.getTxStart();
			stop = gene.getTxEnd();
		}
		else {
			start = exons[exonIndex].getStart();
			stop = exons[exonIndex].getEnd();
		}
		sb.append(start);
		sb.append("-");
		sb.append(stop);
		sb.append("\t");
		//exon index
		sb.append(exonIndex); sb.append("\t");
		//igb link
		sb.append(fetchIGBLink(chr, start, stop)); sb.append("\t");
		//igv link
		sb.append(fetchIGVLink(chr, start, stop)); sb.append("\t");

		return sb.toString();
	}

	public static String fetchIGVLink(String chr, int begin, int stop){
		int start = begin-5000;
		if (start < 0) start = 0;
		int end = stop+5000;
		return "=HYPERLINK(\"http://localhost:60151/goto?locus="+chr+":"+start+ "-" + end+"\",\"IGV\")";
	}

	public String fetchIGBLink(String chr, int begin, int stop){
		int start = begin-5000;
		if (start < 0) start = 0;
		int end = stop+5000;
		return "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid="+chr+"&start="+start+"&end="+end+ "\",\"IGB\")";
	}

	private boolean passesThresholds(GeneExonSample g) {
		if (Math.abs(g.getObsExpLgRto()) < minimumLog2Ratio || Math.abs(g.getResidual()) < significanceThreshold ) return false;
		return true;
	}
	
	private void searchForHetSNPs(){
		//split ges by sample
		GeneExonSample[][] sampleGES = new GeneExonSample[cases.length][numExonsProcessed];
		for (int i=0; i< ges.length; i++){
				sampleGES[ges[i].getSampleIndex()][ges[i].getGlobalExonIndex()] = ges[i];
		}
		//for each sample
		for (int i=0; i< sampleGES.length; i++){
			
			//make a VCF reader
			String name = Misc.removeExtension(caseSampleFileNames[i]) + ".vcf.gz";
			File vcf = new File(vcfDirectory, name);
			VCFFileReader vr = new VCFFileReader(vcf);
			
			GeneExonSample[] sGes = sampleGES[i];
			//for each exon
			for (int j=0; j< sGes.length; j++){
				//only do for deletions
				if (sGes[j].getResidual() > 0) continue;
				
				//define start and stop to search
				int start;
				int stop;
				UCSCGeneLine gene = sGes[j].getGene();
				if (mergeExonCounts){
					start = gene.getTxStart();
					stop = gene.getTxEnd();
				}
				else {
					ExonIntron exon = gene.getExons()[sGes[j].getExonIndex()];
					start = exon.getStart();
					stop = exon.getEnd();
				}
				
				//fetch vcf records
				CloseableIterator<VariantContext> vcfs = vr.query(ges[j].getGene().getChrom(), start, stop);
				StringBuilder sb = new StringBuilder();
				while (vcfs.hasNext()){
					VariantContext vc = vcfs.next();
					if (vc.isFiltered()) continue;
					int position = vc.getStart();
					GenotypesContext gc = vc.getGenotypes();
					Iterator<Genotype> gci = gc.iterator();
					while (gci.hasNext()){
						Genotype g = gci.next();
						if (g.hasGQ() && g.getGQ() < minimumGenotypeQuality) continue;
						sb.append(g.getGenotypeString());
						sb.append("_");
						sb.append(position);
						sb.append(";");
					}
				}
				if (sb.length() !=0) {
					sGes[j].setGenotypePosition(sb.toString().substring(0, sb.length()-1));
System.out.println(sGes[j].toString()+"\t"+sGes[j].getGenotypePosition());
				}
			}
			vr.close();
		}
	}

	private void exportCaseSampleGraphs() {
		//split ges by sample
		GeneExonSample[][] sampleGES = new GeneExonSample[cases.length][numExonsProcessed];
		for (int i=0; i< ges.length; i++){
			sampleGES[ges[i].getSampleIndex()][ges[i].getGlobalExonIndex()] = ges[i];
		}
		//make dirs
		graphDirectory = new File (saveDirectory, "DataTracks");
		graphDirectory.mkdir();
		File resDir = new File (graphDirectory, annoType+"Residual");
		resDir.mkdir();
		File rtoDir = new File (graphDirectory, annoType+"ObsExpLog2Rto");
		rtoDir.mkdir();
		
		//double[] cv = new double[cases.length];
		
		//for each write files
		System.out.println("\nResidual stats:\nDataset\t"+annoType+"sPass\tMean\tMedian\tStdDev\tMin\tMax\t10th\t90th");
		for (int i=0; i< sampleGES.length; i++){
			//make sample specific dirs
			File resSampDir = new File (resDir, "res_"+cases[i].getNameNumber());
			resSampDir.mkdir();
			resSampDir.deleteOnExit();
			File rtoSampDir = new File (rtoDir, "oe_"+cases[i].getNameNumber());
			rtoSampDir.mkdir();
			rtoSampDir.deleteOnExit();
			saveGraphs(cases[i].getNameNumber(), sampleGES[i], resSampDir, rtoSampDir);
			//calc residual stats
			float[] res = getResiduals(sampleGES[i]);
			Arrays.sort(res);
			String stats = Num.statFloatArray(res);
			System.out.println(caseSampleFileNames[i]+ "\t"+passingExonsByCase[i]+ "\t"+ stats);
		}
		
		//double meanCV = Num.mean(cv);
		//System.out.println("\nCoefficient of variation for sample residuals:\n\tMean\t"+meanCV);
		//for (int i=0; i< sampleGES.length; i++){
			//System.out.println(caseSampleFileNames[i]+"\t"+cv[i]);
		//}
		
	}

	private double calculateResidualCV(GeneExonSample[] ges) {
		float[] residuals = new float[ges.length];
		for (int i=0; i< ges.length; i++) residuals[i] = ges[i].getResidual();
		double mean = Num.mean(residuals);
		double std = Num.standardDeviation(residuals, mean);
		return std/mean;
	}
	
	private float[] getResiduals(GeneExonSample[] ges) {
		float[] residuals = new float[ges.length];
		for (int i=0; i< ges.length; i++) residuals[i] = ges[i].getResidual();
		return residuals;
	}

	private void saveGraphs(String sampleName, GeneExonSample[] ges, File resDir, File rtoDir) {
		//walk through data and make a SmoothingWindow[] for each chrom
		String currChr = ges[0].getGene().getChrom();
		ArrayList<SmoothingWindow> smAL = new ArrayList<SmoothingWindow>();

		for (int i=0; i< ges.length; i++){
			//make an sm
			int start;
			int stop;
			UCSCGeneLine gene = ges[i].getGene();
			if (mergeExonCounts){
				start = gene.getTxStart();
				stop = gene.getTxEnd();
			}
			else {
				ExonIntron exon = gene.getExons()[ges[i].getExonIndex()];
				start = exon.getStart();
				stop = exon.getEnd();
			}

			SmoothingWindow sm = new SmoothingWindow(start, stop, new float[]{ges[i].getResidual(), ges[i].getObsExpLgRto()});

			//diff chrom?
			if (ges[i].getGene().getChrom().equals(currChr) == false) {
				//diff chrom, write out old
				Info info = new Info(sampleName, genomeVersion, currChr, ".", 0, null);
				SmoothingWindow[] sms = new SmoothingWindow[smAL.size()];
				smAL.toArray(sms);
				saveStairStepGraph(0, sms, info, resDir, "#FF0000", true); //red
				saveStairStepGraph(1, sms, info, rtoDir, "#FFFF00", true); //yellow

				//reset and set new
				currChr = ges[i].getGene().getChrom();
				smAL.clear();
			}

			smAL.add(sm);
		}
		//process last!
		Info info = new Info(sampleName, genomeVersion, currChr, ".", 0, null);
		SmoothingWindow[] sms = new SmoothingWindow[smAL.size()];
		smAL.toArray(sms);		
		saveStairStepGraph(0, sms, info, resDir, "#FF0000", true); //red
		saveStairStepGraph(1, sms, info, rtoDir, "#FFFF00", true); //yellow
	}

	/**Saves bar heatmap/ stairstep graph files*/
	public void saveStairStepGraph (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color, boolean posNeg){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//save in info
		info.setNotes(map);
		//get heatmap positions and values
		PointData pd;		
		if (posNeg){
			HeatMapMakerPosNeg hm = new HeatMapMakerPosNeg(scoreIndex, 0, 0);
			pd = hm.makeHeatMapPositionValues(sm);
		}
		else {
			HeatMapMaker hm = new HeatMapMaker(scoreIndex,0);
			pd = hm.makeHeatMapPositionValues(sm, false);
		}
		pd.setInfo(info);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
	}

	private void loadSamples(){

		//load cases
		System.out.println("\nLoading case samples...");
		//serialized version present?
		File serCases = new File (saveDirectory, "cases.ser");
		File serGenes = new File (saveDirectory, "genes.ser");
		if (serCases.exists() && serGenes.exists()) {
			cases = (Replica[]) IO.fetchObject(serCases);
			System.out.println("\tWARNING: loading case count data from prior run.  Delete "+serCases+" and restart to load from alignment files.");
			if (cases.length != caseSampleFileNames.length) Misc.printErrAndExit("\nError: the number of cases and case file names differ.  Delete "+serCases+" and restart.\n");
			geneNamesWithMinimumCounts = (LinkedHashSet<String>) IO.fetchObject(serGenes);
		}
		else {
			cases = new Replica[caseSampleFileNames.length];
			//String[] shortNames = new String[cases.length];
			for (int i=0; i< caseSampleFileNames.length; i++){
				File bam = new File(bamDirectory, caseSampleFileNames[i]);
				String name = Misc.removeExtension(caseSampleFileNames[i]);
				//shortNames[i] = name;
				cases[i] = new Replica(name, bam);
				System.out.print("\t"+caseSampleFileNames[i]);
				loadReplica(cases[i]);
				System.out.println("\t"+cases[i].getTotalCounts());
			}
			//shortNames = Misc.trimCommon(shortNames);
			//for (int i=0; i< cases.length; i++) cases[i].setNameNumber(shortNames[i]);
			IO.saveObject(serCases, cases);
			IO.saveObject(serGenes, geneNamesWithMinimumCounts);
		}

		//load controls
		if (controlSampleFileNames != null){
			System.out.println("Loading control samples...");
			String[]  shortNames = new String[controls.length];
			controls = new Replica[this.controlSampleFileNames.length];
			for (int i=0; i< controlSampleFileNames.length; i++){
				File bam = new File(bamDirectory, controlSampleFileNames[i]);
				String name = Misc.removeExtension(controlSampleFileNames[i]);
				shortNames[i] = name;
				controls[i] = new Replica(name, bam);
				System.out.print("\t"+controlSampleFileNames[i]);
				loadReplica(controls[i]);
				System.out.println("\t"+controls[i].getTotalCounts());
			}
			shortNames = Misc.trimCommon(shortNames);
			for (int i=0; i< cases.length; i++) controls[i].setNameNumber(shortNames[i]);
		}

		//any genes with too many reads that were excluded?
		if (flaggedGeneNames.size() !=0) {
			String[] badGeneNames = Misc.hashSetToStringArray(flaggedGeneNames);			
			for (Replica replica: cases) replica.removeFlaggedGenes(badGeneNames);
			if (controlSampleFileNames != null) for (Replica replica: controls) replica.removeFlaggedGenes(badGeneNames);
			//remove them from geneNamesWithMinimumCounts
			for (String baddie: badGeneNames) geneNamesWithMinimumCounts.remove(baddie);
		}

		System.out.println();
		System.out.println(geneNamesWithMinimumCounts.size()+" genes, with >= "+minimumCounts+" and < "+maxAlignmentsDepth+" counts, will be examined for copy number alterations.\n ");
		System.out.println("Skipped genes: "+flaggedGeneNames);
	}

	public void loadBlocks(SAMRecord sam, int chrStartBp, ArrayList<Integer>[] bpNames, boolean[] badBases){

		NameInteger nameIndex = null;
		boolean addIt;
		ArrayList<int[]> blocks = DefinedRegionDifferentialSeq.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
		//add name to each bp
		for (int[] b : blocks){
			int start = b[0] - chrStartBp;
			int stop = b[1] - chrStartBp;
			//need to watch for out of bounds issues, sometimes the length of the chromosome is incorrect in the bam header.
			if (stop > badBases.length) stop = badBases.length;
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

	/**Fetches an old or makes a new Integer to represent the sam read name (e.g. fragment name)*/
	public NameInteger fetchFragmentNameIndex(SAMRecord sam){
		String samReadName = sam.getReadName();
		Integer index;
		Matcher mat = BAD_NAME.matcher(samReadName);
		if (mat.matches()) samReadName = mat.group(1);
		if (workingFragNameIndex.containsKey(samReadName)) index = workingFragNameIndex.get(samReadName);
		else {
			if (sam.getReadNegativeStrandFlag()) {
				index = new Integer(workingFragmentNameIndexMinus--);
			} else {
				index = new Integer(workingFragmentNameIndexPlus++);
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
		for (SAMSequenceRecord sr: seqs) chromLength.put(sr.getSequenceName(), sr.getSequenceLength()+1000);

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

			//unaligned? 
			if (sam.getReadUnmappedFlag()) continue;
			chrom = sam.getReferenceName();

			if (chromGenes.containsKey(chrom) == false) {
				continue;
			}		
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

			//unaligned? 
			if (sam.getReadUnmappedFlag()) continue;

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

		if (chromGenes.containsKey(chrom)) loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);

		reader.close();
		bpNames = null;
		iterator = null;
		seqs = null;
		chromLength = null;
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
						if (bpNames[y] != null) exonReads.addAll(bpNames[y]);
					}
					exonCounts[x] = exonReads.size();
					allReads.addAll(exonReads);
				}
			int numCounts = allReads.size();
			if (numCounts !=0){
				GeneCount tcg = new GeneCount(numCounts, exonCounts);
				geneCounts.put(geneName, tcg);
				replica.setTotalCounts(replica.getTotalCounts() + numCounts);
				if (numCounts >= minimumCounts && geneNamesWithMinimumCounts.contains(geneName) == false) geneNamesWithMinimumCounts.add(geneName);
			}
		}	
		//clean up
		allReads = null;
		exonReads = null;
	}

	/**Just for cases*/
	public void buildCaseGeneExonSamples(){
		try {
			//start counter for all exons to match Alun's code, starts with 1
			int globalExonIndex = 1;
			ArrayList<GeneExonSample> gesAL = new ArrayList<GeneExonSample>();

			//for each gene, must walk through so that these are sorted by chrom
			for (UCSCGeneLine gene: genes){
				//one to analyze?
				String geneName = gene.getDisplayNameThenName();
				if (geneNamesWithMinimumCounts.contains(geneName) == false) continue;
				int numExons = gene.getExons().length;

				//collect exon counts int[sample][exon]
				int[][] caseExonCounts = fetchExonCounts(geneName, numExons, cases);

				//for each exon
				for (int i=0; i< numExons; i++){
					//count total for this exon across all the cases, skip if too few
					int total = 0;
					for (int j=0; j<caseExonCounts.length; j++) total += caseExonCounts[j][i];
					if (total < minimumCounts) continue;
					numExonsProcessed++;

					//for each sample
					for (int j=0; j<caseExonCounts.length; j++){
						GeneExonSample g = new GeneExonSample(gene, globalExonIndex-1, (short)i, (short)j, caseExonCounts[j][i]);
						gesAL.add(g);
					}
					globalExonIndex++;
				}
			}
			//save GES
			ges = new GeneExonSample[gesAL.size()];
			gesAL.toArray(ges);
		}
		catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing out count table.");
		}
	}

	public void chunkGeneExonSamples(){
		//group by global exon index
		ArrayList<GeneExonSample>[] byGlobal = new ArrayList[numExonsProcessed];
		for (int i=0; i< numExonsProcessed; i++) byGlobal[i] = new ArrayList<GeneExonSample>();
		for (int i=0; i< ges.length; i++) byGlobal[ges[i].getGlobalExonIndex()].add(ges[i]);

		//randomize
		Misc.randomize(byGlobal, 0l);

		//chunk
		ArrayList<GeneExonSample>[][] chunks = chunk(byGlobal, minNumInEachChunk);

		//make objects to run on independent threads
		chunkThreads = new PoReDataChunk[chunks.length];
		for (int i=0; i< chunks.length; i++){
			chunkThreads[i] = new PoReDataChunk("batch"+i, chunks[i], this);
		}
	}

	/**Splits an object[] into chunks containing the minNumEach. Any remainder is evenly distributed over the prior.
	 * Note this is by reference, the array is not copied. */
	public static ArrayList<GeneExonSample>[][] chunk (ArrayList<GeneExonSample>[] s, int minNumEach){
		//watch out for cases where the min can't be met
		int numChunks = s.length/minNumEach;
		if (numChunks == 0) return new ArrayList[][]{s};

		double numLeftOver = (double)s.length % (double)minNumEach;

		int[] numInEach = new int[numChunks];
		for (int i=0; i< numChunks; i++) numInEach[i] = minNumEach;

		while (numLeftOver > 0){
			for (int i=0; i< numChunks; i++) {
				numInEach[i]++;
				numLeftOver--;
				if (numLeftOver == 0) break;
			}
		}
		//build chunk array
		ArrayList<GeneExonSample>[][] chunks = new ArrayList[numChunks][];
		int index = 0;
		//for each chunk
		for (int i=0; i< numChunks; i++){
			//create container and fill it
			ArrayList<GeneExonSample>[] sub = new ArrayList[numInEach[i]];
			for (int j=0; j< sub.length; j++) sub[j] = s[index++];
			chunks[i] = sub;
		}
		return chunks;
	}

	/**Just for cases*/
	public void makeGeneExonSamples(){
		//start counter for all exons to match Alun's code, starts with 1
		int globalExonIndex = 1;
		ArrayList<GeneExonSample> gesAL = new ArrayList<GeneExonSample>();
		
		//for each gene, must walk through so that these are sorted by chrom
		for (UCSCGeneLine gene: genes){
			//one to analyze?
			String geneName = gene.getDisplayNameThenName();
			if (geneNamesWithMinimumCounts.contains(geneName) == false) continue;
			
			//merging exon counts?
			int numExons;
			int[][] caseExonCounts;
			if (mergeExonCounts){
				numExons = 1;
				caseExonCounts = fetchMergedExonCounts(geneName, cases);
			}
			else {
				numExons = gene.getExons().length;
				caseExonCounts = fetchExonCounts(geneName, numExons, cases);
			}

			//for each exon
			for (int i=0; i< numExons; i++){
				//count total for this exon across all the cases, skip if too few
				int total = 0;
				for (int j=0; j<caseExonCounts.length; j++) total += caseExonCounts[j][i];
				if (total < minimumCounts) continue;
				numExonsProcessed++;
				//for each sample
				for (int j=0; j<caseExonCounts.length; j++){
					GeneExonSample g = new GeneExonSample(gene, globalExonIndex-1, (short)i, (short)j, caseExonCounts[j][i]);
					gesAL.add(g);
				}
				globalExonIndex++;
			}
		}
		//save GES
		ges = new GeneExonSample[gesAL.size()];
		gesAL.toArray(ges);
	}


	/**For case controls*/
	public void writeCountTableMatrix(){
		countTableFile = new File(saveDirectory, "countTable.txt");

		try {
			//write matrix of name, t,t,t...c,c,c,c...d,d,d, to file for genes with observations
			PrintWriter out = new PrintWriter( new FileWriter(countTableFile));

			//print header
			out.print("#GeneName\tExonIndex\tCoordinates");
			//for each sample
			for (int i=0; i<caseSampleFileNames.length; i++){
				out.print("\t");
				out.print(Misc.removeExtension(caseSampleFileNames[i]));
				if (controls != null){
					out.print("\t");
					out.print(Misc.removeExtension(controlSampleFileNames[i]));
				}
			}
			out.println();
			
			//for each gene, must walk through so that these are sorted by chrom
			for (UCSCGeneLine gene: genes){
				//one to analyze?
				String geneName = gene.getDisplayNameThenName();
				if (geneNamesWithMinimumCounts.contains(geneName) == false) continue;

				ExonIntron[] exons = gene.getExons();
				int numExons = gene.getExons().length;
				String chr= gene.getChrom();

				//collect exon counts int[sample][exon
				int[][] caseExonCounts = fetchExonCounts(geneName, numExons, cases);
				int[][] controlExonCounts = null;
				if (controls != null) controlExonCounts = fetchExonCounts(geneName, numExons, controls);

				//for each exon
				for (int i=0; i< numExons; i++){
					//print geneName and exon index
					StringBuilder sb = new StringBuilder(geneName+"\t"+i+"\t"+chr+":"+exons[i].getStartStopString());
					int total = 0;
					//for each sample
					for (int j=0; j<caseExonCounts.length; j++){
						sb.append("\t");
						sb.append(caseExonCounts[j][i]);
						total += caseExonCounts[j][i];
						if (controls != null){
							sb.append("\t");
							sb.append(controlExonCounts[j][i]);
							total += controlExonCounts[j][i];
						}
					}
					if (total >= minimumCounts) out.println(sb);
				}
			}
			out.close();
		}
		catch (Exception e){
			System.err.println("\nProblem writing out count table.");
			e.printStackTrace();
		}
	}

	private int[][] fetchExonCounts(String geneName, int numExons, Replica[] samples) {
		int[][] exonCounts = new int[samples.length][];
		int[] zero = new int[numExons];
		for (int i=0; i< samples.length; i++){
			GeneCount gc = samples[i].getGeneCounts().get(geneName);
			if (gc == null) exonCounts[i] = zero;
			else exonCounts[i] = gc.getExonCounts();
		}
		return exonCounts;
	}

	/**Returns total gene count.*/
	private int[][] fetchMergedExonCounts(String geneName, Replica[] samples) {
		int[][] exonCounts = new int[samples.length][];
		int[] zero = new int[1];
		for (int i=0; i< samples.length; i++){
			GeneCount gc = samples[i].getGeneCounts().get(geneName);
			if (gc == null) exonCounts[i] = zero;
			else exonCounts[i] = new int[]{gc.getCount()};
		}
		return exonCounts;
	}


	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table, sort by chromosome and position
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(refSeqFile, 0);
		genes = reader.getGeneLines();
		Arrays.sort(genes, new UCSCGeneLineChromComparator());
		
		if (genes == null || genes.length == 0) Misc.printExit("\nProblem loading your USCS gene model table or bed file? No genes/ regions?\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check gene name is unique
		if (reader.uniqueGeneNames() == false) Misc.printExit("\nDuplicate gene names were found in your gene / bed file, these must be unique.\n");
		chromGenes = reader.getChromSpecificGeneLines();
		//exclude sex chromosomes?
		if (excludeSexChromosomes){
			System.out.println("\tExcluding sex chromosomes");
			chromGenes.remove("chrY");
			chromGenes.remove("Y");
			chromGenes.remove("chrX");
			chromGenes.remove("X");
		}
		//load name2Gene
		name2Gene = new HashMap<String, UCSCGeneLine>();
		for (String chr: chromGenes.keySet()){
			for (UCSCGeneLine gl: chromGenes.get(chr)) {
				String name = gl.getDisplayNameThenName();
				name2Gene.put(name, gl);
			}
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PoReCNV(args);
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
					case 's': saveDirectory = new File(args[++i]); break;
					case 'b': bamDirectory = new File(args[++i]); break;
					case 'v': vcfDirectory = new File(args[++i]); break;
					case 'c': caseSampleFileNames = Misc.COMMA.split(args[++i]); break;
					case 'n': controlSampleFileNames = Misc.COMMA.split(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'a': alunRScript = new File(args[++i]); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'g': genomeVersion = args[++i]; break;
					case 'w': mergeExonCounts = true; annoType = "Gene"; break;
					case 'x': excludeSexChromosomes = false; break;
					case 'm': minNumInEachChunk = Integer.parseInt(args[++i]); break;
					case 'e': minimumCounts = Integer.parseInt(args[++i]); break;
					case 'd': maxAlignmentsDepth = Integer.parseInt(args[++i]); break;
					case 'l': minimumLog2Ratio = Float.parseFloat(args[++i]); break;
					case 'p': minimumAdjPVal = Float.parseFloat(args[++i]); break;
					case 't': numberConcurrentThreads = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}

		//check genomeVersion
		if (genomeVersion == null) Misc.printErrAndExit("\nError: please provide a genome version, (e.g. H_sapiens_Mar_2006)\n");

		//look for bam files
		if (bamDirectory == null) Misc.printErrAndExit("\nError: cannot find your directory of bam files?\n");
		File[] bamFiles = IO.extractFiles(bamDirectory, ".bam");
		if (bamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any bam files?\n");
		OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, false);

		//check for case names
		if (caseSampleFileNames == null ) {
			File[] bams = IO.fetchFilesRecursively(bamDirectory, ".bam");
			if (bams == null || bams.length == 0) Misc.printErrAndExit("\nError: cannot find any xxx.bam files to process?!\n");
			caseSampleFileNames = new String[bams.length];
			for (int i=0; i< bams.length; i++) caseSampleFileNames[i] = bams[i].getName();
		}
		else {
			for (int i=0; i< caseSampleFileNames.length; i++){
				File bam = new File(bamDirectory, caseSampleFileNames[i]);
				if (bam.exists() == false) Misc.printErrAndExit("Error: cannot find the "+caseSampleFileNames[i]+" in the bam directory?!\n");
			}
		}
		passingExonsByCase = new int[caseSampleFileNames.length];
		
		//parse vcf file names?
		if (vcfDirectory != null){
			//parse corresponding vcf files
			for (int i=0; i< caseSampleFileNames.length; i++){
				String name = Misc.removeExtension(caseSampleFileNames[i]) + ".vcf.gz";
				File vcf = new File(vcfDirectory, name);
				if (vcf.exists() == false) Misc.printErrAndExit("Error: cannot find "+name+" in the vcf directory?!\n");
				//look for index
				File index = new File(vcfDirectory, name+".tbi");
				if (index.exists() == false) Misc.printErrAndExit("Error: cannot find "+name+".tbi in the vcf directory?!\n");
			}
		}

		//check controls
		/*if (controlSampleFileNames != null){
			//same number?
			if (controlSampleFileNames.length != caseSampleFileNames.length) Misc.printErrAndExit("\nError: the number of case and control sample file names differ?\n");
			for (int i=0; i< controlSampleFileNames.length; i++){
				File bam = new File(bamDirectory, controlSampleFileNames[i]);
				if (bam.exists() == false) Misc.printErrAndExit("Error: cannot find the "+controlSampleFileNames[i]+" in the bam directory?!\n");
			}
		}*/

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdirs();

		//check for R and required libraries, don't need it if they just want the first and last 1/3 count table
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		//find Alun's script
		if (alunRScript == null || alunRScript.canRead()== false) {
			Misc.printErrAndExit("\nError: Cannot find or read Alun's RScript? -> "+alunRScript+"\n");
		}

		//number of threads to use?
		if (numberConcurrentThreads == 0) {
			numberConcurrentThreads = Runtime.getRuntime().availableProcessors();
		}
		
		printParams();

	}	
	
	public void printParams(){
		StringBuilder sb = new StringBuilder("Run parameters:\n");
		sb.append("Results dir\t"); sb.append(saveDirectory); sb.append("\n");
		sb.append("Bam dir    \t"); sb.append(bamDirectory); sb.append("\n");
		sb.append("Annotation \t"); sb.append(refSeqFile); sb.append("\n");
		sb.append("RScript    \t"); sb.append(alunRScript); sb.append("\n");
		sb.append("Genome build\t"); sb.append(genomeVersion); sb.append("\n");
		sb.append("Anno type\t"); sb.append(annoType); sb.append("\n");
		sb.append("Excld sex chrs\t"); sb.append(excludeSexChromosomes); sb.append("\n");
		sb.append("Min counts\t"); sb.append(minimumCounts); sb.append("\n");
		sb.append("Max counts\t"); sb.append(maxAlignmentsDepth); sb.append("\n");
		sb.append("Min log2Rto\t"); sb.append(minimumLog2Ratio); sb.append("\n");
		sb.append("Max adjPVal\t"); sb.append(minimumAdjPVal); sb.append("\n");
		sb.append("Num per batch\t"); sb.append(minNumInEachChunk); sb.append("\n");
		sb.append("Num threads\t"); sb.append(numberConcurrentThreads); sb.append("\n");
		System.out.println(sb);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Po Re CNV: April 2015                             **\n" +
				"**************************************************************************************\n" +
				"Uses poisson regression and Pearson residuals to identify exons or genes whose counts\n"+
				"differ significantly from the fitted value base on all the exon sample counts. This\n"+
				"app wraps an algorithm developed by Alun Thomas.  Data tracks are generated for the\n"+
				"residuals and log2(observed/ fitted counts) as well as detailed spreadsheets. Use for\n"+
				"identifying CNVs in next gen seq datasets with > 10 normal samples.\n"+

				"\nRequired Options:\n"+
				"-s Save directory.\n"+
				"-b BAM file directory, sorted and indexed by coordinate.\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. Tab delimited, see RefSeq Genes\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889\n"+
				"-a Alun Thomas R script file.\n"+
				"-g Genome Version  (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +

				"\nDefault Options:\n"+
				"-c BAM file names for samples to process, comma delimited, no spaces, defaults to all.\n"+
				"-r Full path to R, defaults to /usr/bin/R\n"+
				"-l Minimum log2(obs/exp) for inclusion in the pass spreadsheet, defaults to 0.585\n"+
				"-p Maximum adjusted p-value for inclusion in the pass spreadsheet, defaults to 0.01\n"+
				"-e Minimum all sample exon count for inclusion in analysis, defaults to 20.\n"+
				"-d Max per sample exon alignment depth, defaults to 50000. Exons containing higher\n"+
				"       counts are ignored.\n"+
				"-t Number concurrent threads to run, defaults to the max available to the jvm.\n"+
				"-m Minimum number exons per data chunk, defaults to 1500.\n"+
				"-w Examine whole gene counts for CNVs, defaults to exons.\n"+
				"-x Keep sex chromosomes (X,Y), defaults to removing.\n"+
				//"-v VCF file directory, sorted and indexed. Name.vcf.gz must correspond with Name.bam\n"+
				//"     for matching samples. Use to annotate deletions for possible heterozygous snps. \n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/PoReCNV -s PRCnvResults/  -b BamFiles\n"+
				"       -u hg19EnsGenes.ucsc.gz -a RBambedSource.R  -g H_sapiens_Feb_2009 -d\n\n" +

				"**************************************************************************************\n");

	}


	public File getSaveDirectory() {
		return saveDirectory;
	}


	public boolean isDeleteTempFiles() {
		return deleteTempFiles;
	}


	public File getFullPathToR() {
		return fullPathToR;
	}


	public File getAlunRScript() {
		return alunRScript;
	}


	public float getMinimumLog2Ratio() {
		return minimumLog2Ratio;
	}


	public float getMinimumAdjPVal() {
		return minimumAdjPVal;
	}


	public int getNumExonsProcessed() {
		return numExonsProcessed;
	}
}
