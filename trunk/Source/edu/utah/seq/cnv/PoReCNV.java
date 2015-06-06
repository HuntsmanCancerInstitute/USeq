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
import htsjdk.variant.vcf.VCFFileReader;

/** App to detect CNV's in datasets containing mostly normal genomes.  Wraps Alun Thomas's algorithm.
 * @author Nix
 * */
public class PoReCNV {

	//user defined fields
	private File saveDirectory;
	private File bamDirectory;
	private File vcfDirectory;
	private String[] sampleFileNames;
	private File fullPathToR = new File ("/usr/bin/R");
	private File refSeqFile;
	private int minimumCounts = 20;
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
	private int maxGap = 25000;

	//internal fields
	private int numberAdjacentExonsToScan = 15;
	private Gzipper outPassBed = null;
	private File graphDirectory;
	private UCSCGeneLine[] genes;
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private HashMap<String, UCSCGeneLine> name2Gene;
	public static final Pattern BAD_NAME = Pattern.compile("(.+)/[12]$");
	private LinkedHashSet<String> geneNamesWithMinimumCounts = new LinkedHashSet<String>();
	private Replica[] sampleReplicas = null;
	private PoReDataChunk[] chunkThreads;
	private GeneExonSample[][] sampleGES = null;
	private double[] totalCounts = null;
	private SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

	//container for samples and their exons that pass the minimumCounts threshold
	private GeneExonSample[] ges = null;
	private int numExonsProcessed = 0;
	private int numberPassingExons = 0;
	private int[] passingExonsBySample = null;

	//from R
	private float[] significanceThresholds = null;
	private float minimumSignificanceThreshold;

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

		//find overlapping variants in exons with a deletion?
		//if (vcfDirectory != null) searchForHetSNPs();

		//calculate count zscores
		System.out.println("On target fragment counts:");
		calculateZScores();

		//search for adjacent high scoring exons/ genes
		searchForAdjacents();

		//write out graph files for each sample
		System.out.println("\nSaving graphs and spreadsheets...");
		exportGeneExonSpreadSheets();
		exportSampleGraphs();

		System.out.println("\n"+numExonsProcessed+" "+annoType+"s processed, "+numberPassingExons+" passed thresholds ("+minimumLog2Ratio+" Lg2Rto, "+minimumAdjPVal+" AdjPVal)");

		//convert graphs to useq format
		new Bar2USeq(graphDirectory, true);
	}

	private void setSignificanceThreshold() {
		significanceThresholds = chunkThreads[0].getSignificanceThresholds();
		System.out.println("\nSignificance thresholds for adjacent "+annoType+"s:\n"+ Misc.floatArrayToString(significanceThresholds, ", ")+"\n");
		if (significanceThresholds == null) Misc.printErrAndExit("\nError: did not find significance thresholds for the data trunks?\n");
		minimumSignificanceThreshold = significanceThresholds[significanceThresholds.length-1];
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
			for (int i=0; i< sampleReplicas.length; i++){
				String name = sampleReplicas[i].getNameNumber();
				sb.append("\t"+name+" Lg2(Ob/Ex)\t"+name+" Res\t"+name+" Z\t"+name+" Counts");
			}
			sb.append("\tLog2Rto "+minimumLog2Ratio);
			sb.append("\tResidual "+minimumAdjPVal+" sig threshold +/- for sequential occurances "+Misc.floatArrayToString(significanceThresholds, ","));

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
				//if (passesThresholds(g)) {
				if (g.isPassedThresholds()){
					sampPass.add(sampleReplicas[g.getSampleIndex()].getNameNumber());
					passingExonsBySample[g.getSampleIndex()]++;
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
		if (Math.abs(g.getObsExpLgRto()) < minimumLog2Ratio || Math.abs(g.getResidual()) < significanceThresholds[0] ) return false;
		return true;
	}
	
	private boolean passesThresholds(GeneExonSample first, GeneExonSample second){
		//check obs/exp 
		if ( Math.abs(first.getObsExpLgRto()) < minimumLog2Ratio || Math.abs(second.getObsExpLgRto()) < minimumLog2Ratio) return false;
		
		//check minimum residual
		if ( Math.abs(first.getResidual()) < minimumSignificanceThreshold || Math.abs(second.getResidual()) < minimumSignificanceThreshold) return false;
		
		//check sign
		boolean firstIsPositive = first.getResidual() > 0;
		boolean secondIsPositive = second.getResidual() > 0;
		if (firstIsPositive != secondIsPositive) return false;
		
		//check distance
		int[] firstStartStop = getStartStop(first);
		int[] secondStartStop = getStartStop(second);
		int da = distanceApart(firstStartStop, secondStartStop);
		if (da > maxGap) return false; 
		
		return true;
	}
	
	private boolean passesMinimumThresholds(GeneExonSample first){
		//check obs/exp 
		if ( Math.abs(first.getObsExpLgRto()) < minimumLog2Ratio) return false;
		//check minimum residual
		if ( Math.abs(first.getResidual()) < minimumSignificanceThreshold) return false;
		return true;
	}
	
	private int distanceApart(int[] f, int[] s) {
			//s is left of f
			if (s[1] < f[0]) return f[0] - s[1];
			else if (s[0] > f[1]) return s[0] - f[1];
			//overlap
			return -1;
	}

	private int[] getStartStop(GeneExonSample g){
		UCSCGeneLine gene = g.getGene();
		ExonIntron[] exons = gene.getExons();
		int exonIndex = g.getExonIndex();
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
		return new int[]{start, stop};
	}
	
	private void searchForAdjacents() {

		try{
			outPassBed = new Gzipper(new File(saveDirectory, "results"+annoType+"Pass.bed.gz"));
			outPassBed.println("#Chr\tStart\tStop\tSampleName:lg2Rto(Obs/Exp)_res_z_count;\tNumInGroup\tStrand");

			//for each sample
			for (GeneExonSample[] ges : sampleGES){

				//walk ges and process for each chromosome 
				String currChr = ges[0].getGene().getChrom();
				String sampleName = sampleReplicas[ges[0].getSampleIndex()].getNameNumber();
				ArrayList<GeneExonSample> chrGes = new ArrayList<GeneExonSample>();

				for (int i=0; i< ges.length; i++){
					//diff chrom?
					if (ges[i].getGene().getChrom().equals(currChr) == false) {
						//diff chrom, write out old
						searchChromAdjacents(sampleName, currChr, chrGes);
						//reset and set new
						currChr = ges[i].getGene().getChrom();
						chrGes.clear();
					}
					chrGes.add(ges[i]);
				}
				//process last!
				searchChromAdjacents(sampleName, currChr, chrGes);
			}
			outPassBed.close();

		} catch (Exception e){
			System.err.println("\nERROR: problem finding and writing bed file of adjacent passing CNVs.\n");
			e.printStackTrace();
		}
	}
	
	private void searchChromAdjacents(String sampleName, String chromosome, ArrayList<GeneExonSample> chrGes) throws IOException {
		//walk through ges looking for adjacent exons/genes that exceed threholds
		int size = chrGes.size();
		//build block of passing obs/exp, dist, and minimum sig threshold
		ArrayList<GeneExonSample> block = new ArrayList<GeneExonSample>();
		
		//find first with passing obs/exp,
		GeneExonSample prior = null;
		int i=0;
		for (; i<size; i++){
			prior = chrGes.get(i);
			//pass?
			if (passesMinimumThresholds(prior)){
				block.add(prior);
				break;
			}
		}
		//advance
		i++;
		 
		//scan remainder
		for (; i<size; i++){
			GeneExonSample next = chrGes.get(i);
			//fails thresholds?
			if (passesThresholds(prior, next) == false) {
				//process and clear
				ArrayList<int[]> adjClustIndex = processAdjacentBlock(block);
				if (adjClustIndex != null) printClusters(sampleName, chromosome, block, adjClustIndex);
				//make method
				block.clear();
				//look for next that passes
				//does next pass single thresholds?
				if (passesMinimumThresholds(next)){
					block.add(next);
					prior = next;
				}
				//nope look for it
				else{
					i++;
					for (; i<size; i++){
						prior = chrGes.get(i);
						//pass?
						if (passesMinimumThresholds(prior)){
							block.add(prior);
							break;
						}
					}
				}
			}
			else {
				block.add(next);
				prior = next;
			}
		}
		//process last!
		ArrayList<int[]> adjClustIndex = processAdjacentBlock(block);
		if (adjClustIndex != null) printClusters(sampleName, chromosome, block, adjClustIndex);
	}

	/*This is specific to a particular sample and chromosome*/
	private void printClusters(String sampleName, String chromosome, ArrayList<GeneExonSample> block, ArrayList<int[]> adjClustIndex) throws IOException {
		
		//for each cluster in a sample
		for (int[] indexes : adjClustIndex){
			
			//get first and last for positions
			GeneExonSample first = block.get(indexes[0]);
			GeneExonSample last = block.get(indexes[indexes.length-1]);
			int start = first.getGene().getExons()[first.getExonIndex()].getStart();
			int stop = last.getGene().getExons()[last.getExonIndex()].getEnd();

			//set first
			StringBuilder rtoResCnts = new StringBuilder();
			rtoResCnts.append(first.getDataStringUnderscore());
			first.setPassedThresholds(true);
			//System.out.println("NewCluster\n\t"+fetchGeneInfo(first) + first.getDataString());
			
			//for subsequent
			for (int i=1; i< indexes.length; i++){
				GeneExonSample g = block.get(indexes[i]);
				rtoResCnts.append(";");
				rtoResCnts.append(g.getDataStringUnderscore());
				g.setPassedThresholds(true);
				//System.out.println("\t"+fetchGeneInfo(g) + g.getDataString());
			}
			
			//print out bed line of chrom start stop sampleName;rto_res_cnts;rto_res_cnts..., # adjacent, .
			//System.out.println(chromosome+"\t"+start+"\t"+stop+"\t"+sampleName+ "-"+rtoResCnts+"\t"+indexes.length+"\t.");
			outPassBed.println(chromosome+"\t"+start+"\t"+stop+"\t"+sampleName+ ":"+rtoResCnts+"\t"+indexes.length+"\t.");
			
		}
		
	}

	/**Returns null if none found passing thresholds. Otherwise int[]s of adjacent indexes where all pass thresholds.*/
	private ArrayList<int[]> processAdjacentBlock(ArrayList<GeneExonSample> block) {
		int size = block.size();
		if (size == 0) return null;
		
		ArrayList<int[]> clusteredIndexes = null;
		//OK all of these have passed the distance, obs/exp ratio, and minimum residual
		//just one?
		if (size == 1){
			if (Math.abs(block.get(0).getResidual()) < significanceThresholds[0]) return null;
			clusteredIndexes = new ArrayList<int[]>();
			clusteredIndexes.add(new int[]{0});
		}
		//nope more than one so have to do heavy weight scan
		else {
			//get residual scores
			float[] residuals = new float[block.size()];
			for (int i=0; i< residuals.length; i++) residuals[i] = Math.abs(block.get(i).getResidual());
			//make clusters, returns clusters of one
			clusteredIndexes = findLargestPassingAdjacentBlocks(residuals);
			if (clusteredIndexes.size() == 0) return null;
		}
		return clusteredIndexes;
	}
	
	/**Finds the largest cluster of adjacent residuals.*/
	public ArrayList<int[]> findLargestPassingAdjacentBlocks(float[] residuals){
		int index = significanceThresholds.length - 1;
		if ((residuals.length -1) < index) index = residuals.length-1;
		ArrayList<LinkedHashSet<Integer>> passingHashes = new ArrayList<LinkedHashSet<Integer>>();
		while (true){
			//scan and make runs
			float threshold = significanceThresholds[index];
			int numInARow = index+1;
			//System.out.println("Thres "+threshold+", numInRow "+numInARow);
			//for each position, try to make numInARow
			for (int i=0; i< residuals.length; i++){
				boolean pass = true;
				LinkedHashSet<Integer> passingIndexes = new LinkedHashSet<Integer>();
				for (int j=0; j<numInARow; j++){
					int testingIndex = i+j;
					//too far?
					if (testingIndex>= residuals.length) {
						pass = false;
						break;
					}
					//pass threshold?
					if (residuals[testingIndex] < threshold){
						pass = false;
						break;
					}
					passingIndexes.add(new Integer(testingIndex));
				}
				if (pass){
					//String x = "nada";
					//if (passingIndexes.size() !=0) x= passingIndexes.toString();
					//System.out.println("\t"+i+" "+pass+" "+ x);
					if (passingIndexes.size() > 0) passingHashes.add(passingIndexes);
				}
			}
			index--;
			if (index < 0) break;
		}
		//System.out.println("\nBiggest clusters");

		//need to find the biggest cluster for each index
		HashMap<Integer, HashSet<Integer>> biggestClusters = new HashMap<Integer, HashSet<Integer>>();
		
		//for each cluster
		for (int i=0; i< passingHashes.size(); i++){
			//get the test cluster
			HashSet<Integer> cluster = passingHashes.get(i);
			//for each index
			for (int j=0; j< residuals.length; j++){
				Integer test = new Integer(j);
				//index already found? 
				if (biggestClusters.containsKey(test)) continue;
				//contained in the test cluster?
				if (cluster.contains(test)) biggestClusters.put(test, cluster);
			}
		}
		
		//collapse em
		HashSet<String> collapsedClusterNames = new HashSet<String>();
		ArrayList<int[]> collapsedClusters = new ArrayList<int[]>();
		for (Integer i : biggestClusters.keySet()){
			HashSet<Integer> cluster = biggestClusters.get(i);
			String stringRep = cluster.toString();
			//System.out.println(i+" -> "+biggestClusters.get(i));
			if (collapsedClusterNames.contains(stringRep) == false){
				collapsedClusterNames.add(stringRep);
				collapsedClusters.add(Num.hashSetToInt(cluster));
			}
		}
		return collapsedClusters;
	}
	
	/*private void searchForHetSNPs(){
		
		//for each sample
		for (int i=0; i< sampleGES.length; i++){
			
			//make a VCF reader
			String name = Misc.removeExtension(sampleFileNames[i]) + ".vcf.gz";
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
					//should be just one for most samples
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
					//System.out.println(sGes[j].toString()+"\t"+sGes[j].getGenotypePosition());
				}
			}
			vr.close();
		}
	}*/

	private void exportSampleGraphs() {
		//make dirs
		graphDirectory = new File (saveDirectory, "DataTracks");
		graphDirectory.mkdir();
		File resDir = new File (graphDirectory, annoType+"Residual");
		resDir.mkdir();
		File rtoDir = new File (graphDirectory, annoType+"ObsExpLog2Rto");
		rtoDir.mkdir();
		File zDir = new File (graphDirectory, annoType+"ZScores");
		zDir.mkdir();
		
		//for each write files
		System.out.println("\nResidual stats:\nDataset\t"+annoType+"sPass\tMean\tMedian\tStdDev\tMin\tMax\t10th\t90th");
		for (int i=0; i< sampleGES.length; i++){
			//make sample specific dirs
			File resSampDir = new File (resDir, "res_"+sampleReplicas[i].getNameNumber());
			resSampDir.mkdir();
			resSampDir.deleteOnExit();
			File rtoSampDir = new File (rtoDir, "oe_"+sampleReplicas[i].getNameNumber());
			rtoSampDir.mkdir();
			rtoSampDir.deleteOnExit();
			File zSampDir = new File (zDir, "z_"+sampleReplicas[i].getNameNumber());
			zSampDir.mkdir();
			zSampDir.deleteOnExit();
			saveGraphs(sampleReplicas[i].getNameNumber(), sampleGES[i], resSampDir, rtoSampDir, zSampDir);
			//calc residual stats
			float[] res = getResiduals(sampleGES[i]);
			Arrays.sort(res);
			String stats = Num.statFloatArray(res);
			System.out.println(sampleFileNames[i]+ "\t"+passingExonsBySample[i]+ "\t"+ stats);
		}
	}
	
	private void calculateZScores(){
		//container for every exon
		StandardDeviation[] sd = new StandardDeviation[sampleGES[0].length];
		for (int i=0; i< sd.length; i++) sd[i] = new StandardDeviation();
		
		//for each sample
		for (int i=0; i< sampleGES.length; i++){
			//calc residual stats
			double[] rawCounts = getCounts(sampleGES[i]);
			//normalize
			double total = totalCounts[i];
			for (int j=0; j< rawCounts.length; j++) {
				double normCount = rawCounts[j]/total;
				sd[j].count(normCount);
			}
			System.out.println(Misc.removeExtension("\t"+sampleFileNames[i])+"\t"+(int)totalCounts[i]);
		}
		
		//now calc and set z-score for each exon
		//for each sample
		for (int i=0; i< sampleGES.length; i++){
			double total = totalCounts[i];
			//for each exon
			for (int j=0; j< sd.length; j++){
				double count = sampleGES[i][j].getCount()/total;
				double z = sd[j].getZScore(count);
				sampleGES[i][j].setZscore((float)z);
			}
		}
		
		
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
	
	private double[] getCounts(GeneExonSample[] ges) {
		double[] c = new double[ges.length];
		for (int i=0; i< ges.length; i++) c[i] = ges[i].getCount();
		return c;
	}

	private void saveGraphs(String sampleName, GeneExonSample[] ges, File resDir, File rtoDir, File zDir) {
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

			SmoothingWindow sm = new SmoothingWindow(start, stop, new float[]{ges[i].getResidual(), ges[i].getObsExpLgRto(), ges[i].getZscore()});

			//diff chrom?
			if (ges[i].getGene().getChrom().equals(currChr) == false) {
				//diff chrom, write out old
				Info info = new Info(sampleName, genomeVersion, currChr, ".", 0, null);
				SmoothingWindow[] sms = new SmoothingWindow[smAL.size()];
				smAL.toArray(sms);
				saveStairStepGraph(0, sms, info, resDir, "#FF0000", true); //red
				saveStairStepGraph(1, sms, info, rtoDir, "#0000FF", true); //blue
				saveStairStepGraph(2, sms, info, zDir, "#00FF00", true); //green

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
		saveStairStepGraph(1, sms, info, rtoDir, "#0000FF", true); //blue
		saveStairStepGraph(2, sms, info, zDir, "#00FF00", true); //green
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

		//load samples
		System.out.println("\nLoading samples...");
		//serialized version present?
		File serDir = new File (saveDirectory, "Ser");
		File geneTableInfo = new File (serDir, "geneInfo");

		if (serDir.exists()) {
			System.out.println("\tWARNING: possibly loading sample count data from prior run. This is likely a good thing. If in doubt, delete "+serDir+" and restart to reload from bam files.");
			//check gene table name and number
			if (geneTableInfo.exists() == false) Misc.printErrAndExit("\nError: the geneInfo file does not exist in the Ser directory.  Delete "+serDir+" and restart.\n");
			String[] info = IO.loadFile(geneTableInfo);
			if (info[0].equals(refSeqFile.getName()) == false || Integer.parseInt(info[1])!= genes.length) Misc.printErrAndExit("\nError: the saved gene table name ("+
					info[0]+" vs "+refSeqFile.getName()+") or gene number ("+info[1]+" vs "+genes.length+") don't match.  Delete "+serDir+" and restart.\n");
		}
		else {
			serDir.mkdirs();
			String gi = refSeqFile.getName() +"\n"+genes.length;
			IO.writeString(gi, geneTableInfo);
		}

		sampleReplicas = new Replica[sampleFileNames.length];
		for (int i=0; i< sampleFileNames.length; i++){
			System.out.print("\t"+sampleFileNames[i]);
			//look for saved replica
			File serRep = new File (serDir, sampleFileNames[i]+".ser");
			if (serRep.exists()) {
				sampleReplicas[i] = (Replica)IO.fetchObject(serRep);
				System.out.print("\tser\t");
			}
			else {
				//load from bam
				File bam = new File(bamDirectory, sampleFileNames[i]);
				if (bam.exists() == false) Misc.printErrAndExit("\nCannot find the associated bam file for "+sampleFileNames[i]);
				sampleReplicas[i] = new Replica(Misc.removeExtension(sampleFileNames[i]), bam);
				loadReplica(sampleReplicas[i]);
				//save it
				IO.saveObject(serRep, sampleReplicas[i]);
				System.out.print("\tbam\t");
			}
			System.out.println(sampleReplicas[i].getTotalCounts());
		}

		findGeneNamesWithMinimumCounts();
		
		System.out.println();
		System.out.println(geneNamesWithMinimumCounts.size()+" genes, with >= "+minimumCounts+" counts, will be examined for copy number alterations.");
		
		//container for total counts that hit genes of interest
		totalCounts = new double[sampleReplicas.length];
	}

	private void findGeneNamesWithMinimumCounts() {
		//for each gene, does any sample have enough counts?
		for (UCSCGeneLine gene: genes){
			String geneName = gene.getDisplayNameThenName();
			for (Replica r: sampleReplicas){
				GeneCount sc = r.getGeneCounts().get(geneName);
				if (sc!= null && sc.getCount() >= minimumCounts) {
					geneNamesWithMinimumCounts.add(geneName);
					break;
				}
			}
		}
	}

	public void loadBlocks(SAMRecord sam, int chrStartBp, ArrayList<Integer>[] bpNames){
		NameInteger nameIndex = null;
		ArrayList<int[]> blocks = DefinedRegionDifferentialSeq.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
		//add name to each bp
		for (int[] b : blocks){
			int start = b[0] - chrStartBp;
			int stop = b[1] - chrStartBp;
			//need to watch for out of bounds issues, sometimes the length of the chromosome is incorrect in the bam header.
			if (stop > bpNames.length) stop = bpNames.length;
			for (int i=start; i < stop; i++){
				//never seen before?
				if (bpNames[i] == null) bpNames[i] = new ArrayList<Integer>();
				if (nameIndex == null) nameIndex = fetchFragmentNameIndex(sam);
				bpNames[i].add(nameIndex.index);
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
		try {
			//the idea here is to take a sorted bam file and add the reads to an ArrayList, basically every base should have a ArrayList of reads that overlap that base
			//one can then count the number of overlapping fragments for an exon or a collection of exons by hashing the ArrayList
			//assumes each read (first or second) has the same name

			//make reader
			SamReader reader = factory.open(replica.getBamFile());

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
			while (iterator.hasNext()){
				sam = iterator.next();
				//unaligned? 
				if (sam.getReadUnmappedFlag()) continue;
				chrom = sam.getReferenceName();

				if (chromGenes.containsKey(chrom) == false) continue;		
				priorChroms.add(chrom);
				chrStartBp = sam.getUnclippedStart()-1;
				bpNames = new ArrayList[chromLength.get(chrom) - chrStartBp];
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
			loadBlocks (sam, chrStartBp, bpNames);

			//for each record
			while (iterator.hasNext()){
				sam = iterator.next();

				//unaligned? 
				if (sam.getReadUnmappedFlag()) continue;

				//same chrom?
				if (sam.getReferenceName().equals(chrom)) loadBlocks (sam, chrStartBp, bpNames);
				else {
					//different chrom so time to scan
					loadGeneCounts(replica, bpNames, chrStartBp, chrom);

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
						bpNames = new ArrayList[chromLength.get(chrom) - chrStartBp];
						//load
						loadBlocks (sam, chrStartBp, bpNames);
					}

				}

			}

			if (chromGenes.containsKey(chrom)) loadGeneCounts(replica, bpNames, chrStartBp, chrom);

			reader.close();
			bpNames = null;
			iterator = null;
			seqs = null;
			chromLength = null;

		} catch (IOException e) {
			System.err.println("Problem loading replica "+replica.getBamFile());
			e.printStackTrace();
		}
	}

	private void loadGeneCounts(Replica replica, ArrayList<Integer>[] bpNames, int chrStartBp, String chromosome){
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
			for (int x=0; x< exons.length; x++){
				exonReads.clear();
				int start = exons[x].getStart() - chrStartBp;
				if (start < 0) start = 0;
				int end = exons[x].getEnd() - chrStartBp;
				if (end > lengthBpNames) end = lengthBpNames;
				//for each base in the exon, see if there is a read
				for (int y=start; y< end; y++){
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
			}
		}	
		//clean up
		allReads = null;
		exonReads = null;
	}

	public void buildGeneExonSamples(){
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
				int[][] exonCounts = fetchExonCounts(geneName, numExons, sampleReplicas);

				//for each exon
				for (int i=0; i< numExons; i++){
					//count total for this exon across all the samples, skip if too few
					int total = 0;
					for (int j=0; j<exonCounts.length; j++) total += exonCounts[j][i];
					if (total < minimumCounts) continue;
					numExonsProcessed++;

					//for each sample
					for (int j=0; j<exonCounts.length; j++){
						GeneExonSample g = new GeneExonSample(gene, globalExonIndex-1, (short)i, (short)j, exonCounts[j][i]);
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
			boolean outputSigLevels = (i == 0);
			chunkThreads[i] = new PoReDataChunk("batch"+i, chunks[i], this, outputSigLevels);
		}
	}

	/**Splits an object[] into chunks containing the minNumEach. Any remainder is evenly distributed over the prior.
	 * Note this is by reference, the array is not copied. */
	public static ArrayList<GeneExonSample>[][] chunk (ArrayList<GeneExonSample>[] s, int minNumEach){
		//watch out for samples where the min can't be met
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
			int[][] sampleExonCounts;
			int[][] sampleGeneCounts;
			if (mergeExonCounts){
				numExons = 1;
				sampleExonCounts = fetchMergedExonCounts(geneName, sampleReplicas);
				sampleGeneCounts = sampleExonCounts;
			}
			else {
				numExons = gene.getExons().length;
				sampleExonCounts = fetchExonCounts(geneName, numExons, sampleReplicas);
				//still need to get gene counts
				sampleGeneCounts = fetchMergedExonCounts(geneName, sampleReplicas);
			}
			
			//add gene counts to each sample
			for (int x = 0; x< totalCounts.length; x++) totalCounts[x] += sampleGeneCounts[x][0];

			//for each exon
			for (int i=0; i< numExons; i++){
				//count total for this exon across all the samples, skip if too few
				int total = 0;
				for (int j=0; j<sampleExonCounts.length; j++) total += sampleExonCounts[j][i];
				if (total < minimumCounts) continue;
				numExonsProcessed++;
				//for each sample
				for (int j=0; j<sampleExonCounts.length; j++){
					GeneExonSample g = new GeneExonSample(gene, globalExonIndex-1, (short)i, (short)j, sampleExonCounts[j][i]);
					gesAL.add(g);
				}
				globalExonIndex++;
			}
		}
		//save GES
		ges = new GeneExonSample[gesAL.size()];
		gesAL.toArray(ges);
		
		//split ges by sample for later use
		sampleGES = new GeneExonSample[sampleReplicas.length][numExonsProcessed];
		for (int i=0; i< ges.length; i++){
			sampleGES[ges[i].getSampleIndex()][ges[i].getGlobalExonIndex()] = ges[i];
		}
	}


	/**For case controls
	public void writeCountTableMatrix(){
		countTableFile = new File(saveDirectory, "countTable.txt");

		try {
			//write matrix of name, t,t,t...c,c,c,c...d,d,d, to file for genes with observations
			PrintWriter out = new PrintWriter( new FileWriter(countTableFile));

			//print header
			out.print("#GeneName\tExonIndex\tCoordinates");
			//for each sample
			for (int i=0; i<sampleFileNames.length; i++){
				out.print("\t");
				out.print(Misc.removeExtension(sampleFileNames[i]));
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
				int[][] caseExonCounts = fetchExonCounts(geneName, numExons, sampleReplicas);

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
	}*/

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

	/**Returns total gene count. int[sampleIndex][count but just one value]*/
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
					case 'c': sampleFileNames = Misc.COMMA.split(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'a': alunRScript = new File(args[++i]); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'g': genomeVersion = args[++i]; break;
					case 'w': mergeExonCounts = true; annoType = "Gene"; break;
					case 'x': excludeSexChromosomes = false; break;
					case 'm': minNumInEachChunk = Integer.parseInt(args[++i]); break;
					case 'e': minimumCounts = Integer.parseInt(args[++i]); break;
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
		if (bamFiles == null || bamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any bam files in "+bamDirectory);
		OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, false);
		
		//check for  names
		File serDir = new File (saveDirectory, "Ser");
		if (sampleFileNames == null ) {
			//look for bams
			File[] bams = IO.fetchFilesRecursively(bamDirectory, ".bam");
			HashSet<String> bamNames = new HashSet<String>();
			if (bams != null) for (File f: bams) bamNames.add(f.getName());
			//look for serialized
			if (serDir.exists()){
				File[] serBams = IO.extractFiles(serDir, ".bam.ser");
				if (serBams != null ) for (File f: serBams) bamNames.add(f.getName().replace(".bam.ser", ".bam"));
			}
			//any found?
			if (bamNames.size() == 0) Misc.printErrAndExit("\nError: cannot find any xxx.bam files to process?!\n");
			sampleFileNames = new String[bamNames.size()];
			bamNames.toArray(sampleFileNames);
		}
		else {
			for (int i=0; i< sampleFileNames.length; i++){
				File bam = new File(bamDirectory, sampleFileNames[i]);
				File serBam = new File(serDir, sampleFileNames[i]+".ser");
				if (bam.exists() == false && serBam.exists() == false) Misc.printErrAndExit("Error: cannot find the "+sampleFileNames[i]+" in the bam or ser directories?!\n");
			}
		}
		
		passingExonsBySample = new int[sampleFileNames.length];
		
		//check number
		if (sampleFileNames.length < 10) System.out.println("WARNING: less than 10 samples found! This application performs best with > 10 to fit the glm.\n");

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdirs();

		//check for R 
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
		sb.append("Min log2Rto\t"); sb.append(minimumLog2Ratio); sb.append("\n");
		sb.append("Max adjPVal\t"); sb.append(minimumAdjPVal); sb.append("\n");
		sb.append("Max gap adj\t"); sb.append(maxGap); sb.append("\n");
		sb.append("Num per batch\t"); sb.append(minNumInEachChunk); sb.append("\n");
		sb.append("Num threads\t"); sb.append(numberConcurrentThreads); sb.append("\n");
		System.out.println(sb);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Po Re CNV: June 2015                              **\n" +
				"**************************************************************************************\n" +
				"Uses Poisson regression and Pearson residuals to identify exons or genes whose counts\n"+
				"differ significantly from the fitted value base on all the exon sample counts. This\n"+
				"app wraps an algorithm developed by Alun Thomas.  Data tracks are generated for the\n"+
				"residuals and log2(observed/ fitted counts) as well as detailed spreadsheets. A bed\n"+
				"regions file of merged adjacent exons passing thresholds is also created. Use this\n"+
				"app for identifying CNVs in next gen seq datasets with > 10 normal samples.\n"+

				"\nRequired Options:\n"+
				"-s Save directory.\n"+
				"-b BAM file directory, sorted and indexed by coordinate. One bam per sample. \n"+
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
				"-l Minimum abs(log2(obs/exp)) for inclusion in the pass spreadsheet, defaults to 0.585\n"+
				"-p Maximum adjusted p-value for inclusion in the pass spreadsheet, defaults to 0.01\n"+
				"-e Minimum all sample exon count for inclusion in analysis, defaults to 20.\n"+
				"-d Max per sample exon alignment depth, defaults to 50000. Exons containing higher\n"+
				"       counts are ignored.\n"+
				"-t Number concurrent threads to run, defaults to the max available to the jvm.\n"+
				"-m Minimum number exons per data chunk, defaults to 1500.\n"+
				"-w Examine whole gene counts for CNVs, defaults to exons.\n"+
				"-x Keep sex chromosomes (X,Y), defaults to removing. Don't mix sexes!\n"+
				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/PoReCNV -s PRCnvResults/  -b BamFiles\n"+
				"       -u hg19EnsGenes.ucsc.gz -a RBambedSource.R  -g H_sapiens_Feb_2009 \n\n" +

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

	public int getNumberAdjacentExonsToScan() {
		return numberAdjacentExonsToScan;
	}
}
