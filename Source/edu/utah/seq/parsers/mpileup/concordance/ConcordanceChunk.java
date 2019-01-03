package edu.utah.seq.parsers.mpileup.concordance;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.parsers.mpileup.MpileupLine;
import edu.utah.seq.parsers.mpileup.MpileupSample;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;

public class ConcordanceChunk implements Runnable{
	
	//fields
	private boolean complete = false;
	private boolean failed = false;
	private BamConcordance bcc;
	private Bed[] regions;
	private String chunkName;
	private String cmd = null;
	private File tempBed;
	private File commonSnvBed;
	private Gzipper misMatchBed;
	private int minSnvDP;
	private double minAFForHom;
	private double minAFForMatch;
	private double minAFForHis;
	private int minBaseQuality;
	private int minMapQuality;
	private double maxIndel;
	private Histogram[] afHist;
	private Histogram[] chrXAfHist;
	private Similarity[] similarities;
	private long numMpileupLinesProc = 0;
	private IntervalST<String> commonSnvRegions = null;
	private String currCommonChr = "";
	private TabixReader tabixReader;
	
	public ConcordanceChunk(Bed[] regions, BamConcordance bcc, String chunkName) throws Exception{
		this.bcc = bcc;
		this.regions = regions;
		this.chunkName = chunkName;
		this.minSnvDP = bcc.getMinSnvDP();
		this.minAFForHom = bcc.getMinAFForHom();
		this.minAFForHis = bcc.getMinAFForHis();
		this.minAFForMatch = bcc.getMinAFForMatch();
		this.minBaseQuality = bcc.getMinBaseQuality();
		this.maxIndel = bcc.getMaxIndel();
		this.minMapQuality = bcc.getMinMappingQuality();
		this.commonSnvBed = bcc.getCommonSnvBed();
		if (commonSnvBed != null) tabixReader = new TabixReader(commonSnvBed.toString());
	}
	
	public void run(){
		try {	
			//write out bed
			tempBed = new File (bcc.getTempDirectory(), chunkName+"_temp.bed");
			if (Bed.writeToFile(regions, tempBed) == false) Misc.printErrAndExit("\nERROR: failed to write bed regions to "+tempBed);
			
			//make gzipper
			misMatchBed = new Gzipper(new File (bcc.getTempDirectory(), chunkName+"_MisMatch.bed.gz"));

			//build and execute samtools call
			cmd = bcc.getSamtools()+" mpileup -B -Q "+ minBaseQuality +" -q "+ minMapQuality +" -d 1000000 -f "+bcc.getFasta()+" -l "+tempBed+" "+bcc.getBamNames();
			ProcessBuilder pb = new ProcessBuilder(Misc.WHITESPACE.split(cmd));

			Process proc = pb.start();

			BufferedReader data = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			String line;
			ParsedSample[] parsedSamples = null;

			int counter = 0;
			while ((line=data.readLine())!= null){
				if (line.startsWith("#")){
					System.out.println("\t"+line);
					continue;
				}
				numMpileupLinesProc++;
				
				//parse line
				MpileupLine ml = new MpileupLine(line, minBaseQuality);
				String chr = ml.getChr();
				if (chr == null || ml.getRef().equals("N")) continue;
				boolean chrX = (chr.equals("X") || chr.equals("chrX"));
				
				//load common snv interval tree?
				if (chr.equals(currCommonChr) == false) loadIntervalTree(chr);

				//get samples
				MpileupSample[] samples = ml.getSamples();
				if (samples == null) continue;

				//first set?
				if (afHist == null) {
					makeCounters(samples.length);
					parsedSamples = new ParsedSample[samples.length];
				}

				//check for proper number of samples
				if (samples.length != afHist.length) continue;

				//parse and check samples, this increments the histograms
				boolean passDPIndel = false;
				boolean passAF = false;
				for (int i=0; i<samples.length; i++) {
					parsedSamples[i] = new ParsedSample(samples[i], chrX, i);
					if (parsedSamples[i].passDpIndel) {
						passDPIndel = true;
						if (parsedSamples[i].maxAFIndex!=null && parsedSamples[i].maxAFIndex[0] >= minAFForHom) passAF = true;
					}
				}
				
				//any that pass?
				if (passDPIndel == true && passAF == true){
					//common?
					if (commonSnp(ml.getZeroPos()) == false) {
						//contrast samples
						boolean foundMisMatch = false;
						for (int i=0; i< similarities.length; i++) {
							if (similarities[i].contrast(parsedSamples)) foundMisMatch = true;
						}
						//save it?
						if (foundMisMatch) misMatchBed.println(ml.getBed());
					}
				}
			}

			//finish
			data.close();
			misMatchBed.close();
			if (commonSnvBed != null) tabixReader.close();
			if (numMpileupLinesProc == 0) throw new IOException("Failed to parse any mpileup lines? "+chunkName);
			complete = true;
			failed = false;
			System.out.println(chunkName+ " job complete. "+numMpileupLinesProc+" bases parsed.");
			
			
		} catch (Exception e) {
			System.err.println("Problem executing -> "+cmd);
			e.printStackTrace();
			complete = false;
			failed = true;
		} 
	}
	
	private boolean commonSnp(int zeroPos) {
		if (commonSnvRegions == null) return false;
		Interval1D interval = new Interval1D (zeroPos, zeroPos);
		if (commonSnvRegions.contains(interval)) return true;
		return false;
	}

	class CommonVarChecker {
		boolean checkIfCommon;
		boolean common;
	}
	
	private boolean loadIntervalTree(String chr) {
		try {
			currCommonChr = chr;
			if (commonSnvBed == null) return false;
			commonSnvRegions = new IntervalST<String>();
			TabixReader.Iterator it = tabixReader.query(chr);
			String hit;
			while ((hit = it.next()) != null) {
				String[] t = Misc.TAB.split(hit);
				int start = Integer.parseInt(t[1]);
				int stop = Integer.parseInt(t[2]);
				commonSnvRegions.put(new Interval1D(start, stop-1), chr);
			}
			return true;
		} catch (Exception e){
			return false;
		}
	}

	//a helper classes
	class ParsedSample {
		boolean passDpIndel;
		MpileupSample ms;
		double[] maxAFIndex = null;

		ParsedSample (MpileupSample ms, boolean chrX, int index){
			this.ms = ms;

			//passing?
			if (ms.isPass() == false) passDpIndel = false;
			else{
				//check snv DP
				int snvDp = ms.getReadCoverageForwardBases() + ms.getReadCoverageReverseBases();
				if (snvDp < minSnvDP) passDpIndel = false;
				else {
					//check indel freq
					double indelRto = 1.0 - ((double)snvDp/(double)ms.getReadCoverageAll());
					if (indelRto > maxIndel) passDpIndel = false;
					else {
						passDpIndel = true;
						maxAFIndex = ms.findMaxSnvAFAndIndex();
						if (maxAFIndex[0] >= minAFForHis){
							afHist[index].count(maxAFIndex[0]);
							if (chrX) chrXAfHist[index].count(maxAFIndex[0]);
						}
					}
				}
			}
		}
	}
	
	public void makeCounters(int numSamples) throws Exception{
		int numBams = bcc.getBamFiles().length;
		if (numBams != numSamples) {
			this.failed = true;
			this.complete = false;
			throw new Exception("\nERROR: the number of samples in the mpileup file ("+numSamples+") don't equal the number of bams! "+numBams);
		}
		
		//make histograms
		afHist = new Histogram[numSamples];
		chrXAfHist = new Histogram[numSamples];
		for (int i=0; i< numSamples; i++){
			afHist[i] = new Histogram(minAFForHis,1.01,25);
			chrXAfHist[i] = new Histogram(minAFForHis,1.01,25);
		}

		//make Sims and PCs
		ArrayList<Similarity> al = new ArrayList<Similarity>();
		for (int i=0; i< numSamples; i++){
			for (int j=i+1; j< numSamples; j++){
				al.add(new Similarity(i,j,bcc));
			}
		}
		similarities = new Similarity[al.size()];
		al.toArray(similarities);
	}

	
	//getters and setters
	public String getCmd() {
		return cmd;
	}
	public boolean isComplete() {
		return complete;
	}
	public boolean isFailed() {
		return failed;
	}
	public File getTempBed() {
		return tempBed;
	}

	public Histogram[] getAfHist() {
		return afHist;
	}

	public Histogram[] getChrXAfHist() {
		return chrXAfHist;
	}

	public Similarity[] getSimilarities() {
		return similarities;
	}

	public Gzipper getMisMatchBed() {
		return misMatchBed;
	}

}
