package edu.utah.seq.parsers.mpileup.concordance;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.parsers.jpileup.BamPileupTabixLoaderSingle;
import edu.utah.seq.parsers.jpileup.BaseCount;
import edu.utah.seq.parsers.jpileup.BpileupLine;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;

public class ConcordanceChunkBamPileup implements Runnable{

	//fields
	private boolean complete = false;
	private boolean failed = false;
	private SampleConcordance bcc;
	private Bed[] regions;
	private String chunkName;
	private String cmd = null;
	private File commonSnvBed;
	private Gzipper misMatchBed;
	private int minSnvDP;
	private double minAFForHom;
	private double minAFForHis;
	private double maxIndel;
	private Histogram[] afHist;
	private Histogram[] chrXAfHist;
	private SimilarityBamPileup[] similarities;
	private long numMpileupLinesProc = 0;
	private IntervalST<String> commonSnvRegions = null;
	private String currCommonChr = "";
	private TabixReader tabixReader;
//private int numSimCheck = 0;

	public ConcordanceChunkBamPileup(Bed[] regions, SampleConcordance bcc, String chunkName) throws Exception{
		this.bcc = bcc;
		this.regions = regions;
		this.chunkName = chunkName;
		this.minSnvDP = bcc.getMinSnvDP();
		this.minAFForHom = bcc.getMinAFForHom();
		this.minAFForHis = bcc.getMinAFForHis();
		this.maxIndel = bcc.getMaxIndel();
		this.commonSnvBed = bcc.getCommonSnvBed();
		if (commonSnvBed != null) tabixReader = new TabixReader(commonSnvBed.toString());
	}

	public void run(){
		try {	
			
			//make gzipper
			misMatchBed = new Gzipper(new File (bcc.getTempDirectory(), chunkName+"_MisMatch.bed.gz"));

			//make a bam pileup loader
			BamPileupTabixLoaderSingle loader = new BamPileupTabixLoaderSingle(bcc.getBamPileupFile(), 0);
			
			ParsedSample[] parsedSamples = null;
			
			//for each region
//boolean printMe = false;
			for (Bed region: regions) {
				ArrayList<BpileupLine> lines = loader.fetchBpileupRecords(region);
				
				//for each line
				for (BpileupLine ml: lines) {
					numMpileupLinesProc++;
//printMe = false;
					String chr = ml.getChr();
					if (chr == null || ml.getRef()=='N') continue;
					boolean chrX = (chr.equals("X") || chr.equals("chrX"));
					
//if (chr.equals("chr1") && ml.getZeroPos() > 146983121 && ml.getZeroPos()< 146983141) printMe=true;
//if (printMe) System.out.println(ml.getLine());

					//load common snv interval tree?
					if (chr.equals(currCommonChr) == false) loadIntervalTree(chr);

					//get samples
					BaseCount[] samples = ml.getSamples();
					if (samples == null) continue;
					
					//first set?
					if (afHist == null) {
						makeCounters(samples.length);
						parsedSamples = new ParsedSample[samples.length];
					}

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
					
//if (printMe) System.out.println(passDPIndel+" "+passAF);

					//any ParsedSamples that pass?
					if (passDPIndel == true && passAF == true){
						//common?
						if (commonSnp(ml.getZeroPos()) == false) {
							//contrast samples
							boolean foundMisMatch = false;
//numSimCheck++;							
							for (int i=0; i< similarities.length; i++) {
								if (similarities[i].contrast(parsedSamples)) foundMisMatch = true;
							}
							//save it?
							if (foundMisMatch) misMatchBed.println(ml.getBed());
						}
					}
				}
			}

			//finish
			loader.getTabixReader().close();
			misMatchBed.close();
			if (commonSnvBed != null) tabixReader.close();
			if (numMpileupLinesProc == 0) IO.pl("WARNING: Failed to parse any mpileup lines? "+chunkName);
			complete = true;
			failed = false;
System.out.println("\t"+chunkName+ " job complete. "+numMpileupLinesProc+" bases parsed.");
//System.out.println(chunkName+ " job complete. "+numSimCheck+" sim check.");


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
		boolean passDpIndel = false;
		BaseCount ms;
		double[] maxAFIndex = null;

		ParsedSample (BaseCount ms, boolean chrX, int index){
			this.ms = ms;

			//check snv DP
			int snvDp = (int)ms.getPassingReadCoverageSnv(); 
			if (snvDp < minSnvDP) passDpIndel = false;
			else {
				//check indel freq
				double indelRto = ms.getAlleleFreqIndel();
				if (indelRto > maxIndel) passDpIndel = false;
				else {
					passDpIndel = true;
					//might be null!
					maxAFIndex = ms.findMaxNonReferenceSnvAFAndIndex();					
					if (maxAFIndex != null && maxAFIndex[0] >= minAFForHis){
						afHist[index].count(maxAFIndex[0]);
						if (chrX) chrXAfHist[index].count(maxAFIndex[0]);
					}
				}
			}

		}
	}

	public void makeCounters(int numSamples) throws Exception{

		//make histograms
		afHist = new Histogram[numSamples];
		chrXAfHist = new Histogram[numSamples];
		for (int i=0; i< numSamples; i++){
			afHist[i] = new Histogram(minAFForHis,1.01,25);
			chrXAfHist[i] = new Histogram(minAFForHis,1.01,25);
		}

		//make Sims and PCs
		ArrayList<SimilarityBamPileup> al = new ArrayList<SimilarityBamPileup>();
		for (int i=0; i< numSamples; i++){
			for (int j=i+1; j< numSamples; j++){
				al.add(new SimilarityBamPileup(i,j,bcc));
			}
		}
		similarities = new SimilarityBamPileup[al.size()];
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

	public Histogram[] getAfHist() {
		return afHist;
	}

	public Histogram[] getChrXAfHist() {
		return chrXAfHist;
	}

	public SimilarityBamPileup[] getSimilarities() {
		return similarities;
	}

	public Gzipper getMisMatchBed() {
		return misMatchBed;
	}

}
