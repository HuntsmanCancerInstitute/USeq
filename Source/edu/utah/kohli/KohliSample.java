package edu.utah.kohli;

import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;
import util.gen.Misc;

public class KohliSample {
	
	private String sampleId = null;
	private String type = null;
	private String dateString = null;
	private boolean excludeFromPon = false;
	private boolean germlineSample;
	
	private HashMap<String, AccuGenProbeCounts> accuGenProbeCounts = new HashMap<String, AccuGenProbeCounts>();
	
	//Read Coverage
	private HashMap<String, Double> readCoverageCounts = new HashMap<String, Double>();
	private HashMap<String, Double> scaledReadCoverageCounts = null;
	private double meanRegionCounts = -1.0;
	private double regionCountScalar = -1.0;
	
	public KohliSample (String sampleId, boolean isGermline, String dateString) {
		this.sampleId = sampleId;
		if (isGermline) germlineSample = true;
		else germlineSample = false;
		this.dateString = dateString;
	}
	
	public AccuGenProbeCounts getAccuGenProbeCount(String probeId) {
		return accuGenProbeCounts.get(probeId);
	}
	
	public String toStringProbes(ArrayList<AccuGenProbe> probes) {
		StringBuilder sb = new StringBuilder();
		sb.append(sampleId+"\t"+type+"\t"+dateString+"\n");
		//probe info
		for (AccuGenProbe p: probes) {
			sb.append("\t");
			AccuGenProbeCounts counts = accuGenProbeCounts.get(p.getOriginalInput());
			sb.append(p.getGene()+":"+p.getPos()+":"+counts.toString());
		}
		sb.append("\n");
		
		return sb.toString();
	}
	
	public String toStringRegions(ArrayList<CaptureRegion> regions) {
		StringBuilder sb = new StringBuilder();
		sb.append(sampleId+"\t"+type+"\t"+dateString+"\n");
		//probe info
		for (CaptureRegion p: regions) {
			sb.append("\t");
			Double count = readCoverageCounts.get(p.getOriginalInput());
			sb.append(p.getOriginalInput()+":"+count);
		}
		sb.append("\n");
		return sb.toString();
	}
	
	public String meanRegionCounts(ArrayList<CaptureRegion> regions) {
		//probe info
		double total = 0;
		for (CaptureRegion p: regions) {
			Double count = readCoverageCounts.get(p.getOriginalInput());
			total+= count;
		}
		return sampleId+"\t"+type+"\t"+total/(double)regions.size();
	}
	
	public double getMeanRegionCounts(ArrayList<CaptureRegion> regions) {
		double total = 0;
		for (CaptureRegion p: regions) {
			Double count = readCoverageCounts.get(p.getOriginalInput());
			total+= count;
		}
		meanRegionCounts = total/(double)regions.size();
		return meanRegionCounts;
	}
	
	public double[] getScaledRegionCountArray(ArrayList<CaptureRegion> regions, double targetMean) {
		regionCountScalar = targetMean / getMeanRegionCounts(regions);
		double[] scaled = new double[regions.size()];
		//probe info
		for (int i=0; i< scaled.length; i++) {
			String regionId = regions.get(i).getOriginalInput();
			Double count = readCoverageCounts.get(regionId);
			scaled[i] = count * regionCountScalar;
		}
		return scaled;
	}
	
	public void createScaledRegionCountMap(ArrayList<CaptureRegion> regionsForScalarCalc,  ArrayList<CaptureRegion> regionsToScale, double targetMean) {
		regionCountScalar = targetMean / getMeanRegionCounts(regionsForScalarCalc);
		//IO.pl("Scaling using just one gene for "+this.sampleId+ " numRegionsInGene "+regionsForScalarCalc.size()+" scalar "+regionCountScalar);
		scaledReadCoverageCounts = new HashMap<String, Double>();
		//probe info	
		for (int i=0; i< regionsToScale.size(); i++) {
			String regionId = regionsToScale.get(i).getOriginalInput();
			Double count = readCoverageCounts.get(regionId);
			scaledReadCoverageCounts.put(regionId, count * regionCountScalar);
		}
	}
	
	
	public String toString() {
		return sampleId+"\t"+type+"\t"+dateString;
	}

	public String getSampleId() {
		return sampleId;
	}

	public String getDateString() {
		return dateString;
	}

	public void addAccuGenProbeCounts(AccuGenProbe p, String vcfSampleString) {
		//0,255,255:3150,152
		String[] callsCounts = Misc.COLON.split(vcfSampleString);
		String[] counts = Misc.COMMA.split(callsCounts[1]);
		int ref = Integer.parseInt(counts[0]);
		int alt = Integer.parseInt(counts[1]);
		accuGenProbeCounts.put(p.getOriginalInput(), new AccuGenProbeCounts(ref,alt));
	}
	
	public void addCaptureRegionCount(CaptureRegion p, String count) {
		readCoverageCounts.put(p.getOriginalInput(), Double.parseDouble(count));
	}

	public HashMap<String, Double> getReadCoverageCounts() {
		return readCoverageCounts;
	}

	public HashMap<String, Double> getScaledReadCoverageCounts() {
		return scaledReadCoverageCounts;
	}

	public boolean isExcludeFromPon() {
		return excludeFromPon;
	}

	public void setExcludeFromPon(boolean excludeFromPon) {
		this.excludeFromPon = excludeFromPon;
	}

	public boolean isGermlineSample() {
		return germlineSample;
	}
	
}
