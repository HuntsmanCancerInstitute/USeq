package edu.utah.seq.vcf.cluster;

import java.util.ArrayList;

public class VCFCluster {
	private ArrayList<VCFClusterSample> samples = new ArrayList<VCFClusterSample>();

	public ArrayList<VCFClusterSample> getSamples() {
		return samples;
	}
	
	public void fetchIndexNames(StringBuilder sb){
		sb.append(samples.get(0).getSampleIndex());
		int num = samples.size();
		for (int i=1; i< num; i++){
			sb.append("_");
			sb.append(samples.get(i).getSampleIndex());
		}
	}

	public void fetchSampleNames(StringBuilder sb) {
		sb.append(samples.get(0).getSampleName());
		int num = samples.size();
		for (int i=1; i< num; i++){
			sb.append("_");
			sb.append(samples.get(i).getSampleName());
		}
		
	}
}
