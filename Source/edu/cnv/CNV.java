package edu.cnv;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import trans.misc.GrGraph;
import trans.roc.Gr;
import util.gen.Misc;
import util.gen.Num;
import java.io.*;

public class CNV implements Serializable{
	
	public ArrayList<Gr> grs = new ArrayList<Gr>();
	public GrGraph grGraph = new GrGraph();
	public float medianScore;
	public int start;
	public int stop;	//not included

	public CNV(String chrom){
		grGraph.setChromosome(chrom);
	}

	public void loadGrGraph(){
		Gr[] g = new Gr[grs.size()];
		grs.toArray(g);
		grGraph.loadGrArray(g);
		float[] scores = grGraph.getValues();
		float[] copyScores = new float[scores.length];
		System.arraycopy(scores, 0, copyScores, 0, scores.length);
		Arrays.sort(copyScores);
		medianScore = (float) Num.median(copyScores);
		start = g[0].getPosition();
		stop = g[g.length-1].getPosition()+1;
	}
	
	public String toString(){
		int[] bp = grGraph.getBasePositions();
		StringBuilder sb = new StringBuilder(grGraph.getChromosome());
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(stop);
		sb.append("\t");
		sb.append(stop-start);
		sb.append("\t");
		sb.append(bp.length);
		sb.append("\t");
		sb.append(medianScore);
		sb.append("\t");
		sb.append(Misc.intArrayToString(grGraph.getBasePositions(), ","));
		sb.append("\t");
		sb.append(Misc.floatArrayToString(grGraph.getValues(), ","));
		return sb.toString();
	}
	
	public String toStringBed(String name){
		StringBuilder sb = new StringBuilder(grGraph.getChromosome());
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(stop);
		sb.append("\t");
		sb.append(name);
		sb.append("\t");
		sb.append(medianScore);
		return sb.toString();
	}
	
	/**Splits a CNV[] by chromosome. Assumes CNV[] is sorted by position.*/
	public static HashMap<String, CNV[]> splitByChromsome(CNV[] cnvs){
		ArrayList<CNV> al = new ArrayList<CNV>();
		HashMap<String, CNV[]> map = new HashMap<String, CNV[]>();
		String chromosome = cnvs[0].getGrGraph().getChromosome();
		//for each GenomicRegion
		for (int i=0; i< cnvs.length; i++){
			//same chromosome?
			if (cnvs[i].getGrGraph().getChromosome().equals(chromosome)){
				al.add(cnvs[i]);
			}
			//different chromosome!
			else {
				//set new chrom in map
				CNV[] regions = new CNV[al.size()];
				al.toArray(regions);
				map.put(chromosome, regions);
				//reset params
				al = new ArrayList();
				al.add(cnvs[i]);
				chromosome = cnvs[i].getGrGraph().getChromosome();
			}
		}
		//add last
		CNV[] regions = new CNV[al.size()];
		al.toArray(regions);
		map.put(chromosome, regions);
		return map;
	}

	/**Returns the Gr slice, assumes stop is included*/
	public Gr[] fetchObservations (int start, int stop){
		ArrayList<Gr> subAL = new ArrayList<Gr>();
		int[] pos = grGraph.getBasePositions();
		float[] val = grGraph.getValues();
		for (int i=0; i< pos.length; i++){
			if (pos[i]>= start && pos[i]<= stop) subAL.add(new Gr(pos[i], val[i]));
		}
		Gr[] slice = new Gr[subAL.size()];
		subAL.toArray(slice);
		return slice;
	}
	
	public ArrayList<Gr> getGrs() {
		return grs;
	}

	public void setGrs(ArrayList<Gr> grs) {
		this.grs = grs;
	}

	public GrGraph getGrGraph() {
		return grGraph;
	}

	public void setGrGraph(GrGraph grGraph) {
		this.grGraph = grGraph;
	}

	public float getMedianScore() {
		return medianScore;
	}

	public void setMedianScore(float medianScore) {
		this.medianScore = medianScore;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}
}