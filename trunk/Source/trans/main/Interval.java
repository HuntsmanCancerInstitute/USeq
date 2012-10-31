package trans.main;
import java.io.*;

import util.gen.*;

/**
 * Container for holding information about a potential enriched region: windows, oligos, sub windows, binding peaks....
 */
public class Interval implements Comparable, Serializable{

	//fields
	private String chromosome;
	private int start1stOligo;		//bp start of first oligo
	private int startLastOligo;		//bp start of last oligo, add oligo length if you want total region
	private int sizeOfOligoMinusOne = 24;		
	private int numberOfWindows;
	private Window bestWindow;	
	private SubWindow bestSubWindow;	
	private BindingPeak[] bindingPeaks = null;
	private double sortBy;				
	
	//from LoadIntervalOligoInfo.java
	private Oligo[] oligos = null;
	private File[] celFiles;
	private int numberTreatmentIntensities;
	private int numberControlIntensities;
	
	//from ScoreIntervals.java
	private String sequence; 			//genomic sequence
	private int numberMotifHits;
	private double[] baseScores;
	private int maxCluster; 			//maximum number of hits within a given clusterWindow
	private boolean bestWindowScored;
	
	//from IntervalFilter.java
	private double[] fractionIntersections;
	private String[] regionNames;
	
	public Interval (Window window, int sizeOfOligo){
		chromosome = window.getChromosome();
		start1stOligo = window.getStart1stOligo();
		startLastOligo = window.getStartLastOligo();
		bestWindow = window;
		sizeOfOligoMinusOne = sizeOfOligo -1;
	}
	
	public int compareTo(Object obj){
		Interval other = (Interval)obj;
		//by score
		if (other.sortBy>sortBy) return 1;
		if (other.sortBy<sortBy) return -1;
		return 0;
	}
	
	/**Sums the logged motif hit values using the best AID window or the entire interval. 
	 * This is equivalent to the product of the probabilities.*/
	public double getMotifHitSum(){
		if (baseScores==null) return -1;
		double sum = 0;
		int start;
		int end;
		//sum just the best AID window
		if (bestWindowScored) {
			start = bestWindow.getStart1stOligo()-start1stOligo;
			end = bestWindow.getStartLastOligo()+sizeOfOligoMinusOne-start1stOligo;
		}
		//sum the entire interval
		else {
			start =0;
			end = baseScores.length;
		}
		for (int i=start; i< end; i++) 
			if (baseScores[i]!=0) sum += baseScores[i];
		return sum;
	}

	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start1stOligo);
		sb.append("\t");
		sb.append(startLastOligo+sizeOfOligoMinusOne);
		sb.append("\t");
		sb.append(numberOfWindows);
		sb.append("\t");
		sb.append(bestWindow.stringRep(sizeOfOligoMinusOne));

		if (oligos!=null){
			//print cel files
			sb.append("\t");
			sb.append(celFiles[0].getName());
			int num = celFiles.length;
			for (int i=1; i< num; i++){
				sb.append(",");
				sb.append(celFiles[0].getName());
			}
			sb.append("\t");
			//number treatments
			sb.append(numberTreatmentIntensities);
			sb.append("\t");
			//number controls
			sb.append(numberControlIntensities);
			sb.append("\t");
			//oligos
			num = oligos.length;
			for (int i=0; i< num; i++){
				
			}
		}
		if (sequence!=null){
			sb.append("\t");
			sb.append(maxCluster);
			sb.append("\t");
			sb.append(numberMotifHits);
			sb.append("\t");
			sb.append(sequence);
		}
		return sb.toString();
	}
	
	/**Given a bp start and bp stop as well as a larger Oligo[] returns a sub array of the oligos flanking and including
	 * the start and stop.*/
	public Oligo[] extractSubRegionOligos(int start, int end){
		if (oligos == null) return null;
		int firstOligo=0;
		int lastOligo=0;
		int numOligos = oligos.length;
		//find first
		for (int i=0; i<numOligos; i++){
			if (oligos[i].getStart() == start ){
				firstOligo = i;
				break;
			}
		}
		//find last
		for (int i=firstOligo; i<numOligos; i++){
			if (oligos[i].getStart() == end ){
				lastOligo = i;
				break;
			}
		}
		//make new array
		Oligo[] subOligos = new Oligo[1+lastOligo-firstOligo];
		int last = lastOligo+1;
		int j=0;
		for (int i=firstOligo; i<last; i++){
			subOligos[j++]=oligos[i];
		}
		return subOligos;	
	}
	public String getBestWindowSequence(){
		if (Misc.isEmpty(sequence)) return "";
		//return seq.substring(subStart-intStart,1+subEnd-intStart);		
		return new String(sequence.substring(bestWindow.getStart1stOligo() - start1stOligo, 
				bestWindow.getStartLastOligo()+ sizeOfOligoMinusOne -start1stOligo));
	}
	public String getBestSubWindowSequence(){
		if (Misc.isEmpty(sequence) || bestSubWindow ==null) return "";
		Oligo[] subs = bestSubWindow.getOligos();
		return new String(sequence.substring(subs[0].getStart() - start1stOligo, 
				subs[subs.length-1].getStart()+ sizeOfOligoMinusOne -start1stOligo));		
	}
	public String getChromosome() {
		return chromosome;
	}
	public int getStartLastOligo() {
		return startLastOligo;
	}
	public int getStart1stOligo() {
		return start1stOligo;
	}
	public void setChromosome(String string) {
		chromosome = string;
	}
	public void setStartLastOligo(int i) {
		startLastOligo = i;
	}
	public void setStart1stOligo(int i) {
		start1stOligo = i;
	}
	public int getNumberMotifHits() {
		return numberMotifHits;
	}
	public void setNumberMotifHits(int motifHits) {
		this.numberMotifHits = motifHits;
	}
	public String getSequence() {
		return sequence;
	}
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	public int getMaxCluster() {
		return maxCluster;
	}
	public void setMaxCluster(int maxCluster) {
		this.maxCluster = maxCluster;
	}
	public void setSortBy(double sortBy) {
		this.sortBy = sortBy;
	}
	public double getSortBy() {
		return sortBy;
	}
	public int getNumberOfWindows() {
		return numberOfWindows;
	}
	public void setNumberOfWindows(int numberOfWindows) {
		this.numberOfWindows = numberOfWindows;
	}
	public Oligo[] getOligos() {
		return oligos;
	}
	public void setOligos(Oligo[] oligos) {
		this.oligos = oligos;
	}
	public Window getBestWindow() {
		return bestWindow;
	}
	public void setBestWindow(Window bestWindow) {
		this.bestWindow = bestWindow;
	}
	public File[] getCelFiles() {
		return celFiles;
	}
	public void setCelFiles(File[] celFiles) {
		this.celFiles = celFiles;
	}
	public int getSizeOfOligoMinusOne() {
		return sizeOfOligoMinusOne;
	}
	public int getNumberControlIntensities() {
		return numberControlIntensities;
	}
	public void setNumberControlIntensities(int numberControlIntensities) {
		this.numberControlIntensities = numberControlIntensities;
	}
	public int getNumberTreatmentIntensities() {
		return numberTreatmentIntensities;
	}
	public void setNumberTreatmentIntensities(int numberTreatmentIntensities) {
		this.numberTreatmentIntensities = numberTreatmentIntensities;
	}
	public double[] getBaseScores() {
		return baseScores;
	}
	public void setBaseScores(double[] baseScores) {
		this.baseScores = baseScores;
	}
	public SubWindow getBestSubWindow() {
		return bestSubWindow;
	}
	public void setBestSubWindow(SubWindow bestSubWindow) {
		this.bestSubWindow = bestSubWindow;
	}
	public boolean isBestWindowScored() {
		return bestWindowScored;
	}
	public void setBestWindowScored(boolean bestAIDWindowScored) {
		this.bestWindowScored = bestAIDWindowScored;
	}
	public BindingPeak[] getBindingPeaks() {
		return bindingPeaks;
	}
	public void setBindingPeaks(BindingPeak[] bindingPeaks) {
		this.bindingPeaks = bindingPeaks;
	}

	public double[] getFractionIntersections() {
		return fractionIntersections;
	}

	public void setFractionIntersections(double[] percentIntersections) {
		this.fractionIntersections = percentIntersections;
	}

	public String[] getRegionNames() {
		return regionNames;
	}

	public void setRegionNames(String[] repeatNames) {
		this.regionNames = repeatNames;
	}

	public void setSizeOfOligoMinusOne(int sizeOfOligoMinusOne) {
		this.sizeOfOligoMinusOne = sizeOfOligoMinusOne;
	}
}
