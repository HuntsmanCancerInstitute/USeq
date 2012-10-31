package edu.utah.seq.parsers;

import util.gen.Misc;


public class AlignmentPair {
	
	private ParsedAlignment first;
	private ParsedAlignment second;
	private BaseObservation[] baseObservations;
	
	public AlignmentPair (ParsedAlignment first, ParsedAlignment second){
		this.first = first;
		this.second = second;
		//combine baseObservations
		BaseObservation[] boFirst = first.getBo();
		BaseObservation[] boSecond = second.getBo();
		baseObservations = new BaseObservation[boFirst.length + boSecond.length];
		System.arraycopy(boFirst, 0, baseObservations, 0, boFirst.length);
		System.arraycopy(boSecond, 0, baseObservations, boFirst.length, boSecond.length);
	}
	
	public AlignmentPair (ParsedAlignment first){
		this.first = first;
		baseObservations = first.getBo();
	}
	
	public AlignmentPair (ParsedAlignment first, ParsedAlignment second, BaseObservation[] baseObservations){
		this.first = first;
		this.second = second;
		this.baseObservations = baseObservations;
	}
	/**@return int[]{numberNonConvertedCs, numberConvertedCs}*/
	public int[] getCCounts(){
		int numCon = 0;
		int numNonCon = 0;
		for (BaseObservation bo: baseObservations) {
			if (bo.isConverted()) numCon++;
			else numNonCon++;
		}
		return new int[]{numNonCon, numCon};
	}
	/**@return int[]{numberNonConvertedmCGs, numberConvertedmCGs}*/
	public int[] getmCGCounts(int start, int stop){
		int numCon = 0;
		int numNonCon = 0;
		for (BaseObservation bo: baseObservations) {
			int centerBase = bo.getPosition();
			if (centerBase < start || centerBase >= stop) continue;
			if (bo.getGenSeq().charAt(3) == 'G'){
				if (bo.isConverted()) numCon++;
				else numNonCon++;
			}
		}
		return new int[]{numNonCon, numCon};
	}
	public ParsedAlignment getFirst() {
		return first;
	}

	public ParsedAlignment getSecond() {
		return second;
	}

	public BaseObservation[] getBaseObservations() {
		return baseObservations;
	}
}
