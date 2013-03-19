package edu.utah.seq.vcf;

import java.util.Arrays;


/**Container for base positions and their corresponding vcfRecords.*/
public class VCFLookUp {

	//fields
	private int[] basePosition;
	private VCFRecord[] vcfRecord;
	
	public VCFLookUp(int[] basePosition, VCFRecord[] vcfRecord) {
		this.basePosition = basePosition;
		this.vcfRecord = vcfRecord;
	}
	
	
	
	/**Returns an array of VCFRecord containing the slice defined by the start and stop(excluded). Returns null if none found.*/
	public VCFRecord[] fetchVCFRecords (int startBp, int stopBp){
		int[] indexes = findIndexes (startBp, stopBp);	
		if (indexes == null || indexes[0] == indexes[1]) return null;
		int num = indexes[1] - indexes[0];
		if (num ==0) return null;
		VCFRecord[] vals = new VCFRecord[num];
		int counter = 0;
		for (int i=indexes[0]; i< indexes[1]; i++){
			vals[counter++] = vcfRecord[i];
		}
		return vals ;
	}
	
	/**Returns an array of VCFRecord containing the slice defined by the start and stop(excluded). Returns null if none found.*/
	public VCFRecord[] fetchVCFRecordsDebug (int startBp, int stopBp){
		System.out.println("Looking for "+startBp+" "+stopBp);
		int[] indexes = findIndexes (startBp, stopBp);	
		if (indexes == null || indexes[0] == indexes[1]) return null;
		int num = indexes[1] - indexes[0];
		if (num ==0) return null;
		VCFRecord[] vals = new VCFRecord[num];
		int counter = 0;
		for (int i=indexes[0]; i< indexes[1]; i++){
			vals[counter++] = vcfRecord[i];
			System.out.println("\tFound vcf "+vcfRecord[i].toString());
		}
		return vals ;
	}

	/**Given a start bp (included) and stop bp (not included), returns start (included) and stop (not included) indexes.
	 * May return startIndex = endIndex, therefore nothing found.*/
	public int[] findIndexes(int startBp, int stopBp){
		//find start index, included
		int startIndex = Arrays.binarySearch(basePosition, startBp);

		if (startIndex < 0) {
			startIndex = (startIndex*-1) -1;
		}
		else {
			//find first instance of startBp
			while (true){
				int minOne = startIndex - 1;
				if (minOne < 0) break;
				if (basePosition[minOne] != startBp) break;
				startIndex = minOne;
			}
		}
		//find stop index, not included
		int stopBPMinOne = stopBp-1;
		int endIndex = Arrays.binarySearch(basePosition, stopBPMinOne);		
		if (endIndex < 0) {
			endIndex = (endIndex*-1) -1;

		}
		else {
			//find last instance of endBp
			while (true){
				int addOne = endIndex +1;
				if (addOne >= basePosition.length) break;
				if (basePosition[addOne] != stopBPMinOne) break;
				endIndex = addOne;
			}
			//add one to stop index, it's not included
			endIndex++;
		}
		return new int[]{startIndex, endIndex};
	}

	public int[] getBasePosition() {
		return basePosition;
	}

	public VCFRecord[] getVcfRecord() {
		return vcfRecord;
	}
}
