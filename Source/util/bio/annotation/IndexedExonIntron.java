package util.bio.annotation;

import java.util.ArrayList;

/**For fetching all combinations of transcripts given alternative overlapping exons*/
public class IndexedExonIntron{
	int indexNumber;
	IndexedExonIntron[] indexes;
	ArrayList<ExonIntron> alternativeExons;
	int maxCombinations;
	
	public IndexedExonIntron (int indexNumber, IndexedExonIntron[] indexes, ArrayList<ExonIntron> alternativeExons, int maxCombinations){
		this.indexNumber = indexNumber;
		this.indexes = indexes;
		this.alternativeExons = alternativeExons;
		this.maxCombinations = maxCombinations;
	}
	
	public ArrayList<ArrayList<ExonIntron>> fetchRight(){
		//none to the right?
		if (indexNumber == indexes.length-1) {
			ArrayList<ArrayList<ExonIntron>> toReturn = new ArrayList<ArrayList<ExonIntron>>();
			for (int i=0; i< alternativeExons.size(); i++){
				//create a new AL
				ArrayList<ExonIntron> al = new ArrayList<ExonIntron>();
				//add this
				al.add(alternativeExons.get(i));
				//add toReturn
				toReturn.add(al);
			}
			return toReturn;
		}
		//join with right
		ArrayList<ArrayList<ExonIntron>> rightOverlaps = indexes[indexNumber+1].fetchRight();
		if (rightOverlaps == null || rightOverlaps.size() > maxCombinations) return null;
		
		ArrayList<ArrayList<ExonIntron>> toReturn = new ArrayList<ArrayList<ExonIntron>>();
		//for each one of this overlaps
		for (int i=0; i< alternativeExons.size(); i++){
			for (int j=0; j < rightOverlaps.size(); j++){
				//create a new AL
				ArrayList<ExonIntron> al = new ArrayList<ExonIntron>();
				//add this
				al.add(alternativeExons.get(i));
				//add rights
				al.addAll(rightOverlaps.get(j));
				//add toReturn
				toReturn.add(al);
			}
		}
		return toReturn;
	}

	public ArrayList<ExonIntron> getAlternativeExons() {
		return alternativeExons;
	}
}
