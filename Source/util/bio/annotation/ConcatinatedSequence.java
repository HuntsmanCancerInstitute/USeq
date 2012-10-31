package util.bio.annotation;

import java.util.ArrayList;

public class ConcatinatedSequence {
	//fields
	private SubSequence[] subSequence;
	private String chromosome;
	private String sequence;
	
	//constructor
	public ConcatinatedSequence (ArrayList<SubSequence> fivePrime, ArrayList<SubSequence> threePrime, String chromosome){
		this.chromosome = chromosome;
		//make array
		int numFivePrime = fivePrime.size();
		int numThreePrime = threePrime.size();
		subSequence = new SubSequence[numFivePrime + numThreePrime];
		for (int i=0; i< numFivePrime; i++) subSequence[i] = fivePrime.get(i);
		int counter =0;
		for (int i=numFivePrime; i<subSequence.length; i++) subSequence[i] = threePrime.get(counter++);
	}
	
	public ConcatinatedSequence (ExonIntron[] ei, String chromosome){
		this.chromosome = chromosome;
		subSequence = new SubSequence[ei.length];
		for (int i=0; i< ei.length; i++) subSequence[i] = new SubSequence (ei[i].getSequence(), ei[i].getStart(), ei[i].getEnd());
	}
	
	//methods
	public String getSequence(){
		if (sequence == null){
			StringBuilder sb = new StringBuilder();
			for (SubSequence ss : subSequence) sb.append(ss.getSequence());
			sequence = sb.toString();
		}
		return sequence;
	}
	
	public String getCoordinateString(){
		StringBuilder sb = new StringBuilder();
		sb.append(subSequence[0].getName());
		for (int i=1; i<subSequence.length; i++) {
			sb.append("_");
			sb.append(subSequence[i].getName());
		}
		return sb.toString();
	}
	
	/**Provide something like 'ENSG0003456' will be converted to '>ENSG0003456:chrX:3456-3556_4568-4689_5000-5100'*/
	public String getFastaEntry (String geneName){
		StringBuilder sb = new StringBuilder(">");
		sb.append(geneName);
		sb.append(":");
		sb.append(chromosome);
		sb.append(":");
		
		sb.append(subSequence[0].getName());
		for (int i=1; i<subSequence.length; i++) {
			sb.append("_");
			sb.append(subSequence[i].getName());
		}
		sb.append("\n");
		sb.append(getSequence());
		return sb.toString();
	}
	
	public boolean equals(Object o){
		String seqFirst = getSequence();
		String seqSecond = ((ConcatinatedSequence)o).getSequence();
		return seqFirst.equals(seqSecond);
	}
	
	public int hashCode(){
		return getSequence().hashCode();
	}

}
