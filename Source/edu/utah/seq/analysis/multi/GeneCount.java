package edu.utah.seq.analysis.multi;
import java.io.Serializable;

public class GeneCount implements Serializable{
	//fields
	private int count = 0;
	private int[] exonCounts;
	private static final long serialVersionUID = 1L;
	
	public GeneCount(int count, int[] exonCounts){
		this.count = count;
		this.exonCounts = exonCounts;
	}
	
	/**Calculates the reads per kb per million mapped reads 
	 * # Observed reads in the gene/ bp size of the gene / 1000/ total number reads/ 1000000 */
	public double calculateFPKM(double totalReads, double interrogatedRegionBPSize){
		double millionTotalMappedReads = totalReads/1000000.0;
		double exonicBasesPerKB = interrogatedRegionBPSize/1000;
		double rpkm = ((double)count)/exonicBasesPerKB/millionTotalMappedReads;
		return rpkm;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(count);
		for (int ec: exonCounts){
			sb.append("\t");
			sb.append(ec);
		}
		return sb.toString();
	}
	
	public int getCount() {
		return count;
	}
	public void setCount(int count) {
		this.count = count;
	}
	public int[] getExonCounts() {
		return exonCounts;
	}
	public void setExonCounts(int[] exonCounts) {
		this.exonCounts = exonCounts;
	}
	
}
