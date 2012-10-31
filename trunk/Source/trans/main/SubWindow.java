package trans.main;
import java.io.*;

/**
 * Class to hold information about a subWindow.
 */
public class SubWindow implements Serializable{
	//fields
	private Oligo[] oligos;
	private double medianRatio;
	private int index;				//used internally by FindSubBindingRegions, window index
	
	public SubWindow (int index, double medianRatio){
		this.index = index;
		this.medianRatio = medianRatio;
	}
	public Oligo[] getOligos() {
		return oligos;
	}
	public double getMedianRatio() {
		return medianRatio;
	}
	public void setOligos(Oligo[] oligos) {
		this.oligos = oligos;
	}
	public int getIndex() {
		return index;
	}
}
