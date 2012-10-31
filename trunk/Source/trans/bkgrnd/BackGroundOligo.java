package trans.bkgrnd;
import java.util.*;
import java.io.*;

import util.bio.calc.*;
import util.gen.*;

/**
 * A control oligo spotted multiple times on the array.
 * */
public class BackGroundOligo implements Serializable {
	
	//fields
	private String sequence;
	private int[][] coordinates; //int[][x,y] , the identical sequence is sometimes spotted many times on the same chip
	private ArrayList intensityArrayLists = new ArrayList(); //ArrayList of float[], one per chip

	//contructors
	public BackGroundOligo(String sequence, int[][] coordinates){
		this.sequence = sequence.toLowerCase();
		this.coordinates = coordinates;
	}
	
	//methods
	public String toString(){
		//(seq gc tm 23,36,19;45,46,41;34,38)
		StringBuffer sb = new StringBuffer();
		sb.append(sequence);
		sb.append("\t");
		sb.append(NucleicAcid.calculateFractionGC(sequence));
		sb.append("\t");
		sb.append(NucleicAcid.nearestNeighborOligoTm(sequence, 0.00000005,0.05));
		sb.append("\t");
		sb.append(Misc.floatArrayToString((float[])intensityArrayLists.get(0),","));
		for (int i=1; i< intensityArrayLists.size(); i++){
			sb.append(";");
			sb.append(Misc.floatArrayToString((float[])intensityArrayLists.get(i),","));
		}
		return sb.toString();
	}
	
	//getter setter methods
	public int[][] getCoordinates() {
		return coordinates;
	}
	public ArrayList getIntensityArrayLists() {
		return intensityArrayLists;
	}
	public void setIntensityArrayLists(ArrayList intensityArrayLists) {
		this.intensityArrayLists = intensityArrayLists;
	}
	public String getSequence() {
		return sequence;
	}
	
	

}
