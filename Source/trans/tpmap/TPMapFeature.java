package trans.tpmap;
import java.io.*;
/**
 * Container for info about a particular oligo feature from a tpmap file.
 */
public class TPMapFeature implements Serializable{
	//fields
	private String sequence;
	private String chromosome;
	private int position;
	private short PMX;
	private short PMY;
	
	private short MMX;
	private short MMY;
	private float intensity;
	
	public TPMapFeature (String[] tokens){
		PMX = Short.parseShort(tokens[4]);
		PMY = Short.parseShort(tokens[5]);
		MMX = Short.parseShort(tokens[6]);
		MMY = Short.parseShort(tokens[7]);
	}
	
	public String toString(){
		return PMX+"\t"+PMY+"\t"+MMX+"\t"+MMY;
	}

	public String getChromosome() {
		return chromosome;
	}
	public float getIntensity() {
		return intensity;
	}
	public short getMMX() {
		return MMX;
	}
	public short getMMY() {
		return MMY;
	}
	public short getPMX() {
		return PMX;
	}
	public short getPMY() {
		return PMY;
	}
	public int getPosition() {
		return position;
	}
	public String getSequence() {
		return sequence;
	}
}
