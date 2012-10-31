package trans.tpmap;
import java.io.*;

import util.gen.IO;

/**
 * Minimal info about a feature for combining multiple chips together.
 */
public class IntensityFeature implements Comparable, Serializable{
	//fields
	public static final long serialVersionUID = 1;
	public int start;
	public float intensity;
	
	//constructors
	public IntensityFeature(MapFeature feature, float intensity){
		this.start = feature.start;
		this.intensity = intensity;
	}
	
	//main methods
	/**Sorts by  position.*/
	public int compareTo(Object obj){
		IntensityFeature other = (IntensityFeature)obj;
		//sort by base position
		if (other.start>start) return -1;
		if (other.start<start) return 1;
		return 0;
	}
}
