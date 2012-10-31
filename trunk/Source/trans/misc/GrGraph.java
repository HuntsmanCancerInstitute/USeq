package trans.misc;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

import trans.roc.*;

/**Container for holding info related to a xxx.gr file.*/
public class GrGraph implements Serializable{

	private String chromosome;
	private String genomeVersion;
	private int[] basePositions;
	private float[] values;
	
	public GrGraph (String chromosome, String genomeVersion, int[] basePositions, float[] values){
		this.chromosome = chromosome;
		this.genomeVersion = genomeVersion;
		this.basePositions = basePositions;
		this.values = values;
	}
	
	public GrGraph (String chromosome, String genomeVersion, Gr[] grs){
		this.chromosome = chromosome;
		this.genomeVersion = genomeVersion;
		loadGrArray(grs);
	}
	
	public GrGraph(){}
	
	public void loadGrArray(Gr[] grs){
		basePositions = new int[grs.length];
		values = new float[grs.length];
		for (int i=0; i< grs.length; i++){
			basePositions[i] = grs[i].getPosition();
			values[i]= grs[i].getScore();
		}
	}
	
	/**Loads an sgr file returning chromosome specific GrGraphs.*/
	public static GrGraph[] loadSgrFile(File sgrFile){
		//load
		ArrayList sgrAL = Sgr.loadSgrFile(sgrFile);
		//convert to Sgr[]
		Sgr[] sgr = new Sgr[sgrAL.size()];
		sgrAL.toArray(sgr);
		//sort
		Arrays.sort(sgr, new SgrComparator());
		//split by chromosome
		GrGraph[] grs = Sgr.splitByChrom(sgr);
		sgr = null;
		return grs;
	}
	
	/**Writes to the file base tab value.  
	 * It's strongly recommended to use the chromosome and genomeVersion in the file text!*/
	public boolean writeGrFile(File name){
		try{
			PrintWriter out = new PrintWriter (new FileWriter(name));
			for (int i=0; i<basePositions.length; i++){
				out.print(basePositions[i]);
				out.print("\t");
				out.println(values[i]);
			}
			out.close();
			return true;
		} catch (IOException e){
			e.printStackTrace();
		}
		return false;
	}

	public int[] getBasePositions() {
		return basePositions;
	}

	public void setBasePositions(int[] basePositions) {
		this.basePositions = basePositions;
	}

	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public String getGenomeVersion() {
		return genomeVersion;
	}

	public void setGenomeVersion(String genomeVersion) {
		this.genomeVersion = genomeVersion;
	}

	public float[] getValues() {
		return values;
	}

	public void setValues(float[] values) {
		this.values = values;
	}

}
