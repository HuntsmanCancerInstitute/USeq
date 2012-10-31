package trans.roc;

import trans.misc.GrGraph;
import java.util.*;
import java.io.*;

/**
 * Gr line container.
 */
public class Gr implements Serializable{
	private int position;
	private float score;
	
	public Gr (String line){
		String[] tokens = line.split("\\s");
		position = Integer.parseInt(tokens[0]);
		score = Float.parseFloat(tokens[1]);
	}
	
	public Gr (int position, float score){
		this.position = position;
		this.score = score;
	}
	
	public static Gr[] parseGrGraph(GrGraph gr){
		int[] basePositions = gr.getBasePositions();
		float[] values = gr.getValues();
		Gr[] grs = new Gr[values.length];
		for (int i=0; i< grs.length; i++){
			grs[i]= new Gr (basePositions[i], values[i]);
		}
		return grs;
		
	}
	
	public int getPosition() {
		return position;
	}
	public float getScore() {
		return score;
	}
	public String toString(){
		return position+"\t"+score;
	}
	public void setScore(float score) {
		this.score = score;
	}
	public void setPosition(int position) {
		this.position = position;
	}
}
