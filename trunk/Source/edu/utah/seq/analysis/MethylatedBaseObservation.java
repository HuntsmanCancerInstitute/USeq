package edu.utah.seq.analysis;

import util.gen.Num;

/**Container for holding T vs C methylation observations at a given base position.*/
public class MethylatedBaseObservation {
	
		//fields
		private int position;
		private float nonConT;
		private float conT;
		private float nonConC;
		private float conC;
		
		//constructor
		public MethylatedBaseObservation (int position, float nonConT, float conT, float nonConC, float conC){
			this.position = position;
			this.nonConT = nonConT;
			this.conT = conT;
			this.nonConC = nonConC;
			this.conC = conC;
		}
		
		//methods
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(position +"\tposition\n");
			sb.append(nonConT +"\tnonConT\n");
			sb.append(conT +"\tconT\n");
			sb.append(nonConC +"\tnonConC\n");
			sb.append(conC +"\tconC\n");
			sb.append(((nonConT+1)/(nonConT+conT+2)) +"\ttMethylated\n");
			sb.append(((nonConC+1)/(nonConC+conC+2)) +"\tcMethylated\n");
			sb.append(getDifferentialFractionMethylation() +"\tfraction\n");
			return sb.toString();
		}
		
		public float getFractionMethylatedT(){
			float tMethylated = (nonConT+1)/(nonConT+conT+2);
			return tMethylated;
		}
		
		public float getFractionMethylatedC(){
			float cMethylated = (nonConC+1)/(nonConC+conC+2);
			return cMethylated;
		}
		
		public float getDifferentialFractionMethylation() {
			float tMethylated = (nonConT+1)/(nonConT+conT+2);
			float cMethylated = (nonConC+1)/(nonConC+conC+2);
			return tMethylated/cMethylated;
		}
		
		/**returns nonConT, conT, nonConC, conC*/
		public long[] getCounts() {
			long[] counts = new long[4];
			counts[0] = (long)nonConT;
			counts[1] = (long)conT;
			counts[2] = (long)nonConC;
			counts[3] = (long)conC;
			return counts;
		}
		
		public static int[] fetchPositions (MethylatedBaseObservation[] mbo){
			int[] positions = new int[mbo.length];
			for (int i=0; i< mbo.length; i++) positions[i] = mbo[i].position;
			return positions;
		}
		
		public static float[] fetchLog2DifferentialFractionMethylation (MethylatedBaseObservation[] mbo){
			float[] f = new float[mbo.length];
			for (int i=0; i< mbo.length; i++) f[i] = Num.log2(mbo[i].getDifferentialFractionMethylation());
			return f;
		}

		public int getPosition() {
			return position;
		}

		public void setPosition(int position) {
			this.position = position;
		}

		public float getNonConT() {
			return nonConT;
		}

		public void setNonConT(float nonConT) {
			this.nonConT = nonConT;
		}

		public float getConT() {
			return conT;
		}

		public void setConT(float conT) {
			this.conT = conT;
		}

		public float getNonConC() {
			return nonConC;
		}

		public void setNonConC(float nonConC) {
			this.nonConC = nonConC;
		}

		public float getConC() {
			return conC;
		}

		public void setConC(float conC) {
			this.conC = conC;
		}

	}