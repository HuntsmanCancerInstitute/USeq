package edu.utah.seq.analysis;

import java.util.ArrayList;
import java.util.Arrays;

import edu.utah.seq.data.PointData;
import util.gen.Misc;
import util.gen.Num;

/**Container for holding methylation observations at a given base position.*/
public class MethylatedBaseObservationOneSample {
	
		//fields
		private int position;
		private float nonCon;
		private float con;
		
		
		//constructor
		public MethylatedBaseObservationOneSample (int position, float nonCon, float con){
			this.position = position;
			this.nonCon = nonCon;
			this.con = con;
			
		}
		
		//methods
		public void subsample(int maxReadCoverage) {
			int total = (int)(nonCon+con);
			if (total > maxReadCoverage ){
				boolean[] b = new boolean[total];
				Arrays.fill(b, 0, (int)con, true);
				Misc.randomize(b, 0l);
				nonCon = 0f;
				con = 0f;
				for (int i=0; i< maxReadCoverage; i++){
					if (b[i]) con++;
					else nonCon++;
				}
			}

		}
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(position +"\tposition\n");
			sb.append(nonCon +"\tnonCon\n");
			sb.append(con +"\tcon\n");
			sb.append(((nonCon+1)/(nonCon+con+2)) +"\tMethylated\n");
			
			return sb.toString();
		}
		
		public float getFractionMethylation() {
			return(nonCon+1)/(nonCon+con+2);
		}
		
		public float getFractionMethylationNoAddOne() {
			return(nonCon)/(nonCon+con);
		}
		
		/**returns nonCon, con*/
		public long[] getCounts() {
			long[] counts = new long[2];
			counts[0] = (long)nonCon;
			counts[1] = (long)con;
			return counts;
		}
		
		public static int[] fetchPositions (MethylatedBaseObservationOneSample[] mbo){
			int[] positions = new int[mbo.length];
			for (int i=0; i< mbo.length; i++) positions[i] = mbo[i].position;
			return positions;
		}
		
		/**Collects bases that have a minimum interrogation. Returns null if none found.*/
		public static MethylatedBaseObservationOneSample[] fetchCommonBasesWithMinimumObservations(PointData nonCon, PointData con, int minimumReadCoverage){

			//fetch arrays
			int[] positionsNonCon = nonCon.getPositions();
			float[] readsNonCon = nonCon.getScores();
			int[] positionsCon = con.getPositions();
			float[] readsCon = con.getScores();

			//make containers for graph data
			ArrayList<MethylatedBaseObservationOneSample> ibAL = new ArrayList<MethylatedBaseObservationOneSample>();

			//collect all positions
			int[] allPositions = Num.returnUniques(new int[][]{positionsNonCon, positionsCon});

			//for each position 
			int indexNonCon =0;
			int indexCon =0;
			for (int i=0; i< allPositions.length; i++){
				int testPos = allPositions[i];
				//values for each
				float numNonCon =0;
				float numCon =0;

				//treatment
				//present in nonConT?
				for (int j=indexNonCon; j < positionsNonCon.length; j++){
					//match!
					if (testPos == positionsNonCon[j]){
						numNonCon = readsNonCon[j];
						indexNonCon++;
						break;
					}
					//less than
					if (testPos < positionsNonCon[j]) break;
					//greater than so keep advancing
					indexNonCon = j;
				}
				//present in conT?
				for (int j=indexCon; j < positionsCon.length; j++){
					//match!
					if (testPos == positionsCon[j]){
						numCon = readsCon[j];
						indexCon++;
						break;
					}
					//less than
					if (testPos < positionsCon[j]) break;
					//greater than so keep advancing
					indexCon = j;
				}
				//enough t obs?
				int numTObs = (int)(numCon+numNonCon);
				if (numTObs < minimumReadCoverage) continue;

				//save it
				ibAL.add(new MethylatedBaseObservationOneSample (testPos, numNonCon, numCon));

			}
			MethylatedBaseObservationOneSample[] ibs = null;
			if (ibAL.size() !=0){
				ibs = new MethylatedBaseObservationOneSample[ibAL.size()];
				ibAL.toArray(ibs);
			}
			return ibs;
		}


		public int getPosition() {
			return position;
		}

		public void setPosition(int position) {
			this.position = position;
		}

		public float getNonCon() {
			return nonCon;
		}

		public void setNonCon(float nonCon) {
			this.nonCon = nonCon;
		}

		public float getCon() {
			return con;
		}

		public void setCon(float con) {
			this.con = con;
		}


	}