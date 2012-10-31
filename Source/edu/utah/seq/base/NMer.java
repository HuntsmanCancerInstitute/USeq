package edu.utah.seq.base;

import util.gen.StandardDeviation;

public class NMer {

		//fields
		private StandardDeviation[] sdA;
		private StandardDeviation[] sdC;
		private StandardDeviation[] sdG;
		private StandardDeviation[] sdT;
		private byte nMerSize;
		private float[] meansACGT = null;
		private float[] standardDeviationsACGT = null;

		public NMer(byte nMerSize){
			this.nMerSize = nMerSize;
			sdA = new StandardDeviation[nMerSize];
			sdC = new StandardDeviation[nMerSize];
			sdG = new StandardDeviation[nMerSize];
			sdT = new StandardDeviation[nMerSize];
			for (int i=0; i< nMerSize; i++) sdA[i] = new StandardDeviation();
			for (int i=0; i< nMerSize; i++) sdC[i] = new StandardDeviation();
			for (int i=0; i< nMerSize; i++) sdG[i] = new StandardDeviation();
			for (int i=0; i< nMerSize; i++) sdT[i] = new StandardDeviation();
		}
		
		public void count(int startIndex, int endIndex, short[] aIntensities, short[] cIntensities, short[] gIntensities, short[] tIntensities) throws Exception{
			int index = 0;
			for (int i=startIndex; i< endIndex; i++){
				sdA[index].count(aIntensities[i]);
				sdC[index].count(cIntensities[i]);
				sdG[index].count(gIntensities[i]);
				sdT[index].count(tIntensities[i]);
				index++;
			}
		}
		
		public void loadMeansStandarDeviationsAndNullObjects(){
			//build vector
			meansACGT = new float[4 * nMerSize];
			standardDeviationsACGT = new float[4 * nMerSize];
			//for each position in the nMer
			int counter = 0;
			//A
			for (int i=0; i< nMerSize; i++) {
				meansACGT[counter] = (float)sdA[i].getMean();
				standardDeviationsACGT[counter++] = (float) sdA[i].getStandardDeviation();
			}
			//C
			for (int i=0; i< nMerSize; i++) {
				meansACGT[counter] = (float)sdC[i].getMean();
				standardDeviationsACGT[counter++] = (float) sdC[i].getStandardDeviation();
			}
			//G
			for (int i=0; i< nMerSize; i++) {
				meansACGT[counter] = (float)sdG[i].getMean();
				standardDeviationsACGT[counter++] = (float) sdG[i].getStandardDeviation();
			}
			//T
			for (int i=0; i< nMerSize; i++) {
				meansACGT[counter] = (float)sdT[i].getMean();
				standardDeviationsACGT[counter++] = (float) sdT[i].getStandardDeviation();
			}
			sdA=null;
			sdC=null;
			sdG=null;
			sdT=null;
		}
		
		public double[] fetchACTGMeans(){
			//build vector
			double[] vals = new double[4 * nMerSize];
			//for each position in the nMer
			int counter = 0;
			//A
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdA[i].getMean();
			//C
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdC[i].getMean();
			//G
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdG[i].getMean();
			//T
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdT[i].getMean();
			return vals;
		}
		
		public double[] fetchACTGSig2Noise(){
			//build vector
			double[] vals = new double[4 * nMerSize];
			//for each position in the nMer
			int counter = 0;
			//A
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdA[i].getSignal2Noise();
			//C
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdC[i].getSignal2Noise();
			//G
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdG[i].getSignal2Noise();
			//T
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdT[i].getSignal2Noise();
			return vals;
		}
		
		public double[] fetchACTGStandardDeviations(){
			//build vector
			double[] vals = new double[4 * nMerSize];
			//for each position in the nMer
			int counter = 0;
			//A
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdA[i].getStandardDeviation();
			//C
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdC[i].getStandardDeviation();
			//G
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdG[i].getStandardDeviation();
			//T
			for (int i=0; i< nMerSize; i++) vals[counter++] = sdT[i].getStandardDeviation();
			return vals;
		}
		
		public String getMeanMatrix(){
			StringBuilder sb = new StringBuilder();
			sb.append("mnA");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdA[i].getMean());
			}
			sb.append("\nmnC");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdC[i].getMean());
			}
			sb.append("\nmnG");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdG[i].getMean());
			}
			sb.append("\nmnT");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdT[i].getMean());
			}
			sb.append("\n");
			return sb.toString();
		}
		
		public String getCVMatrix(){
			StringBuilder sb = new StringBuilder();
			sb.append("cvA");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdA[i].getCoefficientOfVariation());
			}
			sb.append("\ncvC");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdC[i].getCoefficientOfVariation());
			}
			sb.append("\ncvG");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdG[i].getCoefficientOfVariation());
			}
			sb.append("\ncvT");
			for (int i=0; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdT[i].getCoefficientOfVariation());
			}
			sb.append("\n");
			return sb.toString();
		}
		

		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(sdA[0].toString());
			sb.append("\t");
			sb.append(sdC[0].toString());
			sb.append("\t");
			sb.append(sdG[0].toString());
			sb.append("\t");
			sb.append(sdT[0].toString());
			//for each base, print ACGT mean, sd
			for (int i=1; i< nMerSize; i++){
				sb.append("\t");
				sb.append(sdA[i].toString());
				sb.append("\t");
				sb.append(sdC[i].toString());
				sb.append("\t");
				sb.append(sdG[i].toString());
				sb.append("\t");
				sb.append(sdT[i].toString());
			}
			return sb.toString();
		}

		public long getNumberObservations(){
			return sdA[0].getNumberObservations();
		}
		
		/**First checks to see if means have been fetched, if not then loads this NMer with these and nulls the StandardDeviation objects.*/
		public float[] getMeansACGT() {
			if (meansACGT == null) loadMeansStandarDeviationsAndNullObjects();
			return meansACGT;
		}

		/**First checks to see if standard deviations have been fetched, if not then loads this NMer with these and nulls the StandardDeviation objects.*/
		public float[] getStandardDeviationsACGT() {
			if (standardDeviationsACGT == null) loadMeansStandarDeviationsAndNullObjects();
			return standardDeviationsACGT;
		}

	
}
