package trans.misc;

/**Container for a ratio and associated %GC content for a particular oligo.*/
public class GCRatio implements Comparable{
		private double ratio;
		private int gc;
		
		public  GCRatio (double ratio, int gc){
			this.ratio = ratio;
			this.gc = gc;
		}
		public int compareTo (Object obj){
			GCRatio other = (GCRatio) obj;
			if (other.gc > gc ) return -1;
			if (other.gc < gc ) return 1;
			return 0;
		}
		public int getGc() {
			return gc;
		}
		public double getRatio() {
			return ratio;
		}
	}



