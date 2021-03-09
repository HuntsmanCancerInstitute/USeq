package edu.utah.seq.vcf.xml.caris;


public class SimpleBed implements Comparable<SimpleBed> {

		//fields
		private String chr;
		private int start;
		private int stop;
		private String name;
		private String score;
		private String strand;

		//constructor
		public SimpleBed (String chr, int start, int stop, String name, String score, String strand){
			this.chr = chr;
			this.start = start;
			this.stop = stop;
			this.name = name;
			this.score = score;
			this.strand = strand;
			if (name == null) this.name = ".";
			if (score == null) this.score = ".";
			if (strand == null) this.strand = ".";
		}
		
		/**Sorts by chromsome, start position, length (smallest to largest).*/
		public int compareTo(SimpleBed otherCoor){
			//sort by chromosome
			int compare = otherCoor.chr.compareTo(chr);
			if (compare !=0) return compare * -1;;
			//sort by start position
			if (start<otherCoor.start) return -1;
			if (start>otherCoor.start) return 1;
			// if same start, sort by length, smaller to larger
			int len = stop-start;
			int otherLen = otherCoor.stop-otherCoor.start;
			if (len<otherLen) return -1;
			if (len>otherLen) return 1;
			return 0;
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder(chr);
			sb.append("\t"); sb.append(new Integer(start).toString());
			sb.append("\t"); sb.append(new Integer(stop).toString());
			sb.append("\t"); sb.append(name);
			sb.append("\t"); sb.append(score);
			sb.append("\t"); sb.append(strand);
			return sb.toString();
		}
	}
