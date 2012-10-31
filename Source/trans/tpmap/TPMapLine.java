package trans.tpmap;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * Minimal tpmap line, chrom, start, nothing else.
 */
public class TPMapLine implements Comparable{
	private String line;
	private String chromosome;
	private String sequence;
	private String orientation = "t";
	private int start;
	private int pmX;
	private int pmY;
	
	/**For parsing from an Agilent data file.*/
	public TPMapLine(String[] items){
		//parse chromosome
		int index = items[12].indexOf(":");
		if (index == -1) {
			start = -1;
			return;
		}
		chromosome = items[12].substring(0,index);
		//sequence
		sequence = items[7];
		//parse start
		int stop = items[12].indexOf("-");
		String number = items[12].substring(index+1, stop);
		start = Integer.parseInt(number);
		
		//parse row and column coordinates, subtract 1 to put in zero based
		pmX = Integer.parseInt(items[2]) - 1;
		pmY = Integer.parseInt(items[3]) - 1;
	}
	
	/**For parsing from an Agilent promoter array data file. Note there isn't a sequence for the oligo!*/
	public TPMapLine(String chrNum, String row, String column, String sequence){
		//parse chromosome
		int index = chrNum.indexOf(":");
		if (index == -1) {
			start = -1;
			return;
		}
		chromosome = chrNum.substring(0,index);
		this.sequence = sequence;
		
		//parse start
		int stop = chrNum.indexOf("-");
		String number = chrNum.substring(index+1, stop);
		start = Integer.parseInt(number);
		
		//parse row and column coordinates, subtract 1 to put in zero based
		pmX = Integer.parseInt(row) - 1;
		pmY = Integer.parseInt(column) - 1;
	}
	
	/**Coordinate are zero based!*/
	public TPMapLine(String chr, int pos, int row, int column, String sequence){
		chromosome = chr;
		start = pos;
		pmX = row;
		pmY = column;
		this.sequence = sequence;
	}
	
	/**Coordinate are zero based!*/
	public TPMapLine(String chr, int pos, int row, int column, String sequence, boolean hybridizesForwardStrand){
		chromosome = chr;
		start = pos;
		pmX = row;
		pmY = column;
		this.sequence = sequence;
		if (hybridizesForwardStrand) orientation = "t";
		else orientation = "f";
	}

	
	/**For parsing from an Nimblegen NDF data file.*/
	public TPMapLine(String[] items, boolean trueFalse){
		//sequence
		sequence = items[5];
		//set orientation
		if (trueFalse) orientation = "t";
		else orientation = "f";
		//parse chromosome
		String chromNumber = items[4].substring(3);
		if (chromNumber.startsWith("0")) chromNumber = chromNumber.substring(1);
		chromosome = "chr"+chromNumber;
		//parse start, subtract to put in zero based coordinates
		start = Integer.parseInt(items[13])-1;
		//start = Integer.parseInt(items[13]);
		//parse row and column coordinates, subtract 1 to put in zero based
		pmX = Integer.parseInt(items[10]) - 1;
		pmY = Integer.parseInt(items[9]) - 1;

	}
	
	public TPMapLine(String line){
		this.line = line;
		String[] tokens = line.split("\\s+");
		sequence = tokens[0];
		orientation = tokens[1];
		chromosome = tokens[2];
		start = Integer.parseInt(tokens[3]);
		pmX = Integer.parseInt(tokens[4]);
		pmY = Integer.parseInt(tokens[5]);
	}
	
	public int compareTo(Object obj){
		TPMapLine other = (TPMapLine)obj;
		//sort by chromosome
		int comp = other.chromosome.compareTo(chromosome);
		if (comp != 0) return comp*-1;
		//sort by base position
		if (other.start>start) return -1;
		if (other.start<start) return 1;
		//sort by pmX
		if (other.pmX>pmX) return -1;
		if (other.pmX<pmX) return 1;
		//sort by pmY
		if (other.pmY>pmY) return -1;
		if (other.pmY<pmY) return 1;
		return 0;
	}
	
	public static boolean saveTPMap(TPMapLine[] sortedTPMap, File tpmapFile){
		try {
			PrintWriter outRes = new PrintWriter(new FileWriter(tpmapFile));
			for (int i=0; i<sortedTPMap.length; i++){
				outRes.println(sortedTPMap[i].getLine());
			}
			outRes.close();
			return true;
		} catch (Exception e){
			return false;
		}
	}
	public String getLine() {
		if (line == null) line = sequence + "\t"+orientation+"\t" + chromosome +"\t"+ start +"\t"+ pmX +"\t"+ pmY;
		return line;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getPmX() {
		return pmX;
	}

	public int getPmY() {
		return pmY;
	}

	public int getStart() {
		return start;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
}
