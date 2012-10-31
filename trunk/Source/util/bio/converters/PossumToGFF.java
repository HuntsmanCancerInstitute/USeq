package util.bio.converters;
import java.io.*;
import java.util.*;
/**
Converts a copied Possum output to GFF annotation, prints to screen, sorted by score.
 */
public class PossumToGFF {
	public static void main(String[] args) {
		new PossumToGFF();
	}
	public PossumToGFF(){
	//read in line by line and print to this format
	//chromosome, source, feature, start, stop, score, strand,frame,attributes/comments
	String line;
	GffLine[] lines = null;
	try {
		ArrayList al = new ArrayList();
		BufferedReader in = new BufferedReader(new FileReader("/Users/nix/Desktop/delme.txt"));
		while ((line = in.readLine()) !=null) {
			line = line.trim();                     //kill whitespace and test if exists
			if (line.length() == 0) continue;       //skip blank lines
			String[] items = line.split(" ");
			
			//o text, 1 start, 3 stop, 4 ori, 5 seq, 6 score(nat log)
			
			al.add( new GffLine(Double.parseDouble(items[6]), items[0]+"\tPossum\t"+items[0]+"\t"+items[1]+"\t"+
				items[3]+"\t"+items[6]+"\t"+items[4]+"\t.\tseq="+items[5]+";"));
		}
		in.close();
		int len = al.size();
		lines = new GffLine[len];
		al.toArray(lines);
		Arrays.sort(lines);
		
		for (int i=0; i<len; i++){
			System.out.println(lines[i].line);
		}
	}
	catch (IOException e) {
		e.printStackTrace();
	}
	
	}
	private class GffLine implements Comparable {
		double score;
		String line;
		public GffLine (double score, String line){
			this.score = score;
			this.line = line;
		}
		public int compareTo(Object obj){
			GffLine other = (GffLine)obj;
			if (other.score> score) return 1;
			if (other.score< score) return -1;
			return 0;
		}
	}
}
