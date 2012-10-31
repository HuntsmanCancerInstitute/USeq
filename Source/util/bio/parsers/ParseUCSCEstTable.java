package util.bio.parsers;
import java.io.*;
import java.util.*;
import util.gen.*;

public class ParseUCSCEstTable {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		File f = new File(args[0]);
		File plus = new File(args[0]+".plus.bed");
		File minus = new File(args[0]+".minus.bed");
		try {
			BufferedReader in = new BufferedReader ( new FileReader(f));
			PrintWriter outPlus = new PrintWriter (new FileWriter(plus));
			PrintWriter outMinus = new PrintWriter (new FileWriter(minus));
			String line;
			in.readLine();
			while ((line = in.readLine()) != null){
				String[] tokens = line.split("\\t");
				String strand = tokens[0];
				String chrom = tokens[2];
				int[] starts = Num.stringArrayToInts(tokens[7], ",");
				int[] sizes = Num.stringArrayToInts(tokens[6], ",");
				PrintWriter out = outPlus;
				if (strand.equals("-")) out = outMinus;
				StringBuilder sb = new StringBuilder();
				for (int i=0; i< starts.length; i++){
					sb.append(chrom);
					sb.append("\t");
					sb.append(starts[i]);
					sb.append("\t");
					sb.append(starts[i]+sizes[i]);
					sb.append("\n");
				}
				out.print(sb.toString());
			}
			outPlus.close();
			outMinus.close();
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

}
