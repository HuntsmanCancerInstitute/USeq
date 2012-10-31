package edu.utah.seq.parsers;
import java.io.*;
import util.gen.*;
import java.util.regex.*;

public class JohnsonToTag {

	public static void main(String[] args) {
		int lengthSeq = 24;
		
		File[] files = IO.extractFiles(new File(args[0]), ".txt");
		for (int i=0; i< files.length; i++){
			try{
				BufferedReader in = new BufferedReader (new FileReader (files[i]));
				PrintWriter out = new PrintWriter (new FileWriter(Misc.removeExtension(files[i].getCanonicalPath())+".tag"));
				String line;
				String[] tokens;
				String strand;
				String chrom;
				int start;
				Pattern pat = Pattern.compile("\\s+");
				Pattern chr = Pattern.compile("(.+):(.+)");
				Matcher mat;
				while ((line = in.readLine()) != null){
					tokens = pat.split(line);
					if (tokens.length != 7) {
						System.out.println("Odd line! "+line);
						continue;
					}
					//strand
					if (tokens[4].equals("R")) strand = "-";
					else strand = "+";
					//chrom start
					mat = chr.matcher(tokens[3]);
					if (mat.matches() == false) {
						System.out.println("Couldn't parse chr:start from "+line);
						continue;
					}
					chrom = mat.group(1);
					start = Integer.parseInt(mat.group(2))-1;
					//print line
					out.println(chrom+"\t"+start+"\t"+(start+lengthSeq)+"\t"+strand);
					
				}
				out.close();
				in.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}

	}

}
