
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import edu.utah.seq.data.*;
import util.gen.*;
import java.util.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class Eland2Fastq {

	public static void main(String[] args) {
		if (args.length ==0) Misc.printErrAndExit("\nEnter a directory containing xxx.sorted.txt.gz files to convert to xxx.sequence.txt.gz");
		File[] files = IO.extractFiles(new File(args[0]));//, ".sorted.txt.gz");
		for (File f: files){
			try {
				BufferedReader in = IO.fetchBufferedReader(f);
				System.out.println(f.getName());
				File seqFile = new File (f.getCanonicalPath().replace("sorted.txt.gz", "sequence.txt"));
				if (seqFile.exists()) Misc.printErrAndExit("\nAlready exits?! "+seqFile);
				PrintWriter out = new PrintWriter( new FileWriter(seqFile));
				String line;
				Pattern tab = Pattern.compile("\\t");
				while ((line = in.readLine()) != null){
					String[] tokens = tab.split(line);
					//make header
					StringBuilder sb = new StringBuilder();
					sb.append(tokens[0]); sb.append (":");
					sb.append(tokens[2]); sb.append (":");
					sb.append(tokens[3]); sb.append (":");
					sb.append(tokens[4]); sb.append (":");
					sb.append(tokens[5]); sb.append ("#0/1");
					//write seq header
					out.print("@");
					out.println(sb.toString());
					//write sequence
					out.println(tokens[8]);
					//write quality header
					out.print("+");
					out.println(sb.toString());
					//write qualities
					out.println(tokens[9]);
				}
				out.close();
				in.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}		

	
	

}
