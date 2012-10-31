package util.bio.parsers;
import java.io.*;
import java.util.*;


/**
 * Parses and extracts a unique set of gene names from a drosophila gene collection file.
 */
public class DrosGeneCollFileParser {

	
	public static void main(String[] args){
		try{
			BufferedReader in = new BufferedReader(new FileReader("/Users/nix/Desktop/DmelR4.0/dgc3.txt"));
			PrintWriter out = new PrintWriter(new FileWriter("/Users/nix/Desktop/dgcB.txt")); 
			HashSet uni = new HashSet();
			String line;
			int counter = 0;
			//skip first line
			System.out.println("Skipping line ->"+in.readLine());
			//
			while ((line = in.readLine())!=null){
				counter++;
				String[] tabbed = line.split("\\t");
				if (tabbed.length == 2){
					String[] tokens = tabbed[1].split(",|:");
					for (int i=0; i<tokens.length; i++){
						uni.add(tokens[i].trim());
					}
				}
			}
			Iterator it = uni.iterator();
			while (it.hasNext()){
				out.println(it.next());
			}
			in.close();
			out.close();
			System.out.println(uni.size()+" Unique gene names");
			System.out.println(counter + " Number of clones");
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	

	
	
}
