package util.bio.parsers;
import java.io.*;
import java.util.regex.*;
import util.gen.*;

public class RipInparanoidResults {

	public static void main(String[] args){
		
		//pull files
		File[] f = IO.extractFiles(new File(args[0]), ".html");
		
		//make patterns
		Pattern ens = Pattern.compile("geneid <b>(ENSG\\d+)</b>");
		Pattern zeb = Pattern.compile(">(ZDB-GENE-\\d+-\\d+)");
		for (int i=0; i<f.length; i++){
			String[] lines = IO.loadFileIntoStringArray(f[i]);
			for (int j=0; j< lines.length; j++){
				Matcher mat = ens.matcher(lines[j]);
				if (mat.find()){
					System.out.print(mat.group(1));
					j++;
					mat = zeb.matcher(lines[j]);
					if (mat.find()) System.out.print("\t"+ mat.group(1));
					else System.out.print("\tNo paralog");
					System.out.println();
				}
			}
			
		}
		
		
	}
	
	
}
