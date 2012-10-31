package util.bio.seq;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.Info;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import edu.utah.seq.parsers.BarParser;
import util.gen.*;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.*;

/**
 * Converts fasta files into bar graph files noting each CpG context.
 */
public class ConvertFasta2GCBarGraph {
	//fields
	private File directory;
	private String genomeVersion;

	public ConvertFasta2GCBarGraph(String[] args){

		processArgs(args);
		System.out.println("Converting...");

		boolean txt = false;
		File[] files = IO.extractFiles(directory);

		for (int i=0; i< files.length; i++){
			System.out.println("\t"+files[i]);
			
			//fetch sequence
			MultiFastaParser fastaParser = new MultiFastaParser();
			fastaParser.parseIt(files[i]);
			if (fastaParser.isFastaFound()== false) Misc.printExit("\nError: fasta not found! Aborting.\n");
			char[] seq = fastaParser.getSeqs()[0].toLowerCase().toCharArray();
			String chrom = fastaParser.getNames()[0];
			String strand = "+";
			
			//scan for CpG contexts
			ArrayList<Point> al = new ArrayList<Point>();
			for (int x=0; x< seq.length; x++){
				if (seq[x] == 'c') {
					int next = x+1;
					if (next == seq.length) break;
					if (seq[next] == 'g') {
						al.add(new Point(x,1));
						x++;
					}
				}
			}
			
			//make pd
			PointData pd = Point.extractPositionScores(al);
			HashMap <String,String> notes = new HashMap <String,String> ();
			notes.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			notes.put(BarParser.STRAND_TAG, strand);
			//make an Info object {
			Info info = new Info(chrom+strand, genomeVersion, chrom, strand, 0, notes);
			info.setNumberObservations(al.size());
			pd.setInfo(info);
			//write to file
			pd.writePointData(directory);
		}

	}

	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': directory = new File(args[i+1]); i++; break;
					case 'v': genomeVersion = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		//check to see if they entered required params
		if (directory==null){
			System.out.println("\nCannot find your directory!\n");
			System.exit(0);
		}
		if (genomeVersion==null){
			System.out.println("\nPlease enter a genome version!\n");
			System.exit(0);
		}
		
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Convert Fasta 2 GC Bar Graphs: April 2011                    **\n" +
				"**************************************************************************************\n" +
				"Converts fasta files into graph files containing a 1 over each C in a CpG context.\n\n"+

				"Required Parameters:\n"+
				"-f Full path name for the directory containing xxx.fasta(.gz/.zip OK).\n" +
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +

				"\n" +
				"Example: java -Xmx4G -jar pathTo/Apps/ConvertFasta2GCBarGraph -f /affy/Fastas/\n" +
				"      -v H_sapiens_Feb_2009\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new ConvertFasta2GCBarGraph(args);
	}

}

