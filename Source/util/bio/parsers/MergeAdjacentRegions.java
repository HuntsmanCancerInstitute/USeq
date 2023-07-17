package util.bio.parsers;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.Region;
import util.apps.FindOverlappingGenes;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Collapses transcripts into consensus gene based on exon overlap.*/
public class MergeAdjacentRegions {

	private File bedFile;
	private File resultsFile;
	private int maxGap = 5000;
	private Gzipper out = null;
	private int totalRegions = 0;
	private int finalRegions = 0;
	private boolean verbose = true;

	/**Merges regions within max gap and tracks number that were merged.  Assumes they do not overlap, thus run MergeRegions on the bed first.*/
	public  MergeAdjacentRegions(String[] args) {
		//process user arguments
		processArgs(args);
		
		//load models
		HashMap<String, Region[]> regions = Bed.parseRegions(bedFile, true);

		//for each chrom strand
		if (verbose) System.out.println("Processing...");
		try {
			out = new Gzipper(resultsFile);
			//for each chromosome
			for (String chrom: regions.keySet()){
				String chr = chrom.substring(0, chrom.length()-1);
				if (verbose) IO.pl("\t"+chr);
				Region[] sortedRegions = regions.get(chrom);
				totalRegions += sortedRegions.length;
				walkRegions(sortedRegions, chr);
			}
			out.close();
		} catch (IOException e) {
			resultsFile.deleteOnExit();
			e.printStackTrace();
		}
		
		if (verbose) IO.pl("\n"+totalRegions+ " regions merged to "+finalRegions);
	}

		
	private void walkRegions(Region[] r, String chr) throws IOException {
		//assuming these don't overlap
		int start = r[0].getStart();
		int stop = r[0].getStop();
		int numRegions = 1;
		for (int i=1; i< r.length; i++){
			int diff = r[i].getStart() - stop;
			if (diff > maxGap){
				//reset
				out.println(chr+"\t"+start+"\t"+stop+"\t.\t"+numRegions+"\t.");
				start = r[i].getStart();
				stop = r[i].getStop();
				numRegions = 1;
				finalRegions++;
			}
			else {
				stop = r[i].getStop();
				numRegions++;
			}
		}
		//complete last
		out.println(chr+"\t"+start+"\t"+stop+"\t.\t"+numRegions+"\t.");
		finalRegions++;
		
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergeAdjacentRegions (args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFile = new File (args[++i]); break;
					case 'r': resultsFile = new File (args[++i]); break;
					case 'm': maxGap = Integer.parseInt(args[++i]); break;
					case 'q': verbose = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your bed file?\n");
		if (resultsFile == null) Misc.printErrAndExit("\nError: please provide a xxx.bed.gz file to save the results.\n");
		
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Merge Adjacent Regions: April 2022                     **\n" +
				"**************************************************************************************\n" +
				"Merges regions within a max bp gap and tracks the number merged.  Regions must not\n"+
				"overlap. Best run the MergeRegions app if in doubt.\n\n" +

				"Options:\n"+
				"-b Path to a bed file of non overlapping regions, xxx.gz/.zip OK.\n"+
				"-r Path for saving the merged xxx.bed.gz file.\n"+
				"-m Max bp gap, defaults to 5000.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MergeAdjacentRegions -b myRegions.bed.zip \n" +
				"     -m 1000 -r mergedRegions.bed.gz \n\n" +

		"**************************************************************************************\n");

	}


}
