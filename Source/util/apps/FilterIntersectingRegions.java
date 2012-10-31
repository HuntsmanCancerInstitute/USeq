package util.apps;
import java.io.*;

import util.bio.annotation.Bed;
import util.gen.*;
import edu.utah.seq.useq.data.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.Region;


/**Takes two bed files. Collapses the first onto each chromosome, for each region in the second, calculates the 
 * fraction intersection and sorts into two files, intersecting and non intersecting.*/
public class FilterIntersectingRegions {

	//fields
	private File toFlatten;
	private HashMap<String,Region[]> regionsForFlattening;
	private File[] toSplit;
	private HashMap<String,RegionScoreText[]> regionsForSplitting;
	private double minimumFractionIntersection = Double.MIN_VALUE;
	private PrintWriter outIntersected;
	private PrintWriter outNonIntersected;
	private double numberIntersecting = 0;
	private double numberNonIntersecting = 0;

	public FilterIntersectingRegions(String[] args) {
		//process args
		processArgs(args);

		//load regions
		System.out.println("Loading regions...");
		regionsForFlattening = Region.parseStartStops(toFlatten, 0, 0, 0);
		
		//for each file to split
		for (File splitMe: toSplit){
			System.out.print("Intersecting "+splitMe.getName());
			
			regionsForSplitting = Bed.parseBedFile(splitMe, true);

			//for each chromosome to split
			Iterator<String> it = regionsForSplitting.keySet().iterator();
			File intersected = new File (splitMe.getParentFile(), Misc.removeExtension(splitMe.getName())+"_Int.bed");
			File nonIntersected = new File (splitMe.getParentFile(), Misc.removeExtension(splitMe.getName())+"_NonInt.bed");
			try {
				outIntersected = new PrintWriter( new FileWriter( intersected));
				outNonIntersected = new PrintWriter( new FileWriter( nonIntersected));
				while (it.hasNext()){
					//get text of chromosome and maxBase
					String chrom = (String)it.next();
					//get SS to split
					RegionScoreText[] ssToSplit = regionsForSplitting.get(chrom);
					//make mask
					boolean[] mask = makeMask(chrom);
					//anything to mask?
					if (mask == null) printNonIntersected(ssToSplit, chrom);
					else {
						RegionScoreText[][] split = intersect(ssToSplit, mask);
						printIntersected(split[0], chrom);
						printNonIntersected(split[1], chrom);
					}
					System.out.print(".");
				}
				outIntersected.close();
				outNonIntersected.close();
			} catch (IOException e){
				e.printStackTrace();
			}	

			//stats
			System.out.println("\nStats:");
			System.out.println((int)numberIntersecting+"\tIntersecting");
			System.out.println((int)numberNonIntersecting+"\tNon Intersecting");
			double fract = numberIntersecting/(numberIntersecting+numberNonIntersecting);
			System.out.println(Num.formatNumber(fract, 3) +"\tFraction intersecting");
			System.out.println();
		}
		System.out.println("\nDone!\n");
	}

	/**Returns intersecting and non intersecting*/
	public RegionScoreText[][] intersect(RegionScoreText[] ss, boolean[] bps){
		ArrayList<Region> intAL = new ArrayList<Region>();
		ArrayList<Region> nonIntAL = new ArrayList<Region>();

		//for each region count number to trues/ masked bases
		int maxBase = bps.length;
		for (int i=0; i< ss.length; i++){
			double numTrue = 0;
			int start = ss[i].getStart();
			int stop = ss[i].getStop();			
			//past the mask?
			if (start >= maxBase) nonIntAL.add(ss[i]);
			else {
				if (stop > bps.length) stop = bps.length;
				for (int j=start; j< stop; j++) if (bps[j]) numTrue++;
				//calc fraction
				double fractInt = numTrue/(double)ss[i].getLength();
				if (fractInt >= minimumFractionIntersection) intAL.add(ss[i]);
				else nonIntAL.add(ss[i]);
			}
		}
		RegionScoreText[] intersecting = new RegionScoreText[intAL.size()];
		intAL.toArray(intersecting);
		RegionScoreText[] nonIntersecting = new RegionScoreText[nonIntAL.size()];
		nonIntAL.toArray(nonIntersecting);
		return new RegionScoreText[][]{intersecting, nonIntersecting};
	}

	public void printNonIntersected(RegionScoreText[] ss, String chromosome) throws IOException{
		if (ss == null) return;
		numberNonIntersecting += ss.length;
		String c = chromosome +"\t";
		for (int i=0; i< ss.length; i++){
			outNonIntersected.println(c + ss[i].getStart()+"\t"+ss[i].getStop()+"\t"+ss[i].getText()+"\t"+ss[i].getScore()+"\t.");
		}
	}

	public void printIntersected(RegionScoreText[] ss, String chromosome) throws IOException{
		if (ss == null) return;
		numberIntersecting += ss.length;
		String c = chromosome +"\t";
		for (int i=0; i< ss.length; i++){
			//chr,start,stop,name,score,strand
			outIntersected.println(c + ss[i].getStart()+"\t"+ss[i].getStop()+"\t"+ss[i].getText()+"\t"+ss[i].getScore()+"\t.");
		}
	}

	public boolean[] makeMask (String chromosome){
		//any mask?
		if (regionsForFlattening.containsKey(chromosome) == false) return null;
		//find max base
		int maxBase = findMaxBase(chromosome);
		//make boolean array to hold whether it's flagged, initially they are all false
		boolean[] bps = new boolean[maxBase+10000];
		//for each RegionScoreText[] scan and throw booleans to true
		Region[] ss = regionsForFlattening.get(chromosome);
		for (int i=0; i<ss.length; i++){
			int start = ss[i].getStart();
			int end = ss[i].getStop();
			for (int j=start; j< end; j++) bps[j] = true;
		}
		return bps;
	}

	public int findMaxBase (String chromosome){
		Region[] ss = regionsForFlattening.get(chromosome);
		int max = ss[ss.length-1].getStop();
		for (int i=0; i< ss.length; i++){
			if (ss[i].getStop()> max) max = ss[i].getStop();
		}
		return max;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FilterIntersectingRegions(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': toFlatten = new File (args[++i]); break;
					case 's': toSplit = IO.extractFiles(new File(args[++i])); break;
					case 'i': minimumFractionIntersection = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (toFlatten == null || toFlatten.canRead() == false) Misc.printErrAndExit("\nCannot find your mask file! "+toFlatten+"\n");
		if (toSplit == null || toSplit[0].canRead() == false) Misc.printErrAndExit("\nCannot find your file(s) to split! "+toSplit+"\n");
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Filter Intersecting Regions: Feb 2011                     **\n" +
				"**************************************************************************************\n" +
				"Flattens the mask regions file and uses it to split the split file into intersecting\n" +
				"and non intersecting regions based on the minimum fraction intersection. Assumes\n" +
				"interbase coordinates.\n" +

				"\nOptions:\n"+
				"-m Full path file text for the masking bed file (tab delim: chr start stop ...).\n"+
				"-s Full path file text for the bed file to split into intersecting and non\n" +
				"        intersecting regions. Can also point to a directory of files to split.\n"+
				"-i Minimum fraction of each split region required to score as an intersection with\n"+
				"        the flattened mask, defaults to 1x10-1074\n"+

				"\nExample: java -Xmx4000M -jar pathTo/Apps/FilterIntersectingRegions -i 0.5\n" +
				"        -m /ArrayDesigns/repMskedDesign.bed -s /ArrayDesigns/novoMskedDesign.bed\n\n" +

		"************************************************************************************\n");
	}

}
