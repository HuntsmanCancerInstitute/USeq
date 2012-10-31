package trans.roc;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.Region;

import util.bio.annotation.*;
import util.bio.seq.Seq;
import util.gen.*;

/**
 * Parses bed files for those that intersect a list of regions
 */
public class ParseBedFilesForParticularRegions {
	private File[] bedFiles;
	private File regionFile;
	private File maskFile;
	private HashMap<String,ArrayListStartStop[]> regions;
	private HashMap<String,Region[]> regionsToSubtract;
	private boolean printSeqReadStats = false;


	public ParseBedFilesForParticularRegions(String[] args){
		processArgs(args);

		//load regions file split by chromosome
		regions = ArrayListStartStop.parseArrayListStartStops(regionFile, 0, 0, false);

		//for each bed file, add lines to regions
		System.out.println("Loading...");
		for (int x=0; x<bedFiles.length; x++){
			Bed[] bedLines = Bed.parseFile(bedFiles[x], 0, 0);
			System.out.println("\t"+bedFiles[x].getName());
			for (int i=0; i< bedLines.length; i++){
				ArrayListStartStop[] ss = regions.get(bedLines[i].getChromosome());
				if (ss == null) continue;
				int start = bedLines[i].getStart();
				int stop = bedLines[i].getStop();
				boolean intersect = false;
				for (int j=0; j< ss.length; j++){
					if (ss[j].intersects(start, stop)) {
						ss[j].getArrayList().add(bedLines[i]);
						intersect = true;
					}
					else if (intersect) break;
				}
			}
		}

		//filter?
		if (maskFile != null) {
			System.out.println("Filtering...");
			regionsToSubtract = Region.parseStartStops(maskFile, 0, 0, 0);
			filterBedLines();
		}

		//print out results
		System.out.println("\nPrinting...");
		if (printSeqReadStats) printResultsSeqStats();
		else printResults();
		System.out.println("\nDone!\n");
	}

	public void filterBedLines(){
		//for each chromosome
		Iterator<String> it = regions.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			ArrayListStartStop[] ss = regions.get(chrom);
			Region[] ssMask = regionsToSubtract.get(chrom);
			if (ssMask == null) continue;
			//for each region see if any of it's bed lines intersect mask
			for (int i=0; i< ss.length; i++){
				ArrayList<Bed> bedLines = ss[i].getArrayList();
				int num = bedLines.size();
				if (num ==0) continue;
				else {
					//for each bed line
					for (int j=0; j< bedLines.size(); j++){
						Bed b = bedLines.get(j);
						int start = b.getStart();
						int stop = b.getStop();
						//for each masked region
						for (int x=0; x< ssMask.length; x++){
							if (ssMask[x].intersects(start, stop)){
								bedLines.remove(j);
								j--;
								break;
							}
						}
					}
				}
			}
		}
	}

	public void printResultsSeqStats(){
		try {
			File results = new File(Misc.removeExtension(regionFile.getCanonicalPath())+"Extracted.bed");
			System.out.println("\t"+results);
			PrintWriter out = new PrintWriter(new FileWriter(results));
			out.println("#\tChromosome\tStart\tStop\tMedianLength\tMeanLength\tFract5'G\tFract5'A\tFract5'T\tFract5'C\tFract5'N");
			//for each chromosome
			Iterator<String> it = regions.keySet().iterator();
			while (it.hasNext()){
				String chrom = it.next();
				ArrayListStartStop[] ss = regions.get(chrom);
				//for each region
				for (int i=0; i< ss.length; i++){
					ArrayList bedLines = ss[i].getArrayList();
					int num = bedLines.size();				
					if (num == 0) out.println("#\t"+chrom+"\t"+ss[i]);
					else {
						//calculate median length
						StringBuilder sb = new StringBuilder();
						int[] lengths = new int[num];
						StringBuilder fivePrimeBases = new StringBuilder();
						for (int j=0; j< num; j++){
							Bed bed = (Bed)bedLines.get(j);
							sb.append(bed); sb.append("\n");
							lengths[j] = bed.getLength();
							fivePrimeBases.append(bed.getName().substring(0,1));
						}
						Arrays.sort(lengths);
						double medianLength = Num.median(lengths);
						double averageLength = Num.mean(lengths);
						double[] fractionGATCN = Seq.calculateFractionBases(fivePrimeBases.toString());
						String strFrac = Num.doubleArrayToString(fractionGATCN, 3, "\t");
						out.println("#\t"+chrom+"\t"+ss[i]+"\t"+Num.formatNumber(medianLength, 2)+"\t"+Num.formatNumber(averageLength, 2)+"\t"+strFrac);
						out.println(sb);
					}
				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void printResults(){
		try {
			File results = new File(Misc.removeExtension(regionFile.getCanonicalPath())+"Extracted.bed");
			System.out.println("\t"+results);
			PrintWriter out = new PrintWriter(new FileWriter(results));
			out.println("#\tChromosome\tStart\tStop");
			//for each chromosome
			Iterator<String> it = regions.keySet().iterator();
			while (it.hasNext()){
				String chrom = it.next();
				ArrayListStartStop[] ss = regions.get(chrom);
				//for each region
				for (int i=0; i< ss.length; i++){
					ArrayList bedLines = ss[i].getArrayList();
					int num = bedLines.size();				
					out.println("#\t"+chrom+"\t"+ss[i]);
					for (int j=0; j< num; j++){
						Bed bed = (Bed)bedLines.get(j);
						out.println(bed);
					}
				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ParseBedFilesForParticularRegions(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File bedDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedDir = new File(args[++i]); break;
					case 'r': regionFile = new File (args[++i]); break;
					case 'm': maskFile = new File (args[++i]); break;
					case 's': printSeqReadStats = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull files
		if (bedDir == null) Misc.printExit("\nCannot find any xxx.bed files to parse?!\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(bedDir,".bed");
		tot[1] = IO.extractFiles(bedDir,".bed.zip");
		tot[2] = IO.extractFiles(bedDir,".bed.gz");
		bedFiles = IO.collapseFileArray(tot);

		if (bedFiles == null || bedFiles.length == 0) Misc.printExit("\nCannot find any xxx.bed files to parse?!\n");
		if (regionFile == null) Misc.printExit("\nCannot find your region file?!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       ParseBedFilesForParticularRegions: May 2009                **\n" +
				"**************************************************************************************\n" +
				"PBFFPRs takes a list of regions uses it to select for intersecting regions from other\n" +
				"lists.  Useful for extracting tag file data for regions of interest.\n\n" +

				"Options:\n"+
				"-r Regions file (tab delimited: chr, start, stop,...) to use in fetching tags.\n"+
				"-m (Optional) Mask file (tab delimited: chr, start, stop,...) to use in excluding tags.\n"+
				"-b Full path file text for the tags xxx.bed(.zip/.gz) file or directory containing\n" +
				"      such to parse. Full bed format required (tab del: chr start stop text\n" +
				"      strand(+/-/.) score)\n"+
				"-s Print sequence stats, assumes the tags xxx.bed use a sequence for their names.\n"+

				"\nExample: java -Xmx1500M -jar pathTo/USeq/Apps/ParseBedFilesForParticularRegions -r\n" +
				"      /BedFiles/ToTile/regions.bed -b /BedFiles/Annotation/exons.bed \n\n" +

		"**************************************************************************************\n");

	}

}
