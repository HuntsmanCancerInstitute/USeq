package util.apps;
import java.io.*;

import trans.anno.GenomicRegion;
import trans.anno.GenomicRegionComparator;
import util.bio.seq.OligoTiler;
import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**Takes several tab delimited region files (chrom, start, stop, ...), creates the union and writes to disk.
 * Assumes the first base is zero, and the stop/stop of a region is excluded (interbase coordinates) .
 * 
 * Enter a text of a directory containing the files.*/
public class MergeRegions {
	
	//fields
	private File directory;
	private File[] regionFiles;
	private File mergedFile;
	private long numberRegions = 0;
	private long numberMergedRegions = 0;
	private int minimumCount = 1;
	private int pad = 0;
	
	public MergeRegions(String[] args) {
		//start clock
		long startTime = System.currentTimeMillis();
		
		//process args
		processArgs(args);
		doWork();
		
		IO.pl("\nMerged "+numberRegions+" -> "+numberMergedRegions);
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public MergeRegions (File[] regionFiles, File mergedFile){
		this.regionFiles = regionFiles;
		this.mergedFile = mergedFile;
		doWork();
	}
	
	public void doWork(){
		//load and sort all of the regions
		GenomicRegion[][] regions = new GenomicRegion[regionFiles.length][];
		GenomicRegionComparator comp = new GenomicRegionComparator();
		for (int i=0; i<regionFiles.length; i++) {
			regions[i] = GenomicRegion.parseRegions(regionFiles[i]);
			if (pad !=0) for (GenomicRegion g: regions[i]) g.pad(pad);
			numberRegions+= regions[i].length;
			Arrays.sort(regions[i], comp);
		}
		
		//find all chromosomes and the maxBase
		TreeMap<String, Integer> map = new TreeMap<String, Integer>();
		for (int i=0; i<regions.length; i++){
			for (int j=0; j< regions[i].length; j++){
				//does chromosome exist
				Integer max = map.get(regions[i][j].getChromosome());
				if (max == null){
					map.put(regions[i][j].getChromosome(), new Integer(regions[i][j].getEnd()));
				}
				//exists, reset maxBase?
				else {
					int num = max.intValue();
					if (regions[i][j].getEnd() > num) map.put(regions[i][j].getChromosome(), new Integer(regions[i][j].getEnd()));
				}
			}
		}
		
		//for each chromosome make merge
		Iterator<String> it = map.keySet().iterator();
		
		try {
			PrintWriter out = new PrintWriter( new FileWriter( mergedFile));
			out.println("# MinOverlap\t"+ minimumCount);
			out.println("# BpPadding\t"+ pad);
			for (File f: regionFiles) out.println("# "+f.getCanonicalPath());
			while (it.hasNext()){
				//get text of chromosome and maxBase
				String chrom = it.next();
				IO.p(chrom+" ");
				int maxBase =map.get(chrom);
				//make boolean array to hold whether it's flagged, initially they are all false
				//boolean[] bps = new boolean[maxBase+1];
				int[] bpsCounts = new int[maxBase+1];
				//for each GenomicRegion[] scan and throw booleans to true
				for (int i=0; i<regions.length; i++){
					boolean found = false;
					for (int j=0; j< regions[i].length; j++){
						//correct chromosome?
						if (regions[i][j].getChromosome().equals(chrom)){
							int stop = regions[i][j].getEnd();
							for (int k=regions[i][j].getStart(); k< stop; k++){
								//bps[k] = true;
								bpsCounts[k]++;
							}
							found = true;
						}
						else if (found) break;
					}
				}
				//print merged chrom
				print(chrom, thresholdCountArray(bpsCounts), out);
			}
			IO.pl();
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
	public boolean[] thresholdCountArray(int[] counts) {
		boolean[] bps = new boolean[counts.length];
		for (int i=0; i< counts.length; i++) if (counts[i]>=minimumCount) bps[i]=true;
		return bps;
	}
	
	public void print(String chromosome, boolean[] bps, PrintWriter out){
		boolean in = false;
		int i=0;
		for (; i<bps.length; i++){
			//start new line?
			if (bps[i] && in==false){
				out.print(chromosome+"\t"+i+"\t");
				in = true;
			}
			//close old?
			else if (in && bps[i] == false){
				out.println(i);
				in = false;
				numberMergedRegions++;
			}
		}
		//close old?
		if (in){
			if (bps[i-1]) out.println(i);
			else out.println(i-1);
			numberMergedRegions++;
		}
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergeRegions(args);
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
					case 'r': directory = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					case 'o': mergedFile = new File (args[++i]); break;
					case 'm': minimumCount = Integer.parseInt(args[++i]); break;
					case 'p': pad = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (directory == null) Misc.printErrAndExit("\nPlease enter a path to a region file or directory containing such.\n");
		regionFiles = IO.extractFiles(directory);
		
		//any mergedFile?
		if (mergedFile == null){
			if (directory.isDirectory() == false){
				mergedFile = new File(directory+".merged");
			}
			else mergedFile = new File(directory,"mergedUnion.txt");
		}
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Merge Regions: April 2020                            **\n" +
				"**************************************************************************************\n" +
				"Flattens tab delimited bed files (chr start stop ...). Assumes interbase coordinates.\n" +
				"Set the -m threshold to restrict the output to bps with that minimum of overlapping\n"+
				"regions.\n"+

				"\nOptions:\n"+
				"-r Path to a region file or directory containing such, xxx.gz/xxx.zip OK.\n"+
				"-m Minimum overlapping bp to save, defaults to 1\n"+
				"-p Pad input regions +/- bps, defaults to 0\n"+
				"-o Path to the merged output file, optional\n"+

				"\nExample: java -Xmx4G -jar pathTo/Apps/MergeRegions -d PassingBedFiles/-m 3 -p 150\n\n" +

		"************************************************************************************\n");
	}

	public void setNumberMergedRegions(int numberMergedRegions) {
		this.numberMergedRegions = numberMergedRegions;
	}
	
}
