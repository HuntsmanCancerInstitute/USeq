package util.apps;
import java.io.*;

import trans.anno.GenomicRegion;
import trans.anno.RegionComparator;
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
	
	public MergeRegions(String[] args) {
		//process args
		processArgs(args);

		//load and sort regions
		GenomicRegion[][] regions = new GenomicRegion[regionFiles.length][];
		RegionComparator comp = new RegionComparator();
		for (int i=0; i<regionFiles.length; i++) {
			regions[i] = GenomicRegion.parseRegions(regionFiles[i]);
			Arrays.sort(regions[i], comp);
		}
		
		//find all chromosomes and the maxBase
		TreeMap map = new TreeMap();
		for (int i=0; i<regions.length; i++){
			for (int j=0; j< regions[i].length; j++){
				//does chromosome exist
				Object obj = map.get(regions[i][j].getChromosome());
				if (obj == null){
					map.put(regions[i][j].getChromosome(), new Integer(regions[i][j].getEnd()));
				}
				//exists, reset maxBase?
				else {
					int num = ((Integer)obj).intValue();
					if (regions[i][j].getEnd() > num) map.put(regions[i][j].getChromosome(), new Integer(regions[i][j].getEnd()));
				}
			}
		}
		
		//for each chromosome make merge
		Iterator it = map.keySet().iterator();
		File mergedFile;
		if (directory.isDirectory() == false){
			mergedFile = new File(directory+".merged");
		}
		else mergedFile = new File(directory,"mergedUnion.txt");
		try {
			PrintWriter out = new PrintWriter( new FileWriter( mergedFile));
			while (it.hasNext()){
				//get text of chromosome and maxBase
				String chrom = (String)it.next();
				int maxBase =((Integer)map.get(chrom)).intValue();
				//make boolean array to hold whether it's flagged, initially they are all false
				boolean[] bps = new boolean[maxBase+1];
				//for each GenomicRegion[] scan and throw booleans to true
				for (int i=0; i<regions.length; i++){
					boolean found = false;
					for (int j=0; j< regions[i].length; j++){
						//correct chromosome?
						if (regions[i][j].getChromosome().equals(chrom)){
							int stop = regions[i][j].getEnd();
							for (int k=regions[i][j].getStart(); k< stop; k++){
								bps[k] = true;
							}
							found = true;
						}
						else if (found) break;
					}
				}
				//print merged chrom
				print(chrom, bps, out);
			}
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}
		
		
	}
	public static void print(String chromosome, boolean[] bps, PrintWriter out){
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
			}
		}
		//close old?
		if (in){
			if (bps[i-1]) out.println(i);
			else out.println(i-1);
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
		File fastaDirectory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': directory = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (directory == null || directory.isDirectory() == false) {
			Misc.printErrAndExit("\nPlease enter a directory containing bed files to merge.\n");
		}
		regionFiles = IO.extractFiles(directory);
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Merge Regions: May 2009                              **\n" +
				"**************************************************************************************\n" +
				"Flattens tab delimited bed files (chr start stop ...). Assumes interbase coordinates.\n" +

				"\nOptions:\n"+
				"-d Directory containing bed files.\n"+

				"\nExample: java -Xmx4000M -jar pathTo/Apps/MergeRegions -d /Anno/TilingDesign/\n\n" +

		"************************************************************************************\n");
	}
	
}
