package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.its.Interval;
import edu.utah.seq.its.IntervalTree;
import edu.utah.seq.useq.data.RegionScoreText;
import edu.utah.seq.useq.data.RegionScoreTextData;
import util.bio.annotation.Bed;
import util.bio.annotation.Coordinate;
import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;

/**Splits variants into those that intersect and those that do not.
 * @author Nix
 * */
public class VCFRegionFilter {

	//user fields
	private File[] vcfFiles;
	private File bedFile;
	private int bpPad = 0;
	private File saveDirectory = null;
	
	//internal fields
	private HashMap<String,IntervalTree<Coordinate>> chrRegionIntervalTrees = null;
	private int numIntRecords = 0;
	private int numNonIntRecords = 0;
	
	
	//constructor
	public VCFRegionFilter(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//build interval tree for bed region file
		createIntervalTrees();

		System.out.println("Name\tNumInt\tNumNonInt");
		for (int i=0; i< vcfFiles.length; i++){
			filterVcf(vcfFiles[i]);
			System.out.println(vcfFiles[i].getName()+"\t"+numIntRecords+"\t"+numNonIntRecords);
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void filterVcf(File vcf) {
		try {
			numIntRecords = 0;
			numNonIntRecords = 0;
			String name = Misc.removeExtension(vcf.getName());
			Gzipper intVcf = new Gzipper( new File (saveDirectory, name + "_int.vcf.gz"));
			Gzipper nonIntVcf = new Gzipper( new File (saveDirectory, name + "_nonInt.vcf.gz"));
			BufferedReader in = IO.fetchBufferedReader(vcf);
			String line;
			String currChrom = "";
			IntervalTree<Coordinate> regions = null;
			//for each line in the file
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#")) {
					intVcf.println(line);
					nonIntVcf.println(line);
				}
				//data line
				else {
					//#CHROM POS ID REF ALT QUAL FILTER INFO ......
					//   0    1   2  3   4   5     6   
					String[] tokens = Misc.TAB.split(line);
					//diff chrom?
					if (tokens[0].equals(currChrom) == false) {
						//attempt to fetch regions for this chrom, might be null
						regions = chrRegionIntervalTrees.get(tokens[0]);
						currChrom = tokens[0];
					}
					//any regions to intersect? Does it intersect?
					if (regions != null) {
						ArrayList<Coordinate> hits = fetchRecords(tokens, regions);
						//no hits so just print out unmodified
						if (hits.size() == 0) {
							nonIntVcf.println(Misc.stringArrayToString(tokens, "\t"));
							numNonIntRecords++;
						}
						else {
							intVcf.println(Misc.stringArrayToString(tokens, "\t"));
							numIntRecords++;
						}
					}
					//nope so just print out as non int
					else {
						nonIntVcf.println(Misc.stringArrayToString(tokens, "\t"));
						numNonIntRecords++;
					}
				}
			}
			in.close();
			intVcf.close();
			nonIntVcf.close();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: filtering "+vcf);
		} 
	}
	
	private ArrayList<Coordinate> fetchRecords(String[] vcfTokens, IntervalTree<Coordinate> regions) {
		//calc start stop to fetch, interbase coordinates
		int size = vcfTokens[4].length();
		int sizeRef = vcfTokens[3].length();
		if (size < sizeRef) size = sizeRef;
		int position = Integer.parseInt(vcfTokens[1]) -1;
		return regions.search(position, position+size+2);
	}
	
	private void createIntervalTrees() {
		//load bed regions 
		Coordinate[] allRegions = Coordinate.parseFile(bedFile, 0, 0);
		Arrays.sort(allRegions);
		HashMap<String,Coordinate[]> chrRegions = Coordinate.splitByChromosome(allRegions);
		//make HashMap of trees
		chrRegionIntervalTrees = new HashMap<String,IntervalTree<Coordinate>>();
		long numRegions = 0;
		for (String chr : chrRegions.keySet()){
			Coordinate[] regions = chrRegions.get(chr);
			numRegions+= regions.length;
			ArrayList<Interval<Coordinate>> ints = new ArrayList<Interval<Coordinate>>();
			for (int i =0; i< regions.length; i++) {		
				int start = regions[i].getStart() - bpPad;
				if (start < 0) start = 0;
				ints.add(new Interval<Coordinate>(start, regions[i].getStop()+bpPad, regions[i]));
			}
			IntervalTree<Coordinate> tree = new IntervalTree<Coordinate>(ints, false);
			chrRegionIntervalTrees.put(chr, tree);
		}
		System.out.println("Loaded "+numRegions+" regions\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFRegionFilter(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFile = new File(args[++i]); break;
					case 'v': forExtraction = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'p': bpPad = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//checkfiles
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: please provide a bed file of regions for vcf intersection.\n");
		
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a vcf file or directory containing such to intersect.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");


		if (saveDirectory != null){
			saveDirectory.mkdirs();
			if (saveDirectory.isDirectory() == false || saveDirectory.exists() == false) Misc.printErrAndExit("\nCannot find or make your save directory?! "+saveDirectory);
		}
		else saveDirectory = vcfFiles[0].getParentFile();
		
	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Region Filter : Dec 2015                        **\n" +
				"**************************************************************************************\n" +
				"Sorts each vcf record (based on its position and max length of the alt or ref)\n" +
				"by intersection with a chr start stop....bed file.\n\n" +

				"Required Params:\n"+
				"-v VCF file or directory containing such (xxx.vcf(.gz/.zip OK)) to parse\n"+
				"-b Bed file of regions (chr start stop ...), interbase coordinates \n"+
				"       (xxx.bed(.gz/.zip OK)) to intersect\n"+
				
				"\nOptional Params:\n"+
				"-s Save directory for the modified vcfs\n"+
				"-p Pad bed start stop by this # bps, defaults to 0\n"+
				"\n"+

				"Example: java -Xmx9G -jar pathTo/USeq/Apps/VCFRegionFilter -v Vcfs/ \n" +
				"       -b /Anno/offTargetRegions.bed.gz -p 10 \n\n"+

		"**************************************************************************************\n");

	}
}
