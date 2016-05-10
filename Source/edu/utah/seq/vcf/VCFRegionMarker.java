package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.its.Interval;
import edu.utah.seq.its.IntervalTree;
import edu.utah.seq.useq.data.RegionScoreText;
import edu.utah.seq.useq.data.RegionScoreTextData;
import util.bio.annotation.Bed;
import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;

/**Modifies the filter field for snv and indel variants falling within the bed file with the bed file name column.
 * Good for flagging variants that hit poor quality or off target regions
 * @author Nix
 * */
public class VCFRegionMarker {

	//user fields
	private File[] vcfFiles;
	private File[] bedFiles;
	private int bpPad = 0;
	private File saveDirectory = null;
	private boolean clearFilter = false;

	//internal fields
	private HashMap<String,IntervalTree<RegionScoreText>> chrRegionIntervalTrees = null;
	private int numRecords = 0;
	private int numModifiedRecords = 0;


	//constructor
	public VCFRegionMarker(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//build interval tree for bed region file
		createIntervalTrees();

		System.out.println("Name\tTotal\tModified");
		for (int i=0; i< vcfFiles.length; i++){
			int[] totalModifiedCounts = markVcf(vcfFiles[i]);
			System.out.println(vcfFiles[i].getName()+"\t"+totalModifiedCounts[0]+"\t"+totalModifiedCounts[1]);
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private int[] markVcf(File vcf) {
		try {
			numRecords = 0;
			numModifiedRecords = 0;
			String name = Misc.removeExtension(vcf.getName());
			Gzipper modVcf = new Gzipper( new File (saveDirectory, name + "_marked.vcf.gz"));
			BufferedReader in = IO.fetchBufferedReader(vcf);
			String line;
			String currChrom = "";
			IntervalTree<RegionScoreText> regions = null;
			//for each line in the file
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#")) modVcf.println(line);
				//data line
				else {
					numRecords++;
					//#CHROM POS ID REF ALT QUAL FILTER INFO ......
					//   0    1   2  3   4   5     6   
					String[] tokens = Misc.TAB.split(line);
					//reset FILTER?
					if (clearFilter) tokens[6] = ".";
					//diff chrom?
					if (tokens[0].equals(currChrom) == false) {
						//attempt to fetch regions for this chrom, might be null
						regions = chrRegionIntervalTrees.get(tokens[0]);
						currChrom = tokens[0];

					}
					//any regions to intersect? Does it intersect?
					if (regions != null) {
						ArrayList<RegionScoreText> hits = fetchRecords(tokens, regions);

						//how many hits?
						int numHits = hits.size();
						String toAdd = null;
						//no hits so just print out unmodified
						if (numHits == 0) modVcf.println(Misc.stringArrayToString(tokens, "\t"));
						//one hit
						else if (numHits == 1) {
							toAdd = hits.get(0).getText();
							//System.out.println("\n"+hits.get(0).toStringUnderscore());
						}
						//more than one
						else toAdd = collapseNames(hits);

						//any to add
						if (toAdd != null){
							numModifiedRecords++;
							//add on to existing?
							if (tokens[6].equals(".") || tokens[6].toUpperCase().equals("PASS")) tokens[6] = toAdd;
							else tokens[6] = tokens[6] +";" + toAdd;
							line = Misc.stringArrayToString(tokens, "\t");
							modVcf.println(line);
							//System.out.println(line);
						}
					}
					//nope so just print out original
					else modVcf.println(Misc.stringArrayToString(tokens, "\t"));
				}
			}
			in.close();
			modVcf.close();
			return new int[]{numRecords, numModifiedRecords};
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: marking "+vcf);
		} 
		return null;
	}

	private ArrayList<RegionScoreText> fetchRecords(String[] vcfTokens, IntervalTree<RegionScoreText> regions) {
		//calc start stop to fetch, interbase coordinates
		int size = vcfTokens[4].length();
		int sizeRef = vcfTokens[3].length();
		if (size < sizeRef) size = sizeRef;
		int position = Integer.parseInt(vcfTokens[1]) -1;
		return regions.search(position, position+size+2);
	}

	private String collapseNames(ArrayList<RegionScoreText> hits) {
		HashSet<String> names = new HashSet<String>();
		int num = hits.size();
		for (int i=0; i< num; i++) names.add(hits.get(i).getText());
		String concat = Misc.hashSetToString(names, ";");
		return concat;
	}

	private void createIntervalTrees() {
		try {
			File bedToParse = bedFiles[0];
			//more than one bed?
			if (bedFiles.length > 1){
				ArrayList<File> beds = new ArrayList<File>();
				for (File f: bedFiles) beds.add(f);
				bedToParse = new File (saveDirectory, "tempMergedBed"+Passwords.createRandowWord(4)+ ".bed");
				bedToParse.deleteOnExit();
				IO.mergeFiles(beds, bedToParse);
			}
			//load bed regions 
			HashMap<String, RegionScoreText[]> chrRegions = Bed.parseBedFile(bedToParse, true, false);
			//make HashMap of trees
			chrRegionIntervalTrees = new HashMap<String,IntervalTree<RegionScoreText>>();
			long numRegions = 0;
			for (String chr : chrRegions.keySet()){
				RegionScoreText[] regions = chrRegions.get(chr);
				numRegions+= regions.length;
				ArrayList<Interval<RegionScoreText>> ints = new ArrayList<Interval<RegionScoreText>>();
				for (int i =0; i< regions.length; i++) {
					int start = regions[i].getStart() - bpPad;
					if (start < 0) start = 0;
					ints.add(new Interval<RegionScoreText>(start, regions[i].getStop()+bpPad, regions[i]));
					//check text
					if (regions[i].getText().length() == 0) {
						Misc.printErrAndExit("\nERROR loading the bed file, each line must contain a text field to use in adding to the FILTER field in the VCF record. See "+regions[i].getBedLine(chr));
					}
				}
				IntervalTree<RegionScoreText> tree = new IntervalTree<RegionScoreText>(ints, false);
				chrRegionIntervalTrees.put(chr, tree);
			}
			System.out.println("Loaded "+numRegions+" regions\n");

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem making interval trees from bed file(s)\n");
		}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFRegionMarker(args);
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
					case 'b': bedFiles = IO.extractFiles(args[++i]); break;
					case 'v': forExtraction = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'p': bpPad = Integer.parseInt(args[++i]); break;
					case 'c': clearFilter = true; break;
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
		if (bedFiles == null || bedFiles.length == 0) Misc.printErrAndExit("\nError: please provide one or more bed files (comma delimited, no spaces) of regions for vcf intersection.\n");

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
				"**                              VCF Region Marker : May 2016                        **\n" +
				"**************************************************************************************\n" +
				"Intersects each vcf record (based on its position and max length of the alt or ref)\n" +
				"with a chr start stop text.... bed file(s). For those that intersect, the bed text is \n"+
				"added to the vcf FILTER field. Multiple hits are concatinated.\n\n" +

				"Required Params:\n"+
				"-v VCF file or directory containing such (xxx.vcf(.gz/.zip OK)) to parse\n"+
				"-b Bed file(s) of regions (minimum chr start stop text), interbase coordinates \n"+
				"       (xxx.bed(.gz/.zip OK)) to intersect, comma delimit multiple files, no spaces.\n"+

				"\nOptional Params:\n"+
				"-s Save directory for the modified vcfs\n"+
				"-p Pad bed start stop by this # bps, defaults to 0\n"+
				"-c Clear starting FILTER field\n"+
				"\n"+

				"Example: java -Xmx9G -jar pathTo/USeq/Apps/VCFRegionMarker -v testHaploCaller.vcf.zip\n" +
				"       -b /Anno/offTargetRegions.bed.gz,/Anno/pseudogene.bed -p 10 -s MarkedVcfs -c \n\n"+

				"**************************************************************************************\n");

	}
}
