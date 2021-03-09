package edu.utah.seq.vcf.sim;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.its.Interval;
import edu.utah.seq.its.IntervalTree;
import edu.utah.seq.query.QueryIndexFileLoader;
import edu.utah.seq.useq.data.RegionScoreText;
import edu.utah.seq.useq.data.RegionScoreTextData;
import edu.utah.seq.vcf.xml.foundation.SimpleVcf;
import util.bio.annotation.Bed;
import util.bio.annotation.Coordinate;
import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;

/**Splits variants into those that intersect and those that do not.
 * @author Nix
 * */
public class VCFSelector {

	//user fields
	private File keyVcfFile;
	private File otherVcfFile;
	private File bedFile;
	private int bpPad = 150;
	private File saveDirectory = null;

	//internal fields
	private HashMap<String,IntervalTree<RegionVcfs>> chrRegionIntervalTrees = null;
	private ArrayList<RegionVcfs> regionVcfs = new ArrayList<RegionVcfs>();
	private Gzipper selectedVcfs = null;
	public static final Pattern DP = Pattern.compile(".+DP=([\\d,]+);.+");
	public int numSelectedVars = 0;


	//constructor
	public VCFSelector(String[] args){
		try {
			//start clock
			long startTime = System.currentTimeMillis();

			//process args
			processArgs(args);

			//build interval tree for collecting vcf records
			IO.pl("Creating interval trees for region vcf intersection...");
			createIntervalTrees();

			//load Interval tree with vcf records
			IO.pl("Intersecting regions with vcf records (# intersect, # don't intersect)...");
			loadTrees(keyVcfFile, true);
			if (otherVcfFile!= null) loadTrees(otherVcfFile, false);
			
			//select variants
			IO.pl("\nSelecting variants...");
			selectVariants();
			IO.pl("\tNumber selected variants "+numSelectedVars);

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");


		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void selectVariants() throws IOException {
		
		for (RegionVcfs rv: regionVcfs) {
			//any intersecting vcfs?
			int numKeyVars = rv.keyVcfs.size();
			int numOtherVars = rv.otherVcfs.size();
			if (numKeyVars == 0 && numOtherVars == 0) continue;
			
			//any key vars?
			ArrayList<SimpleVcf> al = null;
			if (numKeyVars !=0) al = selectBest(rv.keyVcfs);
			
			//any other vars?
			else if (numOtherVars != 0) al = selectBest(rv.otherVcfs);
			
			for (SimpleVcf sv: al) selectedVcfs.println( sv.getOriginalRecord());
			numSelectedVars += al.size();
		}
		selectedVcfs.close();
		
	}

	private ArrayList<SimpleVcf> selectBest(ArrayList<SimpleVcf> vcfs) {
		//just one?
		if (vcfs.size() == 1) return vcfs;

		//calc ave DP and effected bps
		for (SimpleVcf v: vcfs) {
			//DP
			Matcher mat = DP.matcher(v.getInfo());
			double ave = 0;
			if (mat.matches()) {
				double[] numbers = Num.stringArrayToDouble(mat.group(1), ",");
				ave = Num.mean(numbers);
			}
			v.setScore(ave);
			//effected bps
			int[] be = QueryIndexFileLoader.fetchEffectedBpsSingleAlt(v.getPos(), v.getRef(), v.getAlt(), null, true, false);
			v.setPadPos(be[0]);
			v.setPadEnd(be[1]);
		}

		ArrayList<SimpleVcf> al = new ArrayList<SimpleVcf>();

		SimpleVcf left = vcfs.get(0);
		SimpleVcf right = vcfs.get(vcfs.size()-1);
		//add both?
		if ((left.getPadEnd()+ bpPad) < right.getPadPos()) {
			al.add(left);
			al.add(right);
		}
		//find one with max ave depth
		else {
			SimpleVcf maxSV = vcfs.get(0);
			double max = maxSV.getScore();
			for (int i=1; i< vcfs.size(); i++) {
				if (vcfs.get(i).getScore()> max) {
					maxSV = vcfs.get(i);
					max = maxSV.getScore();
				}
			}
			al.add(maxSV);
		}
		return al;
	}

	private void loadTrees(File vcf, boolean keyRecords) {
		try {

			//IO
			String name = Misc.removeExtension(vcf.getName());
			Gzipper intVcf = new Gzipper( new File (saveDirectory, name + "_int.vcf.gz"));
			Gzipper nonIntVcf = new Gzipper( new File (saveDirectory, name + "_nonInt.vcf.gz"));
			BufferedReader in = IO.fetchBufferedReader(vcf);

			String line;
			String currChrom = "";
			int numIntRecords = 0;
			int numNonIntRecords = 0;
			IntervalTree<RegionVcfs> regions = null;

			//for each line in the file
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#")) {
					intVcf.println(line);
					nonIntVcf.println(line);
					if (keyRecords) selectedVcfs.println(line);
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
						ArrayList<RegionVcfs> hits = fetchRegions(tokens, regions);
						//no hits so just print out unmodified
						if (hits == null) {
							nonIntVcf.println(line);
							numNonIntRecords++;
						}
						else {
							intVcf.println(line);
							numIntRecords++;
							if (hits.size() == 1) {
								if (keyRecords) hits.get(0).keyVcfs.add(new SimpleVcf(line, 0));
								else hits.get(0).otherVcfs.add(new SimpleVcf(line, 0));
							}
							else {
								
								throw new Exception("\nFound a vcf record that intersected multiple "+hits.size()+" regions?\n"+line);
							}
						}
					}
					//nope so just print out as non int
					else {
						nonIntVcf.println(line);
						numNonIntRecords++;
					}
				}
			}

			in.close();
			intVcf.close();
			nonIntVcf.close();
			IO.pl("\t"+vcf.getName()+"\t"+numIntRecords+"\t"+numNonIntRecords);

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: loading "+vcf);
		} 
	}

	private ArrayList<RegionVcfs> fetchRegions(String[] vcfTokens, IntervalTree<RegionVcfs> regions) {
		//calc start stop to fetch, interbase coordinates
		int size = vcfTokens[4].length();
		int sizeRef = vcfTokens[3].length();
		if (size < sizeRef) size = sizeRef;
		int start = Integer.parseInt(vcfTokens[1]) -1;
		int stop = start+size+2;
		ArrayList<RegionVcfs> hits = regions.search(start, stop);
		if (hits.size() == 0) return null;
		
		//filter for those regions that entirely contain the vcf record, should only be one
		ArrayList<RegionVcfs> good = new ArrayList<RegionVcfs>();
		for (RegionVcfs r: hits) {
			if (start >= r.coordinate.getStart() && stop <= r.coordinate.getStop()) good.add(r);
		}
		if (good.size()==0) return null;
		return good;
	}

	private void createIntervalTrees() throws FileNotFoundException, IOException {
		//load bed regions 
		Coordinate[] allRegions = Coordinate.parseFile(bedFile, 0, 0);
		Arrays.sort(allRegions);
		HashMap<String,Coordinate[]> chrRegions = Coordinate.splitByChromosome(allRegions);

		//adjust to make padding distance apart
		Gzipper modBed = new Gzipper(new File (saveDirectory, "gapped"+bpPad+ "bps_"+ bedFile.getName()));
		for (Coordinate[] regions: chrRegions.values()) {
			regions = Coordinate.insertGaps(regions, bpPad);
			for (Coordinate c: regions) modBed.println(c);
			chrRegions.put(regions[0].getChromosome(), regions);
		}
		modBed.close();

		//make HashMap of trees
		chrRegionIntervalTrees = new HashMap<String,IntervalTree<RegionVcfs>>();
		long numRegions = 0;
		for (String chr : chrRegions.keySet()){
			Coordinate[] regions = chrRegions.get(chr);
			numRegions+= regions.length;
			ArrayList<Interval<RegionVcfs>> ints = new ArrayList<Interval<RegionVcfs>>();
			for (int i =0; i< regions.length; i++) {	
				RegionVcfs r = new RegionVcfs(regions[i]);
				ints.add(new Interval<RegionVcfs>(regions[i].getStart(), regions[i].getStop(), r));
				regionVcfs.add(r);
			}
			IntervalTree<RegionVcfs> tree = new IntervalTree<RegionVcfs>(ints, false);
			chrRegionIntervalTrees.put(chr, tree);
		}
		System.out.println("\tLoaded "+numRegions+" regions\n");
	}

	private class RegionVcfs {
		Coordinate coordinate;
		ArrayList<SimpleVcf> keyVcfs = new ArrayList<SimpleVcf>();
		ArrayList<SimpleVcf> otherVcfs = new ArrayList<SimpleVcf>();
		
		public RegionVcfs(Coordinate coordinate) {
			this.coordinate = coordinate;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFSelector(args);
	}		


	/**This method will process each argument and assign new variables
	 * @throws Exception 
	 **/
	public void processArgs(String[] args) throws Exception{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFile = new File(args[++i]); break;
					case 'k': keyVcfFile = new File(args[++i]); break;
					case 'o': otherVcfFile = new File(args[++i]); break;
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
		if (keyVcfFile == null || keyVcfFile.canRead() == false) Misc.printExit("\nError: please provide a key vcf file for selection.\n");

		if (saveDirectory == null)  Misc.printErrAndExit("\nCannot find or make your save directory?! "+saveDirectory);
		else {
			saveDirectory.mkdirs();
			if (saveDirectory.isDirectory() == false || saveDirectory.exists() == false) Misc.printErrAndExit("\nCannot find or make your save directory?! "+saveDirectory);
		}
		
		//selected vars
		File selected = new File(saveDirectory, "selectedVcfs.vcf.gz");
		selectedVcfs = new Gzipper(selected);


	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Selector : April 2019                           **\n" +
				"**************************************************************************************\n" +
				"Selects variants for injection with the BamBlaster tool. Prioritizes variants by\n"+
				"position, key vs other, and average read depth.  First run the VCFMpileupAnnotator\n" +
				"and the bam files you plan on modifying.\n"+

				"\nRequired Params:\n"+
				"-b Bed file of regions (chr start stop ...) to use in selecting intersecting vcf \n"+
				"     records, (xxx.bed(.gz/.zip OK)).\n"+
				"-k VCF file (xxx.vcf(.gz/.zip OK)) of key priority variants to attempt to inject\n"+
				"     first, e.g. annotated as pathogenic.\n"+
				"-s Save directory for the modified vcfs\n"+
				
				"\nOptional Params:\n"+
				"-o VCF file (xxx.vcf(.gz/.zip OK)) of other variants to attempt to include when \n"+
				"     priority variants aren't present.\n"+
				"-p BP distance to keep between vcf records, defaults to 150\n"+
				"\n"+

				"Example: java -Xmx9G -jar pathTo/USeq/Apps/VCFSelector -b GBMTargets.bed.gz -p 160 \n" +
				"     -k pathoClinvarSNVs.vcf.gz -o allCosmicClinvarSNVs.vcf.gz \n\n"+

				"**************************************************************************************\n");

	}
}
