package edu.cnv;
import java.io.*;

import trans.anno.*;
import util.gen.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.Region;


import util.bio.annotation.*;
import trans.roc.Gr;


public class IntersectCNVs {

	private File saveDirectory;
	private CNVGroup[] cnvGroups;
	private HashMap<String, CompositCNV[]> upCompositCNVs;
	private HashMap<String, CompositCNV[]> downCompositCNVs;
	private boolean scoreMerge = true;
	private double scoreDivider = 2;
	private LinkedHashSet<String> treatmentFileNames = new LinkedHashSet<String>();
	private LinkedHashSet<String> controlFileNames = new LinkedHashSet<String>();
	private int numberTreatmentFiles;
	private int numberControlFiles;
	private int neighborhood = 10000;
	private File ucscGeneFile;
	private File ucscCytobandFile;
	private boolean noUp = false;
	private boolean noDown = false;
	private String genomeVersion = null;

	public IntersectCNVs(String[] args){
		processArgs(args);

		//make grouping, merge or intersections
		System.out.println("Making composits...");
		if (scoreMerge) makeMergedComposites();
		else makeIntersectedComposites();

		if (noUp && noDown) Misc.printExit("\nNo up or down composite CNVs!\n");

		//intersect each CNVGroup with the targets
		System.out.println("Intersecting composits...");
		intersectCompositeCNVs();

		//load effected genes
		if (ucscGeneFile != null) {
			System.out.println("Loading effected genes...");
			loadEffectedGenes();
		}

		//load effected cytobands
		if (ucscCytobandFile != null) {
			System.out.println("Loading effected cytobands...");
			loadEffectedCytobands();
		}

		//print targets
		printTargets();

		System.out.println("\nDone!");
	}

	/**Loads intersecting genes for each target.*/
	public void loadEffectedGenes(){
		//collect all targets
		ArrayList<CompositCNV> targets = new ArrayList<CompositCNV>();
		Iterator<String> it;
		if (noUp == false){
			it = upCompositCNVs.keySet().iterator();
			while (it.hasNext()){
				CompositCNV[] t = upCompositCNVs.get(it.next());
				for (int i=0; i< t.length; i++) targets.add(t[i]);
			}
		}
		it = downCompositCNVs.keySet().iterator();
		if (noDown == false){
			while (it.hasNext()){
				CompositCNV[] t = downCompositCNVs.get(it.next());
				for (int i=0; i< t.length; i++) targets.add(t[i]);
			}
		}
		CompositCNV[] allTargets = new CompositCNV[targets.size()];
		targets.toArray(allTargets);
		//make coordinates
		Coordinate[] coor = new Coordinate[allTargets.length];
		for (int i=0; i< coor.length; i++) coor[i] = new Coordinate(allTargets[i].getChromosome(), allTargets[i].getStart(), allTargets[i].getStop());
		//find neighboring genes
		FindNeighboringGenes fng = new FindNeighboringGenes(ucscGeneFile, coor, neighborhood, true);
		String[] effGenes = fng.fetchGeneNames();
		//assign genes
		for (int i=0; i<allTargets.length; i++) allTargets[i].setEffectedGenes(effGenes[i]);

	}

	/**Loads intersecting genes for each target.*/
	public void loadEffectedCytobands(){
		//collect all targets
		ArrayList<CompositCNV> targets = new ArrayList<CompositCNV>();
		Iterator<String> it;
		if (noUp == false){
			it = upCompositCNVs.keySet().iterator();
			while (it.hasNext()){
				CompositCNV[] t = upCompositCNVs.get(it.next());
				for (int i=0; i< t.length; i++) targets.add(t[i]);
			}
		}
		it = downCompositCNVs.keySet().iterator();
		if (noDown == false){
			while (it.hasNext()){
				CompositCNV[] t = downCompositCNVs.get(it.next());
				for (int i=0; i< t.length; i++) targets.add(t[i]);
			}
		}
		CompositCNV[] allTargets = new CompositCNV[targets.size()];
		targets.toArray(allTargets);
		//make coordinates
		Coordinate[] coor = new Coordinate[allTargets.length];
		for (int i=0; i< coor.length; i++) coor[i] = new Coordinate(allTargets[i].getChromosome(), allTargets[i].getStart(), allTargets[i].getStop());
		//fetch cytobands
		NamedCoordinate[] namedCoor = NamedCoordinate.parseFile(ucscCytobandFile, 0, 0, 3); 
		Arrays.sort(namedCoor);
		HashMap<String,NamedCoordinate[]> chromNamedCoor = NamedCoordinate.splitByChromosome(namedCoor);
		//for each composite cnv look for intersected
		for (int i=0; i<allTargets.length; i++){
			//any cytobands?
			NamedCoordinate[] cytos = chromNamedCoor.get(coor[i].getChromosome());
			if (cytos!=null){
				ArrayList<String> al = new ArrayList<String>();
				for (int j=0; j< cytos.length; j++){
					if (cytos[j].intersects(coor[i]))  al.add(cytos[j].getName());
				}
				if (al.size()!=0) allTargets[i].setEffectedCytobands(Misc.stringArrayListToString(al, ", "));
			}
		}
	}

	public void printTargets(){
		//for up or down
		//for each target print: rank, chr, start, stop,#Treatment, #Control, median cnv or blank, Effected Genes
		if (noUp == false){
			File upSum = new File(saveDirectory,"upIntSum.xls");
			printResults(upSum, upCompositCNVs);
		}
		if (noDown == false){
			File downSum = new File(saveDirectory,"downIntSum.xls");
			printResults(downSum, downCompositCNVs);
		}
	}

	/**Prints: rank, chr, start, stop,#Treatment, #Control, median cnv or blank, Effected Genes*/
	public void printResults(File file, HashMap<String, CompositCNV[]> compositCNVs){
		try {
			PrintWriter out = new PrintWriter( new FileWriter (file));
			//print header
			out.print("HotLink\tChr\tStart\tStop\t#Treatment\t#Control");
			for (int i=0; i< cnvGroups.length; i++) out.print("\t"+cnvGroups[i].getName());
			if (ucscGeneFile !=null) out.print("\tEffectedGenes");
			if (ucscCytobandFile !=null) out.print("\tEffectedCytobands");
			out.println();

			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			//for each target
			Iterator<String> it = compositCNVs.keySet().iterator();
			int index = 1;
			while (it.hasNext()){
				String chrom = it.next();
				CompositCNV[] targets = compositCNVs.get(chrom);
				for (int i=0; i< targets.length; i++){
					int winStart = targets[i].getStart() - 100000;
					if (winStart < 0) winStart = 0;
					int winEnd = targets[i].getStop() + 100000;
					out.print(url+chrom+"&start="+winStart+"&stop="+winEnd+"\",\""+(index++)+"\")\t");
					out.print(chrom+"\t"+targets[i].getStart()+"\t"+targets[i].getStop()+"\t");
					//for each hit
					ArrayList<CNVHit> hits = targets[i].getHits();				
					//count number treatment and control hits
					int numTreatmentHits = 0;
					int numControlHits = 0;
					for (int x=0; x < hits.size(); x++){
						CNVHit hit = hits.get(x);
						if (treatmentFileNames.contains(hit.getCnvGroup().getName())) numTreatmentHits++;
						else if (controlFileNames.contains(hit.getCnvGroup().getName())) numControlHits++;
					}

					//print vals
					out.print(numTreatmentHits+"\t"+numControlHits+"\t");

					//collect values for all CNVGroups or a particular composite
					String[] vals = new String[numberTreatmentFiles + numberControlFiles];
					for (int x=0; x<vals.length; x++) vals[x]="";
					int counter = 0;
					//run through median value for each treatment
					Iterator<String> it2 = treatmentFileNames.iterator();
					//for each text
					while (it2.hasNext()){
						//get text
						String name = it2.next();
						double median = 0;
						//for each hit
						for (int x=0; x < hits.size(); x++){
							CNVHit hit = hits.get(x);
							//does text match?
							if (hit.getCnvGroup().getName().equals(name)){
								median = hit.getMedianCNV();
								vals[counter] = Num.formatNumber(median, 3);
							}
						}
						counter++;
					}
					//run through median value for each control
					it2 = controlFileNames.iterator();
					//for each text
					while (it2.hasNext()){
						//get text
						String name = it2.next();
						double median = 0;
						//for each hit
						for (int x=0; x < hits.size(); x++){
							CNVHit hit = hits.get(x);
							//does text match?
							if (hit.getCnvGroup().getName().equals(name)){
								median = hit.getMedianCNV();
								vals[counter] = Num.formatNumber(median, 3);
							}
						}
						counter++;
					}
					//print it
					out.print(Misc.stringArrayToString(vals, "\t"));

					//print effected genes and cytobands?
					if (ucscGeneFile !=null) out.print("\t"+targets[i].getEffectedGenes());
					if (ucscCytobandFile !=null) out.print("\t"+targets[i].getEffectedCytobands());
					out.println();
				}
			}

			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void intersectCompositeCNVs(){
		//for each cnvGroup
		for (int i=0; i< cnvGroups.length; i++){
			//for each cnv, intersect with targets
			CNV[] upCNVs = cnvGroups[i].getUpCNVs();
			if (upCNVs != null) intersect(cnvGroups[i], upCNVs, upCompositCNVs);
			CNV[] downCNVs = cnvGroups[i].getDownCNVs();
			if (downCNVs != null) intersect(cnvGroups[i], downCNVs, downCompositCNVs);
		} 
	}

	public void intersect(CNVGroup cnvGroup, CNV[] cnvs, HashMap<String, CompositCNV[]> compositCNVs){
		//split cnvs by chromosome
		HashMap<String, CNV[]> split = CNV.splitByChromsome(cnvs);

		//for each target chromosome
		Iterator<String> it = compositCNVs.keySet().iterator();
		while(it.hasNext()){
			String chrom = it.next();
			CompositCNV[] t = compositCNVs.get(chrom);
			CNV[] c = split.get(chrom);
			if (c == null) continue;
			//intersect and save hits
			//for each compositCNVs
			for (int i=0; i< t.length; i++){
				ArrayList<Float> hits = new ArrayList<Float>();
				//fetch hits
				for (int j=0; j< c.length; j++){
					ArrayList<Float> subHits = intersect(t[i], c[j]);
					if (subHits != null) hits.addAll(subHits);
				}
				//calculate median and set CNVHit?
				if (hits.size() !=0){
					float[] vals = Num.arrayListOfFloatToArray(hits);
					Arrays.sort(vals);
					double median = Num.median(vals);
					t[i].getHits().add(new CNVHit(cnvGroup, median));
				}

			}
		}

	}

	/**Assumes same chromosome.*/
	public ArrayList<Float> intersect (CompositCNV t, CNV c){
		//do they intersect?
		if (c.stop < t.getStart() || c.start > t.getStop()) return null;
		//find intersecting values and take median
		Gr[] slice = c.fetchObservations(t.getStart(), t.getStop());
		ArrayList<Float> vals = new ArrayList<Float>(slice.length);
		for (int i=0; i< slice.length; i++) vals.add(new Float(slice[i].getScore()));
		return vals;
	}

	/**Takes CNVs and merges them to make the CompositCNVs.*/
	public void makeMergedComposites(){
		//make regions
		ArrayList<GenomicRegion[]> upAL = new ArrayList<GenomicRegion[]>();
		ArrayList<GenomicRegion[]> downAL = new ArrayList<GenomicRegion[]>();
		RegionComparator comp = new RegionComparator();
		for (int i=0; i< cnvGroups.length; i++){
			//up regions
			CNV[] up = cnvGroups[i].getUpCNVs();
			if (up != null ) {
				GenomicRegion[] upRegions = new GenomicRegion[up.length];
				for (int j=0; j< up.length; j++){
					upRegions[j] = new GenomicRegion(up[j].getGrGraph().getChromosome(), up[j].getStart(), up[j].getStop(), null);
				}
				Arrays.sort(upRegions,comp);
				upAL.add(upRegions);
			}

			//down regions
			CNV[] down = cnvGroups[i].getDownCNVs();
			if (down != null ) {
				GenomicRegion[] downRegions = new GenomicRegion[down.length];
				for (int j=0; j< down.length; j++){
					downRegions[j] = new GenomicRegion(down[j].getGrGraph().getChromosome(), down[j].getStart(), down[j].getStop(), null);
				}
				Arrays.sort(downRegions,comp);
				downAL.add(downRegions);
			}
		}
		//up
		noUp = true;
		if (upAL.size() !=0){
			noUp = false;
			GenomicRegion[][] up = new GenomicRegion[upAL.size()][];
			for (int i=0; i< upAL.size(); i++) up[i] = upAL.get(i);
			File upMerge = new File (saveDirectory, "upMergedCNVs.bed");
			printMergedRegions (up, upMerge);
			GenomicRegion[] upRegions = GenomicRegion.parseRegions(upMerge);
			HashMap<String, GenomicRegion[]> upChrRegions = GenomicRegion.splitByChromosome(upRegions);
			upCompositCNVs = makeCompositCNVs(upChrRegions);
		}
		//down
		noDown = true;
		if (downAL.size() !=0){
			noDown = false;
			GenomicRegion[][] down = new GenomicRegion[downAL.size()][];
			for (int i=0; i< downAL.size(); i++) down[i] = downAL.get(i);
			File downMerge = new File (saveDirectory, "downMergedCNVs.bed");
			printMergedRegions (down, downMerge);
			GenomicRegion[] downRegions = GenomicRegion.parseRegions(downMerge);
			HashMap<String, GenomicRegion[]> downChrRegions = GenomicRegion.splitByChromosome(downRegions);
			downCompositCNVs = makeCompositCNVs(downChrRegions);
		}
	}


	public HashMap<String, CNV[]> mergeAndSplitCNVs (boolean fetchUp){
		//make merge
		ArrayList<CNV> cnvsAL= new ArrayList<CNV>();
		for (int i=0; i< cnvGroups.length; i++){
			CNV[] c;
			if (fetchUp) c = cnvGroups[i].getUpCNVs();
			else c = cnvGroups[i].getDownCNVs();
			if (c != null) for (int j=0; j< c.length; j++) cnvsAL.add(c[j]);
		}
		CNV[] allCNVs = new CNV[cnvsAL.size()];
		cnvsAL.toArray(allCNVs);
		//sort
		Arrays.sort(allCNVs, new ComparatorCNVChromPosition());
		//split
		return CNV.splitByChromsome(allCNVs);

	}

	public static int findLastBase(CNV[] chromSpecific){
		int max = 0;
		for (int i=0; i< chromSpecific.length; i++){
			if (chromSpecific[i].getStop() > max) max = chromSpecific[i].getStop();
		}
		return max;
	}

	public HashMap<String, GenomicRegion[]> buildChromRegions(HashMap<String, CNV[]> chromMerge){
		HashMap<String, GenomicRegion[]> chromRegions = new HashMap<String, GenomicRegion[]>();
		//for each chromosome
		Iterator<String> it = chromMerge.keySet().iterator();
		while (it.hasNext()){
			//fetch CNVs
			String chrom = it.next();
			CNV[] cnvs = chromMerge.get(chrom);
			//find last base
			int max = findLastBase(cnvs);
			//make array to hold number of hits
			short[] baseHitCount = new short[max];
			//load
			for (int i=0; i< cnvs.length; i++){
				int start = cnvs[i].getStart();
				int stop = cnvs[i].getStop();
				for (int j = start; j< stop; j++) baseHitCount[j]++;
			}
			//set any with less than 2 hits to zero
			for (int i=0; i< baseHitCount.length; i++){
				if (baseHitCount[i] < 2) baseHitCount[i] = 0;
			}
			//fetch start and stops
			Region[] startStops = Region.makeStartStops(baseHitCount);
			if (startStops == null) continue;
			//make regions
			GenomicRegion[] regions = new GenomicRegion[startStops.length];
			for (int i=0; i< startStops.length; i++){
				regions[i] = new GenomicRegion(chrom, startStops[i].getStart(), startStops[i].getStop(), null);
			}
			//add to hash
			chromRegions.put(chrom, regions);
		}
		return chromRegions;
	}



	/**Takes CNVs and intersects them to make the CompositCNVs. Writes intersected regions to bed file.*/
	public void makeIntersectedComposites(){
		noUp = true;
		//up
		HashMap<String, CNV[]> upMerge = mergeAndSplitCNVs(true);
		HashMap<String, GenomicRegion[]> upRegions = buildChromRegions(upMerge);
		if (upRegions.size() != 0) {
			noUp = false;
			File up = new File (saveDirectory, "upIntCNVs.bed");
			GenomicRegion.writeRegions(upRegions, up);
			upCompositCNVs = makeCompositCNVs(upRegions);
		}
		//down
		noDown = true;
		HashMap<String, CNV[]> downMerge = mergeAndSplitCNVs(false);
		HashMap<String, GenomicRegion[]> downRegions = buildChromRegions(downMerge);
		if (downRegions.size() !=0 ){
			noDown = false;
			File down = new File (saveDirectory, "downIntCNVs.bed");
			GenomicRegion.writeRegions(downRegions, down);
			downCompositCNVs = makeCompositCNVs(downRegions);
		}
	}

	public HashMap<String, CompositCNV[]> makeCompositCNVs(HashMap<String, GenomicRegion[]> regions){
		HashMap<String, CompositCNV[]> targets = new HashMap<String, CompositCNV[]>();
		Iterator<String> it = regions.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			GenomicRegion[] r = regions.get(chrom);
			CompositCNV[] t = new CompositCNV[r.length];
			for (int i=0; i< r.length; i++){
				t[i] = new CompositCNV(r[i].getChromosome(), r[i].getStart(), r[i].getEnd());
			}
			targets.put(chrom, t);
		}
		return targets;
	}


	private void printMergedRegions(GenomicRegion[][] regions, File mergedFile){
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
							int stop = regions[i][j].getEnd()+1;
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
	private static void print(String chromosome, boolean[] bps, PrintWriter out){
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
		new IntersectCNVs(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File[] treatmentCNVs = null;
		File[] controlCNVs = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatmentCNVs = IO.extractFiles(args[++i], ".cnv"); break;
					case 'c': controlCNVs = IO.extractFiles(args[++i], ".cnv"); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'g': ucscGeneFile = new File(args[++i]); break;
					case 'y': ucscCytobandFile = new File(args[++i]); break;
					case 'n': neighborhood = Integer.parseInt(args[++i]); break;
					case 'd': scoreDivider = Double.parseDouble(args[++i]); break;
					case 'i': scoreMerge = false; break;
					case 'v': genomeVersion = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for params
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();
		if (treatmentCNVs == null) Misc.printExit("\nError: enter a directory containing treatment xxx.cnv files.\n");
		if (controlCNVs == null) Misc.printExit("\nError: enter a directory containing control xxx.cnv files.\n");
		if (genomeVersion== null) Misc.printExit ("\nPlease enter a genome version, i.e. H_sapiens_Mar_2006\n");

		//collect cnv files
		numberTreatmentFiles = treatmentCNVs.length;
		String[] names = IO.fetchFileNamesNoExtensions(treatmentCNVs);
		for (int i=0; i< names.length; i++) treatmentFileNames.add(names[i]);

		numberControlFiles = controlCNVs.length;
		names = IO.fetchFileNamesNoExtensions(controlCNVs);
		for (int i=0; i< names.length; i++) controlFileNames.add(names[i]);

		//make CNVGroups
		cnvGroups = new CNVGroup[numberTreatmentFiles+numberControlFiles];
		int index = 0;
		for (int i=0; i< numberTreatmentFiles; i++) cnvGroups[index++] = new CNVGroup(treatmentCNVs[i], scoreDivider);
		for (int i=0; i< numberControlFiles; i++) cnvGroups[index++] = new CNVGroup(controlCNVs[i], scoreDivider);

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Intersect CNVs: January 2009                           **\n" +
				"**************************************************************************************\n" +
				"Takes scored lists of CNV from the CNVScanner and identifies the common intersection.\n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-t Full path directory containing treatment xxx.cnv files. See CNVScanner.\n" +
				"-c Full path directory containing control xxx.cnv files. \n" +
				"-v Genome version, i.e. 'H_sapiens_Mar_2006', see IGB or the UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-g Optional, UCSC RefFlat gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables\n"+
				"-y Optional, UCSC cytoband file (tabbed: chr,start,stop,text), full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables\n"+
				"-n Size of neighborhood to use when finding neighboring genes, defaults to 10000\n"+
				"-d Score divider for seperating up vs down cnvs, defaults to 2\n"+
				"-i Score the common intersect (exclusive) instead of the common merge (inclusive).\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/CNVAn/Apps/IntersectCNVs -s /CnvResults/ -t\n" +
				"       /CnvLeuk/ -c /CnvNormal/ -v H_sapiens_May_2004 -g /Anno/hg18RefSeq.txt -g\n" +
				"       /Anno/hg18Cytobands.txt -n 20000 -i\n\n" +

		"**************************************************************************************\n");

	}

}
