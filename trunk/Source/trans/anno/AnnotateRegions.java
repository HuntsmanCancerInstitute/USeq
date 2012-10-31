package trans.anno;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import trans.main.*;
import util.bio.annotation.*;
import util.bio.parsers.gff.*;
import util.gen.*;


/**
 * For annotating Intervals or binding regions with gene information, proximity, neighbors, and bias based on random trials.
 * Specific for dmel release 4.0 gff3.  Each gff3 file needs careful reworking to get it into a form for making gene models.
 * 
 * @see AnnotateRegionsWithGeneList
 */
public class AnnotateRegions {
	
	//fields
	private File gffFile;
	private File picksFile;
	private File dgcFile;
	private GeneGroup[] geneGrps;
	private BindingRegion[] bindingRegions;
	private int numRandomTrials = 0;
	private boolean stripGFF = true;
	private boolean printNumberOfNeighbors = false;
	private int sizeNeighborhood = 10000;
	
	public AnnotateRegions(String[] args){
		processArgs(args);
		System.out.println("\nLaunching Annotator...");
		
		//process gff file retrieving an array of GeneGroup
		System.out.println("\tProcessing GFF file...");
		DmelRel4Extractor ex = new DmelRel4Extractor();
		ex.extract(gffFile, stripGFF);
		ArrayList geneGrpsAL = ex.getGeneGroupArrayList();
		int num = geneGrpsAL.size();
		
		if (stripGFF){
			//run thru ArrayList pitching any without a translation or not a gene
			ArrayList filteredGG = new ArrayList();
			GeneGroup test;
			for (int i=0; i<num; i++){
				test = (GeneGroup)geneGrpsAL.get(i);
				//start with CG?
				if (test.getType().equals("gene")) {
					//has a translation
					if (test.getGeneRep().isTranslationFlag()){
						//reset the orientation
						//if (Math.IEEEremainder(i,2)!=0) test.setOrientation(-1);
						//else test.setOrientation(1);
						filteredGG.add(test);
					}
				}
			}
			geneGrps = new GeneGroup[filteredGG.size()];
			filteredGG.toArray(geneGrps);
		}
		else{
			geneGrps = new GeneGroup[num];
			geneGrpsAL.toArray(geneGrps);
		}
		
		//filter by CG file?
		if (dgcFile!=null){
			HashSet names = IO.loadFileIntoHashSet(dgcFile);
			geneGrps = filterGeneGroups(geneGrps, names);
		}
		System.out.println("\t"+(geneGrps.length+1)+" Filtered Gene Groups!");
		
		//process picks file
		System.out.println("\tProcessing picks/Interval file...");
		bindingRegions = parseIntervalFile(picksFile, sizeNeighborhood);
		if (bindingRegions == null) bindingRegions = parsePicksFile(picksFile, sizeNeighborhood);
		
		//cross to compare
		System.out.println("\tComparing Binding Regions against Gene Groups...");
		compareBindingRegionsVsGeneGrps(geneGrps, bindingRegions);
		
		//sort
		Arrays.sort(bindingRegions);
		
		//print results
		System.out.println("\nBinding Regions:");
		int numBindingRegions = bindingRegions.length;
		for (int i=0; i<numBindingRegions; i++){
			System.out.println(bindingRegions[i]);
			System.out.println("*****************************");
		}
		System.out.println();
		printDistToClosestATGAndTranscript(bindingRegions);
		
		System.out.println("\nBinding Regions Summary Line: rank, chrom, start, stop, dist to closes gene, gene text(s), # neighbors\n");
		System.out.println("Rank #\tChrom\tStart\tEnd\tDist to closest gene\tGene text(s)\t# neighbors within ("+sizeNeighborhood+"bp)");
		for (int i=0; i<numBindingRegions; i++){
			System.out.println(bindingRegions[i].summaryLine());
		}
		
		
		
		//print list of neighbors, cg numbers
		System.out.println("\nList of all genes within "+sizeNeighborhood+"bp of the binding regions.");
		HashSet neighbors = new HashSet(numBindingRegions * 4);
		
		for (int i=0; i<numBindingRegions; i++){
			neighbors.addAll(extractCGNames(bindingRegions[i].getNeighboringGeneGrps()));
		}
		System.out.println(neighbors);
		
		System.out.println("\n");
		
		//print list of immediate neighbors, cg numbers
		System.out.println("List of the closest neighbors. If no genes overlap the binding region/peak, then the closes 5' and the closest 3' genes are " +
				"saved with their distance, provided they are within the "+sizeNeighborhood+"bp neighborhood; otherwise, only the overlapping/ containing genes are saved.");
		HashSet closestNeighbors = new HashSet(numBindingRegions);
		GeneGroup gg;
		ArrayList names;
		boolean foundOverlap;
		for (int i=0; i<numBindingRegions; i++){
			foundOverlap = false;
			names = extractCGNames(bindingRegions[i].getContainingGeneGrps());
			if (names.size()!=0){
				closestNeighbors.addAll(names);
				foundOverlap = true;
			}
			names = extractCGNames(bindingRegions[i].getOverlap3PrimeGeneGrps());
			if (names.size()!=0){
				closestNeighbors.addAll(names);
				foundOverlap = true;
			}
			names = extractCGNames(bindingRegions[i].getOverlap5PrimeGeneGrps());
			if (names.size()!=0){
				closestNeighbors.addAll(names);
				foundOverlap = true;
			}
			if (foundOverlap == false ){
				if (bindingRegions[i].getDistanceTo3PrimeGeneGrp()<=bindingRegions[i].getNeighborhood()){
					gg = bindingRegions[i].getGeneGrp3Prime();
					if(gg!=null) closestNeighbors.add(gg.getName()+"_"+bindingRegions[i].getDistanceTo3PrimeGeneGrp());
				}
				if (bindingRegions[i].getDistanceTo5PrimeGeneGrp()<=bindingRegions[i].getNeighborhood()){
					gg = bindingRegions[i].getGeneGrp5Prime();
					if(gg!=null) closestNeighbors.add(gg.getName()+"_"+bindingRegions[i].getDistanceTo5PrimeGeneGrp());
				}
			}
		}
		System.out.println(closestNeighbors);
		System.out.println("\n"+numBindingRegions+" Number of binding regions.");
		System.out.println(neighbors.size()+" Number of unique neighbors.");
		System.out.println(closestNeighbors.size()+" Number of closest neighbors.");
		
		//count flanking genes
		int[] sides = countGenes (bindingRegions);
		System.out.println(sides[0]+" Number of genes with a binding region on their 5' stop.");
		System.out.println(sides[1]+" Number of genes with a binding region on their 3' stop.");
		System.out.println(sides[2]+" Number of genes with a binding region overlapping their 5' stop.");
		System.out.println(sides[3]+" Number of genes with a binding region overlapping their 3' stop.");
		System.out.println(sides[4]+" Number of genes entirely containing a binding region.");
		System.out.println(sides[5]+" Number of binding regions with neighbors.");
		System.out.println(sides[6]+" Number of binding regions with no neighbors.");
		System.out.println(sides[7]+" Average size of a binding region.");
		System.out.println(sides[8]+" Number of regions in non coding DNA.");
		System.out.println(sides[9]+" Number of regions in coding DNA.");
		System.out.println(sides[10]+" Number of regions that overlap both coding and non coding DNA.\n");
		
		//random trials?
		if (numRandomTrials !=0){
			if (printNumberOfNeighbors) System.out.println("Random trials, number of binding regions with neighbors within "+sizeNeighborhood+"bp of the binding regions.");
			else System.out.println("Random trials, column order as above.");
			//chrom info
			HashMap chromLengths = new HashMap(6);
			chromLengths.put("2L",new Integer(22407834));
			chromLengths.put("2R",new Integer(20766785));
			chromLengths.put("3L",new Integer(23771897));
			chromLengths.put("3R",new Integer(27905053));
			chromLengths.put("4",new Integer(1281640));
			chromLengths.put("X",new Integer(22224390));
			BindingRegion[] random;
			int[] counts;
			for (int i=0; i<numRandomTrials; i++){
				random = makeRandomBindingRegions(bindingRegions, chromLengths, sizeNeighborhood);
				compareBindingRegionsVsGeneGrps(geneGrps, random);
				if (printNumberOfNeighbors){
					System.out.println(countNumberBindingRegionsWithNeighbors(random));
				}
				else{
					counts = countGenes (random);
					System.out.println(counts[0]+"\t"+counts[1]+"\t"+counts[2]+"\t"+counts[3]+"\t"+counts[4]+
							"\t"+counts[5]+"\t"+counts[6]+"\t"+counts[7]+"\t"+counts[8]+"\t"+counts[9]+"\t"+counts[10]);
				}
			}
		}
	}
	
	/**For each binding region this will make another binding region from the same chromosome with the 
	 * same length, yet at a random location.*/
	public static BindingRegion[] makeRandomBindingRegions(BindingRegion[] br, HashMap chromLengths, int sizeNeighborhood){
		int num = br.length;
		BindingRegion[] random = new BindingRegion[num];
		for (int i=0; i<num; i++){
			int end = ((Integer)(chromLengths.get(br[i].getChromosome()))).intValue();
			boolean go = true;
			int start = 0;
			int regionLength = br[i].getEnd()-br[i].getStart();
			while (go){
				end = (int)Math.round(Math.random()*(double)end);
				start = end - regionLength;
				if (start > 0) go = false;
				else end = ((Integer)(chromLengths.get(br[i].getChromosome()))).intValue();
			}
			random[i]= new BindingRegion(i, 0,"chr"+br[i].getChromosome(), start, end, sizeNeighborhood);
		}
		return random;
	}
	
	/**Prints rank, chrom, start, stop, distance to closest ATG, to closest transcript start.*/
	public void printDistToClosestATGAndTranscript(BindingRegion[] br){
		int num = br.length;
		int atg =0;
		int trans = 0;
		for (int i=0; i< num; i++) {
			atg = findDistToClosestATG(br[i]);
			trans = findDistToClosestTranscript(br[i]);
			System.out.println(
					br[i].getRank()+"\t"+
					br[i].getChromosome()+"\t"+
					br[i].getStart()+"\t"+
					br[i].getEnd()+"\t"+
					atg+"\t"+
					trans
			);
		}
	}
	
	/**Finds the distance to the closest ATG translation start site.*/
	public int findDistToClosestATG(BindingRegion br){
		int dist = 1000000000;
		int test;
		//overlap 3'
		test = findDistToClosestATG(br, br.getOverlap3PrimeGeneGrps());
		if (test == 0) return 0;
		if (test<dist) dist = test;
		//overlap 5'
		test = findDistToClosestATG(br, br.getOverlap5PrimeGeneGrps());
		if (test == 0) return 0;
		if (test<dist) dist = test;
		//contained
		test = findDistToClosestATG(br, br.getContainingGeneGrps());
		if (test == 0) return 0;
		if (test<dist) dist = test;		
		//5'
		test = findDistanceToATG(br, br.getGeneGrp5Prime());
		if (test == 0) return 0;
		if (test<dist) dist = test;
		//3'
		test = findDistanceToATG(br, br.getGeneGrp3Prime());
		if (test == 0) return 0;
		if (test<dist) dist = test;
		return dist;
	}
	
	public int findDistToClosestATG(BindingRegion br, ArrayList geneGroups){
		int num = geneGroups.size();
		int dist = 1000000000;
		for (int i=0; i<num; i++){
			GeneGroup gp = (GeneGroup)geneGroups.get(i);
			int testDist = findDistanceToATG(br, gp);		
			if (testDist<dist) dist = testDist;
		}
		return dist;
	}
	
	/**Finds the distance to the closest ATG translation start site.*/
	public int findDistToClosestTranscript(BindingRegion br){
		int dist = 1000000000;
		int test;
		//overlap 3'
		test = findDistToClosestTranscript(br, br.getOverlap3PrimeGeneGrps());
		if (test == 0) return 0;
		if (test<dist) dist = test;
		//overlap 5'
		test = findDistToClosestTranscript(br, br.getOverlap5PrimeGeneGrps());
		if (test == 0) return 0;
		if (test<dist) dist = test;
		//contained
		test = findDistToClosestTranscript(br, br.getContainingGeneGrps());
		if (test == 0) return 0;
		if (test<dist) dist = test;		
		//5'
		GeneGroup gp = br.getGeneGrp5Prime();
		if (gp!=null){
			test = findDistToClosestTranscript(br, gp);
			if (test == 0) return 0;
			if (test<dist) dist = test;
		}
		//3'
		gp = br.getGeneGrp3Prime();
		if (gp!=null){
			test = findDistToClosestTranscript(br, gp);
			if (test == 0) return 0;
			if (test<dist) dist = test;
		}
		return dist;
	}
	
	public int findDistToClosestTranscript(BindingRegion br, ArrayList geneGroups){
		int num = geneGroups.size();
		int dist = 1000000000;
		for (int i=0; i<num; i++){
			GeneGroup gp = (GeneGroup)geneGroups.get(i);
			int testDist = findDistToClosestTranscript(br, gp);
			if (testDist<dist) dist = testDist;
		}
		return dist;
	}
	
	/**Finds the distance to conservative estimate of an ATG, returns 0 if overlaps.*/
	public int findDistanceToATG(BindingRegion br, GeneGroup gp){
		//1 is forward (start stop 5..25), -1 is reverse (start stop 25..5)
		int ori = gp.getOrientation();
		ArrayList codingSegments = gp.getGeneRep().getCodingSegments();
		int num = codingSegments.size();
		//find atg
		int atgPosition;
		if (ori == 1) atgPosition = ((int[])codingSegments.get(0))[0];
		else atgPosition = ((int[])codingSegments.get(num-1))[1];
		
		//System.out.println(gp.getName()+" atgPosition "+atgPosition+" brstart "+br.getStart()+" brstop "+br.getEnd());		
		//contained within
		if (atgPosition>=br.getStart() && atgPosition<= br.getEnd()) return 0;
		//left
		if (atgPosition<br.getStart()) return br.getStart()-atgPosition;
		//right, should work by default
		if (atgPosition>br.getEnd()) return atgPosition- br.getEnd();
		//return atgPosition- br.getEnd();
		System.out.println("problem with findDistanceToATG()");
		return -1;
	}
	
	/**Finds the distance to conservative estimate of start of first exon, returns 0 if overlaps.*/
	public int findDistToClosestTranscript(BindingRegion br, GeneGroup gp){
		//1 is forward (start stop 5..25), -1 is reverse (start stop 25..5)
		int ori = gp.getOrientation();
		int startStop[] = gp.getGeneRep().getFivePrimeNonCodingRegion();
		int start;
		if (ori == 1) start = startStop[0];
		else start = startStop[1];
		//contained within
		if (start>=br.getStart() && start<= br.getEnd()) return 0;
		//left
		if (start<br.getStart()) return br.getStart()-start;
		//right, should work by default
		if (start>br.getEnd()) return start- br.getEnd();
		//return start- br.getEnd();
		System.out.println("problem with findDistanceToStartExon()");
		return -1;
	}
	
	public static int countNumberNeighbors(BindingRegion[] br){
		int num = br.length;
		int numNeighbors =0;
		for (int i=0; i< num; i++) numNeighbors += br[i].getNeighboringGeneGrps().size();
		return numNeighbors;
	}
	
	public static int countNumberBindingRegionsWithNeighbors(BindingRegion[] br){
		int num = br.length;
		int numNeighbors =0;
		for (int i=0; i< num; i++) {
			if (br[i].getNeighboringGeneGrps().size()!=0) numNeighbors++;
		}
		return numNeighbors;
	}
	
	/**Returns the number of genes where the binding region is on the 5' stop and
	 * the number of genes where the binding region is on the 3' stop of the respective gene,
	 * the number of genes that overlap a binding region on their 5' stop and 3' stop, lastly 
	 * the number of regions entirely contained by a gene,
	 * the number of binding regions with neighbors,
	 * the number of regions with no neighbors as defined by the neighborhood,
	 * the number of regions in non coding DNA,
	 * the number of regions in coding DNA,
	 * the number of regions that overlap coding and nonCoding DNA
	 * @return int[9] {num 5', num 3', overlap 5', overlap 3', contained, no neighbors, 
	 * 		non coding, coding, overlap coding and non coding}*/
	public static int[] countGenes (BindingRegion[] br){
		int num = br.length;
		int numFivePrime = 0;
		int numThreePrime = 0;
		int numFiveOverlap = 0;
		int numThreeOverlap = 0;
		int numContained =0;
		int numNoNeighbors =0;
		int numBRsWithNeighbors =0;
		int totalSize =0;
		int numRegionsNonCoding =0;
		int numRegionsCoding =0;
		int numRegionsOverlapCNC =0;
		GeneGroup gg;
		int ori;
		for (int i=0; i< num; i++){
			gg = br[i].getGeneGrp5Prime();
			//5' 3' flanking
			if (gg!=null) {
				ori = gg.getOrientation();
				if (ori == 1) numThreePrime++;
				else numFivePrime++;
			}
			gg= br[i].getGeneGrp3Prime();
			if (gg!=null){
				ori = gg.getOrientation();
				if (ori == -1) numThreePrime++;
				else numFivePrime++;
			}
			//overlaps
			ArrayList over = br[i].getOverlap5PrimeGeneGrps();
			for (int j=0; j<over.size(); j++){
				gg = (GeneGroup)over.get(j);
				ori = gg.getOrientation();
				if (ori ==1 ) numThreeOverlap ++;
				else numFiveOverlap ++;
			}
			over = br[i].getOverlap3PrimeGeneGrps();
			for (int j=0; j<over.size(); j++){
				gg = (GeneGroup)over.get(j);
				ori = gg.getOrientation();
				if (ori == -1 ) numThreeOverlap ++;
				else numFiveOverlap ++;
			}
			//contained
			numContained += br[i].getContainingGeneGrps().size();
			//neighbors
			if (br[i].getNeighboringGeneGrps().size()==0) numNoNeighbors++;
			if (br[i].getNeighboringGeneGrps().size() !=0) numBRsWithNeighbors++;
			//average size
			totalSize += (br[i].getEnd()+1-br[i].getStart());
			//where does the binding region fall coding vs non coding vs overlap both?
			//is it entirely non coding?			
			if (br[i].getOverlap3PrimeGeneGrps().size()==0 && 
					br[i].getOverlap5PrimeGeneGrps().size()==0 && 
					br[i].getContainingGeneGrps().size()==0){
				numRegionsNonCoding++;
			}
			//overlaps or is contained by a gene, now must check coding non coding segments of every overlapping gene group
			else{				
				//assemble set of all overlapping gene groups
				ArrayList allGeneGroups = new ArrayList();
				allGeneGroups.addAll(br[i].getOverlap5PrimeGeneGrps());
				allGeneGroups.addAll(br[i].getOverlap3PrimeGeneGrps());
				allGeneGroups.addAll(br[i].getContainingGeneGrps());			
				//assemble set of all coding and all nonCoding
				ArrayList nonCodingInts = new ArrayList();
				ArrayList codingInts = new ArrayList();
				int numGGs= allGeneGroups.size();
				//run thru each gene group
				for (int j=0; j<numGGs; j++){ 
					gg = (GeneGroup) allGeneGroups.get(j);
					//System.out.println(gg);					
					GeneRep rep = gg.getGeneRep();
					//add nonCoding segments
					nonCodingInts.addAll(rep.getNonCodingSegments());
					if (rep.getFivePrimeNonCodingRegion()!=null)nonCodingInts.add(rep.getFivePrimeNonCodingRegion());
					if (rep.getThreePrimeNonCodingRegion()!=null)nonCodingInts.add(rep.getThreePrimeNonCodingRegion());
					//add coding segments
					codingInts.addAll(rep.getCodingSegments());
				}
				//check for overlap with segments
				boolean overlapsNonCoding = false;
				boolean overlapsCoding = false;
				
				if ( overlap( nonCodingInts, br[i].getStart(), br[i].getEnd() )  ) overlapsNonCoding = true;
				if ( overlap( codingInts, br[i].getStart(), br[i].getEnd() )  ) overlapsCoding = true;
				//increment counters
				if (overlapsNonCoding == true && overlapsCoding == false) {
					numRegionsNonCoding++;
				}
				//also case where non noncoding or coding segments, count as coding, pertains to CRMs
				else if ( (overlapsNonCoding == false && overlapsCoding == true) || (overlapsNonCoding==false && overlapsCoding==false)) {
					numRegionsCoding++;
				}
				else if (overlapsNonCoding == true && overlapsCoding == true) {
					numRegionsOverlapCNC++;
				}
				else{
					System.out.println("Problem with overlap counting in AnnotateRegions!");
					System.exit(1);
				}
				
			}
			
		}
		int aveSizeRegion = Math.round(totalSize/num);
		return new int[]{numFivePrime, numThreePrime, numFiveOverlap, numThreeOverlap, numContained, numBRsWithNeighbors,
				numNoNeighbors,aveSizeRegion, numRegionsNonCoding, numRegionsCoding, numRegionsOverlapCNC};
	}
	
	/**Tests whether any startEnd int[] in the ArrayList of ints ovelaps a region defined by the
	 * startRegion and endRegion. Assumes start is always <= stop.*/
	public static boolean overlap(ArrayList ints, int startRegion, int endRegion){
		int num = ints.size();
		int[]startEndSeg;
		for (int i=0; i<num; i++){
			startEndSeg = (int[])ints.get(i);
			//skip zero and 1bp size segments from 5' and 3' noncoders
			if (startEndSeg[1]-startEndSeg[0]<2) continue;
			//no overlap
			if (startRegion> startEndSeg[1] || endRegion< startEndSeg[0]){}
			//must overlap
			else {
				return true;
			}
		}
		return false;
	}
	
	/**Extracts the names of each gene group returning an ArrayList of Strings.*/
	public static ArrayList extractCGNames(ArrayList geneGroups){
		int numGeneGroups = geneGroups.size();
		ArrayList names = new ArrayList(numGeneGroups);
		GeneGroup gg;
		for (int i=0; i<numGeneGroups; i++){
			gg = (GeneGroup)geneGroups.get(i);
			if (gg!=null) names.add( gg.getName() );
		}
		return names;
	}
	
	/**Does a complete scan, could be optimized.*/
	public static void compareBindingRegionsVsGeneGrps(GeneGroup[] geneGroups, BindingRegion[] bindingRegions){
		int numGeneGrps = geneGroups.length;
		int numBindingRegions = bindingRegions.length;
		//for each binding region, scan it against gene groups
		for (int i=0; i<numBindingRegions; i++){
			for (int j=0; j<numGeneGrps; j++){
				bindingRegions[i].compare(geneGroups[j]);
			}
		}
	}
	
	public static BindingRegion[] parsePicksFile(File picksFile, int sizeNeighborhood){
		BindingRegion[] regions =null;
		try{
			BufferedReader in = new BufferedReader(new FileReader(picksFile));
			String line;
			String[] tokens;
			ArrayList regionsAL = new ArrayList();
			int rank = 0;
			String chromosome;
			int start;
			int stop;
			double score;
			//rank, subWinMedianScore, trimmedMeanPickScore, chrom, start, stop, seq
			//  0           1                     2            3      4     5     6
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.equals("")) continue;
				tokens = line.split("\\s+");
				if (tokens.length<7) continue;
				rank = Integer.parseInt(tokens[0]);
				score = Double.parseDouble(tokens[2]);
				chromosome = tokens[3];
				start = Integer.parseInt(tokens[4]);
				if (tokens.length>5) stop = Integer.parseInt(tokens[5]);
				else stop = start;
				regionsAL.add(new BindingRegion(rank, score,chromosome, start, stop, sizeNeighborhood));
			}
			regions = new BindingRegion[regionsAL.size()];
			regionsAL.toArray(regions);
			
		}catch (IOException e){
			e.printStackTrace();
		}
		return regions;
	}
	
	/**Attempts to fetch a serialized array of Interval[], then sorts/ ranks the intervals by the median
	 * ratio of the sub window.  It then uses it to build an array of BindingRegion.
	 * Will return null if it cannot fetch an Interval[].*/
	public static BindingRegion[] parseIntervalFile(File intervalFile, int sizeNeighborhood){
		//attempt to fetch intervals
		Object a = null;
		try {
			ObjectInputStream in =
				new ObjectInputStream(new FileInputStream(intervalFile));
			a = in.readObject();
			in.close();
		} catch (Exception e) {
			return null;
		}
		Interval[] intervals = (Interval[])a;
		//sort intervals by sub window
		int num = intervals.length;
		for (int i=0; i<num; i++){
			double score = 0;
			if (intervals[i].getBestSubWindow()!=null) score = intervals[i].getBestSubWindow().getMedianRatio();
			intervals[i].setSortBy(score);
		}
		Arrays.sort(intervals);
		//make binding regions
		BindingRegion[] regions = new BindingRegion[num];
		for (int i=0; i<num; i++){
			regions[i] = new BindingRegion(i+1, intervals[i].getSortBy(), intervals[i].getChromosome(), intervals[i].getStart1stOligo(), 
					intervals[i].getStartLastOligo()+ intervals[i].getSizeOfOligoMinusOne(), sizeNeighborhood);
		}
		return regions;
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Annotate Regions:     Jan 2005                          **\n" +
				"**************************************************************************************\n" +
				"Annotates a picks file finding surrounding protein coding genes.\n\n" +
				
				"-g Full path file text for the DmelRel4.0 GFF3 file.\n" +
				"-p Full path file text for the binding region picks or Interval file.\n" +
				"-c Full path file text for the Filtered CG names file (optional)\n"+
				"-b Size of neighborhood in bp, default is 10000 (optional)\n"+
				"-r Number of random trials (optional)\n"+
				"-n Just print number of neighbors for random trials (optional)\n"+
				"-s Skip filtering GFF file\n"+
				"\n" +
				"Example: java -jar pathTo/T2/Apps/AnnotateRegions -g /dmel/dmel_RELEASE4.0.gff3 -p\n" +
				"      /affy/zeste/finalPicks.txt -c /dmel/CGs.txt\n" +
				"\n" +
		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign new varibles*///stripGFF
	public void processArgs(String[] args){
		File directory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'g': gffFile = new File(args[i+1]); i++; break;
					case 'p': picksFile = new File(args[i+1]); i++; break;
					case 'c': dgcFile = new File(args[i+1]); i++; break;
					case 'r': numRandomTrials =Integer.parseInt(args[i+1]); i++; break;
					case 'b': sizeNeighborhood =Integer.parseInt(args[i+1]); i++; break;
					case 's': stripGFF=false; break;
					case 'n': printNumberOfNeighbors=true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		if (gffFile == null || picksFile == null){ //|| dgcFile == null){
			System.out.println("\nPlease enter the required files!\n");
			System.exit(0);
		}
	}
	/**Returns an array of GeneGroup whose names were found in the Hash*/
	public static GeneGroup[] filterGeneGroups(GeneGroup[] genes, HashSet cgNames){
		ArrayList pass = new ArrayList();
		int num = genes.length;
		for (int i=0; i<num; i++){
			if (cgNames.contains(genes[i].getName())) pass.add(genes[i]);
		}
		GeneGroup[] good = new GeneGroup[pass.size()];
		pass.toArray(good);
		return good;
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AnnotateRegions (args);
	}
}
