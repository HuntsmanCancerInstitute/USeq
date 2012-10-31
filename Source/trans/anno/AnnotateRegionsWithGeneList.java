package trans.anno;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import trans.main.*;
import util.bio.annotation.*;
import util.bio.parsers.gff.*;
import util.gen.*;


/**
 * For annotating a picks file (chr, start, stop) with gene information for a particular list of CG names.
 * Specific for dmel release 4.0 gff3.  
 * 
 * @see AnnotateRegions
 */
public class AnnotateRegionsWithGeneList {
	
	//fields
	private File gffFile;
	private File picksFile;
	private File hotGeneListFile;
	private String[] hotGeneList;
	private GeneGroup[] geneGrps;
	private BindingRegion[] bindingRegions;
	private int numRandomTrials = 0;
	private int sizeNeighborhood = 1000;
	private boolean printNumberOfNeighbors = false;
	private boolean printSummaryStats = false;
	private boolean useHotGeneList = true;
	
	public AnnotateRegionsWithGeneList(String[] args){
		processArgs(args);
		System.out.println("\nLaunching...");
		
		//process gff file retrieving an array of GeneGroup
		System.out.println("\tProcessing GFF file...");
		DmelRel4Extractor ex = new DmelRel4Extractor();
		ex.extract(gffFile, true);
		ArrayList geneGrpsAL = ex.getGeneGroupArrayList();
		int num = geneGrpsAL.size();
		
		if (useHotGeneList){
			System.out.println("\tProcessing Hot Gene List...");
			//convert to a hashMap
			HashMap map = new HashMap(num);
			GeneGroup g;
			for (int i=0; i<num; i++){
				g = (GeneGroup)geneGrpsAL.get(i);
				map.put(g.getName(), g);
			}
			//pull out those in hotGeneListFile
			hotGeneList = IO.loadFileIntoStringArray(hotGeneListFile);
			int numHot = hotGeneList.length;
			ArrayList hotGeneGrps = new ArrayList(numHot);
			for (int i=0; i<numHot; i++){
				if (map.containsKey(hotGeneList[i])) hotGeneGrps.add(map.get(hotGeneList[i]));
				else System.out.println("\tWarning, couldn't find '"+hotGeneList[i]+"' in gff derived gene groups.");
			}
			num = hotGeneGrps.size();
			geneGrps = new GeneGroup[num];
			hotGeneGrps.toArray(geneGrps);
		}
		else{
			geneGrps = new GeneGroup[num];
			geneGrpsAL.toArray(geneGrps);
		}
		
		//process picks file
		System.out.println("\tProcessing picks file...");
		bindingRegions = parseIntervalFile(picksFile, sizeNeighborhood);
		if (bindingRegions == null) bindingRegions = parsePicksFile(picksFile, sizeNeighborhood);
		
		//cross to compare
		System.out.println("\tComparing Binding Regions against Gene Groups...");
		compareBindingRegionsVsGeneGrps(geneGrps, bindingRegions);
		
		//sort
		Arrays.sort(bindingRegions);
		
		//print results
		int numBindingRegions = bindingRegions.length;
		if (printSummaryStats == false ){
			System.out.println("\nBinding Regions:");
			
			for (int i=0; i<numBindingRegions; i++){
				System.out.println(bindingRegions[i]);
				System.out.println("*****************************");
			}
			System.out.println();
			//printDistToClosestATGAndTranscript(bindingRegions);
			
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
		}
		System.out.println("\n"+numBindingRegions+" Number of binding regions.");
		System.out.println(geneGrps.length+" Number of gene groups.");
		System.out.println(sizeNeighborhood+" Bp size of neighborhood.\n");
		
		
		//count flanking genes
		int[] sides = countGenes (bindingRegions);
		System.out.println(sides[0]+" Number of genes with a binding region within the neighborhood.");
		System.out.println(sides[1]+" Number of genes with a binding region on their 5' stop.");
		System.out.println(sides[2]+" Number of genes with a binding region on their 3' stop.");
		System.out.println(sides[3]+" Number of genes with a binding region overlapping their 5' stop.");
		System.out.println(sides[4]+" Number of genes with a binding region overlapping their 3' stop.");
		System.out.println(sides[5]+" Number of genes entirely containing a binding region.");
		System.out.println(sides[6]+" Number of binding regions with neighbors within the neighborhood.");
		System.out.println(sides[7]+" Number of binding regions with no neighbors within the neighborhood.");
		System.out.println(sides[8]+" Average size of a binding region.");
		System.out.println(sides[9]+" Number of binding regions in non coding DNA.");
		System.out.println(sides[10]+" Number of binding regions in coding DNA.");
		System.out.println(sides[11]+" Number of binding regions that overlap both coding and non coding DNA.\n");
		
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
			int numColumns = sides.length;
			double[] totals = new double[numColumns]; //total of all for calculating average
			int[] numLess = new int[numColumns]; //num less for p-value
			int[] numGreater = new int[numColumns]; //num more for p-value
			Arrays.fill(totals, 0);
			Arrays.fill(numLess, 0);
			Arrays.fill(numGreater, 0);
			
			for (int i=0; i<numRandomTrials; i++){
				random = makeRandomBindingRegions(bindingRegions, chromLengths, sizeNeighborhood);
				compareBindingRegionsVsGeneGrps(geneGrps, random);
				if (printNumberOfNeighbors){
					System.out.println(countNumberBindingRegionsWithNeighbors(random));
				}
				else{
					counts = countGenes (random);
					for (int j=0; j<numColumns; j++){
						System.out.print("\t"+counts[j]);
						totals[j] += counts[j];
						if (counts[j]>= sides[j]) numGreater[j]++;
						if (counts[j]<= sides[j]) numLess[j]++;
					}
					System.out.println();
					
				}
			}
			
			//print summary lines
			System.out.print("\n\t# Genes w/ region +/- nbrhd \t# Genes w/ region on 5' stop \t#Genes w/ region on 3' stop \t");
			System.out.print("# Genes w/ region overlapping 5' stop \t# Genes w/ region overlapping 3' stop \t# Genes entirely containing a region\t");
			System.out.print("# Regions w/ >=1 genes within nbrhd \t# Regions w/ no genes within nbrhd \tMean size regions \t");
			System.out.print("# Regions in non-coding DNA \t# Regions in coding DNA \t# Regions the overlap coding and non-coding\n");
			//print actual observations
			System.out.print("\nActual\t");
			for (int j=0; j<numColumns; j++){
				System.out.print(sides[j]+"\t");
			}
			//averages
			System.out.print("\nMean Random\t");
			for (int j=0; j<numColumns; j++){
				System.out.print(totals[j]/(double)numRandomTrials+"\t");
			}
			System.out.print("\nRatio\t");
			for (int j=0; j<numColumns; j++){
				if (sides[j]!=0 && totals[j]!=0) System.out.print(Num.formatNumber(((double)sides[j]/(totals[j]/(double)numRandomTrials)),2)+"\t");
				else System.out.print("0\t");
			}
			System.out.println("\n");
			//num less
			System.out.print("Num <= Actual\t");
			for (int j=0; j<numColumns; j++){
				System.out.print(numLess[j]+"/"+numRandomTrials+"\t");
			}
			System.out.println();
			System.out.print("P-value\t");
			for (int j=0; j<numColumns; j++){
				System.out.print(numLess[j]/(double)numRandomTrials+"\t");
			}
			System.out.println();
			
			System.out.println("\n");		
			//num less
			System.out.print("Num >= Actual\t");
			for (int j=0; j<numColumns; j++){
				System.out.print(numGreater[j]+"/"+numRandomTrials+"\t");
			}
			System.out.println();
			System.out.print("P-value\t");
			for (int j=0; j<numColumns; j++){
				System.out.print(numGreater[j]/(double)numRandomTrials+"\t");
			}
			System.out.println();
			System.out.println("\n");
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
		if (br.getGeneGrp5Prime() != null){
			test = findDistanceToATG(br, br.getGeneGrp5Prime());
			if (test == 0) return 0;
			if (test<dist) dist = test;
		}
		//3'
		if (br.getGeneGrp3Prime() != null){
			test = findDistanceToATG(br, br.getGeneGrp3Prime());
			if (test == 0) return 0;
			if (test<dist) dist = test;
		}
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
	
	/**Returns:
	 * the number of genes with one or more binding regions within the neighborhood.
	 * the number of genes where the binding region is on the 5' stop and
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
		int numGenesWithBR = 0;
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
		
		//calc number of genes with a binding region within neighborhood
		HashSet neighbors = new HashSet(num * 4);
		for (int i=0; i<num; i++){
			neighbors.addAll(extractCGNames(br[i].getNeighboringGeneGrps()));
		}
		numGenesWithBR = neighbors.size();
		//run thru binding regions
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
		return new int[]{numGenesWithBR, numFivePrime, numThreePrime, numFiveOverlap, numThreeOverlap, numContained, numBRsWithNeighbors,
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
			int rank = 1;
			String chromosome;
			int start;
			int stop;
			double score;
			//chrom, start, stop, score
			//  0      1      2    3      
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.equals("")) continue;
				tokens = line.split("\\s+");
				if (tokens.length<3) continue;
				chromosome = tokens[0];
				start = Integer.parseInt(tokens[1]);
				if (tokens.length == 3){
					stop = start;
					score = Double.parseDouble(tokens[2]);
				}
				else{
					stop = Integer.parseInt(tokens[2]);
					score = Double.parseDouble(tokens[3]);
				}
				regionsAL.add(new BindingRegion(rank++, score,chromosome, start, stop, sizeNeighborhood));
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
				"**                Gene List Enrichment:     June 2005                               **\n" +
				"**************************************************************************************\n" +
				"Annotates a picks file with a hot gene list file.\n\n" +
				
				"-g Full path file text for the DmelRel4.1 Stripped GFF3 file.\n" +
				"-p Full path file text for the binding region picks or Interval file.\n" +
				"-c Full path file text for the hot CG names file (optional)\n"+
				"-b Size of neighborhood in bp, default is 1000 (optional)\n"+
				"-r Number of random trials (optional)\n"+
				"-n Just print number of neighbors for random trials (optional)\n"+
				"-s Just print the summary statistics (optional)\n\n"+
				
		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
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
					case 'c': hotGeneListFile = new File(args[i+1]); i++; break;
					case 'r': numRandomTrials =Integer.parseInt(args[i+1]); i++; break;
					case 'b': sizeNeighborhood =Integer.parseInt(args[i+1]); i++; break;
					case 'n': printNumberOfNeighbors=true; break;
					case 's': printSummaryStats=true; break;
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
		
		if (hotGeneListFile == null) useHotGeneList = false;
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AnnotateRegionsWithGeneList (args);
	}
}
