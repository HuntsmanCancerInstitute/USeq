package util.bio.parsers;
import java.io.*;
import util.bio.annotation.*;
import util.gen.*;
import java.util.*;


/**Parses a UCSC gene model table file (either refSeq or refFlat), use table options for #hg17.refGene.name	hg17.refLink.product	hg17.geneName.name	hg17.refSeqSummary.summary... example:
 *#refGene.name	refGene.chrom	refGene.strand	refGene.txStart	refGene.txEnd	refGene.cdsStart	refGene.cdsEnd	refGene.exonCount	refGene.exonStarts	refGene.exonEnds	refLink.product	geneName.name	refSeqSummary.summary
 *NM_198576	chr1	+	995569	1031415	995619	1030284	36	995569,997647,1010723,1016111,1016522,1016827,1017258,1018541,1018840,1019125,1019411,1019636,1020463,1020661,1021035,1021266,1021462,1021699,1022122,1022629,1022875,1023078,1023314,1024169,1024538,1024868,1025205,1025535,1025729,1026028,1026555,1026755,1027030,1029055,1029750,1030126,	995820,997909,1010771,1016327,1016747,1017052,1017465,1018760,1019035,1019326,1019560,1019742,1020580,1020826,1021179,1021391,1021568,1022038,1022260,1022757,1022990,1023198,1023668,1024362,1024754,1025098,1025340,1025632,1025894,1026140,1026672,1026948,1027118,1029280,1029854,1031415,	agrin	AGRIN	Agrin is a neuronal aggregating factor that induces the aggregation of ... * */
public class UCSCGeneModelTableReader {
	
	//fields
	private UCSCGeneLine[] geneLines = null;
	private HashMap<String,UCSCGeneLine[]> chromSpecificGeneLines = null;
	
	//test main
	public static void main(String[] args){
		//(new UCSCGeneModelTableReader (new File (args[0]), 1)).readOutHash();
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader (new File (args[0]), 0);
		System.out.println("Num lines "+reader.getGeneLines().length);
		//reader.removeOverlappingExons();
		/*
		UCSCGeneLine[] genes = reader.getGeneLines();
		//print introns
		for (int i=0; i< genes.length; i++){
			if (genes[i].getIntrons()!= null) {
				genes[i].setExons(null);
				System.out.println(genes[i].toStringAll());
			}
		}*/
	}
	
	//constructors
	/**@param numToSubtractFromEnd - UCSC uses interbase numbering, to get to stop inclusive numbering 
	 * you must subtract one from the ends of everything, thus set to 1.*/
	public UCSCGeneModelTableReader (File file, int numToSubtractFromEnd){
		parseGeneTableFile(file, numToSubtractFromEnd);
	}
	public UCSCGeneModelTableReader(){}
	
	
	//methods
	
	public Coordinate[] fetchGeneRegions(){
		Coordinate[] c = new Coordinate[geneLines.length];
		for (int i=0; i< geneLines.length; i++){
			c[i] = new Coordinate(geneLines[i].getChrom(), geneLines[i].getTxStart(), geneLines[i].getTxEnd());
		}
		return c;
	}
	
	/**Switches the exons for introns.  Only saves genes with introns.*/
	public void swapIntronsForExons(){
		ArrayList<UCSCGeneLine> genesWithIntrons = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLines.length; i++){
			ExonIntron[] introns = geneLines[i].getIntrons();
			if (introns != null && introns.length!=0) {
				geneLines[i].setExons(introns);
				genesWithIntrons.add(geneLines[i]);
			}
		}
		geneLines = new UCSCGeneLine[genesWithIntrons.size()];
		genesWithIntrons.toArray(geneLines);
	}
	
	public static byte[] makeBaseCountExonArray(UCSCGeneLine[] genes){
		//find max
		int max = 0;
		for (int i=0; i< genes.length; i++) {
			if (genes[i].getTxEnd() > max) max = genes[i].getTxEnd();
		}
		max++;
		//count hits leaving 0, no hits, 1, one hit, 2 multiple hits
		byte[] counts = new byte[max];
		//for each gene
		for (int i=0; i< genes.length; i++) {
			ExonIntron[] exons = genes[i].getExons();
			//for each exon
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int end = exons[j].getEnd();
				//for each base
				for (int k=start; k< end; k++){
					if (counts[k] !=2) counts[k]++;
				}
			}
		}
		return counts;
	}
	
	public String removeOverlappingExons(){
		StringBuilder eliminatedGenes = new StringBuilder();
		ArrayList<UCSCGeneLine> allGoodGenes = new ArrayList<UCSCGeneLine>();
		//split by chromosome
		splitByChromosome();
		
		//for each chromosome
		Iterator<String> it = chromSpecificGeneLines.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			System.out.print(".");
			UCSCGeneLine[] sub = chromSpecificGeneLines.get(chrom);
			
			//find intersecting exons
			HashMap<String,UCSCGeneLine> intersectingGenes = new HashMap<String,UCSCGeneLine>();
			for (int i=0; i< sub.length; i++){
				int start = sub[i].getTxStart();
				int stop = sub[i].getTxEnd();
				for (int j=i+1; j< sub.length; j++){
					if (sub[j].intersects(start, stop) && ExonIntron.intersect(sub[i].getExons(), sub[j].getExons())) {
						intersectingGenes.put(sub[i].getDisplayNameThenName(), sub[i]);
						intersectingGenes.put(sub[j].getDisplayNameThenName(), sub[j]);
					}
				}
			}
			
			//make array of igs
			int counter = 0;
			UCSCGeneLine[] igs = new UCSCGeneLine[intersectingGenes.size()];
			Iterator<String> names = intersectingGenes.keySet().iterator();
			while (names.hasNext()) igs[counter++] = intersectingGenes.get(names.next());
			
			//container for good genes
			ArrayList<UCSCGeneLine> goodGenes = new ArrayList<UCSCGeneLine>();
			
			//make exon counts array
			byte[] exonHits = makeBaseCountExonArray(igs);
			
			//for each intersecting gene
			for (int i=0; i< igs.length; i++){
				ExonIntron[] exons = igs[i].getExons();
				int[] minMax = ExonIntron.minMax(exons);
				boolean[] falseMask = new boolean[minMax[1]+1];
				Arrays.fill(falseMask, true);
				int numInt = 0;
				int exonBps = 0;
				//for each exon
				for (int j=0; j< exons.length; j++){
					int start = exons[j].getStart();
					int end = exons[j].getEnd();
					exonBps+= exons[j].getLength();
					//for each base, if it's a 1 then flip to false, if it's a 2 the leave true
					for (int k=start; k< end; k++){
						if (exonHits[k] !=2) falseMask[k-minMax[0]] = false;
						else numInt++;
					}
				}
				//retrieve false blocks for A, ends included
				int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(falseMask, 0, 0);
				//convert to interbase coordinate exon introls
				ExonIntron[] ei = new ExonIntron[blocks.length];
				int minPlusOne = minMax[0] +1;
				for (int j=0; j< ei.length; j++){
					ei[j] = new ExonIntron(blocks[j][0]+minMax[0], blocks[j][1]+minPlusOne);
				}
				
				//set exons
				igs[i].setExons(ei);
				//reset boolean mask to all true
				for (int j=0; j< exons.length; j++){
					int start = exons[j].getStart();
					int end = exons[j].getEnd();
					for (int k=start; k< end; k++) falseMask[k] = true;
				}
				//save?
				if (ei != null && ei.length !=0) {
					goodGenes.add(igs[i]);
					double diff = (double)numInt/(double)exonBps;
					if (diff >= 0.5) eliminatedGenes.append (igs[i].getDisplayNameThenName() +"("+Num.formatPercentOneFraction(diff)+"), ");
					//System.out.println(igs[i]);
				}
				else eliminatedGenes.append (igs[i].getDisplayNameThenName() +"(100%), ");
				//sb.append(numInt+"\t"+igs[i].getName()+"\n");
			}
			
			//add non intersecting genes
			for (int i=0; i< sub.length; i++){
				if (intersectingGenes.containsKey(sub[i].getDisplayNameThenName()) == false) {
					goodGenes.add(sub[i]);
					//System.out.println(sub[i]);
				}
			}
			
			//set in hash
			sub = new UCSCGeneLine[goodGenes.size()];
			goodGenes.toArray(sub);
			Arrays.sort(sub, new UCSCGeneLineChromComparator());
			chromSpecificGeneLines.put(chrom, sub);
			//add all
			allGoodGenes.addAll(goodGenes);
		}
		System.out.println();
		//set geneLines
		geneLines = new UCSCGeneLine[allGoodGenes.size()];
		allGoodGenes.toArray(geneLines);
		return eliminatedGenes.toString();
	}
	
	/**Prints all exons to standard out.  Subtracts bpEndMod from start, adds it to stop.*/
	public void printAllExons (int bpEndMod){
		Iterator it = chromSpecificGeneLines.keySet().iterator();
		while (it.hasNext()){
			String chrom = (String) it.next();
			UCSCGeneLine[] sub = (UCSCGeneLine[]) chromSpecificGeneLines.get(chrom);
			for (int i=0; i< sub.length; i++){
				ExonIntron[] exons = sub[i].getExons();
				for (int j=0; j< exons.length; j++){
					int start = exons[j].getStart() - bpEndMod;
					if (start < 0) start = 0;
					int end = exons[j].getEnd()+ bpEndMod;
					System.out.println(chrom+"\t"+start+"\t"+end);
				}
			}
			
		}
	}
	
	public Coordinate[] fetchExons() {
		ArrayList<Coordinate> cAL = new ArrayList<Coordinate>();
		for (int i=0; i< geneLines.length; i++){
			ExonIntron[] exons = geneLines[i].getExons();
			for (int j=0; j< exons.length; j++){
				cAL.add(new Coordinate(geneLines[i].getChrom(), exons[j].getStart(), exons[j].getEnd()));
			}
		}
		Coordinate[] c = new Coordinate[cAL.size()];
		cAL.toArray(c);
		return c;
	}
	
	public void writeGeneTableToFile (File file){
		try{
			PrintWriter out = new PrintWriter(new FileWriter(file));
			for (int i=0; i< geneLines.length; i++) out.println(geneLines[i]);
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**To test your hash.*/
	public void readOutHash(){
		Iterator it = chromSpecificGeneLines.keySet().iterator();
		while (it.hasNext()){
			String chrom = (String) it.next();
			System.out.println(chrom);
			UCSCGeneLine[] sub = (UCSCGeneLine[]) chromSpecificGeneLines.get(chrom);
			System.out.println(sub.length + " "+ sub[0]);
		}
	}
	
	/**Splits a UCSCGeneLine[] by chromosome into a HashMap of chromosome:UCSCGeneLine[].*/
	public void splitByChromosome(){
		if (chromSpecificGeneLines != null || geneLines == null || geneLines.length == 0) return;
		Arrays.sort(geneLines, new UCSCGeneLineChromComparator());
		chromSpecificGeneLines = new HashMap();
		ArrayList al = new ArrayList();
		String currChrom = geneLines[0].getChrom();
		for (int i=0; i< geneLines.length; i++){
			if (geneLines[i].getChrom().equals(currChrom) == false){
				UCSCGeneLine[] sub = new UCSCGeneLine[al.size()];
				al.toArray(sub);
				chromSpecificGeneLines.put(currChrom, sub);
				al.clear();
				currChrom = geneLines[i].getChrom();
			}
			al.add(geneLines[i]);
		}
		//add last to hash
		UCSCGeneLine[] sub = new UCSCGeneLine[al.size()];
		al.toArray(sub);
		chromSpecificGeneLines.put(currChrom, sub);
		al.clear();
	}
	
	/**Splits a UCSCGeneLine[] by chromosome and strand into a HashMap of chromosomeStrand:UCSCGeneLine[].*/
	public HashMap<String,UCSCGeneLine[]> splitByChromosomeAndStrand(){
		//sort by chromStrand, position, short to long
		Arrays.sort(geneLines, new UCSCGeneLineChromStrandComparator());		
		HashMap<String,UCSCGeneLine[]> lines = new HashMap<String,UCSCGeneLine[]>();
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		String currChromStrand = geneLines[0].getChrom()+geneLines[0].getStrand();
		for (int i=0; i< geneLines.length; i++){
			String testChromStrand = geneLines[i].getChrom()+geneLines[i].getStrand();
			if (testChromStrand.equals(currChromStrand) == false){
				UCSCGeneLine[] sub = new UCSCGeneLine[al.size()];
				al.toArray(sub);
				lines.put(currChromStrand, sub);
				al.clear();
				currChromStrand = testChromStrand;
			}
			al.add(geneLines[i]);
		}
		//add last to hash
		UCSCGeneLine[] sub = new UCSCGeneLine[al.size()];
		al.toArray(sub);
		lines.put(currChromStrand, sub);
		al.clear();
		return lines;
	}
	
	/**Make a hash of display name (geneName).*/
	public HashMap<String,UCSCGeneLine> fetchGeneNameHash(){
		HashMap<String,UCSCGeneLine> lines = new HashMap<String,UCSCGeneLine>();
		for (int i=0; i< geneLines.length; i++){
			lines.put(geneLines[i].getDisplayName(), geneLines[i]);
		}
		return lines;
	}

	/**@param numToSubtractFromEnd - UCSC uses interbase numbering, to get to stop inclusive numbering 
	 * you must subtract one from the ends of everything, thus set to 1.*/
	public void parseGeneTableFile(File file, int numToSubtractFromEnd){
		BufferedReader in;
		try{
			in = IO.fetchBufferedReader(file);
			String line;
			ArrayList al = new ArrayList();
			boolean displayNamePresent = false;
			boolean checkForDisplayName = true;
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0 || line.startsWith("#")) continue;
				String[] tokens = line.split("\\t");
				//check for a display text?
				if (checkForDisplayName){
					//look for +,-,or . in tokens[3]
					if (tokens[3].matches("[+-\\.]")) {
						displayNamePresent = true;						
					}
					checkForDisplayName = false;
				}
				if (tokens.length < 10) Misc.printExit("\nError: problem parsing gene table file with line -> "+line);
				al.add(new UCSCGeneLine(tokens, numToSubtractFromEnd,displayNamePresent));
			}
			in.close();
			geneLines = new UCSCGeneLine[al.size()];
			al.toArray(geneLines);
		
		} catch (IOException e){
			e.printStackTrace();
		}
	}

	public HashMap<String,UCSCGeneLine[]> getChromSpecificGeneLines() {
		if (chromSpecificGeneLines == null) splitByChromosome();
		return chromSpecificGeneLines;
	}

	public UCSCGeneLine[] getGeneLines() {
		return geneLines;
	}
	
	public boolean checkStartStopOrder(){
		for (int i=0; i<geneLines.length; i++){
			if (geneLines[i].getTxEnd() < geneLines[i].getTxStart()) {
				System.err.println("Bad Line -> "+geneLines[i]);
				return false;
			}
		}
		return true;
	}

	public void setGeneLines(UCSCGeneLine[] geneLines) {
		this.geneLines = geneLines;
	}

	/**Replaces the sets the Name to what's in the DisplayName (1st column when both are present). Very confusing!.*/
	public void replaceNameWithDisplayName() {
		if (geneLines == null || geneLines[0].getDisplayName() == null) return;
		for (int i=0; i< geneLines.length; i++) geneLines[i].setName(geneLines[i].getDisplayName());
	}

	/**Check to see if gene names are unique.*/
	public boolean uniqueGeneNames() {
		HashSet<String> uni = new HashSet<String>();
		for (UCSCGeneLine l: geneLines) {
			if (uni.contains(l.getDisplayNameThenName())){
				System.err.println("Bad Line -> "+l);
				return false;
			}
			uni.add(l.getDisplayNameThenName());
		}
		return true; 
	}

	/**Checks to see if all of the genes have a strand designation of + or -*/
	public boolean checkStrand() {
		for (UCSCGeneLine l: geneLines) {
			if (l.getStrand().equals("+") == false && l.getStrand().equals("-") == false) {
				System.err.println("Bad Line -> "+ l);
				return false;
			}
		}
		return true;
	}
	
}
