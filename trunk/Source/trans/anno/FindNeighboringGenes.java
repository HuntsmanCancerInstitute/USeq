package trans.anno;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.annotation.*;
import util.bio.parsers.*;
import util.gen.*;


/**
 * Finds closest genes and all genes within a given neighborhood.
 * Needs a UCSC Gene table and a pick list of regions.
 * Distance to a gene is based on the middle of a region to the gene's tss.
 * This is a simplified version of AnnotateRegions.
 */
public class FindNeighboringGenes {

	//fields
	private File ucscGeneTableFile;
	private File workingPicksFile;
	private File[] pickFiles;
	private int sizeNeighborhood = 10000;
	private HashMap chrGeneLines;
	private HashMap chrPicks;
	private CoorGene[][] coorGenes;
	private CoorGene[] allCoorGenes;
	private boolean printAll = true;
	private boolean printOneLine = false;
	private boolean findOverlapping = false; 

	public FindNeighboringGenes(String[] args){
		processArgs(args);
		System.out.println("\nLaunching...");

		//parse ucsc table
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader (ucscGeneTableFile, 1);
		chrGeneLines = reader.getChromSpecificGeneLines();
		System.out.println("\t"+reader.getGeneLines().length+" Gene Models\n");

		//process each picks file
		for (int i=0; i< pickFiles.length; i++){
			workingPicksFile = pickFiles[i];
			System.out.println("\t"+pickFiles[i].getName());
			ScoredCoordinate[] sCorr = ScoredCoordinate.makeRanked(Coordinate.parseFile(workingPicksFile, 0 , 0));
			Arrays.sort(sCorr);
			chrPicks = ScoredCoordinate.splitByChromosome(sCorr);

			findNeighbors();

			printCoorGenes(allCoorGenes);
		}
		//Misc.printArray(fetchGeneNames());

		System.out.println("\nDone!\n");
	}

	public FindNeighboringGenes(File ucscGeneTableFile, Coordinate[] picks, int sizeNeighborhood, boolean findOverlapping){
		this.ucscGeneTableFile = ucscGeneTableFile;
		this.sizeNeighborhood = sizeNeighborhood;
		this.findOverlapping = findOverlapping;
		//make picks
		ScoredCoordinate[] sCorr = ScoredCoordinate.makeRanked(picks);
		Arrays.sort(sCorr);
		chrPicks = ScoredCoordinate.splitByChromosome(sCorr);
		//parse ucsc table
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader (ucscGeneTableFile, 1);
		chrGeneLines = reader.getChromSpecificGeneLines();
		//find neighbors
		findNeighbors();
	}

	public String[] fetchGeneNames(){
		String[] associatedGenes = new String[allCoorGenes.length];
		for (int i=0; i< allCoorGenes.length; i++){
			ArrayList al = allCoorGenes[i].genes;
			DistanceGene[] dg = new DistanceGene[al.size()];
			al.toArray(dg);
			//make a hash to remove dups
			HashSet names = new HashSet();
			for (int j =0; j< dg.length; j++) names.add(dg[j].gene.getName());
			String x = names.toString();
			associatedGenes[i] = x.substring(1, x.length()-1);
		}
		return associatedGenes;
	}

	public void findNeighbors(){

		//intersect ScoredCoordinate with GeneLines
		intersectRegionsWithAnnotation();
		//merge CoorGenes
		int num = Misc.totalLength(coorGenes);
		allCoorGenes = new CoorGene[num];
		int counter = 0;
		for (int i=0; i< coorGenes.length; i++) {
			for (int j=0; j< coorGenes[i].length; j++){
				allCoorGenes[counter++] = coorGenes[i][j];
			}
		}
		//sort by original rank
		CoorGeneComparator comp = new CoorGeneComparator();
		Arrays.sort(allCoorGenes, comp);
	}

	public void printCoorGenes(CoorGene[] all){
		try {
			String trunk = Misc.removeExtension(workingPicksFile.getName());
			PrintWriter out = new PrintWriter( new FileWriter(new File(workingPicksFile.getParent(), trunk+"_FNG.xls")));
			//header
			out.println("Printed below are the region coordinates followed by the closest genes and any genes within " +
			"the neighborhood (Distance to TSS from middle of region, Gene text, Chr, Start, Stop, Strand, Position TSS)\n");
			//for each chromosome
			for (int i=0; i< all.length; i++) {
				//print coordinates
				if (printOneLine) out.print(all[i].coordinates);
				else out.println(all[i].coordinates);
				//make a hash to remove dups
				HashSet names = new HashSet();
				//print closest
				int sizeClosestGenes = all[i].closestGeneTSS.size();
				if (sizeClosestGenes !=0) {
					for (int k=0; k<sizeClosestGenes; k++){
						DistanceGene x = (DistanceGene) all[i].closestGeneTSS.get(k);
						if (printOneLine) out.print("\t"+x.distance+"\t"+x.gene.simpleToString()+"\tClosest");
						else out.println("\t\t\t"+x.distance+"\t"+x.gene.simpleToString()+"\tClosest");
						names.add(x.gene.getName());
					}
				}
				//print others within neighborhood sorted by distance?
				if (printAll){
					ArrayList al = all[i].genes;
					DistanceGene[] dg = new DistanceGene[al.size()];
					al.toArray(dg);
					Arrays.sort(dg);
					for (int k=0; k< dg.length; k++){
						if (names.contains(dg[k].gene.getName())) continue;
						if (printOneLine) out.print("\t"+dg[k].distance+ "\t"+ dg[k].gene.simpleToString());
						else out.println("\t\t\t"+dg[k].distance+ "\t"+ dg[k].gene.simpleToString());
					}
				}
				//close?
				if (printOneLine) out.println();
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void intersectRegionsWithAnnotation(){

		//for each chromosome
		Iterator it = chrPicks.keySet().iterator();
		coorGenes = new CoorGene[chrPicks.keySet().size()][];
		int index = 0;
		while (it.hasNext()){
			String chrom = (String)it.next();
			ScoredCoordinate[] coor = (ScoredCoordinate[]) chrPicks.get(chrom);
			UCSCGeneLine[] genes = (UCSCGeneLine[]) chrGeneLines.get(chrom);
			if (findOverlapping) {
				coorGenes[index++] = overlap (coor, genes);

			}
			else coorGenes[index++] = intersect (coor, genes);
		}
	}

	private CoorGene[] intersect( ScoredCoordinate[] coor, UCSCGeneLine[] genes){
		CoorGene[] coorGenes = new CoorGene[coor.length];
		//for each coor
		for (int i=0; i< coor.length; i++){
			coorGenes[i] = new CoorGene(coor[i]);
			//any genes?
			if (genes == null) continue;
			//for each gene
			int shortestDistance = 1000000000;
			for (int j=0; j< genes.length; j++){
				int distance = Math.abs(coorGenes[i].middle - genes[j].getTss());
				//add to neighborhood?
				DistanceGene dg = null;
				if (distance <= sizeNeighborhood) {
					dg = new DistanceGene (distance, genes[j]);
					coorGenes[i].genes.add(dg);
				}
				//closest gene tss?
				if (distance < shortestDistance){
					shortestDistance = distance;
					if (dg == null){ 
						dg = new DistanceGene (distance, genes[j]); 
					}
					coorGenes[i].closestGeneTSS = new ArrayList();
					coorGenes[i].closestGeneTSS.add(dg);
				}
				else if (distance == shortestDistance){
					if (dg == null){ 
						dg = new DistanceGene (distance, genes[j]); 
					}
					coorGenes[i].closestGeneTSS.add(dg);
				}
			}
		}
		return coorGenes;
	}

	private CoorGene[] overlap( ScoredCoordinate[] coor, UCSCGeneLine[] genes){		
		CoorGene[] coorGenes = new CoorGene[coor.length];
		//for each coor
		for (int i=0; i< coor.length; i++){
			coorGenes[i] = new CoorGene(coor[i]);
			//any genes?
			if (genes == null) continue;
			//for each gene			
			for (int j=0; j< genes.length; j++){
				//do they intersect?
				boolean intersects = coorGenes[i].intersects(genes[j]);
				int distance = Math.abs(coorGenes[i].middle - genes[j].getTss());
				//add to neighborhood?
				DistanceGene dg = null;
				if (intersects) {
					dg = new DistanceGene (distance, genes[j]);
					coorGenes[i].genes.add(dg);
				}

			}
		}
		return coorGenes;
	}

	private class CoorGene {
		private ScoredCoordinate coordinates;
		private int middle;
		int expStart;
		int expStop;
		private ArrayList closestGeneTSS = new ArrayList();
		private  ArrayList genes = new ArrayList();

		public CoorGene (ScoredCoordinate coordinates){
			this.coordinates = coordinates;
			middle = ((coordinates.getStop() - coordinates.getStart())/2) + coordinates.getStart();
			expStart = coordinates.getStart()- sizeNeighborhood;
			expStop = coordinates.getStop() + sizeNeighborhood;
		}

		public boolean intersects (UCSCGeneLine gene) {
			if (gene.getTxEnd() <= expStart || gene.getTxStart()>= expStop) return false;
			return true;
		}
	}

	private class CoorGeneComparator  implements Comparator {
		/**Sorts by ScoredCoordinate*/
		public int compare(Object arg0, Object arg1) {
			ScoredCoordinate first = ((CoorGene)arg0).coordinates;
			ScoredCoordinate second = ((CoorGene)arg1).coordinates;
			if (first.getScore()> second.getScore()) return 1;
			if (first.getScore()< second.getScore()) return -1;
			return 0;
		}

	}


	private class DistanceGene implements Comparable{
		private int distance;
		private UCSCGeneLine gene;
		public DistanceGene (int distance, UCSCGeneLine gene){
			this.distance = distance;
			this.gene = gene;
		}
		public int compareTo(Object obj){
			DistanceGene other = (DistanceGene)obj;
			if (other.distance > distance) return -1;
			if (other.distance < distance) return 1;
			return 0;
		}
	}



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Find Neighboring Genes:   Nov 2008                      **\n" +
				"**************************************************************************************\n" +
				"FNG takes a list of genes in UCSC Gene Table format and intersects them with a list of\n" +
				"regions finding the closest gene to each region as well as all of the genes that fall\n" +
				"within a given neighborhood. Distance is measured from the center of the region to the\n" +
				"transcription start site/ 1st base position in 1st exon. See Tables link under\n" +
				"http://genome.ucsc.edu/ . Note, output coordinates are zero based, stop inclusive.\n\n" +

				"-g Full path file text for a tab delimited UCSC Gene Table (text chrom strand txStart\n" +
				"      txEnd cdsStart cdsEnd exonCount exonStarts exonEnds etc...) .\n" +
				"-p Full path file/directory text for tab delimited region list(s) (chr, start, stop) .\n" +
				"-b Size of neighborhood in bp, default is 10000 \n"+
				"-f Find genes that overlap neighborhood irregardles of distance to TSS.\n"+
				"-c Only print closest genes.\n"+
				"-o Print neighbors on one line.\n"+

				"\n" +
				"Example: java -jar pathTo/T2/Apps/FindNeighboringGenes -g /anno/hg17Ensembl.txt -p\n" +
				"      /affy/p53/finalPicks.txt -b 5000 -c\n" +
				"\n" +
		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign new varibles*///stripGFF
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'g': ucscGeneTableFile = new File(args[i+1]); i++; break;
					case 'p': pickFiles = IO.extractFiles(new File(args[++i])); break;
					case 'b': sizeNeighborhood =Integer.parseInt(args[i+1]); i++; break;
					case 'c': printAll = false; break;
					case 'o': printOneLine = true; break;
					case 'f': findOverlapping = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		if (ucscGeneTableFile == null || pickFiles[0] == null){ 
			Misc.printExit("\nPlease enter the required files!\n");
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FindNeighboringGenes (args);
	}

	public boolean isFindOverlapping() {
		return findOverlapping;
	}

	public void setFindOverlapping(boolean findOverlapping) {
		this.findOverlapping = findOverlapping;
	}
}

