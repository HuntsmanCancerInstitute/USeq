package util.bio.annotation;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import util.gen.*;
import util.bio.parsers.*;

public class AnnotateBedWithGenes {

	private File ucscTableFile;	
	private File bedFile;
	private File resultsFile;
	private HashMap<String,UCSCGeneLine[]> chrGenes;
	private HashMap<String,Bed[]> chrBed;
	private int bpPad = 0;
	private Gzipper out = null;
	private boolean intersectExons = true;
	private int[] coorIndexes = {0,1,2};
	boolean verbose = true;
	boolean addGeneFeature = false;
	private int numIntersectingRegions = 0;
	private int totalNumberRegions = 0;
	private int numIntersecting5UTRs = 0;
	private int numIntersecting3UTRs = 0;
	private int numIntersectingExons = 0;
	private int numIntersectingIntrons = 0;

	public AnnotateBedWithGenes (String[] args){
		try {
			//Parse args
			processArgs(args);

			//parse genes
			parseFiles();

			if (totalNumberRegions != 0) {
				//for each bed chrom
				intersect();
				//sumary stats
				if (verbose)  {
					IO.pl("\n\nIntersection stats:");
					IO.pl("\t"+numIntersectingRegions+"\tNum intersecting regions");
					IO.pl("\t"+(totalNumberRegions-numIntersectingRegions)+"\tNum non intersecting regions");
					if (addGeneFeature) {
						IO.pl("\t"+numIntersecting5UTRs+"\tNum intersecting 5pUTRs");
						IO.pl("\t"+numIntersecting3UTRs+"\tNum intersecting 3pUTRs");
						IO.pl("\t"+numIntersectingExons+"\tNum intersecting Exons");
						IO.pl("\t"+numIntersectingIntrons+"\tNum intersecting Introns");
					}
				}
			}

			out.close();
		} catch (Exception e) {
			System.err.println("\nError annotating "+bedFile.getName());
			if (out != null) out.getGzipFile().delete();
			e.printStackTrace();
		}
	}

	/**Adds genes that intersect onto the Bed.name as name+"\t"+gene1,gene2...*/
	public AnnotateBedWithGenes(File ucscTableFile, Bed[] bed, boolean intersectExons){
		try {
			verbose = false;
			this.intersectExons = intersectExons;

			UCSCGeneModelTableReader tr = new UCSCGeneModelTableReader(ucscTableFile, 0);
			chrGenes = tr.getChromSpecificGeneLines();

			Arrays.sort(bed);
			chrBed = Bed.splitBedByChrom(bed);

			intersect();

		} catch (Exception e) {
			System.err.println("\nError annotating regions!");
			e.printStackTrace();
		}
	}

	private void intersect() throws IOException {
		HashSet<String> toAdd = new HashSet<String>();
		if (verbose)  System.out.print("\n\t");

		//for each chrom of bed regions
		for (String chr: chrBed.keySet()){
			if (verbose)  System.out.print(chr+" ");
			Bed[] bed = chrBed.get(chr);
			UCSCGeneLine[] genes = chrGenes.get(chr);
			if (genes == null) {
				if (verbose) System.err.print("\nWARNING: No genes for chrom "+chr+", skipping.\n\t");
				for (Bed region: bed) {
					region.setName(region.getName()+"\t.");
					if (out != null) out.println(region.getName());
				}
				continue;
			}

			//build interval tree of exons from the genes
			IntervalST<ArrayList<UCSCGeneLine>> exonTree = buildTree(genes);

			//for each bed region
			for (Bed region: bed){
				int start = region.getStart();
				int stop = region.getStop();
				if (start < 0) start = 0;
				//end is included in search so subtract 1
				Iterable<Interval1D> it = exonTree.searchAll(new Interval1D(start, stop-1));
				toAdd.clear();
				for (Interval1D x : it) {
					ArrayList<UCSCGeneLine> g = exonTree.get(x);
					//for each gen line extract wanted info and collapse with the Hash
					for (UCSCGeneLine gl: g){
						if (addGeneFeature) {
							String feature = gl.getIntGeneFeature(start, stop);
							if (feature != null) {
								toAdd.add(gl.getNamesCollapsed(":")+"_"+feature);
								if (feature.equals("Exon")) numIntersectingExons++;
								else if (feature.equals("Intron")) numIntersectingIntrons++;
								else if (feature.equals("5pUTR")) numIntersecting5UTRs++;
								else if (feature.equals("3pUTR")) numIntersecting3UTRs++;
							}
							else toAdd.add(gl.getNamesCollapsed(":"));
						}
						else toAdd.add(gl.getNamesCollapsed(":"));
					}
				}
				//print it
				String anno = ".";
				if (toAdd.size() !=0) {
					anno = Misc.hashSetToString(toAdd, ",");
					numIntersectingRegions++;
				}
				region.setName(region.getName()+"\t"+anno);
				if (out != null) out.println(region.getName());
			}	
		}
	}

	private IntervalST<ArrayList<UCSCGeneLine>> buildTree(UCSCGeneLine[] genes){
		IntervalST<ArrayList<UCSCGeneLine>> exonTree = new IntervalST<ArrayList<UCSCGeneLine>>();
		if (intersectExons){
			//for each gene
			for (UCSCGeneLine gene: genes){
				ExonIntron[] exons = gene.getExons();
				//for each exon
				for (ExonIntron e: exons){
					//the end is included in IntervalST so sub 1 from end
					int start = e.getStart();
					int stop = e.getEnd() -1;
					Interval1D it = new Interval1D(start, stop);
					ArrayList<UCSCGeneLine> al;
					if (exonTree.contains(it)) al = exonTree.get(it);
					else {
						al = new ArrayList<UCSCGeneLine>();
						exonTree.put(it, al);
					}
					al.add(gene);
				}
			}
		}
		else {
			//for each gene
			for (UCSCGeneLine gene: genes){
				//the end is included in IntervalST so sub 1 from end
				int start = gene.getTxStart();
				int stop = gene.getTxEnd() -1;
				Interval1D it = new Interval1D(start, stop);
				ArrayList<UCSCGeneLine> al;
				if (exonTree.contains(it)) al = exonTree.get(it);
				else {
					al = new ArrayList<UCSCGeneLine>();
					exonTree.put(it, al);
				}
				al.add(gene);
			}
		}
		return exonTree;
	}

	private void parseFiles() throws Exception {
		UCSCGeneModelTableReader tr = new UCSCGeneModelTableReader(ucscTableFile, 0);
		chrGenes = tr.getChromSpecificGeneLines();
		if (verbose)  IO.pl(tr.getGeneLines().length+"\tGenes parsed");

		//make writer
		out = new Gzipper(resultsFile);
		//add header if present to output
		String header = IO.fetchHeader(bedFile);
		if (header != null && header.length()!=0) out.print(header);

		//parse bed
		Bed[] regions = Bed.parseFilePutLineInNameNoScoreOrStrand(bedFile, bpPad, (-1*bpPad), coorIndexes);
		totalNumberRegions = regions.length;
		if (verbose)  IO.pl(totalNumberRegions+"\tRegions to intersect");
		if (totalNumberRegions != 0) {
			Arrays.sort(regions);
			chrBed = Bed.splitBedByChrom(regions);
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AnnotateBedWithGenes(args); 
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'u': ucscTableFile = new File (args[++i]); break;
					case 'b': bedFile = new File (args[++i]); break;
					case 'r': resultsFile = new File (args[++i]); break;
					case 'p': bpPad = Integer.parseInt(args[++i]); break;
					case 'g': intersectExons = false; break;
					case 'f': addGeneFeature = true; break;
					case 'i': coorIndexes = Num.parseInts(args[++i], Misc.COMMA); break;
					case 'h': printDocs(); System.exit(0);
					default: IO.pl("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (ucscTableFile == null || ucscTableFile.canRead() == false) Misc.printExit("\nError: cannot find your UCSC table file!\n");
		if (bedFile == null || bedFile.canRead() == false) Misc.printExit("\nError: cannot find your bed file!\n");
		if (resultsFile == null) Misc.printExit("\nError: please provide a file path to save the gzipped results.\n");
	}	

	public static void printDocs(){ 
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Annotate Bed With Genes   Nov 2022                     **\n" +
				"**************************************************************************************\n" +
				"Takes a bed like file and a UCSC gene table, intersects them and adds a new column to\n"+
				"the file with the gene names that intersect the gene exons or regions. \n\n"+

				"Parameters:\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (geneName name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). \n"+
				"-b Bed like file of regions to intersect with genes, gz/zip OK\n"+
				"-i Indexes defining chr start stop columns, defaults to 0,1,2 for bed format.\n"+
				"-r Gzipped results file.\n"+
				"-p Bp padding to expand the bed regions when intersecting with genes.\n"+
				"-g Intersect gene regions with bed, not gene exons.\n"+
				"-f Include gene feature in annotations (5pUTR->3pUTR->All Exons->Intron)\n"+


				"\nExample: java -Xmx2G -jar pathTo/USeq/Apps/AnnotateBedWithGenes -p 100 -g -i 1,2,3\n" +
				"      -b targetRegions.bed -r targetsWithGenes.txt.gz -u hg19EnsGenes.ucsc.gz -f\n"+

				"**************************************************************************************\n");		
	}		
}

