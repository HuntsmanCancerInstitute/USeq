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
	private Gzipper out;

	public AnnotateBedWithGenes (String[] args){
		try {
			//Parse args
			processArgs(args);

			//parse genes
			parseFiles();

			//for each bed chrom
			intersect();

			out.close();
		} catch (Exception e) {
			System.err.println("\nError annotating "+bedFile.getName());
			if (out != null) out.getGzipFile().delete();
			e.printStackTrace();
		}

	}

	private void intersect() throws IOException {
		HashSet<String> toAdd = new HashSet<String>();
		System.out.print("\n\t");
		//for each chrom of bed regions
		for (String chr: chrBed.keySet()){
			System.out.print(chr+" ");
			Bed[] bed = chrBed.get(chr);
			UCSCGeneLine[] genes = chrGenes.get(chr);
			if (genes == null) {
				System.err.print("\nWARNING: No genes for chrom "+chr+", skipping.\n\t");
				for (Bed region: bed) out.println(region.toString()+"\t.");
				continue;
			}

			//build interval tree of exons from the genes
			IntervalST<ArrayList<UCSCGeneLine>> exonTree = buildTree(genes);

			//for each bed region
			for (Bed region: bed){
				int start = region.getStart()-bpPad;
				int stop = region.getStop()+bpPad;
				if (start < 0) start = 0;
				//end is included in search so subtract 1
				Iterable<Interval1D> it = exonTree.searchAll(new Interval1D(start, stop-1));
				toAdd.clear();
				for (Interval1D x : it) {
					ArrayList<UCSCGeneLine> g = exonTree.get(x);
					//for each gen line extract wanted info and collapse with the Hash
					for (UCSCGeneLine gl: g){
						toAdd.add(gl.getNames(":"));
					}
				}
				//print it
				String anno = ".";
				if (toAdd.size() !=0) anno = Misc.hashSetToString(toAdd, ",");
				out.println(region.toString()+"\t"+anno);
			}	
		}
		System.out.println("\n\nDone!");
	}

	private IntervalST<ArrayList<UCSCGeneLine>> buildTree(UCSCGeneLine[] genes){
		IntervalST<ArrayList<UCSCGeneLine>> exonTree = new IntervalST<ArrayList<UCSCGeneLine>>();
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
		return exonTree;
	}

	private void parseFiles() throws Exception {
		UCSCGeneModelTableReader tr = new UCSCGeneModelTableReader(ucscTableFile, 0);
		chrGenes = tr.getChromSpecificGeneLines();
		System.out.println(tr.getGeneLines().length+"\tGenes parsed");

		//parse bed
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		Arrays.sort(regions);
		chrBed = Bed.splitBedByChrom(regions);
		System.out.println(regions.length+"\tRegions to intersect");

		//make writer
		out = new Gzipper(resultsFile);
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
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
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
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
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
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Annotate Bed With Genes   Sept 2017                    **\n" +
				"**************************************************************************************\n" +
				"Takes a bed file and a UCSC gene table, intersects them and adds a new column to the\n"+
				"bed with the gene names with intersecting exons. \n\n"+

				"Parameters:\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (geneName name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). \n"+
				"-b Bed file of regions to intersect with the gene table.\n"+
				"-r Gzipped results bed file.\n"+
				"-p Bp padding to add to the bed coordinates when intersecting with gene exons.\n"+

				"\nExample: java -Xmx2G -jar pathTo/USeq/Apps/AnnotateBedWithGenes -p 100 \n" +
				"      -b targetRegions.bed -r targetRegionsWithGenes.bed.gz -u hg19EnsGenes.ucsc.gz\n"+

				"**************************************************************************************\n");		
	}		



}

