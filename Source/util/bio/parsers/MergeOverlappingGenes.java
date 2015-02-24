package util.bio.parsers;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.apps.FindOverlappingGenes;
import util.bio.annotation.ExonIntron;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Collapses transcripts into consensus gene based on exon overlap.*/
public class MergeOverlappingGenes {

	private File[] ucscGeneFiles;
	private File resultsFile;
	private double minimumFractionOverlap = 0.05;

	/**Download ensembl gene table from UCSC and replace the ENST text column with the name2.*/
	public  MergeOverlappingGenes(String[] args) {
		//process user arguments
		processArgs(args);
		
		//load models
		System.out.println("Parsing gene tables...");
		HashMap<String,UCSCGeneLine[]> chromGenes = loadGenes();

		//for each chrom strand
		int startingGeneCount = 0;
		int endingGeneCount = 0;
		System.out.println("Processing...");
		try {
			Gzipper out = new Gzipper(resultsFile);
			for (String chromStrand: chromGenes.keySet()){
				UCSCGeneLine[] g = chromGenes.get(chromStrand);
				startingGeneCount+= g.length;
				System.out.print("\t"+chromStrand+"\t"+g.length+" -> ");
				ArrayList<UCSCGeneLine> mergedGenes = merge(g, chromStrand);
				for (UCSCGeneLine m: mergedGenes) out.println(m);
				endingGeneCount+= mergedGenes.size();
				System.out.println(mergedGenes.size());
			}
			out.close();
		} catch (IOException e) {
			resultsFile.deleteOnExit();
			e.printStackTrace();
		}
		System.out.println("\n"+startingGeneCount + " genes collapsed to "+endingGeneCount+"\n");
	}

	private ArrayList<UCSCGeneLine> merge(UCSCGeneLine[] genes, String chromStrand) {
		//set tss to index and find min max
		int min = genes[0].getTxStart();
		int max = genes[0].getTxEnd();
		for (int i=0; i< genes.length; i++) {
			genes[i].setTss(i);
			if (genes[i].getTxStart() < min) min = genes[i].getTxStart();
			if (genes[i].getTxEnd() > max) max = genes[i].getTxEnd();
		}
		ArrayList<Integer>[] bpNames = new ArrayList[max-min];
		
		//for each gene set it's unique number at all genic bps
		for (int i=0; i< genes.length; i++) {
			int start = genes[i].getTxStart()-min;
			int stop = genes[i].getTxEnd()-min;
			Integer index = new Integer(genes[i].getTss());
			
			for (int j=start; j< stop; j++){
				ArrayList<Integer> al = bpNames[j];
				if (al == null) {
					al = new ArrayList<Integer>();
					bpNames[j] = al;
				}
				al.add(index);
			}
		}
		
		/*
		System.out.println("\nBlocks! "+min);
		for (int i=0; i< bpNames.length; i++){
			if (bpNames[i] != null){
				int[] geneIndexes = Misc.integerArrayListToIntArray(bpNames[i]);
				System.out.println((i+min) +"\t"+ Misc.intArrayToString(geneIndexes, ","));
			}
		}
		System.exit(0);
		*/
		
		
		//walk across each genic block
		HashSet<Integer> toMerge = new HashSet<Integer>();
		ArrayList<UCSCGeneLine> mergedGenes = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< bpNames.length; i++){
			//null or bp?
			if (bpNames[i] == null){
				//process and clear hash
				if (toMerge.size() !=0) {
					//System.out.println("Processing "+(i+min)+" "+toMerge);
					if (toMerge.size() == 1) {
						mergedGenes.add(genes[toMerge.iterator().next()]);
						toMerge.clear();
					}
					else mergedGenes.addAll(mergeGeneSet(toMerge, genes));
				}
			}
			else toMerge.addAll(bpNames[i]);
		}
		//process last
		if (toMerge.size() !=0) {
			//System.out.println("Processing last "+toMerge);
			if (toMerge.size() == 1) {
				mergedGenes.add(genes[toMerge.iterator().next()]);
				toMerge.clear();
			}
			else mergedGenes.addAll(mergeGeneSet(toMerge, genes));
		}	
		
//System.exit(0);
		
		return mergedGenes;
	}

	private ArrayList<UCSCGeneLine> mergeGeneSet(HashSet<Integer> toMerge, UCSCGeneLine[] genes) {
		int[] indexes = Num.hashSetToInt(toMerge);
		toMerge.clear();
//System.out.print("Indexes ");
//Misc.printArray(indexes);

		//must iteratively loop through until no more merging
		boolean genesMerged = true;
		while (genesMerged){
			genesMerged = false;
			//for each gene
			for (int i=0; i< indexes.length; i++){
				//already merged into another?
				if (genes[indexes[i]] == null) {
					//System.out.println("\n"+i+ " : i TopGene already merged skipping");
				}
				else {
					//System.out.println("\n"+i+" : i Working with TopGene "+genes[indexes[i]].getName());
					//for each subsequent gene
					for (int j=i+1; j< indexes.length; j++){
						//already merged?
						if (genes[indexes[j]] == null) {
							//System.out.println("\t"+j+" : j TestGene already merged skipping ");
						}
						else {
							//System.out.println("\t"+j+" : j Working with testGene "+genes[indexes[j]].getName());
							//calculate max fraction overlap
							double fraction = genes[indexes[i]].getLargestFractionExonicBpIntersection(genes[indexes[j]]);
							//System.out.println("\t\tfraction in "+fraction);
							if (fraction >= minimumFractionOverlap){
								//merge em replacing mergeGene
								//System.out.println("\t\tMerging");
								genes[indexes[i]].mergeWithOtherGene(genes[indexes[j]]);
								genes[indexes[j]] = null;
								genesMerged = true;
							}
							//else System.out.println("\t\tNot merging");
						}
					}
				}
			}
		}
		//return non nulls
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< indexes.length; i++) if (genes[indexes[i]] != null) al.add(genes[indexes[i]]);
		return al;
	}

	private HashMap<String,UCSCGeneLine[]> loadGenes() {
		UCSCGeneLine[][] all = new UCSCGeneLine[ucscGeneFiles.length][];
		int total = 0;
		UCSCGeneModelTableReader reader = null;
		for (int i=0; i< ucscGeneFiles.length; i++){
			reader = new UCSCGeneModelTableReader (ucscGeneFiles[i], 0);
			all[i] = reader.getGeneLines();
			total += all[i].length;
		}
		//make merge
		UCSCGeneLine[] concat = new UCSCGeneLine[total];
		int index = 0;
		for (int i=0; i< all.length; i++){
			for (int j=0; j< all[i].length; j++){
				concat[index++] = all[i][j];
			}
		}
		Arrays.sort(concat, new UCSCGeneLineChromComparator());
		reader.setGeneLines(concat);
		
		return reader.splitByChromosomeAndStrand();
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergeOverlappingGenes (args);
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
					case 'u': ucscGeneFiles = IO.extractFiles(args[++i]); break;
					case 'r': resultsFile = new File (args[++i]); break;
					case 'm': minimumFractionOverlap = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (ucscGeneFiles == null || ucscGeneFiles.length == 0) Misc.printErrAndExit("\nError: cannot find your UCSC formatted gene table(s)?\n");
		if (resultsFile == null) Misc.printErrAndExit("\nError: please provide a file to save the results.\n");

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Merge Overlappng Genes: Feb  2015                      **\n" +
				"**************************************************************************************\n" +
				"Merges transcript models that share exonic bps. Maximizes exons, minimizes introns.\n"+
				"Assumes interbase coordinates.\n\n" +

				"Options:\n"+
				"-u Path to a UCSC RefFlat or RefSeq gene table file or directory with such to merge.\n"+
				"       See http://genome.ucsc.edu/cgi-bin/hgTables, (geneName name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). \n"+
				"-r Path for results file.\n"+
				"-m Minimum fraction exonic bp overlap for merging, defaults to 0.05\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MergeOverlappingGenes -d \n" +
				"      /CufflinkTranscripts/zv9Genes.ucsc.gz -f 0.25 -r merged.ucsc\n\n" +

		"**************************************************************************************\n");

	}


}
