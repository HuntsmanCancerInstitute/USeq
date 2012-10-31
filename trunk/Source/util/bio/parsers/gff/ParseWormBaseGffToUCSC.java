package util.bio.parsers.gff;
import java.io.*;
import java.util.regex.*;

import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;
import java.util.*;

public class ParseWormBaseGffToUCSC{

	private File[] gffFiles;
	private Gff3Feature[] gffs;
	private Gff3Feature[] genes;
	private Gff3Feature[] mRNAs;
	private Gff3Feature[] exons;
	private boolean subtractOneFromEnds = true;
	private String typesToSave = "gene|mRNA|exon|five_prime_UTR|three_prime_UTR";
	private GeneExons[] geneExons;

	public ParseWormBaseGffToUCSC (String[] args){
		processArgs(args);
		
		//for each file
		for (int i=0; i< gffFiles.length; i++){

			//parse file setting types to extract
			parseGff(gffFiles[i]);

			//split to three types gene, mRNA, CDS
			System.out.println("Splitting gff...");
			splitGffToTypes();

			//make and load genes
			System.out.println("Loading genes and exons...");
			makeAndLoadGeneExons();

			//collapse exons into composite
			System.out.println("Merging exons...");
			collapseExons();

			//print in ucsc table format
			System.out.println("Printing UCSC table...");
			printUCSCFormat(gffFiles[i]);
		}
		System.out.println("Done!\n");
	}


	public void splitGffToTypes(){
		//split into genes, mRNAs, and exons/UTRs lines
		ArrayList genesAL = new ArrayList();
		ArrayList mRNAsAL = new ArrayList();
		ArrayList exonsAL = new ArrayList();
		for (int x=0; x< gffs.length; x++){			
			String type = gffs[x].getType();
			if (type.equals("exon") || type.endsWith("prime_UTR")) exonsAL.add(gffs[x]);
			else if (type.equals("mRNA")) mRNAsAL.add(gffs[x]);
			else if (type.equals("gene")) genesAL.add(gffs[x]);
		}
		genes = new Gff3Feature[genesAL.size()];
		genesAL.toArray(genes);
		mRNAs = new Gff3Feature[mRNAsAL.size()];
		mRNAsAL.toArray(mRNAs);
		exons = new Gff3Feature[exonsAL.size()];
		exonsAL.toArray(exons);
	}

	public void parseGff(File file){
		Gff3Parser parser = new Gff3Parser();
		parser.setRegExTypes(typesToSave);
		parser.setRelax(true);
		if (parser.parseIt(file) == false) Misc.printExit("Problem parsing gff! "+file);
		//subtract one from coordinates to put into zero base
		if (subtractOneFromEnds) parser.subtractOneFromFeatures();
		gffs = parser.getFeatures();
	}

	private class GeneExons{
		ArrayList exons = new ArrayList();
		int[][] mergedExons;
		Gff3Feature gene;
		public GeneExons(Gff3Feature gene){
			this.gene=gene;
		}
	}

	public void makeAndLoadGeneExons(){
		geneExons = new GeneExons[genes.length];
		//for each gene
		for (int i=0; i< genes.length; i++){
			//make GeneExons
			geneExons[i] = new GeneExons(genes[i]);
			//find all mRNAs where mRNA.Parent=Gene.ID
			ArrayList mRNAsAL = new ArrayList();
			String id = genes[i].getId();
			for (int j=0; j< mRNAs.length; j++){
				if (id.equals(mRNAs[j].getParent())) mRNAsAL.add(mRNAs[j]);
			}
			//any mRNAs?
			if (mRNAsAL.size() !=0){
				//for each mRNA, find it's assoicated exon
				for (int j=0; j< mRNAsAL.size(); j++){
					Gff3Feature mRNA = (Gff3Feature) mRNAsAL.get(j);
					//parse out id
					String mRNAID = mRNA.getId();
					//for each exon
					for (int k=0; k< exons.length; k++){
						if (mRNAID.equals(exons[k].getParent())) geneExons[i].exons.add(exons[k]);
					}
				}
			}
		}
	}

	public void collapseExons(){
		//for each GeneExons
		for (int i=0; i< geneExons.length; i++){
			if (geneExons[i].exons.size() !=0){
				Gff3Feature[] gfs = new Gff3Feature[geneExons[i].exons.size()];
				geneExons[i].exons.toArray(gfs);
				int[][] collapsedExons = mergeFeatures(gfs);
				geneExons[i].mergedExons = collapsedExons;
			}
		}
	}


	/**Uses a boolean[] to score whether a base is part of an annotation (true) or not (false).
	 * Assumes stop inclusive.
	 * Assumes same chromosome.*/
	public static int[][] mergeFeatures(Gff3Feature[] features){
		//find largest and smallest bases
		int largest = 0;
		int smallest = 1000000000;
		for (int i=0; i< features.length; i++){
			int end = features[i].getEnd();
			if (end > largest) largest = end;
			int start = features[i].getStart();
			if (start < smallest) smallest = start;
		}
		//make boolean to hold
		int length = largest- smallest+1;
		boolean[] bps = new boolean[length];
		Arrays.fill(bps, true);
		//for each feature set booleans to false
		for (int i=0; i< features.length; i++){
			int end = features[i].getEnd() + 1 - smallest;
			int start = features[i].getStart() - smallest;
			for (int j= start; j< end; j++){
				bps[j] = false;
			}
		}
		//find blocks
		int[][] startStop = ExportIntergenicRegions.fetchFalseBlocks(bps,0,0);

		//add smallest to each stop to get back to real coordinates
		for (int j=0; j< startStop.length; j++){
			startStop[j][0] += smallest;
			startStop[j][1] += smallest;
		}

		return startStop;
	}


	/**	Print each cluster in ucsc table format
		refGene.common text refGene.name refGene.chrom refGene.strand refGene.txStart refGene.txEnd refGene.cdsStart refGene.cdsEnd	refGene.exonCount	refGene.exonStarts	refGene.exonEnds
		 p53 NM_198576	chr1	+	995569	1031415	995619	1030284	36	995569,997647,1010723,1016111,1016522,1016827,1017258,1018541,1018840,1019125,1019411,1019636,1020463,1020661,1021035,1021266,1021462,1021699,1022122,1022629,1022875,1023078,1023314,1024169,1024538,1024868,1025205,1025535,1025729,1026028,1026555,1026755,1027030,1029055,1029750,1030126,	995820,997909,1010771,1016327,1016747,1017052,1017465,1018760,1019035,1019326,1019560,1019742,1020580,1020826,1021179,1021391,1021568,1022038,1022260,1022757,1022990,1023198,1023668,1024362,1024754,1025098,1025340,1025632,1025894,1026140,1026672,1026948,1027118,1029280,1029854,1031415,	agrin	AGRIN	Agrin is a neuronal aggregating factor that induces the aggregation of ...
	 */
	public void printUCSCFormat(File gffFile){
		try{
			File ucscFile = new File (gffFile.getParentFile(), Misc.removeExtension(gffFile.getName())+".ucsc");
			PrintWriter out = new PrintWriter (new FileWriter (ucscFile));
			//for each GeneExons
			for (int i=0; i< geneExons.length; i++){
				//collect exon starts and ends, adding one to ends to bring to interbase numbering
				int[][] startStops = geneExons[i].mergedExons;
				int[] starts;
				int[] ends;
				//any startStops?
				if (startStops !=null){
					starts = new int[startStops.length];
					ends = new int[startStops.length];
					for (int j=0; j< startStops.length; j++){
						starts[j] = startStops[j][0];
						ends[j]= startStops[j][1]+1;	//adding one
					}
				}
				//assign one big exon
				else {
					starts = new int[] {geneExons[i].gene.getStart()};
					ends = new int[] {geneExons[i].gene.getEnd()+1};	//adding one
				}
				//what to call common text?
				String commonName;
				String uniqueName = geneExons[i].gene.getId().replaceFirst("Gene:", "");
				if (uniqueName.startsWith("WBGene")) commonName = uniqueName;
				else commonName = geneExons[i].gene.getSource();

				out.println(
						commonName+"\t"+ 
						uniqueName+"\t"+
						geneExons[i].gene.getSeqId()+"\t"+ 
						geneExons[i].gene.getStrand()+"\t"+ 
						geneExons[i].gene.getStart()+"\t"+ 
						//add one to get to interbase numbering
						(geneExons[i].gene.getEnd()+1)+"\t"+ 
						starts[0]+"\t"+
						ends[ends.length-1]+"\t"+
						ends.length+"\t"+
						Misc.intArrayToString(starts, ",")+"\t"+ 
						Misc.intArrayToString(ends, ","));
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ParseWormBaseGffToUCSC(args);
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
					case 'g': gffFiles = IO.extractFiles(new File (args[i+1])); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (gffFiles == null || gffFiles.length ==0) Misc.printExit("\nError: cannot find your gff file(s)!\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Parse WormBase Gff To UCSC    July 2007                    **\n" +
				"**************************************************************************************\n" +
				"Converts a modified wormbase gff3 file to ucsc table format for refseq.br.\n\n"+

				"Parameters:\n"+
				"-g Full path file text for a gff file or directory containing such.\n\n"+

		"**************************************************************************************\n");		
	}
}
