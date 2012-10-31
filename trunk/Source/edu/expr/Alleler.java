package edu.expr;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.parsers.*;
import util.bio.seq.Seq;
import util.gen.*;


/**Takes a table of alleles, fasta sequencing files, and gene models and scores the alleles for effect.
 * Assuming interbase coordinates and UCSC gene models. This is a horrendous app, needs a complete rewrite to incorporate all of the added
 * functionality into a unified framework.*/
public class Alleler {
	
	private File resultsFile;
	private File ucscGeneTableFile = null;
	private File alleleFile = null;
	private HashMap<String, File> fastas;
	private HashMap<String,UCSCGeneLine[]> geneModels;
	private HashMap<String,Allele[]> alleles;
	public HashSet<String> printedTranscripts = new HashSet<String>();
	private int neighborhood = 1000;
	private ArrayList<AffectedGene> affectedGenes = new ArrayList<AffectedGene>();
	private boolean printOnlyDeleteriousAlleles = false;
	private boolean printBedFile = false;
	private boolean collapse = false;
	
	public Alleler(String[] args){
		//process args, fastas, geneModels, and alleles
		processArgs(args);
		
		//scan each gene model for hits
		System.out.println("Scanning...");		
		scanGeneModels();
		
		//print out results
		printResults();
		
		System.out.println("\nDone!\n");
	}
	
	public void printResults(){
		//any hits?
		if (affectedGenes.size() == 0) {
			System.out.println("\nNo allelic hits to gene models.\n");
			return;
		}
		//just baddies?
		if (printOnlyDeleteriousAlleles) System.out.println("\nOnly printing non-synonymous and splice-site affector alleles...");
		else System.out.println("Printing all alleles...");
		try{
			PrintWriter out = new PrintWriter (new FileWriter (resultsFile));
			LinkedHashSet<String> set = new LinkedHashSet<String>();
			if (printBedFile){
				out.println("#Chr\tStart\tStop\tNames_AlleleTblNotes_Effect\tScore\tStrand");
			}
			else {
				out.println("GeneModels:\t"+ucscGeneTableFile.getName()+
					"\nAlleles:\t"+alleleFile.getName()+
					"\n\nResults are printed for each intersecting gene and include:\n" +
					"1) The UCSC gene model (Name(s),Chrom,GeneStart,GeneStop,Strand,TSS)\n" +
					"\t2) The intersecting Alleles (Start,End,Allele,Notes)\n" +
					"\t\t3) Their relative location and potential effects.\n");
			}
			for (int i=0; i< affectedGenes.size(); i++){
				if (printOnlyDeleteriousAlleles) {
					String res;
					if (printBedFile) {
						res = affectedGenes.get(i).toStringDeleteriousBedLine(collapse);
					}
					else res = affectedGenes.get(i).toStringDeleterious(collapse);
					if (res != null) {
						if (set.contains(res) == false){
							out.print(res);
							set.add(res);
						}
					}
				}
				else {
					String toPrint = null;
					if (printBedFile) {
						toPrint = affectedGenes.get(i).toStringBedLine(collapse);
					}
					else toPrint = affectedGenes.get(i).toString();
					if (set.contains(toPrint) == false){
						out.println(toPrint);
						set.add(toPrint);
					}
				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}

	}
	
	/**Associates alleles with gene models.*/
	public void scanGeneModels(){
		//for each chromosome of alleles
		Iterator<String> it = alleles.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			System.out.println("\t"+chrom);
			//Any gene models?
			if (geneModels.containsKey(chrom) == false){
				System.out.println("\tNo gene models for alleles in "+chrom);
				continue;
			}
			
			//get alleles
			Allele[] chromAlleles = alleles.get(chrom);
			
			//get sequence
			MultiFastaParser mfp = null;
			
			//Scan gene models
			UCSCGeneLine[] chromGeneModels = geneModels.get(chrom);
			for (int i=0; i< chromGeneModels.length; i++){
				//filter?
				if (chromGeneModels[i].getTxEnd() < 73412429 && chromGeneModels[i].getTxStart() > 87262584) continue;
				//get start and stop, adjust for neighborhood
				int start = chromGeneModels[i].getTxStart() - neighborhood;
				if (start < 0) start = 0;
				int stop = chromGeneModels[i].getTxEnd() + neighborhood;
				
				//run through each allele
				ArrayList<Allele> hits = new ArrayList<Allele>();
				for (int x=0; x< chromAlleles.length; x++){
					//if it intersects gene then add to AL
					if (chromAlleles[x].intersects(start, stop)) hits.add(chromAlleles[x]);
				}
				
				//any hits? add to effectedGenes
				if (hits.size() !=0){
					//load sequence?
					if (mfp == null){
						File fastaFile = fastas.get(chrom);
						if (fastaFile == null) Misc.printExit("\nError: cannot find fasta file "+chrom+".fasta ?\n");
						mfp = new MultiFastaParser(fastaFile);
					}
					Allele[] a = new Allele[hits.size()];
					hits.toArray(a);
					Arrays.sort(a);
					affectedGenes.add( new AffectedGene (chromGeneModels[i], a, mfp.getSeqs()[0], printedTranscripts));
				}
				
			}
		}
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Alleler(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File fastaDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'u': ucscGeneTableFile = new File(args[++i]); break;
					case 'a': alleleFile = new File(args[++i]); break;
					case 'g': fastaDir = new File(args[++i]); break;
					case 'n': neighborhood = Integer.parseInt(args[++i]); break;
					case 'd': printOnlyDeleteriousAlleles = true; break;
					case 'b': printBedFile = true; break;
					case 'c': collapse = true; break;
					case 'e': Misc.printExit("\n"+alleleTableExample); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (alleleFile == null || alleleFile.exists() == false) Misc.printExit("\nError: cannot find your allele table file. Please complete option -a\n\n"+alleleTableExample);
		if (ucscGeneTableFile == null || ucscGeneTableFile.exists() == false) Misc.printExit("\nError: cannot find your UCSC refseq/refflat gene table file. Please complete option -u\n");
		
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(ucscGeneTableFile, 0);
		reader.splitByChromosome();
		geneModels = reader.getChromSpecificGeneLines();
		
		if (printBedFile) resultsFile = new File (alleleFile.getParentFile(), Misc.removeExtension(alleleFile.getName())+"_Alleler.bed");
		else resultsFile = new File (alleleFile.getParentFile(), Misc.removeExtension(alleleFile.getName())+"_Alleler.txt");
		
		//parse genome files
		if (fastaDir == null) Misc.printExit("\nError: cannot find your fasta genome directory? Please complete option -g");
		fastas = Seq.fetchChromosomeFastaFileHashMap(fastaDir);
		
		//parse alleles
		alleles = Allele.parseAlleles(alleleFile);
	}	
	
	public static final String alleleTableExample =
		"# Allele table example\n"+
		"# Tab delimited: Chr Start Stop ReplaceWith Notes\n"+
		"# Interbase coordinates (1st base is zero, last base is not included, stop >= start)\n"+
		"# Alleles reference the plus strand\n"+
		"# Examples, modifications to start ATG in IL17D (uc001unm.1) Hg18:\n"+
		"\n"+
		"chr13	20176140	20176141	G	Single bp SNP replacing A with a G\n"+
		"chr13 20176140	20176140	T	Single insertion of a T before A\n"+
		"chr13	20176140	20176141		Single bp deletion of the A\n"+
		"chr13	20176140	20177140	CATTCA	1kb deletion with 6bp substitution of CATTCA\n"+
		"chr13	20176140	20176143	C	Replacement of the ATG with a C	\n"+
		"\n";

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Alleler:  Sept 2010                               **\n" +
				"**************************************************************************************\n" +
				"Intersects a list of alleles (SNPs and INDELs) with gene models and returns their\n" +
				"affects on coding sequences and splice-junctions. Assumes interbase coordinates. If\n" +
				"ambiguious bases (ie R,Y,S,W,K,M) are provided the non-reference base is assumed.\n\n" +

				"Options:\n"+
				"-a Full path file text for a table of alleles.\n"+
				"-e Print an example of an allele table.\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables\n"+
				"-g Full path directory text containing fasta files for reference base calling\n" +
				"      (e.g. chr1.fasta, chr5.fasta, ...).\n"+
				"-n Neighborhood to include in intergenic intersection, defaults to 1000\n"+
				"-d Print only non-synonymous and splice affector alleles, defaults to all.\n"+
				"-b Print results in bed format, defaults to detailed report.\n"+
				"-c Collapse multiple hits to the same gene producing the same variant.\n"+
				

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/Alleler -a /APCSeq/apcFam7Alleles.txt\n" +
				"      -u /Anno/ucscKnownGenes.txt -g /Anno/Hg18Fastas/ -n 5000 -d -b\n\n" +

		"**************************************************************************************\n");

	}
	
}
