package edu.utah.seq.vcf.sim;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.annotation.Bed;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Creates snvs and indels for the BAMSurgeon.*/
public class VCFVariantMaker {
	//fields
	private File bedFile;
	private File vcfFile;
	private int numberSNVs = 0;
	private int numberINDELs = 0;
	private int maxIndelSize = 10;
	
	//internal
	private char[] gatc = {'G','A','T','C'};
	private Gzipper vcfOut;
	private Random random = new Random();
	private HashMap<String, Bed[]> chromRegions = null;
	private HashMap<String, File> fastaFiles;
	private String chromSeq = null;
	

	public VCFVariantMaker (String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		try {
			
			//load bed and check length
			pickRegionsToMutate();
			
			//make vcf writer and add header
			vcfOut = new Gzipper(vcfFile);
			vcfOut.println("##fileformat=VCFv4.1");
			vcfOut.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
			
			//for each chrom of region, mutate

			for (String chrom: chromRegions.keySet()){
				//fetch the sequence
				File f = fastaFiles.get(chrom);
				if (f == null) Misc.printErrAndExit("\nERROR: failed to find a fasta reference sequence for "+chrom);
				MultiFastaParser mfp = new MultiFastaParser(f);
				chromSeq = mfp.getSeqs()[0].toUpperCase();
				
				//for each region
				for (Bed bed : chromRegions.get(chrom)){
					if (bed.getName().equals("snv")) makeSnv(bed);
					else makeIndel(bed);
				}
			}
			
			//close gzipper
			vcfOut.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		} 

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void pickRegionsToMutate(){
		Bed[] bed = Bed.parseFile(bedFile, 0, 0);
		int totToMut = numberSNVs + numberINDELs;
		if (totToMut == 0){
			int half = bed.length/2;
			numberSNVs = half;
			numberINDELs = half;
			totToMut= half+half;
		}
		else if (bed.length< totToMut) Misc.printErrAndExit("\nError: More mutations requested than number of regions? Aborting!\n");
		System.out.println("\t"+numberSNVs+"\t#SNVs\n\t"+numberINDELs+"\t#INDELs");
		Misc.randomize(bed, System.currentTimeMillis());
		Bed[] cut = new Bed[totToMut];
		System.arraycopy(bed, 0, cut, 0, totToMut);
		for (int i=0; i< numberSNVs; i++) cut[i].setName("snv");
		for (int i=numberSNVs; i< cut.length; i++) cut[i].setName("indel");
		Arrays.sort(cut);
		chromRegions = Bed.splitBedByChrom(cut);
	}

	private void makeSnv(Bed region) throws Exception{
		int size = region.getLength();
		int rnd = random.nextInt(size);
		int pos = region.getStart()+rnd;
		char refSeq = chromSeq.charAt(pos);
		char varSeq = Seq.pickRandomBase(refSeq);
		//#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
		vcfOut.println(region.getChromosome()+"\t"+(pos+1)+"\t.\t"+refSeq+"\t"+varSeq+"\t.\t.\t.");
	}
	
	private void makeIndel(Bed region) throws Exception{
		int size = region.getLength();
		int pos = random.nextInt(size);
		int start = region.getStart()+pos;
		char refSeq = chromSeq.charAt(start);
		int lengthOfIndel = random.nextInt(maxIndelSize) + 1;

		//insertion 
		if (random.nextBoolean()){
			String varSeq = fetchRandomSeq(lengthOfIndel);
			vcfOut.println(region.getChromosome()+"\t"+(start+1)+"\t.\t"+refSeq+"\t"+refSeq+varSeq+"\t.\t.\t.");
		}
		//deletion
		else {
			int stop = start + lengthOfIndel +1;
			String refSeqToDel = chromSeq.substring(start, stop);
			vcfOut.println(region.getChromosome()+"\t"+(start+1)+"\t.\t"+refSeqToDel+"\t"+refSeq+"\t.\t.\t.");
		}
	}
	
	private String fetchRandomSeq(int length){
		char[] s = new char[length];
		for (int i=0; i< length; i++) s[i] = fetchRandomBase();
		return new String(s);
	}
	
	private char fetchRandomBase(){
		int index = random.nextInt(4);
		return gatc[index];
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFVariantMaker(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File fastaDirectory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFile = new File(args[++i]); break;
					case 's': numberSNVs = Integer.parseInt(args[++i]); break;
					case 'i': numberINDELs = Integer.parseInt(args[++i]); break;
					case 'm': maxIndelSize = Integer.parseInt(args[++i]); break;
					case 'v': vcfFile = new File(args[++i]); break;
					case 'f': fastaDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (bedFile == null || bedFile.canRead()== false) Misc.printErrAndExit("\nPlease enter a bed file of regions to target for mutation.\n");
		if (fastaDirectory == null) Misc.printErrAndExit("\nPlease provide a directory containing fasta files, one per chrom, gz/zip OK.\n");
		//load fastaFiles
		fastaFiles = Seq.fetchChromosomeFastaFileHashMap(fastaDirectory);
		if (fastaFiles.size() == 0) Misc.printErrAndExit("\nFailed to parse any fasta files from your directory? "+fastaDirectory);
		
		if (vcfFile == null) {
			String name = Misc.removeExtension(bedFile.getName());
			vcfFile = new File (bedFile.getParentFile(), name+"_Mut.vcf.gz");
		}
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Vcf Mutant Maker  : June 2015                          **\n" +
				"**************************************************************************************\n" +
				"Generates a VCF file of random SNVs and INDELs for BamBlaster. Only one variant is\n"+
				"made per input region. Be sure to left align this vcf file before using, see\n"+
				"https://github.com/atks/vt .\n" +

				"\nRequired:\n"+
				"-b Full path bed file containing regions in which to create mutations.\n"+
				"-f Path to a directory of chrom specific reference fasta files (gzip/zip OK, \n"+
				"      e.g. 1.fa.gz, 2.fa.gz, X.fa.gz...)\n\n"+

				"Optional:\n" +
				"-s Number of SNVs to make, defaults to 1/2 the number of regions.\n"+
				"-i Number of INDELs to make, ditto.\n"+
				"-m Max indel size, defaults to 10.\n"+
				"-v VCF results file, defaults to a permutation of the bed file.\n"+
				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/VCFMutantMaker -b ~/Sim/highCov.bed \n"+
				"     -f ~/Sim/human_g1k_v37.fasta.gz -m 12 -v ~/Sim/unNormalizedVMM.vcf\n\n"+

				"**************************************************************************************\n");

	}
}
