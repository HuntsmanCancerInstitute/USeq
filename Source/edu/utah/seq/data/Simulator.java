package edu.utah.seq.data;
import java.io.*;

import trans.anno.GenomicRegion;
import trans.anno.RegionComparator;
import trans.main.Interval;
import trans.misc.Util;
import util.gen.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.parsers.*;
import util.bio.seq.*;

public class Simulator {
	//fields
	private int numberSpikeIns = 1000;
	private int numberFragments = 1000;
	private int minSize = 150;
	private int maxSize = 350;
	private int[] fragmentSizeRange;
	//percent error per cycle
	private double[] error = {0.5,0.528,0.556,0.584,0.612,0.64,0.668,0.696,0.724,0.752,0.78,0.808,0.836,0.864,0.892,0.92,0.948,0.976,1.004,1.032,1.06,1.088,1.116,1.144,1.172,1.2,1.228,1.256,1.284,1.312,1.34,1.368,1.396,1.424,1.452,1.48,1.508,1.536,1.564,1.592,1.62,1.648,1.676,1.704,1.732,1.76,1.788,1.816,1.844,1.872,1.9,1.928,1.956,1.984,2.012,2.04,2.068,2.096,2.124,2.152,2.18,2.208,2.236,2.264,2.292,2.32,2.348,2.376,2.404,2.432,2.46,2.488,2.516,2.544,2.572,2.6,2.628,2.656,2.684,2.712,2.74,2.768,2.796,2.824,2.852,2.88,2.908,2.936,2.964,2.992,3.02,3.048,3.076,3.104,3.132,3.16,3.188,3.216,3.244,3.272};
	private int readLength = 26;
	private double maxFractionRepeat = 0.20;
	private double maxFractionN = 0.05;
	private File[] fastaFiles;
	private File resultsDirectory;
	private ArrayList<String> seqReads = new ArrayList();
	private ArrayList<String> tags = new ArrayList();
	private ArrayList<String> lines = new ArrayList();
	private HashMap<String,File> repeats = new HashMap();
	private GenomicRegion[] repeatRegions;
	private File[] repeatFiles;

	/**Constructor*/
	public Simulator(String[] args){
		processArgs(args);
		
		//generate a list of sequence fragments around each point
		System.out.println("Generating reads");
		generateFragments();
		
		//write out results
		System.out.println("Writing files");
		File tagFile = new File (resultsDirectory,"chIPPeaks.bed");
		IO.writeArrayList(tags, tagFile);
		
		File reads = new File (resultsDirectory, "reads.txt");
		IO.writeArrayList(seqReads, reads);
		
		File lineFile = new File (resultsDirectory, "frags.xls");
		IO.writeArrayList(lines, lineFile);
		
		//randomize the reads
		/*
		String[] seqs = new String[seqReads.size()];
		seqReads.toArray(seqs);
		Misc.randomize(seqs, System.currentTimeMillis());
		
		//split into different files
		int stop = numberSpikeIns;
		int length = seqs.length;
		while (true){
			//write reads to file
			ArrayList al = new ArrayList();
			for (int i=0; i<stop; i++) al.add(seqs[i]);
			File reads = new File (resultsDirectory, stop+"_Reads.txt");
			IO.writeArrayList(al, reads);
			//change start
			stop = stop*2;
			if (stop >= length) break;
		}*/
	}

	/**Randomly picks chromosomes for generating points weighted by file size.*/
	public int[] selectHitsToChromosomes(){
		//use file size as proxy for number of bases
		int[] fileSize = new int[fastaFiles.length];
		int minSize = (int)fastaFiles[0].length();
		for (int i=0; i< fileSize.length; i++) {
			fileSize[i] = (int)fastaFiles[i].length();
			if (fileSize[i] < minSize) minSize = fileSize[i];
		}		
		//make start and stop coordinates
		int[][] startStops = new int[fileSize.length][2];
		int oldIndex = 0;
		for (int i=0; i< fileSize.length; i++) {
			int newIndex = (fileSize[i]/minSize) + oldIndex;
			startStops[i] = new int[]{oldIndex,newIndex};
			oldIndex = newIndex;
		}
		//make int[] to track hits
		int[] chromsToHit = new int[fastaFiles.length];
		Random random = new Random();
		for (int i=0; i< numberSpikeIns; i++){
			//pick a random number
			int rndNum = random.nextInt(oldIndex);
			//find the appropriate chromosome
			for (int j=0; j<fastaFiles.length; j++){
				int[] startStop = startStops[j];
				if (rndNum >= startStop[0] && rndNum < startStop[1]) {
					chromsToHit[j]++;
					break;
				}
			}
		}
		return chromsToHit;
	}

	/**Pick points*/
	public void generateFragments(){
		//add header for lines
		lines.add("#chrom\tpos\tstart\tstop\t5'\t3'\tm5'\tm3'\tseq");
		//randomly select chromosomes for making spike ins, weighted by file size
		int[] chromHits = selectHitsToChromosomes();
		//for each chromosome
		for (int i=0; i< chromHits.length; i++){
			//fetch the fasta
			MultiFastaParser mfp = new MultiFastaParser(fastaFiles[i]);			
			//pick random sequences
			pickRandomFragments(mfp.getSeqs()[0], mfp.getNames()[0], chromHits[i]);
		}
	}

	/**Checks repeat content*/
	public boolean lotsOfRepeats(int minStart, int maxEnd){
		if (repeatRegions == null) return false;
		GenomicRegion test = new GenomicRegion(null, minStart, maxEnd-1, null);
		//run thru repeats
		boolean hit = false;
		double totalBpOverlap = 0;
		for (int j=0; j<repeatRegions.length; j++){
			int numBpInt = repeatRegions[j].bpIntersectionSameChromosome(test);
			//overlap
			if (numBpInt < 0 ){
				hit = true;
				totalBpOverlap += -1*numBpInt;
			}
			//no overlap and already hit thus past region
			else if (hit){
				break;
			}
		}
		//calculate fraction int?
		if (hit == false) return false;
		double fractionInt = totalBpOverlap/((double)(maxEnd-minStart));				
		if (fractionInt > maxFractionRepeat) return true;
		return false;
	}
	
	/**Check N content.*/
	public boolean lotsOfNs(String seq){
		char[] bases = seq.toCharArray();
		double numNs = 0;
		for (int i=0; i< bases.length; i++) if (bases[i] == 'n') numNs++;
		double fractN = numNs/((double)bases.length);
		if (fractN >= maxFractionN) return true;
		return false;
	}
	
	/**Picks fragments and reads.*/
	public void pickRandomFragments (String sequence, String chromName, int numberPoints){
		int length = sequence.length();
		int maxSizeFrag = fragmentSizeRange[1] - fragmentSizeRange[0];
		Random rand = new Random();
		System.out.println("\t"+chromName);		
		File repeatFile = repeats.get(chromName);
		if (repeatFile != null) {
			repeatRegions = GenomicRegion.parseRegions(repeatFile);
			//sort by position
			Arrays.sort(repeatRegions, new RegionComparator());
		}
		else repeatRegions = null;
	
		//for each point
		for (int i=0; i< numberPoints; i++){
			//pick a point
			int position = rand.nextInt(length);
			//is it too close to an stop?
			int minStart = position - maxSizeFrag;
			int maxEnd = position + maxSizeFrag;
			if (minStart < 0 || maxEnd >= length || lotsOfRepeats(minStart, maxEnd) || lotsOfNs(sequence.substring(minStart, maxEnd).toLowerCase())) i--;
			//make random fragments around position
			else {
				tags.add(chromName+"\t"+(position - maxSizeFrag)+"\t"+(position+maxSizeFrag));
				//make a fragment and make reads
				for (int x=0; x< numberFragments; x++){
					//pick a size
					int sizeFrag = rand.nextInt(maxSizeFrag)+ fragmentSizeRange[0];
					//pick where to begin
					int div = rand.nextInt(sizeFrag); 
					//set coordinates
					int start = position - div;
					int stop = position + sizeFrag - div;
					//cut frag and assign, be sure to use String constructor to avoid referencing big fasta
					String frag = new String(sequence.substring(start, stop));
					//trim reads off ends
					String a = frag.substring(0, readLength).toLowerCase();
					String b = Seq.reverseComplementDNA(frag).substring(0, readLength).toLowerCase();
					//add error
					String aMod = addError(a);
					String bMod = addError(b);
					seqReads.add(aMod);
					seqReads.add(bMod);
					lines.add(chromName+"\t"+position+"\t"+start+"\t"+stop+"\t"+a+"\t"+b+"\t"+aMod+"\t"+bMod+"\t"+frag);
				}
			}
		}
	}
	
	/**Make sure your read length doesn't exceed the length of the error*/
	public String addError(String read){
		char[] bases = {'g','a','t','c'};
		char[] seq = read.toCharArray();
		//for each base
		Random rand = new Random();
		for (int i=0; i < seq.length; i++){
			//pick a number 100000
			double hit = ((double)rand.nextInt(100000))/1000.0;
			//did it indicate an error?		
			if (hit <= error[i]){
				//modify base
				while (true) {
					int baseIndex = rand.nextInt(4);
					if (seq[i] != bases[baseIndex]) {
						seq[i] = bases[baseIndex];
						break;
					}
				}
			}
		}
		return new String(seq);
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Simulator(args);
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
					case 'f': fastaFiles = IO.extractFiles(new File (args[++i]),"fasta"); break;
					case 'n': numberSpikeIns = Integer.parseInt(args[++i]); break;
					case 'g': numberFragments = Integer.parseInt(args[++i]); break;
					case 's': minSize = Integer.parseInt(args[++i]); break;
					case 'x': maxSize = Integer.parseInt(args[++i]); break;
					case 'r': resultsDirectory = new File (args[++i]); break;
					case 'b': repeatFiles = IO.extractFiles(new File(args[++i]), "bed"); break;
					case 'l': readLength = Integer.parseInt(args[++i]); break;
					case 'e': error = Num.stringArrayToDouble(args[++i],",");break;
					
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		fragmentSizeRange = new int[]{minSize,maxSize};
		if (error.length<readLength)Misc.printErrAndExit("Error array is less than size of read?!");
		if (fastaFiles == null || fastaFiles.length == 0) Misc.printErrAndExit("Cannot find your xxx.fasta files?!");
		if (resultsDirectory == null) Misc.printErrAndExit("Enter a directory where you would like to save the results.");
		resultsDirectory.mkdir();
		if (repeatFiles == null || repeatFiles.length ==0) Misc.printErrAndExit("Cannot find your xxx.bed file containing repeats?");
		for (int i=0; i< repeatFiles.length; i++) repeats.put(Misc.removeExtension(repeatFiles[i].getName()), repeatFiles[i]);
	}	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Simulator: Nov 2008                                **\n" +
				"**************************************************************************************\n" +
				"Generates chIP-seq simulated sequences for aligning to a reference genome.\n\n" +

				"-f Directory containing xxx.fasta files with genomic sequence. File names should\n" +
				"     represent chromosome names (e.g. chr1.fasta, chrY.fasta...)\n" +
				"-r Results directory\n" +
				"-b Bed file containing repeat locations (e.g. RepeatMasker.bed)\n" +
				"-n Number of spike-ins, defaults to 1000\n" +
				"-g Number of random fragments to generate for each spike-in, defaults to 1000\n" +
				"-s Minimum size of a fragment, defaults to 150\n" +
				"-x Maximum size of a fragment, defaults to 350\n" +
				"-l Length of read, defaults to 26\n" +
				"-e Comma delimited text of per base % error rates, defaults to 0.5,0.528,0.556,...\n" +
				

				"\nExample: java -Xmx1500M -jar pathTo/USeq/Apps/Simulator -f /Hg18/Fastas -r /Spikes/\n" +
				"    -b /Hg18/Repeats/repMsker.bed -l 36\n\n" +

		"**************************************************************************************\n");
	}

}
