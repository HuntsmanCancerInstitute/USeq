package util.bio.seq;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.useq.data.Region;


import util.gen.*;
import util.bio.parsers.*;

/**Tiles oligos across a list of genomic regions.*/
public class OligoTiler {

	//fields
	private File bedFile;
	private File resultsSeqFileForward;
	private File resultsSeqFileReverse;
	private File resultsBedFile;
	private HashMap<String,File> fastaFiles;
	private HashMap<String,Region[]> regions;
	private int threePrimeOffSet = 10;
	private int spacerSize = 25;
	private double halfSpacerSize;
	private int minimumRegionSize = 20;
	private int oligoSize = 50;
	private int oligoSizeMinusBarcode;
	private double halfOligoSize;
	private int totalNumberOligos = 0;
	private long totalBPTiled = 0;
	private boolean tileCpGs = false;
	private Pattern cpg = Pattern.compile("cg",Pattern.CASE_INSENSITIVE);
	private Pattern nonGATC = Pattern.compile("[^gatc]",Pattern.CASE_INSENSITIVE);
	private String currentSequence;
	private String currentChromosome;
	private int maxGap = 8;
	private boolean printNameSeqFile = true;
	private boolean alternateFR = true;
	private boolean barcode = false;
	private String barcodeSequence = "ccgatacgtcg";
	private PrintWriter outSeqForward;
	private PrintWriter outSeqReverse;
	private PrintWriter outBed;

	//constructor
	public OligoTiler(String[] args){
		processArgs(args);

		//load regions
		regions = Region.parseStartStops(bedFile, 0, 0, 0); 

		//for each chromosome
		try {
			//make print writers
			outSeqForward = new PrintWriter( new FileWriter (resultsSeqFileForward));
			outSeqReverse = new PrintWriter( new FileWriter (resultsSeqFileReverse));
			outBed = new PrintWriter( new FileWriter (resultsBedFile));

			Iterator<String> it = regions.keySet().iterator();
			while (it.hasNext()){
				currentChromosome = it.next();
				//get sequence
				File fasta = fastaFiles.get(currentChromosome);
				if (fasta == null || fasta.canRead()== false) Misc.printExit("\nCannot find "+currentChromosome+ ".fasta in your fasta directory, aborting");
				MultiFastaParser mfp = new MultiFastaParser(fasta);
				currentSequence = mfp.getSeqs()[0];

				//fetch centers
				int[] centers = fetchAllCenters();
				if (centers == null) continue;

				//print results
				if (printNameSeqFile) printTxtBed(centers);
				else printFastaBed(centers);

			}
			outSeqForward.close();
			outSeqReverse.close();
			outBed.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		System.out.println(totalNumberOligos+"\tNumber Oligos");
		System.out.println(totalBPTiled+"\tTotal bp of tiled regions");
		double bpT = totalBPTiled;
		double tO = totalNumberOligos;
		double bpPerOligo = bpT/tO;

		System.out.println(Num.formatNumber(bpPerOligo, 3)+"\tBP tiled per oligo");


	}

	/**Prints a fasta file and bed file for the tiled oligos*/
	public void printFastaBed (int[] centers){
		//number passing oligos
		int numOligos = 0;
		int seqLengthMinusOne = currentSequence.length()-1;
		boolean printForward = true;
		//for each center pull oligo
		for (int i=0; i< centers.length; i++){
			int start = (int)Math.round(centers[i]-halfOligoSize);
			int stop = (int)Math.round(centers[i]+halfOligoSize);

			//watch out for out of bounds start and stop
			if (start < 0 || stop > seqLengthMinusOne) continue;

			//check for non GATC characters in critical region
			String intSeq = currentSequence.substring(start,stop);
			Matcher badBase = nonGATC.matcher(intSeq);
			if (badBase.find()) continue;


			//collect forward oligo sequence			
			int a = start;
			if (a<0) continue;
			int b = stop+threePrimeOffSet;
			if (b > seqLengthMinusOne)continue;
			String forward = currentSequence.substring(a,b);
			forward = nonGATC.matcher(forward).replaceAll("a");
			String name = currentChromosome+":"+start+"-"+stop+"\n";

			//collect reverse 
			a = start-threePrimeOffSet;
			if (a<0) continue;
			b = stop;
			if (b > seqLengthMinusOne) continue;
			String reverse = currentSequence.substring(a, b);
			reverse = Seq.reverseComplementDNA(reverse);
			reverse = nonGATC.matcher(reverse).replaceAll("a");

			totalNumberOligos++;

			//add barcode?
			if (barcode) {
				forward = forward.substring(0, oligoSizeMinusBarcode)+barcodeSequence;
				reverse = reverse.substring(0, oligoSizeMinusBarcode)+barcodeSequence;
			}

			//print
			outBed.println(currentChromosome+"\t"+start+"\t"+stop+"\t"+intSeq);
			if (alternateFR){
				if (printForward){
					outSeqForward.println(">F_"+name+forward);
					outSeqReverse.println(">R_"+name+reverse);
					printForward = false;
				}
				else {
					outSeqReverse.println(">F_"+name+forward);
					outSeqForward.println(">R_"+name+reverse);
					printForward = true;
				}
			}
			else {
				outSeqForward.println(">F_"+name+forward);
				outSeqReverse.println(">R_"+name+reverse);
			}

		}
	}

	/**Prints a txt file and bed file for the tiled oligos*/
	public void printTxtBed (int[] centers){

		//number passing oligos
		int seqLengthMinusOne = currentSequence.length()-1;
		boolean printForward = true;
		//for each center pull oligo
		for (int i=0; i< centers.length; i++){
			int start = (int)Math.round(centers[i]-halfOligoSize);
			int stop = (int)Math.round(centers[i]+halfOligoSize);

			//check for non GATC characters in critical region
			String intSeq = currentSequence.substring(start,stop);
			Matcher badBase = nonGATC.matcher(intSeq);
			if (badBase.find()) continue;

			totalNumberOligos++;
			//collect forward oligo sequence
			int a = start;
			if (a<0) a=0;
			int b = stop+threePrimeOffSet;
			if (b > seqLengthMinusOne) b= seqLengthMinusOne;
			String forward = currentSequence.substring(a,b);
			forward = nonGATC.matcher(forward).replaceAll("a");
			String name = currentChromosome+":"+start+"-"+stop+"\t";

			//collect reverse 
			a = start-threePrimeOffSet;
			if (a<0) a=0;
			b = stop;
			if (b > seqLengthMinusOne) b= seqLengthMinusOne;
			String reverse = currentSequence.substring(a, b);
			reverse = Seq.reverseComplementDNA(reverse);
			reverse = nonGATC.matcher(reverse).replaceAll("a");

			//add barcode?
			if (barcode) {
				forward = forward.substring(0, oligoSizeMinusBarcode)+barcodeSequence;
				reverse = reverse.substring(0, oligoSizeMinusBarcode)+barcodeSequence;
			}

			//print
			outBed.println(currentChromosome+"\t"+start+"\t"+stop+"\t"+intSeq);
			if (alternateFR){
				if (printForward){
					outSeqForward.println("F_"+name+forward);
					outSeqReverse.println("R_"+name+reverse);
					printForward = false;
				}
				else {
					outSeqReverse.println("F_"+name+forward);
					outSeqForward.println("R_"+name+reverse);
					printForward = true;
				}
			}
			else {
				outSeqForward.println("F_"+name+forward);
				outSeqReverse.println("R_"+name+reverse);
			}
		}
	}

	/**For each region finds oligo centers.*/
	public int[] fetchAllCenters(){
		Region[] ss = regions.get(currentChromosome);
		ArrayList<Integer> centers = new ArrayList<Integer>();
		//for each region find centers
		if (tileCpGs) for (int i=0; i< ss.length; i++) centers.addAll(findCpGCenters(ss[i]));
		else for (int i=0; i< ss.length; i++) centers.addAll(findCenters(ss[i]));
		if (centers.size() !=0) return Num.arrayListOfIntegerToInts(centers);
		return null;

	}

	/**Within a given region finds oligo centers.*/
	public ArrayList<Integer> findCenters(Region ss){
		double regionSize = ss.getLength();
		ArrayList<Integer> centers = new ArrayList<Integer>();
		if (regionSize < minimumRegionSize){
			//System.out.println("Too small, skipping "+ss);
		}
		else if (regionSize <= spacerSize) {
			totalBPTiled+= regionSize;
			int halfSize = (int) Math.round(regionSize/2);
			//System.out.println("Smaller "+(ss.getStart()+halfSize));
			centers.add(new Integer(ss.getStart()+halfSize));
		}
		
		else {

			totalBPTiled+= regionSize;
			double halfRegionSize = regionSize/2.0;
			//System.out.println("Bigger");
			double numOligos = (int)(regionSize/spacerSize);
			if (numOligos == 1){
				centers.add(new Integer((int)Math.round(ss.getStart() + halfRegionSize)));
			}
			else if (numOligos <= 3){
				
				//add start
				int center = (int)Math.round(ss.getStart() + halfOligoSize);
				centers.add(new Integer(center));
				//add middle?
				if (numOligos == 3) {
					int testCenter = (int)Math.round(ss.getStart() + halfRegionSize);
					if (testCenter != center) {
						centers.add(new Integer(testCenter));
						center = testCenter;
					}
				}
				//add stop
				int end = new Integer((int)Math.round(ss.getStop() - halfOligoSize));
				if (end > 0 && end != center) centers.add(end);
			}
			//for greater than three
			else {
				//add start
				int center = (int)Math.round(ss.getStart() + halfOligoSize);
				centers.add(new Integer(center));
				double numberRemainingOligos = numOligos-2;
				double spacer = (regionSize-oligoSize)/ (numberRemainingOligos+1);
				//add first
				double first = (double)ss.getStart() + halfOligoSize + spacer;
				int firstInt = (int)Math.round(first);
				if (firstInt != center){
					centers.add(new Integer(firstInt));
					center = firstInt;
				}
				//add others
				for (int i=1; i< numberRemainingOligos; i++){
					first+= spacer;
					firstInt = (int)Math.round(first);
					if (firstInt != center){
						centers.add(new Integer(firstInt));
						center = firstInt;
					}
				}
				//add stop
				int end = new Integer((int)Math.round(ss.getStop() - halfOligoSize));
				if (end > 0 && end!= center) centers.add(end);
			}
		}
		return centers;
	}

	/**Places a center over each c in a CpG.*/
	public ArrayList<Integer> findCpGCenters(Region ss){
		double regionSize = ss.getLength();
		ArrayList<Integer> centers = new ArrayList<Integer>();
		if (regionSize >= minimumRegionSize){
			totalBPTiled+= regionSize;
			int start = ss.getStart();
			if (start >= currentSequence.length()) {
				Misc.printErrAndExit("Error: start of region exceeds end of sequence for "+ss.toString());
			}
			if (start < 0) start = 0;
			int end = ss.getStop();
			if (end >= currentSequence.length()) {
				Misc.printErrAndExit("Error: end of region exceeds end of sequence for "+ss.toString());
			}
			String subSeq = currentSequence.substring(start, end);
			Matcher mat = cpg.matcher(subSeq);
			while (mat.find()){
				//find first
				int first = mat.start();
				//look ahead
				int old = first;
				while (mat.find()){
					int next = mat.start();
					//keep going?
					if (next-first > maxGap){
						//set center
						int center = first+ (int)Math.round(((double)(old-first))/2) + ss.getStart();
						centers.add(new Integer(center));
						//reset
						first = next;
						old = next;
					}
					else {
						old = next;
					}
				}
				//set last
				int center;
				if (old!=first) center = first+ (int)Math.round(((double)(old-first))/2) + ss.getStart();
				else center = first + ss.getStart();
				centers.add(new Integer(center));
			}
		}
		return centers;
	}




	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new OligoTiler(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		File fastaDirectory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fastaDirectory = new File (args[++i]); break;
					case 'r': bedFile = new File (args[++i]); break;
					case 's': spacerSize = Integer.parseInt(args[++i]); break;
					case 'g': maxGap = Integer.parseInt(args[++i]); break;
					case 'o': oligoSize = Integer.parseInt(args[++i]); break;
					case 't': threePrimeOffSet = Integer.parseInt(args[++i]); break;
					case 'm': minimumRegionSize = Integer.parseInt(args[++i]); break;
					case 'b': barcode = true; break;
					case 'c': tileCpGs = true; break;
					case 'e': alternateFR = false; break;
					case 'a': printNameSeqFile = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//check params
		if (bedFile == null || bedFile.canRead() == false) Misc.printExit("\nCannot find or read your region/ bed file?!\n");
		if (fastaDirectory == null) Misc.printExit("\nCannot find your directory containing chromosome specific xxx.fasta(.zip/.gz) files?\n");

		//make fastas
		fastaFiles = Seq.fetchChromosomeFastaFileHashMap(fastaDirectory);

		//calculate sizes
		halfSpacerSize = ((double)spacerSize)/2;
		halfOligoSize = ((double)oligoSize)/2;
		oligoSizeMinusBarcode = oligoSize +threePrimeOffSet - barcodeSequence.length();

		//set extension for sequence files
		String ext;
		if (printNameSeqFile) ext = ".txt";
		else ext = ".fasta";

		//set gap
		String gap;
		if (tileCpGs) gap = maxGap+"CpG"+minimumRegionSize+"MinOT";
		else gap = spacerSize+"bp"+minimumRegionSize+"MinOT";

		//set For Rev
		String alt1;
		String alt2;
		if (alternateFR){
			alt1 = "_ForRev1_";
			alt2 = "_ForRev2_";
		}
		else {
			alt1 = "_For";
			alt2 = "_Rev";
		}

		//make results files
		String name = Misc.removeExtension(bedFile.toString());
		resultsSeqFileForward = new File (name+ alt1 +gap+ext);
		resultsSeqFileReverse = new File (name+ alt2 +gap+ext);
		resultsBedFile = new File (name+"_"+gap+".bed");


		//print notes
		if (tileCpGs) {
			System.out.println("Tiling CpGs");
			System.out.println(maxGap+ "\tMax gap between adjacent CpGs");
		}
		else System.out.println(spacerSize+ "\tSpacing");
		System.out.println(oligoSize+ "\tEffective Oligo Size");
		System.out.println(threePrimeOffSet+"\t3' offset");
		System.out.println((int)minimumRegionSize+ "\tMinimum region size");
		System.out.println(printNameSeqFile+"\tPrint text sequence files for Agilent eArray upload");
		System.out.println(alternateFR+"\tSplit results by alternating strand");
		System.out.println(barcode+"\tReplace 3' with barcode");
		System.out.println("Results Forward->\t"+resultsSeqFileForward);
		System.out.println("Results Reverse->\t"+resultsSeqFileReverse);
		System.out.println("Results Bed->\t"+resultsBedFile);
		System.out.println();

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Oligo Tiler: Oct 2009                               **\n" +
				"**************************************************************************************\n" +
				"OT tiles oligos across genomic regions returning their forward and reverse sequences.\n" +
				"Won't tile oligos with non GATC characters, case insensitive. Replaces non GATC chars\n" +
				"in offset regions with 'a'. Note, the defaults are set for generating a 60 mer Agilent\n" +
				"specific tiling microarray design where the first 10bp of the 3' stop are buried in the\n" +
				"matrix and the effective oligo length is 50bp. Adjust accordingly for other platforms.\n\n" +

				"Options:\n"+
				"-f Fasta file directory, should contain chromosome specific xxx.fasta files.\n" +
				"-r Regions file to tile (tab delimited: chr start stop ...) interbase coordinates.\n"+
				"-o Effective oligo size, defaults to 50.\n"+
				"-s Spacing to place oligos, defaults to 25.\n"+
				"-t Three prime offset, defaults to 10.\n"+
				"-m Minimum size of region to tile, defaults to 20.\n"+
				"-a Print oligo FASTA instead of an Agilent eArray text seq formatted results.\n"+
				"-c Tile CpG (spacing not used, see max gap option).\n"+
				"-g Max gap between adjacent CpGs to include in same oligo, defaults to 8.\n"+
				"-e Split export files by strand instead of alternating strand.\n"+
				"-b Replace 3' stop of oligos with the human 11-nullomer 'ccgatacgtcg'. The first\n" +
				"       ~10bp don't contribute to hybridization on Agilent arrays.\n"+

				"\n"+

				"Example: java -Xmx4000M -jar pathTo/Apps/OligoTiler -s 40 -f /Genomes/Hg18/Fastas/ \n" +
				"     -r /Designs/cancerArray.bed -p -a \n\n" +

		"************************************************************************************\n");

	}
}
