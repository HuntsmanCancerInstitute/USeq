package util.bio.seq;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;
import util.bio.parsers.*;
import util.bio.annotation.*;

/**Generates splices and transcripts from a gene table.
 * @author davidnix*/
public class MakeTranscriptome {

	//fields
	private File ucscTableFile = null;
	private File fastaFileDirectory = null;
	private HashMap<String,File> fastaFiles;
	private HashMap<String,UCSCGeneLine[]> chromSpecificGeneLines;
	private HashMap<String, ArrayList<UCSCGeneLine>> geneNameModelHM;
	private String chromosome;
	private String chrSequence;
	private UCSCGeneLine[] models;
	private Gzipper outSplices = null;
	private Gzipper outTranscripts = null;
	private int radius = 0;
	private int maxNumberSpicesPerTranscript = 100000;
	private double maxMinutesThreshold = 10;
	private int maxExonThreshold = 35;
	private HashMap<String, StringBuffer> skippedGenes = new HashMap<String, StringBuffer>();
	private HashSet<String> spliceCoordinates = new HashSet<String>(1000000);
	private boolean skipRepeatedSplices = false;
	private long numberSplices = 0;
	private long numberTranscripts = 0;

	public MakeTranscriptome(String[] args){
		try {

			long startTime = System.currentTimeMillis();

			//process user args and load data
			processArgs(args);

			File spliceFasta = new File(ucscTableFile.getParent(), Misc.removeExtension(ucscTableFile.getName())+"Rad"+radius+"Num"+(int)(maxNumberSpicesPerTranscript/1000)+"kMin"+(int)maxMinutesThreshold+"Splices.fasta.gz");
			File transcriptFasta = new File(ucscTableFile.getParent(), Misc.removeExtension(ucscTableFile.getName())+"Rad"+radius+"Num"+(int)(maxNumberSpicesPerTranscript/1000)+"kMin"+ (int)maxMinutesThreshold+"Transcripts.fasta.gz");
			System.out.println("\nSaving fasta files to "+ucscTableFile.getParent()+"\n");
			outSplices = new Gzipper (spliceFasta);
			outTranscripts = new Gzipper (transcriptFasta);


			//for each chromosome of gene models
			System.out.println("Processing...");
			Iterator<String> it = chromSpecificGeneLines.keySet().iterator();
			while (it.hasNext()){
				chromosome = it.next();
				System.out.print("\t" + chromosome+" ");

				//fetch the sequence
				MultiFastaParser mfp = new MultiFastaParser(fastaFiles.get(chromosome));
				chrSequence = mfp.getSeqs()[0];

				//for each gene line
				models = chromSpecificGeneLines.get(chromosome);

				//clear coordinates
				if (skipRepeatedSplices) spliceCoordinates.clear();

				//make transcripts
				addTranscripts();

				//make splice junctions
				addSpliceJunctions();

				System.out.println();
			}
			
			outSplices.close();
			outTranscripts.close();

			//any skipped genes?
			if (skippedGenes.size() !=0){
				StringBuilder genes = new StringBuilder();
				for (String name: skippedGenes.keySet()){
					genes.append(name);
					genes.append(":");
					genes.append(skippedGenes.get(name).toString());
					genes.append(" ");
				}
				System.out.println("\nWARNING: splice junctions for the following genes were only partially made. They exceeded the maximum number (N) or " +
						"minutes (M) thresholds. Only those present in the transcripts are guarenteed to be present in the fasta file.\n\n"+genes);
				//double fracSkipped = ((double)skippedGenes.size())/ ((double)geneNameModelHM.size());
				//System.out.println("\t"+Num.formatPercentOneFraction(fracSkipped)+" fraction genes with partial theoretical splice junctions. ");
			}

			System.out.println("\n"+numberTranscripts+"\t number transcripts created");
			System.out.println(numberSplices+"\t number splices created");


			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" min\n");

		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public void addSkippedGene (String gene, String note){
		StringBuffer sb;
		if (skippedGenes.containsKey(gene)) {
			sb = skippedGenes.get(gene);
			sb.append(",");
		}
		else {
			sb = new StringBuffer();
			skippedGenes.put(gene, sb);
		}
		sb.append(note);
		
	}

	private class SpliceMakerThread extends Thread {

		private String geneName;
		private HashSet<ConcatinatedSequence> uniqueSplices;
		private volatile boolean stop = false;

		public SpliceMakerThread (String geneName, HashSet<ConcatinatedSequence> uniqueSplices){
			this.geneName = geneName;
			this.uniqueSplices = uniqueSplices;
		}

		public void run(){
			ArrayList<UCSCGeneLine> transcripts = geneNameModelHM.get(geneName);

			//add known splices without expansion
			makeKnownSpliceJunctions(transcripts, uniqueSplices);
			if (stop) return;

			//add splices from each transcript model, this is needed to avoid cases where an exon intersects two others and leads to a collapse in the  merge
			UCSCGeneLine nextTranscript = transcripts.get(0);
			boolean go = true;
			long numMade = makeAllSpliceJunctionsWithCoordinates(nextTranscript.getExons(), nextTranscript.getChrom(), uniqueSplices);
			if (numMade != -1) {
				addSkippedGene(geneName, "1N");
				go = false;
			}
			if (stop) return;
			//more than one transcript? make merge 
			if (transcripts.size() > 1 && go){
				UCSCGeneLine mergeTranscript = nextTranscript.partialCloneForSpliceJunctionGeneration();
				for (int i=1; i< transcripts.size(); i++){
					nextTranscript = transcripts.get(i);
					//add splices, watch out for threshold
					numMade = makeAllSpliceJunctionsWithCoordinates(nextTranscript.getExons(), nextTranscript.getChrom(), uniqueSplices);
					if (numMade != -1) {
						addSkippedGene(geneName, "2N");
						go = false;
						break;
					}
					if (stop) return;
					//merge em with the running merged
					ExonIntron[] merged = ExonIntron.merge(mergeTranscript.getExons(), nextTranscript.getExons());
					mergeTranscript.setExons(merged);
					//reset tx start and stop
					if (mergeTranscript.getTxStart() > nextTranscript.getTxStart()) mergeTranscript.setTxStart(nextTranscript.getTxStart());
					if (mergeTranscript.getTxEnd() < nextTranscript.getTxEnd()) mergeTranscript.setTxEnd(nextTranscript.getTxEnd());
					//reset cds start stop
					if (mergeTranscript.getCdsStart() > nextTranscript.getCdsStart()) mergeTranscript.setCdsStart(nextTranscript.getCdsStart());
					if (mergeTranscript.getCdsEnd() < nextTranscript.getCdsEnd()) mergeTranscript.setCdsEnd(nextTranscript.getCdsEnd());
				}

				if (go) {
					//create container to generate multiple merged transcripts
					ExonIntron[] mergedExons = mergeTranscript.getExons();
					//load with sequence and this array
					for (int i=0; i< mergedExons.length; i++){
						mergedExons[i].setSequence(new String (chrSequence.substring(mergedExons[i].getStart(), mergedExons[i].getEnd())));
						mergedExons[i].setExonsIntrons(mergedExons);
						mergedExons[i].setIndex(i);
						mergedExons[i].setMaxNumberSplices(maxNumberSpicesPerTranscript);
					}
					//collect overlapping exons 
					IndexedExonIntron[] indexedExons = new IndexedExonIntron[mergedExons.length];
					//add merged exons
					for (int i=0; i< mergedExons.length; i++){
						ArrayList<ExonIntron> al = new ArrayList<ExonIntron>();
						al.add(mergedExons[i]);
						indexedExons[i] = new IndexedExonIntron(i, indexedExons,al, maxNumberSpicesPerTranscript );
					}
					//add alternative overlapping exons
					//for each transcript
					for (int i=0; i< transcripts.size(); i++){
						if (stop) return;
						//fetch exons
						ExonIntron[] exons = transcripts.get(i).getExons();
						//see which intersect
						for (int j=0; j< indexedExons.length; j++){
							//fetch alternatives for this index
							ArrayList<ExonIntron> alternativeExons = indexedExons[j].getAlternativeExons();
							//for each exon in the transcript
							for (ExonIntron transcriptExon: exons){
								//does it intersect the merged?
								if (mergedExons[j].intersects(transcriptExon.getStart(), transcriptExon.getEnd())){
									//has it already been added to the alternatives
									boolean notFound = true;
									for (ExonIntron alt : alternativeExons){
										if (transcriptExon.compareTo(alt) == 0) {
											notFound = false;
											break;
										}
									}
									//add if not found
									if (notFound) {
										alternativeExons.add(transcriptExon);
									}
								}
							}
						}
					}

					//fetch merged transcripts
					ArrayList<ArrayList<ExonIntron>> mergedTranscripts = indexedExons[0].fetchRight();
					if (mergedTranscripts == null || stop) {
						if (mergedTranscripts == null) {
							addSkippedGene(geneName, "eN");
						}
						return;
					}

					//for each one make splice junctions
					for (int i=0; i< mergedTranscripts.size(); i++) {
						ExonIntron[] ei = ExonIntron.arrayList2Array(mergedTranscripts.get(i));
						//System.out.println(geneName+"\ttrans"+i+"\t"+ExonIntron.toString(ei));
						//load with sequence and this array
						for (int j=0; j< ei.length; j++){
							ei[j].setSequence(new String( chrSequence.substring(ei[j].getStart(), ei[j].getEnd())));
							ei[j].setExonsIntrons(ei);
							ei[j].setIndex(j);
							ei[j].setMaxNumberSplices(maxNumberSpicesPerTranscript);
						}
						numMade = makeAllSpliceJunctionsWithCoordinates(ei, mergeTranscript.getChrom(), uniqueSplices);
						if (numMade != -1) {
							addSkippedGene(geneName, "3N");
							go=false;
							break;
						}
						if (stop) return;
					}
				}
			}
		}
	}

	public int maxExons (String geneName){
		ArrayList<UCSCGeneLine> transcripts = geneNameModelHM.get(geneName);
		int max = 0;
		for (UCSCGeneLine line : transcripts){
			int numExons = line.getExons().length;
			if (numExons > max) max = numExons;
		}
		return max;
	}


	public void makeKnownSpliceJunctions(ArrayList<UCSCGeneLine> transcripts, HashSet<ConcatinatedSequence> uniqueSplices){
		String chromosome = transcripts.get(0).getChrom();

		//for each transcript
		for (UCSCGeneLine transcript : transcripts){
			//fetch exons
			ExonIntron[] exons = transcript.getExons();
			//any junctions?
			if (exons.length == 1) continue;
			//for each exon pair
			for (int j=0; j< exons.length; j++){
				//at end?
				int nextIndex = j+1;
				if (nextIndex == exons.length) break;
				//make ConcatinatedSequence and add to hash
				ConcatinatedSequence cs = new ConcatinatedSequence (exons[j].fetch5PrimeSubSeqNoExpand(radius), exons[nextIndex].fetch3PrimeSubSeqNoExpand(radius), chromosome);
				uniqueSplices.add(cs);
			}
		}
	}

	public void makeSplicesNoThread(String geneName) throws IOException{
		ArrayList<UCSCGeneLine> transcripts = geneNameModelHM.get(geneName);
		HashSet<ConcatinatedSequence> uniqueSplices = new HashSet<ConcatinatedSequence>();

		//add known splices without expansion
		makeKnownSpliceJunctions(transcripts, uniqueSplices);

		//add splices from each transcript model, this is needed to avoid cases where an exon intersects two others and leads to a collapse in the  merge
		UCSCGeneLine nextTranscript = transcripts.get(0);
		boolean go = true;
		long numMade = makeAllSpliceJunctionsWithCoordinates(nextTranscript.getExons(), nextTranscript.getChrom(), uniqueSplices);
		if (numMade != -1) {
			addSkippedGene(geneName, "4N");
			go = false;
		}

		//more than one transcript? make merge 
		if (transcripts.size() > 1 && go){
			UCSCGeneLine mergeTranscript = nextTranscript.partialCloneForSpliceJunctionGeneration();
			for (int i=1; i< transcripts.size(); i++){
				nextTranscript = transcripts.get(i);
				//add splices, watch out for threshold
				numMade = makeAllSpliceJunctionsWithCoordinates(nextTranscript.getExons(), nextTranscript.getChrom(), uniqueSplices);
				if (numMade != -1) {
					addSkippedGene(geneName, "5N");
					go = false;
					break;
				}
				//merge em with the running merged
				ExonIntron[] merged = ExonIntron.merge(mergeTranscript.getExons(), nextTranscript.getExons());
				mergeTranscript.setExons(merged);
				//reset tx start and stop
				if (mergeTranscript.getTxStart() > nextTranscript.getTxStart()) mergeTranscript.setTxStart(nextTranscript.getTxStart());
				if (mergeTranscript.getTxEnd() < nextTranscript.getTxEnd()) mergeTranscript.setTxEnd(nextTranscript.getTxEnd());
				//reset cds start stop
				if (mergeTranscript.getCdsStart() > nextTranscript.getCdsStart()) mergeTranscript.setCdsStart(nextTranscript.getCdsStart());
				if (mergeTranscript.getCdsEnd() < nextTranscript.getCdsEnd()) mergeTranscript.setCdsEnd(nextTranscript.getCdsEnd());
			}

			if (go) {
				//create container to generate multiple merged transcripts
				ExonIntron[] mergedExons = mergeTranscript.getExons();
				//load with sequence and this array
				for (int i=0; i< mergedExons.length; i++){
					mergedExons[i].setSequence(new String (chrSequence.substring(mergedExons[i].getStart(), mergedExons[i].getEnd())));
					mergedExons[i].setExonsIntrons(mergedExons);
					mergedExons[i].setIndex(i);
					mergedExons[i].setMaxNumberSplices(maxNumberSpicesPerTranscript);
				}
				//collect overlapping exons 
				IndexedExonIntron[] indexedExons = new IndexedExonIntron[mergedExons.length];
				//add merged exons
				for (int i=0; i< mergedExons.length; i++){
					ArrayList<ExonIntron> al = new ArrayList<ExonIntron>();
					al.add(mergedExons[i]);
					indexedExons[i] = new IndexedExonIntron(i, indexedExons,al, maxNumberSpicesPerTranscript);
				}
				//add alternative overlapping exons
				//for each transcript
				for (int i=0; i< transcripts.size(); i++){
					//fetch exons
					ExonIntron[] exons = transcripts.get(i).getExons();
					//see which intersect
					for (int j=0; j< indexedExons.length; j++){
						//fetch alternatives for this index
						ArrayList<ExonIntron> alternativeExons = indexedExons[j].getAlternativeExons();
						//for each exon in the transcript
						for (ExonIntron transcriptExon: exons){
							//does it intersect the merged?
							if (mergedExons[j].intersects(transcriptExon.getStart(), transcriptExon.getEnd())){
								//has it already been added to the alternatives
								boolean notFound = true;
								for (ExonIntron alt : alternativeExons){
									if (transcriptExon.compareTo(alt) == 0) {
										notFound = false;
										break;
									}
								}
								//add if not found
								if (notFound) {
									alternativeExons.add(transcriptExon);
								}
							}
						}
					}
				}

				//fetch merged transcripts
				ArrayList<ArrayList<ExonIntron>> mergedTranscripts = indexedExons[0].fetchRight();
				if (mergedTranscripts != null) { 
					//for each one make splice junctions
					for (int i=0; i< mergedTranscripts.size(); i++) {
						ExonIntron[] ei = ExonIntron.arrayList2Array(mergedTranscripts.get(i));
						//load with sequence and this array
						for (int j=0; j< ei.length; j++){
							ei[j].setSequence(new String( chrSequence.substring(ei[j].getStart(), ei[j].getEnd())));
							ei[j].setExonsIntrons(ei);
							ei[j].setIndex(j);
							ei[j].setMaxNumberSplices(maxNumberSpicesPerTranscript);
						}
						numMade = makeAllSpliceJunctionsWithCoordinates(ei, mergeTranscript.getChrom(), uniqueSplices);
						if (numMade != -1) {
							addSkippedGene(geneName, "N");
							go=false;
							break;
						}
					}
				}
				else addSkippedGene(geneName, "eN"); 

			}
		}

		//scan for junctions entirely contained in another bigger one
		removeSmallSplices(uniqueSplices);

		//print out splices, this might be partial
		for (ConcatinatedSequence seq : uniqueSplices) {
			if (skipRepeatedSplices) {
				String coor = seq.getCoordinateString();
				if (spliceCoordinates.contains(coor)) continue;
				else spliceCoordinates.add(coor);
			}
			outSplices.println(seq.getFastaEntry(geneName));
			numberSplices++;
		}
	}


	public void addSpliceJunctions() throws IOException{
		try {
			//for each gene
			int counter = 0;
			for (String geneName : geneNameModelHM.keySet()){
				if (++counter == 100) {
					System.out.print(".");
					counter = 0;
				}
				//how many exons?
				int maxExons = maxExons(geneName);
				//use a timed thread? this is slow so best to avoid
				if (maxExons > maxExonThreshold) {
					//start a thread and timer, this is needed because some splices cause recursion so deep as to never return.
					HashSet<ConcatinatedSequence> uniqueSplices = new HashSet<ConcatinatedSequence>();
					SpliceMakerThread smt = new SpliceMakerThread(geneName, uniqueSplices);
					smt.start();
					long startTime = System.currentTimeMillis();
					while (smt.isAlive()){
						Thread.sleep(100);
						double min = ((double)(System.currentTimeMillis() -startTime)/60000.0);				
						//expired?
						if (min > maxMinutesThreshold){
							smt.stop = true;
							smt.interrupt();
							addSkippedGene(geneName, "M");
							break;
						}
					}

					removeSmallSplices(uniqueSplices);

					//print out splices
					for (ConcatinatedSequence seq : uniqueSplices) {
						if (skipRepeatedSplices) {
							String coor = seq.getCoordinateString();
							if (spliceCoordinates.contains(coor)) continue;
							else spliceCoordinates.add(coor);
						}
						outSplices.println(seq.getFastaEntry(geneName));
						numberSplices++;
					}



				}
				//nope just run it directly
				else makeSplicesNoThread(geneName);
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

	}



	public void removeSmallSplices(HashSet<ConcatinatedSequence> uniqueSplices){
		//scan for junctions entirely contained in another bigger one
		ArrayList<ConcatinatedSequence> toRemove = new ArrayList<ConcatinatedSequence>();
		for (ConcatinatedSequence seq1 : uniqueSplices){
			int lengthSeq1 = seq1.getSequence().length();
			for (ConcatinatedSequence seq2 : uniqueSplices){
				//make sure they are different lengths
				if (lengthSeq1 != seq2.getSequence().length()){
					//is it part of if
					if(seq2.getSequence().contains(seq1.getSequence())) {
						toRemove.add(seq1);
						break;
					}
				}
			}
		}
		//remove them
		for (ConcatinatedSequence s : toRemove) uniqueSplices.remove(s);
	}


	/**Makes all linear combinations of splice junctions given and array of non overlapping ExonIntrons.
	 * Will reach through additional exons if the first exon is too short.
	 * Set reverseComplement = true for negative stranded annotations
	 * @return -1 is ok, all splices made, otherwise the number is the number that would have been made but failed the maxSplicesPerTranscript*/
	public long makeAllSpliceJunctionsWithCoordinates(ExonIntron[] ei, String chromosome, HashSet<ConcatinatedSequence> uniqueSplices){
		//make all pairwise splice junctions
		for (int i=0; i< ei.length; i++){
			for (int j=i+1; j< ei.length; j++){
				//fetch seqs
				ArrayList<SubSequence>[] fivePrimeSeqs = ei[i].fetch5PrimeSubSeq(radius);
				if (fivePrimeSeqs == null) return 0;
				ArrayList<SubSequence>[] threePrimeSeqs = ei[j].fetch3PrimeSubSeq(radius);
				if (threePrimeSeqs == null) return 0;
				//join em
				for (ArrayList<SubSequence> fpSeq: fivePrimeSeqs){
					for (ArrayList<SubSequence> tpSeq: threePrimeSeqs){
						//make container
						ConcatinatedSequence csj = new ConcatinatedSequence(fpSeq, tpSeq, chromosome);
						uniqueSplices.add(csj);
					}
				}
				//check for too many
				if (uniqueSplices.size() >= maxNumberSpicesPerTranscript) return uniqueSplices.size();
			}
		}
		return -1;
	}


	public void addTranscripts(){
		try{
			geneNameModelHM = new HashMap<String, ArrayList<UCSCGeneLine>>();
			for (int i=0; i< models.length; i++){
				//fetch exons
				ExonIntron[] exons = models[i].getExons();
				//make text and strand
				String geneName = models[i].getDisplayName();
				String transcriptName = models[i].getName();
				if (geneName == null || transcriptName == null) Misc.printErrAndExit("\nAborting, missing gene and or transcript name from -> "+models[i]);
				//for each exon
				for (int j=0; j< exons.length; j++){
					//check end
					int end = exons[j].getEnd();
					if (end > chrSequence.length()) Misc.printErrAndExit("\nAborting, exon end is past the end of the sequence?! -> "+models[i]);
					String exonSeq = new String (chrSequence.substring(exons[j].getStart(), end));
					//save info for splices
					exons[j].setSequence(exonSeq);
					exons[j].setIndex(j);
					exons[j].setExonsIntrons(exons);
					exons[j].setMaxNumberSplices(maxNumberSpicesPerTranscript);
				}
				//make composite
				ConcatinatedSequence cs = new ConcatinatedSequence(exons, models[i].getChrom());
				//add
				outTranscripts.println(cs.getFastaEntry(geneName+ ":"+ transcriptName));
				numberTranscripts++;

				//save to hash map for splice junctions
				ArrayList<UCSCGeneLine> al = geneNameModelHM.get(geneName);
				if (al == null){
					al = new ArrayList<UCSCGeneLine>();
					geneNameModelHM.put(geneName, al);
				}
				al.add(models[i]);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MakeTranscriptome(args);
	}		

	/**This method will process each argument and assign new varibles*/
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
					case 'f': fastaFileDirectory = new File (args[++i]); break;
					case 'u': ucscTableFile = new File (args[++i]); break;
					case 'r': radius = Integer.parseInt(args[++i]);break;
					case 'n': maxNumberSpicesPerTranscript = Integer.parseInt(args[++i]);break;
					case 'm': maxMinutesThreshold = Double.parseDouble(args[++i]);break;
					case 's': skipRepeatedSplices = true; break;
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
		if (radius ==0) Misc.printErrAndExit("\nPlease enter a sequence length radius.\n");
		if (ucscTableFile == null || ucscTableFile.canRead() == false) Misc.printExit("\nCannot find or read your UCSC gene table?!\n");
		if (fastaFileDirectory == null || fastaFileDirectory.isDirectory() == false) Misc.printExit("\nCannot find your directory of chromosome specific xxx.fasta files?\n");

		//load UCSC table
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(ucscTableFile,0);
		chromSpecificGeneLines = reader.getChromSpecificGeneLines();

		//load fastaFiles
		fastaFiles = Seq.fetchChromosomeFastaFileHashMap(fastaFileDirectory);

		//check that all chroms are there in the fasta directory
		Iterator<String> it = chromSpecificGeneLines.keySet().iterator();
		StringBuilder missingChromosomes = new StringBuilder();
		while (it.hasNext()){
			String chr = it.next();
			if (fastaFiles.containsKey(chr) == false) {
				missingChromosomes.append(" ");
				missingChromosomes.append(chr);
			}
		}
		if (missingChromosomes.length() !=0) Misc.printExit("\nCould not find fasta sequence files for :"+missingChromosomes+"! \nEither delete the offending lines from your UCSC gene table or put an appropriate chromosme in the fasta directory.\n");

		//print args
		System.out.println("Transcript table\t"+ucscTableFile);
		System.out.println("Fasta directory\t"+fastaFileDirectory);
		System.out.println("Splice radius\t"+radius);
		System.out.println("Max # splices per gene\t"+maxNumberSpicesPerTranscript);
		System.out.println("Max minutes per splice\t"+(int)maxMinutesThreshold);
		System.out.println("Skip duplicate splice junctions\t"+skipRepeatedSplices);

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Make Transcriptome:  June 2012                         **\n" +
				"**************************************************************************************\n" +
				"Takes a UCSC ref flat table of transcripts and generates two multi fasta files of\n" +
				"transcripts and splices (known and theoretical). All possible unique splice junctions\n" +
				"are created given the exons from each gene's transcripts. In some cases this is\n" +
				"computationally intractable and theoretical splices from these are not complete.\n" +
				"Read through occurs with small exons to the next up or downstream so keep the sequence\n" +
				"length radius to a minimum to reduce the number of junctions. Overlapping exons are\n" +
				"assumed to be mutually exclusive. All sequence is from the plus genomic stand, no\n" +
				"reverse complementation. Interbase coordinates. This app can take a very long time to\n" +
				"run. Break up gene table by chromosome and run on a cluster. \n\n" +

				"To incorporate additional splice-junctions, add a new annotation line containing two\n" +
				"exons representing the junction to the table. If needed, set the -s option to skip\n" +
				"duplicates. \n\n"+

				"Options:\n"+
				"-f Fasta file directory, one per chromosome (e.g. chrX.fasta or chrX.fa, .gz/.zip OK)\n"+
				"-u UCSC RefFlat gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (geneName transcriptName chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 ENST00000329454 chr1 + \n" +
				"       16203317 16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 .\n"+
				"-r Sequence length radius. Set to the read length - 4bp.\n" +
				"-n Max number splices per transcript, defaults to 100000.\n"+
				"-m Max minutes to process each gene's splices before interrupting, defaults to 10.\n"+
				"-s Skip subsequent occurrences of splices with the same coordinates. Memory intensive.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MakeTranscriptome -f /Genomes/Hg18/Fastas/\n" +
				"      -u /Anno/Hg18/ensemblGenes.txt.ucsc -r 46 -s \n\n" +

		"************************************************************************************\n");

	}
}
