package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SamSubsampler{
	//user defined fields
	private File[] samFiles;
	private File saveDirectory;
	private float minimumPosteriorProbability = 13;
	private float maximumAlignmentScore = 300;
	private String adapter = "chrAdapt";
	private String phiX = "chrPhiX";
	private String lambda = "chrLamb";
	private int numberOfAlignmentsToPrint = 0;
	private double fractionAlignmentsToPrint = 0;
	private boolean sortFinal = false;
	private boolean applyFilters = true;
	private boolean verbose = true;

	//internal fields
	private Gzipper[] gzippers = null;
	private Random random = new Random();
	private int numberChunks = 100;
	private String samHeader;

	//alignment counts for sam files
	private long numberAlignmentsFailingQualityScore = 0;
	private long numberAlignmentsFailingAlignmentScore = 0;
	private long numberControlAlignments = 0;
	private long numberAlignmentsFailingQC = 0;
	private long numberAlignmentsUnmapped = 0;
	private long numberPassingAlignments = 0;



	//constructors
	public SamSubsampler(String[] args){
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		doWork();
		
		//write out final
		long toSave = numberPassingAlignments;
		File f = new File (saveDirectory, "randomized"+toSave+".sam.gz");
		if (fractionAlignmentsToPrint != 0) numberOfAlignmentsToPrint = (int)Math.round(fractionAlignmentsToPrint * (double) numberPassingAlignments); 
		printRandomAlignments(toSave, f);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}





	private void doWork() {
		if (verbose) System.out.println("Filtering and splitting alignments into "+numberChunks+" chunks...");

		//make Gzippers
		makeGzippers();

		//for each file, filter, split into numberChunks chunks
		for (int i=0; i< samFiles.length; i++){
			if (verbose) System.out.print("\t"+samFiles[i]);
			if (parseWorkingSAMFile(samFiles[i]) == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
		}

		//close Gzippers
		closeGzippers();

		//Alignment filtering stats
		if (verbose){
			double total = numberAlignmentsFailingQualityScore + numberAlignmentsFailingAlignmentScore + numberControlAlignments + numberAlignmentsFailingQC + numberAlignmentsUnmapped + numberPassingAlignments;
			System.out.println("\nFiltering statistics for "+(int)total+" alignments:");
			System.out.println(numberAlignmentsFailingQualityScore +"\tFailed mapping quality score ("+minimumPosteriorProbability+")");
			System.out.println(numberAlignmentsFailingAlignmentScore +"\tFailed alignment score ("+maximumAlignmentScore+")");
			System.out.println(numberControlAlignments +"\tAligned to chrPhiX*, chrAdapt*, or chrLamb*");
			System.out.println(numberAlignmentsFailingQC +"\tFailed vendor QC");
			System.out.println(numberAlignmentsUnmapped +"\tAre unmapped\n");
			System.out.println(numberPassingAlignments +"\tPassed filters ("+Num.formatPercentOneFraction(((double)numberPassingAlignments)/ total)+")");
		}
		
		//randomize each
		randomizeChunks();

		//ok now ready to write out
	}





	public SamSubsampler(File bam, boolean verbose) {
		this.verbose = verbose;
		applyFilters = false;
		sortFinal = false;
		samFiles = new File[]{bam};
		doWork();
		//now call printRandomAlignments()
	}



	

	public void printRandomAlignments(long toSave, File outputGzSamFile) {
		try {
			//do they want a diff number than all?
			if (numberOfAlignmentsToPrint != 0){
				//too many requested
				if (numberOfAlignmentsToPrint> numberPassingAlignments){
					Misc.printErrAndExit("\nERROR: The number of passing alignments is less than the requested number to save!  Saving only "+numberPassingAlignments+" alignments.\n");
				}
				else toSave = numberOfAlignmentsToPrint;
			}
			//make readers
			BufferedReader[] readers = new BufferedReader[numberChunks];
			for (int i=0; i<readers.length; i++) readers[i] = IO.fetchBufferedReader(gzippers[i].getGzipFile());

			if (verbose) System.out.println("Writing randomized alignments...");
			

			Gzipper out = new Gzipper(outputGzSamFile);
			out.println(samHeader);

			for (int i=0; i< toSave; i++){
				//pick random reader
				while (true){
					int num = random.nextInt(readers.length);
					BufferedReader test = readers[num];
					String line = test.readLine();
					if (line == null){
						//close old
						test.close();
						//remove it from the readers
						BufferedReader[] newReaders = new BufferedReader[readers.length -1];
						for (int j=0; j< num; j++) newReaders[j] = readers[j];
						for (int j=num+1; j< readers.length; j++) newReaders[j-1] = readers[j];
						readers = newReaders;
					}
					else {
						out.println(line);
						break;
					}
				}
			}

			//close final sam
			out.close();

			//close readers
			for (int i=0; i<readers.length; i++) if (readers[i]!= null) readers[i].close();

			//sort and write bam?
			if (sortFinal){
				if (verbose) System.out.println("Sorting and indexing...");
				File bam = new File (saveDirectory, "randomized"+toSave+".bam");
				//sort and convert to BAM
				new PicardSortSam (outputGzSamFile, bam);
				outputGzSamFile.delete();
			}

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: problem writing final merged sam file\n");
		} 
	}





	private void randomizeChunks() {
		try {
			if (verbose) System.out.print("\nRandomizing chunks");
			for (int i=0; i< numberChunks; i++){
				if (verbose) System.out.print(".");
				//load chunk
				File ori = gzippers[i].getGzipFile();
				String[] lines = IO.loadFile(ori);
				//randomize and print
				Misc.randomize(lines, System.currentTimeMillis());
				File rnd = new File(ori.getParentFile(), i+"rand.sam.gz");
				rnd.deleteOnExit();
				gzippers[i] = new Gzipper(rnd);
				gzippers[i].println(lines);
				gzippers[i].close();
				//delete temp file
				ori.delete();
			}
			if (verbose) System.out.println();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError writting randomized chunk.\n");
		} 
	}





	private void makeGzippers() {
		try {
			gzippers = new Gzipper[numberChunks];
			File tempDir = new File (saveDirectory, "TempFiles");
			tempDir.mkdirs();
			tempDir.deleteOnExit();
			for (int i=0; i< gzippers.length; i++) {
				File t = new File (tempDir, i+"_temp.sam.gz");
				t.deleteOnExit();
				gzippers[i] = new Gzipper(t);
			}
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem making temporary gzippers\n");
		} 
	}

	private void closeGzippers() {
		try {
			for (int i=0; i< gzippers.length; i++) gzippers[i].close();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem closing temporary gzippers\n");
		} 
	}

	private Gzipper fetchRandomGzipper(){
		int randInt = random.nextInt(numberChunks);
		return gzippers[randInt];
	}





	public boolean parseWorkingSAMFile(File workingFile){
		try{
			//make reader
			SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(workingFile);
			SAMRecordIterator iterator = reader.iterator();
			if (samHeader == null) samHeader = reader.getFileHeader().getSAMString().trim();

			int counter =0;
			int numBadLines = 0;
			//for each record
			while (iterator.hasNext()){
				SAMRecord samRecord = iterator.next();

				//print status blip
				if (++counter == 1000000){
					if (verbose) System.out.print(".");
					counter = 0;
				}

				//this is a bit inefficient but gives absolute control on the sam data
				SamAlignment sa;
				String samLine = samRecord.getSAMString().trim();
				try {
					sa = new SamAlignment(samLine, false);
				} catch (Exception e) {
					if (verbose) System.out.println("\nSkipping malformed sam alignment ->\n"+samRecord.getSAMString()+"\n"+e.getMessage());
					if (numBadLines++ > 1000) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
					continue;
				}

				if (applyFilters){
					//is it aligned?
					if (sa.isUnmapped()) {
						numberAlignmentsUnmapped++;
						continue;
					}

					//does it pass the vendor qc?
					if (sa.failedQC()) {
						numberAlignmentsFailingQC++;
						continue;
					}

					//skip phiX, adapter, lambda
					String chr = sa.getReferenceSequence();
					if (chr.startsWith(phiX) || chr.startsWith(adapter) || chr.startsWith(lambda)) {
						numberControlAlignments++;
						continue;
					}

					//does it pass the scores threshold?
					if (sa.getAlignmentScore() > maximumAlignmentScore) {
						numberAlignmentsFailingAlignmentScore++;
						continue;
					}
					if (sa.getMappingQuality() < minimumPosteriorProbability) {
						numberAlignmentsFailingQualityScore++;
						continue;
					}
				}

				//increment counter
				numberPassingAlignments++;

				//save it to a random gzipper
				fetchRandomGzipper().println(samLine);

			}
			reader.close();
			if (verbose) System.out.println();
		} catch (Exception e){
			System.err.println("\nError parsing Novoalign file or writing split binary chromosome files.\nToo many open files? Too many chromosomes? " +
					"If so then login as root and set the default higher using the ulimit command (e.g. ulimit -n 10000)\n");
			e.printStackTrace();
			return false;
		}
		return true;
	}






	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamSubsampler(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': forExtraction = new File(args[++i]); break;
					case 'r': saveDirectory = new File(args[++i]); saveDirectory.mkdir(); break;
					case 'x': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'q': minimumPosteriorProbability = Float.parseFloat(args[++i]); break;
					case 'n': numberOfAlignmentsToPrint = Integer.parseInt(args[++i]); break;
					case 'f': fractionAlignmentsToPrint = Double.parseDouble(args[++i]); break;
					case 's': sortFinal = true; break;
					case 'b': applyFilters = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");
		samFiles = IO.collapseFileArray(tot);
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) file(s)!\n");

		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printErrAndExit("\nError: cannot find your alignment file(s)!\n");
		if (saveDirectory == null || saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nPlease enter a directory to use in saving your results.\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              SamSubsampler: June 2016                            **\n" +
				"**************************************************************************************\n" +
				"Filters, randomizes, subsamples and sorts sam/bam alignments, doesn't keep pairs.\n" +

				"\nOptions:\n"+
				"-a Alignment file or directory containing SAM/BAM (xxx.sam(.zip/.gz OK) or xxx.bam).\n" +
				"      Multiple files are merged.\n" +
				"-r Results directory.\n"+

				"\nDefault Options:\n"+
				"-n Number of alignments to print, defaults to all passing thresholds.\n"+
				"-f Fraction alignments to keep, defaults to 1.\n"+
				"-s Sort and index output alignments.\n"+
				"-x Maximum alignment score. Defaults to 300, smaller numbers are more stringent.\n"+
				"-q Minimum mapping quality score. Defaults to 13, bigger numbers are more stringent.\n" +
				"      For RNASeq data, set this to 0.\n" +
				"-b Bypass all filters and thresholds.\n"+

				"\nExample: java -Xmx25G -jar pathToUSeq/Apps/SamSubsampler -x 240 -q 20 -a\n" +
				"      /Novo/Run7/ -s /Novo/Run7/SR -f 0.05 \n\n" +


				"**************************************************************************************\n");

	}

}
