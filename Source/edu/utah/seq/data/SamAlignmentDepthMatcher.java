package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.ValidationStringency;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SamAlignmentDepthMatcher{
	//user defined fields
	private File toMatchBam;
	private File toSubSampleBam;
	private File finalMatchedSam;
	private File targetRegionsFile;
	
	//internal
	private Bed[] regions;
	private int numberThreads = 0;
	private File tempDir;
	private SamAlignmentDepthLoader[] loaders;
	private SamReaderFactory readerFactory;
	private ArrayList<File> samsToConcat = new ArrayList<File>();
	private double numPassing = 0;

	//constructors
	public SamAlignmentDepthMatcher(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			makeCheckBams();

			parseBams();

			System.out.println("\n\nConcatinating final unsorted gzipped sam...");
			IO.concatinateFiles(samsToConcat, finalMatchedSam);
			
			double fracPass = numPassing/ (double)regions.length;
			String fp = Num.formatNumber(fracPass, 3);
			IO.pl("\n"+fp+" of regions were successfully matched (>= 0.95 of target alignments extracted)");
			
			IO.deleteDirectory(tempDir);
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
			System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void makeCheckBams() throws IOException {
		readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

		//check that they are sorted by coordinate
		SamReader bamReader = readerFactory.open(toMatchBam);
		if (bamReader.getFileHeader().getSortOrder().equals(SortOrder.coordinate) == false || bamReader.hasIndex() == false) throw new IOException ("\nError: you bam file must be sorted by coordinate");
		
		bamReader.close();
		
		bamReader = readerFactory.open(toSubSampleBam);
		if (bamReader.getFileHeader().getSortOrder().equals(SortOrder.coordinate) == false || bamReader.hasIndex() == false) throw new IOException ("\nError: you bam file must be sorted by coordinate");
		
		SAMFileHeader h = bamReader.getFileHeader();
		String[] shLines = Misc.RETURN.split(h.getTextHeader()); //this gets the original header, setting SamSort.unsorted doesn't change it. Ugg.
		bamReader.close();
		for (int i=0; i< shLines.length; i++){
			if (shLines[i].contains("SO:coordinate")){
				shLines[i] = shLines[i].replace("SO:coordinate", "SO:unsorted");
				break;
			}
		}
		
		File header = new File(tempDir, "samHeader.sam.gz");
		Gzipper hOut = new Gzipper(header);
		for (String line: shLines) hOut.println(line);
		hOut.close();
		samsToConcat.add(header);
	}

	private void parseBams() throws Exception{

		//chunk regions
		regions = Bed.parseFile(targetRegionsFile, 0, 0);
		Misc.randomize(regions, System.currentTimeMillis());
		int numRegiosPerChunk = (int)Math.round((double)regions.length / (double)numberThreads);
		Bed[][] chunks = Bed.chunk(regions, numRegiosPerChunk);
		System.out.print(regions.length+" regions to parse with "+ numberThreads+ " matchers");
		
		//make workers
		loaders = new SamAlignmentDepthLoader[chunks.length];
		for (int i=0; i< loaders.length; i++) loaders[i] = new SamAlignmentDepthLoader(i, chunks[i], this);

		//load mpileup data
		ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
		for (SamAlignmentDepthLoader l: loaders) executor.execute(l);
		executor.shutdown();
		while (!executor.isTerminated()) {}

		//check loaders and fetch gzipped sams
		for (SamAlignmentDepthLoader l: loaders) {
			if (l.isFailed()) throw new IOException("ERROR: SamReadDepthLoader issue! \n"+l.getErrorMessage());	
			samsToConcat.add(l.getSam());
			numPassing += l.getNumPassing();
		}
		
		Arrays.sort(regions);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamAlignmentDepthMatcher(args);
	}		

	/**This method will process each argument and assign new varibles
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': toMatchBam = new File(args[++i]); break;
					case 's': toSubSampleBam = new File(args[++i]); break;
					case 'o': finalMatchedSam = new File(args[++i]); break;
					case 'b': targetRegionsFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (toMatchBam == null || toMatchBam.exists() == false || toMatchBam.getName().endsWith(".bam") == false) {
			Misc.printErrAndExit("\nPlease provide a xxx.bam file to use in matching alignment depths.");
		}
		if (toSubSampleBam == null || toSubSampleBam.exists() == false || toSubSampleBam.getName().endsWith(".bam") == false) {
			Misc.printErrAndExit("\nPlease provide a xxx.bam file to subsample and write out.");
		}
		if (finalMatchedSam == null || finalMatchedSam.getName().endsWith(".sam.gz") == false) {
			Misc.printErrAndExit("\nPlease provide a xxx.sam.gz file to use in writing out the extracted alignments.");
		}
		if (targetRegionsFile == null || targetRegionsFile.exists() == false ) {
			Misc.printErrAndExit("\nPlease provide a bed like regions file to match alignments over.");
		}
		tempDir = new File(finalMatchedSam.getCanonicalFile().getParentFile(), "TempSADM_"+Misc.getRandomString(6));
		tempDir.mkdirs();
		
		//threads to use
		numberThreads = Runtime.getRuntime().availableProcessors() - 1;
	}	
	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            SamAlignmentDepthMatcher: Oct 2018                    **\n" +
				"**************************************************************************************\n" +
				"Performs a region by region alignment depth subsampling to output a sam file with\n"+
				"matching alignment depths. Alignments are extracted, mates matched, randomized, then\n"+
				"saved to match the count from the other file. Uses all threads available.\n"+

				"\nRequired options:\n"+
				"-m Bam file, coordinate sorted with index, to use in calculating the alignment depths\n"+
				"      to match.\n" +
				"-s Bam file, coordinate sorted with index, to subsample to match -m and save. Be sure\n"+
				"      this is significantly bigger than -m .\n"+
				"-o Gzipped sam file to write the unsorted alignments from -s, must end in xxx.sam.gz\n"+
				"-b Bed file of regions in which to match alignment depths, xxx.bed(.gz/.zip OK).\n"+

				"\nExample: java -Xmx25G -jar pathToUSeq/Apps/SamAlignmentDepthMatcher -m tumor.bam\n" +
				"      -s bigMockTumor.bam -o matchedMock.sam.gz -b exomeCaptureTargets.bed.gz\n\n" +


		        "**************************************************************************************\n");

	}

	public File getToMatchBam() {
		return toMatchBam;
	}

	public File getToSubSampleBam() {
		return toSubSampleBam;
	}

	public File getTempDir() {
		return tempDir;
	}






	public SamReaderFactory getReaderFactory() {
		return readerFactory;
	}

}
