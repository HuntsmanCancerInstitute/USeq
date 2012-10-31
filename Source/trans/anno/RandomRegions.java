package trans.anno;
import java.io.*;

import util.bio.annotation.*;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;

/**Use to generate random genomicRegions of a defined length from a file of genomicRegions (chrom, start, stop).
 * No attention is paid to the chromosome so if you want chromosome matching then split the genomicRegions file by 
 * chromosome and make a seperate RandomRegions object.*/
public class RandomRegions {
	
	//fields
	private GenomicRegion[] genomicRegions;
	private Random randomGenerator;
	private int divider = 10;	//used to convert the length of a region into array pointers, set to 1 for no length bias, higher for faster processing
	private File genomicSequenceDirectory = null;
	private double fractionGCTolerance = 0.025;		
	
	//constructors
	public RandomRegions(File regionsFile){
		//make genomicRegions
		makeRegions(regionsFile);
		//make random number generator
		randomGenerator = new Random();
	}
	public RandomRegions(File regionsFile, File genomicSequenceDirectory){
		this.genomicSequenceDirectory = genomicSequenceDirectory;	
		//make genomicRegions
		makeRegions(regionsFile);
		//make random number generator
		randomGenerator = new Random();
	}
	
	public static void main(String[] args){
		RandomRegions rr = new RandomRegions(new File(args[0]));
		int size = 100;
		int number = 10;
		int[][] randomRegions = rr.fetchRandomCoordinates(size, number);
		System.out.println("\nRandom Regions:");
		for (int i=0; i< randomRegions.length; i++){
			System.out.println(randomRegions[i][0]+"\t"+randomRegions[i][1]);
		}
	}
	
	/**Looks for a binary version of the file xxx.corr, if found loads binary, otherwise it parses the txt file 
	 * and then writes the binary for future use.*/
	public static GenomicRegion[] loadWriteBinaryCoordinates(File regionsFile){
		if (regionsFile.getName().endsWith(".corr")){
			return loadBinaryCoordinatesAsRegionArray(regionsFile);
		}
		File binary = new File (IO.getFullPathName(regionsFile)+".corr");
		if (binary.exists()) return loadBinaryCoordinatesAsRegionArray(binary);
		else {
			GenomicRegion[] rs = GenomicRegion.parseRegions(regionsFile);
			System.out.println("\tSaving binary version of genomicRegions file -> "+binary);
			saveRegionsAsCoordinatesArray(binary, rs);
			return rs;
		}
	}
	
	/**Reads binary chrom start stop file into GenomicRegion[]*/
	public static GenomicRegion[] loadBinaryCoordinatesAsRegionArray(File file){
		Coordinate[] c = Coordinate.readBinary(file);
		GenomicRegion[] rs = new GenomicRegion[c.length];
		String chrom = Misc.removeExtension(file.getName());
		for (int i=0; i< c.length; i++){
			rs[i] = new GenomicRegion(chrom, c[i].getStart(), c[i].getStop(), null);
		}
		c = null;
		return rs;
	}
	
	/**Writes the chrom start stop of a GenomicRegion[] as a binary file.*/
	public static void saveRegionsAsCoordinatesArray(File file, GenomicRegion[] r){
		Coordinate[]  startStop = new Coordinate[r.length];
		for (int i=0; i< r.length; i++){
			startStop[i] = new Coordinate(r[i].getChromosome(), r[i].getStart(), r[i].getEnd());
		}
		Coordinate.writeBinary(startStop, file);
		startStop = null;
	}
	
	/**Given a desired length of sequence, returns multiple random region start stop coordinates.
	 * @return int[number of random genomicRegions][start stop], will return null if random genomicRegions of a given length could not be found.*/
	public int[][] fetchRandomCoordinates (int length, int number){
		int[][] rrs = new int[number][2];
		for (int i=0; i< number; i++){
			//randomly select a region of the correct size
			GenomicRegion r = findRandomRegion(length);
			if (r==null) return null;
			int[] startStop = findRandomSegment(r, length);
			rrs[i] = startStop;
		}
		return rrs;
	}
	
	/**Given a desired length of sequence, returns multiple random region start stop coordinates.
	 * @return int[number of random genomicRegions][start stop], will return null if random genomicRegions of a given length could not be found.*/
	public int[][] fetchRandomCoordinates (int length, int number, double fractionGCContent){
		int[][] rrs = new int[number][2];
		for (int i=0; i< number; i++){
			//attempt 10K times to find a length and GC matched random region
			for (int x =0; x<10000; x++){
				//randomly select a region of the correct size
				GenomicRegion r = findRandomRegion(length);
				if (r==null) return null;
				int[] startStop = findRandomSegment(r, length);
				//calc gc
				double total = 0;
				boolean[] gc = r.getGcContent();
				int stop = startStop[1]-r.getStart()+1;
				for (int y=startStop[0]-r.getStart(); y<stop; y++){
					if (gc[y]) total++;
				}
				double perGC = total/(double)length;
				double min = perGC - fractionGCTolerance;
				double max = perGC + fractionGCTolerance;
				//if gc within tolerances then break and get next one.
				if (fractionGCContent >= min && fractionGCContent <= max) {
					rrs[i] = startStop;
					break;
				}
			}
			//found a matched gc frag?
			if (rrs[i] == null) return null;
		}
		return rrs;
	}
	
	/**Returns a random segment of the appropriate length within the region.
	 * Assumes the GenomicRegion is larger than the length.*/
	public int[] findRandomSegment(GenomicRegion genomicRegion, int length){
		int[] startStop = new int[2];
		int size = genomicRegion.getEnd()-genomicRegion.getStart() - length + 1;
		startStop[0] = randomGenerator.nextInt(size) + genomicRegion.getStart();
		startStop[1] = startStop[0]+ length -1;
		return startStop;
	}
	
	/**Attempt to find a random GenomicRegion of the appropriate minimum length.
	 * Will return null after genomicRegions.length attempts. */
	public GenomicRegion findRandomRegion(int minimumLength){
		for (int x =0; x<genomicRegions.length; x++){
			int index = randomGenerator.nextInt(genomicRegions.length);
			int size = genomicRegions[index].getEnd() - genomicRegions[index].getStart();
			if (size >= minimumLength) return genomicRegions[index];
		}
		return null;
	}
	
	/**Converts the GenomicRegion[] sequences into boolean[]s, everything not g or c are recorded as false.
	 * If stop exceeds the length of the chromosome, gc trunkated.*/
	public boolean[][] fetchGCContent(GenomicRegion[] r){
		//load chromosome sequence
		String chromosomeSequence = null;
		//binary?
		File binarySeq = new File (genomicSequenceDirectory,r[0].getChromosome()+".binarySeq");
		if (binarySeq.exists()) chromosomeSequence = Seq.readBinarySequence(binarySeq);
		else {
			File chromFastaFile = new File (genomicSequenceDirectory,r[0].getChromosome()+".fasta");
			if (chromFastaFile.exists() ==  false )  chromFastaFile = new File (genomicSequenceDirectory,r[0].getChromosome()+".fasta.gz");
			if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,r[0].getChromosome()+".fa");
			if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,r[0].getChromosome()+".fa.gz");
			if (chromFastaFile.exists() ==  false ) Misc.printExit("Error: Cannot find genomic sequence file -> "+chromFastaFile +" Aborting.");
			MultiFastaParser fastaParser = new MultiFastaParser();
			fastaParser.parseIt(chromFastaFile);
			chromosomeSequence = fastaParser.getSeqs()[0];
		}
		boolean[][] gcContent = new boolean[r.length][];
		
		for (int i=0; i<r.length; i++) {
			int end = r[i].getEnd();
			int start = r[i].getStart();
			if (end > chromosomeSequence.length() || start > chromosomeSequence.length()) {
				int size = end-start;
				end = chromosomeSequence.length();
				start = end - size;
			}
			String seq = chromosomeSequence.substring(start, end).toLowerCase();
			char[] sequence = seq.toCharArray();
			gcContent[i] = new boolean[sequence.length];
			for (int j=0; j< sequence.length; j++){
				if (sequence[j] == 'g' || sequence[j] == 'c') gcContent[i][j] = true;
				else gcContent[i][j] = false;
			}
		}

		chromosomeSequence = null;
		return gcContent;
	}
	

	
	/**Parses the genomicRegions file but then adds a pointer per base pair length/ divider to
	 * the GenomicRegion[] for every region.  This is to prevent bias in randomly selecting
	 * genomicRegions.*/
	public void makeRegions(File regionsFile){
		//parse file
		System.out.println("\tLoading interrogated region coordinates for "+regionsFile.getName()+" ...");
		GenomicRegion[] rs = GenomicRegion.parseRegions(regionsFile);
		
		//check starts
		for (int i=0; i<rs.length; i++) {
			if (rs[i].getStart() < 0) rs[i].setStart(0);
		}
		
		//load genomicRegions with gc content?
		if (genomicSequenceDirectory != null){
			System.out.println("\tLoading gc content...");	
			boolean[][] gcContent = fetchGCContent(rs);
			for (int i=0; i<rs.length; i++) rs[i].setGcContent(gcContent[i]);
		}
		
		//make array of pointers, one for each base/ divider
		genomicRegions = new GenomicRegion[countLengths(rs)];
		int index = 0;
		for (int i=0; i<rs.length; i++){
			int length = (rs[i].getEnd() - rs[i].getStart())/divider;
			if (length < 1) length = 1;
			for (int j=0; j< length; j++) {
				genomicRegions[index] = rs[i];
				index++;
			}
		}
	}
	
	public int countLengths(GenomicRegion[] rs){
		int length = 0;
		for (int i=0; i< rs.length; i++){
			int size = (rs[i].getEnd() - rs[i].getStart())/divider;
			if (size < 1) size = 1;
			length += size;
		}
		return length;
	}
	
	public void setGenomicSequenceDirectory(File genomicSequenceDirectory) {
		this.genomicSequenceDirectory = genomicSequenceDirectory;
	}
	
	
	
}
