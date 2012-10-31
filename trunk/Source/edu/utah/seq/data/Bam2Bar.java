package edu.utah.seq.data;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.BarParser;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import util.gen.*;

/**Per base read coverage
 * @author Nix
 * */
public class Bam2Bar {

	//fields
	private File[] bamFiles;
	private File bamFile;
	private File saveDirectory;
	private String versionedGenome;
	private boolean makeRelativeTracks = true;
	private float scalar = 0;
	private Pattern cigarSub = Pattern.compile("(\\d+)([MSDHN])");
	private Pattern supportedCigars = Pattern.compile(".*[^\\dMIDSHN].*");
	private double totalNumberAlignments = 0;
	private boolean stranded = false;
	private boolean workingWithNegativeStrand = false;

	//per chromosome fields
	private String chromosome;
	private Point point = null;
	private ArrayList<Point> points = new ArrayList<Point>();
	private int maxNumberBases;

	/**For stand alone app.*/
	public Bam2Bar(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);

		//for each bam file
		System.out.println("Processing:");
		for (int i=0; i< bamFiles.length; i++){
			bamFile = bamFiles[i];
			System.out.print("\t"+bamFile.getName()+" ");

			//make coverage track
			makeCoverageTrack();
			if (stranded) {
				workingWithNegativeStrand = true;
				makeCoverageTrack();
			}

			//save for scaling or write it out?
			if (makeRelativeTracks) scaleCoverageTrack();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public SAMRecord fetchFirstAlignment(SAMRecordIterator it){
		SAMRecord priorAlignment = null;
		if (stranded){
			while (true){
				priorAlignment = it.next();
				if (priorAlignment.getReadUnmappedFlag() == false && priorAlignment.getReadNegativeStrandFlag() == workingWithNegativeStrand) break;
				else priorAlignment = null;
			}
		}
		else {
			while (true){
				priorAlignment = it.next();
				if (priorAlignment.getReadUnmappedFlag() == false) break;
				else priorAlignment = null;
			}
		}
		return priorAlignment;
	}

	public void makeCoverageTrack(){		
		//make a reader
		SAMFileReader reader = new SAMFileReader(bamFile);
		//is it sorted
		if (reader.getFileHeader().getSortOrder().compareTo(SortOrder.coordinate) !=0) Misc.printErrAndExit("\nPlease sort your bam files by coordinate. Note SAMtools sort is broken, use Picard's SortSam app.\n");
		SAMRecordIterator it = reader.iterator();

		//fetch first alignment
		SAMRecord priorAlignment = fetchFirstAlignment(it);
		if (priorAlignment == null) Misc.printErrAndExit("\nNo aligned reads found in "+bamFile+" !?\n");

		chromosome = priorAlignment.getReferenceName();
		System.out.print("\n"+chromosome);		
		ArrayList<SAMRecord> alignments = new ArrayList<SAMRecord>();
		alignments.add(priorAlignment);
		int maxEnd = priorAlignment.getUnclippedEnd();

		//for each record
		while (it.hasNext()){
			SAMRecord nextAlignment = it.next();
			//pass flags?
			if (nextAlignment.getReadUnmappedFlag() == true || nextAlignment.getReadFailsVendorQualityCheckFlag() == true) continue;
			//stranded?
			if (stranded && priorAlignment.getReadNegativeStrandFlag() == workingWithNegativeStrand) continue;
			//different chromosome
			if (chromosome.equals(nextAlignment.getReferenceName()) == false){
				processAlignments(alignments);
				closeChromosome(nextAlignment.getReferenceName());
				maxEnd = nextAlignment.getUnclippedEnd();
			}
			//non intersecting and not adjacent
			else if (maxEnd < nextAlignment.getUnclippedStart()-1){
				processAlignments(alignments);
			}
			alignments.add(nextAlignment);
			priorAlignment = nextAlignment;
			if (nextAlignment.getUnclippedEnd() > maxEnd) maxEnd = nextAlignment.getUnclippedEnd();
		}	

		//process and close last chromosome
		processAlignments(alignments);
		closeChromosome(null);

		reader.close();
	}

	public static void printSam(SAMRecord sam){
		System.out.println(sam+" "+sam.getCigarString());
		System.out.println(sam.getAlignmentStart()+" - "+sam.getAlignmentEnd());
		System.out.println(sam.getUnclippedStart()+" - "+sam.getUnclippedEnd());

	}

	public void processAlignments(ArrayList<SAMRecord> alignments){
		totalNumberAlignments+= alignments.size();

		//make base hit count array
		int firstBase = alignments.get(0).getUnclippedStart() -1;
		int lastBase = -1;
		for (SAMRecord sam : alignments){
			int last = sam.getAlignmentEnd();
			if (last > lastBase) lastBase = last;
		}

		int[] baseCount = new int[lastBase-firstBase];

		//for each read in block
		for (SAMRecord sam : alignments){

			//check to see if cigar contains any unsupported characters
			String cigar = sam.getCigarString();
			checkCigar(cigar, sam);

			//get start relative to firstBase
			int start = sam.getUnclippedStart() -1 - firstBase;

			//for each block
			Matcher mat = cigarSub.matcher(cigar);
			while (mat.find()){
				String call = mat.group(2);
				int numberBases = Integer.parseInt(mat.group(1));
				//a match
				if (call.equals("M")) {
					for (int i = 0; i< numberBases; i++) baseCount[start++]++;
					if (numberBases > maxNumberBases) maxNumberBases = numberBases;
				}
				//just advance for all but insertions which should be skipped via the failure to match
				else start += numberBases;
			}

		}

		//save block, last base is included!

		//find first non zero point point
		int hits = 0;
		int basePosition = 0;
		int i = 0;
		for (; i< baseCount.length; i++){
			if (baseCount[i] !=0){
				hits = baseCount[i];
				basePosition = firstBase +i;
				break;
			}
		}


		//set beginning zero point?
		if (basePosition !=0) {
			int basePositionMinusOne = basePosition-1;
			int priorPosition = -1;
			if (point != null) priorPosition = point.getPosition();
			if (priorPosition != basePositionMinusOne) {
				point = new Point(basePositionMinusOne, 0);
				points.add(point);
			}
		}
		//set first data point
		point = new Point(basePosition,hits);
		points.add(point);

		//for each subsequent
		for (int m=i+1; m<baseCount.length; m++){
			//add Point only if hits differs
			if (baseCount[m] != hits){
				//add old if not the same position as prior
				if (point.getPosition() != basePosition){
					point = new Point(basePosition,hits);
					points.add(point);
				}			
				//set new
				basePosition = m + firstBase;
				hits = baseCount[m];
				point = new Point(basePosition,hits);
				points.add(point);	
			}
			else basePosition = m + firstBase;
		}

		//add last
		if (point.position != basePosition) {
			point = new Point(basePosition,hits);
			points.add(point);
		}

		//set ending point
		point = new Point(basePosition+1,0);
		points.add(point);

		//clear the alignments
		alignments.clear();
	}

	public void checkCigar(String cigar, SAMRecord sam){
		Matcher mat = supportedCigars.matcher(cigar);
		if (mat.matches()) Misc.printErrAndExit("\nUnsupported cigar string in -> \n"+sam.toString());
	}


	public void closeChromosome(String nextChromosome){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = new HashMap<String,String>();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color 
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#0066FF"); //light blue
		//what's the source
		map.put(BarParser.SOURCE_TAG, bamFile.getName());

		if (makeRelativeTracks) map.put(BarParser.DESCRIPTION_TAG, "Relative read coverage track (e.g. #Reads/TotalReads/1M).");
		else map.put(BarParser.DESCRIPTION_TAG, "Read coverage track.");
		//save in info
		Info info = new Info ("ReadCoverageTrack_"+bamFile.getName(), versionedGenome,  chromosome, ".", maxNumberBases, map);

		//write into PointData
		Point[] pts = new Point[points.size()];
		points.toArray(pts);

		PointData hitTrack = Point.extractPositionScores(pts);
		hitTrack.setInfo(info);

		//calculate interrogated region hits?
		//if (regions != null) calculateInterrogatedRegionCoverage(hitTrack);

		saveDirectory = new File (bamFile.getParentFile(), "ReadCoverage_"+Misc.removeExtension(bamFile.getName()));
		saveDirectory.mkdir();
		hitTrack.writePointData(saveDirectory);

		//set next chromosome
		chromosome = nextChromosome;
		if (nextChromosome != null) System.out.print(" "+ chromosome);	
		else System.out.println("\n");

		//reset points
		points.clear();
	}



	public void scaleCoverageTrack(){
		//calculate scalar
		scalar = (float)(totalNumberAlignments/ 1000000.0);
		//System.out.println("\t"+scalar+"\tScalar (#Matching BPs/1,000,000)");

		//for each chromosome load and scale
		HashMap<String, ArrayList<PointData>> pdAL = PointData.fetchPointData(saveDirectory);

		for (ArrayList<PointData> al : pdAL.values()){
			PointData pd = al.get(0);
			float[] scores = pd.getScores();
			for (int i=0; i< scores.length; i++){
				if (scores[i] !=0) scores[i] = scores[i]/scalar;
			}
			pd.writePointData(saveDirectory);
		}

		totalNumberAlignments = 0;
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Bam2Bar(args);
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
					case 'b': bamFiles = IO.extractFiles(args[++i], ".bam"); break;
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 'r': makeRelativeTracks = false; break;
					case 's': stranded = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//genome version?
		if (versionedGenome == null) Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");

		//look for point directories
		if (bamFiles == null || bamFiles.length == 0) Misc.printExit("\nError: cannot find your bam files?\n");


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Bam 2 Bar : July 2011                             **\n" +
				"**************************************************************************************\n" +
				"Depreciated use Sam2USeq\n" +
				"\n" +
				"Generates per base read depth stair-step xxx.bar graph files for visualization in IGB\n" +
				"from BAM alignment files. By default, values are scaled per million mapped reads.\n\n" +

				"Options:\n"+
				"-b Full path to a BAM file or directory containing such.  These must be sorted by\n" +
				"      coordinate. Multiple files will be processed independently. \n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006, D_rerio_Jul_2010), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-r Don't scale graph values. Leave as actual read counts. \n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/Bam2Bar -b /Data/BamFiles/ -r\n\n"+

		"**************************************************************************************\n");

	}		

}
