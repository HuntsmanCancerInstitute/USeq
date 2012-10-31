package trans.main;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;
import util.bio.seq.*;

/**
 * Loads arrays of Interval[] with oligo and sequence information needed for other applications.
 * Uses chromosome specific files.
 */
public class LoadChipSetIntervalOligoInfo {
	//fields
	private File[] treatmentDirectories = null;
	private File[] controlDirectories = null;
	private File[] chromOligoPositions = null;
	private File oligoPositionDirectory = null;

	private File[] intervalFiles;
	private Interval[] intervals;
	private int[] positions;
	private float[] intensities; 
	private boolean loadSequence = false;
	private File indexedSeqDirectory;
	private FetchIndexedSequence seqFetcher;

	private boolean notExact = false; //forces exact matching of oligo positions, keep false unless your using made up intervals.

	public LoadChipSetIntervalOligoInfo(String[] args){

		//process user params
		processArgs(args);

		//combine cel files
		int numTreatments = treatmentDirectories.length;
		int numControls = controlDirectories.length;
		int totalDirectories = numTreatments+numControls;
		File[] intensityDirectories = new File[totalDirectories];	
		
		for (int i=0; i<numTreatments; i++) intensityDirectories[i] = treatmentDirectories[i];
		int index = 0;
		for (int i=numTreatments; i<totalDirectories; i++) intensityDirectories[i] = controlDirectories[index++];		
		
		//instantiate the seq fetcher
		if (loadSequence) seqFetcher = new FetchIndexedSequence (indexedSeqDirectory,false);
		
		System.out.println("\nLaunching Loader ...");
		
		try {
			//for each Interval[] file
			for (int x=0; x<intervalFiles.length; x++){
				intervals = (Interval[])IO.fetchObject(intervalFiles[x]);
				System.out.println("Processing interval file "+intervalFiles[x]);				
				//sort intervals by chromosome
				Arrays.sort(intervals, new IntervalComparator());
				String chromosome = "";
				
				//for each cel file directory
				for (int i=0; i<totalDirectories; i++){
					System.out.println("\nProcessing "+intensityDirectories[i]);
					boolean intensitiesLoaded = false;
					
					//for each interval
					for (int j=0; j<intervals.length; j++){
						Oligo[] oligos = intervals[j].getOligos();
						
						//check if new chromosome
						if (chromosome.equals(intervals[j].getChromosome())==false){
							chromosome = intervals[j].getChromosome();
							//only load if first pass thru interval array
							if (oligos==null) {
								loadOligoPositions(chromosome);
								//set indexedSequence[]
								if (loadSequence) seqFetcher.setIndexedSequence(chromosome);
							}
							//throw intensity loaded flag
							intensitiesLoaded = false;

						}
						
						//load intensities?
						if (intensitiesLoaded == false){
							intensities = loadIntensities (intensityDirectories[i], chromosome);
							intensitiesLoaded = true;
						}
						
						//if first, set celFile directories in each interval, sort of a hack to keep track of where intensity values came from
						if (i==0){
							intervals[j].setCelFiles(intensityDirectories);
							intervals[j].setNumberTreatmentIntensities(numTreatments);
							intervals[j].setNumberControlIntensities(numControls);
						}

						//make oligos if they don't already exist
						if (oligos==null) {
							//find start index
							int[] indexes = Num.findStartStopIndexes (positions, intervals[j].getStart1stOligo(), intervals[j].getStartLastOligo(), notExact);
							//find number of oligos
							int numberOligos = indexes[1]-indexes[0]+1;
							Oligo[] newOligos = new Oligo[numberOligos];
							//make array to hold start stop positions for fetching oligo sequence
							int[][] startStop = new int[numberOligos][2];
							int counter =0;
							for (int k= indexes[0]; k<=indexes[1]; k++){
								ArrayList al = new ArrayList();
								al.add(new Float(intensities[k]));
								//original int start = positions[k]-1;
								//mod 
								int start = positions[k];
								if (start == -1) System.out.println("Start = -1");
								int stop = positions[k]+intervals[j].getSizeOfOligoMinusOne();
								startStop[counter] = new int[] {start, stop};
								newOligos[counter++] = new Oligo(k, positions[k], 0, al, null);
							}	
							//fetch and load oligo sequences
							if (loadSequence) {
								String[] seqs = seqFetcher.fetchSequences(startStop);
								for (int k=0; k< numberOligos; k++) {
									newOligos[k].setSequence(seqs[k]);
								}
							}
							//set in interval
							intervals[j].setOligos(newOligos);
						}
						//oligos exist thus just add new intensities
						else {
							for (int k=0; k<oligos.length; k++){							
								ArrayList intAL = oligos[k].getIntensities();
								intAL.add(new Float(intensities[oligos[k].getIndex()]));
							}
						}
						//if no sequence present, add sequence for interval 
						if (loadSequence && Misc.isEmpty(intervals[j].getSequence())){	
							int[][] startStop = {{intervals[j].getStart1stOligo(), intervals[j].getSizeOfOligoMinusOne()+intervals[j].getStartLastOligo()-1}};
							//String sequence = new String(chromSequence.substring(intervals[j].getStart1stOligo(), intervals[j].getSizeOfOligoMinusOne()+intervals[j].getStartLastOligo()-1));
							intervals[j].setSequence(seqFetcher.fetchSequences(startStop)[0]);
						}
					}
				}
				//save intervals
				File save = new File (intervalFiles[x].getCanonicalPath()+"Ld");
				System.out.println("\tSaving loaded Interval file "+save.getCanonicalPath());
				IO.saveObject(save, intervals);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		if (args.length<5){
			printDocs();
			System.exit(0);
		}
		new LoadChipSetIntervalOligoInfo(args);
	}
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		//comma delimited list
		String treatment = null;
		String control = null;
		File intervalDirectory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatment = args[i+1]; i++; break;
					case 'c': control = args[i+1]; i++; break;
					case 'i': intervalDirectory = new File (args[i+1]); i++; break;
					case 's': indexedSeqDirectory = new File (args[i+1]); i++; loadSequence= true; break;
					case 'o': oligoPositionDirectory = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test);
				}
			}
		}
		
		//look for required parameters
		if (treatment == null || control == null || oligoPositionDirectory == null || intervalDirectory == null){
			Misc.printExit("\nPlease complete one or more of the following required parameters: -t, -c, -o, -s, or -i .\n");
		}
		//get intervals
		intervalFiles = IO.extractFiles(intervalDirectory);
		if (intervalFiles.length == 0) Misc.printExit("\nProblem extracting intervals to process -> "+intervalDirectory);

		//parse treatments
		treatmentDirectories = IO.extractFiles(treatment);
		if (treatmentDirectories == null) Misc.printExit("\nProblem parsing treatment directories -> "+treatment);
		
		//parse control
		controlDirectories = IO.extractFiles(control);
		
		//parse oligo positions
		chromOligoPositions = IO.extractFiles(oligoPositionDirectory);
		if (chromOligoPositions == null) Misc.printExit("\nProblem parsing oligo positions -> "+oligoPositionDirectory);
		
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                 Load Chip Set Interval Oligo Info : April 2007                    **\n" +
				"**************************************************************************************\n" +
				"LCSIOI fetches oligo intensities, start positions, and their sequence for a given\n" +
				"array of intervals using the following parameters:\n\n" +
				
				"-o The 'OligoPositions' directory, full path, generated by ChipSetCelProcessor.\n" +
				"-s Load sequence information, provide the full path directory text containing indexed\n" +
				"       genomic fasta sequences. Run the IndexFastas app on your chromosome split seqs.\n" +
				"-t Treatment chip set directories, full path, comma delimited, no spaces.\n" +
				"-c Control chip set directories, full path, comma delimited, no spaces.\n" +
				"-i The full path Interval[] file text or directory containing Interval[] files.\n\n" +
				
				"Example: java -Xmx256M -jar pathTo/T2/Apps/LoadChipSetIntervalOligoInfo -b\n" +
				"      /affy/675bpTPMapFiles -s /seq/dmel/ -t /affy/tCels -c /affy/cCels \n\n" +
				
		"**************************************************************************************\n");		
	}	
	
	/**Loads the oligo positions.*/
	public void loadOligoPositions(String chrom){
		File pos = new File(oligoPositionDirectory, chrom);
		if (pos.exists() == false) Misc.printExit("\nCannot find the oligo positions for chromosome "+chrom+" in "+oligoPositionDirectory);
		positions = (int[])IO.fetchObject(pos);
	}
	
	/**Loads the oligo intensities.*/
	public static float[] loadIntensities(File directory, String chrom){
		File ints = new File(directory, chrom);			
		if (ints.exists() == false) Misc.printExit("\nCannot find the intensities for chromosome "+chrom+" in "+directory);
		return (float[])IO.fetchObject(ints);
	}

}
