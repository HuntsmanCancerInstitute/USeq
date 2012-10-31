package trans.main;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import trans.tpmap.MapSplitter;
import util.bio.parsers.*;
import util.gen.*;



/**
 * Loads arrays of Interval[] with oligo and sequence information needed for other applications.
 */
public class LoadIntervalOligoInfo {
	//fields
	private File[] intervalFiles;
	private Interval[] intervals;
	private int[] positions;
	private int[] matches;
	private float[] intensities;
	private File tpmapFile;  
	private File infoFile; 
	private File[] treatmentCelFiles;
	private File[] controlCelFiles;
	private File chromSeqDirectory;
	private String chromosomeSequence;
	private int chromeSeqLength;
	private MultiFastaParser fastaParser = new MultiFastaParser();
	
	
	public LoadIntervalOligoInfo(String[] args){
		//process user params
		processArgs(args);
		//get .tpmap file info
		ArrayList info = (ArrayList)IO.fetchObject(infoFile);
		//combine cel files
		int numTreatmentCelFiles = treatmentCelFiles.length;
		int numControlCelFiles = controlCelFiles.length;
		int numCelFiles = numTreatmentCelFiles+numControlCelFiles;
		File[] celFiles = new File[numCelFiles];	
		
		for (int i=0; i<numTreatmentCelFiles; i++) celFiles[i] = treatmentCelFiles[i];
		int index = 0;
		for (int i=numTreatmentCelFiles; i<numCelFiles; i++) celFiles[i] = controlCelFiles[index++];		

		String uniqueId;
		String celFilePathName;
		System.out.println("\nLaunching Loader ...");
		
		try {
			//for each Interval[] file
			for (int x=0; x<intervalFiles.length; x++){
				intervals = (Interval[])IO.fetchObject(intervalFiles[x]);
				//sort intervals
				sortIntervalsByChromosome(intervals);
				
				//for each cel file
				for (int i=0; i<numCelFiles; i++){
					String chromosome = "";
					//fetch float[] of normalized transformed intensity values 
					System.out.println("\tProcessing "+celFiles[i]);
					intensities = (float[])IO.fetchObject(celFiles[i]);
				
					//break both apart by chromosome using the tpmapInfo ArrayList, save to disk,
					celFilePathName = celFiles[i].getCanonicalPath();
					uniqueId = MapSplitter.breakSaveIntensityValues(info, intensities, celFilePathName);
					intensities = null;
					
					//for each interval
					int numIntervals = intervals.length;
					ArrayList intAL;
					int[] indexes;
					int numberOligos;
					int counter;
					Oligo[] oligos=null;
					
					for (int j=0; j<numIntervals; j++){
						//if first interval set # celFiles in each interval, sort of a hack to keep track of where intensity values came from
						if (i==0){
							intervals[j].setCelFiles(celFiles);
							intervals[j].setNumberTreatmentIntensities(numTreatmentCelFiles);
							intervals[j].setNumberControlIntensities(numControlCelFiles);
						}
						oligos = intervals[j].getOligos();
						
						//check if new chromosome, if so load intensities and positions and sequences
						if (chromosome.equals(intervals[j].getChromosome())==false){
							chromosome = intervals[j].getChromosome();
							intensities = (float[])IO.fetchObject(new File(celFilePathName+chromosome+uniqueId));
							//only load if first pass
							if (oligos==null){
								chromosomeSequence = null;
								positions = null;
								matches = null;
								fastaParser.resetFields();
								setGenomicSequenceParams(chromosome);
								positions = (int[])IO.fetchObject(new File(tpmapFile+chromosome));
								//any matches?
								File m = new File (tpmapFile+chromosome+"Matches");
								if (m.exists()){
									matches = (int[])IO.fetchObject(m);
								}
								else matches = new int[positions.length];
							}
						}
						//make oligos if they don't already exist
						if (oligos==null) {
							//find start index
							indexes = Num.findStartStopIndexes (positions, intervals[j].getStart1stOligo(), intervals[j].getStartLastOligo(), false);
							//find number of oligos
							numberOligos = indexes[1]-indexes[0]+1;							
							Oligo[] newOligos = new Oligo[numberOligos];
							counter =0;
							numberOligos++;
							for (int k= indexes[0]; k<=indexes[1]; k++){
								ArrayList al = new ArrayList();
								al.add(new Float(intensities[k]));
								int start = positions[k]-1;
								int stop = positions[k]+intervals[j].getSizeOfOligoMinusOne();
								if (stop>chromeSeqLength){
									IO.deleteFiles(celFiles[i].getParent(), uniqueId);
									Misc.printErrAndExit("\nError: one of your oligo coordinates exceeds the length of " + chromosome+
											". Check that you are using the exact same genome build in the tpmap and xxx.fasta files.\n");
								}
								String oligoSequence = new String(chromosomeSequence.substring(start, stop));
								newOligos[counter++] = new Oligo(k, positions[k], matches[k], al, oligoSequence);
							}						
							intervals[j].setOligos(newOligos);
						}
						//oligos exist thus just add new intensities
						else {
							numberOligos = oligos.length;
							for (int k=0; k<numberOligos; k++){							
								intAL = oligos[k].getIntensities();							
								intAL.add(new Float(intensities[oligos[k].getIndex()]));
							}
						}
						//if no sequence present, add sequence for interval
						//by some lucky fluke the Start and End return the exact sequence from 3.1 and 4.0, be sure to check
						//only holds when using the affy tpmap coordinates, when using to blast or mummer mapped coordinates subtract one 

						if (Misc.isEmpty(intervals[j].getSequence())){			
							String sequence = new String(chromosomeSequence.substring(intervals[j].getStart1stOligo()-1, 
									intervals[j].getSizeOfOligoMinusOne()+intervals[j].getStartLastOligo()-1));
							intervals[j].setSequence(sequence);
						}
					}
					//remove split intensity tmp files
					IO.deleteFiles(celFiles[i].getParent(), uniqueId);
				}
				//save intervals
				File save = new File (intervalFiles[x].getCanonicalPath()+"Ld");
				System.out.println("\tSaving oligo loaded Interval file "+save.getCanonicalPath());
				IO.saveObject(save, intervals);
			}
		} catch (Exception e){
			System.out.println("\nSomething is wrong with one of your cel or interval files");
			e.printStackTrace();
		}
		
	}
	
	/**Loads up the fasta file for a particular chromosome and sets some params.*/
	public void setGenomicSequenceParams(String chromosome){
		File chromFastaFile = new File (chromSeqDirectory,chromosome+".fasta");
		fastaParser.parseIt(chromFastaFile);
		chromosomeSequence = fastaParser.getSeqs()[0];
		chromeSeqLength = chromosomeSequence.length();
	}
	
	public static void main(String[] args) {
		if (args.length<5){
			printDocs();
			System.exit(0);
		}
		new LoadIntervalOligoInfo(args);
	}
	/**This method will process each argument and assign new fields*/
	public void processArgs(String[] args){
		File treatmentCelDirectory = null;
		File controlCelDirectory = null;
		File directory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatmentCelDirectory = new File (args[i+1]); i++; break;
					case 'c': controlCelDirectory = new File (args[i+1]); i++; break;
					case 'i': directory = new File (args[i+1]); i++; break;
					case 'b': tpmapFile = new File(args[i+1],"tpmap.fa");infoFile = new File (args[i+1],"tpmap.faInfo");i++; break;
					case 's': chromSeqDirectory = new File (args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		//get intervals
		if (directory == null || directory.exists() == false){
			System.out.println("\nCould not find your interval file or directory?!\n");
			System.exit(1);
		}
		intervalFiles = IO.extractFiles(directory);
		
		//check directories and files
		StringBuffer sb = new StringBuffer();
		if (treatmentCelDirectory == null || treatmentCelDirectory.isDirectory()==false) sb.append("Cannot find your treatment directory.\n");
		if (controlCelDirectory == null || controlCelDirectory.isDirectory()==false) sb.append("Cannot find your control directory.\n");
		if (tpmapFile == null || tpmapFile.exists() ==false) sb.append("Cannot find your TPMapFiles directory.\n");
		if (chromSeqDirectory == null || chromSeqDirectory.isDirectory() == false) sb.append("Cannot find your split genomic sequence file directory.\n");
		if (Misc.isNotEmpty(sb.toString())){
			System.out.println(sb);
			System.exit(1);
		}
		
		//extract processed cel files
		treatmentCelFiles = IO.extractFilesReturnFiles(treatmentCelDirectory, ".celp");
		controlCelFiles = IO.extractFilesReturnFiles(controlCelDirectory, ".celp");
		Arrays.sort(treatmentCelFiles);
		Arrays.sort(controlCelFiles);
		
		//check that cel files were found
		if (treatmentCelFiles.length==0 || controlCelFiles.length == 0){
			System.out.println("\nNo .celp files were found in your treatment or control directories'!\n");
			System.exit(1);
		}
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Load Interval Oligo Info: Dec 2005                          **\n" +
				"**************************************************************************************\n" +
				"LIOI fetches oligo intensities, start positions, and their sequence for a given\n" +
				"array of intervals using the following parameters:\n\n" +
				
				"-b The full path 'xxxTPMapFiles' directory text generated by the TPMapProcessor\n" +
				"-s The full path directory text containing the split genomic fasta files.\n" +
				"-t The full path directory text containing the serialized float[] treatment file(s)\n" +
				"-c The full path directory text containing the serialized float[] control file(s)\n" +
				"-i The full path Interval[] file text or directory containing Interval[] files.\n\n" +
				
				"Example: java -Xmx256M -jar pathTo/T2/Apps/LoadIntervalOligoInfo -b\n" +
				"      /affy/675bpTPMapFiles -s /seq/dmel/ -t /affy/tCels -c /affy/cCels \n\n" +
				
		"**************************************************************************************\n");		
	}	
	
	/**Ugly, ugly, ugly, use a Comparator!*/
	public static void sortIntervalsByChromosome(Interval[] intervals){
		//run thru and find text of all chromosomes
		HashSet chroms = new HashSet();
		int num = intervals.length;
		for (int i=0; i<num; i++){
			chroms.add(intervals[i].getChromosome());
		}
		//convert to array
		String[] chromNames = new String[chroms.size()];
		Iterator it = chroms.iterator();
		int counter = 0;
		while (it.hasNext()){
			chromNames[counter++] = (String)it.next();
		}
		//sort
		Arrays.sort(chromNames);
		//make hashMap
		int numChroms = chromNames.length;
		HashMap chromNumber = new HashMap();
		for (int i=0; i< numChroms; i++){
			chromNumber.put(chromNames[i], new Integer(i));
		}
		//run through and assign number
		for (int i=0; i<num; i++){
			Integer number = (Integer)chromNumber.get(intervals[i].getChromosome());
			intervals[i].setSortBy(number.intValue());
		}
		//sort intervals
		Arrays.sort(intervals);
	}
}
