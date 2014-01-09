package edu.utah.seq.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Misc;
import util.gen.Num;

public class MicrosatelliteCounter {
	private File mergedFastq = null;
	private File primerFile = null;
	private File outputDirectory = null;
	private String sampleName = "";
	private ArrayList<MsPrimer> primers = new ArrayList<MsPrimer>();
	private ArrayList<Pattern> pattern1 = new ArrayList<Pattern>();
	private ArrayList<Pattern> pattern2 = new ArrayList<Pattern>();
	private ArrayList<RepeatContainer> repContainers = new ArrayList<RepeatContainer>();
	private boolean bothPrimers = false;
	

	public MicrosatelliteCounter(String[] args) {
		this.parseArgs(args);
		this.run();
	}
	
	public void run() {
		//Read in primer infomration from reference file
		System.out.println("Parsing Primers");
		this.parsePrimers();
		
		//initialize count containers
		int[] primerCountBoth = new int[primers.size()];
		int[] primerCount1 = new int[primers.size()];
		int[] primerCount2 = new int[primers.size()];
		int noneCount = 0;
		int diffCount = 0;
		int totalCount = 0;
		
		//initialize counter
		int counter = 0;
		
	    //Read through merged fastq file
		System.out.println("Reading Fastq");
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.mergedFastq));
			
			String temp = null;
			
			
			while ((temp = br.readLine()) != null ) {
				counter++;
				//The sequence is found on the second line of the fastq file
				if (counter % 4 == 2) {
					totalCount++;
					//Store standard and reverse complement sequences
					String seq1 = temp;
					String seq2 = this.reverseComplement(temp);
					String usedSequence = ""; //The sequence in the 'forward' orientation with primers lopped off
					
					//Primer indexes for primer1 and primer2
					int index1 = -1; 
					int index2 = -1;
					int usedIndex = -1;
					int tempIndex = -1;
					
					
					//Orientation of primer
					boolean reversed = false;
					
					//Run through all combinations of primer / sequence orientiation to determine the origin of the sequence
					//and whether it should be flipped.
				    tempIndex = this.findPrimersInSequence(seq1, pattern1); //Check for primer1 at start of the sequence: standard orientation
				   
				    if (tempIndex != -1) {
				    	index1 = tempIndex;
				    
				    } else {
				    	tempIndex = this.findPrimersInSequence(seq2, pattern1); //Check for primer1 at the start of the rc sequence: reverse orient
				    	
				    	if (tempIndex != -1) {
				    		index1 = tempIndex;
				    		reversed = true;
				    	}
				    }
				    
				    tempIndex = this.findPrimersInSequence(seq1,pattern2); //Check for primer2 at the start of the sequence: reverse orient
				    
				    if (tempIndex != -1) {
				    	index2 = tempIndex;
				    	reversed = true;
				    } else {
				    	tempIndex = this.findPrimersInSequence(seq2,pattern2); //Check for primer2 at the start of the rc sequence: standard
				    	
				    	if (tempIndex != -1) {
				    		index2 = tempIndex;
				    	}
				    }
				    
				    //grab the sequence based on orientation
				    if (reversed) {
				    	usedSequence = seq2;
				    } else {
				    	usedSequence = seq1;
				    }
				    
				    //Check for inconsistencies and determine sequence origination.
				    if (index1 == -1) {
				    	if (index2 == -1) {
				    		noneCount += 1;
				    	} else {
				    		primerCount2[index2] += 1;
				    		usedIndex = index2;
				    		usedSequence = usedSequence.substring(0,usedSequence.length()-primers.get(usedIndex).getPrimer2().length());
				    	}
				    } else if (index2 == -1) {
				    	primerCount1[index1] += 1;
				    	usedIndex = index1;
				    	usedSequence = usedSequence.substring(primers.get(usedIndex).getPrimer1().length());
				    	
				    } else {
				    	if (index1 == index2) {
				    		primerCountBoth[index1] += 1;
				    		usedIndex = index1;
				    		usedSequence = usedSequence.substring(primers.get(usedIndex).getPrimer1().length(),usedSequence.length()-primers.get(usedIndex).getPrimer2().length());
				    	} else {
				    		diffCount += 1;
				    	}
				    }
				    
				    //If both primers are required, bail if either have -1 index
				    if (this.bothPrimers) {
				    	if (index1 == -1 || index2 == -1 || usedIndex == -1) {
				    		continue;
				    	}
				    } else if (usedIndex == -1) { //If only one primer is required, only bail if usedIndex is -1
				    	continue;
				    }
				    
				    //Store information about repeat
				    int maxPos = -1;
				    int maxCount = 0;
				    
				    String repeatToCheck = primers.get(usedIndex).getRepeat();
				    
				    for (int i=0; i<usedSequence.length()-repeatToCheck.length();i++) {
				    	//Create variable
				    	int tempPos = i;
				    	int tempCount = 0;
				    	int repeatLength = repeatToCheck.length();
				    	
				    	//Stop checking if a longer repeat can't be generated
				    	int remainingBases = usedSequence.length() - tempPos;
				    	int maxRemainingRepeats = remainingBases / repeatLength;
				    	if (maxRemainingRepeats < maxCount) {
				    		break;
				    	}
				    	
				    	//Check if repeat matches expected
				    	String repeat = usedSequence.substring(tempPos,tempPos+repeatLength);
				    	if (!repeat.equals(repeatToCheck)) {
				    		continue;
				    	}
				    	tempCount++; //If it matches, start counting at 1
				    	
				    	//Check remaining sequence for additional repeat units
				    	int jumps = (usedSequence.length() - tempPos) / repeatLength; //How many possible repeats remain
				    	if (jumps > 1) {
				    		for (int j=1;j<jumps;j++) { //Go through each potential repeat 
				    			int tempPos2 = tempPos + (j * repeatLength); //Calculate potential repeat position
				    			repeat = usedSequence.substring(tempPos2,tempPos2+repeatLength);
				    			
				    			if (repeat.equals(repeatToCheck)) {
						    		tempCount++;
						    	} else {
						    		break;
						    	}
				    		}
				    	}
				    	
				    	//If current repeat run is the longest, store
				    	if (tempCount > maxCount) {
				    		maxCount = tempCount;
				    		maxPos = tempPos;
				    	}
				    	
				    	//Set the sequence position to end of repeat
				    	i = tempPos + (maxCount * repeatLength);
				    }
				    	
				    //If primer1 is missing and the repeat starts at 0, don't record
				    if (maxPos == 0 && index1 == -1) {
				    	repContainers.get(usedIndex).foundAmbiguous();
				    }
				    //If primer2 is missing and the reapeat ends at the last position, don't record
				    else if (maxCount*repeatToCheck.length()+repeatToCheck.length()-1+maxPos == usedSequence.length() && index2 == -1) {
				    	repContainers.get(usedIndex).foundAmbiguous();
				    	
				    } 
				    //If the repeat couldn't be found, skip
				    else if (maxCount == 0) {
				    	repContainers.get(usedIndex).foundNothing();
		
				    } 
				    //Write the repeat
				    else {
				    	repContainers.get(usedIndex).addSeq(maxPos, maxCount);
				    }
				}
				
			}
			
			br.close();
			
			System.out.println("Write repeat data to file");
			this.printRepeatFile();
			System.out.println("Write primer data to file");
			this.printPrimerFile(primerCountBoth, primerCount1, primerCount2, noneCount, diffCount, totalCount);
			System.out.println("Finished");
			
		} catch (FileNotFoundException fnfe) {
			System.out.println("Could not find the merged fastq file: " + this.mergedFastq.getName());
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Error reading merged fastq file: " + this.mergedFastq.getName());
			System.exit(1);
		}
		
	}
	
	public void printRepeatFile() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(this.outputDirectory + "/" + this.sampleName + ".repeat.txt"));

			for (int i=0; i<repContainers.size();i++) {
				bw.write(String.format("%s\n",primers.get(i).getName()));
				bw.write(repContainers.get(i).getRepeatOutput());
				bw.write("\n********************\n\n");
			}
			
			bw.close();
		} catch (FileNotFoundException fnfe) {
			System.out.println("Could not find the repeat output file: " + this.sampleName + ".repeat.txt");
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Error writing to the repeat output file: " + this.sampleName + ".repeat.txt");
			System.exit(1);
		}
	}
	
	public void printPrimerFile(int[] primerCountBoth, int[] primerCount1, int[] primerCount2, int noneCount, int diffCount, int totalCount) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(this.outputDirectory + "/" + this.sampleName + ".primer.txt"));

			bw.write(String.format("Total Reads: %d\n",totalCount));
			bw.write(String.format("No primers identified: %d (%.2f%%)\n",noneCount, (float)noneCount / totalCount * 100));
			bw.write(String.format("Identified primers didn't match: %d (%.2f%%)\n\n", diffCount, (float)diffCount / totalCount * 100));
			
			int sumPrimer1 = Num.sumIntArray(primerCount1);
			bw.write(String.format("Primer 1 only %d (%.2f%%)\n", sumPrimer1, (float)sumPrimer1 / totalCount * 100 ));
			for (int i=0; i<primers.size();i++) {
				bw.write(String.format("%s %d (%.2f%%)\n", primers.get(i).getName(), primerCount1[i], (float)primerCount1[i] / sumPrimer1 * 100));
			}
			bw.write("\n");
			
			int sumPrimer2 = Num.sumIntArray(primerCount2);
			bw.write(String.format("Primer 2 only %d (%.2f%%)\n", sumPrimer2, (float)sumPrimer2 / totalCount * 100 ));
			for (int i=0; i<primers.size();i++) {
				bw.write(String.format("%s %d (%.2f%%)\n", primers.get(i).getName(), primerCount2[i], (float)primerCount2[i] / sumPrimer2 * 100));
			}
			bw.write("\n");
			
			int sumBothPrimer = Num.sumIntArray(primerCountBoth);
			bw.write(String.format("Both primers found %d (%.2f%%)\n", sumBothPrimer, (float)sumBothPrimer / totalCount * 100 ));
			for (int i=0; i<primers.size();i++) {
				bw.write(String.format("%s %d (%.2f%%)\n", primers.get(i).getName(), primerCountBoth[i], (float)primerCountBoth[i] / sumBothPrimer * 100));
			}
			bw.write("\n");
			
			bw.close();
		} catch (FileNotFoundException fnfe) {
			System.out.println("Could not find the primer output file: " + this.sampleName + ".primer.txt");
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Error writing to the primer output file: " + this.sampleName + ".primer.txt");
			System.exit(1);
		}
	}
	
	public int findPrimersInSequence(String sequence, ArrayList<Pattern> patternSet) {
		int index = -1;
		for (int i=0; i<this.primers.size();i++) {
			Pattern pattern = patternSet.get(i);
			Matcher m = pattern.matcher(sequence);
			if (m.matches()) {
				index = i;
				break;
			}	
		}
		return index;
	}
	
	
	/* This method parses a primer reference file. The information is stored in an ArrayList of MsPrimer.  This method also creates patterns for
	 * primer1 and primer2 for use in matching.
	 */
	public void parsePrimers() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.primerFile));
			String temp = null;
			
			while ((temp = br.readLine()) != null) {
				String[] entries = temp.split("\t");
				if (entries.length != 4) {
					System.out.println("There aren't enough entries in the primer reference file. The format should be NAME\tPRIMER1\tPRIMER2\tREPEAT\n");
					System.exit(1);
				}
				
				System.out.println("Found: " + entries[0]);
				
				MsPrimer newPrimer = new MsPrimer(entries[0],entries[1],entries[2],entries[3]);
				primers.add(newPrimer);
				
				Pattern p1 = Pattern.compile("^" + entries[1] + ".+");
				pattern1.add(p1);
				
				Pattern p2 = Pattern.compile("^" + entries[2] + ".+");
				pattern2.add(p2);
				
				RepeatContainer newRC = new RepeatContainer(entries[3], entries[0], this.sampleName, this.outputDirectory.toString());
				repContainers.add(newRC);
			}
			
			br.close();
		} catch (FileNotFoundException fnfe) {
			System.out.println("Could not find primer file, exiting: " + this.primerFile.getName());
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Error reading primer file, exiting: " + this.primerFile.getName()); 
		}
	}
	
	public String reverseComplement(String sequence) {
		HashMap<Character,Character> rcHash = new HashMap<Character,Character>();
		rcHash.put('A','T');
		rcHash.put('T','A');
		rcHash.put('C','G');
		rcHash.put('G','C');
		rcHash.put('N','N');
		
		StringBuffer compSeq = new StringBuffer("");
		
		for (int i=0; i < sequence.length();i++) {
			compSeq.append(rcHash.get(sequence.charAt(i)));
		}
		
		String revCompSeq = compSeq.reverse().toString();
		
		return revCompSeq;
	}

	public static void main(String[] args) {
		new MicrosatelliteCounter(args);
	}
	
	private void parseArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': mergedFastq = new File(args[++i]); break;
					case 'p': primerFile = new File(args[++i]); break;
					case 'n': sampleName = args[++i]; break;
					case 'd': outputDirectory = new File(args[++i]); break;
					case 'b': bothPrimers = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		
		if (mergedFastq == null) {
			exitMessage("Merged fastq file not specified, exiting.");
		}
		
		if (primerFile == null) {
			exitMessage("Primer reference file not specified, exiting.");
		}
		
		if (outputDirectory == null) {
			exitMessage("Output directory not specified, exiting.");
		}
		
		if (!outputDirectory.exists()) {
			exitMessage("Ouput directory does not exist, exiting.");
		}
	 
		if (!outputDirectory.isDirectory()) {
			exitMessage("Specified output directory is not a directory, exiting.");
		}
		
	}
	
	private class RepeatContainer {
		private String repeat = null; //Repeat
		private ArrayList<Integer> repeatPositionList = null; //List of observed repeat start positions
		private ArrayList<Integer> repeatCountList = null; //List of observed repeat counts
		private int sequenceCount = 0;
		private int nothing = 0;
		private int ambiguous = 0;
		private String primer = null;
		private String sampleName = null;
		private String directoryName = null;

		
		/* This class handles sample-wide repeat information */
		public RepeatContainer(String repeat, String primer, String sampleName, String directoryName) {
			this.repeat = repeat;
			this.primer = primer;
			this.directoryName = directoryName;
			this.sampleName = sampleName;
			repeatPositionList = new ArrayList<Integer>();
			repeatCountList = new ArrayList<Integer>();
		}
		
		public void foundNothing() {
			this.nothing++;
		}
		
		public void foundAmbiguous() {
			this.ambiguous++;
		}
		
		/* This method adds repeat information for a given sequence */
		public void addSeq(int position, int count) {
			this.repeatPositionList.add(position);
			this.repeatCountList.add(count);
			this.sequenceCount += 1;
		}
		
		/* This method prints out repeat information.  Currently text-only */
		public String getRepeatOutput() {
			if (this.repeatCountList.size() == 0) {
				return "";
			} else {
				String repeatData = this.createHistogram(this.repeatCountList, "NumRepeatUnits", "Observations", "repeat");
				String posData = this.createHistogram(this.repeatPositionList, "PositionInSequence", "Observations","position");
				
				
				return String.format("%s\t%s\t%s\t%s\n\n%s\n\n%s\n", this.repeat,
						Integer.toString(this.sequenceCount), Integer.toString(this.nothing), Integer.toString(this.ambiguous),repeatData,
						posData);
			}
		}
		
		/* This method takes a list of observed repeat counts, or repeat positions and makes a histogram.*/
		private String createHistogram(ArrayList<Integer> rawCounts, String xaxis, String yaxis, String type) {
			//Determine min/max values
			int min = Integer.MAX_VALUE;
			int max = Integer.MIN_VALUE;
			for (Integer i: rawCounts) {
				if (i < min) {
					min = i;
				} else if (i > max) {
					max = i;
				}
			}
			
			//Create and fill observation array
			int[] histArray = new int[max+1];
			
		    for (Integer i: rawCounts) {
		    	histArray[i] += 1;
		    }
		    String histString = "";
		    
		    for (int i=min; i<max+1; i++) {
		    	histString += "\t" + Integer.toString(histArray[i]);
		    }
		    histString = histString.substring(1);
		    histArray = Arrays.copyOfRange(histArray, min, max+1);
		    
		    
		    //Create and fill count array
		    String numberString = "";
		    int[] intArray = new int[max+1-min];
		    int p=0;
		    for (int i=min;i<max+1;i++) {
		    	numberString += "\t" + Integer.toString(i);
		    	intArray[p] = i;
		    	p++;
		    }
		    numberString = numberString.substring(1);
		    
		    //Create histogram
		    this.createRhistogram(histArray, intArray, xaxis, yaxis,type);
		    
		    //Return string
		    return String.format(xaxis + ":\t%s\n" + yaxis + ":\t%s\n", numberString, histString);
		   
		}
		
		
		
		private void createRhistogram(int[] observations, int[] position, String xaxis, String yaxis, String type) {
			String filename = this.directoryName + "/" + this.sampleName + "." + this.primer + "." + type;
			File script = new File(filename + ".R");
			File data = new File(filename + ".txt");
			File image = new File(filename + ".png");
			
			try {
				
				//Write data file
				BufferedWriter bw = new BufferedWriter(new FileWriter(data));
				bw.write(xaxis + "\t" + yaxis + "\n");
				for (int i=0; i<observations.length; i++) {
					bw.write(position[i] + "\t" + observations[i] + "\n");
				}
				bw.close();
				
				//Write Rscript
				BufferedWriter bw2 = new BufferedWriter(new FileWriter(script));
				bw2.write("library(ggplot2)\n");
				bw2.write("d <- read.table(\"" + data.toString() + "\",header=TRUE)\n");
				bw2.write("ggplot(d,aes(" + xaxis + "," + yaxis + ")) + geom_bar(stat=\"identity\") + ggtitle(\"" + this.sampleName + " " + this.primer + 
						": " + xaxis + " vs " + yaxis + "\")\n");
				bw2.write("ggsave(\"" + image.getAbsolutePath().toString() + "\")\n");
				bw2.close();
				
				//Run Script
				try {
					ProcessBuilder pb = new ProcessBuilder("Rscript",script.getAbsoluteFile().toString());
					Process p = pb.start();
					
					int val = p.waitFor();
					
					if (val != 0) {
						System.out.println("Error running R script");
						BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
						String line2 = null;
						while((line2 = br2.readLine()) != null) {
							System.out.println(line2);
						}
						System.exit(1);
					}
					
				} catch (IOException ioex) {
					System.out.println("IO Exception while trying run RScript.");
					System.exit(1);
				} catch (InterruptedException ieex) {
					System.out.println("Process was interrupted while trying to run RScript");
					System.exit(1);
				}
				
				//Delete input files
				script.delete();
				data.delete();
				
			} catch (FileNotFoundException fnfe) {
				System.out.println("Could not find file: " + data.toString());
				System.exit(1);
			} catch (IOException ioex) {
				System.out.println("Error writing to file: " + data.toString());
				System.exit(1);
			}
		}
		
	}
	
	private class MsPrimer {
		private String name;
		private String primer1;
		private String primer2;
		private String repeat;
		
		public MsPrimer(String name, String primer1, String primer2, String repeat) {
			this.name = name;
			this.primer1 = primer1;
			this.primer2 = primer2;
			this.repeat = repeat;
		}
		
		public String getName() {
			return this.name;
		}
		
		public String getPrimer1() {
			return this.primer1;
		}
		
		public String getPrimer2() {
			return this.primer2;
		}
		
		public String getRepeat() {
			return this.repeat;
		}
	}
	
	private void exitMessage(String message) {
		printDocs();
		System.out.println(message);
		System.exit(0);
	}
	
	private static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Microsatellite Counter: Jan 2014                           **\n" +
				"**************************************************************************************\n" +
				"MicrosatelliteCounter identifies and counts microsatellite repeats in MiSeq fastq \n" +
				"files. This iteration of the software requires you to specify the primers used in the\n" +
				"sequencing project.  It will automatically find the most likely microsatellite by \n" +
				"looking at all possible repeats of length 1 through length 10 and finding the longest\n" + 
				"repeat by length, not repeat unit.  There are two output files generated, the first \n" + 
				"lists primer statistics (currently only reads with both primers are used), the \n" +
				"second lists repeat data.  Note that the input file are fastq sequence that were \n" +
				"merged using a program like PEAR\n\n" +
				
				"\nRequired Arguments:\n\n"+
				"-f Merged fastq file. Path to merged fastq file. We currently suggest using PEAR to \n" +
				"       merge fastq sequences.\n" +
				"-p Primer file.  Path to primer reference file.  This file lists each primer used in \n" +
				"       in the sequencing project in the format NAME<tab>PRIMER1<tab>PRIMER2.\n" +
				"-n Sample name.  Sample name.  This string will be appended to the output files names.\n" +
				"-d Directory. Output directory. Output files will be written to this directory\n" +
				"-b Require both primers.  Both primers must be identified in order to more forward \n" +
				"       with the analysis.\n" +
				
				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MicrosatelliteCounter -f Merged.fastq \n" +
				"      -r PrimerReference.txt -p 10511X1.primer.txt -o 10511X1.repeat.txt\n\n" +

				"**************************************************************************************\n");

	}
	
	
	
	

}
