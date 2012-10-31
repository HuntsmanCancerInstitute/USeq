package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.data.*;
import util.gen.*;
import java.util.*;

import util.bio.annotation.*;
import util.bio.seq.Seq;

/**
 * Use for filtering identical alignments.
 * @author david.nix@hci.utah.edu 
 **/
public class FilterDuplicateAlignments {
	//fields
	private File[] dataFiles;
	private int chromosomeColumnIndex = -1;
	private int positionColumnIndex = -1;
	private int sequenceColumnIndex = -1;
	private int qualityColumnIndex = -1;
	private int strandColumnIndex = -1;
	private File saveDirectory;
	private File workingFile;
	private File tempDirectory;
	private int maxReadCount;
	private int totalNumberAlignments = 0;
	private int totalNumberFilteredAlignments = 0;
	private HashMap <String, DataOutputStream> chromOut = new HashMap <String, DataOutputStream>();
	private Pattern tab = Pattern.compile("\\t");
	private ArrayList<File> tempChrData = new ArrayList<File>();
	private boolean saveUniques = true;
	private boolean skipSpliceJunctions = true;
	private Pattern spliceJunction = Pattern.compile(".+_\\d+_\\d+");
	private Pattern negativeStrand = Pattern.compile("[-Rr]");

	//constructors
	public FilterDuplicateAlignments(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		//for each file, parse by chromosome and save to disk
		System.out.println("Splitting by chromosome...");
		parseFiles();

		//close the writers
		closeWriters();

		//load, sort, make point chromPointData, and save
		System.out.print("Filtering split data");
		filter();
		System.out.println();

		//cleanup
		IO.deleteDirectory(tempDirectory);

		//stats
		System.out.println("Stats:");
		System.out.println(totalNumberAlignments+"\tTotal number alignments");
		System.out.println(totalNumberFilteredAlignments+"\tTotal number filtered alignments");
		double fract = (double)totalNumberFilteredAlignments/(double)totalNumberAlignments;
		System.out.println(Num.formatNumber(fract, 3)+"\tFraction passing filters");

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void parseFiles(){
		for (int i=0; i< dataFiles.length; i++){
			//set working objects and parse tag file text
			workingFile = dataFiles[i];
			System.out.print("\t"+workingFile);

			//split file to chromosome strand specific temp files
			boolean parsed = parseWorkingFile(); 
			if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
			System.out.println();
		}
	}

	/**Closes writers.*/
	public void closeWriters(){
		try{
			Iterator<DataOutputStream> it = chromOut.values().iterator();
			while (it.hasNext()) it.next().close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/**Filters the split chromosome file by the max number of reads.*/
	public void filter(){
		try {
			File filteredAlignments = new File (saveDirectory, "filteredAlignments.txt");
			PrintWriter out = new PrintWriter ( new FileWriter( filteredAlignments));
			//for each composite data file
			for (int i=0; i< tempChrData.size(); i++){
				File chromDataFile = tempChrData.get(i);
				System.out.print(".");
				//get data and sort
				PositionLine[] pl = PositionLine.loadBinary(chromDataFile);
				Arrays.sort(pl);
				totalNumberAlignments += pl.length;
				//count and subsample
				pl = countAndSubSample(pl);
				totalNumberFilteredAlignments += pl.length;
				//write to file
				for (int j=0; j< pl.length; j++) out.println(pl[j].getLine());
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		System.out.println();
	}
	
	/**Counts the number of identical alignments and finds the best one based on the sum of the quality values.*/
	public PositionLine[] countAndSubSample(PositionLine[] sortedPL){
		Random random = new Random();
		ArrayList<PositionLine> al = new ArrayList<PositionLine>();
		ArrayList<PositionLine> block = new ArrayList<PositionLine>();
		block.add(sortedPL[0]);
		int currentPosition =sortedPL[0].getPosition();
		for (int i=1; i< sortedPL.length; i++){
			int testPosition = sortedPL[i].getPosition();
			if (currentPosition == testPosition) block.add(sortedPL[i]);
			else {
				int numReads = block.size();
				//just one read, then save
				if (numReads == 1) al.addAll(block);
				//save best?
				else if (qualityColumnIndex != -1){
					String bestLine = block.get(0).getLine();
					int bestScore = Num.sumIntArray(Seq.convertScores( tab.split(bestLine)[qualityColumnIndex]));
					for (int x=1; x< numReads; x++){
						String line = block.get(x).getLine();
						int totalScore = Num.sumIntArray(Seq.convertScores( tab.split(line)[qualityColumnIndex]));
						if (totalScore > bestScore) {
							bestLine = line;
							bestScore = totalScore;
						}
					}
					al.add(new PositionLine(currentPosition,bestLine));
				}
				//saving uniques?
				else if (saveUniques){					
					//make a hash map, load with sequence
					HashMap<String,ArrayList<String>> hash = new HashMap<String, ArrayList<String>>();
					for (int x=0; x< numReads; x++){
						String line = block.get(x).getLine();
						String seq = tab.split(line)[sequenceColumnIndex];
						if (hash.containsKey(seq)){
							hash.get(seq).add(line);
						}
						else {
							ArrayList<String> alx = new ArrayList<String>();
							alx.add(line);
							hash.put(seq, alx);
						}
					}
					//iterate through the hash randomly selecting one of each sequence
					Iterator<String> it = hash.keySet().iterator();
					while (it.hasNext()){
						ArrayList<String> lines = hash.get(it.next());
						int rand = random.nextInt(lines.size());
						PositionLine pl = new PositionLine(currentPosition,lines.get(rand));
						al.add(pl);
					}
				}
				//saving set number?
				else if (numReads > maxReadCount){
					//subsample
					PositionLine[] pl = new PositionLine[block.size()];
					block.toArray(pl);
					Misc.randomize(pl, 0);
					for (int j=0; j< maxReadCount; j++) al.add(pl[j]);
				}
				else al.addAll(block);
				//reset
				block.clear();
				block.add(sortedPL[i]);
				currentPosition = testPosition;
			}
		}

		//add last block
		int numReads = block.size();
		//just one read, then save
		if (numReads == 1) al.addAll(block);
		//save best?
		else if (qualityColumnIndex != -1){
			String bestLine = block.get(0).getLine();
			int bestScore = Num.sumIntArray(Seq.convertScores( tab.split(bestLine)[qualityColumnIndex]));
			for (int x=1; x< numReads; x++){
				String line = block.get(x).getLine();
				int totalScore = Num.sumIntArray(Seq.convertScores( tab.split(line)[qualityColumnIndex]));
				if (totalScore > bestScore) {
					bestLine = line;
					bestScore = totalScore;
				}
			}
			al.add(new PositionLine(currentPosition,bestLine));
		}
		//saving uniques?
		else if (saveUniques){					
			//make a hash map, load with sequence
			HashMap<String,ArrayList<String>> hash = new HashMap<String, ArrayList<String>>();
			for (int x=0; x< numReads; x++){
				String line = block.get(x).getLine();
				String seq = tab.split(line)[sequenceColumnIndex];
				if (hash.containsKey(seq)){
					hash.get(seq).add(line);
				}
				else {
					ArrayList<String> alx = new ArrayList<String>();
					alx.add(line);
					hash.put(seq, alx);
				}
			}
			//iterate through the hash randomly selecting one of each sequence
			Iterator<String> it = hash.keySet().iterator();
			while (it.hasNext()){
				ArrayList<String> lines = hash.get(it.next());
				int rand = random.nextInt(lines.size());
				PositionLine pl = new PositionLine(currentPosition,lines.get(rand));
				al.add(pl);
			}
		}
		//saving set number?
		else if (numReads > maxReadCount){
			//subsample
			PositionLine[] pl = new PositionLine[block.size()];
			block.toArray(pl);
			Misc.randomize(pl, 0);
			for (int j=0; j< maxReadCount; j++) al.add(pl[j]);
		}
		else al.addAll(block);

		//return final
		sortedPL = new PositionLine[al.size()];
		al.toArray(sortedPL);
		return sortedPL;
	}




	/**Splits a file by chromosome.*/
	public boolean parseWorkingFile(){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			String line;
			String[] tokens = null;
			int counter = 0;
			String currentChrom = "";
			DataOutputStream dos = null;
			while ((line = in.readLine()) !=null){
				try {
					if (line.startsWith("#")) continue;
					if (line.contains("chrAdapter")) continue;
					tokens = tab.split(line);

					//parse chromosome
					String chromosome = tokens[chromosomeColumnIndex];
					//check for splice junction
					if (skipSpliceJunctions && spliceJunction.matcher(chromosome).matches()) continue;
					//parse strand
					String strand = tokens[strandColumnIndex];
					String chromStrand = chromosome+strand;
					//parse position
					int position = Integer.parseInt(tokens[positionColumnIndex]);
					//if - strand then add seq length to set position as the 5' end
					if (negativeStrand.matcher(strand).matches()){
						position += tokens[sequenceColumnIndex].length();
					}
					//get PrintWriter
					if (currentChrom.equals(chromStrand) == false){
						currentChrom = chromStrand;
						if (chromOut.containsKey(chromStrand)) dos = chromOut.get(chromStrand);
						else {
							File f = new File(tempDirectory, chromStrand);
							tempChrData.add(f);
							dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f)));
							chromOut.put(chromStrand, dos);
						}
					}
					//save data
					//position
					dos.writeInt(position);
					//length of line
					dos.writeInt(line.length());
					//line
					dos.writeBytes(line);
					//status
					if (++counter == 25000){
						System.out.print(".");
						counter = 0;
					}
				} catch (Exception e){
					System.out.println("\nProblem parsing line -> "+line +" Skipping!\n");
					e.printStackTrace();
				}
			}
			in.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FilterDuplicateAlignments(args);
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
					case 'f': dataFiles = IO.extractFiles(new File(args[++i])); break;
					case 'c': chromosomeColumnIndex = Integer.parseInt(args[++i]); break;
					case 'p': positionColumnIndex = Integer.parseInt(args[++i]); break;
					case 's': sequenceColumnIndex = Integer.parseInt(args[++i]); break;
					case 'b': qualityColumnIndex = Integer.parseInt(args[++i]); break;
					case 't': strandColumnIndex = Integer.parseInt(args[++i]); break;
					case 'm': maxReadCount = Integer.parseInt(args[++i]); saveUniques = false; break;
					case 'r': saveDirectory = new File(args[++i]); break;
					case 'j': skipSpliceJunctions = false; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check params
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printErrAndExit("\nError: cannot find your alignment file(s)!\n");
		if (chromosomeColumnIndex == -1 || positionColumnIndex == -1 || sequenceColumnIndex== -1 || strandColumnIndex == -1) Misc.printErrAndExit("\nPlease enter indexes (0 based) for the chromosome, position, sequence, and strand columns.\n");
		if (saveDirectory == null ) Misc.printErrAndExit("\nPlease provide a directory to use in saving the results.\n");
		saveDirectory.mkdir();
		tempDirectory = new File (saveDirectory, "TempDir_"+Passwords.createRandowWord(7));
		tempDirectory.mkdir();

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          FilterDuplicateAlignments: Mar 2010                     **\n" +
				"**************************************************************************************\n" +
				"Filters alignments for potential amplification bias by randomly selecting X alignments\n" +
				"from those with the same chromosome, position, and strand. Can also filter for the\n" +
				"best unique alignment based on read score. Column indexes start with 0.\n"+

				"\nOptions:\n"+
				"-f Full path file/ directory text containing tab delimited alignments.\n"+
				"-r Full path directory for saving the results.\n"+
				"-c Alignment chromosome column index.\n"+
				"-p Alignment position column index, assumes this is always referenced to the + strand\n"+
				"-s Alignment sequence column index.\n"+
				"-t Strand column index.\n"+
				"-m Save a max number of identical alignments, choose number, defaults to random\n" +
				"        unique sequences.\n" +
				"-b Save only the best alignment per start postion, defined by total score.  Indicate\n" +
				"        which column contains the quality ascii text.\n"+
				"-j Include splice junction chromosomes in filtering (e.g. chr7_101267544_101272320).\n" +
				"        Defaults to removing them. (Only keep for RNA-Seq datasets.)\n"+
				
				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/FilterDuplicateAlignments -f\n" +
				"     /Novoalign/Run7/ -s /Novoalign/Run7/DupFiltered/ -c 7 -p 8 -s 2 -b 3 -t 9\n\n" +
				
				"	  Use -c 10 -p 12 -s 8 -t 13 -b 9  for ELAND sorted or export alignments.\n\n" +

		"**************************************************************************************\n");

	}	

}
