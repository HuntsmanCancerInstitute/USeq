package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;

/**Parses a Novoalign paired alignment txt file into 12 column bed file format.
 * @author david.nix@hci.utah.edu 
 **/
public class NovoalignPairedParser{
	//fields
	private File[] dataFiles;
	private File workingFile;
	private Pattern tab = Pattern.compile("\\t+");
	private int radiusSpliceJunction = 34;
	private int averageFragmentLength = 110;
	public static final Pattern COMMENT = Pattern.compile("^\\s*#.*");
	public static final Pattern QC = Pattern.compile(".+QC$");
	private Pattern adapter = Pattern.compile(".*chrAdapt.*");
	//pulls chromosome start stop for splice junctions
	public static final Pattern SPLICE = Pattern.compile(">(.+)_(\\d+)_(\\d+)$");
	private boolean saveHalfMatches = true;
	private PrintWriter out = null;
	private int maxSizePair = 100000;

	//indexes for columns
	private int readNameIndex = 0;
	private int sequenceIndex = 2;
	private int qualityIndex = 3;
	private int posteriorProbabilityIndex = 6;
	private int chromosomeIndex = 7;
	private int positionIndex = 8;
	private int strandIndex = 9;
	private int numberOddPairs =0;
	private PrintWriter oddPairsOut;

	//constructors
	public NovoalignPairedParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		System.out.println("\nParsing and filtering...");

		//for each file, parse and save to disk	
		for (int i=0; i< dataFiles.length; i++){
			//set working objects and parse tag file text
			workingFile = dataFiles[i];
			System.out.print("\t"+workingFile);
			//make PrintWriter
			String name = workingFile.getName();
			name = name.replace(".gz", "");
			name = name.replace(".zip", "");
			name = Misc.removeExtension(name);
			File bed = new File (workingFile.getParentFile(), name+"_NPP.bed");
			File oddPairs = new File (workingFile.getParentFile(), name+"_NPPOddPairs.txt");
			try {
				out = new PrintWriter(new FileWriter (bed));
				oddPairsOut = new PrintWriter(new FileWriter (oddPairs));
				//split file to chromosome strand specific temp files
				boolean  parsed = parseWorkingFile(); 
				if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
				
				out.close();
				oddPairsOut.close();
				//odd pairs
				System.out.println("\tNumber odd pairs "+numberOddPairs);
				numberOddPairs = 0;
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public boolean parseWorkingFile(){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			String line;
			String[] tokens = null;
			int counter =0;
			String workingReadName = "";
			HashMap<String, String[]> lines = new HashMap<String, String[]>();
			while ((line = in.readLine()) !=null){
				//print status blip
				if (++counter == 25000){
					System.out.print(".");
					counter = 0;
				}
				//comment line
				if (COMMENT.matcher(line).matches()) continue;
				//qc failed?
				if (QC.matcher(line).matches()) continue;
				//adapter chrom?
				if (adapter.matcher(line).matches()) continue;
				//check number of rows
				tokens = tab.split(line);
				if (tokens.length < 5) {
					//System.err.println("\nToo few columns, skipping -> "+line);
					continue;
				}
				//check size of read against size of quality score
				if (tokens[sequenceIndex].length() != tokens[qualityIndex].length()) {
					//System.err.println("\nSeq length != Qual length, skipping -> "+line);
					continue;
				}
				//check read name
				String testReadName = tokens[0].substring(0, tokens[0].length()-2);
				if (testReadName.equals(workingReadName)) lines.put(line, tokens);
				else {
					processGroup(lines);
					lines.clear();
					lines.put(line, tokens);
					workingReadName = testReadName;
				}
			}
			//process last
			processGroup(lines);
			System.out.println();
		} catch (Exception e){
			e.printStackTrace();
		}
		return true;
	}

	public String makeBedLineFromHalfMatch(String[][] pair){
		//get match
		String[] match = pair[0];
		if (pair[0].length < 13) match = pair[1];
		//set strand
		String strand = "+";
		if (match[strandIndex].equals("F") == false) strand = "-";
		//splice junction?
		String chromosome;
		int start = 0;
		int stop = 0;
		String name;
		int blockStart;
		Matcher mat = SPLICE.matcher(match[chromosomeIndex]);
		if (mat.matches()){
			chromosome = mat.group(1);
			//splice position - radius + actual start of read
			if (strand.equals("+")) {
				start = Integer.parseInt(mat.group(2)) - radiusSpliceJunction + Integer.parseInt(match[positionIndex]) -1;
				stop = start + averageFragmentLength;
				blockStart = 0;
			}
			else {
				stop = Integer.parseInt(mat.group(3)) - radiusSpliceJunction + Integer.parseInt(match[positionIndex]) + match[sequenceIndex].length() -1;
				start = stop - averageFragmentLength;
				blockStart = (stop-start-match[sequenceIndex].length());
			}
			name = "HalfMatchSpliceJunction_"+match[chromosomeIndex];
		}
		//parse position
		else {  
			chromosome = match[chromosomeIndex].substring(1);
			name = "HalfMatch";
			start = Integer.parseInt(match[positionIndex]) -1;
			blockStart = 0;
			if (match[strandIndex].equals("F")) stop = start + averageFragmentLength;
			else {
				stop = start+ match[sequenceIndex].length();
				start = stop - averageFragmentLength;
				blockStart = (stop-start-match[sequenceIndex].length());
			}
		}
		
		//watch out for negative starts! on chrM
		if (start <= 0) return null;
		//check length
		if ((stop-start)> maxSizePair) return null;
		int length = match[sequenceIndex].length();
		return chromosome+"\t"+start+"\t"+stop+"\t"+ name+ "\t"+ match[posteriorProbabilityIndex]+"\t"+strand+"\t"+start+"\t"+stop+"\t0\t1\t"+length+"\t"+blockStart;
	}

	public String makeBedLineFromSameChromosomeMatch(String[][] pair){
		//find F-R
		String[] left;
		String[] right;
		if (pair[0][strandIndex].equals("F")) {
			left = pair[0];
			right = pair[1];
		}
		else {
			left = pair[1];
			right = pair[0];
		}
		//parse start stop chrom
		int start = Integer.parseInt(left[positionIndex])-1;
		int stop = Integer.parseInt(right[positionIndex])+ right[sequenceIndex].length()-1;
		String chrom = left[chromosomeIndex].substring(1);
		
		
		
		//check length
		int size = stop - start;
		if (size > maxSizePair || size < left[sequenceIndex].length()) {
			oddPairsOut.println("Skipping Odd Match Pair Length "+size);
			oddPairsOut.println("L:\t"+Misc.stringArrayToString(left, "\t"));
			oddPairsOut.println("R:\t"+Misc.stringArrayToString(right, "\t"));
			numberOddPairs++;
			return null;
		}
		//save
		return chrom+"\t"+start+"\t"+stop+"\tPaired\t"+left[posteriorProbabilityIndex]+"\t+\t"+start+"\t"+stop+"\t0\t2\t"+left[sequenceIndex].length()+","+right[sequenceIndex].length()+"\t0,"+(stop-start-right[sequenceIndex].length());

	}

	public String makeBedLineFromPairedSpliceMatch(String[][] pair){
		//find F-R
		String[] left;
		String[] right;
		if (pair[0][strandIndex].equals("F")) {
			left = pair[0];
			right = pair[1];
		}
		else {
			left = pair[1];
			right = pair[0];
		}
		//parse chromosome
		Matcher matL = SPLICE.matcher(left[chromosomeIndex]);
		matL.matches();
		String chrom = matL.group(1);

		//parse start stop chrom
		int start = Integer.parseInt(matL.group(2)) - radiusSpliceJunction + Integer.parseInt(left[positionIndex]) -1;
		Matcher matR = SPLICE.matcher(right[chromosomeIndex]);
		matR.matches();
		int stop = Integer.parseInt(matR.group(3)) - radiusSpliceJunction + Integer.parseInt(right[positionIndex]) + right[sequenceIndex].length() -1;
		
		//check length
		int size = stop - start;
		if (size > maxSizePair || size < left[sequenceIndex].length()) {
			oddPairsOut.println("Skipping Odd Paired Splice Length "+size);
			oddPairsOut.println("L:\t"+Misc.stringArrayToString(left, "\t"));
			oddPairsOut.println("R:\t"+Misc.stringArrayToString(right, "\t"));
			numberOddPairs++;
			return null;
		}
		//save
		return chrom+"\t"+start+"\t"+stop+"\tPairedSplices\t"+left[posteriorProbabilityIndex]+"\t+\t"+start+"\t"+stop+"\t0\t2\t"+left[sequenceIndex].length()+","+right[sequenceIndex].length()+"\t0,"+(stop-start-right[sequenceIndex].length());
	}

	public String makeBedLineFromSingleSpliceMatch(String[][] pair){
		//find F-R
		String[] left;
		String[] right;
		if (pair[0][strandIndex].equals("F")) {
			left = pair[0];
			right = pair[1];
		}
		else {
			left = pair[1];
			right = pair[0];
		}
		//get start, stop, and chrom
		int start;
		int stop;
		String chrom;
		//attempt to match to splice
		Matcher mat = SPLICE.matcher(left[chromosomeIndex]);
		if (mat.matches()){
			chrom = mat.group(1);
			start = Integer.parseInt(mat.group(2)) - radiusSpliceJunction + Integer.parseInt(left[positionIndex]) -1; 
			stop = Integer.parseInt(right[positionIndex])+ right[sequenceIndex].length()-1; 
		}
		else {
			start = Integer.parseInt(left[positionIndex]) -1;
			chrom = left[chromosomeIndex].substring(1);
			Matcher matR = SPLICE.matcher(right[chromosomeIndex]);
			matR.matches();
			stop = Integer.parseInt(matR.group(3)) - radiusSpliceJunction + Integer.parseInt(right[positionIndex]) + right[sequenceIndex].length()-1;
		}
		
		//check length
		int size = stop - start;
		if (size > maxSizePair || size < left[sequenceIndex].length()) {
			oddPairsOut.println("Skipping Odd Single Splice Pair Length "+size);
			oddPairsOut.println("L:\t"+Misc.stringArrayToString(left, "\t"));
			oddPairsOut.println("R:\t"+Misc.stringArrayToString(right, "\t"));
			numberOddPairs++;
			return null;
		}
		
		//save
		return chrom+"\t"+start+"\t"+stop+"\tSingleSplice\t"+left[posteriorProbabilityIndex]+"\t+\t"+start+"\t"+stop+"\t0\t2\t"+left[sequenceIndex].length()+","+right[sequenceIndex].length()+"\t0,"+(stop-start-right[sequenceIndex].length());

	}


	public void processGroup(HashMap<String, String[]> lines){
		int numLines = lines.size();
		//too few lines
		if (numLines < 2) return;
		//System.out.println("\nNewGroup:");
		//convert to String[]
		String[][] alignments = new String[numLines][];
		Iterator<String>it = lines.keySet().iterator();
		int index =0;
		while (it.hasNext()){
			String a = it.next();
			//System.out.println(a);
			alignments[index++] = lines.get(a);
		}
		//pair
		if (numLines == 2){
			processAlignmentPair(alignments);
		}
		//multiple alignments
		else {
			//System.out.println("\tMultiple alignments, skipping");
		}

	}

	public void processAlignmentPair (String[][] pair){
		//watch out for duplicates
		if (pair[0][readNameIndex].equals(pair[1][readNameIndex])) {
			//System.out.println("\tRepetitive");
			return;
		}
		//both no match
		if (pair[0].length < 13 && pair[1].length <13){
			//System.out.println("\tBoth NM");
			return;
		}
		//one or other no match
		if (saveHalfMatches && (pair[0].length < 13 || pair[1].length <13)){
			//System.out.println("\tOne NM");
			String line = makeBedLineFromHalfMatch(pair);
			if (line!=null) out.println(line);
		}
		//both lines hit same chromosome, tokens[chromosomeIndex]
		else if (pair[0][chromosomeIndex].equals(pair[1][chromosomeIndex])){
			//System.out.println("\tSame chromosome");
			String line = makeBedLineFromSameChromosomeMatch(pair);
			if (line!=null) out.println(line);
		}
		//different chroms?
		else {
			//splices?
			Matcher mat0 = SPLICE.matcher(pair[0][chromosomeIndex]);
			Matcher mat1 = SPLICE.matcher(pair[1][chromosomeIndex]);
			boolean mat0Matches = mat0.matches();
			boolean mat1Matches = mat1.matches();
			if (mat0Matches || mat1Matches){
				//both match?
				if (mat0Matches && mat1Matches){
					//System.out.println("\tBoth match splices!");
					String chrom0 = mat0.group(1);
					String chrom1 = mat1.group(1);
					//same chrom?
					if (chrom0.equals(chrom1)){
						//System.out.println("\tSame chrom");
						String line = makeBedLineFromPairedSpliceMatch(pair);
						if (line!=null) out.println(line);
					}
					//else //System.out.println("\tDiff chroms");
				}
				//only one matches
				else {
					//System.out.println("\tOne splice match");
					String chrom0; 
					String chrom1; 
					if (mat0Matches) chrom0 = ">"+mat0.group(1);
					else chrom0 = pair[0][chromosomeIndex];
					if (mat1Matches) chrom1 = ">"+mat1.group(1);
					else chrom1 = pair[1][chromosomeIndex];
					//same chrom?
					if (chrom0.equals(chrom1)){
						//System.out.println("\tSame chrom");
						String line = makeBedLineFromSingleSpliceMatch(pair);
						if (line!=null) out.println(line);
					}
					//else System.out.println("\tDiff chroms");
				}
			}
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new NovoalignPairedParser(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[i+1]); i++; break;
					case 'e': saveHalfMatches = false; break;
					case 'm': maxSizePair = Integer.parseInt(args[++i]); break;
					case 's': radiusSpliceJunction = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull files
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".txt");
		tot[1] = IO.extractFiles(forExtraction,".txt.zip");
		tot[2] = IO.extractFiles(forExtraction,".txt.gz");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.txt(.zip/.gz) file(s)!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Novoalign Paired Parser: January 2009                     **\n" +
				"**************************************************************************************\n" +
				"Parses Novoalign paired alignment files xxx.txt(.zip/.gz) into xxx.bed format.\n" +

				"\nOptions:\n"+
				"-f The full path directory/file text of your Novoalign xxx.txt(.zip or .gz) file(s).\n" +
				"-e Exclude half matches with a high quality unmatched pair, defaults to keeping them.\n"+
				"-m Maximum size for paired reads mapping to the same chromosome, defaults to 100000.\n"+
				"-s Splice junction radius, defaults to 34. See the MakeSpliceJunctionFasta app.\n"+


				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/NovoalignPairedParser -f /Novo/Run7/\n" +
				" \n\n" +

		"**************************************************************************************\n");

	}	

}
