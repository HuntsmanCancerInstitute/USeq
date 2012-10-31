package trans.tpmap;
import java.util.*;
import util.bio.seq.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;
import util.gen.*;

/**
 * Application to map 1lq files to fasta files.
 * Uses Mummer to map oligos to the fasta files.
 * Reverses the oligo sequences in the 1lq file prior to mapping with Mummer. 
 */
public class MummerMapper {
	
	//fields
	private File mummer = new File ("/nfs/transcriptome/software/noarch/T2/64Bit_MUMmer3.18/mummer");
	private File orig1lq;
	private File filtered1lq;
	private String resultsFileName;
	private String finalResultsFileName;
	private String mummerResultsLine = "";
	private int numLinesInFiltered1lq = 0;
	private File genomeDirectory;
	private File[] genomeFastas;
	private File mummerResults;
	private File fasta1lq;
	private Pattern mummerPattern = Pattern.compile("> (\\d+).*"); //for matching mummer lines "> 22345 misc"
	private Pattern colonComma = Pattern.compile(":|,"); //for stripping problem chars
	private Pattern lengthMummerMatch = Pattern.compile("\\d+$"); //for extracting last number 'chr5 39 25 25'
	private Pattern extract1lqOligoLength = Pattern.compile("^\\d+"); //for extracting first number
	private Pattern destype = null;
	private String destypeString = null;
	private boolean readMummer = true;
	private boolean preG7TPMap = true;
	private int numberRows = 0;
	private int minimumOligoSize = 25;
	private boolean printMMCoordinates = true;
	private boolean saveToScratch = false;
	private boolean replaceOriWithHits = false;
	private boolean matchForwardComplement = false;	//if both false then matches forward and reverse
	private boolean matchReverseComplement = false;
	private boolean reverseComplementOligos = false;	//this flips the orientation of the 1lq file, S <-> AS
	private boolean reverseOligos = true; //leave true for 1lq files since their seq is listed 3'- 5'
	private boolean firstMummerResultsRead = true;
	private int maxNumberOfMatches = 100;
	private int numberExactMatches;
	private int numberNoMatches;
	private int number11OrMoreMatches;
	private int number2Matches;
	private int number3Matches;
	private int number4Matches;
	private int number5Matches;
	private int number6Matches;
	private int number7Matches;
	private int number8Matches;
	private int number9Matches;
	private int number10Matches;
	/**String given to naming the control sequences.*/
	public static final String controlChromosomeName = "chrCtrls";
	private boolean controlChromosomePresent = false;
	private boolean removeCrossMatchControlProbes = false;
	
	
	//main constructor
	public MummerMapper(String[] args){		
		//process args and load files
		processArgs(args);
		System.out.println("\n"+(maxNumberOfMatches-1)+"  Maximum number of exact matches ");
		System.out.println(minimumOligoSize+ "  Minimum size of oligos to map ");
		if (preG7TPMap == false){
			System.out.println("Making a rotated, 7G scanner specific tpmap...\n");
		}
		else System.out.println("Making an unrotated, pre 7G scanner tpmap");
		//what kind of matching
		System.out.print("Matching ");
		if (matchForwardComplement) System.out.print("forward ");
		else {
			if (matchReverseComplement) System.out.print("reverse ");
			else System.out.print("both forward and reverse ");
		}
		System.out.println("complement sequences.");
		
		//reversing oligo sequences?
		if (reverseComplementOligos) System.out.println("Reverse complementing 1lq oligo sequences.\n");
		if (reverseOligos == false) System.out.println("Not reversing oligo sequences. This would be inappropriate for real 1lq files!\n");
		
		//filter and make fasta file from 1lq
		System.out.println("\nMaking FASTA file from 1lq file...\n");
		makeFastaAndFilter();
		
		//run mummer to find perfect matches on each genomeFastas file
		boolean go = true;
		mummerResults = new File (resultsFileName+ new Random().nextInt(10000)+"TempMummerMapper.mummerResults");
		for (int i=0; i<genomeFastas.length; i++){
			//skip directories
			if (genomeFastas[i].isDirectory()) continue;
			System.out.println("Running MUMMER on "+genomeFastas[i].getName()+"...");
			//run mummer to generate a mummerResults file
			go = runMummerExactMatcher(genomeFastas[i]);
			//append results onto probe exporter file
			if (go) go = appendMummerHits();
			if (go == false) Misc.printExit("\nFatal error while running mummer!\n");
			//is it a control chromosome?
			if (genomeFastas[i].getName().startsWith(controlChromosomeName)) controlChromosomePresent = true;
		}
		//make tpmap
		System.out.println("\nMaking tpmap...\n");
		buildTPMap();
		
		//remove temp files
		fasta1lq.delete();
		mummerResults.delete();
		filtered1lq.delete();
		
		//print summary statistics
		System.out.println();
		System.out.println("Number of (destype(s) ="+destypeString+") oligos with: ");
		System.out.println("\tNo exact matches\t"+numberNoMatches);
		System.out.println("\t1 exact match\t"   + numberExactMatches);
		System.out.println("\t2 exact matches\t" + number2Matches);
		System.out.println("\t3 exact matches\t" + number3Matches);
		System.out.println("\t4 exact matches\t" + number4Matches);
		System.out.println("\t5 exact matches\t" + number5Matches);
		System.out.println("\t6 exact matches\t" + number6Matches);
		System.out.println("\t7 exact matches\t" + number7Matches);
		System.out.println("\t8 exact matches\t" + number8Matches);
		System.out.println("\t9 exact matches\t" + number9Matches);
		System.out.println("\t10 exact matches\t" + number10Matches);
		System.out.println("\t11 or > exact matches\t" + number11OrMoreMatches);
		System.out.println("\tTotal Matching Oligos\t" + (numberExactMatches + number2Matches + number3Matches + number4Matches + number5Matches +
				number6Matches + number7Matches + number8Matches + number9Matches + number10Matches + number11OrMoreMatches));
	}
	
	/**Runs through mummer results and appends any hits to filtered1lq */
	public boolean appendMummerHits(){
		try{
			BufferedReader pe = new BufferedReader (new FileReader (filtered1lq));
			BufferedReader mummer = new BufferedReader (new FileReader (mummerResults));
			File modified1lq = new File (resultsFileName+ new Random().nextInt(10000)+"TempMummerMapper.mod1lq");
			PrintWriter tpmap = new PrintWriter (new FileWriter (modified1lq));
			
			//run through lines and fetch matching mummer results, set first read to true to load initial mummer result line.
			String probeLine;
			firstMummerResultsRead = true;
			//run thru mummer results
			while (readMummer) {
				//read in the filtered1lq line and parse length from front
				probeLine = pe.readLine();
				Matcher mat = extract1lqOligoLength.matcher(probeLine);
				if (mat.find() == false) Misc.printExit("\nError: problem extracting length of oligo from modified 1lq file.\n");
				int lengthOligo = Integer.parseInt(mat.group());
				//fetch the possible matches from the mummer results
				ArrayList lines = fetchSameNumberLines(mummer, lengthOligo);				
				//append and save
				if (lines.size() !=0) tpmap.println(probeLine+":"+lines);
				else tpmap.println(probeLine);
			}
			
			//check to see if anymore probeExporter lines
			probeLine = pe.readLine();
			if (probeLine !=null){
				System.out.println("\nError: Additional 1lq lines found!  Mummer results don't match. Mapping failed!!!\n");
				return false;
			}
			//close handles
			pe.close();
			mummer.close();
			tpmap.close();
			//rename tempMod to original probe exporter file, this kills the original!
			modified1lq.renameTo(filtered1lq);
			//reset readMummer
			readMummer = true;
			
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	/**Builds a tpmap file from a mummer results file*/
	public void buildTPMap(){
		
		try{
			int numRowsMinusOne = numberRows - 1;
			//make readers and writers 
			BufferedReader pe = new BufferedReader (new FileReader (filtered1lq));
			String line;
			String[] tokens;
			PrintWriter noMatchOut = new PrintWriter(new FileWriter(Misc.removeExtension(finalResultsFileName)+".noMatch"));
		
			
			//run thru modified 1lq file, place tpmaplines in hash set to collapse
			ArrayList linesAL = new ArrayList(numLinesInFiltered1lq);		
			while ((line = pe.readLine()) != null) {
				//split on : or a comma
				tokens = line.split(":|,");
				//number of matches?
				int num = tokens.length;
				
				//score type of matches
				if (num == 2) numberExactMatches++;
				else if (num == 1) numberNoMatches++;
				else if (num == 3) number2Matches++;
				else if (num == 4) number3Matches++;
				else if (num == 5) number4Matches++;
				else if (num == 6) number5Matches++;
				else if (num == 7) number6Matches++;
				else if (num == 8) number7Matches++;
				else if (num == 9) number8Matches++;
				else if (num == 10) number9Matches++;
				else if (num == 11) number10Matches++;
				else number11OrMoreMatches++;
				
				//matches
				if (num >1 && num <= maxNumberOfMatches) {
					//length x y seq misc.... : [(mummerResults) chrom, start, startPosInOligo_1orLength, lengthOfMatch]:[...]:[...]
					//25 0 0 CAGCAGTTCTACGATGGCAAGTCCT -3 control default_at -1 25 0 ! ! ! ! 0 0 0:
					//[chr1 1 1 25, chr5 39 25 25, chrX 74 25 25]:
					//[chrY 5 1 25, chrZ 80 25 25, chr22 44 25 25]
					//make bpmap line 626		205	chr21	9929431	AACAAAAAAAGAACCTCACCATGAA		[chr1         1         1        25]
					//				  0       1   2         3              4                     0            1        2         3
					//extract coordinates and sequence
					String[] peTokens = tokens[0].split("\\s+");					
					String pmX = peTokens[1];
					String pmY = peTokens[2];
					String mmX = "";
					String mmY = "";
					
					//rotate 1lq coordinates for 7G cel files
					if (preG7TPMap == false){					
						int transPmX = Integer.parseInt(pmY);
						int transPmY = numRowsMinusOne - Integer.parseInt(pmX);
						//mm coordinates? Y is the same X differes if the oligo in question is a 111 (PM) or 113 (MM)
						int transMmY = transPmY;
						int transMmX;
						//mm oligo?, mm off to left
						if (peTokens[4].indexOf("113") != -1) transMmX = transPmX -1;
						//mm is off to the right
						else transMmX = transPmX +1;
						
						//reset
						pmX = transPmX+"";
						pmY = transPmY+"";
						mmX = transMmX+"";
						mmY = transMmY+"";
					}
					//normal pre 7G
					else {
						mmX = pmX;
						int parsedPmY = Integer.parseInt(pmY);
						//mm oligo?, mm on top
						if (peTokens[4].indexOf("113") != -1) mmY = (parsedPmY -1)+"";
						//pm oligo there fore mm on bottom
						else mmY = (parsedPmY +1)+"";
					}
					
					//reverse sequence
					String seq; 
					if (reverseComplementOligos) seq = Seq.complementDNA(peTokens[3]);
					else if (reverseOligos) seq = new StringBuffer(peTokens[3]).reverse().toString();
					else seq = peTokens[3];
					
					//drop control probes that match non control genomic sequences?
					boolean printControlProbes = true;
					if (removeCrossMatchControlProbes && controlChromosomePresent){
						for (int i=1; i<num; i++){
							if (tokens[i].indexOf(controlChromosomeName) == -1){
								printControlProbes = false;								
								break;
							}
						}
					}
					
					
					//for each match, create a new tpmap line and place in hash to clear duplicates, is this necessary?
					for (int i=1; i<num; i++){						
						String match = tokens[i].replaceAll("\\[","");
						String[] mummerTokens = match.trim().split("\\s+");
						//make chr text
						String chr = mummerTokens[0];
						//print control probe only if true
						if (printControlProbes == false && chr.equals(controlChromosomeName)) continue;
						//start, subtract one since Affy is in zero coordinates
						int start = Integer.parseInt(mummerTokens[1]) -1;
						//forward orientation? true or false
						String forward = "t";	
						if (replaceOriWithHits) forward = (num -1) + ""; 
						else if (mummerTokens[2].equals("1")) forward = "f";
						//tpmap		sequence t chrom start pmX pmY mmX mmY mystery#1
						String tpLine;
						if (printMMCoordinates) tpLine = seq +"\t"+ forward +"\t"+ chr +"\t"+ start +"\t"+ pmX +"\t"+ pmY +"\t"+ mmX +"\t"+ mmY+"\t"+1;
						else tpLine = seq +"\t"+ forward +"\t"+ chr +"\t"+ start + "\t"+ pmX +"\t"+ pmY +"\t"+1;
						//add to ArrayList
						linesAL.add(new TPMapLine(tpLine));					
					}				
				}
				//print no exact matches
				if (num == 1){
					noMatchOut.println(line);
				}
			}
			//close handles
			pe.close();
			noMatchOut.close();
			
			//convert ArrayList
			int numLines = linesAL.size();
			TPMapLine[] lines = new TPMapLine[numLines];
			linesAL.toArray(lines);
			
			//sort lines
			Arrays.sort(lines);
			
			//print lines to make actual tpmap file
			PrintWriter out = new PrintWriter(new FileWriter(Misc.removeExtension(finalResultsFileName)+".tpmap"));
			
			//add header
			out.println("#seq_group_name\n#version");
			for (int i=0; i<numLines; i++){
				out.println(lines[i].getLine());
			}
			out.close();
			
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	/**Loads the mummer results for one oligo into an ArrayList.  Checks each mummer hit to ensure it is equal to the oligo length.
	 * No partial matches.*/
	public ArrayList fetchSameNumberLines (BufferedReader lines, int oligoLength){
		ArrayList al = new ArrayList();
		try {
			boolean go = true;
			Matcher mat = null;
			String lineNumber = null;
			while (go){
				//read in first line?
				if (firstMummerResultsRead) {
					//set to not read in again
					firstMummerResultsRead = false;
					//load it
					mummerResultsLine = lines.readLine();
				}
				//if null return, stop of file reached
				if (mummerResultsLine == null){
					readMummer = false;
					return al;
				}
				//nope, result line found
				mat = mummerPattern.matcher(mummerResultsLine);
				//does result begin with a '>' indicating a new oligo? 
				if (mat.matches()){
					String newLineNumber = mat.group(1);
					//has the lineNumber been set yet?
					if (lineNumber == null) lineNumber = newLineNumber;
					//is it the same number? if not return matches.
					else if (lineNumber.equals(newLineNumber)==false) return al;
				}
				//must be a mummer match or possibly a blank line
				else {
					mummerResultsLine = mummerResultsLine.trim();
					if (mummerResultsLine.length()!=0){
						//check length of match
						Matcher matcher = lengthMummerMatch.matcher(mummerResultsLine);
						if (matcher.find()){
							int length = Integer.parseInt(matcher.group());
							if (length == oligoLength){
								al.add(mummerResultsLine);
							}
						}
						else Misc.printExit("Error: problem extracting length of mummer match from stop of line "+mummerResultsLine);
					}
				}
				//read next line
				mummerResultsLine = lines.readLine();
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		return al;
	}
	
	/**Runs MUMMER to find perfect matches*/
	public boolean runMummerExactMatcher(File genomeFastaFile){
		try {
			//make command line, -b looks forward and reverse, -r reverse complement, nothing is forward complement
			String[] commandArray;
			if (matchForwardComplement) commandArray = new String[]{mummer.getCanonicalPath(),"-n", "-l", minimumOligoSize+"", "-F", "-maxmatch", genomeFastaFile.getCanonicalPath(),fasta1lq.getCanonicalPath()};
			else {
				String match = "-b";
				if (matchReverseComplement) match = "-r";
				commandArray = new String[]{mummer.getCanonicalPath(),"-n", "-l", minimumOligoSize+"", match,"-c", "-F", "-maxmatch", genomeFastaFile.getCanonicalPath(),fasta1lq.getCanonicalPath()};
			}
			
			System.out.println ("Mummer command ->"+Misc.stringArrayToString(commandArray," "));
			Runtime rt = Runtime.getRuntime();
			//rt.traceInstructions(true); //for debugging
			//rt.traceMethodCalls(true); //for debugging
			Process p = rt.exec(commandArray);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			//BufferedReader data = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			PrintWriter out = new PrintWriter(new FileWriter(mummerResults));
			String line;
			while ((line = data.readLine()) != null){
				//System.out.println(line);
				out.println(line);
			}
			data.close();
			out.close();
			return true;
		}catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}
	
	/**Makes a fasta file from a 1lq file, reverses the sequence of each oligo, drops control probes*/
	public void makeFastaAndFilter(){
		//fields
		String line;
		String[] tokens;
		int counter = 1;
		try{
			//make files
			filtered1lq = new File (resultsFileName+ new Random().nextInt(10000)+ "TempMummerMapper.filtered1lq"); 
			fasta1lq = new File (resultsFileName+ new Random().nextInt(10000)+"TempMummerMapper.fasta");
			PrintWriter outFasta = new PrintWriter(new FileWriter(fasta1lq));
			PrintWriter outFiltered1lq = new PrintWriter(new FileWriter(filtered1lq));
			BufferedReader in = new BufferedReader(new FileReader(orig1lq));
			String reversed;
			//skip header
			while (true){
				line = in.readLine();
				if (line.startsWith("X")) break;
			}
			//write sequences
			//filter for particular destypes?
			if (destype != null){
				Matcher mat;
				while ((line = in.readLine()) !=null) {
					//X       Y       Seq     DESTYPE
					//0       1        2       3
					line = line.trim();
					if (line.length() == 0) continue;
					tokens = line.split("\\s+");
					mat = destype.matcher(tokens[3]);
					if (tokens[2].length() >= minimumOligoSize && ( mat.matches() )){
						//remove : and , from line to prevent improper splitting when making tpmap
						Matcher strip = colonComma.matcher(line);
						line = strip.replaceAll("_");
						//print 1lq line prepended by length of sequence
						outFiltered1lq.println(tokens[2].length()+"\t"+line);
						numLinesInFiltered1lq++;
						//reversing oligo sequence?
						if (reverseComplementOligos) reversed = Seq.complementDNA(tokens[2]);
						else if (reverseOligos) reversed = new StringBuffer(tokens[2]).reverse().toString();
						else reversed = tokens[2];
						outFasta.println(">"+counter+"\n"+reversed);
					}
					counter++;
				}
			}
			//no filtering add everything
			else {
				while ((line = in.readLine()) !=null) {
					//X       Y       Seq     etc
					//0       1        2     
					line = line.trim();
					if (line.length() == 0) continue;
					tokens = line.split("\\s+");
					if (tokens[2].length() >= minimumOligoSize ){
						//remove : and , from line to prevent improper splitting when making tpmap
						Matcher strip = colonComma.matcher(line);
						line = strip.replaceAll("_");
						//print 1lq line prepended by length of sequence
						outFiltered1lq.println(tokens[2].length()+"\t"+line);
						numLinesInFiltered1lq++;
						//reversing oligo sequence!!!
						if (reverseComplementOligos) reversed = Seq.complementDNA(tokens[2]);
						else if (reverseOligos) reversed = new StringBuffer(tokens[2]).reverse().toString();
						else reversed = tokens[2];
						outFasta.println(">"+counter+"\n"+reversed);
					}
					//else System.out.println("\tSkipping -> "+line);
					counter++;
				}
			}
			in.close();
			outFiltered1lq.close();
			outFasta.close();
		}catch(Exception e){
			e.printStackTrace();
			Misc.printExit("\nError: Problem parsing 1lq file!\n");
		}
	}
	
	public static void main (String[] args){
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		long start = System.currentTimeMillis();
		new MummerMapper(args);
		int elapse = (int)(System.currentTimeMillis()-start)/1000;
		System.out.println("\nFinished "+elapse+" seconds");
	}
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File resultsDirectory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': orig1lq = new File(args[i+1]); i++; break;
					case 'g': genomeDirectory = new File(args[i+1]); i++; break;
					case 'm': mummer = new File(args[i+1]); i++; break;
					case 'r': resultsDirectory = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					case 'n': preG7TPMap = false; numberRows = Integer.parseInt(args[i+1]); i++; break;
					case 'd': printMMCoordinates = false; break;
					case 'o': replaceOriWithHits = true; break;
					case 'f': destypeString = args[i+1]; i++; break;
					case 's': saveToScratch = true; break;
					case 'e': removeCrossMatchControlProbes = true; break;
					case 'a': matchForwardComplement = true; matchReverseComplement = false; break;
					case 'b': matchReverseComplement = true; matchForwardComplement = false;break;
					case 'c': reverseComplementOligos = true; break;
					case 'i': reverseOligos = false; break;
					case 'x': maxNumberOfMatches = Integer.parseInt(args[i+1]); i++; break;
					case 'y': minimumOligoSize = Integer.parseInt(args[i+1]); i++; break;
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}		
		
		//check to see if they entered required params
		String error = null;
		if (orig1lq == null || orig1lq.canRead()== false) error = "-p 1lq file";
		else if (genomeDirectory == null || genomeDirectory.canRead()==false) error = "-g genome sequence file(s)";
		else if (mummer == null || mummer.canRead() == false) error = "-m mummer application";
		else if (resultsDirectory == null || resultsDirectory.isDirectory()==false) error = "-r results directory";
		if (error != null) Misc.printExit("\nCannot find or read your "+error+"!\n");
		
		//fetch genome files
		genomeFastas = IO.extractFiles(genomeDirectory);
		resultsFileName = IO.getFullPathName(resultsDirectory)+File.separator+orig1lq.getName();
		//using scratch?
		finalResultsFileName = resultsFileName;
		if (saveToScratch){
			resultsFileName = "/scratch/"+orig1lq.getName();
		}
		//filterForDestypes?
		if (destypeString != null){
			String[] types = destypeString.split(",");
			destypeString = Misc.stringArrayToString(types,"|");
			destype = Pattern.compile(destypeString);
		}
		//increment numberOfMatches by one
		maxNumberOfMatches++;
		
		
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Mummer Mapper: April 2007                            **\n" +
				"**************************************************************************************\n" +
				"MM uses MUMMER to generate a 7G or pre 7G tpmap from a standard, unrotated, 1lq file\n" +
				"using genomic fasta file(s). Note, DO NOT use a rotated 7G specific 1lq file! Oligos\n" +
				"in the 1lq file are reversed prior to mapping. MUMMER matches are case-insensitive but,\n" +
				"restricted to GATC, thus use Ns or Xs to mask your genomic fasta files. For large\n" +
				"fasta's compile and run MUMMER on a 64 bit machine. Coordinate are zero based and\n" +
				"relative to the forward (+) strand. In addition to the genomic fasta files, it may be\n" +
				"usefull to include one named "+controlChromosomeName+ " and populated with control sequences\n" +
				"found on the array (e.g. bacterial, Arabidopsis, known regions that don't change).\n" +
				"These can be used by the TiMAT2 CelProcessor app for defined region median scaling.\n" +
				"See MUMMER from http://sourceforge.net/projects/mummer .\n\n" +
				
				"Parameters:\n"+
				"    -p Full path file text for the 1lq file.\n" +
				"    -r Full path directory file text to save the results.\n"+
				"    -g Full path file or directory text containing genomic (multi) fasta file(s).\n" +
				"    -m Full path file text for the mummer application. Defaults to\n"+
				"          /nfs/transcriptome/software/noarch/T2/64Bit_MUMmer3.18/mummer\n"+
				"    -x Maximum number of exact matches allowed per oligo, defaults to 100, set to 1\n" +
				"         for unique oligos.\n"+
				"    -y Minimum size of oligos to map, others will be excluded. Defaults to 25.\n"+
				"    -a Match forward complement, defaults to both.\n"+
				"    -b Match reverse complement, defaults to both.\n"+
				"    -c Reverse complement the 1lq oligo sequences. This switches the orientation of\n" +
				"         the 1lq file S <-> AS. (Use of -b -c is the same as -a.)\n"+
				"    -i Don't reverse the oligo sequence, by default oligo sequences are reversed\n"+
				"         since seqs in 1lq format are listed 3' to 5'.\n"+
				"    -o Replace the orientation column with the number of exact matches for a\n" +
				"         particular oligo in the genome. Doing so enables graphing of matches when\n" +
				"         generating interval plots in TiMAT2 (recommended).\n"+
				"    -s Save temp files to '/scratch/'.\n"+
				"    -e Remove control probes that also match non-control fasta sequences.\n"+
				"    -f Enter a comma delimited list, no spaces, of particular DESTYPES to map, others\n" +
				"         will be skipped. Using '-111,111' is recommended! (PM= -111,111; MM= -113,113;\n" +
				"         'manufacture/ gridding controls= 0,1,-1,-3,3,-4,4,-6,6').\n"+
				"    -d Don't print mismatch coordinates (ie they don't exist).\n"+
				"    -n Make a new 7G scanner tpmap, enter the number of rows/ columns in a cel file,\n" +
				"         default is to make an old pre-7G tpmap.\n\n"+
		
				"MM probes are assumed to be directly below their PM counter-part on the array. If\n" +
				"    mapping MM oligos (not recommended) their 'MM but really PM' counter-part is\n" +
				"    assumed to be directly above.\n"+
				"To make an antisense stranded tpmap from an antisense 1lq file use option -b. To\n"+
				"    make a sense stranded tpmap from an antisense 1lq file use options -c -a.\n"+
				"To make a mock 1lq file, create a tab delimited text file with the following:\n"+
				"    X  Y  Seq  Destype  ... Be sure to have this header immediately preceed your map\n"+
				"    data.  The Destype column is optional but can be used to filter your mock 1lq\n"+
				"    for particular lines.  If needed, use the -i flag to prevent reversing of your\n" +
				"    sequences. Rember that X and Y are zero based. The PrintSelectColumns app can be\n" +
				"    used to manipulate your map file to make a mock 1lq.\n\n"+
				
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/MummerMapper -p /affy/human.1lq -f\n" +
				"    -111,111 -g /affy/human/NCBIv34/ -m\n" +
				"    /nfs/transcriptome/software/noarch/T2/64Bit_MUMmer3.18/mummer -r /affy/results/\n" +
				"    -s -x 10 -n 2560 -y 22\n" +
				"\n" +
		"**************************************************************************************\n");
	}
}
