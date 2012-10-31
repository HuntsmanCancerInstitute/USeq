
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.bio.parsers.*;
import edu.utah.seq.data.*;
import util.gen.*;

import java.util.*;

/**Parses an Eland Extended xxx.export.txt or xxx.sorted.txt for sequencing information.
 * Generates track data representing the sum of the base quality scores for every position, GATC
 * as well as a consensus track (1-fraction consensus)
 * Positions are in interbase coordinates (0 start, stop excluded).
 * @author david.nix@hci.utah.edu 
 **/
public class ElandSequenceParser {
	//fields
	private File[] dataFiles;
	private File saveDirectory;
	private File tempDirectory;
	private String versionedGenome;
	private HashMap <String, ArrayList<File>> chromFile = new HashMap <String, ArrayList<File>> ();
	private float minimumAlignmentScore = 13;
	private float minimumConsensusScore = 60;
	private HashMap<String, File> fastas;
	private static final Pattern splitPat = Pattern.compile("\\s");
	private int numFailQC = 0;
	private int numNoAlign = 0;
	private int numPoorAlign = 0;
	private int numPass = 0;



	//constructors
	public ElandSequenceParser(String[] args){
		processArgs(args);


		//for each  file, split by chromosome
		System.out.println("Spliting by chromosome...");
		for (int i=0; i< dataFiles.length; i++) {
			System.out.println("\t"+dataFiles[i].getName());
			splitByChromosome(dataFiles[i]);
		}

		//stats
		int total = numPass+numFailQC+numNoAlign+numPoorAlign;
		System.out.println("\nRead Statistics:");
		System.out.println("\t"+total+"\tTotal\n");
		System.out.println("\t"+formatFraction(numFailQC,total,2)+"\tFail QC");
		System.out.println("\t"+formatFraction(numNoAlign,total,2)+"\tNo alignment");
		System.out.println("\t"+formatFraction(numPoorAlign,total,2)+"\tPoor alignment (score < "+minimumAlignmentScore+")");
		System.out.println("\t"+formatFraction(numPass,total,2)+"\tAligned");

		//for each chromosome
		System.out.println("\nParsing...");
		Iterator<String> it = chromFile.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			System.out.println("\t"+chrom);
			//load chromosome
			File fasta = fastas.get(chrom);
			if (fasta == null) {
				Misc.printExit("\nError: cannot find appropriate fasta file? Where's "+chrom+".fasta ?");
			}
//System.out.println("\t\tMFP "+Misc.fetchUsedMemory());
			MultiFastaParser mfp = new MultiFastaParser(fasta);
			//make int arrays to hold gatc counts
			int[][] gatc = new int[4][mfp.getSeqs()[0].length()+100];
			//parse data loading gatc counts
//System.out.println("\t\tChromData "+Misc.fetchUsedMemory());
			parseChromData(chrom, gatc);
			//write out bars
//System.out.println("\t\tWriteGATC "+Misc.fetchUsedMemory());
			writeGATCBars(chrom, gatc);
			//score and write consensus
//System.out.println("\t\tConsensus "+Misc.fetchUsedMemory());
			float[] fractionConsensus = scoreConsensus(gatc, mfp.getSeqs()[0].toUpperCase(), minimumConsensusScore);
//System.out.println("\t\tWriteConsensus "+Misc.fetchUsedMemory());
			writeConsensusBars(chrom, fractionConsensus);
		}

		//remove temp directory
		IO.deleteDirectory(tempDirectory);

		System.out.println("\nDone!\n");
	}

	public float[] scoreConsensus(int[][] gatc, String seq, float minTotalToScore){
		float[] fract = new float[gatc[0].length];
		for (int i=0; i< fract.length; i++){
			//sum
			float total =0;
			for (int j=0; j< gatc.length; j++) total += gatc[j][i];
			if (total < minTotalToScore) {
				fract[i] = -1.0f;
				continue;
			}
			//score position
			if (seq.charAt(i) == 'G') fract[i] = ((float)gatc[0][i])/total;
			else if (seq.charAt(i) == 'A') fract[i] = ((float)gatc[1][i])/total;
			else if (seq.charAt(i) == 'T') fract[i] = ((float)gatc[2][i])/total;
			else if (seq.charAt(i) == 'C') fract[i] = ((float)gatc[3][i])/total;
			//non gatc in reference
		}
		return fract;
	}

	public void parseChromData (String chromosome, int[][] gatc){
		ArrayList<File> files = chromFile.get(chromosome);
		for (int i=0; i< files.size(); i++){
			parseChromFile (files.get(i), gatc);
		}
	}

	public void parseChromFile(File file, int[][] gatc){
		try{
			BufferedReader in = new BufferedReader (new FileReader (file));
			String line;
			//read in file and print out to tempDirectory
			while ((line = in.readLine()) != null) new ElandExportLine(line, gatc);
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**Returns number space (percent%)*/
	public static String formatFraction(int number, int total, int numberDecimals){
		double fraction = 100.0 *(double)number / (double)total ;
		String fractionString = Num.formatNumber(fraction, numberDecimals);
		return fractionString+"%\t"+number;
	}

	public void splitByChromosome(File dataFile){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(dataFile);
			//make writers for those that pass qc but don't align and those that don't align well
			String name = dataFile.getName().substring(0,dataFile.getName().indexOf(".txt"));
			PrintWriter noAlign = new PrintWriter (new FileWriter( new File (saveDirectory, name+"_NoAlign.txt")));
			PrintWriter poorAlign = new PrintWriter (new FileWriter( new File (saveDirectory, name+"_PoorAlign.txt")));
			HashMap<String,PrintWriter> pws = new HashMap<String,PrintWriter>();
			HashMap<String,File> files = new HashMap<String,File>();
			String line;
			String[] tokens;
			String workingChromosome = "";
			PrintWriter workingPw = null;
			//read in file and print out to tempDirectory
			while ((line = in.readLine()) != null){
				tokens = splitPat.split(line);
				if (tokens.length !=22) {
					IO.deleteDirectory(tempDirectory);
					Misc.printExit("Error: line does not contain enough columns -> \n\t"+line+"\n\tFrom file "+dataFile);
				}
				//trim tokens
				Misc.trim(tokens);
				//worth parsing? 
				//pass QC?
				if (tokens[21].equals("N")) {
					numFailQC++;
					continue;
				}
				//was aligned?
				if (tokens[12].equals("")){
					noAlign.println(line);
					numNoAlign++;
					continue;
				}
				//passed alignment score?
				if (Float.parseFloat(tokens[15]) < minimumAlignmentScore) {
					poorAlign.println(line);
					numPoorAlign++;
					continue;
				}
				//parse the chromosome
				int index = tokens[10].indexOf('.');
				//anything found?
				if (index == -1) {
					IO.deleteDirectory(tempDirectory);
					Misc.printExit("Error: missing '.' (e.g. chrX.fasta) from chromosome text cannot parse! See -> \n\t"+line+"\n\tFrom file "+dataFile);
				}
				//parse it
				numPass++;
				String chromosome = tokens[10].substring(0, index);
				//get print writer
				if (chromosome.equals(workingChromosome) == false){
					workingChromosome = chromosome;
					//does pw exist?
					if (pws.containsKey(chromosome)) workingPw = pws.get(chromosome);
					else {
						//make new file
						File chrFile = new File (tempDirectory, chromosome+"_"+dataFile.getName());
						files.put(chromosome, chrFile);
						//make PW
						workingPw = new PrintWriter (chrFile);
						pws.put(chromosome, workingPw);
					}
				}
				//print line
				workingPw.println(Misc.stringArrayToString(tokens, "\t"));
			}
			in.close();
			noAlign.close();
			poorAlign.close();
			//close pws and add files to master hash
			Iterator<String> it = pws.keySet().iterator();
			while (it.hasNext()) {
				String chr = it.next();
				pws.get(chr).close();
				ArrayList<File> fileAL;
				if (chromFile.containsKey(chr)) fileAL = chromFile.get(chr);
				else {
					fileAL = new ArrayList<File>();
					chromFile.put(chr, fileAL);
				}
				fileAL.add(files.get(chr));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void writeConsensusBars (String chrom, float[] frac){
		File consensusDir = new File (saveDirectory, "Consensus");
		if (consensusDir.exists() == false) consensusDir.mkdir();
		//make Point[]
		ArrayList<Point> pointsAL = new ArrayList<Point>();
		for (int x=0; x< frac.length; x++) {
			if (frac[x] != -1.0f) {
				pointsAL.add(new Point(x,1.0f - frac[x]));
			}
		}
		if (pointsAL.size() == 0) return;
		Point[] points = new Point[pointsAL.size()];
		pointsAL.toArray(points);
		//make notes
		HashMap <String,String> notes = new HashMap <String,String> ();
		notes.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		notes.put(BarParser.SOURCE_TAG, IO.concatinateFileFullPathNames(dataFiles, ","));
		notes.put(BarParser.UNIT_TAG, "Eland extended base quality score sum");
		notes.put(BarParser.DESCRIPTION_TAG, "Fraction consensus track. Generated by running the ElandSeqSumParser on Solexa ELAND file(s), interbase coordinates");
		//make an Info object  public Info (String text, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
		Info info = new Info(chrom, versionedGenome, chrom, ".", 0, notes);
		//make pd
		PointData pd = Point.extractPositionScores(points);			
		pd.setInfo(info);
		//write to file
		File barFile = new File (consensusDir, info.getChromosome()+".bar");
		pd.writePointDataToFile(barFile);
		//cleanup
		pd = null;
		pointsAL = null;
		points = null;
		System.out.print(".");

	}



	public void writeGATCBars (String chrom, int[][]gatc){
		File GATC = new File (saveDirectory, "GATC");
		if (GATC.exists() == false) GATC.mkdir();
		String[] bases = {"G", "A", "T", "C"};
		//for each base
		for (int i=0; i< 4; i++){
			//make Point[]
			ArrayList<Point> pointsAL = new ArrayList<Point>();
			for (int x=0; x< gatc[i].length; x++) if (gatc[i][x] !=0) pointsAL.add(new Point(x,gatc[i][x]));
			if (pointsAL.size() == 0) continue;
			Point[] points = new Point[pointsAL.size()];
			pointsAL.toArray(points);
			//make notes
			HashMap <String,String> notes = new HashMap <String,String> ();
			notes.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			notes.put(BarParser.SOURCE_TAG, IO.concatinateFileFullPathNames(dataFiles, ","));
			notes.put(BarParser.UNIT_TAG, "Eland extended base quality score sum");
			notes.put(BarParser.DESCRIPTION_TAG, bases[i]+" Track. Generated by running the ElandSeqSumParser on Solexa ELAND file(s), interbase coordinates");
			//make an Info object  public Info (String text, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
			Info info = new Info(chrom, versionedGenome, chrom, ".", 0, notes);
			//make pd
			PointData pd = Point.extractPositionScores(points);			
			pd.setInfo(info);
			//write to file
			File barFile = new File (GATC, bases[i]+"_"+info.getChromosome()+".bar");
			pd.writePointDataToFile(barFile);
			//System.out.println("\n**********************");
			//System.out.println(pd.getInfo());
			//System.out.println("\tNum Obs "+pd.getNumberObservations());
			//Misc.printArray(pd.getPositions());
			//Misc.printArray(pd.getScores());
			//cleanup
			pd = null;
			pointsAL = null;
			points = null;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ElandSequenceParser(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		File[] genomeFiles = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[i+1]); i++; break;
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 'g': genomeFiles = IO.extractFiles(new File(args[i+1]), ".fasta"); i++; break;
					case 'a': minimumAlignmentScore = Float.parseFloat(args[i+1]); i++; break;
					case 'c': minimumConsensusScore = Float.parseFloat(args[i+1]); i++; break;
					case 'r': saveDirectory = new File (args[i+1]); i++; break;
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
		File[][] tot = new File[6][];
		tot[0] = IO.extractFiles(forExtraction,"_export.txt");
		tot[1] = IO.extractFiles(forExtraction,"_export.txt.zip");
		tot[2] = IO.extractFiles(forExtraction,"_sorted.txt");
		tot[3] = IO.extractFiles(forExtraction,"_sorted.txt.zip");
		tot[4] = IO.extractFiles(forExtraction,"_export.txt.gz");
		tot[5] = IO.extractFiles(forExtraction,"_sorted.txt.gz");
		dataFiles = IO.collapseFileArray(tot);
		
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx_sorted.txt(.zip/.gz) or xxx_export.txt(.zip/.gz) file(s)!\n");
		if (versionedGenome == null) Misc.printExit("\nPlease enter a genome version recognized by UCSC, see http://genome.ucsc.edu/FAQ/FAQreleases.\n");
		if (saveDirectory == null) saveDirectory = dataFiles[0].getParentFile();
		else if (saveDirectory.exists() == false) saveDirectory.mkdir();
		//parse genome files
		if (genomeFiles == null || genomeFiles.length ==0) Misc.printExit("\nError: cannot find any xxx.fasta genome files? Please complete option -g");
		fastas = new HashMap<String, File>();
		for (int i=0; i< genomeFiles.length; i++){
			int index = genomeFiles[i].getName().lastIndexOf('.');
			//anything found?
			if (index == -1) Misc.printExit("\nError: cannot parse chromosome text from genome file text, no '.'? "+genomeFiles[i]);
			String chromosome = genomeFiles[i].getName().substring(0, index);
			fastas.put (chromosome, genomeFiles[i]);
		}
		//make temp directory
		tempDirectory = new File(saveDirectory,"TempFiles_"+Passwords.createRandowWord(7));
		tempDirectory.mkdir();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Eland Sequence Parser: March 2009                        **\n" +
				"**************************************************************************************\n" +
				"Parses sequence information from Eland Extended alignment summary files. For every\n" +
				"base, sums the quality scores generating a G, A, T, and C track xxx.bar file for \n" +
				"visualization in IGB.  Also generates a consensus track (1-fraction consensus) for\n" +
				"each base.\n\n" +

				"-f The full path directory/file text of your xxx_export.txt(.zip/.gz) or\n" +
				"      xxx_sorted.txt(.zip/.gz) file(s).\n" +
				"-r Full path directory text for saving the results.\n"+
				"-g Full path directory text containing fasta files for reference base calling\n" +
				"      (e.g. chr1.fasta, chr5.fasta).\n"+
				"-v Versioned Genome (ie hg18, dm2, ce2, mm8), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-a Minimum aligment score, -10Log10(p-value), defaults to 13.\n"+
				"-c Minimum consensus score, -10Log10(p-value), defaults to 60.\n"+


				"\n\nExample: java -Xmx1500M -jar pathToUSeq/Apps/ElandSequenceParser -v hg18 -c 90\n" +
				"      -f /data/ExportFiles/ -r /data/Results -g /genomes/Hg18Fastas \n\n" +

		"**************************************************************************************\n");

	}	

}
