
package edu.utah.seq.parsers.mpileup;
import java.io.*;
import java.util.Random;
import java.util.regex.*;

import util.gen.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class MpileupRandomizer{
	//fields
	private File mpileupFile;
	private File scrambledFile;
	private int minReadDepth = 10;
	private int minNumSamples = 3;
	private long numPrintedLines = 0;
	private long numSkippedLines = 0;
	private long numBlocks = 0;
	private int minGap = 125;
	private String priorChrom = "";
	private int priorPos = 0;
	private int[] sampleIndexes = null;
	private Random random = new Random();
	private PrintWriter out;
	

	//constructors
	public MpileupRandomizer(String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			
			processArgs(args);
			IO.p("Thresholds:");
			IO.p(minReadDepth+"\tMinimum sample read depth");
			IO.p(minNumSamples+"\tMinimum num passing samples to save");
			IO.p(minGap+"\tMinimum bp gap to trigger sample reordering\n");
			
			parse();

			IO.p(numPrintedLines+"\tLines saved\n"+numSkippedLines+"\tLines failing thresholds\n"+numBlocks+"\tRandomized blocks");

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.p("\nDone! "+Math.round(diffTime)+" Min\n");
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public void checkSamples(String line) throws Exception{
	
		String[] fields = Misc.TAB.split(line);
		
		//skip Ns
		if (fields[2].equals("N")) {
			numSkippedLines++;
			priorChrom = fields[0];
			priorPos = Integer.parseInt(fields[1]);
			return;
		}
		
		//watch out for lines where last is 0 obs and no bases or quals so skip it
		int numFields = fields.length;
		int fieldsToParse = numFields-3;
		if (fieldsToParse % 3 !=0) {
			numFields--;
			fieldsToParse--;
			//other junkers
			if (fieldsToParse % 3 !=0) {
				System.err.println("Malformed, skipping: "+line);
				numSkippedLines++;
				priorChrom = fields[0];
				priorPos = Integer.parseInt(fields[1]);
				return;
			}
		}
		
		//parse samples
		int numSamples = fieldsToParse/3;
		
		//create sampleIndexes?
		if (sampleIndexes == null) createSampleIndexes(numSamples);
		
		//check num samples
		if (numSamples != sampleIndexes.length) throw new Exception ("Different number samples found in this line, aborting. \n"+line);
		
		String[] samples = new String[numSamples];
		numSamples = 0;
		int numPassing = 0;
		for (int i=3; i< numFields; i+=3) {
			samples[numSamples++] = new String(fields[i]+"\t"+fields[i+1]+"\t"+ fields[i+2]);
			int depth = Integer.parseInt(fields[i]);
			if (depth >= minReadDepth) numPassing++;
		}
		
		//scramble it?
		boolean scramble = false;
		int currPos = Integer.parseInt(fields[1]);
		//diff chroms?
		if (priorChrom.equals(fields[0]) == false) {
			scramble = true;
			priorChrom = fields[0];
		}
		else {
			//same chrom enough of a gap?
			int diff = currPos - priorPos;
			if (diff >= minGap) scramble = true;
		}
		priorPos = currPos;
		
		if (scramble) {
			Misc.randomize(sampleIndexes, random);
			numBlocks++;
		}
		
		//print it?
		if (numPassing >= minNumSamples){
		
			//print chr pos ref
			StringBuilder sb = new StringBuilder();
			sb.append(fields[0]);
			sb.append("\t");
			sb.append(fields[1]);
			sb.append("\t");
			sb.append(fields[2]);
		
			//print scrambled samples
			for (int i=0; i< sampleIndexes.length; i++){
				sb.append("\t");
				sb.append(samples[sampleIndexes[i]]);
			}
			out.println(sb);
			numPrintedLines++;
		}
		else numSkippedLines++;
	}

	private void createSampleIndexes(int numSamples) {
		sampleIndexes = new int[numSamples];
		for (int i=0; i< numSamples; i++) sampleIndexes[i] = i;
	}

	/**The goal here is to identify bases in non problematic regions (good read depth, low indels, low poor qual bps, ref not N) two good bp flanks, so working xxTxx */
	private void parse() {
		try{
			BufferedReader in = IO.fetchBufferedReader(mpileupFile);
			String line;
			out = new PrintWriter( new FileWriter(scrambledFile));

			while ((line=in.readLine())!= null){
				if (line.startsWith("#") || line.trim().length() == 0) {
					out.println(line);
				}
				else checkSamples(line);
			}
			in.close();
			out.close();
		}catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem paring mpileup file! Aborting.");
		}

	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.p("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': mpileupFile = new File(args[++i]); break;
					case 'r': minReadDepth = Integer.parseInt(args[++i]); break;
					case 's': minNumSamples = Integer.parseInt(args[++i]); break;
					case 'g': minGap = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (mpileupFile == null || mpileupFile.exists() == false) Misc.printErrAndExit("\nERROR: Can't find your mpileup file? "+mpileupFile);
		
		//make output file
		String name = mpileupFile.getName();
		if (name.endsWith(".gz")) name = name.replaceFirst(".gz", "");
		name = name +"_DP"+minReadDepth+"MS"+minNumSamples +".txt";
		scrambledFile = new File(mpileupFile.getParentFile(), name);

		
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MpileupRandomizer(args);
	}


	public static void printDocs(){
		IO.p("\n" +
				"**************************************************************************************\n" +
				"**                             Mpileup Randomizer: March 2018                       **\n" +
				"**************************************************************************************\n" +
				"Upon finding a gap in the coverage, the sample order is randomized and maintained. Use\n"+
				"this app to 'de-identify' a multi sample mpileup file while maintaining INDEL blocks.\n"+

				"\nRequired Options:\n"+
				"-m Path to a Samtools mpileup file (gz/zip OK).\n"+

				"\nDefault Options:\n"+
				"-r Minimum read depth to pass a sample, default 10\n"+
				"-s Minimum number of samples that must pass to save line, default 3\n"+
				"-g Minimum gap, defaults to 75\n"+

				"\nExample: java -Xmx4G -jar pathToUSeq/Apps/MpileupRandomizer -m normExo.mpileup.gz\n" +
				"     -o randomized.mpileup.gz -r 20 -s 4 \n" +

				"**************************************************************************************\n");

	}
}
