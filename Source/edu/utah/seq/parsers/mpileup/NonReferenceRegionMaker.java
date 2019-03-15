
package edu.utah.seq.parsers.mpileup;
import java.io.*;
import java.util.regex.*;

import util.gen.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class NonReferenceRegionMaker{
	//fields
	private File mpileupFile;
	private File bedFile;
	private int minReadCoverage = 10;
	private int minNonRefCount = 3;
	private double minNonRefAlleleFreq = 0.05;
	private int minBaseQuality = 13;

	//internal
	private long numParsedBps = 0;
	private long numPassReadDeptBps = 0;
	private long numSavedBps = 0;
	private Gzipper outBed = null;
	private Histogram histAFs = new Histogram(0.0, 1.01, 101);

	//constructors
	public NonReferenceRegionMaker(String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);
			printSettings();

			IO.p("\nParsing");
			parse();
			outBed.close();

			IO.pl("\nStats:\n"+numParsedBps+"\tParsed mpileup records/ BPs");
			IO.pl(numPassReadDeptBps+"\tBPs passing read depth");
			IO.pl(numSavedBps+"\tBPs passing allele frequency and saved to file\n");
			IO.pl("Histogram of non zero AFs:");
			histAFs.printScaledHistogram();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Min\n");
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	private void parse() {
		try{
			BufferedReader in;
			if (mpileupFile != null) in = IO.fetchBufferedReader(mpileupFile);
			else in = new BufferedReader( new InputStreamReader(System.in) );
			
			String line;
			int counter = 0;
			while ((line=in.readLine())!= null){
				if (line.startsWith("#") || line.trim().length() == 0) continue;
				if (counter++ > 1000000) {
					IO.p(".");
					counter = 0;
				}
				MpileupLine ml = new MpileupLine(line, minBaseQuality);
				numParsedBps++;
				
				if (ml.getChr() == null) continue;
				if (ml.getSamples().length != 1) Misc.printErrAndExit("\nNot just one sample? See "+line);
				saveCounts(ml);
			}
			IO.pl();
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem paring mpileup file! Aborting."+mpileupFile);
		}

	}

	private void saveCounts(MpileupLine ml) throws IOException {
		MpileupSample s = ml.getSamples()[0];
		if (s.getMpileupReadCount() < minReadCoverage) return;
		numPassReadDeptBps++;
		
		int snvCount = s.getNonRefBaseCounts();
		int insertions = s.getInsertions();
		int deletions = s.getDeletions();
		if ((snvCount+insertions+deletions) < minNonRefCount) return;

		double nonRefAF = s.getAlleleFreqNonRefPlusIndels();
		if (nonRefAF > 0.0) histAFs.count(nonRefAF);
		if (nonRefAF < minNonRefAlleleFreq) return;

		//save it
		numSavedBps++;
		outBed.printTab(ml.getChr());
		outBed.printTab(Integer.toString(ml.getZeroPos()));
		if (insertions == 0) outBed.printTab(Integer.toString(ml.getZeroPos()+1));
		else outBed.printTab(Integer.toString(ml.getZeroPos()+2));

		//Read Depth
		outBed.printTab(Integer.toString(s.getReadCoverageAll()));

		//AF
		outBed.printTab(Num.formatNumber(nonRefAF, 4));

		outBed.println("\t.");
	}


	/**This method will process each argument and assign new variables
	 * @throws IOException 
	 * @throws FileNotFoundException */
	public void processArgs(String[] args) throws FileNotFoundException, IOException{
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': mpileupFile = new File(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'r': minReadCoverage = Integer.parseInt(args[++i]); break;
					case 'c': minNonRefCount = Integer.parseInt(args[++i]); break;
					case 'a': minNonRefAlleleFreq = Double.parseDouble(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (mpileupFile != null && mpileupFile.exists() == false) Misc.printErrAndExit("\nERROR: Can't find your mpileup file? "+mpileupFile);
		if (bedFile == null ) Misc.printErrAndExit("\nERROR: Please provide a file path for writing the output bed file, should end in xxx.bed.gz "+bedFile);

		outBed = new Gzipper(bedFile);
		outBed.println("#Chr\tStart\tStop\tReadDepth\tNonRefAlleleFreq\t.");
	}

	public void printSettings(){
		if (mpileupFile != null) IO.pl("Parsing\t"+mpileupFile);
		else IO.pl("Parsing\tStandard In");
		IO.pl("Non Ref Bps\t"+bedFile);
		IO.pl(minReadCoverage+"\tMin Mpileup Read Coverage");
		IO.pl(minNonRefCount+"\tMinimum Non Ref Count");
		IO.pl(minNonRefAlleleFreq+"\tMinimum Non Ref Allele Freq");
		IO.pl(minBaseQuality+"\tMin Base Quality");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new NonReferenceRegionMaker(args);
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Non Reference Region Maker: Jan 2018                   **\n" +
				"**************************************************************************************\n" +
				"NRRM scans a single sample mpileup file looking for non reference base pairs. If these\n"+
				"pass read depth, allele frequency, and non ref base count thresholds, the base is\n"+
				"written to a bed file. BPs with insertions are saved as a 2 BP region. Run MergeRegions\n"+
				"or MergeAdjacentRegions to join proximal non ref BPs. \n"+

				"\nOptions:\n"+
				"-m Provide a path to a single sample samtools mpileup file or pipe mpileup output.\n"+
				"-b Path to write the bed file output, should end in xxx.bed.gz\n"+
				"-r Minimum read depth, 10\n"+
				"-a Minimum non reference allelic frequency (SNVs + INDELS), default 0.05\n"+
				"-c Minimum non reference base count, default 3\n"+
				"-q Minimum base quality for inclusion in AF calculation, default 10\n"+


				"\nExample: samtools mpileup -B -d 1000000 -f $faIndex -l $bed $bam | java\n"+
				"     -Xmx4G -jar pathToUSeq/Apps/NonReferenceRegionMaker -q 13 -r 20 -a 0.025 -b \n"+
				"     0.025normExoNonRefMask.bed.gz -c 4 \n" +

				"**************************************************************************************\n");

	}
}
