package edu.utah.seq.parsers.mpileup;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**@author Nix
 * Calculates the mean and std dev of the poor quality bases for each 5mer context.  Used to put a z-score on the background. 
 * 
 * MUST use samtools 1.8, the only param you should change is the --min-MQ, any others and st will generate garbage
 * 
 samtools=/uufs/chpc.utah.edu/common/HIPAA/u0028003/BioApps/Samtools/1.8/bin/samtools
fasta=/uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa
bed=/scratch/mammoth/serial/u0028003/Underhill/SpikeBuild/VcfToInject/Bed/HSV1_GBM_IDT_Probes_Hg38Pad150bps_91K.bed
$samtools mpileup \
--no-BAQ \
--count-orphans \
--max-depth 100000 \
--max-idepth 100000 \
--gap-frac 0.001 \
--per-sample-mF \
--ignore-RG \
--min-MQ 13 \
--fasta-ref $fasta \
--ff UNMAP,SECONDARY,QCFAIL \
-l $bed \
-o 15352X1_0.01.mpileup \
15352X1_0.01.bam

 * */
public class MpileupStats {

	//user defined fields
	private File mpileupFile;
	private File statsFile;
	private int minBaseQuality = 20;
	private int minDepth = 1000;

	//internal fields
	private BufferedReader in = null;
	//4^5=1024 4^3=64, unique k-mers of 5 and 3 length
	private TreeMap<String, ArrayList<Float>> fiveMerNFreq = new TreeMap<String, ArrayList<Float>>();
	private HashMap<String, double[]> seqStats = new HashMap<String, double[]>();

	//constructor
	public MpileupStats(String[] args){
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);
			
			//create IO
			if (mpileupFile != null) in = IO.fetchBufferedReader(mpileupFile);
			else in = new BufferedReader(new InputStreamReader(System.in));

			parse();
			
			printMeanStdDev();
			
			IO.saveObject(statsFile, seqStats);

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("Problem encountered, aborting!");
		} finally {
			try {
				if (in != null) in.close();
			} catch (IOException e) {}}
	}


	private void printMeanStdDev() {
		IO.pl("Seq\tObs\tMean\tStd");
		for (String seq : fiveMerNFreq.keySet()) {
			ArrayList<Float> nAFs = fiveMerNFreq.get(seq);
			float[] n = Num.arrayListOfFloatToArray(nAFs);
			double mean = Num.mean(n);
			double std = Num.standardDeviation(n, mean);
			IO.pl(seq +"\t"+ n.length + "\t" + mean+ "\t" +std);
			seqStats.put(seq, new double[] {mean,std});
		}
		
	}


	public void parse() throws Exception{

		String line = null;
		
		ArrayList<MpileupSample> block = new ArrayList<MpileupSample>();
		while ((line = in.readLine())!= null) {
			if (line.startsWith("#") || line.trim().length() == 0) continue;
			else {
				MpileupLine ml = new MpileupLine(line, minBaseQuality);
				MpileupSample currentSample = ml.getSamples()[0];
				if (block.size() ==0) block.add(currentSample);
				else {
					//check chrom and position
					MpileupLine priorML = block.get(block.size()-1).getRecord();
					if (priorML.getChr().equals(ml.getChr()) == false || priorML.getZeroPos() != (ml.getZeroPos()-1)) {
						processBlock(block);
						block.clear();
						block.add(currentSample);
					}
					else block.add(currentSample);
				}
			}
		}
		//process last		
		processBlock(block);
	}

	private void processBlock(ArrayList<MpileupSample> block) throws IOException {

		//fetch N AFs
		float[] nAFs = new float[block.size()];
		String[] bases = new String[nAFs.length];
		for (int i=0; i< nAFs.length; i++) {
			MpileupSample ms = block.get(i);

			//pass read depth?
			if ((ms.getReadCoverageAll()+ms.getPoorQualBases()) < minDepth) nAFs[i] = -1f;
			else nAFs[i] = (float)ms.getAlleleFreqNs();
			bases[i]  = ms.getRecord().getRef();
		}

		//5mer
		int windowSize = 5;
		int nAFIndex = 2;
		int num = 1+ nAFs.length - windowSize;
		if (num > 0 ) {
			for (int i=0; i<num; i++){
				//does center base pass read depth
				float nAF = nAFs[i+nAFIndex];
				if (nAF != -1) {
					int stop = i+windowSize;
					StringBuilder seq = new StringBuilder();
					for (int j=i; j< stop; j++) seq.append(bases[j]);
					String s = seq.toString();
					
					ArrayList<Float> al = fiveMerNFreq.get(s);
					if (al == null) {
						al = new ArrayList<Float>();
						fiveMerNFreq.put(s, al);
					}
					al.add(nAF);
				}
			}
		}

	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MpileupStats(args);
	}		



	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': mpileupFile = new File(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'd': minDepth = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (mpileupFile == null || mpileupFile.canRead() == false) Misc.printExit("\nError: cannot find or read your mpileup file\n");
		String name = Misc.removeExtension(mpileupFile.getName());
		statsFile = new File(mpileupFile.getParent(), name+ "."+minBaseQuality +".mstats.obj");
		
	}	




	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                    Mpileup Stats                                 **\n" +
				"**************************************************************************************\n" +
				"Scans a samtools v1.8 mpileup file calculating 5mer base context info to aid in\n"+
				"flagging FPs. Produces an xxx.mstats.obj file for use in other USeq apps.  Example: \n"+
				
				"\nsamtools1.8 mpileup --no-BAQ --count-orphans --max-depth 100000 --max-idepth 100000 \\\n"+
				"--gap-frac 0.001 --per-sample-mF --ignore-RG --min-MQ 13 --fasta-ref /Hg38/hs38DH.fa \\\n"+
				"--ff UNMAP,SECONDARY,QCFAIL -l regionsToScan.bed -o 15352X1.mpileup 15352X1.bam \n"+
				
				"\nRequired:\n"+
				"-m Path to a mpileup file, gz/zip OK.\n" +

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MpileupStats -m tumorA73.mpileup.gz\n\n"+

				"**************************************************************************************\n");
	}

	//getters and setters
	public int getMinBaseQuality() {
		return minBaseQuality;
	}

}
