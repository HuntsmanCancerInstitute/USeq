package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;

import edu.utah.seq.parsers.mpileup.MpileupSample;
import edu.utah.seq.parsers.mpileup.MpileupTabixLoader;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class GraphiteZScore {
	
	private File vcfIn;
	private File graphiteExe;
	private File indexedFasta;
	private File[] normalBams;
	private File graphiteVcf;
	private String optionalBashCmds = null;
	private double minimumSampleReadCount;
	private int minimumNumSamples;
	private double maximumSampleAF;
	private boolean deleteTemp = false;
	private HashMap<String, GraphiteResult> vcfGraphiteResults = new HashMap<String, GraphiteResult>();

	/**The vcfToProcess must be stripped of any FORMAT and sample fields, only CHROM thru INFO.
	 * @throws IOException */
	public GraphiteZScore(File vcfToProcess, File graphiteExe, File[] normalBams, File indexedFasta, String optionalBashCmds, int minimumSampleReadCount, int minimumNumSamples, double maximumSampleAF) throws Exception {
		this.vcfIn = vcfToProcess;
		this.graphiteExe = graphiteExe;
		this.normalBams = normalBams;
		this.indexedFasta = indexedFasta;
		this.optionalBashCmds = optionalBashCmds;
		this.minimumSampleReadCount = minimumSampleReadCount;
		this.minimumNumSamples = minimumNumSamples;
		this.maximumSampleAF = maximumSampleAF;
		doWork();
	}

	private void doWork() throws Exception {
		//graphite -f indexFasta -v vcf -o GResults -t 20 -b bam1 -b bam2 -b bam3
		String name = Misc.removeExtension(vcfIn.getName());		
		File resultsDir = new File(vcfIn.getParentFile(), "GraphiteTemp_"+name);
		resultsDir.mkdirs();
		graphiteVcf = new File(resultsDir, vcfIn.getName());

		executeGraphite(resultsDir);

		parseGraphiteResults();
		
		if (deleteTemp) {
			graphiteVcf.delete();
			IO.deleteDirectory(resultsDir);
		}
		
	}
	
	private void parseGraphiteResults() throws Exception{
		
		BufferedReader in = IO.fetchBufferedReader(graphiteVcf);
		String line;
		String[] t;
		String[] format;
		String[] sample;
		String[] dp4;
		ArrayList<Double> afs = new ArrayList<Double>();
		while ((line = in.readLine())!= null){
			if (line.startsWith("#")) continue;
			t = Misc.TAB.split(line);
			 
			format = Misc.COLON.split(t[8]);
			//find DP4_NFP
			int index =0;
			for (;index<format.length; index++){
				if (format[index].equals("DP4_NFP")) break;
			}
			//for each sample
			afs.clear();
			for (int i=9; i< t.length; i++){
				sample = Misc.COLON.split(t[i]);
				dp4 = Misc.COMMA.split(sample[index]);
				if (dp4.length!=4) throw new IOException("Failed to parse sample info for "+line);
				double fRef = Double.parseDouble(dp4[0]);
				double rRef = Double.parseDouble(dp4[1]);
				double fAlt = Double.parseDouble(dp4[2]);
				double rAlt = Double.parseDouble(dp4[3]);
				double total = fRef+rRef+fAlt+rAlt;
				if (total < minimumSampleReadCount) continue;
				double sampleAF = (fAlt+rAlt)/total;
				if (sampleAF < maximumSampleAF) afs.add(sampleAF);
			}
			String oriVcf = combineFirst7(t);
			
			//calculate z score?
			if (afs.size() < minimumNumSamples) {
				vcfGraphiteResults.put(oriVcf, null);
			}
			else {
				double[] values = Num.arrayListOfDoubleToArray(afs);
				double mean = Num.mean(values);
				double stndDev = Num.standardDeviation(values,mean);
				Matcher mat = MpileupTabixLoader.AF.matcher(t[7]);
				if (mat.matches() == false) throw new IOException ("Failed to parse AF= number from the INFO field in this graphite processed variant:\n"+line);
				double freq = Double.parseDouble(mat.group(1));
				double z = (freq-mean)/stndDev;
				if (Double.isInfinite(z)) z = VCFBackgroundChecker.zscoreForInfinity;
				vcfGraphiteResults.put(oriVcf, new GraphiteResult(z, values));
			}
		}
		in.close();
	}
	
	private String combineFirst7(String[] t) {
		StringBuilder sb = new StringBuilder(t[0]);
		for (int i=1; i< 8; i++){
			sb.append("\t");
			sb.append(t[i]);
		}
		return sb.toString();
	}

	public class GraphiteResult {
		private double[] bkAFs;
		private double zScore;
		
		public GraphiteResult(double zScore, double[] bkAFs){
			this.bkAFs = bkAFs;
			this.zScore = zScore;
		}
		public double[] getBkAFs() {
			return bkAFs;
		}
		public double getzScore() {
			return zScore;
		}
	}

	private void executeGraphite(File resultsDir) throws IOException {
		StringBuilder bs = new StringBuilder();
		
		//any bash cmds? e.g.  module load gcc/4.9.2
		if (optionalBashCmds != null) {
			bs.append(optionalBashCmds);
			bs.append("\n");
		}
		bs.append(graphiteExe.getCanonicalPath());
		bs.append(" -f ");
		bs.append(indexedFasta.getCanonicalPath());
		bs.append(" -v ");
		bs.append(vcfIn.getCanonicalPath());
		bs.append(" -o ");
		bs.append(resultsDir.getCanonicalPath());
		bs.append(" -t ");
		bs.append(Runtime.getRuntime().availableProcessors());
		for (File b: normalBams){
			bs.append(" -b ");
			bs.append(b.getCanonicalPath());
		}
		bs.append("\n");		
		
		String[] messages = IO.executeShellScript(bs.toString(), resultsDir);
		if (messages == null || messages.length > 0) {
			if (messages != null) throw new IOException("\nProblem executing\n"+bs.toString()+"\n"+Misc.stringArrayToString(messages, "\n"));
			throw new IOException("\nProblem executing\n"+bs.toString()+"\n");
		}		
	}

	public HashMap<String, GraphiteResult> getVcfGraphiteResults() {
		return vcfGraphiteResults;
	}

	
}
