package edu.utah.seq.analysis.tele;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import util.bio.annotation.ExonIntron;
import util.bio.parsers.UCSCGeneLine;
import util.gen.Gzipper;
import util.gen.Misc;
import util.gen.Num;

public class TeleTranscript {
	private UCSCGeneLine transcript;
	private TeleStats treatment;
	private TeleStats control;
	private int cDNALength = 0;
	
	public TeleTranscript (UCSCGeneLine transcript, int numberExonicTreatmentAlignments, int numberExonicControlAlignments){
		this.transcript = transcript;
		treatment = new TeleStats();
		control = new TeleStats();
		treatment.setNumberExonicAlignments(numberExonicTreatmentAlignments);
		control.setNumberExonicAlignments(numberExonicControlAlignments);
	}
	
	public String toStringTabLine(double windowSize, double fractionLengthBackground){
		StringBuilder sb = new StringBuilder();
		//igb hyper link
		sb.append("=HYPERLINK(\"Genes/");
		sb.append(transcript.getDisplayName());
		sb.append("/");
		sb.append(transcript.getName());
		sb.append("_");
		sb.append(transcript.getDisplayName());
		sb.append("_exonic.png\",\"");
		sb.append(transcript.getName());
		sb.append("\")");
		sb.append("\t");
		//length cDNA
		sb.append(cDNALength); sb.append("\t");
		//5' best window start index
		sb.append(treatment.getFivePrimeWindowIndex()); sb.append("\t");
		//add median info
		treatment.addMedianInfo(sb); sb.append("\t");
		control.addMedianInfo(sb);sb.append("\t");
		double ratio = treatment.getLog2MedianSkew() - control.getLog2MedianSkew();
		sb.append(ratio);sb.append("\t");
		//add count info
		treatment.addCountInfo(sb, windowSize, cDNALength, fractionLengthBackground); sb.append("\t");
		control.addCountInfo(sb, windowSize, cDNALength, fractionLengthBackground);sb.append("\t");
		ratio = treatment.getLog2CountSkew(windowSize, cDNALength, fractionLengthBackground) - control.getLog2CountSkew(windowSize, cDNALength, fractionLengthBackground);
		sb.append(ratio);sb.append("\t");
		//pAdj count skew
		sb.append(treatment.getpAdjSkewedReadCount());
		return sb.toString();
	}
	
	public void calculateSkewStats(float fractionLengthForBackground, int windowSize){
		loadMedians(fractionLengthForBackground, windowSize, treatment, -1);
		loadMedians(fractionLengthForBackground, windowSize, control, (int)(treatment.getFivePrimeWindowIndex()));
	}
	
	/**Calculates the mean read coverage for the 5' and 3' ends using the windowSize for the 5'.*/
	public void loadMedians(float fractionLengthForBackground, int windowSize, TeleStats ts, int index){
		//create spliced exonic read coverage transcript
		float[] exonicReadCoverage = fetchSplicedReadCoverage(ts);
		
//System.out.println(transcript.getNames("_"));
//Misc.printArray(tTransRC);
		//calculate background
		float[] backgroundMedianIndex = calculateMedianBackground(fractionLengthForBackground, exonicReadCoverage);
		//calc best 5' window and window index? or used one found in prior scan?
		float median5Window;
		int bestIndex;
		if (index >= 0) {
			median5Window = calculateFixedWindowMedian(exonicReadCoverage, windowSize, index);
			bestIndex = index;
//System.out.println("Fixed: "+median5Window+" "+bestIndex+"\n"+Misc.floatArrayToString(exonicReadCoverage, "\t"));
		}
		else {
			float[] s = scan5PrimeEnd(exonicReadCoverage, windowSize, fractionLengthForBackground);
			median5Window = s[0];
			bestIndex = (int) s[1];
		}
		ts.setFivePrimeMedianBaseCoverage(median5Window);
		ts.setFivePrimeWindowIndex(bestIndex);
		ts.setThreePrimeMedianBaseCoverage(backgroundMedianIndex[0]);
		ts.setReadCountsForBestAndBackground(windowSize, (int)backgroundMedianIndex[1]);
//System.out.println("Final 5' and 3' "+ts.getFivePrimeMedianBaseCoverage()+"  "+ts.getThreePrimeMedianBaseCoverage());
//System.out.println("Best 5' window indexs "+ts.getFivePrimeWindowIndex()+" - "+(ts.getFivePrimeWindowIndex()+windowSize)+"\n");
	}

	private float calculateFixedWindowMedian(float[] rc, int windowSize, int index) {
		float[] slice = new float[windowSize];
		System.arraycopy(rc, index, slice, 0, windowSize);
		Arrays.sort(slice);		
		return (float)Num.median(slice);
	}
	
	/**Returns maxMedian and the bestIndex */
	private float[] scan5PrimeEnd (float[] rc, int windowSize, float fractionLengthForBackground){
		int num = Math.round( (float)rc.length * fractionLengthForBackground );
		int stop = rc.length - num;
		float[] toScan = new float[stop];
		System.arraycopy(rc, 0, toScan, 0, stop);
		float[] medians = Num.windowMedianScores(toScan, windowSize);
		int bestIndex = 0;
		float maxMedian = medians[0];
		for (int i=1; i< medians.length; i++){
			if (medians[i]> maxMedian){
				maxMedian = medians[i];
				bestIndex = i;
			}
		}
		return new float[]{maxMedian, (float)bestIndex};
	}

	/**Returns the median of the slice and the start index.*/
	private float[] calculateMedianBackground(float fractionLengthForBackground,float[] readCoverage) {
		int num = Math.round( (float)readCoverage.length * fractionLengthForBackground );
		int startIndex = readCoverage.length - num;
		float[] slice = new float[num];
		System.arraycopy(readCoverage, startIndex, slice, 0, num);
		Arrays.sort(slice);
		float median = (float)Num.median(slice);
//System.out.println("Background calcs: bpNum of end "+num+"\tStartIndex "+startIndex+"\tLen "+readCoverage.length+"\t");
		return new float[]{median, (float)startIndex};
	}

	private float[] fetchSplicedReadCoverage(TeleStats ts) {
		//calc length of cDNA
		cDNALength = transcript.getTotalExonicBasePairs();
		float[] rc = new float[cDNALength];
		ArrayList<String>[] rcNamesExonic = new ArrayList[cDNALength];
		float[] rcGenic = ts.getBaseCoverage();
		ArrayList<String>[] rcNamesGenic = ts.getBaseCoverageNames();
		ExonIntron[] exons = transcript.getExons();
		int firstBase = transcript.getTxStart();
		int rcIndex = 0;
		//for each exon
		for (int i=0; i< exons.length; i++){
			int start = exons[i].getStart() - firstBase;
			int end = exons[i].getEnd() - firstBase;
			for (int j = start; j < end; j++){
				rc[rcIndex] = rcGenic[j];
				rcNamesExonic[rcIndex] = rcNamesGenic[j];
				rcIndex++;
			}
		}
		//invert it for - strand genes
		if (transcript.isMinusStrand()) {
			Num.invertArray(rc);
			Misc.invertArray(rcNamesExonic);
		}

		ts.setExonicBaseCoverage(rc);
		ts.setExonicBaseCoverageNames(rcNamesExonic);
		return rc;
	}

	public String toString(float minimumBaseReadCoverage, boolean printGraphs){
		StringBuilder sb = new StringBuilder();
		sb.append(transcript.toString());
		sb.append("\n");
		sb.append(treatment.getNumberExonicAlignments()+"\t# exonic alignments for treatment");
		sb.append("\n");
		sb.append(control.getNumberExonicAlignments()+"\t# exonic alignments for control");
		sb.append("\n");
		sb.append(Num.intArrayToString(treatment.getTypes(), "\t"));
		sb.append("\tTreatment type count (Exon, SplicedExon, Intron, MisSpliced)");
		sb.append("\n");
		sb.append(Num.intArrayToString(control.getTypes(), "\t"));
		sb.append("\tControl type count (Exon, SplicedExon, Intron, MisSpliced)");
		if (printGraphs) {
			sb.append("\nTreatment base coverage\n");
			sb.append(Num.floatArrayToString(treatment.getBaseCoverage(), "\t"));
			sb.append("\nControl base coverage\n");
			sb.append(Num.floatArrayToString(control.getBaseCoverage(), "\t"));
			sb.append("\nNormalized Log2Rtos\n");
			sb.append(Num.floatArrayToString(getScaledRatios(minimumBaseReadCoverage, false), "\t"));
			sb.append("\n");
			sb.append("\nTreatment misSplice graph\n");
			sb.append(Num.floatArrayToString(getMisSpliceGraph(true), "\t"));
			sb.append("\n");
			sb.append("\nControl misSplice graph\n");
			sb.append(Num.floatArrayToString(getMisSpliceGraph(false), "\t"));
		}
		sb.append("\n");
		return sb.toString();
	}
	
	public void printSgrGraphs(File directory, float minimumBaseReadCoverage){
		File saveDir = new File (directory, transcript.getDisplayName());
		saveDir.mkdirs();
		
		File f;
		
		//treatment bc with introns
		f = new File (saveDir, transcript.getNames("_")+"_TreatmentBC.sgr.gz");
		writeSgr(f, treatment.getBaseCoverage(), true);
		//control bc with introns
		f = new File (saveDir, transcript.getNames("_")+"_ControlBC.sgr.gz");
		writeSgr(f, control.getBaseCoverage(), true);

		//log2Rto with introns
		f = new File (saveDir, transcript.getNames("_")+"_NormLog2Rto.sgr.gz");
		writeSgr(f, getScaledRatios(minimumBaseReadCoverage, false), true);
		
		//misSplice Treatment
		f = new File (saveDir, transcript.getNames("_")+"_TreatmentMisSplice.sgr.gz");
		float[] s = getMisSpliceGraph(true);
		if (s!= null) writeSgr(f, s, true);
		
		//misSplice Control
		f = new File (saveDir, transcript.getNames("_")+"_ControlMisSplice.sgr.gz");
		s = getMisSpliceGraph(false);
		if (s!= null) writeSgr(f, s, true);
		
		//treatment exonic bc
		f = new File (saveDir, transcript.getNames("_")+"_TreatmentExonBC.sgr.gz");
		writeSgr(f, treatment.getExonicBaseCoverage(), false);
		//control exonic bc
		f = new File (saveDir, transcript.getNames("_")+"_ControlExonBC.sgr.gz");
		writeSgr(f, control.getExonicBaseCoverage(), false);
		//log2Rto with introns
		f = new File (saveDir, transcript.getNames("_")+"_ExonNormLog2Rto.sgr.gz");
		writeSgr(f, getScaledRatios(minimumBaseReadCoverage, true), false);
	}

	private void writeSgr(File f, float[] vals, boolean setFirstBase) {
		try {
			Gzipper out = new Gzipper(f);
			int firstBase = 0;
			String chr;
			if (setFirstBase) {
				firstBase = transcript.getTxStart();
				chr = transcript.getChrom() +"\t";
			}
			else {
				chr = transcript.getNames("_")+"\t";
			}
			for (int i=0; i< vals.length; i++){
				if (vals[i] != 0) {
					int pos = firstBase + i;
					out.println(chr+ pos+ "\t"+ vals[i]);
				}
			}
			out.close();
		} catch (Exception e) {
			System.err.println("Problem writing "+f);
			e.printStackTrace();
		} 
	}
	
	public void printExonicGGPlot(File directory, StringBuilder rScript) {
		try {
			File saveDir = new File (directory, transcript.getDisplayName());
			saveDir.mkdirs();

			//generate table
			File table = new File (saveDir, transcript.getNames("_")+"_ggplot2Data.txt");
			table.deleteOnExit();
			PrintWriter out = new PrintWriter( new FileWriter (table));
			out.println("Group\tPosition\tRelativeCoverage");
			
			//find max
			float[] t = treatment.getExonicBaseCoverage();
			float[] c = control.getExonicBaseCoverage();
			float maxT = Num.findHighestFloat(t);
			float maxC = Num.findHighestFloat(c);
			//print
			for (int i=0; i< t.length; i++){
				if (t[i] == 0.0f) continue;
				float val = t[i]/maxT;
				out.println("T\t"+i+"\t"+val);
			}
			for (int i=0; i< c.length; i++){
				if (c[i] == 0.0f) continue;
				float val = c[i]/maxC;
				out.println("C\t"+i+"\t"+val);
			}
			//append R script
			File png = new File (saveDir, transcript.getNames("_")+"_Exonic.png");
			rScript.append("png('"+png+"', height=400, width=1200)\n");
			rScript.append("dfn = read.table(header=T, file='"+table+"')\n");
			rScript.append("ggplot(data=dfn, aes(x=Position, y=RelativeCoverage, group=Group, colour=Group)) + geom_line() + geom_point()\n");
			rScript.append("dev.off()\n");
			
			
			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public float[] getMisSpliceGraph(boolean isTreatment){	
		float[] counts = new float[treatment.getBaseCoverage().length];
		int[] pos;
		if (isTreatment) pos = treatment.getMisSplicePositions();
		else pos = control.getMisSplicePositions();	
		if (pos == null) return null;
		for (int i=0; i< pos.length; i++) counts[pos[i]]++;
		return counts;
	}

	//calculates log2(numT+1/numC+1) for bases passing the minimumBaseReadCoverage threshold
	public float[] getScaledRatios(float minimumBaseReadCoverage, boolean useExonicRC) {
		//scalar to multiply C or divide T
		float scalar = (float)(treatment.getNumberExonicAlignments()+1) / (float)(control.getNumberExonicAlignments()+1);
		//for each base
		float[] tbc;
		float[] cbc;
		float[] ratios;
		if (useExonicRC){
			tbc = treatment.getExonicBaseCoverage();
			cbc = control.getExonicBaseCoverage();
			ratios = new float[tbc.length];
		}
		else {
			tbc = treatment.getBaseCoverage();
			cbc = control.getBaseCoverage();
			ratios = new float[tbc.length];
		}
		for (int i=0; i< ratios.length; i++){
			//fails base read coverage? if so the defaults to zero
			if ((tbc[i] + cbc[i]) < minimumBaseReadCoverage) continue;
			//divide T?
			if (tbc[i] >= cbc[i]) {
				float t = tbc[i]/ scalar;
				ratios[i] = Num.log2( (t+1.0f)/ (cbc[i]+1.0f));
			}
			else {
				float c = cbc[i] * scalar;
				ratios[i] = Num.log2( (tbc[i]+1.0f)/ (c+1.0f));
			}
		}
		return ratios;
	}

	public UCSCGeneLine getTranscript() {
		return transcript;
	}

	public void setTranscript(UCSCGeneLine transcript) {
		this.transcript = transcript;
	}

	public TeleStats getTreatment() {
		return treatment;
	}

	public void setTreatment(TeleStats treatment) {
		this.treatment = treatment;
	}

	public TeleStats getControl() {
		return control;
	}

	public void setControl(TeleStats control) {
		this.control = control;
	}

	

}
