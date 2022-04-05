package edu.utah.seq.qc;

import java.io.BufferedReader;
import java.io.File;
import org.json.JSONArray;
import org.json.JSONObject;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class DnaSample {

	//Cutadapt -> Align -> mark dup -> filt
	private boolean alignParsed = false;
	private long numberFastqRead1And2;

	//Markdup
	private boolean dupParsed = false;
	private long numberUnfilteredAlignments;
	private long numberFilteredAlignments;
	private long numberDuplicateAlignments;

	//Readcov
	private boolean readCovParsed = false;
	private double meanOnTargetCoverage;
	private int coverageAt090OfTargetBps;
	private int coverageAt095OfTargetBps;
	private float[] fractionTargetBpsWithIndexedCoverage;

	public DnaSample(File alignLog, File dupLog, File readCovJson) throws Exception {
		parseAlignLog(alignLog);
		parseDupLog(dupLog);
		parseReadCoverageJson(readCovJson);
	}
	
	private void parseReadCoverageJson(File readCovJson) throws Exception {
		if (readCovJson == null) return;
		JSONObject jo = new JSONObject(IO.loadFile(readCovJson, " ", true));
		meanOnTargetCoverage = Double.parseDouble(jo.getString("meanCoverage"));
		coverageAt095OfTargetBps = jo.getInt("coverageAt0.95OfTargetBps");
		coverageAt090OfTargetBps = jo.getInt("coverageAt0.90OfTargetBps");
		JSONArray ja = jo.getJSONArray("fractionTargetBpsWithIndexedCoverage");
		int num = ja.length();
		fractionTargetBpsWithIndexedCoverage = new float[num];
		for (int i=0; i< num; i++) fractionTargetBpsWithIndexedCoverage[i] = (float)ja.getDouble(i);		
		readCovParsed = true;
	}

	private void parseDupLog(File dupLog) throws Exception {
		if (dupLog == null) return;
		String line;
		BufferedReader in = IO.fetchBufferedReader(dupLog);
		while ((line = in.readLine())!=null) {
			if (line.startsWith("READ:")) {
				String[] tokens = Misc.COLON.split(line);
				String num = tokens[1].trim();
				numberUnfilteredAlignments = Long.parseLong(num);
			}
			else if (line.startsWith("WRITTEN:")) {
				String[] tokens = Misc.COLON.split(line);
				String num = tokens[1].trim();
				numberFilteredAlignments = Long.parseLong(num);
			}
			else if (line.startsWith("DUPLICATE TOTAL:")) {
				String[] tokens = Misc.COLON.split(line);
				String num = tokens[1].trim();
				numberDuplicateAlignments = Long.parseLong(num);
				in.close();
				dupParsed = true;
				return;
			}
		}
	}

	private void parseAlignLog(File alignLog) throws Exception {
		if (alignLog == null) return;
		String line;
		BufferedReader in = IO.fetchBufferedReader(alignLog);
		while ((line = in.readLine())!=null) {
			if (line.startsWith("Total read pairs processed:")) {
				String[] tokens = Misc.COLON.split(line);
				String num = tokens[1].trim();
				num = Misc.COMMA.matcher(num).replaceAll("");
				numberFastqRead1And2 = Long.parseLong(num) * 2;
				in.close();
				alignParsed = true;
				return;
			}
		}
		
	}
	
	/**At what index does the fraction hit or fall below 0.25?*/
	public int whenHit25ReadCoverage() {
		for (int i=0; i< fractionTargetBpsWithIndexedCoverage.length; i++){
			if (fractionTargetBpsWithIndexedCoverage[i] <= 0.25) return i;
		}
		return fractionTargetBpsWithIndexedCoverage.length-1;
	}
	
	public static String getHeader() {
		StringBuilder sb = new StringBuilder();
		sb.append("# Fastq Reads"); //numberFastqRead1And2;
		sb.append("\t# Unfiltered Alignments"); //numberUnfilteredAlignments
		sb.append("\tFraction Duplicate"); //numberDuplicateAlignments/numberUnfilteredAlignments
		sb.append("\tFraction Saved"); //numberFilteredAlignments/numberUnfilteredAlignments
		sb.append("\t# Saved Alignments"); //numberFilteredAlignments

		sb.append("\tMean on Target Coverage"); //meanOnTargetCoverage
		sb.append("\tCoverage at 0.9 of Target BPs"); //coverageAt090OfTargetBps
		sb.append("\tCoverage at 0.95 of Target BPs"); //coverageAt095OfTargetBps
		return sb.toString();
	}

	public String getTabbedLine(String type) {
		StringBuilder sb = new StringBuilder();
		sb.append("\t");
		sb.append(type); sb.append("\t");
		
		if (alignParsed) {
			sb.append(Num.formatNumber(numberFastqRead1And2, 0)); 
			sb.append("\t");
		}
		else sb.append("NA\t");
		
		
		if (dupParsed) {
			sb.append(Num.formatNumber(numberUnfilteredAlignments,0)); sb.append("\t");
			double fracDup = (double)numberDuplicateAlignments/(double)numberUnfilteredAlignments; 
			sb.append(Num.formatNumber(fracDup, 2)); sb.append("\t");
			double fracSaved = (double)numberFilteredAlignments/(double)numberUnfilteredAlignments; 
			sb.append(Num.formatNumber(fracSaved, 2)); sb.append("\t");
			sb.append(Num.formatNumber(numberFilteredAlignments,0)); sb.append("\t");
		}
		else  sb.append("NA\tNA\tNA\tNA\t");
		
		if (readCovParsed) {
			sb.append(Num.formatNumber(meanOnTargetCoverage,1)); sb.append("\t");
			sb.append(coverageAt090OfTargetBps); sb.append("\t");
			sb.append(coverageAt095OfTargetBps); sb.append("\t");
		}
		else  sb.append("NA\tNA\tNA\t");
		
		return sb.toString();
	}

	public float[] getFractionTargetBpsWithIndexedCoverage() {
		return fractionTargetBpsWithIndexedCoverage;
	}

}
