package edu.utah.seq.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;

import com.eclipsesource.json.JsonArray;
import com.eclipsesource.json.JsonObject;
import com.eclipsesource.json.JsonValue;
import util.gen.IO;
import util.gen.Misc;

public class SampleQC {

	//fields
	String sampleName;
	
	//Fastq
	private boolean fastqParsed = false;
	private long numberFastqReads;

	//SAE
	private boolean saeParsed = false;
	private long numberUnfilteredAlignments;
	private double fractionAlignmentsPassQCScoreFilters;
	private double fractionOnTargetAndPassQCScoreFilters;
	private double estimatedFractionDuplicateAlignments;
	private double mappingQualityThreshold;
	private double alignmentScoreThreshold;
	private boolean divideAlignmentScoreByCigarM;

	//MPA
	private boolean mpaParsed = false;
	private double fractionOverlappingBpsInPairedReads;
	private double meanInsertSize;
	private long numberPassingBps;
	private double fractionPassingQ20bps;
	private double fractionPassingQ30bps;

	//S2U
	private boolean s2uParsed = false;
	private double meanOnTargetCoverage;
	private int minimumCoverageThreshold;
	private double coverageAt095OfTargetBps;
	private float[] fractionTargetBpsWithIndexedCoverage;
	private String[] lowCoverageRegions; 
	private String[] exonicMedianPerBpCoverage;
	private long numberLowCoverageBps;

	//constructor
	public SampleQC( String sampleName){
		this.sampleName = sampleName;
	}

	@SuppressWarnings("deprecation")
	public void loadJson(File f, String type){

		try {
			//use own buffered reader to support gzipped files
			BufferedReader in = IO.fetchBufferedReader(f);
			JsonObject jo = JsonObject.readFrom(in);
			if (type.equals("fastq")){
				fastqParsed = true;
				//total number of fastq reads
				numberFastqReads = jo.getLong("numberFastqReads", -1);
			}
			else if (type.equals("sae")){
				saeParsed = true;
				//approx number of paired fastq off machine, plus some supplementary and secondary as well those that don't align or align's very poorly
				numberUnfilteredAlignments = jo.getLong("numberUnfilteredAlignments", -1);
				//fraction that actually pass basic quality and alignment uniqueness filters, ideally > 0.8
				fractionAlignmentsPassQCScoreFilters = jo.getDouble("fractionAlignmentsPassQCScoreFilters", -1);
				//fraction that also are on target, this represents the usable paired and often overlapping alignments, ideally > 0.75
				fractionOnTargetAndPassQCScoreFilters = jo.getDouble("fractionOnTargetAndPassQCScoreFilters", -1);
				//fraction of these usable alignments that are also marked as duplicate, all but one will be tossed at a later step, ideally < 0.2
				estimatedFractionDuplicateAlignments = jo.getDouble("estimatedFractionDuplicateAlignments", -1);
				//MQ threshold,
				mappingQualityThreshold = jo.getDouble("mappingQualityThreshold", -1);
				//AS threshold see boolean it might be divided the the length of the M's in the CIGAR String
				alignmentScoreThreshold = jo.getDouble("alignmentScoreThreshold", -1);
				divideAlignmentScoreByCigarM = jo.getBoolean("divideAlignmentScoreByCigarM", true);
			}
			else if (type.equals("mpa")){
				mpaParsed = true;
				//mean insert size, ideally big enough to have < 0.1 overlap
				meanInsertSize = jo.getDouble("meanInsertSize", -1);
				//fraction of bps that overlap each other in paired alignments, ideally < 0.1
				fractionOverlappingBpsInPairedReads = jo.getDouble("fractionOverlappingBpsInPairedReads", -1);
				//total bps or yield after collapsing overlapping alignments
				numberPassingBps = jo.getLong("numberPassingBps", -1);
				//fraction of the total with good quality, ideally > 0.95
				fractionPassingQ20bps = jo.getDouble("fractionPassingQ20bps", -1);
				//fraction of the total with very good quality, ideally > 0.9
				fractionPassingQ30bps = jo.getDouble("fractionPassingQ30bps", -1);
			}
			else if (type.equals("s2u")){
				s2uParsed = true;
				//the traditional view of coverage, mean per bp coverage over captured regions, note with all of these stats, overlapping alignments have been merged so no double counting
				meanOnTargetCoverage= jo.getDouble("meanOnTargetCoverage", -1);
				//this is the fraction of target bps with >=0x, >=1x, >=2x .....
				fractionTargetBpsWithIndexedCoverage = parseFloatArray(jo.get("fractionTargetBpsWithIndexedCoverage"));
				//the coverage # found at 0.95 of target bps, calculate by asking what fraction of target bps have at minimum, 0x 1x, 2x, 3x, ... stop when fraction hits 0.95.
				coverageAt095OfTargetBps = jo.getDouble("coverageAt0.95OfTargetBps", -1);
				//this is the median of the per bp coverage across each region supplied to S2U
				exonicMedianPerBpCoverage = parseStringArray(jo.get("exonicMedianPerBpCoverage"));
				//the threshold, typically 20-30x for germline, 300x for somatic
				minimumCoverageThreshold = jo.getInt("minimumCoverageThreshold", -1);
				//bps that fail the minimum coverage threshold
				lowCoverageRegions = parseStringArray(jo.get("lowCoverageRegions"));
				//number of target bps failing the coverage threshold
				numberLowCoverageBps = jo.getLong("numberLowCoverageBps", -1);
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR reading json file "+f);
		}
	}
	
	public String fetchTabbedLine() {
		ArrayList<String> al = new ArrayList<String>();
		al.add(sampleName);
		if (fastqParsed) al.add(new Long(numberFastqReads).toString());
		if (saeParsed){
			al.add(new Long(numberUnfilteredAlignments).toString());
			al.add(new Double(fractionAlignmentsPassQCScoreFilters).toString());
			al.add(new Double(fractionOnTargetAndPassQCScoreFilters).toString());
			al.add(new Double(estimatedFractionDuplicateAlignments).toString());
		}
		if (mpaParsed){
			al.add(new Double(meanInsertSize).toString());
			al.add(new Double(fractionOverlappingBpsInPairedReads).toString());
			al.add(new Long(numberPassingBps).toString());
			al.add(new Double(fractionPassingQ20bps).toString());
			al.add(new Double(fractionPassingQ30bps).toString());
		}
		if (s2uParsed){
			al.add(new Double(meanOnTargetCoverage).toString());
			al.add(new Double(coverageAt095OfTargetBps).toString());
			al.add(new Long(numberLowCoverageBps).toString());
		}
		return Misc.stringArrayListToString(al, "\t");
	}
	
	public String fetchTabbedHeader(){
		ArrayList<String> al = new ArrayList<String>();
		al.add("Sample Name");
		if (fastqParsed) al.add("# Fastq Reads");
		if (saeParsed){
			al.add("# Unfiltered Alignments");
			al.add("Fraction Passing QC");
			al.add("Fraction On Target");
			al.add("Fraction Duplicate");
		}
		if (mpaParsed){
			al.add("Mean Insert Size");
			al.add("Fraction Overlapping BPs");
			al.add("# Unique BPs");
			al.add("Fraction Q20 BPs");
			al.add("Fraction Q30 BPs");
		}
		if (s2uParsed){
			al.add("Mean on Target Coverage");
			al.add("Coverage at 0.95 of Target BPs");
			al.add("# Low Coverage Target BPs");
		}
		return Misc.stringArrayListToString(al, "\t");
	}
	
	public void appendHtmlColumns(StringBuilder sb) {
		//sb.append("	data.addColumn('string', '');\n");
		//sb.append("	data.addColumn('number', '');\n");
		//sb.append("	data.addColumn('boolean', '');\n");
		
		sb.append("	data.addColumn('string', 'Sample Name');\n");
		if (fastqParsed) sb.append("	data.addColumn('number', '# Fastq Reads');\n"); 
		if (saeParsed){
			sb.append("	data.addColumn('number', '# Unfiltered Alignments');\n");
			sb.append("	data.addColumn('number', 'Fraction Passing QC');\n");
			sb.append("	data.addColumn('number', 'Fraction On Target');\n");
			sb.append("	data.addColumn('number', 'Fraction Duplicate');\n");
		}
		if (mpaParsed){
			sb.append("	data.addColumn('number', 'Mean Insert Size');\n");
			sb.append("	data.addColumn('number', 'Fraction Overlapping BPs');\n");
			sb.append("	data.addColumn('number', '# Unique BPs');\n");
			sb.append("	data.addColumn('number', 'Fraction Q20 BPs');\n");
			sb.append("	data.addColumn('number', 'Fraction Q30 BPs');\n");
		}
		if (s2uParsed){
			sb.append("	data.addColumn('number', 'Mean on Target Coverage');\n");
			sb.append("	data.addColumn('number', 'Coverage at 0.95 of Target BPs');\n");
			sb.append("	data.addColumn('number', '# Low Coverage Target BPs');\n");
		}
	}
	
	public void appendHtmlDataRow(StringBuilder sb, boolean skipComma) {
		//sb.append("		['Mike',  100.3, true],\n");
		//sb.append("		['Jim',   {v:8000,   f: '$8,000'},  false],\n");
		//sb.append("		['Alice', {v: 12500, f: '$12,500'}, true],\n");
		//sb.append("		['Bob',   {v: 7000,  f: '$7,000'},  true]\n");
		ArrayList<String> al = new ArrayList<String>();
		int numDec = 3;
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(numDec);
		
		al.add("'"+sampleName+"'");
		if (fastqParsed) al.add(new Long(numberFastqReads).toString());
		if (saeParsed){
			al.add(new Long(numberUnfilteredAlignments).toString());
			al.add(f.format(fractionAlignmentsPassQCScoreFilters) );
			al.add(f.format(fractionOnTargetAndPassQCScoreFilters) );
			al.add(f.format(estimatedFractionDuplicateAlignments) );
		}
		if (mpaParsed){
			al.add(f.format(meanInsertSize) );
			al.add(f.format(fractionOverlappingBpsInPairedReads) );
			al.add(new Long(numberPassingBps).toString());
			al.add(f.format(fractionPassingQ20bps) );
			al.add(f.format(fractionPassingQ30bps) );
		}
		if (s2uParsed){
			al.add(f.format(meanOnTargetCoverage) );
			al.add(f.format(coverageAt095OfTargetBps) );
			al.add(new Long(numberLowCoverageBps).toString());
		}
		sb.append("\t\t[");
		sb.append(Misc.stringArrayListToString(al, ","));
		if (skipComma) sb.append("]\n");
		else sb.append("],\n");
	}
	
	
	public String fetchThresholds(String d, String r){
		ArrayList<String> al = new ArrayList<String>();
		if (saeParsed){
			al.add(mappingQualityThreshold+ d+ "MQ Threshold - Minimum mapping quality threshold, typically 13 or 20.");
			al.add(alignmentScoreThreshold+ d+ "AS Threshold - Alignment score threshold.");
			al.add(divideAlignmentScoreByCigarM+ d+ "AS/ CIGAR M - Was the given AS score divided by the CIGAR M length? This is needed for bwa but not novoalign to normalize for the length of the alignment.");
		}
		if (s2uParsed){
			al.add(minimumCoverageThreshold+ d+ "Coverage Threshold - Minimum unique observationcoverage threshold for counting as a failed BP.");
		}
		return Misc.stringArrayListToString(al, r);
	}
	
	/**begin, divide, end*/
	public String fetchDescriptions(String b, String d, String e){
		ArrayList<String> al = new ArrayList<String>();
		al.add(b+ "Sample Name"+ d+ "Name parsed from the json.gz files.");
		if (fastqParsed) al.add(b+ "# Fastq Reads"+ d+ "Total number of fastq reads passed to the aligner, not pairs of reads, all.");
		if (saeParsed){
			al.add(b+ "# Unfiltered Alignments"+ d+ "Number of unfiltered alignments.");
			al.add(b+ "Fraction Passing QC"+ d+ "Fraction of primary alignments passing vendor QC, MQ, and AS scores. This is a good measure of the fraction of quality alignments coming from a sample.");
			al.add(b+ "Fraction On Target"+ d+ "Fraction of these quality alignments that are also on target.");
			al.add(b+ "Fraction Duplicate"+ d+ "Fraction of these on target, quality alignments, that are marked as duplicate.");
		}
		if (mpaParsed){
			al.add(b+ "Mean Insert Size"+ d+ "Mean insert size for paired reads.");
			al.add(b+ "Fraction Overlapping BPs"+ d+ "Fraction BPs in paired alignments that overlap.");
			al.add(b+ "# Unique BPs"+ d+ "Total number of unique, duplicate free bps. This is a good estimation of the usable yield of data for a given samples.");
			al.add(b+ "Fraction Q20 BPs"+ d+ "Fraction of these unique BPs with good quality.");
			al.add(b+ "Fraction Q30 BPs"+ d+ "Fraction of these unique BPs with very good quality.");
		}
		if (s2uParsed){
			al.add(b+ "Mean on Target Coverage"+ d+ "Traditional measure of coverage over target BPs.");
			al.add(b+ "Coverage at 0.95 of Target BPs"+ d+ "Better measure of coverage, calculated by asking what fraction of target BPs have 0x, 1x, 2x or more coverage. Stop when it hits 0.95.");
			al.add(b+ "# Low Coverage Target BPs"+ d+ "Number of target BPs with less than the minimum coverage threshold.");
		}
		return Misc.stringArrayListToString(al, e);
	}
	
	public String[] parseStringArray(JsonValue jv){
		JsonArray ja = jv.asArray();
		int size = ja.size();
		String[] values = new String[size];
		for (int i=0; i< size; i++) values[i] = ja.get(i).asString();
		return values;
	}
	
	public float[] parseFloatArray(JsonValue jv){
		JsonArray ja = jv.asArray();
		int size = ja.size();
		float[] values = new float[size];
		for (int i=0; i< size; i++) values[i] = ja.get(i).asFloat();
		return values;
	}
	
	/**At what index does the fraction hit or fall below 0.25?*/
	public int whenHit25ReadCoverage() {
		for (int i=0; i< fractionTargetBpsWithIndexedCoverage.length; i++){
			if (fractionTargetBpsWithIndexedCoverage[i] <= 0.25) return i;
		}
		return fractionTargetBpsWithIndexedCoverage.length-1;
	}

	public boolean checkThresholds(SampleQC s) {
		if (s.getMappingQualityThreshold() != this.mappingQualityThreshold || s.getAlignmentScoreThreshold() != this.alignmentScoreThreshold || 
				s.divideAlignmentScoreByCigarM != this.divideAlignmentScoreByCigarM) return false;
		return true;
	}
	
	public boolean checkJsonFiles(SampleQC s) {
		if (s.isFastqParsed() != this.fastqParsed || s.isSaeParsed() != this.saeParsed || s.isMpaParsed() != this.mpaParsed || 
				s.isS2uParsed() != this.s2uParsed) return false;
		return true;
	}
	
	public String getSampleName() {
		return sampleName;
	}
	public boolean isSaeParsed() {
		return saeParsed;
	}
	public long getNumberUnfilteredAlignments() {
		return numberUnfilteredAlignments;
	}
	public double getFractionOnTargetAndPassQCScoreFilters() {
		return fractionOnTargetAndPassQCScoreFilters;
	}
	public double getEstimatedFractionDuplicateAlignments() {
		return estimatedFractionDuplicateAlignments;
	}
	public boolean isMpaParsed() {
		return mpaParsed;
	}
	public double getFractionOverlappingBpsInPairedReads() {
		return fractionOverlappingBpsInPairedReads;
	}
	public double getMeanInsertSize() {
		return meanInsertSize;
	}
	public long getNumberPassingBps() {
		return numberPassingBps;
	}
	public double getFractionPassingQ20bps() {
		return fractionPassingQ20bps;
	}
	public double getFractionPassingQ30bps() {
		return fractionPassingQ30bps;
	}
	public boolean isS2uParsed() {
		return s2uParsed;
	}
	public double getMeanOnTargetCoverage() {
		return meanOnTargetCoverage;
	}
	public int getMinimumCoverageThreshold() {
		return minimumCoverageThreshold;
	}
	public double getCoverageAt095OfTargetBps() {
		return coverageAt095OfTargetBps;
	}
	public float[] getFractionTargetBpsWithIndexedCoverage() {
		return fractionTargetBpsWithIndexedCoverage;
	}
	public String[] getLowCoverageRegions() {
		return lowCoverageRegions;
	}
	public String[] getExonicMedianPerBpCoverage() {
		return exonicMedianPerBpCoverage;
	}
	public double getFractionAlignmentsPassQCScoreFilters() {
		return fractionAlignmentsPassQCScoreFilters;
	}
	public double getMappingQualityThreshold() {
		return mappingQualityThreshold;
	}
	public double getAlignmentScoreThreshold() {
		return alignmentScoreThreshold;
	}
	public boolean isDivideAlignmentScoreByCigarM() {
		return divideAlignmentScoreByCigarM;
	}

	public boolean isFastqParsed() {
		return fastqParsed;
	}












}