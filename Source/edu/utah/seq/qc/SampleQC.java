package edu.utah.seq.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;

import org.json.JSONArray;
import org.json.JSONObject;

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
	private String targetRegionsFileNameSAE;
	
	//DUP
	private boolean dupParsed = false;
	//will also load the estimatedFractionDuplicateAlignments

	//MPA
	private boolean mpaParsed = false;
	private String mpaBamFileName;
	private double fractionOverlappingBpsInPairedReads;
	private double meanInsertSize;
	private long numberPassingBps;
	private double fractionPassingQ20bps;
	private double fractionPassingQ30bps;

	//S2U
	private boolean s2uParsed = false;
	private double meanOnTargetCoverage;
	private int minimumCoverageThreshold;
	private double coverageAt090OfTargetBps;
	private double coverageAt095OfTargetBps;
	private float[] fractionTargetBpsWithIndexedCoverage;
	private String[] lowCoverageRegions; 
	private String[] exonicMedianPerBpCoverage;
	private long numberLowCoverageBps;
	private String targetRegionsFileNameS2U;
	
	//AvatarInfo
	JSONObject avatarInfo = null;
	String diagnosis = null;
	boolean parseReadCoverageStats = true;
	boolean swapExomeForDNA = false;

	//constructor
	public SampleQC( String sampleName, boolean parseReadCoverageStats, boolean swapExomeForDNA){
		this.sampleName = sampleName;
		this.parseReadCoverageStats = parseReadCoverageStats;
		this.swapExomeForDNA = swapExomeForDNA;
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
				//only load if not set, better to get from the dup json
				if (dupParsed == false) estimatedFractionDuplicateAlignments = jo.getDouble("estimatedFractionDuplicateAlignments", -1);
				//MQ threshold,
				mappingQualityThreshold = jo.getDouble("mappingQualityThreshold", -1);
				//AS threshold see boolean it might be divided the the length of the M's in the CIGAR String
				alignmentScoreThreshold = jo.getDouble("alignmentScoreThreshold", -1);
				divideAlignmentScoreByCigarM = jo.getBoolean("divideAlignmentScoreByCigarM", true);
				//target file
				targetRegionsFileNameSAE = jo.getString("targetRegionsFileName", "notFound");
			}
			else if (type.equals("dup")){				
				dupParsed = true;
				//fraction duplicate reads, should overwrite anyting from sae
				estimatedFractionDuplicateAlignments = jo.getDouble("fractionDuplicateAlignments", -1.0);
				if (estimatedFractionDuplicateAlignments == -1.0) estimatedFractionDuplicateAlignments = jo.getDouble("estimatedFractionDuplicateAlignments", -1.0);
			}
			else if (type.equals("mpa")){
				mpaParsed = true;
				//bam file name that was parsed, fix issue with prior naming?
				mpaBamFileName = jo.getString("bamFileName", "notFound");
				if (swapExomeForDNA) mpaBamFileName = mpaBamFileName.replace("Exome", "DNA");
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
				if (parseReadCoverageStats) {
					//this is the fraction of target bps with >=0x, >=1x, >=2x .....
					fractionTargetBpsWithIndexedCoverage = parseFloatArray(jo.get("fractionTargetBpsWithIndexedCoverage"));
					//this is the median of the per bp coverage across each region supplied to S2U
					exonicMedianPerBpCoverage = parseStringArray(jo.get("exonicMedianPerBpCoverage"));
					//bps that fail the minimum coverage threshold
					lowCoverageRegions = parseStringArray(jo.get("lowCoverageRegions"));
				}
				//the coverage # found at 0.90 of target bps, calculate by asking what fraction of target bps have at minimum, 0x 1x, 2x, 3x, ... stop when fraction hits 0.90.
				coverageAt090OfTargetBps = jo.getDouble("coverageAt0.90OfTargetBps", -1);
				//the coverage # found at 0.95 of target bps, calculate by asking what fraction of target bps have at minimum, 0x 1x, 2x, 3x, ... stop when fraction hits 0.95.
				coverageAt095OfTargetBps = jo.getDouble("coverageAt0.95OfTargetBps", -1);
				//the threshold, typically 20-30x for germline, 300x for somatic
				minimumCoverageThreshold = jo.getInt("minimumCoverageThreshold", -1);
				//number of target bps failing the coverage threshold
				//numberLowCoverageBps = jo.getLong("numberLowCoverageBps", -1);
				//name of targets
				targetRegionsFileNameS2U = jo.getString("targetRegionsFileName", "notFound");
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR reading json file "+f);
		}
	}
	
	/**Loads all of the sample diagnosis.*/
	public String getDiagnosis(){
		if (diagnosis != null && diagnosis.equals("NA") == false) return diagnosis;
		TreeSet<String> d = new TreeSet<String>();
		try {
			JSONArray ja = avatarInfo.getJSONArray("Samples");
			int numSamples = ja.length();
			for (int i=0; i< numSamples; i++){
				JSONObject jo = ja.getJSONObject(i);
				if (jo.has("Diagnosis")) d.add(jo.getString("Diagnosis"));
			}
		}
		catch (Exception e){}

		if (d.size() == 0) diagnosis = "NA";
		else diagnosis = Misc.treeSetToString(d, ",");
		return diagnosis;
	}

	public String fetchTabbedLine(BamConcordanceQC bc) {
		try {
			ArrayList<String> al = new ArrayList<String>();
			//add sample name
			al.add(sampleName);
			//add Analysis ID?
			if (avatarInfo != null){
				al.add(getDiagnosis());
				//Add analysis ID
				String analysisId = null;
				try {
					analysisId = avatarInfo.getString("AnalysisId");
				} catch (Exception e){}
				if (analysisId == null) analysisId = "NA";
				al.add(analysisId);
			}
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
				al.add(new Double(coverageAt090OfTargetBps).toString());
				al.add(new Double(coverageAt095OfTargetBps).toString());
				//al.add(new Long(numberLowCoverageBps).toString());
			}
			String genderCheck = null;
			if (bc != null){
				al.add(bc.getSimilarity());
				genderCheck = bc.getGenderCheck();
				al.add(genderCheck);
			}
			else {
				al.add("NA");
				al.add("NA");
			}
			//add gender?
			String gender = null;
			if (avatarInfo != null){
				try {
					gender = avatarInfo.getString("Gender");
				} catch (Exception e){}
				
				if (gender == null) al.add("NA");
				else if (genderCheck == null) al.add(gender);
				else {
					//ok both genderCheck and gender are present
					if (gender.equals("F")){
						if (genderCheck.contains(" MALE"))al.add(gender+" FAIL");
						else al.add(gender+" PASS");
					}
					else if (gender.equals("M")) {
						if (genderCheck.contains(" FEMALE"))al.add(gender+" FAIL");
						else al.add(gender+" PASS");
					}
					else al.add(gender);
				}
			}
			
			
			return Misc.stringArrayListToString(al, "\t");
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing sample "+sampleName);
		}
		return null;
	}
	
	public String fetchTabbedHeader(boolean includeAvatarInfo, boolean includeBC){
		ArrayList<String> al = new ArrayList<String>();
		al.add("Sample Name");
		if (includeAvatarInfo) {
			al.add("Diagnosis");
			al.add("AnalysisId");
		}
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
			al.add("Coverage at 0.9 of Target BPs");
			al.add("Coverage at 0.95 of Target BPs");
			//al.add("# Low Coverage Target BPs");
		}
		if (includeBC){
			al.add("Similarities Fwd Rev");
			al.add("Het/Hom All ChrX Log2(All/ChrX)");
		}
		if (includeAvatarInfo){
			al.add("Gender Check");
		}
		return Misc.stringArrayListToString(al, "\t");
	}
	
	public void appendHtmlColumns(StringBuilder sb, boolean addBC) {
		sb.append("	data.addColumn('string', 'Sample Name');\n");
		if (avatarInfo != null) {
			sb.append("	data.addColumn('string', 'Diagnosis');\n");
			sb.append("	data.addColumn('string', 'AnalysisId');\n");
		}
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
			sb.append("	data.addColumn('number', 'Coverage at 0.9 of Target BPs');\n");
			sb.append("	data.addColumn('number', 'Coverage at 0.95 of Target BPs');\n");
		}
		if (addBC){
			sb.append("	data.addColumn('string', 'Similarities Fwd Rev');\n");
			sb.append("	data.addColumn('string', 'Het/Hom All ChrX Log2(All/ChrX)');\n");
		}
		if (avatarInfo != null) sb.append("	data.addColumn('string', 'Gender Check');\n");
	}
	
	
	
	public void appendHtmlDataRow(StringBuilder sb, boolean skipComma, HashMap<String, BamConcordanceQC> bamFileNameBCR) throws IOException {
		//sb.append("		['Mike',  100.3, true],\n");
		//sb.append("		['Jim',   {v:8000,   f: '$8,000'},  false],\n");
		//sb.append("		['Alice', {v: 12500, f: '$12,500'}, true],\n");
		//sb.append("		['Bob',   {v: 7000,  f: '$7,000'},  true]\n");
		ArrayList<String> al = new ArrayList<String>();
		int numDec = 3;
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(numDec);
		f.setGroupingUsed(false); //no commas!
		
		al.add("'"+sampleName+"'");
		
		if (avatarInfo != null) {
			al.add("'"+getDiagnosis()+"'");
			//Add Analysis ID
			String aId = null;
			try {
				aId = avatarInfo.getString("AnalysisId");
			} catch (Exception e){}
			if (aId == null) aId = "NA";
			al.add("'"+aId+"'");
			
		}
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
			al.add(f.format(coverageAt090OfTargetBps) );
			al.add(f.format(coverageAt095OfTargetBps) );
			//al.add(new Long(numberLowCoverageBps).toString());
		}
		if (bamFileNameBCR != null){
			BamConcordanceQC bc = bamFileNameBCR.get(mpaBamFileName);
			if (bc == null) {
				al.add("'NA'");
				al.add("'NA'");
			}
			else {
				al.add("'"+bc.getSimilarity()+"'");
				al.add("'"+bc.getGenderCheck()+"'");
			}
		}
		if (avatarInfo != null) {
			String gender = null;
			try {
				gender = avatarInfo.getString("Gender");
			} catch (Exception e){}
			if (gender == null) gender = "NA";
			al.add("'"+gender+"'");
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
	public String fetchDescriptions(String b, String d, String e, boolean addBC){
		ArrayList<String> al = new ArrayList<String>();
		al.add(b+ "Sample Name"+ d+ "Name parsed from the json.gz files.");
		if (avatarInfo != null) {
			al.add(b+ "Diagnosis"+ d+ "Submitted diagnosis, typically tumor type.");
			al.add(b+ "AnalysisId"+ d+ "GNomEx Analysis ID.");
		}
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
			al.add(b+ "Coverage at 0.9 of Target BPs"+ d+ "Better measure of coverage, calculated by asking what fraction of target BPs have 0x, 1x, 2x or more coverage. Stop when it hits 0.9.");
			al.add(b+ "Coverage at 0.95 of Target BPs"+ d+ "Better measure of coverage, calculated by asking what fraction of target BPs have 0x, 1x, 2x or more coverage. Stop when it hits 0.95.");
			//al.add(b+ "# Low Coverage Target BPs"+ d+ "Number of target BPs with less than the minimum coverage threshold.");
		}
		if (addBC){
			al.add(b+ "Similarities Fwd Rev"+ d+ "BamConcordance similarity estimates for each sample set in the forward and reverse direction.");
			al.add(b+ "Het/Hom All ChrX Log2(All/ChrX)"+ d+ "BamConcordance gender check for each sample, ratio of Het/Hom for All, chrX, and log2(All/ChrX).");
		}
		if (avatarInfo != null) al.add(b+ "Gender"+ d+ "Reported Gender in GNomEx");
		return Misc.stringArrayListToString(al, e);
	}
	
	public static String[] parseStringArray(JsonValue jv){
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
		if (parseReadCoverageStats == false) return -1; 
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
	public double getCoverageAt090OfTargetBps() {
		return coverageAt090OfTargetBps;
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
	public String getTargetRegionsFileNameS2U() {
		return targetRegionsFileNameS2U;
	}
	public String getTargetRegionsFileNameSAE() {
		return targetRegionsFileNameSAE;
	}
	public String getMpaBamFileName() {
		return mpaBamFileName;
	}
	public void setAvatarInfo(JSONObject avatarInfo) {
		this.avatarInfo = avatarInfo;
	}














}