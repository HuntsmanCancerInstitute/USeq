package edu.utah.seq.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeSet;
import org.json.JSONArray;
import org.json.JSONObject;
import util.gen.IO;
import util.gen.Misc;

public class SampleQC2 {

	//fields
	private String sampleName;
	private DnaSample normalDnaSample = null;
	private DnaSample tumorDnaSample = null;
	
	//Sample Concordance
	private SampleConcordanceQC2 sampleConcordance = null;
	
	//ClinInfo
	private String diagnosis = "NA";
	private String analysisId = "NA";
	private String gender = "NA";

	//constructor
	public SampleQC2(String name, File normAlignLog, File normDupLog, File normReadCovJson, File tumAlignLog, File tumDupLog, File tumReadCovJson, File sampConcJson, File clinInfoJson, File clinInfoTxt) throws Exception {
		sampleName = name;
		//load DNA stats
		if (normAlignLog != null || normDupLog != null || normReadCovJson != null) normalDnaSample = new DnaSample(normAlignLog, normDupLog, normReadCovJson);
		if (tumAlignLog != null || tumDupLog != null || tumReadCovJson != null) tumorDnaSample = new DnaSample(tumAlignLog, tumDupLog, tumReadCovJson);
		//load sample concordance
		if (sampConcJson != null) sampleConcordance = new SampleConcordanceQC2(sampConcJson);
		//parse clinical info
		parseClinical(clinInfoJson, clinInfoTxt);
	}

	/**Loads all of the sample diagnosis.*/
	public void parseClinical (File json, File txt) throws Exception {

		// parse from Avatar?
		if (json != null) {
			JSONObject jo = new JSONObject(IO.loadFile(json, " ", true));
			if (jo.has("AnalysisId")) analysisId = jo.getString("AnalysisId");
			else analysisId = "NA";
			if (jo.has("Gender")) gender = jo.getString("Gender");
			else gender = "NA";

			JSONArray ja = jo.getJSONArray("Samples");
			int numSamples = ja.length();
			TreeSet<String> d = new TreeSet<String>();
			for (int i=0; i< numSamples; i++){
				JSONObject sample = ja.getJSONObject(i);
				if (sample.has("Diagnosis")) d.add(sample.getString("Diagnosis"));
			}
			if (d.size() == 0) diagnosis = "NA";
			else diagnosis = Misc.treeSetToString(d, ",");
		}
		else if (txt != null) setGender(txt);
	}

	private void setGender(File txt) throws IOException {	
		BufferedReader in = IO.fetchBufferedReader(txt);
		String line = null;
		while ((line = in.readLine()) !=null) {
			line = line.toLowerCase();
			if (line.contains("gender")) {
				if (line.contains("female")) {
					gender = "FEMALE";
					in.close();
					return;
				}
				else if (line.contains("male")) {
					gender = "MALE";
					in.close();
					return;
				}
			}
		}
		in.close();
	}

	public String fetchTabbedLine() {
		
			ArrayList<String> preAl = new ArrayList<String>();
			//pre
			preAl.add(sampleName);
			preAl.add(diagnosis);
			preAl.add(analysisId);
			preAl.add(gender);
			String pre = Misc.stringArrayListToString(preAl, "\t");
			
			//post
			String post = null;
			if (sampleConcordance!=null) post = sampleConcordance.getTabbedLine();
			else post = SampleConcordanceQC2.getDefaultTabbedLine();
			
			
			//look for tumor
			ArrayList<String> tumorLineAl = null;
			if (tumorDnaSample != null) {
				tumorLineAl = new ArrayList<String>();
				tumorLineAl.add(pre);
				tumorLineAl.add(tumorDnaSample.getTabbedLine("tumor"));
				tumorLineAl.add(post);
			}
			
			//look for normal
			ArrayList<String> normalLineAl = null;
			if (normalDnaSample != null) {
				normalLineAl = new ArrayList<String>();
				normalLineAl.add(pre);
				normalLineAl.add(normalDnaSample.getTabbedLine("normal"));
				normalLineAl.add(post);
			}
			
			ArrayList<String> finalLines = new ArrayList<String>();
			if (tumorLineAl != null) {
				finalLines.addAll(tumorLineAl);
				finalLines.add("\n");
			}
			if (normalLineAl != null) {
				finalLines.addAll(normalLineAl);
				finalLines.add("\n");
			}
			
			return Misc.stringArrayListToString(finalLines, "");
	}
	
	public static String fetchTabbedHeader(){
		ArrayList<String> al = new ArrayList<String>();
		al.add("Sample Name"); //Job Dir Name
		al.add("Diagnosis"); //ALL or Ovarian
		al.add("Analysis Id"); //A1234
		al.add("Reported Gender"); // Male or Female or NA
		al.add("Sample Type"); // Tumor or Normal

		al.add(DnaSample.getHeader());
		al.add(SampleConcordanceQC2.getHeader());

		return Misc.stringArrayListToString(al, "\t");
	}
	

	
	/**begin, divide, end*/
	public static String fetchDescriptions(String b, String d, String e){
		ArrayList<String> al = new ArrayList<String>();
		al.add(b+ "Sample Name"+ d+ "Name parsed from the root job directory");
		al.add(b+ "Diagnosis"+ d+ "Submitted diagnosis, typically tumor type.");
		al.add(b+ "Analysis Id"+ d+ "GNomEx Analysis ID.");
		al.add(b+ "Reported Gender"+ d+ "Gender reported in the job folder.");
		al.add(b+ "# Fastq Reads"+ d+ "Total number of fastq reads passed to the aligner, not pairs of reads, all.");
		al.add(b+ "# Unfiltered Alignments"+ d+ "Number of unfiltered alignments.");
		al.add(b+ "Fraction Duplicate"+ d+ "Fraction of quality alignments that are marked as duplicate.");
		al.add(b+ "Fraction Passing QC"+ d+ "Fraction of alignments passing all filters and written to file.");
		al.add(b+ "# Filtered Alignments"+ d+ "Number of alignments passing all filters and written to file.");
		al.add(b+ "Mean on Target Coverage"+ d+ "Traditional measure of coverage over target BPs.");
		al.add(b+ "Coverage at 0.9 of Target BPs"+ d+ "Better measure of coverage, calculated by asking what fraction of target BPs have 0x, 1x, 2x or more coverage. Stop when it hits 0.9.");
		al.add(b+ "Coverage at 0.95 of Target BPs"+ d+ "Better measure of coverage, calculated by asking what fraction of target BPs have 0x, 1x, 2x or more coverage. Stop when it hits 0.95.");
		al.add(b+ "Similarities Fwd Rev"+ d+ "SampleConcordance similarity estimates for each sample set in the forward and reverse direction.");
		al.add(b+ "Concordance Check"+ d+ "Do the samples pass homozygous variant concordance");
		al.add(b+ "Het/Hom All ChrX Log2(All/ChrX)"+ d+ "SampleConcordance gender check for each sample, ratio of Het/Hom for All, chrX, and log2(All/ChrX).");
		al.add(b+ "Gender Call"+ d+ "Estimated gender.");
		al.add(b+ "Gender Check"+ d+ "Does the gender estimate match the reported gender.");
		return Misc.stringArrayListToString(al, e);
	}
	
	public String getSampleName() {
		return sampleName;
	}

	public DnaSample getNormalDnaSample() {
		return normalDnaSample;
	}

	public DnaSample getTumorDnaSample() {
		return tumorDnaSample;
	}


}
