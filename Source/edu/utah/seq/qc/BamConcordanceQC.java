package edu.utah.seq.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import com.eclipsesource.json.JsonObject;
import util.gen.IO;
import util.gen.Misc;

public class BamConcordanceQC {
	
	private String[] sampleNames = null;
	private String[] bamFileNames = null;
	private String[] similarities = null;
	private String[] genderChecks = null;
	private String similarity;
	private String genderCheck;
	
	public BamConcordanceQC(File bcj, HashMap<String, BamConcordanceQC> bamFileNameBCR) throws IOException {
		
		//load json objects
		BufferedReader in = IO.fetchBufferedReader(bcj);
		JsonObject jo = JsonObject.readFrom(in);
		sampleNames = SampleQC.parseStringArray(jo.get("sampleNames"));
		bamFileNames = SampleQC.parseStringArray(jo.get("bamFileNames"));
		similarities = SampleQC.parseStringArray(jo.get("similarities"));
		genderChecks = SampleQC.parseStringArray(jo.get("genderChecks"));
		in.close();
		
		//check that all were found
		if (sampleNames == null || bamFileNames == null || similarities == null || genderChecks == null) throw new IOException ("Failed to parse the four required fields from "+bcj);

		//load hash
		for (String bfn: bamFileNames){
			if (bamFileNameBCR.containsKey(bfn)) throw new IOException("This file name with associated bam concordance has already been seen. Are your bam file names unique? Look for "+bfn);
			bamFileNameBCR.put(bfn, this);
		}

		//make similarity and gender strings
		similarity = Misc.stringArrayToString(similarities, "; ");
		genderCheck = Misc.stringArrayToString(genderChecks, "; ");
	}
	public String getSimilarity() {
		return similarity;
	}
	public String getGenderCheck() {
		return genderCheck;
	}

}
