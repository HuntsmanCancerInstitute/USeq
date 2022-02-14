package edu.utah.seq.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import org.json.JSONObject;

import com.eclipsesource.json.JsonObject;
import util.gen.IO;
import util.gen.Misc;

public class SampleConcordanceQC2 {
	
	private String[] sampleNames = null;
	private String[] bamFileNames = null;
	private String[] similarityStats = null;
	private String[] genderStats = null;
	private String concordanceCheck = null; //PASS FAIL
	private String genderCheck = null; //PASS FAIL
	private String genderCall = null;  //MALE FEMALE
	
	public SampleConcordanceQC2(File json) throws IOException {
		//load json objects
		BufferedReader in = IO.fetchBufferedReader(json);
		JsonObject jo = JsonObject.readFrom(in);
		sampleNames = SampleQC.parseStringArray(jo.get("sampleNames"));
		bamFileNames = SampleQC.parseStringArray(jo.get("bamFileNames"));
		similarityStats = SampleQC.parseStringArray(jo.get("similarities"));
		genderStats = SampleQC.parseStringArray(jo.get("genderChecks"));
		concordanceCheck = jo.get("concordanceCheck").asString();
		genderCall = jo.get("genderCall").asString();
		genderCheck = jo.get("genderComparison").asString();
		in.close();
	}
	public static String getHeader() {
		return "Similarities Fwd Rev\tConcordance Check\tHet/Hom All ChrX Log2(All/ChrX)\tGender Call\tGender Check";
	}
	public String getTabbedLine() {
		StringBuilder sb = new StringBuilder(Misc.stringArrayToString(similarityStats, ", "));
		sb.append("\t");
		sb.append(concordanceCheck);
		sb.append("\t");
		sb.append(Misc.stringArrayToString(genderStats, ", "));
		sb.append("\t");
		sb.append(genderCall);
		sb.append("\t");
		sb.append(genderCheck);
		return sb.toString();
	}
	public static String getDefaultTabbedLine() {
		return "NA\tNA\tNA\tNA\tNA";
	}

}
