package edu.utah.seq.query;

import java.io.File;
import java.io.IOException;
import java.text.StringCharacterIterator;
import java.util.ArrayList;
import java.util.HashMap;

import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.Bed;
import util.gen.Misc;

public class TabixQuery {
	
	public static void main (String[] args){
		try {
			TabixQuery tq = new TabixQuery("chrX", new RegionScoreText(1234, 1239, 0.4f, "descript of this region"));
			ArrayList<String> data = new ArrayList<String>();
			String dirty = "20\t18142564\trs3828013\tA\tG\t.\tPASS\tDB;ECNT=1;HCNT=16;MAX_ED=.;MIN_ED=.;NLOD=39.13;TLOD=452.16";
			data.add(dirty);
			data.add("21\t28142564\trs3828013\tA\tG\t.\tPASS\tDB;ECNT=1;HCNT=16;MAX_ED=.;MIN_ED=.;NLOD=39.13;TLOD=452.16");
			tq.getSourceResults().put(new File("/Users/u0028003/Desktop/IE/Data/chr20_1_3ConMut2.raw.vcf.gz"), null);
			
			ArrayList<String> data2 = new ArrayList<String>();
			data2.add("21\t28142564\trs3828013\tA\tG\t.\tPASS\tDB;ECNT=1;HCNT=16;MAX_ED=.;MIN_ED=.;NLOD=39.13;TLOD=452.16");
			tq.getSourceResults().put(new File("/Users/u0028003/Desktop/IE/Data/chr20_1_3ConMut2.pass.vcf.gz"), data);
			
			System.err.println(tq.toJson(true, "/Users/u0028003/Desktop/"));
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	//fields
	private String chr;
	private int start; //interbase	
	private int stop; //interbase
	private float score;
	private String info = null;
	
	//just for vcf, need pos since the start and stop might be padded
	private String pos = null;
	private String ref = null;
	private String[] alts = null;
	
	private HashMap<File, ArrayList<String>> sourceResults = new HashMap<File, ArrayList<String>>();
	
	//constructors
	/**For bed file region data, interbase coordinates assumed!
	 * @throws IOException */
	public TabixQuery(String chr, RegionScoreText r) throws IOException {
		this.chr = chr;
		start = r.getStart();
		stop = r.getStop();
		score = r.getScore();
		info = r.getText();
		//check coor
		if (stop <= start) throw new IOException("\nERROR: the stop is not > the start for TabixQuery "+getInterbaseCoordinates());
	}

	public TabixQuery(Bed bed) {
		chr = bed.getChromosome();
		start = bed.getStart();
		stop = bed.getStop();
		score = (float) bed.getScore();
		info = bed.getName();
	}
	
	/**Good to call if working with vcf data so that an exact match can be made*/
	public void parseVcf(){
		String[] t = Misc.TAB.split(info);
		pos = t[1];
		ref = t[3];
		alts = Misc.COMMA.split(t[4]);
	}
	
	/**Comparses these vcf args against the TQ for exact matching, note only one of the alts must match, not all.*/
	public boolean compareVcf(String otherPos, String otherRef, String[] otherAlts){
		if (pos.equals(otherPos) == false) return false;
		if (ref.equals(otherRef) == false) return false;
		//one of the alts needs to match
		for (String thisAlt: alts){
			for (String thatAlt: otherAlts){
				if (thisAlt.equals(thatAlt)) return true;
			}
		}
		return false;
	}

	/*Need to synchronize this since multiple threads could be adding results simultaneously.*/
	public synchronized void addResults(File source, ArrayList<String> results){
		ArrayList<String> al = sourceResults.get(source);
		if (al == null) {
			al = new ArrayList<String>();
			sourceResults.put(source, al);
		}
		al.addAll(results);
	}

	public String getTabixCoordinates() {
		return chr+":"+(start+1)+"-"+stop;
	}
	
	public String getInterbaseCoordinates() {
		return chr+":"+start+"-"+stop;
	}
	
	public String toJson(boolean includeScore, String pathToTrimmedFile){
		StringBuilder sb = new StringBuilder();
		sb.append("{\"chr\": \""); sb.append(chr); 
		sb.append("\", \"start\": "); sb.append(start); 
		sb.append(", \"stop\": "); sb.append(stop);
		if (info.length()!=0) {
			sb.append(", \"info\": \""); sb.append(info); sb.append("\"");
		}
		if (includeScore) {
			sb.append(", \"score\": "); sb.append(score);
		}
		if (sourceResults.size()!=0) {
			sb.append(",\n\"hits\": [\n");
			appendHits(sb, pathToTrimmedFile);
			sb.append("\t]");
		}
		sb.append("}");
		return sb.toString();
	}
	
	public void appendHits(StringBuilder sb, String pathToTrimmedFile){
		int numHits = sourceResults.size();
		int hitCounter = 0;
		for (File source: sourceResults.keySet()){
			String trimmedName = source.toString().replaceFirst(pathToTrimmedFile, "");
			//any results? 
			ArrayList<String> results = sourceResults.get(source);
			if (results == null || results.size() ==0){
				sb.append("\t{\"source\": \""); sb.append(trimmedName); sb.append("\"}");
				if (++hitCounter != numHits) sb.append(",\n");
				else sb.append("\n");
			}
			else {
				sb.append("\t{\"source\": \""); sb.append(trimmedName); sb.append("\",\n\t \"data\": [");

				int numData = sourceResults.get(source).size();
				int dataCounter = 0;
				for (String data: sourceResults.get(source)){
					//replace tabs with \t, " with '
					String clean = data.replace("\t", "\\t");
					clean = clean.replace("\"", "'");
					sb.append("\n\t\t\""); sb.append(clean); sb.append("\"");
					if (++dataCounter != numData) sb.append(", ");
				}
				if (++hitCounter != numHits) sb.append("]},\n");
				else sb.append("]}\n");
			}
		}
	}
	
	public static String getInterbaseCoordinates(ArrayList<TabixQuery> al){
		StringBuilder sb = new StringBuilder();
		sb.append(al.get(0).getInterbaseCoordinates());
		for (int i=1; i< al.size(); i++){
			sb.append(",");
			sb.append(al.get(i).getInterbaseCoordinates());
		}
		return sb.toString();
	}

	public HashMap<File, ArrayList<String>> getSourceResults() {
		return sourceResults;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public String getInfo() {
		return info;
	}
	
	
	
}
