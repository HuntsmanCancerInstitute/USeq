package util.gen;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.json.*;



public class Delme {

	public static void main(String[] args) throws IOException {
		/*HashMap<String, String> nameStats = loadStats(new File("/Users/u0028003/Downloads/OvarBamConParsing/stats.txt"));
		
		File[] jsonFiles = IO.extractFiles(new File("/Users/u0028003/Downloads/OvarBamConParsing/Jsons"), ".json");
		 HashSet<String> noPrint = addJsons(jsonFiles, nameStats);
		 */

		test();
		 
	}
	
	private static void test() {
		String[] out = {"a","b","c"};
		for (String s: out) {
			IO.pl("Looking at "+s+" "+s.length());
			if (s.startsWith("b")) {
				if (s.length() == 1) {
					IO.pl("\tbreaking");
					break;
				}
				return;
			}
		}
		IO.pl("break worked");
	}

	private static HashMap<String, String> loadStats(File file) throws IOException {
		HashMap<String, String> nameStats = new HashMap<String,String>();
		BufferedReader in = IO.fetchBufferedReader(file);
		String[] fields;
		String line;
		while ((line = in.readLine())!= null) {
			line = line.trim();
			if (line.length()==0) continue;
			fields = Misc.TAB.split(line);
			nameStats.put(fields[0], line);
		}
		in.close();
		return nameStats;
	}

	private static HashSet<String> addJsons(File[] jsonFiles, HashMap<String, String> oriStats) throws IOException {
		
		
		PrintWriter out = new PrintWriter (new FileWriter( new File("/Users/u0028003/Downloads/OvarBamConParsing/merged.xls")));
		HashSet<String> noPrint = new HashSet<String>();
		
		for (File jf: jsonFiles) {
			String patientId = jf.getName().replace("_BamConcordance.json", "");
			IO.pl(patientId);
			BufferedReader is = IO.fetchBufferedReader(jf);
	        
	        JSONTokener tokener = new JSONTokener(is);
	        JSONObject object = new JSONObject(tokener);
	        
	        String concordanceStats = object.getJSONArray("similarities").getString(0);
	        String concordanceComparisonResult = object.getString("concordanceCheck");
	        String genderStats = object.getJSONArray("genderChecks").getString(0)+","+object.getJSONArray("genderChecks").getString(1);
	        String gender = object.getString("genderCall");
	        String genderComparisonResult = object.getString("genderComparison");
	        String stats = patientId+"\t"+ concordanceStats+"\t"+ concordanceComparisonResult+"\t"+ genderStats+"\t"+ gender+"\t"+ genderComparisonResult;

	        //for each sample name
	        JSONArray courses = object.getJSONArray("sampleNames");
	        for (int i = 0; i < courses.length(); i++) {
	        	String name = courses.getString(i);
	        	name = name.replace("_NormalDNA_Hg38_final.cram", "");
	        	name = name.replace("_Hg38_final.cram", "");
	            out.println(oriStats.get(name) +"\t"+name+"\t"+ stats );
	            noPrint.add(name);
	        }
	        is.close();
		}
		
		Iterator<String> it = oriStats.keySet().iterator();
		while (it.hasNext()) {
			String name = it.next();
			if (noPrint.contains(name)) continue;
			out.println(oriStats.get(name));
		}
		
		out.close();
		return noPrint;
	}
	

}
