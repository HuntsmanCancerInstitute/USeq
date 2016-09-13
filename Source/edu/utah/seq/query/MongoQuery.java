package edu.utah.seq.query;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import com.mongodb.DBObject;
import edu.utah.seq.useq.data.RegionScoreText;

public class MongoQuery {
	
	//fields
	private String chr;
	private int start; //interbase	
	private int stop; //interbase
	private HashMap<String, ArrayList<DBObject>> sourceResults;
	
	//constructors
	/**For bed file region data, interbase coordinates assumed!
	 * @throws IOException */
	public MongoQuery(String chr, RegionScoreText r) throws IOException {
		start = r.getStart();
		stop = r.getStop();
		this.chr = chr;
		//check coor
		if (stop <= start) throw new IOException("\nERROR: the stop is not > the start for MongoQuery "+getInterbaseCoordinates());
	}

	public String getOneBasedCoordinates() {
		return chr+":"+(start+1)+"-"+stop;
	}
	
	public String getInterbaseCoordinates() {
		return chr+":"+start+"-"+stop;
	}
	
	public static String getInterbaseCoordinates(ArrayList<MongoQuery> al){
		StringBuilder sb = new StringBuilder();
		sb.append(al.get(0).getInterbaseCoordinates());
		for (int i=1; i< al.size(); i++){
			sb.append(",");
			sb.append(al.get(i).getInterbaseCoordinates());
		}
		return sb.toString();
	}

	public HashMap<String, ArrayList<DBObject>> getSourceResults() {
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

	public void setSourceResults(HashMap<String, ArrayList<DBObject>> sourceResults) {
		this.sourceResults = sourceResults;
	}
	
	
	
}
