package edu.utah.seq.query.brokenMongo;

import com.mongodb.BasicDBObjectBuilder;
import com.mongodb.DBObject;

import util.gen.Misc;

/**Simple container for one or more bed records with the same source, chr, start */
public class MongoBed {
	
	//fields
	private String _id; 
	private String[] records;
	
	public MongoBed(String source, String[] records){
		//chr start stop ...
		this.records = records;
		String[] fields = Misc.TAB.split(records[0]);
		setId(source, fields);
	}
	
	public DBObject createDBObject(){
		BasicDBObjectBuilder d = BasicDBObjectBuilder.start();
		d.append("_id", _id);
		for (String bed: records) d.append("bed", bed);
		return d.get();
	}
	
	/** Source\tchr\tstart */
	public void setId(String source, String[] bedFields){
		StringBuilder sb = new StringBuilder (source);
		sb.append("\t");
		sb.append(bedFields[0]);
		sb.append("\t");
		sb.append(Integer.parseInt(bedFields[1]));
		_id = sb.toString();
	}

	public String get_id() {
		return _id;
	}

	public String[] getRecords() {
		return records;
	}
}
