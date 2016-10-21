package edu.utah.seq.query.brokenMongo;

import com.mongodb.BasicDBObjectBuilder;
import com.mongodb.DBObject;

import util.gen.Misc;

/**Container for one vcf record. Interbase coordinates. 
 * Start and Stop define the region of effected bps for SNV and DEL and flanking bps for INS*/
public class MongoVcf {
	
	//fields
	private String source;
	private String chr;
	private int start;
	private int stop;
	private String record;
	
	public MongoVcf(String source, int[] startStop, String record){
		this.source = source;
		this.record = record;
		
		//CHROM POS ID REF ALT QUAL FILTER INFO SAMPLE1 2 3 ...
		String[] fields = Misc.TAB.split(record);
		chr = fields[0];
		start = startStop[0]; 
		stop = startStop[1];
	}
	
	public DBObject createDBObject(){
		BasicDBObjectBuilder d = BasicDBObjectBuilder.start();
		d.append("source", source);
		d.append("chr", chr);
		d.append("start", start);
		d.append("stop", stop);
		d.append("record", record);
		return d.get();
	}

	public String getSource() {
		return source;
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

	public String getRecord() {
		return record;
	}
}
