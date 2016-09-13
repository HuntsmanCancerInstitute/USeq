package edu.utah.seq.query;

import java.util.ArrayList;
import java.util.HashMap;
import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import util.gen.Misc;

public class MongoLoader implements Runnable{

	//fields
	private ArrayList<MongoQuery> toQuery;
	private DBCollection collection; 


	public MongoLoader (DBCollection collection, ArrayList<MongoQuery> toQuery) {
		this.toQuery = toQuery;
		this.collection = collection;
	}

	/**This walks through each query, fetching the sources to query, executes each and then loads the results back into the query.
	 * Assumes no other threads will be working with these queries until all the mongo executions are complete!*/
	public void run() {
		int size = toQuery.size();
		
		//for each query
		for (int i=0; i< size; i++){
			MongoQuery mq = toQuery.get(i);
			//container for source: results
			HashMap<String, ArrayList<DBObject>> sourceResults = mq.getSourceResults();

			//for each source
			for (String source : sourceResults.keySet()){
				BasicDBObject q = new BasicDBObject();
				//exact match
				q.put("source", source);
				q.put("chr", mq.getChr());
				//range query, might want to tune this by bracketing the start and stop?
				q.put("start", new BasicDBObject("$lt", mq.getStop()));
				q.put("stop", new BasicDBObject("$gt", mq.getStart()));

				//there should be results every time
				ArrayList<DBObject> results = new ArrayList<DBObject>();
				DBCursor cursor = collection.find(q);
				while(cursor.hasNext())  results.add(cursor.next());
				if (results.size() == 0) Misc.printErrAndExit("\nError: query returned no results for "+mq.getInterbaseCoordinates()+" "+source+"\n"+q);
				sourceResults.put(source, results);
			}

			//save it, not thread safe so don't modify this object
			mq.setSourceResults(sourceResults);
		}
	}

	public ArrayList<MongoQuery> getToQuery() {
		return toQuery;
	}
}
