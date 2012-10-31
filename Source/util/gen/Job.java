package util.gen;
import java.util.*;

public class Job {
	
	//fields
	private String jobId;					//Job ID
	private String userName;				//Username
	private String queue;					//Queue
	private String jobName;				//Jobname
	private String sessionId;				//SessID
	private String numberRequestedNodes;	//NDS
	private String numberRequestedTasks;	//TSK
	private String requestedMemory;		//Req'd Memory
	private String requestedCPUWallTime;	//Req'd Time
	private String state;					//S
	private String cpuWallTimeUsed;		//Elap Time
	
	//constructor
	public Job (String[] lines){
		//attempt to parse jobId
		jobId = parseJobId(lines); 
		//if sucessfull parse everything else
		if (jobId != null){
			HashMap hash = loadHash(lines);
		}
	}
	
	
	
	/**Converts a String[] of key = values into a HashMap. Skips lines that don't
	 * contain a ' = '.*/
	public static HashMap loadHash(String[] lines){
		HashMap hash = new HashMap(25);
		for (int i=0; i<lines.length; i++){
			lines[i]= lines[i].trim();
			int index = lines[i].indexOf(" = ");
			if (index != -1){
				String key = lines[i].substring(0,index-3);
				String value = lines[i].substring(index);
				if (hash.containsKey(key)) System.err.println("Error in Job.loadHash(): key already exists -> "+lines[i]);
				hash.put(key, value);
			}
		}
		return hash;
	}
	
	/**Attempts to parse 'Job Id:', returns null if not found.*/
	public static String parseJobId(String[] lines) {
		for (int i=0; i<lines.length; i++){
			int index = lines[i].indexOf("Job Id: ");
			if (index != -1){
				return lines[i].substring(index);
			}
		}
		return null;
	}
	
}
