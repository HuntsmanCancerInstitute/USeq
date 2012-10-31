package util.bio.wrappers;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import util.gen.IO;
import util.gen.Passwords;

public class ClusterJob {
	
	//fields
	private String jobID;
	private File shellScript;
	private File results;
	private File gzippedResults;
	private File data1;
	private File data2;
	private File log;
	private String name;
	private CHPCAligner chpcAligner;
	private ArrayList<String> script;
	private boolean filterForChrLines = false;
	
	public void writeShellScript() {
		String rndName = "CHPCAlignerDeleteMe_"+Passwords.createRandowWord(7);

		script = new ArrayList<String>();
		//add header
		script.add(fetchPBSHeader());
		//name
		script.add("echo '"+name+"'");
		//start time
		script.add("echo -n 'Start:\t'; date +%s");
		//change umask so others can delete the files
		script.add("umask 000");
		//attempt a cleanup, this likely won't work due to permissions
		script.add("rm -rf /scratch/local/*");
		//make tmp dir
		script.add("mkdir /scratch/local/"+rndName+"/");
		//copy over genome index
		script.add("cp -f "+chpcAligner.getGenomeIndex() +" /scratch/local/"+rndName+"/");
		//uncompress it?
		//if (chpcAligner.getGenomeIndex().getName().endsWith(".gz")){
			//script.add("gunzip -f "+chpcAligner.getGenomeIndex() +" /scratch/local/"+rndName+"/");
		//}
		//copy over the data
		script.add("mkdir /scratch/local/"+rndName+"/Data1/ /scratch/local/"+rndName+"/Data2/");
		script.add("cp -f "+data1 +" /scratch/local/"+rndName+"/Data1/");
		if (data2 != null) script.add("cp -f "+data2 +" /scratch/local/"+rndName+"/Data2/");
		//add alignment command
		script.add(fetchNovoAlignmentCommand(rndName));
		//force compress results, deletes any that already exist
		script.add("gzip -f /scratch/local/"+rndName+"/"+results.getName());
		//move results off node 
		script.add("mv -f /scratch/local/"+rndName+"/"+results.getName()+".gz "+gzippedResults);
		//do a clean up 
		script.add("rm -rf /scratch/local/"+rndName+"/");
		//end time
		script.add("echo -n 'End:\t'; date +%s");
		//write out
		IO.writeArrayList(script, shellScript);
		//change permissions so it can be executed
		shellScript.setExecutable(true);
	}
	
	/**Returns the minutes of run time or -1 if log file not present or -2 if the log file didn't contain start and end times*/
	public int getRunTime(){
		int runTime = -1;
		if (log.exists()){
			HashMap<String,String> lines = IO.loadFileIntoHashMap(log);
			if (lines == null || lines.containsKey("Start:")== false || lines.containsKey("End:")==false) runTime = -2;
			else {
				int startSec = Integer.parseInt(lines.get("Start:"));
				int endSec = Integer.parseInt(lines.get("End:"));
				runTime = endSec- startSec;
			}
		}
		return runTime/60;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(name);
		sb.append("\t");
		sb.append(jobID);
		sb.append("\t");
		sb.append(log.getName());
		sb.append("\t");
		sb.append(gzippedResults.getName());
		sb.append("\t");
		sb.append(data1.getName());
		sb.append("\t");
		if (data2 != null){
			sb.append(data2.getName());
			sb.append("\t");
		}
		sb.append(shellScript.getName());
		return sb.toString();
	}
	
	public boolean isGzippedResultFileOK(){
		if (gzippedResults.exists() == false) return false;
		if (gzippedResults.length() < 5000) return false;
		return true;
	}
	
	public boolean isLogFilePresent(){
		return log.exists();
	}

	public String fetchNovoAlignmentCommand(String dirName){
		String cmd = null;
		//for novo
		if (chpcAligner.getAlignerName().contains("novoalign")){
			String data = "/scratch/local/"+dirName+"/Data1/"+data1.getName();
			//paired?
			if (data2 != null) data = data +" /scratch/local/"+dirName+"/Data2/"+data2.getName();
			cmd = 
				chpcAligner.getAlignerApp() + " "+
				chpcAligner.getAlignerParams() + " "+
				" -d /scratch/local/"+dirName+"/"+ chpcAligner.getGenomeIndex().getName() +" "+
				" -f "+data;
			String stripHeaders = "";
			if (chpcAligner.isStripSAMSQHeaders()) stripHeaders = " | grep -v ^@SQ ";
				
			String stripNonChr = "";
			if (filterForChrLines) stripNonChr =" | grep chr ";
		
			cmd = cmd + stripHeaders + stripNonChr +" > /scratch/local/"+dirName+"/"+ results.getName()+" \n";
		}
		else {
			chpcAligner.printError("Unsupported aligner "+chpcAligner.getAlignerName());
		}
		return cmd;
	}

	public String fetchPBSHeader(){
		//add ann account?
		String account = "";
		if (chpcAligner.getChpcAccount() != null) account = "#PBS -A "+chpcAligner.getChpcAccount()+"\n";

		String h = 
			"#PBS -l nodes=1:ppn="+chpcAligner.getNumberCPUs()+",walltime="+chpcAligner.getHrsWallTime()+":00:00 \n"+
			"#PBS -m a \n"+ 
			"#PBS -M "+chpcAligner.getEmail()+"\n"+ 
			"#PBS -N "+name+"\n"+  
			"#PBS -j oe \n"+ 
			account+
			"#PBS -o "+ log +"\n";
		return h;
	}
	
	
	public File getShellScript() {
		return shellScript;
	}
	public void setShellScript(File shellScript) {
		this.shellScript = shellScript;
	}
	public File getResults() {
		return results;
	}
	public void setResults(File results) {
		this.results = results;
		gzippedResults = new File (results+".gz");
	}
	public File getData1() {
		return data1;
	}
	public void setData1(File data) {
		this.data1 = data;
	}
	public File getLog() {
		return log;
	}
	public void setLog(File log) {
		this.log = log;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}

	public CHPCAligner getChpcAligner() {
		return chpcAligner;
	}

	public void setChpcAligner(CHPCAligner chpcAligner) {
		this.chpcAligner = chpcAligner;
	}

	public ArrayList<String> getScript() {
		return script;
	}

	public void setJobID(String jobID) {
		this.jobID = jobID;
	}

	public String getJobID() {
		return jobID;
	}

	public File getGzippedResults() {
		return gzippedResults;
	}

	public File getData2() {
		return data2;
	}

	public void setData2(File data2) {
		this.data2 = data2;
	}

	public boolean isFilterForChrLines() {
		return filterForChrLines;
	}

	public void setFilterForChrLines(boolean filterForChrLines) {
		this.filterForChrLines = filterForChrLines;
	}

}
