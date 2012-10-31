package util.apps;
import java.io.*;
import java.util.*;
import util.gen.IO;
import util.gen.Misc;
/**
 * Fires command lines on clusters through qsub. Logs commands. Directs results into a daily directory.
 */
public class JQSub {
	//fields
	private File paramFile;
	private File logFile;
	private File logFileDirectory;
	//defaults
	private String pathToJava32Bit = "/nfs/linux/pkg/java_pkg/jdk1.4/bin/java";
	private String pathToJava64Bit = "/nfs/linux/pkg/java_pkg/jdk1.5.0-amd64/bin/java";
	private File tempDirectory = new File ("/home/dnix/JQSub/");
	private String queue32Bit = "tome-n@pbs1.cluster.ev.affymetrix.com";
	private String queue64Bit = "general-n@torque1.cluster.ev.affymetrix.com";
	private String pathToBioTools = "/home/dnix/Code/BioTools/";
	private String walltime = "72:00:00";
	
	public JQSub(String[] args){
		//load param file from TiMAT directory
		loadParamFile();
		//initialize params
		int num = args.length;
		StringBuffer command = new StringBuffer();
		StringBuffer commandsToSave = new StringBuffer();
		String qsubParams = "-q "+queue32Bit+" -l nodes=1,walltime="+walltime+" ";
		//does logFileDirectory exist?
		if (logFileDirectory.exists()==false) {
			logFileDirectory.mkdir();
		}
		//run through args
		for (int i=0; i<num; i++){
			//look for return and java
			if (args[i].equalsIgnoreCase("return")) {
				command.append(";");
				commandsToSave.append("; ");
			}
			else if (args[i].equalsIgnoreCase("\\>")) {
				command.append(">");
				commandsToSave.append("> ");
			}
			//look for wall time
			else if (args[i].equalsIgnoreCase("-w")){
				i++;
				double fractionalTime = Double.parseDouble(args[i]);
				walltime = Misc.getFormattedTimeFromFraction(fractionalTime);
			}
			else if (args[i].equals("java")) {
				command.append("cd "+pathToBioTools+"; ");
				commandsToSave.append("java ");
				//look for memory
				if (args[i+1].startsWith("-Xmx")){
					i++;
					commandsToSave.append(args[i]+ " ");
					//extract memory requirement
					int memory = Integer.parseInt(args[i].replaceAll("[^\\d]",""));
					if (memory <= 1500) {
						command.append(pathToJava32Bit+" " +args[i]+" ");
						qsubParams = "-q "+queue32Bit+" -l nodes=1,walltime="+walltime+" ";
					}
					else if (memory <= 4000) {
						command.append(pathToJava64Bit+" " +args[i]+" ");
						qsubParams = "-q "+queue64Bit+" -l nodes=1,walltime="+walltime+" ";
					}
					else {
						command.append(pathToJava64Bit+" " +args[i]+" ");
						qsubParams = "-q "+queue64Bit+" -l nodes=1:ppn=1:MEM16G,walltime="+walltime+" ";
					}
				}
				else {
					command.append(pathToJava32Bit+" -Xmx1500M "); 
				}
			}
			else {
				command.append(args[i]); 
				command.append(" ");
				commandsToSave.append(args[i]);
				commandsToSave.append(" ");
			}
		}
		//System.out.println("Executing:\n"+qsubParams+"\n"+command.toString()+"\n");
		String results  = IO.fireQSub(qsubParams, command.toString(), logFileDirectory);
		//System.out.println("\nResults: "+results+"\n");
		String line = results+"\t"+commandsToSave +"\n";
		System.out.print("Submitted: "+line);
		if (logFile.exists()) IO.writeStringAppend(line, logFile);
		else {
			String header = "Name\tJobID\tCommand\n"; 
			IO.writeString(header + line, logFile);
		}
		
	}
	
	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}
		new JQSub(args);
	}
	
	public void loadParamFile(){
		
		File homeDir = new File (System.getProperty("user.home"));
		paramFile = new File(homeDir,".paramFileJQSub.txt");
		if (paramFile.exists()){
			LinkedHashMap map = IO.loadKeyValueFile(paramFile);
			if (map.containsKey("pathToJava32Bit")) pathToJava32Bit = (String)map.get("pathToJava32Bit");
			if (map.containsKey("pathToJava64Bit")) pathToJava64Bit = (String)map.get("pathToJava64Bit");
			if (map.containsKey("tempDirectory")) tempDirectory = new File((String)map.get("tempDirectory"));
			if (map.containsKey("queue32Bit")) queue32Bit = (String)map.get("queue32Bit");
			if (map.containsKey("queue64Bit")) queue64Bit = (String)map.get("queue64Bit"); 
			if (map.containsKey("pathToBioTools")) pathToBioTools = (String)map.get("pathToBioTools");
		}
		else IO.writeString(defaults, paramFile);
		
		logFile = new File (tempDirectory, Misc.getDateNoSpaces()+"LogFile.xls");
		logFileDirectory = new File (tempDirectory, Misc.getDateNoSpaces()+"Results");
	}
	
	public String defaults = 
		"//parameters for using JQSub, use full paths, spaces are ignored\n"+
		"\npathToJava32Bit = "+pathToJava32Bit+
		"\npathToJava64Bit = "+pathToJava64Bit+
		"\ntempDirectory = "+tempDirectory+
		"\nqueue32Bit = "+queue32Bit+
		"\nqueue64Bit = "+queue64Bit+
		"\npathToBioTools = "+pathToBioTools; 
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 JQSub: April 2007                                **\n" +
				"**************************************************************************************\n" +
				"JQSub executes a given command line on the cluster logging your submission and results\n" +
				"in a temp directory designated by the ~/.paramFileJQSub.txt in your home directory.\n" +
				"Use the word 'return' to separate multiple shell commands. JQSub is java friendly and\n" +
				"will direct your job to the appropriate cluster depending on your -Xmx parameter.\n\n" +
				
				"Options:\n"+
				"-w Wall time in hours, defaults to 72.  Set to 10-15% longer than expected runtime\n"+
				"     for your application.  Otherwise it gets the chop.\n\n"+

				"Example: java util/apps/JQSub -w 1.75 java -Xmx4000M trans/main/ScanChip\n\n" +
				
				"**************************************************************************************\n");
	}

	
}

