package edu.utah.pysano.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.pysano.utils.Logger;

public class CommandPostProcess extends Command {
	public static ArrayList<ArrayList<Command>> checkForInputs(File projectDirectory, HashMap<String, String> properties, Logger logFile, 
			LinkedHashMap<String,Boolean> analysisSteps, String cluster, String email, Integer wallTime, boolean suppress) {
		//Create input matching patterns
		Pattern pattern1 = Pattern.compile("(.+?).raw.bam");
		Pattern pattern2 = Pattern.compile("(.+?).raw.bai");
		
		//Check final directory first, then project directory for matches 
		File alignmentDirectory = new File(projectDirectory,"Alignments");
		HashMap<String,File> matchingFiles1 = new HashMap<String, File>();
		HashMap<String,File> matchingFiles2 = new HashMap<String, File>();
		if (alignmentDirectory.exists()) {
			logFile.writeInfoMessage("Looking for matches in: " + alignmentDirectory.getAbsolutePath());
			File[] finalDirContents = alignmentDirectory.listFiles();
			matchingFiles1 = matchFiles(pattern1,finalDirContents);
			matchingFiles2 = matchFiles(pattern2,finalDirContents);
		}
		if (matchingFiles1.isEmpty() || matchingFiles2.isEmpty()) {
			logFile.writeInfoMessage("No matching files in " + alignmentDirectory.getAbsolutePath() + ", checking for matches: " + projectDirectory.getAbsolutePath());
			File[] projectDirContents = projectDirectory.listFiles();
			matchingFiles1 = matchFiles(pattern1,projectDirContents);
			matchingFiles2 = matchFiles(pattern2,projectDirContents);
		}
	
		//Find samples with both files and initialize command
		ArrayList<Command> commandList = new ArrayList<Command>();
		for (String key: matchingFiles1.keySet()) {
			if (matchingFiles2.containsKey(key)) { 
				//Required files are found, determine names
				String sampleName = key;
				
				ArrayList<File> inputFiles = new ArrayList<File>();
				inputFiles.add(matchingFiles1.get(key));
				inputFiles.add(matchingFiles2.get(key));
				
				CommandPostProcess cpp = new CommandPostProcess(inputFiles, projectDirectory, logFile, sampleName, properties, cluster,email,wallTime,suppress);
		
				commandList.add(cpp);
			}
		}	
		
		if (commandList.isEmpty()) {
			logFile.writeErrorMessage("[CommandPostProcess] No matching files found for post processing, exiting.",false);
			System.exit(1);
		}
		
		//Determine if command finished, needs cleaning or needs to be run
		ArrayList<Command> toRun = new ArrayList<Command>();
		ArrayList<Command> toClean = new ArrayList<Command>();
		for (Command c: commandList) {
			if (c.doesFinalExist()) {
				//Finished!
			} else if (c.doesWorkingExist()) {
				//Just cleanup
				toClean.add(c);
			} else {
				//everything
				toClean.add(c);
				toRun.add(c);
			}
		}
		
		//Return list of files to run/clean
		ArrayList<ArrayList<Command>> returnList = new ArrayList<ArrayList<Command>>();
		returnList.add(toRun);
		returnList.add(toClean);
		
		return returnList;
	}

	//Constructor for stand alone command
	public CommandPostProcess(ArrayList<File> inputFiles, File projectDirectory, Logger logFile, String sampleName, HashMap<String, String> properties, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		initialize(projectDirectory, sampleName, properties);
		addInputFiles(inputFiles);
		
		for (File inputFile: inputFiles) {
			addProtect(inputFile);
		}
		
		//setup delete
		for (File deleteFile: inputFiles) {
			addDelete(deleteFile);
		}
		
		
	}

	//Constructor for combined command (don't check for input files pre-run)
	public CommandPostProcess(File projectDirectory, Logger logFile, String sampleName, HashMap<String, String> properties, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		initialize(projectDirectory, sampleName, properties);
	}
	
	private void initialize(File projectDirectory, String sampleName, HashMap<String, String> properties) {
		//Setup basic information
		setGlobalProperties(properties);
		File finalDirectory = new File(projectDirectory,"Alignments");
		File jobDirectory = new File(projectDirectory,"JOB_postprocess_" + sampleName);
		setFinalDirectory(finalDirectory);
		setJobDirectory(jobDirectory);
		setSampleName(sampleName);
		
		//Setup output files
		addOutputFile(new File(finalDirectory,sampleName + ".final.bam"));
		addOutputFile(new File(finalDirectory,sampleName + ".final.bai"));
		
		//Setup skip files (these files will be gone by this point, so shoudn't be checked
		addSkipFile(new File(finalDirectory,sampleName + ".raw.bam"));
		addSkipFile(new File(finalDirectory,sampleName + ".raw.bai"));
		
		//Setup properties
		addLocalProperty("NAME", sampleName);
		
		//Setup template
		String version = properties.get("VERSION");
		setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/align_postprocess." + version + ".txt"));
		
		File commandFile = new File(jobDirectory,"cmd.txt");
		addProtect(commandFile);
	}
}
