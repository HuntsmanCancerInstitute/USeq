package edu.utah.pysano.commands;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;

import edu.utah.pysano.utils.Logger;

public class CommandMetrics extends Command {
	public static ArrayList<ArrayList<Command>> checkForInputs(File projectDirectory, HashMap<String, String> properties, Logger logFile, String studyName, File targetFile, 
			LinkedHashMap<String,Boolean> analysisSteps, String cluster, String email, Integer wallTime, boolean suppress) {
		//Create input matching patterns
		Pattern pattern1 = Pattern.compile("(.+?).final.bam");
		Pattern pattern2 = Pattern.compile("(.+?).final.bai");
		
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
				
				CommandMetrics cmt = new CommandMetrics(inputFiles, projectDirectory, logFile, studyName, sampleName, properties, targetFile, cluster, email, wallTime, suppress);
		
				commandList.add(cmt);
			}
		}	
		
		if (commandList.isEmpty()) {
			logFile.writeErrorMessage("[CommandMetrics] No matching files found for metrics calling, exiting.",false);
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
	public CommandMetrics(ArrayList<File> inputFiles, File projectDirectory,  Logger logFile, String studyName, String sampleName, HashMap<String, String> properties, File targetFile, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		initialize(projectDirectory, sampleName, properties, studyName, targetFile);
		addInputFiles(inputFiles);
		for (File inputFile: inputFiles) {
			addProtect(inputFile);
		}
		
		for (File deleteFile: inputFiles) {
			addDelete(deleteFile);
		}
		
	}

	//Constructor for combined command (don't check for input files pre-run)
	public CommandMetrics(File projectDirectory, Logger logFile, String studyName, String sampleName, HashMap<String, String> properties, File targetFile, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		initialize(projectDirectory, sampleName, properties, studyName, targetFile);
	}

	private void initialize(File projectDirectory, String sampleName, HashMap<String, String> properties, String studyName, File targetFile) {
		//Setup basic information
		setGlobalProperties(properties);
		File finalDirectory = new File(projectDirectory,"Metrics");
		File jobDirectory = new File(projectDirectory,"JOB_metrics_" + sampleName);
		setFinalDirectory(finalDirectory);
		setJobDirectory(jobDirectory);
		setSampleName(sampleName);
		
		//setup external
		if (targetFile == null) {
			targetFile = new File(properties.get("TARGET_DEFAULT"));
		}
		addExternalFile(targetFile);
		
		//setup delete
		addDelete(targetFile);
		
		//Setup output files
		File dictFile = new File(finalDirectory,sampleName + ".dict.txt");
		File imageDir = new File(finalDirectory,"images");
		File excelFile = new File(finalDirectory,studyName + ".xlsx");
		
		addOutputFile(dictFile);
		addOutputFile(imageDir);
		addConditionalSkip(excelFile, dictFile);
		
		//Setup properties
		addLocalProperty("NAME", sampleName);
		addLocalProperty("TARGETS", targetFile.getName());

		//Setup template
		String version = properties.get("VERSION");
		setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/metrics." + version + ".txt"));


		//setup protect
		addProtect(new File(jobDirectory,"cmd.txt"));
		addProtect(targetFile);

	}
	
	public static void mergeMetrics(File projectDirectory, HashMap<String,String> properties, String studyName, Logger logFile) {
		try {
			
			File metricsDir = new File(projectDirectory,"Metrics");
			File finalFile = new File(metricsDir,studyName + ".xlsx");
			if (finalFile.exists()) {
				return;
			}
			
			
			ProcessBuilder pb = new ProcessBuilder("java","-Xmx10g","-jar",properties.get("USEQ_PATH_LOCAL") + "/Apps/MergeExonMetrics","-f",metricsDir.getAbsolutePath(),"-o",studyName);
			Process p = pb.start();
			
			
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[CommandMetrics] Could not merge metrics files",true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[CommandMetrics] IO Exception while trying to merge metrics",true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[CommandMetrics] Process was interrupted while trying to merge metrics",true);
			System.exit(1);
		}
	}
}
