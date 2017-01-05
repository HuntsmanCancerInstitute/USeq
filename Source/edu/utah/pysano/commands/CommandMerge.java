package edu.utah.pysano.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;

import edu.utah.pysano.utils.Logger;

public class CommandMerge extends Command {
	public static ArrayList<ArrayList<Command>> checkForInputs(File projectDirectory, HashMap<String, String> properties, Logger logFile, String studyName, File targetFile, 
			LinkedHashMap<String,Boolean> analysisSteps, String cluster, String email, Integer wallTime, boolean suppress) {
		Pattern pattern1 = Pattern.compile("(.+?).raw.bam");
		Pattern pattern2 = Pattern.compile("(.+?).raw.bai");
		
		//Check final directory first, then project directory for matches 
		File alignmentDirectory = new File(projectDirectory,"Alignments");
		HashMap<String,File> matchingFiles1 = new HashMap<String, File>();
		HashMap<String,File> matchingFiles2 = new HashMap<String, File>();
		if (alignmentDirectory.exists()) {
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
	
		
		//Groups input files by sample names
		HashMap<String,ArrayList<File>> sampleNameToInputs = new HashMap<String,ArrayList<File>>();
		HashMap<String,ArrayList<String>> sampleNameToSampleIDs = new HashMap<String,ArrayList<String>>();
			
		//Group samples
		for (String key: matchingFiles1.keySet()) {
			if (matchingFiles2.containsKey(key)) { 
				//Required files are found, determine names
				String sampleID = key;
				String sampleSubID = key;
				String sampleName = null;
				
				String[] splitOnName = key.split("__");
				
				if (splitOnName.length > 1) {
					sampleName = splitOnName[0];
					sampleSubID = splitOnName[1];
				} 
				
				String[] idParts = sampleSubID.split("_",2);
				
				if (idParts.length != 2) {
			    	logFile.writeErrorMessage("[CommandMerge] File name improperly formed. Should have an underscore separating the sample name and the flowcell name: SAMPLE_FLOWCELL_1.fastq and SAMPLE_FLOWCELL_2.fastq",false);
			    	System.exit(1);
			    } 
				
				if (sampleName == null) {
					sampleName = idParts[0];
				}
				
				if (!sampleNameToInputs.containsKey(sampleName)) {
					sampleNameToInputs.put(sampleName, new ArrayList<File>());
					sampleNameToSampleIDs.put(sampleName, new ArrayList<String>());
				}
				sampleNameToInputs.get(sampleName).add(matchingFiles1.get(key));
				sampleNameToInputs.get(sampleName).add(matchingFiles2.get(key));
				sampleNameToSampleIDs.get(sampleName).add(sampleID);
			
			}
		}	
	
		//set this to true if at least one merge found
		boolean foundMerged = false;
	
		//Initialize commandList
		ArrayList<Command> commandList = new ArrayList<Command>();
		
		//Iterate through each sampleName and create commands as needed
		for (String sampleName: sampleNameToInputs.keySet()) {
			ArrayList<File> inputList = sampleNameToInputs.get(sampleName);
			if (inputList.size() > 2) {
				CommandMerge cm = new CommandMerge(inputList, projectDirectory, logFile, sampleName, sampleNameToSampleIDs.get(sampleName), properties, cluster, email, wallTime, suppress);
				commandList.add(cm);
				foundMerged = true;
				//If no merging, add downstream steps if full run!
				if (analysisSteps.containsKey("postprocess")) {
					CommandPostProcess cpp = new CommandPostProcess(projectDirectory, logFile, cm.getSampleName(), properties, cluster, email, wallTime, suppress);
					cm.merge(cpp);	
				}
				if (analysisSteps.containsKey("metrics")) {
					CommandMetrics cmt = new CommandMetrics(projectDirectory, logFile, studyName, cm.getSampleName(), properties, targetFile, cluster, email, wallTime, suppress);
					cm.merge(cmt);
				}
				if (analysisSteps.containsKey("haplotype")) {
					CommandHaplotype ch = new CommandHaplotype(projectDirectory, logFile, cm.getSampleName(), properties, targetFile, cluster, email,wallTime,suppress);
					cm.merge(ch);
				}
			}
		}
		
		//If nothing can be merged throw a warning message!
		if (commandList.isEmpty()) {
			logFile.writeErrorMessage("[CommandMerge] No matching files found for merging, exiting.",false);
			System.exit(1);
		}
		
		if (!foundMerged) {
			logFile.writeWarningMessage("[CommandMerge] Could not find any data to merge, this step should not have been called!  Please contact core with this message!");
		}
		
		//Update analysis steps
		if (analysisSteps.containsKey("postprocess")) {
			analysisSteps.put("postprocess",false);	
		}
		if (analysisSteps.containsKey("metrics")) {
			analysisSteps.put("metrics",false);	
		}
		if (analysisSteps.containsKey("haplotype")) {
			analysisSteps.put("haplotype",false);	
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
		
		ArrayList<ArrayList<Command>> returnList = new ArrayList<ArrayList<Command>>();
		returnList.add(toRun);
		returnList.add(toClean);
		
		return returnList;
	}
	
	
	public CommandMerge(ArrayList<File> inputFiles, File projectDirectory, Logger logFile, String sampleName, ArrayList<String> sampleIdList, HashMap<String, String> properties, String cluster, String email, 
			Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		
		//Setup basic information
		setGlobalProperties(properties);
		File finalDirectory = new File(projectDirectory,"Alignments");
		File jobDirectory = new File(projectDirectory,"JOB_merge_" + sampleName);
		setFinalDirectory(finalDirectory);
		setJobDirectory(jobDirectory);
		setSampleName(sampleName);
		
		//setup input files
		addInputFiles(inputFiles);
		
		for (File inputFile: inputFiles) {
			addProtect(inputFile);
		}
		
		//setup delete
		for (File deleteFile: inputFiles) {
			addDelete(deleteFile);
		}
		
		//Setup output files
		addOutputFile(new File(finalDirectory,sampleName + ".raw.bam"));
		addOutputFile(new File(finalDirectory,sampleName + ".raw.bai"));
		
		//Mark input files to be skipped and deleted
		for (String sampleID: sampleIdList) {
			addSkipFile(new File(finalDirectory,sampleID + ".raw.bam"));
			addSkipFile(new File(finalDirectory,sampleID + ".raw.bai"));
			addDelete(new File(finalDirectory,sampleID + ".raw.bam"));
			addDelete(new File(finalDirectory,sampleID + ".raw.bai"));
		}
		
		
		//Setup properties
		StringBuffer inputFileList = new StringBuffer();
		for (File inputFile: inputFiles) {
			if (inputFile.getName().endsWith(".bam")) {
				inputFileList.append(inputFile.getName() + " ");
			}
			
		}
		
		addLocalProperty("OUTPUT", sampleName + ".raw.bam");
		addLocalProperty("INPUT_LIST", inputFileList.toString());
		
		
		//Setup template
		String version = properties.get("VERSION");
		setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/align_merge." + version + ".txt"));
		
		addProtect(new File(jobDirectory,"cmd.txt"));
		
	}
}
