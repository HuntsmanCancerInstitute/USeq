package edu.utah.pysano.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;

import edu.utah.pysano.utils.Logger;

public class CommandAlign extends Command {
	
	public static ArrayList<ArrayList<Command>> checkForInputs(File projectDirectory, HashMap<String, String> properties, Logger logFile, String studyName,
			boolean isNovo, boolean validateFastq, File targetFile, LinkedHashMap<String,Boolean> analysisSteps, String cluster, String email, Integer wallTime, boolean suppress) {
		Pattern pattern1 = Pattern.compile("(.+?)(_1|_R1|_001)\\..+\\.gz$");
		Pattern pattern2 = Pattern.compile("(.+?)(_2|_R2|_002)\\..+\\.gz$");
		
		//Get file listing of root directory
		File[] projectDirContents = projectDirectory.listFiles();
		HashMap<String,File> matchingFiles1 = matchFiles(pattern1,projectDirContents);
		HashMap<String,File> matchingFiles2 = matchFiles(pattern2,projectDirContents);
		
		//Tracks sample names to sampleIDs
		HashMap<String,Integer> sampleNames = new HashMap<String,Integer>();
			
		//Find samples with both files and initialize command
		ArrayList<Command> commandList = new ArrayList<Command>();
		for (String key: matchingFiles1.keySet()) {
			if (matchingFiles2.containsKey(key)) { 
				//Required files are found, determine names
				String sampleID = key;
				String sampleSubID = key;
				String library = "NA";
				String flowcell = "NA";
				String sampleName = null;
				
				String[] splitOnName = key.split("__");
				
				if (splitOnName.length > 1) {
					sampleName = splitOnName[0];
					sampleSubID = splitOnName[1];
				} 
				
				String[] idParts = sampleSubID.split("_",2);
				
				if (idParts.length != 2) {
			    	logFile.writeErrorMessage("[CommandAlignNovo] File name improperly formed. Should have an underscore separating the sample name and the flowcell name: SAMPLE_FLOWCELL_1.fastq and SAMPLE_FLOWCELL_2.fastq",false);
			    	System.exit(1);
			    } 
				
				library = idParts[0];
				flowcell= idParts[1];
				if (sampleName == null) {
					sampleName = library;
				}
				
				if (validateFastq) {
					validateFastqSet(matchingFiles1.get(key),matchingFiles2.get(key), logFile);
				}
				
				ArrayList<File> inputFiles = new ArrayList<File>();
				inputFiles.add(matchingFiles1.get(key));
				inputFiles.add(matchingFiles2.get(key));
				
				CommandAlign can = new CommandAlign(inputFiles, projectDirectory, logFile, sampleName, sampleID, flowcell, library, properties, isNovo, cluster, email,wallTime,suppress);
		
				commandList.add(can);
				
				//Keep track of samples
				if (!sampleNames.containsKey(sampleName)) {
					sampleNames.put(sampleName, 0);
				}
				sampleNames.put(sampleName, sampleNames.get(sampleName) + 1);
			}
		}	
		
		//Write message if nothing found
		if (commandList.isEmpty()) {
			logFile.writeErrorMessage("[CommandAlign] No matching files found for alignment, moving to next step.",false);
			System.exit(1);
		}
		
		//Figure out if any merging is needed!
		boolean mergeNeeded = false;
		for (Command c: commandList) {
			if (sampleNames.get(c.getSampleName()) > 1) {
				mergeNeeded = true;
			} else {
				CommandRename cr = new CommandRename(projectDirectory, logFile, c.getSampleName(), c.getSampleId(), properties, cluster, email, wallTime, suppress);
				c.merge(cr);
				//If no merging, add downstream steps if full run!
				if (analysisSteps.containsKey("postprocess")) {
					CommandPostProcess cpp = new CommandPostProcess(projectDirectory, logFile, c.getSampleName(), properties, cluster, email, wallTime, suppress);
					c.merge(cpp);	
				}
				if (analysisSteps.containsKey("metrics")) {
					CommandMetrics cmt = new CommandMetrics(projectDirectory, logFile, studyName, c.getSampleName(), properties, targetFile, cluster, email, wallTime, suppress);
					c.merge(cmt);
				}
				if (analysisSteps.containsKey("haplotype")) {
					CommandHaplotype ch = new CommandHaplotype(projectDirectory, logFile, c.getSampleName(), properties, targetFile, cluster, email,wallTime,suppress);
					c.merge(ch);
				}
			}
		}
		
		
		//Make changes to run order depending on run/merge type
		if (!mergeNeeded) {
			analysisSteps.put("merge", false);
			//If full, modify steps
			if (analysisSteps.containsKey("postprocess")) {
				analysisSteps.put("postprocess",false);	
			}
			if (analysisSteps.containsKey("metrics")) {
				analysisSteps.put("metrics",false);	
			}
			if (analysisSteps.containsKey("haplotype")) {
				analysisSteps.put("haplotype",false);	
			}
			
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

	public CommandAlign(ArrayList<File> inputFiles, File projectDirectory, Logger logFile, String sampleName, String sampleID, String flowcell, String library, HashMap<String, String> properties, 
			boolean isNovo, String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		//Setup basic information
		setGlobalProperties(properties);
		
		File finalDirectory = new File(projectDirectory,"Alignments");
		File jobDirectory = new File(projectDirectory,"JOB_align_" + sampleID);
		setFinalDirectory(finalDirectory);
		setJobDirectory(jobDirectory);
		setSampleID(sampleID);
		setSampleName(sampleName);
		
		//setup input files
		addInputFiles(inputFiles);
		
		for (File inputFile: inputFiles) {
			addProtect(inputFile);
		}
		
		//setup delete
		for (File deleteFile: inputFiles) {
			addProtect(deleteFile);
			addDelete(deleteFile);
		}
		
		//Setup output files
		addOutputFile(new File(finalDirectory,sampleID + ".raw.bam"));
		addOutputFile(new File(finalDirectory,sampleID + ".raw.bai"));
		
		//Setup properties
		addLocalProperty("NAME", sampleID);
		addLocalProperty("LIBRARY", library);
		addLocalProperty("SAMPLE", sampleName);
		addLocalProperty("FLOWCELL", flowcell);
		
		//Setup template
		String version = properties.get("VERSION");
		if (isNovo) {
			setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/align_novo." + version + ".txt"));
		} else {
			setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/align_bwa." + version + ".txt"));
		}
		
		File commandFile = new File(jobDirectory,"cmd.txt");
		addProtect(commandFile);
		
	}
	

}
		
		



