package edu.utah.pysano.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;

import edu.utah.pysano.utils.Logger;

public class CommandHaplotype extends Command {
	public static ArrayList<ArrayList<Command>> checkForInputs(File projectDirectory, HashMap<String, String> properties, Logger logFile, File targetFile, 
			LinkedHashMap<String,Boolean> analysisSteps, String cluster, String email, Integer wallTime, boolean suppress) {
		//Create input matching patterns
		Pattern pattern1 = Pattern.compile("(.+?).final.bam");
		Pattern pattern2 = Pattern.compile("(.+?).final.bai");
		
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
		
		//Find samples with both files and initialize command
		ArrayList<Command> commandList = new ArrayList<Command>();
		for (String key: matchingFiles1.keySet()) {
			if (matchingFiles2.containsKey(key)) { 
				//Required files are found, determine names
				String sampleName = key;
				
				ArrayList<File> inputFiles = new ArrayList<File>();
				inputFiles.add(matchingFiles1.get(key));
				inputFiles.add(matchingFiles2.get(key));
				
				CommandHaplotype ch = new CommandHaplotype(inputFiles, projectDirectory, logFile, sampleName, properties, targetFile,cluster, email,wallTime, suppress);
		
				commandList.add(ch);
			}
		}
		
		//Write message if nothing found
		if (commandList.isEmpty()) {
			logFile.writeErrorMessage("[CommandHaplotype] No matching files found for halotype calling, exiting.",false);
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
	public CommandHaplotype(ArrayList<File> inputFiles, File projectDirectory, Logger logFile, String sampleName, HashMap<String, String> properties, File targetFile, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		initialize(projectDirectory, sampleName, properties, targetFile);
		addInputFiles(inputFiles);
		for (File inputFile: inputFiles) {
			addProtect(inputFile);
		}
		for (File deleteFile: inputFiles) {
			addDelete(deleteFile);
		}
		
	}

	//Constructor for combined command (don't check for input files pre-run)
	public CommandHaplotype(File projectDirectory, Logger logFile, String sampleName, HashMap<String, String> properties, File targetFile, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		initialize(projectDirectory, sampleName, properties, targetFile);
	}

	private void initialize(File projectDirectory, String sampleName, HashMap<String, String> properties, File targetFile) {
		//Setup basic information
		setGlobalProperties(properties);
		File finalDirectory = new File(projectDirectory,"Variants/GVCF");
		File jobDirectory = new File(projectDirectory,"JOB_haplotype_" + sampleName);
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
		addOutputFile(new File(finalDirectory,sampleName + ".genomic.vcf.gz"));
		addOutputFile(new File(finalDirectory,sampleName + ".genomic.vcf.gz.tbi"));

		//Setup properties
		addLocalProperty("SAMPLE", sampleName);
		addLocalProperty("TARGETS", "-L " + targetFile.getName());

		//Setup template
		String version = properties.get("VERSION");
		setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/variant_haplotype." + version + ".txt"));

		File commandFile = new File(jobDirectory,"cmd.txt");
		
		//setup protect
		addProtect(commandFile);
		addProtect(targetFile);
	}
}
