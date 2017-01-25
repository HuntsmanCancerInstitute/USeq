package edu.utah.pysano.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Pattern;

import edu.utah.pysano.utils.Logger;

public class CommandVariantMerge extends Command {
	
	public static ArrayList<ArrayList<Command>> checkForInputs(File projectDirectory, HashMap<String, String> properties, Logger logFile, File targetFile, String studyName,  boolean isVQSR,
			boolean use1KGenomes, boolean callGenotypes, LinkedHashMap<String,Boolean> analysisSteps, String cluster, String email, Integer wallTime, boolean suppress) {
		//Create input matching patterns
		Pattern pattern1 = Pattern.compile("(.+?).genomic.vcf.gz");
		Pattern pattern2 = Pattern.compile("(.+?).genomic.vcf.gz.tbi");
		
		//Check final directory first, then project directory for matches 
		File variantDirectory = new File(projectDirectory,"Variants/GVCF");
		HashMap<String,File> matchingFiles1 = new HashMap<String, File>();
		HashMap<String,File> matchingFiles2 = new HashMap<String, File>();
		if (variantDirectory.exists()) {
			File[] finalDirContents = variantDirectory.listFiles();
			matchingFiles1 = matchFiles(pattern1,finalDirContents);
			matchingFiles2 = matchFiles(pattern2,finalDirContents);
		}
		if (matchingFiles1.isEmpty() || matchingFiles2.isEmpty()) {
			logFile.writeInfoMessage("No matching files in " + variantDirectory.getAbsolutePath() + ", checking for matches: " + projectDirectory.getAbsolutePath());
			File[] projectDirContents = projectDirectory.listFiles();
			matchingFiles1 = matchFiles(pattern1,projectDirContents);
			matchingFiles2 = matchFiles(pattern2,projectDirContents);
		}
	
		//Find samples with both files and initialize command
		ArrayList<File> inputFileList = new ArrayList<File>();
		for (String key: matchingFiles1.keySet()) {
			if (matchingFiles2.containsKey(key)) { 
				//Required files are found, determine names				
				inputFileList.add(matchingFiles1.get(key));
				inputFileList.add(matchingFiles2.get(key));
			}
		}
		
		ArrayList<Command> commandList = new ArrayList<Command>();
		if (inputFileList.isEmpty()) { //Nothing found
			logFile.writeErrorMessage("[CommandVariantMerge] No matching files found for variant merging, exiting.",false);
			System.exit(1);
		}
		
		CommandVariantMerge cvm = new CommandVariantMerge(inputFileList, projectDirectory, logFile, studyName, properties, cluster, email, wallTime, suppress);
		commandList.add(cvm);
		if (callGenotypes) {
			CommandGenotype cvg = new CommandGenotype(projectDirectory, properties, logFile, targetFile, studyName, isVQSR, use1KGenomes, inputFileList.size()/2, cluster, email, wallTime, suppress);
			cvm.merge(cvg);
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
	public CommandVariantMerge(ArrayList<File> inputFiles, File projectDirectory, Logger logFile, String studyName, HashMap<String, String> properties, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		//Setup basic information
		setGlobalProperties(properties);
		File finalDirectory = new File(projectDirectory,"Variants/GVCF");
		File jobDirectory = new File(projectDirectory,"JOB_genotype_" + studyName);
		setFinalDirectory(finalDirectory);
		setJobDirectory(jobDirectory);
		
		//Setup input files
		addInputFiles(inputFiles);
		
		//setup protect
		for (File inputFile: inputFiles) {
			addProtect(inputFile);
		}
		
		//setup delete
		for (File deleteFile: inputFiles) {
			addDelete(deleteFile);
		}

		//Setup output files
		addOutputFile(new File(finalDirectory,studyName + ".genomic.vcf.gz"));
		addOutputFile(new File(finalDirectory,studyName + ".genomic.vcf.gz.tbi"));
		
		String version = properties.get("VERSION");

		//Setup properties
		addLocalProperty("STUDY", studyName);
		
		if (inputFiles.size() == 2) {
			for (File f: inputFiles) {
				if (f.getName().endsWith(".vcf.gz")) {
					addLocalProperty("VCF_FILE", f.getName());
				} else {
					addLocalProperty("VCF_INDEX", f.getName());
				}
			}
			setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/variant_rename." + version + ".txt"));
		} else {
			StringBuilder gvcfString = new StringBuilder("");
			for (File input: inputFiles) {
				if (input.getName().endsWith("vcf.gz")) {
					gvcfString.append(" --variant " + input.getName());
					addDelete(input);
				}
			}
			addLocalProperty("VARIANT_LIST", gvcfString.toString());
			setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/variant_merge." + version + ".txt"));
		}
	
		addProtect(new File(jobDirectory,"cmd.txt"));
		
		
	}
	
}
