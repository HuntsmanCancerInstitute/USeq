package edu.utah.pysano.commands;

import java.io.File;
import java.util.HashMap;

import edu.utah.pysano.utils.Logger;

public class CommandRename extends Command {
	
	public CommandRename(File projectDirectory, Logger logFile, String sampleName, String sampleID, HashMap<String, String> properties, 
			String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		
		//Setup basic information
		setGlobalProperties(properties);
		File finalDirectory = new File(projectDirectory,"Alignments");
		File jobDirectory = new File(projectDirectory,"JOB_align_" + sampleID);
		setFinalDirectory(finalDirectory);
		setJobDirectory(jobDirectory);
		setSampleID(sampleID);
		setSampleName(sampleName);
		
		//Setup properties
		addLocalProperty("SAMPLE_ID", sampleID);
		addLocalProperty("NAME", sampleName);
		
		//Setup blocking
		addSkipFile(new File(finalDirectory,sampleID + ".raw.bam"));
		addSkipFile(new File(finalDirectory,sampleID + ".raw.bai"));
		
		//Setup template
		String version = properties.get("VERSION");
		setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/align_rename." + version + ".txt"));
		
		addProtect(new File(jobDirectory,"cmd.txt"));

	}
}
