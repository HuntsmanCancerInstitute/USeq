package edu.utah.pysano.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import edu.utah.pysano.utils.Logger;

public class CommandGenotype extends Command {
	
	//Constructor for combined command (don't check for input files pre-run)
	public CommandGenotype(File projectDirectory, HashMap<String, String> properties, Logger logFile, File targetFile, String studyName, boolean isVQSR,
			boolean use1KGenomes, int totalSamples, String cluster, String email, Integer wallTime, boolean suppress) {
		super(projectDirectory, cluster,email,wallTime,suppress, logFile);
		
		//Setup basic information
		setGlobalProperties(properties);
		File finalDirectory = new File(projectDirectory,"Variants/VCF");
		File jobDirectory = new File(projectDirectory,"JOB_genotype_" + studyName);
		setFinalDirectory(finalDirectory);
		setJobDirectory(jobDirectory);
		this.logFile = logFile;
		
		//setup external
		if (targetFile == null) {
			targetFile = new File(properties.get("TARGET_DEFAULT"));
		}
		addExternalFile(targetFile);
		
		//setup delete
		addDelete(targetFile);
		
		
		
		//Setup output files
		addOutputFile(new File(finalDirectory,studyName + ".raw.vcf.gz"));
		addOutputFile(new File(finalDirectory,studyName + ".raw.vcf.gz.tbi"));
		addOutputFile(new File(finalDirectory,studyName + ".filterFieldSetAll.vcf.gz"));
		addOutputFile(new File(finalDirectory,studyName + ".filterFieldSetAll.vcf.gz.tbi"));
		addOutputFile(new File(finalDirectory,studyName + ".filterFieldSetPassing.vcf.gz"));
		addOutputFile(new File(finalDirectory,studyName + ".filterFieldSetPassing.vcf.gz.tbi"));
		
		
		//Setup properties
		addLocalProperty("STUDY", studyName);
		
		if (use1KGenomes) {
			logFile.writeInfoMessage("[CommandGenotype] Adding 1K Genome background files.");
			ArrayList<File> backgroundFiles = this.get1KGenomes(properties);
			StringBuilder backgroundString = new StringBuilder("");
			
			for (File f: backgroundFiles) {
				if (f.getAbsolutePath().endsWith("vcf.gz")) {
					backgroundString.append(" --variant " + f.getName() + " ");
					
					totalSamples++;
				}
				addExternalFile(f);
				addProtect(f);
				addDelete(f);
			}
		
			addLocalProperty("BACKGROUND",backgroundString.toString());
		} else {
			addLocalProperty("BACKGROUND","");
		}
		
		addLocalProperty("TARGETS", "-L " + targetFile.getName());
		if (totalSamples > 20) {
			addLocalProperty("THREADS","-nt \\$NCPU");
			addLocalProperty("INBREED","-an InbreedingCoeff");
		} else {
			addLocalProperty("THREADS","");
			addLocalProperty("INBREED","");
		}

		//Setup template
		String version = properties.get("VERSION");
		if (isVQSR) {
			setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/variant_vqsr." + version + ".txt"));
		} else {
			setTemplateFile(new File(properties.get("TEMPLATE_PATH"),"templates/variant_raw." + version + ".txt"));
		}
		

		File commandFile = new File(jobDirectory,"cmd.txt");
		
		//setup protect
		addProtect(commandFile);
		addProtect(targetFile);
		
	}
	
	private ArrayList<File> get1KGenomes(HashMap<String,String> properties) {
		File oneKDir = new File(properties.get("BACKGROUND_PATH"));
		if (oneKDir.canRead()) {
			File[] gvcfs = oneKDir.listFiles();
			Pattern gvcfPattern1 = Pattern.compile("(.+?).vcf.gz");
			Pattern gvcfPattern2 = Pattern.compile("(.+?).vcf.gz.tbi");
			HashMap<String,File> matchingFiles1 = matchFiles(gvcfPattern1,gvcfs);
			HashMap<String,File> matchingFiles2 = matchFiles(gvcfPattern2,gvcfs);
			
			//Find samples with both files and initialize command
			ArrayList<File> inputFileList = new ArrayList<File>();
			for (String key: matchingFiles1.keySet()) {
				if (matchingFiles2.containsKey(key)) { 
					//Required files are found, determine names				
					inputFileList.add(matchingFiles1.get(key));
					inputFileList.add(matchingFiles2.get(key));
				}
			}
			
			return inputFileList;
		} else {
			return null;
		}
	}
}
