package edu.utah.tomato.model.command;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.tomato.daemon.TFThread;
import edu.utah.tomato.daemon.TFThreadDaemon;
import edu.utah.tomato.model.TFCommand;
import edu.utah.tomato.model.TFFileObject;
import edu.utah.tomato.model.TFSampleInfo;
import edu.utah.tomato.util.TFConstants;
import edu.utah.tomato.util.TFLogger;
import edu.utah.tomato.model.TFMatchObject;





public class TFCommandExomeVariant extends TFCommand {
	
	private String splitType = "none";
	private String study = null;
	//private String[] hg19Chrom = {"chr1","chr2","chr3","chr4","chr5","chr6",
	//		"chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
	//		"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
	//		"chrX","chrY","chrM"};
	private String[] hg19Chrom = {"1","2","3","4","5","6","7","8","9","10",
			"11","12","13","14","15","16","17","18","19","20","21","22",
			"X","Y","MT"};
	private File targetFile;
	private boolean use1KGenomes;
	private boolean deleteReducedBams;
	private File finalGvcfDirectory = null;
	private File finalVcfDirectory = null;
	
	public TFCommandExomeVariant(ArrayList<File> templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, boolean deleteReducedBam, boolean isFull, boolean use1KGenomes,
			String study, File targetFile,
			HashMap<String,String> properties, String cluster) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress, isFull, properties, cluster);
		this.study = study;
		this.targetFile = targetFile;
		this.use1KGenomes = use1KGenomes;
		this.deleteReducedBams = deleteReducedBam;	
		this.finalDirectory = new File(this.rootDirectory,"Variants");
		this.finalGvcfDirectory = new File(this.finalDirectory,"GVCF");
		this.finalVcfDirectory = new File(this.finalDirectory,"VCF");
		this.jobDirectory = new File(this.finalDirectory,"Jobs");
	
	}
	
	@Override
	public ArrayList<TFSampleInfo> run(ArrayList<TFSampleInfo> sampleList) {
		TFThread.setFailCount(0);
		
		//Check for samples.
		if (sampleList.size() == 0) {
			sampleList = this.findPrereqsNew(sampleList);
		} else {
			sampleList = this.findPrereqsExisting(sampleList);
		}
		
	
		HashMap<String, Integer> commandMap = new HashMap<String,Integer>();
		commandMap.put("ugp_variant_haplotype",1);
		commandMap.put("ugp_variant_merge",2);
		commandMap.put("ugp_vqsr",3);
		commandMap.put("ugp_raw",3);
		commandMap.put("core_variant", 4);

		for (int i=0; i<this.templateFiles.size();i++) {
			String filename = this.templateFiles.get(i).getName();
			String[] fileNameParts = filename.split("\\.");
			String name = fileNameParts[0];
			if (!commandMap.containsKey(name)) {
				logFile.writeErrorMessage("[TFExomeVariantUgp] Command step " + name + " is not recognized",true);
				System.exit(1);
			}

			if (commandMap.get(name) == 1) {
				this.runHaplotypeCaller(sampleList,i);
			} else if (commandMap.get(name) == 2) {
				this.runGvcfMerge(sampleList,i);
			} else if (commandMap.get(name) == 3) {
				this.runGenotyping(sampleList, i);
			} else if (commandMap.get(name) == 4) {
				this.runUnified(sampleList, i);
			} else {
				logFile.writeErrorMessage("[TFExomeVariantUgp] Command step " + name + " is not assigned to a valid method.",true);
				System.exit(1);
			}
		}

		return sampleList;
	}
	
	@Override
	protected ArrayList<TFSampleInfo> validateSampleSet(ArrayList<TFSampleInfo> sampleList) {
		ArrayList<TFSampleInfo> validSamples = new ArrayList<TFSampleInfo>();
		boolean valid = true;
		for (TFSampleInfo tfsi: sampleList) {
			if (tfsi.finalFileExists(TFConstants.FILE_FINAL_BAM) && tfsi.finalFileExists(TFConstants.FILE_FINAL_BAI)) {
				//OK!
			} else {
				valid = false;
			}
			
			if (valid) {
				validSamples.add(tfsi);
			}
		}
		return validSamples;
	}
	
	@Override
	protected ArrayList<TFSampleInfo> findPrereqsExisting(ArrayList<TFSampleInfo> sampleList) {
		//Create patterns of interest
		ArrayList<TFMatchObject> dependantPatterns = new ArrayList<TFMatchObject>();
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAM,Pattern.compile("(.+?)\\.final\\.bam$"),TFConstants.PREFIX_SAMPLENAME));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAI,Pattern.compile("(.+?)\\.final\\.bai$"),TFConstants.PREFIX_SAMPLENAME));
		
		File rootAlignments = new File(this.rootDirectory,"Alignments");
		File processedAlignments = new File(rootAlignments,"ProcessedAlignments");
		if (processedAlignments.exists()) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsExisting(sampleList,processedAlignments, dependantPatterns);
			if (foundSampleList.size() > 0) {
				if (foundSampleList.size() == sampleList.size()) {
					sampleList = foundSampleList;
					logFile.writeInfoMessage(String.format("[TFExomeVariant] Found %d of %d samples in the ProcessedAlignments directory",foundSampleList.size(),sampleList.size()));
				} else {
					logFile.writeErrorMessage(String.format("[TFExomeVariant] Found fewer than expected samples (%d of %d) in the ProcessedAlignments directory, exiting",foundSampleList.size(),sampleList.size()), false);
					System.exit(1);
				} 
			} else {
				logFile.writeErrorMessage("[TFExomeVariant] Did not find any potential samples in the ProcessedAlignments directory, exiting",false);
				System.exit(1);
			}
		} else {
			logFile.writeErrorMessage("[TFExomeVariant] Could not find the ProcessedAlignments directory, exiting",true);
			System.exit(1);
		}
		
		return sampleList;
	}

	@Override
	protected ArrayList<TFSampleInfo> findPrereqsNew(ArrayList<TFSampleInfo> sampleList) {
		ArrayList<TFMatchObject> masterPatterns = new ArrayList<TFMatchObject>();
		
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAM,Pattern.compile("(.+?)\\.final\\.bam$"),TFConstants.PREFIX_SAMPLENAME));
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAI,Pattern.compile("(.+?)\\.final\\.bai$"),TFConstants.PREFIX_SAMPLENAME));
		
		File rootAlignments = new File(this.rootDirectory,"Alignments");
		File processedAlignments = new File(rootAlignments,"ProcessedAlignments");
		if (processedAlignments.exists()) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsNew(processedAlignments, masterPatterns, new ArrayList<TFMatchObject>());
			
			if (foundSampleList.size() > 0) {
				sampleList = foundSampleList;
				logFile.writeInfoMessage(String.format("[TFExomeVariant] Found %d samples in the ProcessedAlignments directory.",foundSampleList.size()));
			} else {
				logFile.writeInfoMessage("[TFExomeVariant] Did not find any potential samples in the ProcessedAlignments directory, checking run directory");
			}
		}
		
		if (sampleList.size() == 0) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsNew(this.rootDirectory, masterPatterns, new ArrayList<TFMatchObject>());
			if (foundSampleList.size() > 0) {
				sampleList = foundSampleList;
				logFile.writeInfoMessage(String.format("[TFExomeVariant] Found %d samples in the run directory directory.",foundSampleList.size()));
			} else {
				logFile.writeErrorMessage("[TFExomeVariant] Did not find any potential samples in the run directory, exiting",false);
				System.exit(1);
			}
		}
		
		return sampleList;
	}
	
	private void runUnified(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		File workingDir = new File(this.rootDirectory,"JOB_" + this.study + "_variant");
		
		TFFileObject tfoVcf = new TFFileObject(this.study + ".raw.vcf.gz",this.finalDirectory,workingDir);
		TFFileObject tfoVcfIdx = new TFFileObject(this.study + ".raw.vcf.gz.tbi",this.finalDirectory,workingDir);
		
		if (tfoVcf.doesFinalExist() && tfoVcfIdx.doesFinalExist())  {
			return;
		}
			
		
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<File> runDirectoryList = new HashSet<File>();
		ArrayList<File> rawVcfList = new ArrayList<File>();
		
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		
		if (!tfoVcf.doesWorkingExist() || !tfoVcfIdx.doesWorkingExist()) {
			this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(templateIdx).getName(), 1, this.jobs);
			this.daemon.start();
			
			if (!workingDir.exists()) {
				workingDir.mkdir();
			}
			runDirectoryList.add(workingDir);
			
		    int counter  = 1;
			ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
			ArrayList<String> bamList = new ArrayList<String>(); //bam names 
			
			for (String sampleName: samples.keySet()) {
				TFSampleInfo tfsi = samples.get(sampleName).get(0);
				
				TFFileObject tfoFinalBam = tfsi.getFileObject(TFConstants.FILE_FINAL_BAM);
				TFFileObject tfoFinalBai = tfsi.getFileObject(TFConstants.FILE_FINAL_BAI);
				
				File fileReduceBam = tfoFinalBam.createDestForFileObject(workingDir);
				File fileReduceBai = tfoFinalBai.createDestForFileObject(workingDir);
				
				this.createLink(tfoFinalBam.getFinalPath(), fileReduceBam);
				this.createLink(tfoFinalBai.getFinalPath(), fileReduceBai);
				
				//Add bam/bai to preserve
				protectList.add(fileReduceBam);
				protectList.add(fileReduceBai);
				deleteList.add(fileReduceBam);
				deleteList.add(fileReduceBai);
				
				//Add bam name to list
				bamList.add(tfoFinalBam.getFileName());
			}
			
			//Create bam String
			String bamString = "";
			for (String b: bamList) {
				bamString += " -I " + b;
			}
			
			//Create replacement tokens
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("STUDY", this.study);
			replacements.put("BAM_LIST", bamString);
			
			if (targetFile == null) {
				replacements.put("TARGETS","");
			} else {
				File localTarget = new File(workingDir,targetFile.getName());
				replacements.put("TARGETS","-L " + targetFile.getName());
				this.cpFile(targetFile, localTarget);
				protectList.add(localTarget);
			}
			replacements.putAll(this.properties);
			
			//Create cmd.txt file
			File cmdFile = new File(workingDir,"cmd.txt");
			protectList.add(cmdFile);
			this.createCmd(replacements, cmdFile, templateIdx);
			
			TFThread thread = new TFThread(workingDir,this.failmax, counter, this.heartbeat, protectList, this.logFile);
			this.daemon.addJob(thread);
			
			//Wait for command to finish
			try {
				this.daemon.join();
				if (this.daemon.getFailed()) {
					System.exit(1);
				}
				Thread.sleep(5000);
			} catch (InterruptedException ie) {
				logFile.writeErrorMessage("[TFExomeVariant] Daemon interrupted",true);
				System.exit(1);
			}
			
			File fullVcf = new File(workingDir,this.study + ".raw.vcf.gz");
			File fullVcfIdx = new File(workingDir,this.study + ".raw.vcf.gz.tbi");
			
			rawVcfList.add(fullVcf);
			deleteList.add(fullVcf);
			deleteList.add(fullVcfIdx);
		} 
		
		this.finalDirectory.mkdirs();
		this.jobDirectory.mkdir();
				
		logFile.writeInfoMessage("[TFExomeVariant] Moving vcf files");

		this.moveFile(tfoVcf.getWorkingPath(),tfoVcf.getFinalPath());
		this.moveFile(tfoVcfIdx.getWorkingPath(),tfoVcfIdx.getFinalPath());
		
		this.cleanup(runDirectoryList, deleteList);		
		
	}
	
	private void runHaplotypeCaller(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		/**************************************************************************
		 * Run haplotype caller
		 ************************************************************************/
		logFile.writeInfoMessage("[TFExomeVariant] Running Haplotype Caller");

		
		//Containers for bam-splitting processing 
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<String> samplesToRun = new HashSet<String>();
		HashSet<String> samplesToPostProcess = new HashSet<String>();
		
		/**************************************************************************
		 * Group by sample
		 ************************************************************************/
		
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		
		/**************************************************************************
		 * Determine sample stage
		 ************************************************************************/
		
		for (String sampleName: samples.keySet()) {
			File workingDir = new File(this.rootDirectory,"JOB_" + sampleName + "_haplotype"); //Create working directory
			
			for (TFSampleInfo si: samples.get(sampleName)) {
				
				
				//Create output objects for each sample
				TFFileObject tfoGvcf = new TFFileObject(si.getSampleName() + ".genomic.vcf.gz",this.finalGvcfDirectory,workingDir);
				TFFileObject tfoGvcfIdx = new TFFileObject(si.getSampleName() + ".genomic.vcf.gz.tbi",this.finalGvcfDirectory,workingDir);
				TFFileObject tfoStudyGvcf = new TFFileObject(this.study + ".genomic.vcf.gz",this.finalGvcfDirectory,workingDir);
				TFFileObject tfoStudyGvcfIdx = new TFFileObject(this.study + ".genomic.vcf.gz.tbi",this.finalGvcfDirectory,workingDir);
			
				si.setFileObject(TFConstants.FILE_GVCF, tfoGvcf);
				si.setFileObject(TFConstants.FILE_GVCF_IDX, tfoGvcfIdx);
				si.setFileObject(TFConstants.FILE_GVCF_STUDY, tfoStudyGvcf);
				si.setFileObject(TFConstants.FILE_GVCF_STUDY_IDX, tfoStudyGvcfIdx);

				
				//Check for the existance of output files
				if ((!tfoGvcf.doesFinalExist() || !tfoGvcfIdx.doesFinalExist()) && (!tfoStudyGvcf.doesFinalExist() || !tfoStudyGvcfIdx.doesFinalExist())) {
					if (!si.getFileObject(TFConstants.FILE_FINAL_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_FINAL_BAI).doesFinalExist()) {
						logFile.writeErrorMessage("[TFExomeVariant] Final alignment files not found by TF", true);
						System.exit(1);
					}
					
					if (!tfoGvcf.doesWorkingExist() || !tfoGvcf.doesWorkingExist()) {
						samplesToRun.add(sampleName);
					} else {
						samplesToPostProcess.add(sampleName);
					}
				}
			}
		}
		
		
		/**************************************************************************
		 * Process samples
		 ************************************************************************/
		
		if (samplesToRun.size() > 0) {
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(templateIdx).getName(),samplesToRun.size(),this.jobs);
			this.daemon.start();
			
			int counter = 1;
			
			for(String sampleName: samplesToRun) {
				samplesToPostProcess.add(sampleName); //Add samples to post-processing list
				
				ArrayList<File> protectList = new ArrayList<File>();
				
				//Get representative sample and create input/output objects
				TFSampleInfo repSI = samples.get(sampleName).get(0);
				TFFileObject tfoGvcf = repSI.getFileObject(TFConstants.FILE_GVCF);
				TFFileObject tfoFinalBam = repSI.getFileObject(TFConstants.FILE_FINAL_BAM);
				TFFileObject tfoFinalBai = repSI.getFileObject(TFConstants.FILE_FINAL_BAI);
				
				//Create run directory
				File runDirectory = tfoGvcf.getWorkingDirectory();
				if (runDirectory.exists()) {
					this.deleteFolder(runDirectory);
				}
				runDirectory.mkdir();
				
				//Create working files
				File fileFinalBam = tfoFinalBam.createDestForFileObject(runDirectory);
				File fileFinalBai = tfoFinalBai.createDestForFileObject(runDirectory);
				
				//Create links
				this.createLink(tfoFinalBam.getFinalPath(), fileFinalBam);
				this.createLink(tfoFinalBai.getFinalPath(), fileFinalBai);
				
				//create properties
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("SAMPLE", sampleName);
				if (targetFile == null) {
					replacements.put("TARGETS","");
				} else {
					File localTarget = new File(runDirectory,targetFile.getName());
					replacements.put("TARGETS","-L " + targetFile.getName());
					this.cpFile(targetFile, localTarget);
					protectList.add(localTarget);
				}
				replacements.putAll(this.properties);
				
				//create files
				File cmdFile = new File(runDirectory,"cmd.txt");
				this.createCmd(replacements,cmdFile,templateIdx);
				
				deleteList.add(fileFinalBam);
				deleteList.add(fileFinalBai);
				protectList.add(fileFinalBam);
				protectList.add(fileFinalBai);
				protectList.add(cmdFile);
	
				//Run job
				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
				this.daemon.addJob(thread);
				counter ++;
			}
				
			//Wait for command to finish
			try {
				this.daemon.join();
				Thread.sleep(5000);
				if (this.daemon.getFailed()) {
					System.exit(1);
				}
			} catch (InterruptedException ie) {
				logFile.writeErrorMessage("[TFExomeAlignUgp] Daemon interrupted",true);
				System.exit(1);
			}
		}
		
		/**************************************************************************
		 * Post-Process Samples
		 ************************************************************************/
		
		if (samplesToPostProcess.size() > 0) {
			//Create final directories
			this.finalDirectory.mkdir();
			this.finalGvcfDirectory.mkdir();
			this.jobDirectory.mkdir();
			
			HashSet<File> runDirectoryList = new HashSet<File>();
			
			for (String sampleName: samplesToPostProcess) {
				 //Get representative sample
				TFSampleInfo repSI = samples.get(sampleName).get(0);
				
				//Get final bam file objects
			    TFFileObject tfoGvcf = repSI.getFileObject(TFConstants.FILE_GVCF);
			    TFFileObject tfoGvcfIdx = repSI.getFileObject(TFConstants.FILE_GVCF_IDX);
			    
			    File runDirectory = tfoGvcf.getWorkingDirectory();
			    runDirectoryList.add(runDirectory);
			    
			    File fileGvcf = tfoGvcf.createDestForFileObject(runDirectory);
			    File fileGvcfIdx = tfoGvcfIdx.createDestForFileObject(runDirectory);
			   
			    this.moveFile(fileGvcf, tfoGvcf.getFinalPath());
			    this.moveFile(fileGvcfIdx, tfoGvcfIdx.getFinalPath());
			}
			
			this.cleanup(runDirectoryList, deleteList);
		}
	}
	
	private void runGvcfMerge(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		logFile.writeInfoMessage("[TFExomeVariant] Running GVCF Merge");
		File workingDir = new File(this.rootDirectory,"JOB_" + this.study + "_variantmerge");
		
		
		TFFileObject tfoGvcfStudy = new TFFileObject(this.study + ".genomic.vcf.gz",this.finalGvcfDirectory,workingDir);
		TFFileObject tfoGvcfStudyIdx = new TFFileObject(this.study + ".genomic.vcf.gz.tbi",this.finalGvcfDirectory,workingDir);
		
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<File> runDirectoryList = new HashSet<File>();
		
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		if (tfoGvcfStudy.doesFinalExist() && tfoGvcfStudyIdx.doesFinalExist())  {
			return;
		}
		
		for (TFSampleInfo si: sampleList) {
			if (!si.getFileObject(TFConstants.FILE_GVCF).doesFinalExist() || !si.getFileObject(TFConstants.FILE_GVCF_IDX).doesFinalExist()) {
				logFile.writeErrorMessage("[TFExomeVariant] Sample GVCF files can't be found by the GVCF merge step", true);
				System.exit(1);
			}
			
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		
		
		if (samples.size() == 1 ) {
			if (workingDir.exists()) {
				this.deleteFolder(workingDir);
			}
			workingDir.mkdir();
			//Only iterate over one
			for (String sampleName: samples.keySet()) {
				TFSampleInfo tfsi = samples.get(sampleName).get(0);
				TFFileObject sampleGvcf = tfsi.getFileObject(TFConstants.FILE_GVCF);
				TFFileObject sampleGvcfIdx = tfsi.getFileObject(TFConstants.FILE_GVCF_IDX);
				TFFileObject studyGvcf = tfsi.getFileObject(TFConstants.FILE_GVCF_STUDY);
				TFFileObject studyGvcfIdx = tfsi.getFileObject(TFConstants.FILE_GVCF_STUDY_IDX);
				
				this.moveFile(sampleGvcf.getFinalPath(), studyGvcf.getWorkingPath());
				this.moveFile(sampleGvcfIdx.getFinalPath(), studyGvcfIdx.getWorkingPath());
			}
			
			workingDir.delete();
		}
		
		
		if (!tfoGvcfStudy.doesWorkingExist() || !tfoGvcfStudyIdx.doesWorkingExist()) {
			this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(templateIdx).getName(), 1, this.jobs);
			this.daemon.start();
			
			runDirectoryList.add(workingDir);

			if (workingDir.exists()) {
				this.deleteFolder(workingDir);
			}
			workingDir.mkdir();
			
		    int counter  = 1;
			ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
			ArrayList<String> gvcfList = new ArrayList<String>(); //bam names 
			
			for (String sampleName: samples.keySet()) {
				TFSampleInfo tfsi = samples.get(sampleName).get(0);
				
				TFFileObject tfoGvcf = tfsi.getFileObject(TFConstants.FILE_GVCF);
				TFFileObject tfoGvcfIdx = tfsi.getFileObject(TFConstants.FILE_GVCF_IDX);
				
				
				
				File fileGvcf = tfoGvcf.createDestForFileObject(workingDir);
				File fileGvcfIdx = tfoGvcfIdx.createDestForFileObject(workingDir);
				
				this.createLink(tfoGvcf.getFinalPath(), fileGvcf);
				this.createLink(tfoGvcfIdx.getFinalPath(), fileGvcfIdx);
				
				//Add bam/bai to preserve
				protectList.add(fileGvcf);
				protectList.add(fileGvcfIdx);
				deleteList.add(fileGvcf);
				deleteList.add(fileGvcfIdx);
				
				//Add bam name to list
				gvcfList.add(tfoGvcf.getFileName());
			}
			
			//Create bam String
			String gvcfString = "";
			for (String g: gvcfList) {
				gvcfString += " --variant " + g;
			}
			
			//Create replacement tokens
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("STUDY", this.study);
			replacements.put("VARIANT_LIST", gvcfString);
			
			replacements.putAll(this.properties);
			
			//Create cmd.txt file
			File cmdFile = new File(workingDir,"cmd.txt");
			protectList.add(cmdFile);
			this.createCmd(replacements, cmdFile, templateIdx);
			
			TFThread thread = new TFThread(workingDir,this.failmax, counter, this.heartbeat, protectList, this.logFile);
			this.daemon.addJob(thread);
			
			//Wait for command to finish
			try {
				this.daemon.join();
				if (this.daemon.getFailed()) {
					System.exit(1);
				}
				Thread.sleep(5000);
			} catch (InterruptedException ie) {
				logFile.writeErrorMessage("[TFExomeVariant] Daemon interrupted",true);
				System.exit(1);
			}
		}
			
		//Move the results
		this.finalDirectory.mkdirs();
		this.jobDirectory.mkdir();
		this.finalGvcfDirectory.mkdir();
		
		logFile.writeInfoMessage("[TFExomeVariant] Moving study gvcf file");
		
		this.moveFile(tfoGvcfStudy.getWorkingPath(),tfoGvcfStudy.getFinalPath());
		this.moveFile(tfoGvcfStudyIdx.getWorkingPath(),tfoGvcfStudyIdx.getFinalPath());
		
		this.cleanup(runDirectoryList, deleteList);	
	}
	
	
	private void runGenotyping(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		logFile.writeInfoMessage("[TFExomeVariant] Running Genotyper");
		File workingDir = new File(this.rootDirectory,"JOB_" + this.study + "_genotype");
		
		boolean isFinished = false;
		boolean isCleanup = false;
		TFFileObject tfoRawVcf = new TFFileObject(this.study + ".raw.vcf.gz",this.finalVcfDirectory,workingDir);
		TFFileObject tfoRawVcfIdx = new TFFileObject(this.study + ".raw.vcf.gz.tbi",this.finalVcfDirectory,workingDir);
		TFFileObject tfoFilterVcf = new TFFileObject(this.study + ".filterFieldSetAll.vcf.gz",this.finalVcfDirectory,workingDir);
		TFFileObject tfoFilterVcfIdx = new TFFileObject(this.study + ".filterFieldSetAll.vcf.gz.tbi",this.finalVcfDirectory,workingDir);
		TFFileObject tfoPassingVcf = new TFFileObject(this.study + ".filterFieldSetPassing.vcf.gz",this.finalVcfDirectory,workingDir);
		TFFileObject tfoPassingVcfIdx = new TFFileObject(this.study + ".filterFieldSetPassing.vcf.gz.tbi",this.finalVcfDirectory,workingDir);
		
		if (tfoRawVcf.doesFinalExist() && tfoRawVcfIdx.doesFinalExist() 
				&& tfoFilterVcf.doesFinalExist() && tfoFilterVcfIdx.doesFinalExist() 
				&& tfoPassingVcf.doesFinalExist() && tfoPassingVcfIdx.doesFinalExist() ) {
			isFinished = true;
		} else if (tfoRawVcf.doesWorkingExist() && tfoRawVcfIdx.doesWorkingExist() 
				&& tfoFilterVcf.doesWorkingExist() && tfoFilterVcfIdx.doesWorkingExist() 
				&& tfoPassingVcf.doesWorkingExist() && tfoPassingVcfIdx.doesWorkingExist()) {
			isCleanup = true;
		} else {
			for (TFSampleInfo tfsi: sampleList) {
				tfsi.setFileObject(TFConstants.FILE_VCF_RAW, tfoRawVcf);
				tfsi.setFileObject(TFConstants.FILE_VCF_RAW_IDX,tfoRawVcfIdx);
				tfsi.setFileObject(TFConstants.FILE_VCF_FILTER,tfoFilterVcf);
				tfsi.setFileObject(TFConstants.FILE_VCF_FILTER_IDX, tfoFilterVcfIdx);
				tfsi.setFileObject(TFConstants.FILE_VCF_PASSING,tfoPassingVcf);
				tfsi.setFileObject(TFConstants.FILE_VCF_PASSING_IDX, tfoPassingVcfIdx);
			}
		}
		
		TFSampleInfo repSample = sampleList.get(0);
		
		TFFileObject tfoGvcf = repSample.getFileObject(TFConstants.FILE_GVCF_STUDY);
		TFFileObject tfoGvcfIdx = repSample.getFileObject(TFConstants.FILE_GVCF_STUDY_IDX);
		
		if (!tfoGvcf.doesFinalExist() || !tfoGvcfIdx.doesFinalExist()) {
			logFile.writeErrorMessage("[TFExomeVariant] Study GVCF file can't be found by the genotyping step", true);
			System.exit(1);
		}
		
		
		if (!isFinished) {
			HashSet<File> deleteList = new HashSet<File>();
			HashSet<File> runDirectoryList = new HashSet<File>();
			if (!isCleanup) {
				if (workingDir.exists()) {
					this.deleteFolder(workingDir);
				}
				workingDir.mkdir();
				
				int counter  = 1;
				ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
				
				
				File fileGvcf = tfoGvcf.createDestForFileObject(workingDir);
				File fileGvcfIdx = tfoGvcfIdx.createDestForFileObject(workingDir);
				
				this.createLink(tfoGvcf.getFinalPath(), fileGvcf);
				this.createLink(tfoGvcfIdx.getFinalPath(), fileGvcfIdx);
				
				protectList.add(fileGvcf);
				protectList.add(fileGvcfIdx);
				deleteList.add(fileGvcf);
				deleteList.add(fileGvcfIdx);
				
				//Add merged to runDirectoryList
				runDirectoryList.add(workingDir);
				
				//Initialize daemon
				this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(templateIdx).getName(), 1, this.jobs);
				this.daemon.start();
				
				//Create replacement tokens
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("STUDY", this.study);
				
				if (this.use1KGenomes) {
					logFile.writeInfoMessage("[TFExomeVariant] Adding 1K Genome background files.");
					ArrayList<TFSampleInfo> background = this.get1KGenomes();
					String backgroundString = "";
					for (TFSampleInfo tfsi: background) {
						TFFileObject backVcf = tfsi.getFileObject(TFConstants.FILE_GVCF);
						TFFileObject backVcfIdx = tfsi.getFileObject(TFConstants.FILE_GVCF_IDX);
						File backFileGvcf = backVcf.createDestForFileObject(workingDir);
						File backFileGvcfIdx = backVcfIdx.createDestForFileObject(workingDir);
						this.createLink(backVcf.getFinalPath(), backFileGvcf);
						this.createLink(backVcfIdx.getFinalPath(), backFileGvcfIdx);
						
						protectList.add(backFileGvcf);
						protectList.add(backFileGvcfIdx);
						deleteList.add(backFileGvcf);
						deleteList.add(backFileGvcfIdx);
						
						backgroundString += " --variant " + backFileGvcf.getName() + " ";
					}
					replacements.put("BACKGROUND",backgroundString);
					
				} else {
					replacements.put("BACKGROUND","");
				}
				
				replacements.putAll(this.properties);
				
				//Create cmd.txt file
				File cmdFile = new File(workingDir,"cmd.txt");
				protectList.add(cmdFile);
				this.createCmd(replacements, cmdFile, templateIdx);
				
				TFThread thread = new TFThread(workingDir,this.failmax, counter, this.heartbeat, protectList, this.logFile);
				this.daemon.addJob(thread);
				
				//Wait for command to finish
				try {
					this.daemon.join();
					if (this.daemon.getFailed()) {
						System.exit(1);
					}
					Thread.sleep(5000);
				} catch (InterruptedException ie) {
					logFile.writeErrorMessage("[TFExomeVariant] Daemon interrupted",true);
					System.exit(1);
				}
			}
			
			//Move output files
			this.finalVcfDirectory.mkdir();
			logFile.writeInfoMessage("[TFExomeVariant] Moving vcf files");
			this.moveFile(tfoRawVcf.getWorkingPath(), tfoRawVcf.getFinalPath());
			this.moveFile(tfoRawVcfIdx.getWorkingPath(), tfoRawVcfIdx.getFinalPath());
			this.moveFile(tfoFilterVcf.getWorkingPath(), tfoFilterVcf.getFinalPath());
			this.moveFile(tfoFilterVcfIdx.getWorkingPath(), tfoFilterVcfIdx.getFinalPath());
			this.moveFile(tfoPassingVcf.getWorkingPath(), tfoPassingVcf.getFinalPath());
			this.moveFile(tfoPassingVcfIdx.getWorkingPath(), tfoPassingVcfIdx.getFinalPath());
			
			this.cleanup(runDirectoryList, deleteList);
		}
	}

	public ArrayList<TFSampleInfo> get1KGenomes() {
		File oneKDir = new File(this.properties.get("BACKGROUND_PATH"));
		if (oneKDir.canRead()) {
			File[] gvcfs = oneKDir.listFiles();
			ArrayList<TFSampleInfo> oneKTFList = new ArrayList<TFSampleInfo>();
			Pattern gvcfPattern = Pattern.compile("(.+?).vcf.gz$");
			for (File gvcf: gvcfs) {
				String name = gvcf.getName();
				Matcher m = gvcfPattern.matcher(name);
				if (m.matches()) {
					String sampleName = m.group(1);
					TFSampleInfo nsi = new TFSampleInfo(sampleName,sampleName,"",this.logFile);
					TFFileObject tfoGvcf = new TFFileObject(sampleName + ".vcf.gz",oneKDir,oneKDir);
					TFFileObject tfoGvcfIdx = new TFFileObject(sampleName + ".vcf.gz.tbi",oneKDir,oneKDir);
			
					nsi.setFileObject(TFConstants.FILE_GVCF, tfoGvcf);
					nsi.setFileObject(TFConstants.FILE_GVCF_IDX, tfoGvcfIdx);
				
					oneKTFList.add(nsi);
				}
				
			}
			return oneKTFList;
		} else {
			return null;
		}
	}
}
