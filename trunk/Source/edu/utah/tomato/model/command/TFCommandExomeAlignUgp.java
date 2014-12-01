package edu.utah.tomato.model.command;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import edu.utah.tomato.daemon.TFThread;
import edu.utah.tomato.daemon.TFThreadDaemon;
import edu.utah.tomato.model.TFCommand;
import edu.utah.tomato.model.TFFileObject;
import edu.utah.tomato.model.TFMatchObject;
import edu.utah.tomato.model.TFSampleInfo;
import edu.utah.tomato.util.TFConstants;
import edu.utah.tomato.util.TFLogger;


public class TFCommandExomeAlignUgp extends TFCommand {
	private boolean validateFastq = false;

	public TFCommandExomeAlignUgp(ArrayList<File> templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, boolean isFull, boolean validateFastq, 
			HashMap<String,String> properties, String cluster) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress, isFull, properties,cluster);
		// TODO Auto-generated constructor stub
		File subDirectory = new File(this.rootDirectory,"Alignments");
		this.finalDirectory = new File(subDirectory,"ProcessedAlignments");
		this.jobDirectory = new File(subDirectory,"Jobs");
		this.validateFastq = validateFastq;
	}
	
	@Override
	public ArrayList<TFSampleInfo> run(ArrayList<TFSampleInfo> sampleList) {
		TFThread.setFailCount(0);
		
		if (sampleList.size() > 0) {
			logFile.writeErrorMessage("[TFExomeAlignUgp] Command modules upstream of this one aren't currently supported.",true);
			System.exit(1);
		} else {
			sampleList = this.findPrereqsNew(sampleList);
		}
		
		if (this.validateFastq) {
			for (TFSampleInfo fi: sampleList) {
				validateFileSet(fi);
			}
		}
		
		HashMap<String, Integer> commandMap = new HashMap<String,Integer>();
		commandMap.put("core_nov",0);
		commandMap.put("core_bwa",0);
		commandMap.put("ugp_nov",0);
		commandMap.put("ugp_bwa",0);
		commandMap.put("common_merge", 1);
		commandMap.put("core_dedup", 2);
		commandMap.put("ugp_realign", 2);
		
		for (int i=0; i<this.templateFiles.size();i++) {
			String filename = this.templateFiles.get(i).getName();
			String[] fileNameParts = filename.split("\\.");
			String name = fileNameParts[0];
			if (!commandMap.containsKey(name)) {
				logFile.writeErrorMessage("[TFExomeAlignUgp] Command step " + name + " is not recognized",true);
				System.exit(1);
			}
			
			if (commandMap.get(name) == 0) {
				this.alignment(sampleList,i);
			} else if (commandMap.get(name) == 1) {
				this.mergeSampleLanes(sampleList,i);
			} else if (commandMap.get(name) == 2) {
				this.postProcessSample(sampleList, i);
			} else {
				logFile.writeErrorMessage("[TFExomeAlignUgp] Command step " + name + " is not assigned to a valid method.",true);
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
			if (tfsi.finalFileExists(TFConstants.FILE_FASTQ1) && tfsi.finalFileExists(TFConstants.FILE_FASTQ2)) {
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
		//There are currently no command modules upstream of this one, so leaving this method blank for now.
		return null;
	}
	
	@Override
	protected ArrayList<TFSampleInfo> findPrereqsNew(ArrayList<TFSampleInfo> sampleList) {
		
		//Create patterns of interest
		ArrayList<TFMatchObject> dependantPatterns = new ArrayList<TFMatchObject>(); //No dependent patterns.
		ArrayList<TFMatchObject> masterPatterns = new ArrayList<TFMatchObject>();
		
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_FASTQ1, Pattern.compile("(.+?)(_1|_R1|_001)\\..+\\.gz$"), TFConstants.PREFIX_SAMPLEID));
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_FASTQ2, Pattern.compile("(.+?)(_2|_R2|_002)\\..+\\.gz$"), TFConstants.PREFIX_SAMPLEID));
		
		ArrayList<TFSampleInfo> foundSampleList = this.findPatternsNew(this.rootDirectory, masterPatterns, dependantPatterns);
		
		if (foundSampleList.size() > 0) {
			sampleList = foundSampleList;
			logFile.writeInfoMessage(String.format("[TFExomeAlignUgp] Found %d samples in the run directory.",foundSampleList.size()));
		} else {
			logFile.writeErrorMessage("[TFExomeAlignUgp] Did not find any matching samples in the run directory.",false);
			System.exit(1);
		}
		
		return sampleList;

	}
	
	
	public void alignment(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		/******************************************************************************************
		 *  Standard Alignment
		 *****************************************************************************************/
		logFile.writeInfoMessage("[TFExomeAlignUgp] Starting sample alignment");
		
		//Create project-specific storage
		HashSet<File> deleteList = new HashSet<File>(); //list of files to delete at the end of the run
		HashSet<TFSampleInfo> samplesToRun = new HashSet<TFSampleInfo>();
		HashSet<TFSampleInfo> samplesToPostProcess = new HashSet<TFSampleInfo>();
		
		
		
		for (TFSampleInfo si: sampleList) {
			//Create working directory
			File workingDir = new File(this.rootDirectory,"JOB_" + si.getSampleID() + "_align");
			
			//Create raw bam objects
			TFFileObject tfoBam = new TFFileObject(si.getSampleID() + ".raw.bam",this.finalDirectory,workingDir);
			TFFileObject tfoBai = new TFFileObject(si.getSampleID() + ".raw.bai",this.finalDirectory,workingDir);
			TFFileObject tfoSampleBam = new TFFileObject(si.getSampleName() + ".raw.bam",this.finalDirectory,workingDir);
			TFFileObject tfoSampleBai = new TFFileObject(si.getSampleName() + ".raw.bai",this.finalDirectory,workingDir);
			
			si.setFileObject(TFConstants.FILE_BAM,tfoBam);
			si.setFileObject(TFConstants.FILE_BAI,tfoBai);
			
			
			//Check to see if the output files already exist
			if ((tfoBam.doesFinalExist() && tfoBai.doesFinalExist()) || (tfoSampleBam.doesFinalExist() && tfoSampleBai.doesFinalExist())) {
				//OK
			} else{
				//Make sure files are initialized
				if (!si.getFileObject(TFConstants.FILE_FASTQ1).doesFinalExist() || !si.getFileObject(TFConstants.FILE_FASTQ2).doesFinalExist()) {
					logFile.writeErrorMessage("[TFExomeAlignUgp] Your fastq files weren't found by TF",true);
					System.exit(1);
				}
				
				if (!tfoBam.doesWorkingExist() || !tfoBai.doesWorkingExist()) {
					samplesToRun.add(si); 
				} else {
					samplesToPostProcess.add(si);
				}
			}
		}
		
		if (samplesToRun.size() > 0) {
			int counter = 1;
			this.daemon = new TFThreadDaemon(this.logFile,this.commandString,samplesToRun.size(),this.jobs);
			this.daemon.start();
			
			for (TFSampleInfo si: samplesToRun) {
				//Create sample-specific storage
				ArrayList<File> protectList = new ArrayList<File>(); //don't remove these files on job re-submission cleanup
				
				//Get necessary file objects.
				TFFileObject tfoBam = si.getFileObject(TFConstants.FILE_BAM);
				TFFileObject tfoFastq1 = si.getFileObject(TFConstants.FILE_FASTQ1);
				TFFileObject tfoFastq2 = si.getFileObject(TFConstants.FILE_FASTQ2);
			
				//Create run directory
				samplesToPostProcess.add(si);
				File runDirectory = tfoBam.getWorkingDirectory();
				if (runDirectory.exists()) {
					this.deleteFolder(runDirectory);
				}
				runDirectory.mkdir();
				
				//Create files
				File cmdFile = new File(runDirectory,"cmd.txt");
				File fastq1 = tfoFastq1.createDestForFileObject(runDirectory,new String[]{"fastq"},new String[]{"txt"});
				File fastq2 = tfoFastq2.createDestForFileObject(runDirectory,new String[]{"fastq"},new String[]{"txt"});
				
				//Link necessary files
				this.createLink(tfoFastq1.getFinalPath(),fastq1);
				this.createLink(tfoFastq2.getFinalPath(),fastq2);
				
				//Create replacement hash
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("NAME", si.getSampleID());
				replacements.put("LIBRARY",si.getSampleName());
				replacements.put("SAMPLE",si.getSampleName());
				replacements.put("FLOWCELL",si.getPuID());
				if (si.isQual64()) {
					replacements.put("QUAL_FLAG", "-I");
				} else {
					replacements.put("QUAL_FLAG","");
				}
				
				//Add properties to replacement hash
				replacements.putAll(this.properties);
				
				//Create cmd.txt file
				this.createCmd(replacements,cmdFile,templateIdx);
				
				//Mark files for deletion or cleanup protection
				protectList.add(fastq1);
				protectList.add(fastq2);
				protectList.add(cmdFile);
				
				deleteList.add(fastq1);
				deleteList.add(fastq2);
				
				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
				
				//this.taskList.add(thread);
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
		
		
		if (samplesToPostProcess.size() > 0) {
			
			//Create final directories
			this.finalDirectory.mkdirs();
			this.jobDirectory.mkdir();
			
			HashSet<File> runDirectoryList = new HashSet<File>();
			
			for (TFSampleInfo si: samplesToPostProcess) {
				//Get sample objects
				TFFileObject tfoBam = si.getFileObject(TFConstants.FILE_BAM);
				TFFileObject tfoBai = si.getFileObject(TFConstants.FILE_BAI);
				
				//Add for cleanup
				runDirectoryList.add(tfoBam.getWorkingDirectory());
				
				this.moveFile(tfoBam.getWorkingPath(), tfoBam.getFinalPath());
				this.moveFile(tfoBai.getWorkingPath(), tfoBai.getFinalPath());
			}
			
			this.cleanup(runDirectoryList, deleteList);	
		}
		
		
	}
	
	
	

	
	public void mergeSampleLanes(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		//This method expects FILE_SPLIT_BAM
		//This method creates FILE_SAMPLE_BAM
		
		/**************************************************************************
		 *  Merge Sample Bams across lanes
		 ************************************************************************/
		
		logFile.writeInfoMessage("[TFExomeAlignUgp] Merging samples across lanes");
		
		//Create project-specific storage
		HashSet<File> deleteList = new HashSet<File>(); //List of files to delete at the end of the run
		HashSet<String> samplesToRun = new HashSet<String>();
		HashSet<String> samplesToPostProcess = new HashSet<String>();
		HashSet<String> samplesToRename = new HashSet<String>();
		
		/**************************************************************************
		 * Group by sample
		 **************************************************************************/
		
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		/**************************************************************************
		 * Determine stage of each sample
		 **************************************************************************/
		
		for (String sampleName: samples.keySet()) {
			//Create working Directory
			File workingDir = new File(this.rootDirectory,"JOB_" + sampleName + "_merge");
			
			boolean skipFile = true;
			ArrayList<TFSampleInfo> al = samples.get(sampleName);
			if (al.size() > 1) {
				skipFile = false;
			}
			
			for (TFSampleInfo si: samples.get(sampleName)) {
				//Create output objects for each sample
				TFFileObject tfoSampleBam = new TFFileObject(si.getSampleName() + ".raw.bam",this.finalDirectory,workingDir);
				TFFileObject tfoSampleBai = new TFFileObject(si.getSampleName() + ".raw.bai",this.finalDirectory,workingDir);
				
				//Check for existence of downstream output files.
				if (tfoSampleBam.doesFinalExist() && tfoSampleBai.doesFinalExist()) {
					si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoSampleBam);
					si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoSampleBai);
				} else {
					si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoSampleBam);
					si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoSampleBai);
					
					if (!si.getFileObject(TFConstants.FILE_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_BAI).doesFinalExist()) {
						logFile.writeErrorMessage("[TFExomeAlignUgp] Raw alignment files weren't found by TF.", true);
						System.exit(1);
					}
					
					if (skipFile) {
						samplesToRename.add(sampleName);
					} else if (!tfoSampleBam.doesWorkingExist() || !tfoSampleBai.doesWorkingExist()) {
						samplesToRun.add(sampleName);
					} else {
						samplesToPostProcess.add(sampleName);
					}
				}
			}
		}
		
		/**************************************************************************
		 * Rename bam file
		 **************************************************************************/
		
		if (samplesToRename.size() > 0) {
			for (String sampleName: samplesToRename) {
				TFSampleInfo tfsi = samples.get(sampleName).get(0);
					
				this.moveFile(tfsi.getFileObject(TFConstants.FILE_BAM).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAM).getFinalPath());
				this.moveFile(tfsi.getFileObject(TFConstants.FILE_BAI).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAI).getFinalPath());
			}
		}
		
		/**************************************************************************
		 * Process Samples
		 **************************************************************************/
		
		if (samplesToRun.size() > 0) {
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(templateIdx).getName(),samplesToRun.size(),this.jobs);
			this.daemon.start();
		
			int counter=1;
			for(String sampleName: samplesToRun) {
				//add samples for post-processing
				samplesToPostProcess.add(sampleName);
				
				//Get representative sample
				ArrayList<TFSampleInfo> al = samples.get(sampleName);
				TFSampleInfo repSI = al.get(0);
				TFFileObject tfoSampleBam = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAM);
				
				//Create run directory
				File runDirectory = tfoSampleBam.getWorkingDirectory();
				if (runDirectory.exists()) {
					this.deleteFolder(runDirectory);
				}
				runDirectory.mkdir();
				
				ArrayList<File> protectList = new ArrayList<File>();
				String mergeList = "";
				
				//Prep working directory.
				for (TFSampleInfo tfsi: samples.get(sampleName)) {
					//Create input objects
					TFFileObject tfoLaneBam = tfsi.getFileObject(TFConstants.FILE_BAM);
					TFFileObject tfoLaneBai = tfsi.getFileObject(TFConstants.FILE_BAI);
					
					//Create local versions of input files
					File fileLaneBam = tfoLaneBam.createDestForFileObject(runDirectory);
					File fileLaneBai = tfoLaneBai.createDestForFileObject(runDirectory);
					
					deleteList.add(fileLaneBam);
					deleteList.add(fileLaneBai);
					protectList.add(fileLaneBam);
					protectList.add(fileLaneBai);
					
					mergeList += " INPUT=" + fileLaneBam.getName();
					
					//Create links
					this.createLink(tfoLaneBam.getFinalPath(), fileLaneBam);
					this.createLink(tfoLaneBai.getFinalPath(), fileLaneBai);
				}
				
				//create properties list	
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("INPUT_LIST", mergeList);
				replacements.put("OUTPUT", tfoSampleBam.getFileName());
				replacements.putAll(this.properties);
				
				//create files
				File cmdFile = new File(runDirectory,"cmd.txt");
				this.createCmd(replacements,cmdFile,templateIdx);
				protectList.add(cmdFile);
				
				//Start job
				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
				this.daemon.addJob(thread);
				counter++;
	
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
			
			
			for(String sampleName: samplesToRun) {
				for (TFSampleInfo tfsi: samples.get(sampleName)) {
					//Create input objects
					TFFileObject tfoLaneBam = tfsi.getFileObject(TFConstants.FILE_BAM);
					TFFileObject tfoLaneBai = tfsi.getFileObject(TFConstants.FILE_BAI);
					
					this.deleteFile(tfoLaneBam.getFinalPath());
					this.deleteFile(tfoLaneBai.getFinalPath());
				}
			}
			
			
			
			
		}
		
		
		/**************************************************************************
		 * Post-Process
		 ************************************************************************/
		
		if (samplesToPostProcess.size() > 0) {
			//Create final directories
			this.finalDirectory.mkdir();
			this.jobDirectory.mkdir();
			
			HashSet<File> runDirectoryList = new HashSet<File>();
			for (String sampleName: samplesToPostProcess) {				
				//Get representative si
				TFSampleInfo repSI = samples.get(sampleName).get(0);
				
				TFFileObject tfoSampleBam = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAM);
				TFFileObject tfoSampleBai = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAI);
				
				File workingDir = tfoSampleBam.getWorkingDirectory();
				runDirectoryList.add(workingDir);
				
				File fileSampleBam = tfoSampleBam.createDestForFileObject(workingDir);
				File fileSampleBai = tfoSampleBai.createDestForFileObject(workingDir);
				
				this.moveFile(fileSampleBam, tfoSampleBam.getFinalPath());
				this.moveFile(fileSampleBai, tfoSampleBai.getFinalPath());
			}
			
			this.cleanup(runDirectoryList, deleteList);
		}
	}
	
	
	
	

	
	
	
	public void postProcessSample(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		//This method expects FILE_SAMPLE_BAM
		//This method creates FILE_SAMPLE_REALIGN
				
		/**************************************************************************
		 * Realign merged sample
		 ************************************************************************/
		
		logFile.writeInfoMessage("[TFExomeAlignUgp] Postprocessing merged samples");
		
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
			File workingDir = new File(this.rootDirectory,"JOB_" + sampleName + "_postprocess"); //Create working directory
			
			for (TFSampleInfo si: samples.get(sampleName)) {
				
				
				//Create output objects for each sample
				TFFileObject tfoFinalBam = new TFFileObject(si.getSampleName() + ".final.bam",this.finalDirectory,workingDir);
				TFFileObject tfoFinalBai = new TFFileObject(si.getSampleName() + ".final.bai",this.finalDirectory,workingDir);
				si.setFileObject(TFConstants.FILE_FINAL_BAM, tfoFinalBam);
				si.setFileObject(TFConstants.FILE_FINAL_BAI, tfoFinalBai);
				
				//Check for the existance of output files
				if (!tfoFinalBam.doesFinalExist() || !tfoFinalBai.doesFinalExist()) {
					if (!si.getFileObject(TFConstants.FILE_SAMPLE_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_SAMPLE_BAI).doesFinalExist()) {
						logFile.writeErrorMessage("[TFExomeAlignUgp] Merged sample files weren't found by TF", true);
						System.exit(1);
					}
					
					if (!tfoFinalBam.doesWorkingExist() || !tfoFinalBai.doesWorkingExist()) {
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
				TFFileObject tfoFinalBam = repSI.getFileObject(TFConstants.FILE_FINAL_BAM);
				TFFileObject tfoSampleBam = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAM);
				TFFileObject tfoSampleBai = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAI);
				
				//Create run directory
				File runDirectory = tfoFinalBam.getWorkingDirectory();
				if (runDirectory.exists()) {
					this.deleteFolder(runDirectory);
				}
				runDirectory.mkdir();
				
				//Create working files
				File fileSampleBam = tfoSampleBam.createDestForFileObject(runDirectory);
				File fileSampleBai = tfoSampleBai.createDestForFileObject(runDirectory);
				
				//Create links
				this.createLink(tfoSampleBam.getFinalPath(), fileSampleBam);
				this.createLink(tfoSampleBai.getFinalPath(), fileSampleBai);
				
				//create properties
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("NAME", sampleName);
				replacements.putAll(this.properties);
				
				//create files
				File cmdFile = new File(runDirectory,"cmd.txt");
				this.createCmd(replacements,cmdFile,templateIdx);
				
				deleteList.add(fileSampleBam);
				deleteList.add(fileSampleBai);
				protectList.add(fileSampleBam);
				protectList.add(fileSampleBai);
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
			this.jobDirectory.mkdir();
			
			HashSet<File> runDirectoryList = new HashSet<File>();
			
			for (String sampleName: samplesToPostProcess) {
				 //Get representative sample
				TFSampleInfo repSI = samples.get(sampleName).get(0);
				
				//Get final bam file objects
			    TFFileObject tfoFinalBam = repSI.getFileObject(TFConstants.FILE_FINAL_BAM);
			    TFFileObject tfoFinalBai = repSI.getFileObject(TFConstants.FILE_FINAL_BAI);
			    
			    File runDirectory = tfoFinalBam.getWorkingDirectory();
			    runDirectoryList.add(runDirectory);
			    
			    File fileFinalBam = tfoFinalBam.createDestForFileObject(runDirectory);
			    File fileFinalBai = tfoFinalBai.createDestForFileObject(runDirectory);
			    
			    this.moveFile(fileFinalBam, tfoFinalBam.getFinalPath());
			    this.moveFile(fileFinalBai, tfoFinalBai.getFinalPath());
			}
			
			this.cleanup(runDirectoryList, deleteList);
		}
	}
	
}
