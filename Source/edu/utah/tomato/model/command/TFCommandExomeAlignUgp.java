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
			Integer jobs, boolean suppress, boolean isFull, boolean validateFastq, HashMap<String,String> properties) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress, isFull, properties);
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
		
		this.alignment(sampleList,0);
		this.mergeSampleLanes(sampleList,1);
		this.realignSample(sampleList,2);
		this.reduceSample(sampleList,3);
		
//		this.alignment(sampleList,0);
//		this.realignLane(sampleList,1);
//		this.splitLaneBams(sampleList,2);
//		this.mergeSampleLanes(sampleList,3);
//		this.realignSample(sampleList,3);
//		this.reduceSample(sampleList,4);
		
		
		
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
		
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_FASTQ1, Pattern.compile("(.+?)_1\\..+\\.gz$"), TFConstants.PREFIX_SAMPLEID));
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_FASTQ2, Pattern.compile("(.+?)_2\\..+\\.gz$"), TFConstants.PREFIX_SAMPLEID));
		
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
			TFFileObject tfoLaneBam = new TFFileObject(si.getSampleID() + ".realign.bam",this.finalDirectory,workingDir);
			TFFileObject tfoLaneBai = new TFFileObject(si.getSampleID() + ".realign.bai",this.finalDirectory,workingDir);
			TFFileObject tfoFinalBam = new TFFileObject(si.getSampleName() + ".final.bam",this.finalDirectory,workingDir);
			TFFileObject tfoFinalBai = new TFFileObject(si.getSampleName() + ".final.bai",this.finalDirectory,workingDir);
			si.setFileObject(TFConstants.FILE_BAM,tfoBam);
			si.setFileObject(TFConstants.FILE_BAI,tfoBai);
			si.setFileObject(TFConstants.FILE_LANE_BAM,tfoLaneBam);
			si.setFileObject(TFConstants.FILE_LANE_BAI,tfoLaneBai);
			
			//Check to see if the output files already exist
			if ((tfoBam.doesFinalExist() && tfoBai.doesFinalExist()) && ((tfoLaneBam.doesFinalExist() && !tfoLaneBai.doesFinalExist()) ||
					(tfoFinalBam.doesFinalExist() && tfoFinalBai.doesFinalExist()))) {
				//OK
			} else{
				//Make sure files are initialized
				if (!si.getFileObject(TFConstants.FILE_FASTQ1).doesFinalExist() || !si.getFileObject(TFConstants.FILE_FASTQ2).doesFinalExist()) {
					logFile.writeErrorMessage("[TFExomeAlignUgp] Your fastq files weren't found by TF",true);
					System.exit(1);
				}
				
				if (!tfoBam.doesWorkingExist() || !tfoBai.doesWorkingExist() || !tfoLaneBam.doesWorkingExist() || !tfoLaneBai.doesWorkingExist()) {
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
				TFFileObject tfoLaneBam = si.getFileObject(TFConstants.FILE_LANE_BAM);
				TFFileObject tfoLaneBai = si.getFileObject(TFConstants.FILE_LANE_BAI);
				
				//Add for cleanup
				runDirectoryList.add(tfoBam.getWorkingDirectory());
				
				this.moveFile(tfoBam.getWorkingPath(), tfoBam.getFinalPath());
				this.moveFile(tfoBai.getWorkingPath(), tfoBai.getFinalPath());
				this.moveFile(tfoLaneBam.getWorkingPath(), tfoLaneBam.getFinalPath());
				this.moveFile(tfoLaneBai.getWorkingPath(), tfoLaneBai.getFinalPath());
				
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
				TFFileObject tfoFinalBam = new TFFileObject(si.getSampleName() + ".final.bam",this.finalDirectory,workingDir);
				TFFileObject tfoFinalBai = new TFFileObject(si.getSampleName() + ".final.bai",this.finalDirectory,workingDir);
				TFFileObject tfoSampleBam = new TFFileObject(si.getSampleName() + ".sample.bam",this.finalDirectory,workingDir);
				TFFileObject tfoSampleBai = new TFFileObject(si.getSampleName() + ".sample.bai",this.finalDirectory,workingDir);
				
				//Check for existence of downstream output files.
				if (tfoFinalBam.doesFinalExist() && tfoFinalBai.doesFinalExist()) {
					si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoFinalBam);
					si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoFinalBai);
				} else if (tfoSampleBam.doesFinalExist() && tfoSampleBai.doesFinalExist()) {
					si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoSampleBam);
					si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoSampleBai);
				} else {
					si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoSampleBam);
					si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoSampleBai);
					
					if (!si.getFileObject(TFConstants.FILE_LANE_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_LANE_BAI).doesFinalExist()) {
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
					
				this.moveFile(tfsi.getFileObject(TFConstants.FILE_LANE_BAM).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAM).getFinalPath());
				this.moveFile(tfsi.getFileObject(TFConstants.FILE_LANE_BAI).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAI).getFinalPath());
			}
		}
		
		/**************************************************************************
		 * Process Samples
		 **************************************************************************/
		
		if (samplesToRun.size() > 0) {
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(templateIdx).getName(),samplesToRun.size(),this.jobs);
			this.daemon.start();
		
			int counter=1;
			for(String sampleName: samples.keySet()) {
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
					TFFileObject tfoLaneBam = tfsi.getFileObject(TFConstants.FILE_LANE_BAM);
					TFFileObject tfoLaneBai = tfsi.getFileObject(TFConstants.FILE_LANE_BAI);
					
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
	
	
	
	

	
	
	
	public void realignSample(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		//This method expects FILE_SAMPLE_BAM
		//This method creates FILE_SAMPLE_REALIGN
				
		/**************************************************************************
		 * Realign merged sample
		 ************************************************************************/
		
		logFile.writeInfoMessage("[TFExomeAlignUgp] Realigning merged samples");
		
		//Containers for bam-splitting processing 
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<String> samplesToRun = new HashSet<String>();
		HashSet<String> samplesToPostProcess = new HashSet<String>();
		HashSet<String> samplesToRename = new HashSet<String>();
		
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
			File workingDir = new File(this.rootDirectory,"JOB_" + sampleName + "_realign"); //Create working directory
			
			boolean skipFile = true;
			ArrayList<TFSampleInfo> al  = samples.get(sampleName);
			if (al.size() > 1) {
				skipFile = false;
			}
			
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
					if (skipFile) {
						samplesToRename.add(sampleName);
					} else if (!tfoFinalBam.doesWorkingExist() || !tfoFinalBai.doesWorkingExist()) {
						samplesToRun.add(sampleName);
					} else {
						samplesToPostProcess.add(sampleName);
					}
				}
			}
		}
		
		/**************************************************************************
		 * Rename samples
		 ************************************************************************/
		
		if ( samplesToRename.size() > 0 ) {
			for (String sampleName: samplesToRename) {
				TFSampleInfo tfsi = samples.get(sampleName).get(0);
				
				this.moveFile(tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAM).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_FINAL_BAM).getFinalPath());
				this.moveFile(tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAI).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_FINAL_BAI).getFinalPath());
			}
		}
		
		/**************************************************************************
		 * Process samples
		 ************************************************************************/
		
		if (samplesToRun.size() > 0) {
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(templateIdx).getName(),samplesToRun.size(),this.jobs);
			this.daemon.start();
			
			int counter = 1;
			
			for(String sampleName: samples.keySet()) {
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
			    
			    //Delete merged sample bam
			    this.deleteFile(repSI.getFileObject(TFConstants.FILE_SAMPLE_BAM).getFinalPath());
			    this.deleteFile(repSI.getFileObject(TFConstants.FILE_SAMPLE_BAI).getFinalPath());
			}
			
			this.cleanup(runDirectoryList, deleteList);
		}
	}
	
	public void reduceSample(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		//This method expects SAMPLE_REALIGN_BAM
		//This method creates SAMPLE_REDUCE_BAM
		
		/**************************************************************************
		 * Create reduced bams ( for variant calling )
		 ************************************************************************/
		
		logFile.writeInfoMessage("[TFExomeAlignUgp] Reducing sample bams");
		
		//Containers for bam-splitting processing 
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<String> samplesToRun = new HashSet<String>();
		HashSet<String> samplesToPostProcess = new HashSet<String>();

		//Group samples by Lane
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		/**************************************************************************
		 * Determine Sample Stage
		 ************************************************************************/
		
		for (String sampleName: samples.keySet()) {
			//Create working directory
			File workingDir = new File(this.rootDirectory,"JOB_" + sampleName + "_reduce");
			
			for (TFSampleInfo tfsi: samples.get(sampleName)) {
				
				
				//Create output objects for each sample
				TFFileObject tfoReduceBam = new TFFileObject(tfsi.getSampleName() + ".reduce.bam",this.finalDirectory,workingDir);
				TFFileObject tfoReduceBai = new TFFileObject(tfsi.getSampleName() + ".reduce.bai",this.finalDirectory,workingDir);
				tfsi.setFileObject(TFConstants.FILE_REDUCE_BAM, tfoReduceBam);
				tfsi.setFileObject(TFConstants.FILE_REDUCE_BAI, tfoReduceBai);
				
				//Check for existance of output files
				if (!tfoReduceBam.doesFinalExist() || !tfoReduceBai.doesFinalExist()) {
					if (!tfsi.getFileObject(TFConstants.FILE_FINAL_BAM).doesFinalExist() || !tfsi.getFileObject(TFConstants.FILE_FINAL_BAI).doesFinalExist()) {
						logFile.writeErrorMessage("[TFExomeAlignUgp] Final bam files weren't found by TF", true);
						System.exit(1);
					}
					
					if (!tfoReduceBam.doesWorkingExist() || !tfoReduceBai.doesWorkingExist()) {
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
			
			for(String sampleName: samples.keySet()) {
				//Get representative sample
				TFSampleInfo siRep = samples.get(sampleName).get(0);
				
				//Create input/output objecst
				TFFileObject tfoReduceBam = siRep.getFileObject(TFConstants.FILE_REDUCE_BAM);
				TFFileObject tfoFinalBam = siRep.getFileObject(TFConstants.FILE_FINAL_BAM);
				TFFileObject tfoFinalBai = siRep.getFileObject(TFConstants.FILE_FINAL_BAI);
				
				//Create job directory
				File runDirectory = tfoReduceBam.getWorkingDirectory();
				if (runDirectory.exists()) {
					this.deleteFolder(runDirectory);
				}
				runDirectory.mkdir();
				samplesToPostProcess.add(sampleName);
				
				//Create files
				File fileFinalBam = tfoFinalBam.createDestForFileObject(runDirectory);
				File fileFinalBai = tfoFinalBai.createDestForFileObject(runDirectory);
				
				//Create links
				this.createLink(tfoFinalBam.getFinalPath(),fileFinalBam);
				this.createLink(tfoFinalBai.getFinalPath(),fileFinalBai);
				
				//Create run-specific protection
				ArrayList<File> protectList = new ArrayList<File>();
				
				//create properties list
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("NAME", sampleName);
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
			this.jobDirectory.mkdir();
			
			HashSet<File> runDirectoryList = new HashSet<File>();
			
			for (String sampleName: samplesToPostProcess) {
				TFSampleInfo repSI = samples.get(sampleName).get(0); //get representative sample
				
				TFFileObject tfoReduceBam = repSI.getFileObject(TFConstants.FILE_REDUCE_BAM);
				TFFileObject tfoReduceBai = repSI.getFileObject(TFConstants.FILE_REDUCE_BAI);
				
				File workingDir = tfoReduceBam.getWorkingDirectory();
				runDirectoryList.add(workingDir);
				
				File fileReduceBam = tfoReduceBam.createDestForFileObject(workingDir);
				File fileReduceBai = tfoReduceBai.createDestForFileObject(workingDir);
				
				this.moveFile(fileReduceBam, tfoReduceBam.getFinalPath());
				this.moveFile(fileReduceBai, tfoReduceBai.getFinalPath());
				
			}
			
			this.cleanup(runDirectoryList, deleteList);
		}
	}
	

/*
 * This commented out code was used when we did lane-level (vs library-level) recalibration.  It wasn't wrong, but not as efficient and more
 * complicated.
 */
//	public void mergeSampleLanes(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
//	//This method expects FILE_SPLIT_BAM
//	//This method creates FILE_SAMPLE_BAM
//	
//	/**************************************************************************
//	 *  Merge Sample Bams across lanes
//	 ************************************************************************/
//	
//	logFile.writeInfoMessage("[TFExomeAlignUgp] Merging samples across lanes");
//	
//	//Create project-specific storage
//	HashSet<File> deleteList = new HashSet<File>(); //List of files to delete at the end of the run
//	HashSet<String> samplesToRun = new HashSet<String>();
//	HashSet<String> samplesToPostProcess = new HashSet<String>();
//	HashSet<String> samplesToRename = new HashSet<String>();
//	
//	/**************************************************************************
//	 * Group by sample
//	 **************************************************************************/
//	
//	HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
//	
//	for (TFSampleInfo si: sampleList) {
//		if (!samples.containsKey(si.getSampleName())) {
//			samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
//		}
//		samples.get(si.getSampleName()).add(si);
//	}
//	
//	/**************************************************************************
//	 * Determine stage of each sample
//	 **************************************************************************/
//	
//	for (String sampleName: samples.keySet()) {
//		//Create working Directory
//		File workingDir = new File(this.rootDirectory,"JOB_" + sampleName + "_merge");
//		
//		boolean skipFile = true;
//		ArrayList<TFSampleInfo> al = samples.get(sampleName);
//		if (al.size() > 1) {
//			skipFile = false;
//		}
//		
//		for (TFSampleInfo si: samples.get(sampleName)) {
//			if (!si.getFileObject(TFConstants.FILE_SPLIT_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_SPLIT_BAI).doesFinalExist()) {
//				logFile.writeErrorMessage("[TFExomeAlignUgp] Split sample bam files weren't found by TF", true);
//				System.exit(1);
//			}
//			
//			
//			
//			//Create output objects for each sample
//			TFFileObject tfoFinalBam = new TFFileObject(si.getSampleName() + ".final.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoFinalBai = new TFFileObject(si.getSampleName() + ".final.bai",this.finalDirectory,workingDir);
//			TFFileObject tfoSampleBam = new TFFileObject(si.getSampleName() + ".sample.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoSampleBai = new TFFileObject(si.getSampleName() + ".sample.bai",this.finalDirectory,workingDir);
//			
//			//Check for existence of downstream output files.
//			if (tfoFinalBam.doesFinalExist() && tfoFinalBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoFinalBam);
//				si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoFinalBai);
//			} else if (tfoSampleBam.doesFinalExist() && tfoSampleBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoSampleBam);
//				si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoSampleBai);
//			} else {
//				si.setFileObject(TFConstants.FILE_SAMPLE_BAM,tfoSampleBam);
//				si.setFileObject(TFConstants.FILE_SAMPLE_BAI,tfoSampleBai);
//				
//				if (skipFile) {
//					samplesToRename.add(sampleName);
//				} else if (!tfoSampleBam.doesWorkingExist() || !tfoSampleBai.doesWorkingExist()) {
//					samplesToRun.add(sampleName);
//				} else {
//					samplesToPostProcess.add(sampleName);
//				}
//			}
//		}
//	}
//	
//	/**************************************************************************
//	 * Rename bam file
//	 **************************************************************************/
//	
//	if (samplesToRename.size() > 0) {
//		for (String sampleName: samplesToRename) {
//			TFSampleInfo tfsi = samples.get(sampleName).get(0);
//				
//			this.moveFile(tfsi.getFileObject(TFConstants.FILE_SPLIT_BAM).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAM).getFinalPath());
//			this.moveFile(tfsi.getFileObject(TFConstants.FILE_SPLIT_BAI).getFinalPath(), tfsi.getFileObject(TFConstants.FILE_SAMPLE_BAI).getFinalPath());
//		}
//	}
//	
//	/**************************************************************************
//	 * Process Samples
//	 **************************************************************************/
//	
//	if (samplesToRun.size() > 0) {
//		this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(3).getName(),samplesToRun.size(),this.jobs);
//		this.daemon.start();
//	
//		int counter=1;
//		for(String sampleName: samples.keySet()) {
//			//add samples for post-processing
//			samplesToPostProcess.add(sampleName);
//			
//			//Get representative sample
//			ArrayList<TFSampleInfo> al = samples.get(sampleName);
//			TFSampleInfo repSI = al.get(0);
//			TFFileObject tfoSampleBam = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAM);
//			
//			//Create run directory
//			File runDirectory = tfoSampleBam.getWorkingDirectory();
//			if (runDirectory.exists()) {
//				this.deleteFolder(runDirectory);
//			}
//			runDirectory.mkdir();
//			
//			ArrayList<File> protectList = new ArrayList<File>();
//			String mergeList = "";
//			
//			//Prep working directory.
//			for (TFSampleInfo tfsi: samples.get(sampleName)) {
//				//Create input objects
//				TFFileObject tfoSplitBam = tfsi.getFileObject(TFConstants.FILE_SPLIT_BAM);
//				TFFileObject tfoSplitBai = tfsi.getFileObject(TFConstants.FILE_SPLIT_BAI);
//				
//				//Create local versions of input files
//				File fileSplitBam = tfoSplitBam.createDestForFileObject(runDirectory);
//				File fileSplitBai = tfoSplitBai.createDestForFileObject(runDirectory);
//				
//				deleteList.add(fileSplitBam);
//				deleteList.add(fileSplitBai);
//				protectList.add(fileSplitBam);
//				protectList.add(fileSplitBai);
//				
//				mergeList += " INPUT=" + fileSplitBam.getName();
//				
//				//Create links
//				this.createLink(tfoSplitBam.getFinalPath(), fileSplitBam);
//				this.createLink(tfoSplitBai.getFinalPath(), fileSplitBai);
//			}
//			
//			//create properties list	
//			HashMap<String,String> replacements = new HashMap<String,String>();
//			replacements.put("INPUT_LIST", mergeList);
//			replacements.put("OUTPUT", tfoSampleBam.getFileName());
//			replacements.putAll(this.properties);
//			
//			//create files
//			File cmdFile = new File(runDirectory,"cmd.txt");
//			this.createCmd(replacements,cmdFile,templateIdx);
//			protectList.add(cmdFile);
//			
//			//Start job
//			TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
//			this.daemon.addJob(thread);
//			counter++;
//
//		}
//			
//		//Wait for command to finish
//		try {
//			this.daemon.join();
//			Thread.sleep(5000);
//			if (this.daemon.getFailed()) {
//				System.exit(1);
//			}
//		} catch (InterruptedException ie) {
//			logFile.writeErrorMessage("[TFExomeAlignUgp] Daemon interrupted",true);
//			System.exit(1);
//		}
//		
//		
//	}
//	
//	
//	/**************************************************************************
//	 * Post-Process
//	 ************************************************************************/
//	
//	if (samplesToPostProcess.size() > 0) {
//		//Create final directories
//		this.finalDirectory.mkdir();
//		this.jobDirectory.mkdir();
//		
//		HashSet<File> runDirectoryList = new HashSet<File>();
//		for (String sampleName: samplesToPostProcess) {				
//			//Get representative si
//			TFSampleInfo repSI = samples.get(sampleName).get(0);
//			
//			TFFileObject tfoSampleBam = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAM);
//			TFFileObject tfoSampleBai = repSI.getFileObject(TFConstants.FILE_SAMPLE_BAI);
//			
//			File workingDir = tfoSampleBam.getWorkingDirectory();
//			runDirectoryList.add(workingDir);
//			
//			File fileSampleBam = tfoSampleBam.createDestForFileObject(workingDir);
//			File fileSampleBai = tfoSampleBai.createDestForFileObject(workingDir);
//			
//			this.moveFile(fileSampleBam, tfoSampleBam.getFinalPath());
//			this.moveFile(fileSampleBai, tfoSampleBai.getFinalPath());
//		}
//		
//		this.cleanup(runDirectoryList, deleteList);
//	}
//}
	
	
	
//	public void realignLane(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
//	/*
//	 * This method expects the input files TFConstants.FILE_BAM
//	 * This method creates the output files TFContants.REALIGN_BAM
//	 */
//	
//	/**************************************************************************
//	 *  Realign and recalibrate by lane
//	 ************************************************************************/
//	
//	logFile.writeInfoMessage("[TFExomeAlignUgp] Running realignment and recalibration on each lane");
//	
//	
//	/**************************************************************************
//	 * Group raw alignments by lane
//	 **************************************************************************/
//	
//	//Group samples by Lane
//	HashMap<String,ArrayList<TFSampleInfo>> lanes = new HashMap<String,ArrayList<TFSampleInfo>>();
//	for (TFSampleInfo si: sampleList) {
//		if (!lanes.containsKey(si.getPuID())) {
//			lanes.put(si.getPuID(), new ArrayList<TFSampleInfo>());
//		}
//		lanes.get(si.getPuID()).add(si);
//	}
//	
//	/**************************************************************************
//	 * Determine stage of each lane
//	 **************************************************************************/
//	
//	//Create project-specific storage
//	HashSet<File> deleteList = new HashSet<File>(); //list of files to delete at the end of the run
//	HashSet<String> lanesToRun = new HashSet<String>();
//	HashSet<String> lanesToPostProcess = new HashSet<String>();
//			
//	//Iterate through each lane/sample.  Make sure pre-requisites are met.  Check to see if job already finished.
//	for (String laneName: lanes.keySet()) {
//		//Create working directory
//		File workingDir = new File(this.rootDirectory,"JOB_" + laneName + "_realign");
//		
//		for (TFSampleInfo si: lanes.get(laneName)) {
//			if (!si.getFileObject(TFConstants.FILE_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_BAI).doesFinalExist()) {
//				logFile.writeErrorMessage("[TFExomeAlignUgp] Raw bam files weren't found by TF", true);
//				System.exit(1);
//			}
//			
//			//Create output objects for each sample
//			TFFileObject tfoLaneBam = new TFFileObject(si.getPuID() + ".realign.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoLaneBai = new TFFileObject(si.getPuID() + ".realign.bai",this.finalDirectory,workingDir);
//			TFFileObject tfoSplitBam = new TFFileObject(si.getSampleID() + ".split.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoSplitBai = new TFFileObject(si.getSampleID() + ".split.bai",this.finalDirectory,workingDir);
//			TFFileObject tfoFinalBam = new TFFileObject(si.getSampleName() + ".final.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoFinalBai = new TFFileObject(si.getSampleName() + ".final.bai",this.finalDirectory,workingDir);
//			TFFileObject tfoSampleBam = new TFFileObject(si.getSampleName() + ".sample.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoSampleBai = new TFFileObject(si.getSampleName() + ".sample.bai",this.finalDirectory,workingDir);
//			
//			//Check for existence of downstream output files.
//			if (tfoLaneBam.doesFinalExist() && tfoLaneBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_LANE_BAM,tfoLaneBam);
//				si.setFileObject(TFConstants.FILE_LANE_BAI,tfoLaneBai);
//			} else if (tfoSplitBam.doesFinalExist() && tfoSplitBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_LANE_BAM,tfoSplitBam);
//				si.setFileObject(TFConstants.FILE_LANE_BAI,tfoSplitBai);
//			} else if (tfoFinalBam.doesFinalExist() && tfoFinalBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_LANE_BAM,tfoFinalBam);
//				si.setFileObject(TFConstants.FILE_LANE_BAI,tfoFinalBai);
//			} else if (tfoSampleBam.doesFinalExist() && tfoSampleBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_LANE_BAM,tfoSampleBam);
//				si.setFileObject(TFConstants.FILE_LANE_BAI,tfoSampleBai);
//			} else {
//				si.setFileObject(TFConstants.FILE_LANE_BAM,tfoLaneBam);
//				si.setFileObject(TFConstants.FILE_LANE_BAI,tfoLaneBai);
//				if (!tfoLaneBam.doesWorkingExist() || !tfoLaneBam.doesWorkingExist()) {
//					lanesToRun.add(laneName);
//				} else {
//					lanesToPostProcess.add(laneName);
//				}
//			}
//		}	
//	}
//	
//	
//	/**************************************************************************
//	 * Process files
//	 **************************************************************************/
//	
//	if (lanesToRun.size() > 0) {
//		int counter = 1;
//		this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(1).getName(),lanesToRun.size(),this.jobs);
//		this.daemon.start();
//		
//		for(String laneName: lanes.keySet()) {
//			//Create lane-specific storage
//			ArrayList<File> protectList = new ArrayList<File>(); //don't remove these files on job re-submission cleanup
//		
//			//create run directory
//			File runDirectory = new File(this.rootDirectory,"JOB_" + laneName + "_realign");
//			if (runDirectory.exists()) {
//				this.deleteFolder(runDirectory);
//			}
//			runDirectory.mkdir();
//			
//			//Add lane to postprocess list
//			lanesToPostProcess.add(laneName);
//			
//			//Build dup input list and grab 
//			String dupList = "";
//			ArrayList<TFSampleInfo> laneBams = lanes.get(laneName);
//			
//			
//			//Create links in working directory to raw bam files
//			for (TFSampleInfo si: laneBams) {
//				File rawBamLocal = si.getFileObject(TFConstants.FILE_BAM).createDestForFileObject(runDirectory);
//				File rawBaiLocal = si.getFileObject(TFConstants.FILE_BAI).createDestForFileObject(runDirectory);
//				
//				protectList.add(rawBamLocal);
//				protectList.add(rawBaiLocal);
//				deleteList.add(rawBamLocal);
//				deleteList.add(rawBaiLocal);
//				
//				this.createLink(si.getFileObject(TFConstants.FILE_BAM).getFinalPath(),rawBamLocal);
//				this.createLink(si.getFileObject(TFConstants.FILE_BAI).getFinalPath(),rawBaiLocal);
//				
//				dupList += " INPUT=" + rawBamLocal.getName();
//		
//			}
//			
//			
//			//create properties list
//			HashMap<String,String> replacements = new HashMap<String,String>();
//			replacements.put("NAME", laneName);
//			replacements.put("DUP_LIST", dupList);
//			replacements.putAll(this.properties);
//			
//			//Create cmd.txt file
//			File cmdFile = new File(runDirectory,"cmd.txt");
//			this.createCmd(replacements,cmdFile,templateIdx);
//			protectList.add(cmdFile);
//			
//			TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
//			
//			this.daemon.addJob(thread);
//			
//			counter ++;
//			
//		}
//		
//		//Wait for command to finish
//		try {
//			this.daemon.join();
//			Thread.sleep(5000);
//			if (this.daemon.getFailed()) {
//				System.exit(1);
//			}
//		} catch (InterruptedException ie) {
//			logFile.writeErrorMessage("[TFExomeAlignUgp] Daemon interrupted",true);
//			System.exit(1);
//		}
//		
//		
//	} //end lane processing
//	
//	
//	/**************************************************************************
//	 * Post-Process files
//	 **************************************************************************/
//	
//	if (lanesToPostProcess.size() > 0) {
//		//Create final directories
//		this.finalDirectory.mkdir();
//		this.jobDirectory.mkdir();
//		
//		HashSet<File> runDirectoryList = new HashSet<File>();
//		
//		
//		for(String laneName: lanesToPostProcess) {			
//			//Get representative sample object
//			TFFileObject realignBam = lanes.get(laneName).get(0).getFileObject(TFConstants.FILE_LANE_BAM);
//			TFFileObject realignBai = lanes.get(laneName).get(0).getFileObject(TFConstants.FILE_LANE_BAI);
//			
//			//Move results
//			this.moveFile(realignBam.getWorkingPath(),realignBam.getFinalPath());
//			this.moveFile(realignBai.getWorkingPath(),realignBai.getFinalPath());
//			
//
//			runDirectoryList.add(realignBam.getWorkingDirectory());
//		}
//		
//		this.cleanup(runDirectoryList, deleteList);
//	}
//}

	
	

	

//public void splitLaneBams(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
//	/*
//	 * This method expects FILE_LANE_BAM
//	 * This method creates FILE_SPLIT_BAM
//	 */
//	
//	/**************************************************************************
//	 *  Split lane bams into sample bams
//	 ************************************************************************/
//	
//	logFile.writeInfoMessage("[TFExomeAlignUgp] Splitting lane bams into sample bams");
//	
//	/**************************************************************************
//	 * Group raw alignments by lane
//	 **************************************************************************/
//	
//	//Group samples by Lane
//	HashMap<String,ArrayList<TFSampleInfo>> lanes = new HashMap<String,ArrayList<TFSampleInfo>>();
//	for (TFSampleInfo si: sampleList) {
//		if (!lanes.containsKey(si.getPuID())) {
//			lanes.put(si.getPuID(), new ArrayList<TFSampleInfo>());
//		}
//		lanes.get(si.getPuID()).add(si);
//	}
//
//	
//	/**************************************************************************
//	 * Determine stage of each lane
//	 **************************************************************************/
//	
//	//Create project-specific storage
//	HashSet<File> deleteList = new HashSet<File>(); //list of files to delete at the end of the run
//	HashSet<String> lanesToRun = new HashSet<String>();
//	HashSet<String> lanesToPostProcess = new HashSet<String>();
//	HashSet<String> lanesToMove = new HashSet<String>();
//
//			
//	//Iterate through each lane/sample.  Make sure pre-requisites are met.  Check to see if job already finished.
//	for (String laneName: lanes.keySet()) {
//		//Create working directory
//		File workingDir = new File(this.rootDirectory,"JOB_" + laneName + "_split");
//		
//		boolean skipFile = true;
//		
//		//Determine number of daemon jobs to run
//		ArrayList<TFSampleInfo> al = lanes.get(laneName);
//		if (al.size() > 1) {
//			skipFile = false;
//		}
//		
//		for (TFSampleInfo si: lanes.get(laneName)) {
//			if (!si.getFileObject(TFConstants.FILE_LANE_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_LANE_BAI).doesFinalExist()) {
//				logFile.writeErrorMessage("[TFExomeAlignUgp] Lane-level bam files weren't found by TF", true);
//				System.exit(1);
//			}
//			
//			//Create output objects for each sample
//			TFFileObject tfoSplitBam = new TFFileObject(si.getSampleID() + ".split.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoSplitBai = new TFFileObject(si.getSampleID() + ".split.bai",this.finalDirectory,workingDir);
//			TFFileObject tfoFinalBam = new TFFileObject(si.getSampleName() + ".final.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoFinalBai = new TFFileObject(si.getSampleName() + ".final.bai",this.finalDirectory,workingDir);
//			TFFileObject tfoSampleBam = new TFFileObject(si.getSampleName() + ".sample.bam",this.finalDirectory,workingDir);
//			TFFileObject tfoSampleBai = new TFFileObject(si.getSampleName() + ".sample.bai",this.finalDirectory,workingDir);
//			
//			//Check for existence of downstream output files.
//			if (tfoSplitBam.doesFinalExist() && tfoSplitBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_SPLIT_BAM,tfoSplitBam);
//				si.setFileObject(TFConstants.FILE_SPLIT_BAI,tfoSplitBai);
//			} else if (tfoFinalBam.doesFinalExist() && tfoFinalBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_SPLIT_BAM,tfoFinalBam);
//				si.setFileObject(TFConstants.FILE_SPLIT_BAI,tfoFinalBai);
//			} else if (tfoSampleBam.doesFinalExist() && tfoSampleBai.doesFinalExist()) {
//				si.setFileObject(TFConstants.FILE_SPLIT_BAM,tfoSampleBam);
//				si.setFileObject(TFConstants.FILE_SPLIT_BAI,tfoSampleBai);
//			} else {
//				si.setFileObject(TFConstants.FILE_SPLIT_BAM,tfoSplitBam);
//				si.setFileObject(TFConstants.FILE_SPLIT_BAI,tfoSplitBai);
//				
//				if (skipFile) {
//					lanesToMove.add(laneName);
//				} else if (!tfoSplitBam.doesWorkingExist() || !tfoSplitBam.doesWorkingExist()) {
//					lanesToRun.add(laneName);
//				} else {
//					lanesToPostProcess.add(laneName);
//				}
//			}
//		}	
//	}
//	
//	/**************************************************************************
//	 * Simply rename the files, if there isn't more than one sample per lane
//	 **************************************************************************/
//	if (lanesToMove.size() > 0) {
//		for (String laneName: lanes.keySet()) {
//			TFSampleInfo tfsi = lanes.get(laneName).get(0);
//
//			this.moveFile(tfsi.getFileObject(TFConstants.FILE_LANE_BAM).getFinalPath(),tfsi.getFileObject(TFConstants.FILE_SPLIT_BAM).getFinalPath());
//			this.moveFile(tfsi.getFileObject(TFConstants.FILE_LANE_BAI).getFinalPath(),tfsi.getFileObject(TFConstants.FILE_SPLIT_BAI).getFinalPath());
//		}
//	}
//	
//	
//	/**************************************************************************
//	 * Process lanes
//	 **************************************************************************/
//	
//	if (lanesToRun.size() > 0) {
//			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(2).getName(),lanesToRun.size(),this.jobs);
//			this.daemon.start();
//			
//			int counter=1;
//		
//			for(String laneName: lanes.keySet()) {
//				//Add lanes for post-processing
//				lanesToPostProcess.add(laneName);
//				
//				ArrayList<TFSampleInfo> al = lanes.get(laneName);
//				
//				//Get representative sample
//				TFSampleInfo repSI = al.get(0);
//				
//				//Create FileObjects
//				TFFileObject tfoLaneBam = repSI.getFileObject(TFConstants.FILE_LANE_BAM);
//				TFFileObject tfoLaneBai = repSI.getFileObject(TFConstants.FILE_LANE_BAI);
//				TFFileObject tfoSplitBam = repSI.getFileObject(TFConstants.FILE_SPLIT_BAM);
//				
//				//Create run directory
//				File runDirectory =  tfoSplitBam.getWorkingDirectory();
//				if (runDirectory.exists()) {
//					this.deleteFolder(runDirectory);
//				}
//				runDirectory.mkdir();
//				
//				//Create local versions of input files
//				File fileLaneBam = tfoLaneBam.createDestForFileObject(runDirectory);
//				File fileLaneBai = tfoLaneBai.createDestForFileObject(runDirectory);
//				
//				//Create link
//				this.createLink(tfoLaneBam.getFinalPath(), fileLaneBam);
//				this.createLink(tfoLaneBai.getFinalPath(), fileLaneBai);
//				
//				//Create sample-specific storage
//				ArrayList<File> protectList = new ArrayList<File>();
//								
//				String indexing = "";
//				
//				for (TFSampleInfo tfsi: al) {
//					//runningListDirectory.add(tfsi);
//					indexing += " " + tfsi.getSampleName();
//				}
//				
//				//runningList.add(runningListDirectory);
//				indexing = indexing.substring(1);
//				
//				//create properties list
//				HashMap<String,String> replacements = new HashMap<String,String>();
//				replacements.put("INPUT_BAM", fileLaneBam.getName());
//				replacements.put("INPUT_LIST", indexing);
//				replacements.putAll(this.properties);
//				
//				//create cmd file
//				File cmdFile = new File(runDirectory,"cmd.txt");
//				this.createCmd(replacements,cmdFile,templateIdx);
//				
//			
//				//Protect files
//				protectList.add(cmdFile);
//				protectList.add(fileLaneBam);
//				protectList.add(fileLaneBai);
//				deleteList.add(fileLaneBam);
//				deleteList.add(fileLaneBai);
//				
//				//Submit job
//				
//				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
//				this.daemon.addJob(thread);
//				counter ++;
//						
//			}
//				
//			//Wait for command to finish
//			try {
//				this.daemon.join();
//				Thread.sleep(5000);
//				if (this.daemon.getFailed()) {
//					System.exit(1);
//				}
//			} catch (InterruptedException ie) {
//				logFile.writeErrorMessage("[TFExomeAlignUgp] Daemon interrupted",true);
//				System.exit(1);
//			}
//		
//	}
//	
//	/**************************************************************************
//	 * Post-Process lanes
//	 **************************************************************************/
//	
//	if (lanesToPostProcess.size() > 0) {
//		//Create final directories
//		this.finalDirectory.mkdir();
//		this.jobDirectory.mkdir();
//		
//		HashSet<File> runDirectoryList = new HashSet<File>();
//		
//		for (String laneName: lanesToPostProcess) {
//		
//			for (TFSampleInfo tfsi: lanes.get(laneName)) {
//				TFFileObject siSplitBam = tfsi.getFileObject(TFConstants.FILE_SPLIT_BAM);
//				TFFileObject siSplitBai = tfsi.getFileObject(TFConstants.FILE_SPLIT_BAI);
//				TFFileObject siLaneBam = tfsi.getFileObject(TFConstants.FILE_LANE_BAM);
//				TFFileObject siLaneBai = tfsi.getFileObject(TFConstants.FILE_LANE_BAI);
//				
//				File fileSplitBam = siSplitBam.createDestForFileObject(siSplitBam.getWorkingDirectory(),new String[]{tfsi.getSampleID() + "."},new String[]{tfsi.getSampleName() + "."});
//				File fileSplitBai = siSplitBai.createDestForFileObject(siSplitBam.getWorkingDirectory(),new String[]{tfsi.getSampleID() + "."},new String[]{tfsi.getSampleName() + "."});
//				
//				this.moveFile(fileSplitBam, siSplitBam.getFinalPath());
//				this.moveFile(fileSplitBai, siSplitBai.getFinalPath());
//				
//				this.deleteFile(siLaneBam.getFinalPath());
//				this.deleteFile(siLaneBai.getFinalPath());
//				
//				runDirectoryList.add(siSplitBam.getWorkingDirectory());
//			}
//		}
//		
//		this.cleanup(runDirectoryList, deleteList);
//	}
//
//}
	
	
}
