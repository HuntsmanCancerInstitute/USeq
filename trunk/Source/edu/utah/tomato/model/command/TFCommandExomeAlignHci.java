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

public class TFCommandExomeAlignHci extends TFCommand {
	private boolean validateFastq = false;
	
	public TFCommandExomeAlignHci(ArrayList<File> templateFile, File rootDirectory,
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
			logFile.writeErrorMessage("[TFExomeAlignHci] Command modules upstream of this one aren't currently supported.",true);
			System.exit(1);
		} else {
			sampleList = this.findPrereqsNew(sampleList);
		}
		
		if (this.validateFastq) {
			for (TFSampleInfo fi: sampleList) {
				validateFileSet(fi);
			}
		}
		
		this.align(sampleList, 0);
		this.mergeSampleLanes(sampleList, 1);
		this.reduceSample(sampleList, 2);
		
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
			logFile.writeInfoMessage(String.format("[TFExomeAlignHci] Found %d samples in the run directory.",foundSampleList.size()));
		} else {
			logFile.writeErrorMessage("[TFExomeAlignHci] Did not find any matching samples in the run directory.",false);
			System.exit(1);
		}
		return sampleList;
	}
	
	
	
	private void align(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		/******************************************************************************************
		 *  HCI Alignment parameters
		 *****************************************************************************************/
		logFile.writeInfoMessage("[TFExomeAlignHci] Starting sample alignment");
		
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<TFSampleInfo> samplesToRun = new HashSet<TFSampleInfo>();
		HashSet<TFSampleInfo> samplesToPostProcess = new HashSet<TFSampleInfo>();
		
		/******************************************************************************************
		 *  Check sample stage
		 *****************************************************************************************/
		
		for (TFSampleInfo tfsi: sampleList) {
			if (!tfsi.getFileObject(TFConstants.FILE_FASTQ1).doesFinalExist() || !tfsi.getFileObject(TFConstants.FILE_FASTQ2).doesFinalExist()) {
				logFile.writeErrorMessage("[TFExomeAlignHci] Fastq files could not be found by TF,  exiting",true);
				System.exit(1);
			}
			
			File runDirectory = new File(this.rootDirectory,"JOB_" + tfsi.getSampleID() + "_align");
			File mergeDirectory = new File(this.rootDirectory,"JOB_" + tfsi.getSampleName() + "_merge");
			
			//Create file objects
			TFFileObject tfoIdFinalBam = new TFFileObject(tfsi.getSampleID() + ".final.bam",this.finalDirectory,runDirectory);
			TFFileObject tfoIdFinalBai = new TFFileObject(tfsi.getSampleID() + ".final.bai",this.finalDirectory,runDirectory);
			TFFileObject tfoRawBam = new TFFileObject(tfsi.getSampleID() + ".raw.bam",this.finalDirectory,runDirectory);
			TFFileObject tfoRawBai = new TFFileObject(tfsi.getSampleID() + ".raw.bai",this.finalDirectory,runDirectory);
			TFFileObject tfoFinalBam = new TFFileObject(tfsi.getSampleName() + ".final.bam",this.finalDirectory,mergeDirectory);
			TFFileObject tfoFinalBai = new TFFileObject(tfsi.getSampleName() + ".final.bai",this.finalDirectory,mergeDirectory);
			tfsi.setFileObject(TFConstants.FILE_FINAL_BAM, tfoFinalBam);
			tfsi.setFileObject(TFConstants.FILE_FINAL_BAI, tfoFinalBai);
			tfsi.setFileObject(TFConstants.FILE_BAM, tfoRawBam);
			tfsi.setFileObject(TFConstants.FILE_BAI, tfoRawBai);
			tfsi.setFileObject(TFConstants.FILE_ID_FINAL_BAM,tfoIdFinalBam);
			tfsi.setFileObject(TFConstants.FILE_ID_FINAL_BAI,tfoIdFinalBai);
			
			//(found && found) && ((found && found) || (missing && missing))
			
			if (tfoIdFinalBam.doesFinalExist() && tfoIdFinalBai.doesFinalExist() && 
					tfoRawBam.doesFinalExist() && tfoRawBai.doesFinalExist()) {
			} else if (tfoFinalBam.doesFinalExist() && tfoFinalBai.doesFinalExist() && 
					tfoRawBam.doesFinalExist() && tfoRawBai.doesFinalExist()) {
				tfsi.setFileObject(TFConstants.FILE_ID_FINAL_BAM, tfoFinalBam);
				tfsi.setFileObject(TFConstants.FILE_ID_FINAL_BAI, tfoFinalBai);
			} else if ((tfoRawBam.doesWorkingExist() && tfoRawBai.doesWorkingExist()) && 
						((tfoIdFinalBam.doesWorkingExist() && tfoIdFinalBai.doesWorkingExist()))) {
				samplesToPostProcess.add(tfsi);
			} else {
				samplesToRun.add(tfsi);
			}
		}
		
		/******************************************************************************************
		 *  Process samples
		 *****************************************************************************************/
		
		if (samplesToRun.size() > 0) {
			int counter = 1;
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(templateIdx).getName(),samplesToRun.size(),this.jobs);
			this.daemon.start();
			
			
			for (TFSampleInfo si: samplesToRun) {
				samplesToPostProcess.add(si);
				
				//Create sample-specific storage
				ArrayList<File> protectList = new ArrayList<File>(); //don't remove these files on job re-submission cleanup
			
				TFFileObject tfoFastq1 = si.getFileObject(TFConstants.FILE_FASTQ1);
				TFFileObject tfoFastq2 = si.getFileObject(TFConstants.FILE_FASTQ2);
				TFFileObject tfoRawBam = si.getFileObject(TFConstants.FILE_BAM);
				
				//Get working
				File workingDir = tfoRawBam.getWorkingDirectory();
				workingDir.mkdir();
				
				File fileFastq1 = tfoFastq1.createDestForFileObject(workingDir,new String[]{"fastq"},new String[]{"txt"});
				File fileFastq2 = tfoFastq2.createDestForFileObject(workingDir,new String[]{"fastq"},new String[]{"txt"});
				
				this.createLink(tfoFastq1.getFinalPath(), fileFastq1);
				this.createLink(tfoFastq2.getFinalPath(), fileFastq2);
				
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
				

				//Create files
				File cmdFile = new File(workingDir,"cmd.txt");
				this.createCmd(replacements,cmdFile,templateIdx);
				
				//Mark files for deletion or cleanup protection
				protectList.add(fileFastq1);
				protectList.add(fileFastq2);
				protectList.add(cmdFile);
				deleteList.add(fileFastq1);
				deleteList.add(fileFastq2);
			
				//Run job
				TFThread thread = new TFThread(workingDir,this.failmax, counter, this.heartbeat, protectList, this.logFile);
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
				logFile.writeErrorMessage("[TFExomeAlignHci] Daemon interrupted",true);
				System.exit(1);
			}
			
		}
		
		/******************************************************************************************
		 *  Post-Process samples
		 *****************************************************************************************/
		
		if (samplesToPostProcess.size() > 0) {
			
			this.finalDirectory.mkdirs();
			this.jobDirectory.mkdir();
			
			HashSet<File> workingDirectoryList = new HashSet<File>();
			for (TFSampleInfo tfsi: samplesToPostProcess) {
				//Output file objects
				TFFileObject tfoRawBam = tfsi.getFileObject(TFConstants.FILE_BAM);
				TFFileObject tfoRawBai = tfsi.getFileObject(TFConstants.FILE_BAI);
				TFFileObject tfoIdFinalBam = tfsi.getFileObject(TFConstants.FILE_ID_FINAL_BAM);
				TFFileObject tfoIdFinalBai = tfsi.getFileObject(TFConstants.FILE_ID_FINAL_BAI);
				
				workingDirectoryList.add(tfoIdFinalBam.getWorkingDirectory());
				
			    //move files
				this.moveFile(tfoRawBam.getWorkingPath(),tfoRawBam.getFinalPath());
				this.moveFile(tfoRawBai.getWorkingPath(),tfoRawBai.getFinalPath());
				this.moveFile(tfoIdFinalBam.getWorkingPath(),tfoIdFinalBam.getFinalPath());
				this.moveFile(tfoIdFinalBai.getWorkingPath(),tfoIdFinalBai.getFinalPath());
			}
			
			this.cleanup(workingDirectoryList, deleteList);
		}
		
		
		
	}
	
	private void mergeSampleLanes(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
		//This method expects FILE_SPLIT_BAM
		//This method creates FILE_SAMPLE_BAM
		
		/**************************************************************************
		 *  Merge Sample Bams across lanes
		 ************************************************************************/
		
		logFile.writeInfoMessage("[TFExomeAlignHci] Merging samples across lanes");
		
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
			boolean skipFile = true;
			ArrayList<TFSampleInfo> al = samples.get(sampleName);
			if (al.size() > 1) {
				skipFile = false;
			}
			
			for (TFSampleInfo si: samples.get(sampleName)) {
				if (!si.getFileObject(TFConstants.FILE_ID_FINAL_BAM).doesFinalExist() || !si.getFileObject(TFConstants.FILE_ID_FINAL_BAI).doesFinalExist()) {
					logFile.writeErrorMessage("[TFExomeAlignHci] Final/Reduced bam/bai files cannot be found, exiting", true);
					System.exit(1);
				}
				
				//Create output objects for each sample
				TFFileObject tfoFinalBam = si.getFileObject(TFConstants.FILE_FINAL_BAM);
				TFFileObject tfoFinalBai = si.getFileObject(TFConstants.FILE_FINAL_BAI);
	
				
				
				//Check for the existence of output files
				if (tfoFinalBam.doesFinalExist() && tfoFinalBai.doesFinalExist() ) {
					//Skip
				} else if (skipFile) {
					samplesToRename.add(sampleName);
				} else if (!tfoFinalBam.doesWorkingExist() || !tfoFinalBai.doesWorkingExist()) {
					samplesToRun.add(sampleName);
				} else {
					samplesToPostProcess.add(sampleName);
				}
			}
		}
		
		/**************************************************************************
		 * Rename bam file
		 **************************************************************************/
		
		if (samplesToRename.size() > 0) {
			for (String sampleName: samplesToRename) {
				TFSampleInfo repSi = samples.get(sampleName).get(0);
				
				this.moveFile(repSi.getFileObject(TFConstants.FILE_ID_FINAL_BAM).getFinalPath(), repSi.getFileObject(TFConstants.FILE_FINAL_BAM).getFinalPath());
				this.moveFile(repSi.getFileObject(TFConstants.FILE_ID_FINAL_BAI).getFinalPath(), repSi.getFileObject(TFConstants.FILE_FINAL_BAI).getFinalPath());
				
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
				
				TFFileObject tfoFinalSampleBam = repSI.getFileObject(TFConstants.FILE_FINAL_BAM);
				
				
				//Create run directory
				File runDirectory = tfoFinalSampleBam.getWorkingDirectory();
				runDirectory.mkdir();
				
				ArrayList<File> protectList = new ArrayList<File>();
				String mergeList = "";
				
				
				//Prep working directory.
				for (TFSampleInfo tfsi: samples.get(sampleName)) {
					//Create input objects
					TFFileObject tfoFinalBam = tfsi.getFileObject(TFConstants.FILE_ID_FINAL_BAM);
					TFFileObject tfoFinalBai = tfsi.getFileObject(TFConstants.FILE_ID_FINAL_BAI);
					
					
					//Create local versions of input files
					File fileFinalBam = tfoFinalBam.createDestForFileObject(runDirectory);
					File fileFinalBai = tfoFinalBai.createDestForFileObject(runDirectory);
					
					
					deleteList.add(fileFinalBam);
					deleteList.add(fileFinalBai);
					protectList.add(fileFinalBam);
					protectList.add(fileFinalBai);
					
					
					mergeList += " INPUT=" + fileFinalBam.getName();
					
					//Create links
					this.createLink(tfoFinalBam.getFinalPath(), fileFinalBam);
					this.createLink(tfoFinalBai.getFinalPath(), fileFinalBai);

				}
				
				//create properties list	
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("INPUT_LIST", mergeList);
				replacements.put("OUTPUT", tfoFinalSampleBam.getFileName());

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
				logFile.writeErrorMessage("Daemon interrupted",true);
				System.exit(1);
			}
		}
		
		
		/**************************************************************************
		 * Post-Process
		 ************************************************************************/
		
		if (samplesToPostProcess.size() > 0) {
			this.finalDirectory.mkdirs();
			this.jobDirectory.mkdir();
			
			HashSet<File> runDirectoryList = new HashSet<File>();
			for (String sampleName: samplesToPostProcess) {				
				//Get representative si
				TFSampleInfo repSI = samples.get(sampleName).get(0);

				TFFileObject tfoFinalBam = repSI.getFileObject(TFConstants.FILE_FINAL_BAM);
				TFFileObject tfoFinalBai = repSI.getFileObject(TFConstants.FILE_FINAL_BAI);
				
				runDirectoryList.add(tfoFinalBam.getWorkingDirectory());
				
				this.moveFile(tfoFinalBam.getWorkingPath(), tfoFinalBam.getFinalPath());
				this.moveFile(tfoFinalBai.getWorkingPath(), tfoFinalBai.getFinalPath());
				
				for (TFSampleInfo tfsi: samples.get(sampleName)) {
					//Detele lane-level bam, replace pointer with merged
					this.deleteFile(tfsi.getFileObject(TFConstants.FILE_ID_FINAL_BAM).getFinalPath());
					this.deleteFile(tfsi.getFileObject(TFConstants.FILE_ID_FINAL_BAI).getFinalPath());
				}
			}
			this.cleanup(runDirectoryList, deleteList);
		}
	}
	
	
	
	private void reduceSample(ArrayList<TFSampleInfo> sampleList, int templateIdx) {
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
				if (!tfsi.getFileObject(TFConstants.FILE_FINAL_BAM).doesFinalExist() || !tfsi.getFileObject(TFConstants.FILE_FINAL_BAI).doesFinalExist()) {
					logFile.writeErrorMessage("[TFExomeAlignUgp] Final bam files weren't found by TF", true);
					System.exit(1);
				}
				
				//Create output objects for each sample
				TFFileObject tfoReduceBam = new TFFileObject(tfsi.getSampleName() + ".reduce.bam",this.finalDirectory,workingDir);
				TFFileObject tfoReduceBai = new TFFileObject(tfsi.getSampleName() + ".reduce.bai",this.finalDirectory,workingDir);
				tfsi.setFileObject(TFConstants.FILE_REDUCE_BAM, tfoReduceBam);
				tfsi.setFileObject(TFConstants.FILE_REDUCE_BAI, tfoReduceBai);
				
				//Check for existance of output files
				if (!tfoReduceBam.doesFinalExist() || !tfoReduceBai.doesFinalExist()) {
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
//			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(templateIdx).getName(),samplesToRun.size(),this.jobs);
//			this.daemon.start();
//			int counter = 1;
//			
//			for(String sampleName: samples.keySet()) {
//				//Get representative sample
//				TFSampleInfo siRep = samples.get(sampleName).get(0);
//				
//				//Create input/output objecst
//				TFFileObject tfoReduceBam = siRep.getFileObject(TFConstants.FILE_REDUCE_BAM);
//				TFFileObject tfoFinalBam = siRep.getFileObject(TFConstants.FILE_FINAL_BAM);
//				TFFileObject tfoFinalBai = siRep.getFileObject(TFConstants.FILE_FINAL_BAI);
//				
//				//Create job directory
//				File runDirectory = tfoReduceBam.getWorkingDirectory();
//				if (runDirectory.exists()) {
//					this.deleteFolder(runDirectory);
//				}
//				runDirectory.mkdir();
//				samplesToPostProcess.add(sampleName);
//				
//				//Create files
//				File fileFinalBam = tfoFinalBam.createDestForFileObject(runDirectory);
//				File fileFinalBai = tfoFinalBai.createDestForFileObject(runDirectory);
//				
//				//Create links
//				this.createLink(tfoFinalBam.getFinalPath(),fileFinalBam);
//				this.createLink(tfoFinalBai.getFinalPath(),fileFinalBai);
//				
//				//Create run-specific protection
//				ArrayList<File> protectList = new ArrayList<File>();
//				
//				//create properties list
//				HashMap<String,String> replacements = new HashMap<String,String>();
//				replacements.put("NAME", sampleName);
//				replacements.putAll(this.properties);
//				
//				//create files
//				File cmdFile = new File(runDirectory,"cmd.txt");
//				this.createCmd(replacements,cmdFile,templateIdx);
//				
//				deleteList.add(fileFinalBam);
//				deleteList.add(fileFinalBai);
//				protectList.add(fileFinalBam);
//				protectList.add(fileFinalBai);
//				protectList.add(cmdFile);
//				
//				//Run job
//				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
//				this.daemon.addJob(thread);
//				counter ++;		
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
			
			for(String sampleName: samples.keySet()) {
				//Get representative sample
				TFSampleInfo siRep = samples.get(sampleName).get(0);
				
				//Create input/output objecst
				TFFileObject tfoFinalBam = siRep.getFileObject(TFConstants.FILE_FINAL_BAM);
				TFFileObject tfoFinalBai = siRep.getFileObject(TFConstants.FILE_FINAL_BAI);
				TFFileObject tfoReduceBam = siRep.getFileObject(TFConstants.FILE_REDUCE_BAM);
				TFFileObject tfoReduceBai = siRep.getFileObject(TFConstants.FILE_REDUCE_BAI);
				
				File runDirectory = tfoReduceBam.getWorkingDirectory();
				if (runDirectory.exists()) {
					this.deleteFolder(runDirectory);
				}
				runDirectory.mkdir();
				samplesToPostProcess.add(sampleName);
			
				this.cpFile(tfoFinalBam.getFinalPath(), tfoReduceBam.getWorkingPath());
				this.cpFile(tfoFinalBai.getFinalPath(), tfoReduceBai.getWorkingPath());
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

	
}

