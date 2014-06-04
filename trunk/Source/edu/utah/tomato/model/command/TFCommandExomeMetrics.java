package edu.utah.tomato.model.command;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
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




public class TFCommandExomeMetrics extends TFCommand {
	private File targetFile;
	private String studyName;
	private boolean deleteMetricsBam;
	
	
	public TFCommandExomeMetrics(ArrayList<File> templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, boolean deleteMetricsBam, boolean isFull, String studyName, File targetFile, HashMap<String,String> properties) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress, isFull, properties);
		this.targetFile = targetFile;
		this.studyName = studyName;
		this.deleteMetricsBam = deleteMetricsBam;
	
		this.finalDirectory = new File(this.rootDirectory,"Metrics");
		this.jobDirectory = new File(this.finalDirectory,"Jobs");
		// TODO Auto-generated constructor stub
	}
	
	@Override
	protected ArrayList<TFSampleInfo> findPrereqsExisting(ArrayList<TFSampleInfo> sampleList) {
		//Create patterns of interest
		ArrayList<TFMatchObject> dependantPatterns = new ArrayList<TFMatchObject>();
		
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_BAM,Pattern.compile("(.+?)(\\.raw\\.bam$)"),TFConstants.PREFIX_SAMPLEID));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_BAI,Pattern.compile("(.+?)(\\.raw\\.bai$)"),TFConstants.PREFIX_SAMPLEID));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_SPLIT_BAM,Pattern.compile("(.+?)(\\.split.bam$)"),TFConstants.PREFIX_SAMPLEID));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_SPLIT_BAI,Pattern.compile("(.+?)(\\.split.bai$)"),TFConstants.PREFIX_SAMPLEID));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAM,Pattern.compile("(.+?)(\\.final\\.bam$)"),TFConstants.PREFIX_SAMPLENAME));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAI,Pattern.compile("(.+?)(\\.final\\.bai$)"),TFConstants.PREFIX_SAMPLENAME));
		
		//Grab all files in the ProcessedAlignments Directory.
		File rootAlignments = new File(this.rootDirectory,"Alignments");
		File processedAlignments = new File(rootAlignments,"ProcessedAlignments");
		if (processedAlignments.exists()) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsExisting(sampleList, processedAlignments, dependantPatterns);
			
			if (foundSampleList.size() > 0) {
				if (foundSampleList.size() == sampleList.size()) {
					sampleList = foundSampleList;
					logFile.writeInfoMessage(String.format("[TFExomeMetrics] Found %d of %d samples in the ProcessedAlignments directory.",foundSampleList.size(),sampleList.size()));
				} else {
					logFile.writeErrorMessage(String.format("[TFExomeMetrics] Found fewer than expected samples (%d of %d) in the ProcessedAlignments directory, exiting",foundSampleList.size(),sampleList.size()), false);
					System.exit(1);
				}
			} else {
				logFile.writeErrorMessage("[TFExomeMetrics] Did not find any potential samples in the ProcessedAlignments directory, exiting",false);
				System.exit(1);
			}
		} else {
			logFile.writeErrorMessage("[TFExomeMetrics] Could not find ProcessedAlignments directory, exiting",false);
			System.exit(1);
		}
		return sampleList;
	}
	
	@Override
	protected ArrayList<TFSampleInfo> findPrereqsNew(ArrayList<TFSampleInfo> sampleList) {
		//Create patterns of interest
		ArrayList<TFMatchObject> masterPatterns = new ArrayList<TFMatchObject>();
		ArrayList<TFMatchObject> dependantPatterns = new ArrayList<TFMatchObject>();
		
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_BAM,Pattern.compile("(.+?)(\\.raw\\.bam$)"),TFConstants.PREFIX_SAMPLEID));
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_BAI,Pattern.compile("(.+?)(\\.raw\\.bai$)"),TFConstants.PREFIX_SAMPLEID));

		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_SPLIT_BAM,Pattern.compile("(.+?)(\\.split.bam$)"),TFConstants.PREFIX_SAMPLEID));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_SPLIT_BAI,Pattern.compile("(.+?)(\\.split.bai$)"),TFConstants.PREFIX_SAMPLEID));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAM,Pattern.compile("(.+?)(\\.final\\.bam$)"),TFConstants.PREFIX_SAMPLENAME));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_FINAL_BAI,Pattern.compile("(.+?)(\\.final\\.bai$)"),TFConstants.PREFIX_SAMPLENAME));
		
		//Grab all files in the ProcessedAlignments Directory.
		File rootAlignments = new File(this.rootDirectory,"Alignments");
		File processedAlignments = new File(rootAlignments,"ProcessedAlignments");
		if (processedAlignments.exists()) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsNew(processedAlignments, masterPatterns, dependantPatterns);
			
			if (foundSampleList.size() > 0) {
				sampleList = foundSampleList;
				logFile.writeInfoMessage(String.format("[TFExomeMetrics] Found %d samples in the ProcessedAlignments directory.",foundSampleList.size()));
			} else {
				logFile.writeInfoMessage("[TFExomeMetrics] Did not find any potential samples in the ProcessedAlignments directory, checking run directory");
			}
		}
		
		if (sampleList.size() == 0) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsNew(this.rootDirectory, masterPatterns, dependantPatterns);
			if (foundSampleList.size() > 0) {
				sampleList = foundSampleList;
				logFile.writeInfoMessage(String.format("[TFExomeMetrics] Found %d samples in the run directory directory.",foundSampleList));
			} else {
				logFile.writeErrorMessage("[TFExomeMetrics] Did not find any potential samples in the run directory, exiting",false);
				System.exit(1);
			}
		}
		
		return sampleList;
	}
	
	/**This method goes through each sample object and checks if necessary prerequisites were found.
	 * Note that a whole set WON'T fail if one fails.  This method only checks one directory at time, 
	 * which means that samples can't be spread across two locations (i.e. final results directory 
	 * and run directory), which seems reasonable to me.  
	 *
	 * If everything is found, an ArrayList of starting files are returned.  The starting
	 * files might be .split.bam or .final.bam, depending if samples were split across lanes.  Note
	 * there can be a mix
	 * 
	 */
	@Override
	protected ArrayList<TFSampleInfo> validateSampleSet(ArrayList<TFSampleInfo> sampleList) {
		ArrayList<TFSampleInfo> validSamples = new ArrayList<TFSampleInfo>();
	
		boolean valid = true;
		
		for (TFSampleInfo tfsi: sampleList) {				
			if (tfsi.finalFileExists(TFConstants.FILE_SPLIT_BAM) && tfsi.finalFileExists(TFConstants.FILE_SPLIT_BAI)) {
				//OK!
			} else if (tfsi.finalFileExists(TFConstants.FILE_FINAL_BAM) && tfsi.finalFileExists(TFConstants.FILE_FINAL_BAI)) {
				tfsi.setFileObject(TFConstants.FILE_SPLIT_BAM, tfsi.getFileObject(TFConstants.FILE_FINAL_BAM));
				tfsi.setFileObject(TFConstants.FILE_SPLIT_BAI, tfsi.getFileObject(TFConstants.FILE_FINAL_BAI));
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
	public ArrayList<TFSampleInfo> run(ArrayList<TFSampleInfo> sampleList) {
		TFThread.setFailCount(0);
		
		if (sampleList.size() == 0) {
			sampleList = this.findPrereqsNew(sampleList);
		} else {
			sampleList = this.findPrereqsExisting(sampleList);
		}
		
		this.generateMetrics(sampleList);
		
//		if (this.deleteMetricsBam) {
//			for (TFSampleInfo tfsi: sampleList) {
//				if (tfsi.fileExists(TFConstants.FILE_FINAL_BAM)) {
//					if (tfsi.getFileObject(TFConstants.FILE_SPLIT_BAM).getFinalPath() != tfsi.getFileObject(TFConstants.FILE_FINAL_BAM).getFinalPath()) {
//						tfsi.getFileObject(TFConstants.FILE_SPLIT_BAM).getFinalPath().delete();
//						tfsi.getFileObject(TFConstants.FILE_SPLIT_BAI).getFinalPath().delete();
//					}
//				} else {
//					tfsi.getFileObject(TFConstants.FILE_SPLIT_BAM).getFinalPath().delete();
//					tfsi.getFileObject(TFConstants.FILE_SPLIT_BAI).getFinalPath().delete();
//				}
//			}
//		}
		
		return sampleList;
	}
	
	
	protected void mergeMetrics(File metricsDir) {
		try {
			
			ProcessBuilder pb = new ProcessBuilder("java","-Xmx10g","-jar",properties.get("USEQ_PATH_LOCAL") + "/Apps/MergeExonMetrics","-f",metricsDir.getAbsolutePath(),"-o",studyName);
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[TFExomeMetrics] Could not merge metrics files",true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[TFExomeMetrics] IO Exception while trying to merge metrics",true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[TFExomeMetrics] Process was interrupted while trying to merge metrics",true);
			System.exit(1);
		}
	}
	
	private void generateMetrics(ArrayList<TFSampleInfo> sampleList) {
		ArrayList<TFSampleInfo> samplesToRun = new ArrayList<TFSampleInfo>();
		ArrayList<TFSampleInfo> samplesToPostProcess = new ArrayList<TFSampleInfo>();
		HashSet<File> deleteList = new HashSet<File>();
		
		//Check sample stage
		for (TFSampleInfo si: sampleList) {
			File workingDir = new File(this.rootDirectory,"JOB_" + si.getSampleID() + "_metrics");
			
			//Create output objects for each sample
			TFFileObject tfoMetricsDict = new TFFileObject(si.getSampleID() + ".dict.txt",this.finalDirectory,workingDir);
			TFFileObject tfoFinalMetrics = new TFFileObject(this.studyName + ".xlsx",this.finalDirectory,this.finalDirectory);
			si.setFileObject(TFConstants.FILE_METRICS, tfoMetricsDict);
			
			if (!tfoFinalMetrics.doesFinalExist()) { //Final output exists, do nothing
				if (!tfoMetricsDict.doesWorkingExist()) { 
					samplesToRun.add(si);
				} else {
					samplesToPostProcess.add(si);
				}
			}
		}
		
		//Process samples
		if (samplesToRun.size() > 0) {
			//Initialize daemon-specific variables
			int counter = 1;
			this.daemon = new TFThreadDaemon(this.logFile,this.commandString,samplesToRun.size(),this.jobs);
			this.daemon.start();
			
	
			for (TFSampleInfo si: samplesToRun) {
				samplesToPostProcess.add(si);
				
				//Create sample-specific storage
				ArrayList<File> protectList = new ArrayList<File>();
				
				TFFileObject tfoMetricsDict = si.getFileObject(TFConstants.FILE_METRICS);

				//Create run directory
				File runDirectory = tfoMetricsDict.getWorkingDirectory();
				runDirectory.mkdir();
				
				
				//Get input objects
				TFFileObject tfoRawBam = si.getFileObject(TFConstants.FILE_BAM);
				TFFileObject tfoRawBai = si.getFileObject(TFConstants.FILE_BAI);
				TFFileObject tfoFinalBam = si.getFileObject(TFConstants.FILE_SPLIT_BAM);
				TFFileObject tfoFinalBai = si.getFileObject(TFConstants.FILE_SPLIT_BAI);
				
				//Create local versions of the input files
				File fileRawBam = tfoRawBam.createDestForFileObject(runDirectory);
				File fileRawBai = tfoRawBai.createDestForFileObject(runDirectory);
				
				File fileFinalBam = tfoFinalBam.createDestForFileObject(runDirectory,new String[]{"split",si.getSampleName() + "."},new String[]{"final",si.getSampleID() + "."});
				File fileFinalBai = tfoFinalBai.createDestForFileObject(runDirectory,new String[]{"split",si.getSampleName() + "."},new String[]{"final",si.getSampleID() + "."});
				
				//Create links
				this.createLink(tfoRawBam.getFinalPath(), fileRawBam);
				this.createLink(tfoRawBai.getFinalPath(), fileRawBai);
				this.createLink(tfoFinalBam.getFinalPath(), fileFinalBam);
				this.createLink(tfoFinalBai.getFinalPath(), fileFinalBai);
				
				//Create replacement hash
				HashMap<String,String> replacements = new HashMap<String,String>();
				
				//Copy target region to directory
				if (this.targetFile == null) {
					this.targetFile = new File(properties.get("TARGET_DEFAULT"));
				}
				
				File destFile = new File(runDirectory,targetFile.getName());
				this.cpFile(targetFile, destFile);
		
				
				replacements.put("NAME", si.getSampleID());
				replacements.put("TARGETS", destFile.getName());
				replacements.put("FILE",targetFile.getName());
				replacements.putAll(this.properties);
				
				
				//Create cmd.txt file
				File cmdFile = new File(runDirectory,"cmd.txt");
				this.createCmd(replacements, cmdFile,0);
				
			
				//Mark files for deletion of cleanup protection
				protectList.add(fileRawBam);
				protectList.add(fileRawBai);
				protectList.add(fileFinalBam);
				protectList.add(fileFinalBai);
				protectList.add(cmdFile);
				protectList.add(destFile);
				deleteList.add(fileRawBam);
				deleteList.add(fileRawBai);
				deleteList.add(fileFinalBam);
				deleteList.add(fileFinalBai);
							
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
				logFile.writeErrorMessage("[TFExomeMetrics] Daemon interrupted",true);
				System.exit(1);
			}		
		}
		
		
		if (samplesToPostProcess.size() > 0) {
			//Create final directories
			this.finalDirectory.mkdirs();
			this.jobDirectory.mkdir();
			HashSet<File> runDirectoryList = new HashSet<File>();
			
			for (TFSampleInfo tfsi: samplesToPostProcess) {
				TFFileObject tfoMetricsDict = tfsi.getFileObject(TFConstants.FILE_METRICS);
				
				File workingDir = tfoMetricsDict.getWorkingDirectory();
				runDirectoryList.add(workingDir);
				
				File fileMetricsDict = tfoMetricsDict.createDestForFileObject(workingDir);
				File fileImageDir = new File(workingDir,"images");
				
				this.moveFile(fileMetricsDict, tfoMetricsDict.getFinalPath());
				this.cpFile(fileImageDir,this.finalDirectory);
			
			}
			
			//Merge metrics data
			mergeMetrics(this.finalDirectory);
			this.moveFile(new File(studyName + ".xlsx"), this.finalDirectory);
			
			
			//cleanup
			this.cleanup(runDirectoryList, deleteList);
		}
	}
	
}



	



