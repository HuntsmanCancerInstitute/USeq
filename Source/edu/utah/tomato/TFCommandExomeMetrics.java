package edu.utah.tomato;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;


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
		// TODO Auto-generated constructor stub
	}

	@Override
	public void run(ArrayList<TFSampleInfo> sampleList) {
		//Initialize daemon-specific variables
		int counter = 1;
		this.daemon = new TFThreadDaemon(this.logFile,this.commandString,sampleList.size(),this.jobs);
		this.daemon.start();
		
		//Create project-specific storage
		ArrayList<File> runDirectoryList = new ArrayList<File>();
		ArrayList<File> deleteList = new ArrayList<File>();
		
		
		for (TFSampleInfo si: sampleList) {
			//Create sample-specific storage
			ArrayList<File> protectList = new ArrayList<File>();

			//Create run directory
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleID() + "_metrics");
			runDirectoryList.add(runDirectory);
			runDirectory.mkdir();
			
			//Create  files
			File mbam = new File(runDirectory,si.getSampleID() + ".mate.bam");
			File mbai = new File(runDirectory,si.getSampleID() + ".mate.bai");
			File rbam = new File(runDirectory,si.getSampleID() + ".split.lane.bam");
			File rbai = new File(runDirectory,si.getSampleID() + ".split.lane.bai");
			
			File cmdFile = new File(runDirectory,"cmd.txt");
			
			//Link necessary files
			this.createLink(si.getFile(TFConstants.FILE_BAM),mbam);
			this.createLink(si.getFile(TFConstants.FILE_BAI),mbai);
			this.createLink(si.getFile(TFConstants.FILE_SPLIT_LANE_BAM),rbam);
			this.createLink(si.getFile(TFConstants.FILE_SPLIT_LANE_BAI),rbai);
			
			//Create replacement hash
			HashMap<String,String> replacements = new HashMap<String,String>();
			
			//Copy target region to directory
			if (this.targetFile == null) {
				this.targetFile = new File("/home/u0855942/tim_genomes/hg19.ccds.20131019.bed");
			}
			
			File destFile = new File(runDirectory,targetFile.getName());
			this.cpFile(targetFile, destFile);
	
			
			replacements.put("NAME", si.getSampleID());
			replacements.put("TARGETS", destFile.getName());
			replacements.put("FILE",targetFile.getName());
			replacements.putAll(this.properties);
			
			
			//Create cmd.txt file
			this.createCmd(replacements, cmdFile,0);
			
		
			//Mark files for deletion of cleanup protection
			protectList.add(rbam);
			protectList.add(rbai);
			protectList.add(mbam);
			protectList.add(mbai);
			protectList.add(cmdFile);
			protectList.add(destFile);
			deleteList.add(mbam);
			deleteList.add(mbai);
			deleteList.add(rbam);
			deleteList.add(rbai);
						
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
			logFile.writeErrorMessage("Daemon interrupted",true);
			System.exit(1);
		}
		
		//Create output files
		File metricsDir = new File(this.rootDirectory,"Metrics");
		File jobDir = new File(metricsDir,"Jobs");
		if (metricsDir.exists()) {
			IO.deleteDirectory(metricsDir);
		}
		File imageDir = new File(metricsDir,"images");
		metricsDir.mkdir();
		imageDir.mkdir();
		jobDir.mkdir();
		
		for (int i=0; i<sampleList.size();i++) {
			//get sample-specific info
			File runDirectory = runDirectoryList.get(i);
			File runImageDirectory = new File(runDirectory,"images");
			TFSampleInfo si = sampleList.get(i);
			

			//Make destination files
			File metricsFile = new File(metricsDir,si.getSampleID() + ".dict.txt");
			
			
			//Move results
			this.moveFile(new File(runDirectory,si.getSampleID() + ".dict.txt"), metricsFile);
			this.cpFile(runImageDirectory, metricsDir);

			
			//Add files to si object
			si.setFile(TFConstants.FILE_METRICS, metricsFile);
		}
		
		//Merge metrics data
		mergeMetrics(metricsDir);
		
		this.moveFile(new File(studyName + ".xlsx"), metricsDir);
		
		//clean up unwanted files
		for (File df: deleteList) {
			this.deleteFile(df);
		}
		
		//Move JOB directories
		for (File rd: runDirectoryList) {
			File existDir = new File(jobDir,rd.getName());
			if (existDir.exists()) {
				deleteFolder(existDir);
			}
			this.moveFile(rd, jobDir);
		}
		
		//Clean up metrics files if specified
		this.logFile.writeInfoMessage("Deleting raw alignments");
		if (this.deleteMetricsBam && this.isFull) {
			for (TFSampleInfo si: sampleList) {
				if (si.getFile(TFConstants.FILE_BAM).exists()) {
					si.getFile(TFConstants.FILE_BAM).delete();
				}
				if (si.getFile(TFConstants.FILE_BAI).exists()) {
					si.getFile(TFConstants.FILE_BAI).delete();
				}
				if (si.getFile(TFConstants.FILE_SPLIT_LANE_BAM).exists()) {
					si.getFile(TFConstants.FILE_SPLIT_LANE_BAM).delete();
				}
				if (si.getFile(TFConstants.FILE_SPLIT_LANE_BAI).exists()) {
					si.getFile(TFConstants.FILE_SPLIT_LANE_BAI).delete();
				}
			}
		}
		
        //FIN
	}
	
	
	protected void mergeMetrics(File metricsDir) {
		try {
			
			ProcessBuilder pb = new ProcessBuilder("java","-Xmx2g","-jar","/home/u0855942/tim_scripts/USeq_vcf/Apps/MergeExonMetrics","-f",metricsDir.getAbsolutePath(),"-o",studyName);
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[mergeMetrics] Could not merge metrics files",true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[mergeMetrics] IO Exception while trying to merge metrics",true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[mergeMetrics] Process was interrupted while trying to merge metrics",true);
			System.exit(1);
		}
	}
	
}



	



