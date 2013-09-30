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
	
	public TFCommandExomeMetrics(File templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, String studyName, File targetFile) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress);
		this.targetFile = targetFile;
		this.studyName = studyName;
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
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleName() + "_metrics");
			runDirectoryList.add(runDirectory);
			runDirectory.mkdir();
			
			//Create  files
			File bam = new File(runDirectory,si.getSampleName() + ".bam");
			File bai = new File(runDirectory,si.getSampleName() + ".bai");
			File sam = new File(runDirectory,si.getSampleName() + ".sam.gz");
			File cmdFile = new File(runDirectory,"cmd.txt");
			
			//Link necessary files
			this.createLink(si.getFile(TFConstants.FILE_BAM),bam);
			this.createLink(si.getFile(TFConstants.FILE_BAI),bai);
			this.createLink(si.getFile(TFConstants.FILE_SAM),sam);
			
			//Create replacement hash
			HashMap<String,String> replacements = new HashMap<String,String>();
			
			//Copy target region to directory
			if (this.targetFile == null) {
				this.targetFile = new File("/home/u0855942/tim_genomes/hg19.ccds.20130301.bed");
			}
			File destFile = new File(runDirectory,targetFile.getName());
			this.cpFile(targetFile, destFile);
			
			
			replacements.put("NAME", si.getSampleName());
			replacements.put("TARGETS", destFile.getName());
			replacements.put("FILE",targetFile.getName());
			
			
			//Create cmd.txt file
			this.createCmd(replacements, cmdFile);
			
		
			//Mark files for deletion of cleanup protection
			protectList.add(bam);
			protectList.add(sam);
			protectList.add(bai);
			protectList.add(cmdFile);
			protectList.add(destFile);
			deleteList.add(bam);
			deleteList.add(bai);
			deleteList.add(sam);
						
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
			File metricsFile = new File(metricsDir,si.getSampleName() + ".dict.txt");
			
			
			//Move results
			this.moveFile(new File(runDirectory,si.getSampleName() + ".dict.txt"), metricsFile);
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
			this.moveFile(rd, jobDir);
		}
		

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



	



