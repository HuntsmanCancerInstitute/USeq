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
		int counter = 1;
		this.daemon = new TFThreadDaemon(this.logFile,this.commandString,sampleList.size(),this.jobs);
		this.daemon.start();
		for (TFSampleInfo si: sampleList) {

			//Create run directory
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleName() + "_metrics");
			runDirectory.mkdir();
			
			
			//Create command file
			File cmdFile = new File(runDirectory,"cmd.txt");
			
			//Create keepers files
			File bam = new File(runDirectory,si.getSampleName() + ".bam");
			File bai = new File(runDirectory,si.getSampleName() + ".bai");
			File sam = new File(runDirectory,si.getSampleName() + ".sam.gz");
			
			ArrayList<File> keepers = new ArrayList<File>();
			keepers.add(bam);
			keepers.add(sam);
			keepers.add(bai);
			keepers.add(cmdFile);
			
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
			keepers.add(destFile);
			
			replacements.put("NAME", si.getSampleName());
			replacements.put("TARGETS", destFile.getName());
			replacements.put("FILE",targetFile.getName());
			
			
			//Create cmd.txt file
			this.createCmd(replacements, cmdFile);
			
			TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, keepers, this.logFile);
			
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
		
		File metricsDir = new File(this.rootDirectory,"metrics");
		if (metricsDir.exists()) {
			IO.deleteDirectory(metricsDir);
		}
		File imageDir = new File(metricsDir,"images");
		metricsDir.mkdir();
		imageDir.mkdir();
		
		for (TFSampleInfo si: sampleList) {
			//Make destination directories
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleName() + "_metrics");
			
			File runImageDirectory = new File(runDirectory,"images");
			
			//Make destination files
			File metricsFile = new File(metricsDir,si.getSampleName() + ".html");
			
			
			//Move results
			this.moveFile(new File(runDirectory,si.getSampleName() + ".html"), metricsFile);
			this.cpFile(runImageDirectory, metricsDir);

			
			//Add files to si object
			si.setFile(TFConstants.FILE_METRICS, metricsFile);
		}
		
		mergeMetrics(metricsDir);

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



	



