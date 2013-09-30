package edu.utah.tomato;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

public class TFCommandBisulfite extends TFCommand {
	
	private String genome;

	public TFCommandBisulfite(File templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress) {
		super(templateFile, rootDirectory, commandString, commandType, logFile,
				email, wallTime, heartbeat, failmax, jobs, suppress);
		
	}

	@Override
	public void run(ArrayList<TFSampleInfo> sampleList) {
		int counter = 1;
		this.daemon = new TFThreadDaemon(this.logFile,this.commandString,sampleList.size(),this.jobs);
		this.daemon.start();
		
		
		for (TFSampleInfo si: sampleList) {
			//set split number
			int splitNumber = 100000000;

			//Create run directory
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleName() + "_bisulfite");
			runDirectory.mkdir();
			
			//Create output files
			File fastq1 = new File(runDirectory,si.getSampleName() + "_1.txt.gz");
			File fastq2 = new File(runDirectory,si.getSampleName() + "_2.txt.gz");
			
			//Link necessary files
			this.createLink(si.getFile(TFConstants.FILE_FASTQ1),fastq1);
			this.createLink(si.getFile(TFConstants.FILE_FASTQ2),fastq2);
			
			//Calculate split file count
			int recordCount = si.getRecordCount() * 4;
			int fileCount = recordCount / splitNumber + 1;
			
			//split files
			this.fileSplitter(fastq1, splitNumber);
			this.fileSplitter(fastq2, splitNumber);
			
			//Create run directories and run stuff
			for (int i=1;i<=fileCount;i++) {
				File dir = new File(runDirectory,String.valueOf(i));
				File sf1 = new File(runDirectory,String.valueOf(i) + "_" + si.getSampleName() + "_1.txt.gz");
				File sf2 = new File(runDirectory,String.valueOf(i) + "_" + si.getSampleName() + "_2.txt.gz");
				
				this.moveFile(sf1, dir);
				this.moveFile(sf2, dir);
				
				File cmdFile = new File(dir,"cmd.txt");
				
				//store files
				ArrayList<File> keepers = new ArrayList<File>();
				keepers.add(sf1);
				keepers.add(sf2);
				keepers.add(cmdFile);
				
				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, keepers, this.logFile);
				
				//this.taskList.add(thread);
				this.daemon.addJob(thread);
				
				counter ++;
				
				
			}
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
		
		for (TFSampleInfo si: sampleList) {
			//Make destination directories
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleName() + "_align");
			File processedDir = new File(this.rootDirectory,"processed_alignments");
			File rawDir = new File(this.rootDirectory,"raw_alignments");
			processedDir.mkdir();
			rawDir.mkdir();
			
			//Make destination files
			File bamFile = new File(processedDir,si.getSampleName() + ".bam");
			File baiFile = new File(processedDir,si.getSampleName() + ".bai");
			File samFile = new File(rawDir,si.getSampleName() + ".sam.gz");
			
			//Move results
			this.moveFile(new File(runDirectory,si.getSampleName() + ".bam"), bamFile);
			this.moveFile(new File(runDirectory,si.getSampleName() + ".bai"), baiFile);
			this.moveFile(new File(runDirectory,si.getSampleName() + ".sam.gz"), samFile);
			
			//Add files to si object
			si.setFile(TFConstants.FILE_BAM, bamFile);
			si.setFile(TFConstants.FILE_BAI, baiFile);
			si.setFile(TFConstants.FILE_SAM, samFile);
		}

	}
	
	private void fileSplitter(File source, Integer count) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[fileSplitter] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			ProcessBuilder pb = new ProcessBuilder("java ","-Xmx10g","-jar","/home/BioApps/USeq/Apps/FileSplitter",
					"-f",source.getAbsolutePath(),"-n",String.valueOf(count),"-g");

			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[fileSplitter] System could not split your file " + source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[fileSplitter] IO Exception while trying to split your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[fileSplitter] Process was interrupted while trying to split your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}

}
