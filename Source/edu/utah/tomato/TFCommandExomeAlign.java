package edu.utah.tomato;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

public class TFCommandExomeAlign extends TFCommand {

	
	public TFCommandExomeAlign(File templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void run(ArrayList<TFSampleInfo> sampleList) {
		int counter = 1;
		this.daemon = new TFThreadDaemon(this.logFile,this.commandString,sampleList.size(),this.jobs);
		this.daemon.start();
		for (TFSampleInfo si: sampleList) {
		
			//Create run directory
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleName() + "_align");
			runDirectory.mkdir();
			
			//Create command file
			File cmdFile = new File(runDirectory,"cmd.txt");
			
			//Create output files
			File fastq1 = new File(runDirectory,si.getSampleName() + "_1.txt.gz");
			File fastq2 = new File(runDirectory,si.getSampleName() + "_2.txt.gz");
			ArrayList<File> keepers = new ArrayList<File>();
			keepers.add(fastq1);
			keepers.add(fastq2);
			keepers.add(cmdFile);
			
			//Link necessary files
			this.createLink(si.getFile(TFConstants.FILE_FASTQ1),fastq1);
			this.createLink(si.getFile(TFConstants.FILE_FASTQ2),fastq2);
			
			//Create replacement hash
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("NAME", si.getSampleName());
			if (si.isQual64()) {
				replacements.put("QUAL_FLAG", "-I");
			} else {
				replacements.put("QUAL_FLAG","");
			}
			
			//Create cmd.txt file
			this.createCmd(replacements,cmdFile);
			
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

	

}
