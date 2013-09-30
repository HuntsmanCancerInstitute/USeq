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
		
		//Create project-specific storage
		ArrayList<File> runDirectoryList = new ArrayList<File>(); //List of the run directories
		ArrayList<File> deleteList = new ArrayList<File>(); //list of files to delete at the end of the run
		
		for (TFSampleInfo si: sampleList) {
			//Create sample-specific storage
			ArrayList<File> protectList = new ArrayList<File>(); //don't remove these files on job re-submission cleanup
		
			//Create run directory
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleName() + "_align"); //tomato job directory
			runDirectoryList.add(runDirectory); //store directory information
			runDirectory.mkdir();
			
			//Create files
			File cmdFile = new File(runDirectory,"cmd.txt");
			File fastq1 = new File(runDirectory,si.getSampleName() + "_1.txt.gz");
			File fastq2 = new File(runDirectory,si.getSampleName() + "_2.txt.gz");
			
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
			logFile.writeErrorMessage("Daemon interrupted",true);
			System.exit(1);
		}
		
		//Create final directories
		File alignDir = new File(this.rootDirectory,"Alignments");
		File processedDir = new File(alignDir,"ProcessedAlignments");
		File rawDir = new File(alignDir,"RawAlignments");
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		rawDir.mkdir();
		jobDir.mkdir();
		
		
		for (int i=0;i<sampleList.size();i++) {
			//Get sample-specific info
			File runDirectory = runDirectoryList.get(i);
			TFSampleInfo si = sampleList.get(i);
			
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
		
		//clean up unwanted files
		for (File df: deleteList) {
			this.deleteFile(df);
		}
		
		//Move JOB directories
		for (File rd: runDirectoryList) {
			this.moveFile(rd, jobDir);
		}
		
		//Fin

	}

	

}
