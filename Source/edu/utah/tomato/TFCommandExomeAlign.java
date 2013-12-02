package edu.utah.tomato;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;


public class TFCommandExomeAlign extends TFCommand {

	
	public TFCommandExomeAlign(ArrayList<File> templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, boolean isFull, HashMap<String,String> properties) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress, isFull, properties);
		// TODO Auto-generated constructor stub
	}
	
	
	public void alignment(ArrayList<TFSampleInfo> sampleList) {
		/******************************************************************************************
		 *  Standard Alignment
		 *****************************************************************************************/
		logFile.writeInfoMessage("Starting sample alignment");
		
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
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleID() + "_align"); //tomato job directory
			runDirectoryList.add(runDirectory); //store directory information
			runDirectory.mkdir();
			
			//Create files
			File cmdFile = new File(runDirectory,"cmd.txt");
			File fastq1 = new File(runDirectory,si.getSampleID() + "_1.txt.gz");
			File fastq2 = new File(runDirectory,si.getSampleID() + "_2.txt.gz");
			
			//Link necessary files
			this.createLink(si.getFile(TFConstants.FILE_FASTQ1),fastq1);
			this.createLink(si.getFile(TFConstants.FILE_FASTQ2),fastq2);
			
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
			this.createCmd(replacements,cmdFile,0);
			
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
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		jobDir.mkdir();
		
		
		for (int i=0;i<sampleList.size();i++) {
			//Get sample-specific info
			File runDirectory = runDirectoryList.get(i);
			TFSampleInfo si = sampleList.get(i);
			
			//Make destination files
			File bamFile = new File(processedDir,si.getSampleID() + ".raw.bam");
			File baiFile = new File(processedDir,si.getSampleID() + ".raw.bai");
			
			//Move results
			this.moveFile(new File(runDirectory,si.getSampleID() + ".raw.bam"), bamFile);
			this.moveFile(new File(runDirectory,si.getSampleID() + ".raw.bai"), baiFile);
			
			//Add files to si object
			si.setFile(TFConstants.FILE_BAM, bamFile);
			si.setFile(TFConstants.FILE_BAI, baiFile);
		}
		
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
			
	}
	
	public void realignLane(ArrayList<TFSampleInfo> sampleList) {
		logFile.writeInfoMessage("Running realignment and recalibration on each lane");
		
		/**************************************************************************
		 * Group raw alignments by lane
		 **************************************************************************/
		
		//Group samples by Lane
		HashMap<String,ArrayList<TFSampleInfo>> lanes = new HashMap<String,ArrayList<TFSampleInfo>>();
		for (TFSampleInfo si: sampleList) {
			if (!lanes.containsKey(si.getPuID())) {
				lanes.put(si.getPuID(), new ArrayList<TFSampleInfo>());
			}
			lanes.get(si.getPuID()).add(si);
		}
		
		/**************************************************************************
		 *  Realign and recalibrate by lane
		 ************************************************************************/
		
		int counter = 1;
		this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(1).getName(),lanes.size(),this.jobs);
		this.daemon.start();
		
		//Containers for lane-level processing 
		ArrayList<File> runDirectoryList = new ArrayList<File>();
		ArrayList<File >deleteList = new ArrayList<File>();
		ArrayList<String> laneList = new ArrayList<String>();
		
		for(String laneName: lanes.keySet()) {
			
			//create run directory
			File runDirectory = new File(this.rootDirectory,"JOB_" + laneName + "_realign");
			runDirectoryList.add(runDirectory);
			runDirectory.mkdir();
			
			//Create lane-specific storage
			ArrayList<File> protectList = new ArrayList<File>(); //don't remove these files on job re-submission cleanup
			
			//Build dup input list and grab 
			String dupList = "";
			
			ArrayList<TFSampleInfo> laneBams = lanes.get(laneName);
			for (TFSampleInfo si: laneBams) {
				
				
				File sBam = new File(runDirectory,si.getSampleID() + ".raw.bam");
				File sBai = new File(runDirectory,si.getSampleID() + ".raw.bai");
				
				protectList.add(sBam);
				protectList.add(sBai);
				deleteList.add(sBam);
				deleteList.add(sBai);
				
				this.createLink(si.getFile(TFConstants.FILE_BAM),sBam);
				this.createLink(si.getFile(TFConstants.FILE_BAM),sBai);
				
				dupList += " INPUT=" + sBam.getName();
		
			}
			laneList.add(laneName);
			
			//create files
			File cmdFile = new File(runDirectory,"cmd.txt");
			protectList.add(cmdFile);
			
			//create properties list
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("NAME", laneName);
			replacements.put("DUP_LIST", dupList);
			replacements.putAll(this.properties);
			
			//Create cmd.txt file
			this.createCmd(replacements,cmdFile,1);
			
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
		
		File alignDir = new File(this.rootDirectory,"Alignments");
		File processedDir = new File(alignDir,"ProcessedAlignments");
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		jobDir.mkdir();
		
		for(int i=0;i<laneList.size();i++) {
			//Get lane-specific info
			File runDirectory = runDirectoryList.get(i);
			String laneName = laneList.get(i);
			
			//Make destination Files
			File bamFile = new File(processedDir,laneName + ".realign.bam");
			File baiFile = new File(processedDir,laneName + ".realign.bai");
			
			//Move Results
			this.moveFile(new File(runDirectory,laneName + ".realign.bam"), bamFile);
			this.moveFile(new File(runDirectory,laneName + ".realign.bai"), baiFile);
			
		}
		
		for (TFSampleInfo si: sampleList) {
			si.setFile(TFConstants.FILE_LANE_BAM, new File(processedDir,si.getPuID() + ".realign.bam"));
			si.setFile(TFConstants.FILE_LANE_BAI, new File(processedDir,si.getPuID() + ".realign.bai"));
		}
		
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
	}
		

	public void splitLaneBams(ArrayList<TFSampleInfo> sampleList) {
		
		
		//Group samples by Lane
		HashMap<String,ArrayList<TFSampleInfo>> lanes = new HashMap<String,ArrayList<TFSampleInfo>>();
		int numberJobs = 0;
		for (TFSampleInfo si: sampleList) {
			if (!lanes.containsKey(si.getPuID())) {
				lanes.put(si.getPuID(), new ArrayList<TFSampleInfo>());
			} 
			lanes.get(si.getPuID()).add(si);
		}
		
		for(String laneName: lanes.keySet()) {
			ArrayList<TFSampleInfo> al = lanes.get(laneName);
			if (al.size() > 1) {
				numberJobs += 1;
			}
		}
		
		if (numberJobs == 0) {
			logFile.writeInfoMessage("No lanes to split, moving on");
		} else {
			logFile.writeInfoMessage("Splitting lane bams into sample bams");
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(2).getName(),numberJobs,this.jobs);
			this.daemon.start();
		}
		
		
		/**************************************************************************
		 * Split alignments back into readgroups
		 ************************************************************************/
		
		int counter =1;
	
		
		//Containers for bam-splitting processing 
		ArrayList<File> runDirectoryList = new ArrayList<File>();
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<ArrayList<TFSampleInfo>> runningList = new ArrayList<ArrayList<TFSampleInfo>>();
		
		
		//Create necessary directories
		File alignDir = new File(this.rootDirectory,"Alignments");
		File processedDir = new File(alignDir,"ProcessedAlignments");
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		jobDir.mkdir();
		
		for(String laneName: lanes.keySet()) {
			ArrayList<TFSampleInfo> al = lanes.get(laneName);
			if (al.size() > 1) {
				TFSampleInfo si = al.get(0);
				
				ArrayList<File> protectList = new ArrayList<File>();
				
				File runDirectory = new File(this.rootDirectory,"JOB_" + laneName + "_split");
				runDirectoryList.add(runDirectory);
				ArrayList<TFSampleInfo> runningListDirectory = new ArrayList<TFSampleInfo>();
				
				runDirectory.mkdir();
				
				String indexing = "";
				
				for (TFSampleInfo tfsi: al) {
					runningListDirectory.add(tfsi);
					
					indexing += " " + tfsi.getSampleName();
				}
				
				runningList.add(runningListDirectory);
				indexing = indexing.substring(1);
				
				//create files
				File cmdFile = new File(runDirectory,"cmd.txt");
				protectList.add(cmdFile);
				
				File sBam = new File(runDirectory,si.getPuID() + ".realign.bam");
				File sBai = new File(runDirectory,si.getPuID() + ".realign.bai");
				
				
				protectList.add(sBam);
				protectList.add(sBai);
				deleteList.add(sBam);
				deleteList.add(sBai);
				
				this.createLink(si.getFile(TFConstants.FILE_LANE_BAM), sBam);
				this.createLink(si.getFile(TFConstants.FILE_LANE_BAI), sBai);
				
				//create properties list
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("INPUT_BAM", sBam.getName());
				replacements.put("INPUT_LIST", indexing);
				
				replacements.putAll(this.properties);
				this.createCmd(replacements,cmdFile,2);
				
				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
				
				this.daemon.addJob(thread);
				
				counter ++;
					
			
			} else {
				TFSampleInfo tfsi = al.get(0);
				File laneBam = new File(processedDir,tfsi.getSampleID() + ".split.bam");
				File laneBai = new File(processedDir,tfsi.getSampleID() + ".split.bai");
				
				tfsi.setFile(TFConstants.FILE_SPLIT_LANE_BAM, laneBam);
				tfsi.setFile(TFConstants.FILE_SPLIT_LANE_BAI, laneBai);
				
				this.moveFile(tfsi.getFile(TFConstants.FILE_LANE_BAM),laneBam);
				this.moveFile(tfsi.getFile(TFConstants.FILE_LANE_BAI),laneBai);
				
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
		
		
		for(int i=0;i<runDirectoryList.size();i++) {
			//Get lane-specific info
			File runDirectory = runDirectoryList.get(i);
			ArrayList<TFSampleInfo> tfsil = runningList.get(i);
			

			for (TFSampleInfo tfsi: tfsil) {
				//Make destination Files
				File bamFile = new File(processedDir,tfsi.getSampleID() + ".split.bam");
				File baiFile = new File(processedDir,tfsi.getSampleID() + ".split.bai");
				
				//Set sampleInfoObject
				tfsi.setFile(TFConstants.FILE_SPLIT_LANE_BAM, bamFile);
				tfsi.setFile(TFConstants.FILE_SPLIT_LANE_BAI, baiFile);
				
				//Move Results
				this.moveFile(new File(runDirectory,tfsi.getSampleName() + ".bam"), bamFile);
				this.moveFile(new File(runDirectory,tfsi.getSampleName() + ".bai"), baiFile);
			}
			
		}
		
		
		
		
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
		
		//Delete lane-level bams
		for (TFSampleInfo tfsi: sampleList) {
			if (tfsi.getFile(TFConstants.FILE_LANE_BAM).exists()) {
				tfsi.getFile(TFConstants.FILE_LANE_BAM).delete();
			}
			if (tfsi.getFile(TFConstants.FILE_LANE_BAI).exists()) {
				tfsi.getFile(TFConstants.FILE_LANE_BAI).delete();
			}
		}
	}
	
	public void mergeSampleLanes(ArrayList<TFSampleInfo> sampleList) {
		
		
		//Group samples by Lane
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		int numberJobs = 0;
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		for(String laneName: samples.keySet()) {
			ArrayList<TFSampleInfo> al = samples.get(laneName);
			if (al.size() > 1) {
				numberJobs += 1;
			}
		}
		

		if (numberJobs == 0) {
			logFile.writeInfoMessage("No samples to merge, moving on");
		} else {
			logFile.writeInfoMessage("Merging samples across lanes");
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(3).getName(),numberJobs,this.jobs);
			this.daemon.start();
		}
		
		/**************************************************************************
		 * Split alignments back into readgroups
		 ************************************************************************/
		
		int counter = 1;
		
		//Containers for bam-splitting processing 
		ArrayList<File> runDirectoryList = new ArrayList<File>();
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<String> runningNames = new ArrayList<String>();
		
		//Create necessary directories
		File alignDir = new File(this.rootDirectory,"Alignments");
		File processedDir = new File(alignDir,"ProcessedAlignments");
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		jobDir.mkdir();
		
		for(String sampleName: samples.keySet()) {
			ArrayList<TFSampleInfo> al = samples.get(sampleName);
			if (al.size() > 1) {
				//Create job directory
				File runDirectory = new File(this.rootDirectory,"JOB_" + al.get(0).getSampleName() + "_merge");
				runDirectoryList.add(runDirectory);
				runDirectory.mkdir();
				
				//Create run-specific protection
				ArrayList<File> protectList = new ArrayList<File>();
				
				//Add name
				runningNames.add(sampleName);
				
				String mergeList = "";
				
				for (TFSampleInfo si: samples.get(sampleName)) {
					File sBam = new File(runDirectory,si.getSampleID() + ".split.bam");
					File sBai = new File(runDirectory,si.getSampleID() + ".split.bai");

					deleteList.add(sBam);
					deleteList.add(sBai);
					protectList.add(sBam);
					protectList.add(sBai);
					
					this.createLink(si.getFile(TFConstants.FILE_SPLIT_LANE_BAM), sBam);
					this.createLink(si.getFile(TFConstants.FILE_SPLIT_LANE_BAI), sBai);
					
					mergeList += " INPUT=" + sBam.getName();
					
				}
					
				//create files
				File cmdFile = new File(runDirectory,"cmd.txt");
				protectList.add(cmdFile);
				
				String oBam = sampleName + ".sample.bam";
				
				
				//create properties list
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("INPUT_LIST", mergeList);
				replacements.put("OUTPUT", oBam);
				replacements.putAll(this.properties);
				this.createCmd(replacements,cmdFile,3);
				
				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
				
				this.daemon.addJob(thread);
				
				counter ++;
					
			} else {
				TFSampleInfo tfsi = al.get(0);
				File sampleBam = new File(processedDir,tfsi.getSampleName() + ".sample.bam");
				File sampleBai = new File(processedDir,tfsi.getSampleName() + ".sample.bai");
				
				tfsi.setFile(TFConstants.FILE_SAMPLE_BAM, sampleBam);
				tfsi.setFile(TFConstants.FILE_SAMPLE_BAI, sampleBai);
				
				this.cpFile(tfsi.getFile(TFConstants.FILE_SPLIT_LANE_BAM),sampleBam);
				this.cpFile(tfsi.getFile(TFConstants.FILE_SPLIT_LANE_BAI),sampleBai);
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
		
		
		for(int i=0;i<runningNames.size();i++) {
			//Get lane-specific info
			File runDirectory = runDirectoryList.get(i);
			String sampleName = runningNames.get(i);
			
			//Make destination Files
			File bamFile = new File(processedDir,sampleName + ".sample.bam");
			File baiFile = new File(processedDir,sampleName + ".sample.bai");
			
			//Set sampleInfoObject
			for (TFSampleInfo tfsi: sampleList) {
				if (tfsi.getSampleName().equals(sampleName)) {
					tfsi.setFile(TFConstants.FILE_SAMPLE_BAM, bamFile);
					tfsi.setFile(TFConstants.FILE_SAMPLE_BAI, baiFile);
				}
			}
			
			//Move Results
			this.moveFile(new File(runDirectory,sampleName + ".sample.bam"), bamFile);
			this.moveFile(new File(runDirectory,sampleName + ".sample.bai"), baiFile);
			
			
			
		}
		
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
	}
	
	public void realignSample(ArrayList<TFSampleInfo> sampleList) {
		
		//Group samples by Lane
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		int numberJobs = 0;
		
		for(String laneName: samples.keySet()) {
			ArrayList<TFSampleInfo> al = samples.get(laneName);
			if (al.size() > 1) {
				numberJobs += 1;
			}
		}
		
		
		if (numberJobs == 0) {
			logFile.writeInfoMessage("No merged samples to realign, moving on");
		} else {
			logFile.writeInfoMessage("Realigning merged samples");
			this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(4).getName(),numberJobs,this.jobs);
			this.daemon.start();
		}
		
		
		/**************************************************************************
		 * Split alignments back into readgroups
		 ************************************************************************/
		
		int counter = 1;
		
		
		//Containers for bam-splitting processing 
		ArrayList<File> runDirectoryList = new ArrayList<File>();
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<String> runningNames = new ArrayList<String>();
		
		//Create necessary directories
		File alignDir = new File(this.rootDirectory,"Alignments");
		File processedDir = new File(alignDir,"ProcessedAlignments");
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		jobDir.mkdir();
		
		for(String sampleName: samples.keySet()) {
			ArrayList<TFSampleInfo> al = samples.get(sampleName);
			if (al.size() > 1) {
				TFSampleInfo tfsi = al.get(0);
				
				//Create job directory
				File runDirectory = new File(this.rootDirectory,"JOB_" + tfsi.getSampleName() + "_realign");
				runDirectoryList.add(runDirectory);
				runDirectory.mkdir();
				
				//Create run-specific protection
				ArrayList<File> protectList = new ArrayList<File>();
				
				//Add name
				runningNames.add(sampleName);
				
				File sBam = new File(runDirectory,tfsi.getSampleName() + ".sample.bam");
				File sBai = new File(runDirectory,tfsi.getSampleName() + ".sample.bai");
				
				deleteList.add(sBam);
				deleteList.add(sBai);
				protectList.add(sBam);
				protectList.add(sBai);
				
				this.createLink(tfsi.getFile(TFConstants.FILE_SAMPLE_BAM), sBam);
				this.createLink(tfsi.getFile(TFConstants.FILE_SAMPLE_BAI), sBai);
						
				//create files
				File cmdFile = new File(runDirectory,"cmd.txt");
				protectList.add(cmdFile);
				
				//create properties list
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("NAME", sampleName);
				replacements.putAll(this.properties);
				this.createCmd(replacements,cmdFile,4);
				
				TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
				
				this.daemon.addJob(thread);
				
				counter ++;
						
			} else {
				TFSampleInfo tfsi = al.get(0);
				File realignBam = new File(processedDir,tfsi.getSampleName() + ".final.bam");
				File realignBai = new File(processedDir,tfsi.getSampleName() + ".final.bai");
				
				tfsi.setFile(TFConstants.FILE_REALIGN_SAMPLE_BAM, realignBam);
				tfsi.setFile(TFConstants.FILE_REALIGN_SAMPLE_BAI, realignBai);
				
				this.moveFile(tfsi.getFile(TFConstants.FILE_SAMPLE_BAM),realignBam);
				this.moveFile(tfsi.getFile(TFConstants.FILE_SAMPLE_BAI),realignBai);
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
		
		
		for(int i=0;i<runningNames.size();i++) {
			//Get lane-specific info
			File runDirectory = runDirectoryList.get(i);
			String sampleName = runningNames.get(i);
			
			//Make destination Files
			File bamFile = new File(processedDir,sampleName + ".final.bam");
			File baiFile = new File(processedDir,sampleName + ".final.bai");
			
			//Set sampleInfoObject
			for (TFSampleInfo tfsi: sampleList) {
				if (tfsi.getSampleName().equals(sampleName)) {
					tfsi.setFile(TFConstants.FILE_REALIGN_SAMPLE_BAM, bamFile);
					tfsi.setFile(TFConstants.FILE_REALIGN_SAMPLE_BAI, baiFile);
				}
			}
			
			//Move Results
			this.moveFile(new File(runDirectory,sampleName + ".final.bam"), bamFile);
			this.moveFile(new File(runDirectory,sampleName + ".final.bai"), baiFile);
			
		}
		
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
		
		//Delete lane-level bams
		for (TFSampleInfo tfsi: sampleList) {
			if (tfsi.getFile(TFConstants.FILE_SAMPLE_BAM).exists()) {
				tfsi.getFile(TFConstants.FILE_SAMPLE_BAM).delete();
			}
			if (tfsi.getFile(TFConstants.FILE_SAMPLE_BAI).exists()) {
				tfsi.getFile(TFConstants.FILE_SAMPLE_BAI).delete();
			}
		}
	}
	
	public void reduceSample(ArrayList<TFSampleInfo> sampleList) {
		logFile.writeInfoMessage("Reducing sample bams");
		//Group samples by Lane
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		int numberJobs = samples.size();

		
		/**************************************************************************
		 * Split alignments back into readgroups
		 ************************************************************************/
		
		int counter = 1;
		this.daemon = new TFThreadDaemon(this.logFile,this.templateFiles.get(5).getName(),numberJobs,this.jobs);
		this.daemon.start();
		
		//Containers for bam-splitting processing 
		ArrayList<File> runDirectoryList = new ArrayList<File>();
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<String> runningNames = new ArrayList<String>();
		
		//Create necessary directories
		File alignDir = new File(this.rootDirectory,"Alignments");
		File processedDir = new File(alignDir,"ProcessedAlignments");
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		jobDir.mkdir();
		
		for(String sampleName: samples.keySet()) {
			ArrayList<TFSampleInfo> al = samples.get(sampleName);
			TFSampleInfo tfsi = al.get(0);
			
			//Create job directory
			File runDirectory = new File(this.rootDirectory,"JOB_" + tfsi.getSampleName() + "_reduce");
			runDirectoryList.add(runDirectory);
			runDirectory.mkdir();
			
			//Create run-specific protection
			ArrayList<File> protectList = new ArrayList<File>();
			
			//Add name
			runningNames.add(sampleName);
			
			File sBam = new File(runDirectory,tfsi.getSampleName() + ".final.bam");
			File sBai = new File(runDirectory,tfsi.getSampleName() + ".final.bai");
			
			deleteList.add(sBam);
			deleteList.add(sBai);
			protectList.add(sBam);
			protectList.add(sBai);
			
			this.createLink(tfsi.getFile(TFConstants.FILE_SPLIT_LANE_BAM), sBam);
			this.createLink(tfsi.getFile(TFConstants.FILE_SPLIT_LANE_BAI), sBai);
					
			//create files
			File cmdFile = new File(runDirectory,"cmd.txt");
			protectList.add(cmdFile);
			
			//create properties list
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("NAME", sampleName);
			replacements.putAll(this.properties);
			this.createCmd(replacements,cmdFile,5);
			
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
			logFile.writeErrorMessage("Daemon interrupted",true);
			System.exit(1);
		}
		
		
		for(int i=0;i<runningNames.size();i++) {
			//Get lane-specific info
			File runDirectory = runDirectoryList.get(i);
			String sampleName = runningNames.get(i);
			
			//Make destination Files
			File bamFile = new File(processedDir,sampleName + ".reduce.bam");
			File baiFile = new File(processedDir,sampleName + ".reduce.bai");
			
			//Set sampleInfoObject
			for (TFSampleInfo tfsi: sampleList) {
				if (tfsi.getSampleName().equals(sampleName)) {
					tfsi.setFile(TFConstants.FILE_REDUCE_BAM, bamFile);
					tfsi.setFile(TFConstants.FILE_REDUCE_BAI, baiFile);
				}
			}
			
			//Move Results
			this.moveFile(new File(runDirectory,sampleName + ".reduce.bam"), bamFile);
			this.moveFile(new File(runDirectory,sampleName + ".reduce.bai"), baiFile);
			
		}
		
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
	}
	
	public void simpleAlign(ArrayList<TFSampleInfo> sampleList) {
		/******************************************************************************************
		 *  Standard Alignment
		 *****************************************************************************************/
		logFile.writeInfoMessage("Starting sample alignment");
		
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
			File runDirectory = new File(this.rootDirectory,"JOB_" + si.getSampleID() + "_align"); //tomato job directory
			runDirectoryList.add(runDirectory); //store directory information
			runDirectory.mkdir();
			
			//Create files
			File cmdFile = new File(runDirectory,"cmd.txt");
			File fastq1 = new File(runDirectory,si.getSampleID() + "_1.txt.gz");
			File fastq2 = new File(runDirectory,si.getSampleID() + "_2.txt.gz");
			
			//Link necessary files
			this.createLink(si.getFile(TFConstants.FILE_FASTQ1),fastq1);
			this.createLink(si.getFile(TFConstants.FILE_FASTQ2),fastq2);
			
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
			this.createCmd(replacements,cmdFile,0);
			
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
		File jobDir = new File(alignDir,"Jobs");
		alignDir.mkdir();
		processedDir.mkdir();
		jobDir.mkdir();
		
		
		for (int i=0;i<sampleList.size();i++) {
			//Get sample-specific info
			File runDirectory = runDirectoryList.get(i);
			TFSampleInfo si = sampleList.get(i);
			
			//Make destination files
			File rawBamFile = new File(processedDir,si.getSampleID() + ".raw.bam");
			File rawBaiFile = new File(processedDir,si.getSampleID() + ".raw.bai");
			File finalBamFile = new File(processedDir,si.getSampleID() + ".final.bam");
			File finalBaiFile = new File(processedDir,si.getSampleID() + ".final.bai");
			File reduceBamFile = new File(processedDir,si.getSampleID() + ".reduce.bam");
			File reduceBaiFile = new File(processedDir,si.getSampleID() + ".reduce.bai");
			
			//Move results
			this.moveFile(new File(runDirectory,si.getSampleID() + ".raw.bam"), rawBamFile);
			this.moveFile(new File(runDirectory,si.getSampleID() + ".raw.bai"), rawBaiFile);
			this.moveFile(new File(runDirectory,si.getSampleID() + ".final.bam"), finalBamFile);
			this.moveFile(new File(runDirectory,si.getSampleID() + ".final.bai"), finalBaiFile);
			this.moveFile(new File(runDirectory,si.getSampleID() + ".reduce.bam"), reduceBamFile);
			this.moveFile(new File(runDirectory,si.getSampleID() + ".reduce.bai"), reduceBaiFile);
			
			
			//Add files to si object
			si.setFile(TFConstants.FILE_BAM, rawBamFile);
			si.setFile(TFConstants.FILE_BAI, rawBaiFile);
			si.setFile(TFConstants.FILE_SPLIT_LANE_BAM, finalBamFile);
			si.setFile(TFConstants.FILE_SPLIT_LANE_BAI, finalBaiFile);
			si.setFile(TFConstants.FILE_REDUCE_BAM, reduceBamFile);
			si.setFile(TFConstants.FILE_REDUCE_BAI, reduceBaiFile);
		}
		
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
	}
	
	
	
	

	@Override
	public void run(ArrayList<TFSampleInfo> sampleList) {
		String commandName = (this.commandString.split(":"))[0];
		if (commandName.equals("exome_align_best")) {
			this.simpleAlign(sampleList);
		} else {
		
			this.alignment(sampleList);
			this.realignLane(sampleList);
			this.splitLaneBams(sampleList);
			this.mergeSampleLanes(sampleList);
			this.realignSample(sampleList);
			this.reduceSample(sampleList);
		}
		
		
		

	}
	
	public void check(ArrayList<TFSampleInfo> sampleList) {
		for (TFSampleInfo si: sampleList) {
			System.out.println("Alignment bam " + si.getFile(TFConstants.FILE_BAM));
			System.out.println("Lane realign bam " + si.getFile(TFConstants.FILE_LANE_BAM));
			System.out.println("Split bam " + si.getFile(TFConstants.FILE_SPLIT_LANE_BAM));
			System.out.println("Merge bam " + si.getFile(TFConstants.FILE_SAMPLE_BAM));
			System.out.println("Sample realign bam " + si.getFile(TFConstants.FILE_REALIGN_SAMPLE_BAM));
			System.out.println("reduce bam " + si.getFile(TFConstants.FILE_REDUCE_BAM));
		}
	}
	
	
	
}
