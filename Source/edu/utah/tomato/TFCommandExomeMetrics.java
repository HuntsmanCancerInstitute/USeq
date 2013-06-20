package edu.utah.tomato;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

public class TFCommandExomeMetrics extends TFCommand {

	public TFCommandExomeMetrics(File templateFile, File rootDirectory,
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
			File runDirectory = new File(this.rootDirectory,si.getSampleName() + "_metrics");
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
			
			//Create replacemnet hash
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("NAME", si.getSampleName());
			
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
		
		for (TFSampleInfo si: sampleList) {
			//Make destination directories
			File runDirectory = new File(this.rootDirectory,si.getSampleName() + "_metrics");
			File metricsDir = new File(this.rootDirectory,"metrics");
			metricsDir.mkdir();
			
			
			//Make destination files
			File metricsFile = new File(metricsDir,si.getSampleName() + ".pdf");
			
			
			//Move results
			this.moveFile(new File(runDirectory,si.getSampleName() + ".pdf"), metricsFile);

			
			//Add files to si object
			si.setFile(TFConstants.FILE_METRICS, metricsFile);
			
		}

	}
}



	



