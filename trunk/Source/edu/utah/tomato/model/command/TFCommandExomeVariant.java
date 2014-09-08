package edu.utah.tomato.model.command;

import java.io.BufferedInputStream;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.tomato.daemon.TFThread;
import edu.utah.tomato.daemon.TFThreadDaemon;
import edu.utah.tomato.model.TFCommand;
import edu.utah.tomato.model.TFFileObject;
import edu.utah.tomato.model.TFSampleInfo;
import edu.utah.tomato.util.TFConstants;
import edu.utah.tomato.util.TFLogger;
import edu.utah.tomato.model.TFMatchObject;





public class TFCommandExomeVariant extends TFCommand {
	
	private String splitType = "none";
	private String study = null;
	//private String[] hg19Chrom = {"chr1","chr2","chr3","chr4","chr5","chr6",
	//		"chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
	//		"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
	//		"chrX","chrY","chrM"};
	private String[] hg19Chrom = {"1","2","3","4","5","6","7","8","9","10",
			"11","12","13","14","15","16","17","18","19","20","21","22",
			"X","Y","MT"};
	private File targetFile;
	private double coveragePerChunk = 1000000000d;
	private boolean use1KGenomes;
	private boolean deleteReducedBams;
	private boolean variantBest = false;
	
	public TFCommandExomeVariant(ArrayList<File> templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, boolean deleteReducedBam, boolean isFull, boolean use1KGenomes,
			String study, String splitType, File targetFile,
			HashMap<String,String> properties) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress, isFull, properties);
		this.splitType = splitType;
		this.study = study;
		this.targetFile = targetFile;
		this.use1KGenomes = use1KGenomes;
		this.deleteReducedBams = deleteReducedBam;	
		this.finalDirectory = new File(this.rootDirectory,"Variants");
		this.jobDirectory = new File(this.finalDirectory,"Jobs");
		
		//exome variant best does not have 
		if (this.commandType.equals("exome_variant_best")) {
			this.variantBest = true;
		}
	}
	
	@Override
	public ArrayList<TFSampleInfo> run(ArrayList<TFSampleInfo> sampleList) {
		TFThread.setFailCount(0);
		
		//Check for samples.
		if (sampleList.size() == 0) {
			sampleList = this.findPrereqsNew(sampleList);
		} else {
			sampleList = this.findPrereqsExisting(sampleList);
		}
		
	
		//Use the final list object.
		if (this.use1KGenomes) {
			logFile.writeInfoMessage("[TFExomeVariant] Attempting to load 1K genomes samples");
			ArrayList<TFSampleInfo> onekList = get1KGenomes();
			if (onekList == null) {
				logFile.writeWarningMessage("[TFExomeVariant] You don't have the permissions to access the 1K genomes archive, moving on without them");
			} else {
				sampleList.addAll(onekList);
			}
		}
		
		File runDirectory = new File(this.rootDirectory,"JOB_" + this.study + "_variant");
		File mergedDirectory = new File(runDirectory,"merged_calls");
		boolean isFinished = false;
		
		//Check for the existence of the output file
		if (!this.variantBest) {
			TFFileObject tfoRawVcf = new TFFileObject(this.study + ".raw.vcf.gz",this.finalDirectory,mergedDirectory);
			TFFileObject tfoRawVcfIdx = new TFFileObject(this.study + ".raw.vcf.gz.tbi",this.finalDirectory,mergedDirectory);
			TFFileObject tfoFilterVcf = new TFFileObject(this.study + ".filterFieldSetAll.vcf.gz",this.finalDirectory,mergedDirectory);
			TFFileObject tfoFilterVcfIdx = new TFFileObject(this.study + ".filterFieldSetAll.vcf.gz.tbi",this.finalDirectory,mergedDirectory);
			TFFileObject tfoPassingVcf = new TFFileObject(this.study + ".filterFieldSetPassing.vcf.gz",this.finalDirectory,mergedDirectory);
			TFFileObject tfoPassingVcfIdx = new TFFileObject(this.study + ".filterFieldSetPassing.vcf.gz.tbi",this.finalDirectory,mergedDirectory);
			
			if (tfoRawVcf.doesFinalExist() && tfoRawVcfIdx.doesFinalExist() 
					&& tfoFilterVcf.doesFinalExist() && tfoFilterVcfIdx.doesFinalExist() 
					&& tfoPassingVcf.doesFinalExist() && tfoPassingVcfIdx.doesFinalExist() ) {
				isFinished = true;
			} else {
				for (TFSampleInfo tfsi: sampleList) {
					tfsi.setFileObject(TFConstants.FILE_VCF_RAW, tfoRawVcf);
					tfsi.setFileObject(TFConstants.FILE_VCF_RAW_IDX,tfoRawVcfIdx);
					tfsi.setFileObject(TFConstants.FILE_VCF_FILTER,tfoFilterVcf);
					tfsi.setFileObject(TFConstants.FILE_VCF_FILTER_IDX, tfoFilterVcfIdx);
					tfsi.setFileObject(TFConstants.FILE_VCF_PASSING,tfoPassingVcf);
					tfsi.setFileObject(TFConstants.FILE_VCF_PASSING_IDX, tfoPassingVcfIdx);
				}
			}
		} else {
			TFFileObject tfoRawVcf = new TFFileObject(this.study + ".raw.vcf.gz",this.finalDirectory,mergedDirectory);
			TFFileObject tfoRawVcfIdx = new TFFileObject(this.study + ".raw.vcf.gz.tbi",this.finalDirectory,mergedDirectory);
			
			if (tfoRawVcf.doesFinalExist() && tfoRawVcfIdx.doesFinalExist()) {
				isFinished = true;
			} else {
				for (TFSampleInfo tfsi: sampleList) {
					tfsi.setFileObject(TFConstants.FILE_VCF_RAW, tfoRawVcf);
					tfsi.setFileObject(TFConstants.FILE_VCF_RAW_IDX,tfoRawVcfIdx);
				}
			}
		}
		
		if (!isFinished) {
			if (!runDirectory.exists()) {
				runDirectory.mkdir();
			}
			
			//Create Links to bam files
			for (TFSampleInfo tfsi: sampleList) {
				TFFileObject tfoReduceBam = tfsi.getFileObject(TFConstants.FILE_REDUCE_BAM);
				TFFileObject tfoReduceBai = tfsi.getFileObject(TFConstants.FILE_REDUCE_BAI);
				
				File fileReduceBam = tfoReduceBam.createDestForFileObject(runDirectory);
				File fileReduceBai = tfoReduceBai.createDestForFileObject(runDirectory);
				
				this.createLink(tfoReduceBam.getFinalPath(), fileReduceBam);
				this.createLink(tfoReduceBai.getFinalPath(), fileReduceBai);
			}
			
			if (this.targetFile == null) {
				logFile.writeInfoMessage("[TFExomeVariant] Calling variations across all regions");
			} else {
				logFile.writeInfoMessage("[TFExomeVariant] Calling variations using: " + targetFile.getAbsolutePath());
			}
			
			if (splitType.equals("chrom")) {
				runByChrom(sampleList,runDirectory);
			} else if (splitType.equals("chunk")) {
				runByChunk(sampleList,runDirectory);
			} else if (splitType.equals("none")) {
				runFull(sampleList,runDirectory);
			} else {
				logFile.writeErrorMessage("[TFExomeVariant] Don't recognize run type: " + splitType, true);
				System.exit(1);
			}
			
			//Delete job container directory
			if (runDirectory.exists()) {
				this.deleteFolder(runDirectory);
			}
		}
		
		
//		//Clean up reduced files if specified
//		if (this.deleteReducedBams && this.isFull) {
//			this.logFile.writeInfoMessage("Deleting reduced alignments");
//			for (TFSampleInfo si: sampleList) {
//				this.deleteFile(si.getFileObject(TFConstants.FILE_REDUCE_BAM).getFinalPath());
//				this.deleteFile(si.getFileObject(TFConstants.FILE_REDUCE_BAI).getFinalPath());
//			}
//		}
		
		return sampleList;
	}
	
	@Override
	protected ArrayList<TFSampleInfo> validateSampleSet(ArrayList<TFSampleInfo> sampleList) {
		ArrayList<TFSampleInfo> validSamples = new ArrayList<TFSampleInfo>();
		boolean valid = true;
		for (TFSampleInfo tfsi: sampleList) {
			if (tfsi.finalFileExists(TFConstants.FILE_REDUCE_BAM) && tfsi.finalFileExists(TFConstants.FILE_REDUCE_BAI)) {
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
		//Create patterns of interest
		ArrayList<TFMatchObject> dependantPatterns = new ArrayList<TFMatchObject>();
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_REDUCE_BAM,Pattern.compile("(.+?)\\.reduce\\.bam$"),TFConstants.PREFIX_SAMPLENAME));
		dependantPatterns.add(new TFMatchObject(TFConstants.FILE_REDUCE_BAI,Pattern.compile("(.+?)\\.reduce\\.bai$"),TFConstants.PREFIX_SAMPLENAME));
		
		File rootAlignments = new File(this.rootDirectory,"Alignments");
		File processedAlignments = new File(rootAlignments,"ProcessedAlignments");
		if (processedAlignments.exists()) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsExisting(sampleList,processedAlignments, dependantPatterns);
			if (foundSampleList.size() > 0) {
				if (foundSampleList.size() == sampleList.size()) {
					sampleList = foundSampleList;
					logFile.writeInfoMessage(String.format("[TFExomeVariant] Found %d of %d samples in the ProcessedAlignments directory",foundSampleList.size(),sampleList.size()));
				} else {
					logFile.writeErrorMessage(String.format("[TFExomeVariant] Found fewer than expected samples (%d of %d) in the ProcessedAlignments directory, exiting",foundSampleList.size(),sampleList.size()), false);
					System.exit(1);
				} 
			} else {
				logFile.writeErrorMessage("[TFExomeVariant] Did not find any potential samples in the ProcessedAlignments directory, exiting",false);
				System.exit(1);
			}
		} else {
			logFile.writeErrorMessage("[TFExomeVariant] Could not find the ProcessedAlignments directory, exiting",true);
			System.exit(1);
		}
		
		return sampleList;
	}

	@Override
	protected ArrayList<TFSampleInfo> findPrereqsNew(ArrayList<TFSampleInfo> sampleList) {
		ArrayList<TFMatchObject> masterPatterns = new ArrayList<TFMatchObject>();
		
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_REDUCE_BAM,Pattern.compile("(.+?)\\.reduce\\.bam$"),TFConstants.PREFIX_SAMPLENAME));
		masterPatterns.add(new TFMatchObject(TFConstants.FILE_REDUCE_BAI,Pattern.compile("(.+?)\\.reduce\\.bai$"),TFConstants.PREFIX_SAMPLENAME));
		
		File rootAlignments = new File(this.rootDirectory,"Alignments");
		File processedAlignments = new File(rootAlignments,"ProcessedAlignments");
		if (processedAlignments.exists()) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsNew(processedAlignments, masterPatterns, new ArrayList<TFMatchObject>());
			
			if (foundSampleList.size() > 0) {
				sampleList = foundSampleList;
				logFile.writeInfoMessage(String.format("[TFExomeVariant] Found %d samples in the ProcessedAlignments directory.",foundSampleList.size()));
			} else {
				logFile.writeInfoMessage("[TFExomeVariant] Did not find any potential samples in the ProcessedAlignments directory, checking run directory");
			}
		}
		
		if (sampleList.size() == 0) {
			ArrayList<TFSampleInfo> foundSampleList = this.findPatternsNew(this.rootDirectory, masterPatterns, new ArrayList<TFMatchObject>());
			if (foundSampleList.size() > 0) {
				sampleList = foundSampleList;
				logFile.writeInfoMessage(String.format("[TFExomeVariant] Found %d samples in the run directory directory.",foundSampleList.size()));
			} else {
				logFile.writeErrorMessage("[TFExomeVariant] Did not find any potential samples in the run directory, exiting",false);
				System.exit(1);
			}
		}
		
		return sampleList;
	}
	
	

	
	private void runByChrom(ArrayList<TFSampleInfo> sampleList, File mainJobDir) {
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<File> runDirectoryList = new HashSet<File>();
		ArrayList<File> fullVcfList = new ArrayList<File>();
		
		//Determine which chromosomes need to be run
		ArrayList<String> chromsToRun = new ArrayList<String>();
		
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
	
		for (String chrom: this.hg19Chrom) {
			File chromRunDir = new File(mainJobDir,"JOB_" + chrom + "_variant");
			File chromVcf = new File(chromRunDir,this.study + ".raw.vcf.gz");
			if (!chromVcf.exists()) {
				chromsToRun.add(chrom);
				runDirectoryList.add(chromRunDir);
			} else {
				fullVcfList.add(chromVcf);
			}
		}
		
		
		if (chromsToRun.size() > 0) {
			this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(0).getName(), chromsToRun.size(), this.jobs);
			this.daemon.start();
			int counter = 1;
			
			ArrayList<BedEntry> bedFile = this.getTargetRegions();
			
			for (String chrom: chromsToRun) {
				
				//Check if there is space on the cluster for another job
				boolean firstPass = true;
				while(true) {
					if (firstPass && this.daemon.isMaxed()) {
						logFile.writeInfoMessage("[TFExomeVariant] Pausing splitting, waiting for servers to catch up");
						firstPass = false;
					} else if (!this.daemon.isMaxed()) {
						logFile.writeInfoMessage("[TFExomeVariant] Resuming splitting");
						break;
					}
					
					try {
						Thread.sleep(5000);
					} catch (InterruptedException e) {
						logFile.writeErrorMessage("[TFExomeVariant] Interupted while waiting to resume bam splitting, exiting",true);
						System.exit(1);
					}
				}
				
				//Create chromosome specific directory
				File chromRunDir = new File(mainJobDir,"JOB_" + chrom + "_variant");
				if (chromRunDir.exists()) {
					this.deleteFolder(chromRunDir);
				}
				chromRunDir.mkdir();
				runDirectoryList.add(chromRunDir);
				
				//Create interval file
				File chromIntervals = new File(chromRunDir,"interval.bed");
				ArrayList<BedEntry> intervalBed = new ArrayList<BedEntry>();
				for (BedEntry interval: bedFile) {
					if (interval.getChrom().equals(chrom)) {
						intervalBed.add(interval);
					}
				}
				
				if (intervalBed.size() == 0) {
					logFile.writeInfoMessage("[TFExomeVariant] There are no targets for chromosome: " + chrom + " moving on.");
					this.daemon.decrementFinishedJobs();
					continue;
				}
				
				this.writeTargetRegions(intervalBed, chromIntervals);
				
				//Create chrom specific lists
				ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
				ArrayList<String> bamList = new ArrayList<String>(); //bam names 
				ArrayList<File> deleteEarly = new ArrayList<File>();
				
				for (String sampleName: samples.keySet()) {
					TFSampleInfo si = samples.get(sampleName).get(0);
					
					//Create file for original bams/split bams/output files
					File bam = si.getFileObject(TFConstants.FILE_REDUCE_BAM).createDestForFileObject(mainJobDir);
					//File bam = new File(mainJobDir,si.getSampleName() + ".reduce.bam");
					File chromBam = new File(chromRunDir,si.getSampleName() + "." + chrom + ".bam");
					File chromBai = new File(chromRunDir,si.getSampleName() + "." + chrom + ".bam.bai");
				
					//Run split
					this.logFile.writeInfoMessage("[TFExomeVariant] Splitting bam file (" + si.getSampleName() + ") by chromosome (" + chrom + ")");
					this.splitByChrom(chrom, bam, chromBam);
					
					//Add bam/bai to preserve
					protectList.add(chromBam);
					protectList.add(chromBai);
					deleteEarly.add(chromBam);
					deleteEarly.add(chromBai);
					
					//Add bam name to list
					bamList.add(si.getSampleName() + "." + chrom + ".bam");
				}
				
				//Create file for original bams/split bams/output files
				File fullVcf = new File(chromRunDir,this.study + ".raw.vcf.gz");
				File fullVcfIdx = new File(chromRunDir,this.study + ".raw.vcf.gz.tbi");
		
				//Add files to list
				deleteList.add(fullVcf);
				deleteList.add(fullVcfIdx);
				fullVcfList.add(fullVcf);
				
				//Create bam String
				String bamString = "";
				for (String b: bamList) {
					bamString += " -I " + b;
				}
				bamString = bamString.substring(1);
				
				//Create replacement tokens
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("STUDY", this.study);
				replacements.put("BAM_LIST", bamString);
				if (targetFile == null) {
					replacements.put("TARGETS","");
				} else {
					//File localTarget = new File(chromRunDir,targetFile.getName());
					replacements.put("TARGETS","-L " + chromIntervals.getName());
					//this.cpFile(targetFile, localTarget);
					deleteList.add(chromIntervals);
					protectList.add(chromIntervals);
				}
				replacements.putAll(this.properties);
				
				//Create cmd.txt file
				File cmdFile = new File(chromRunDir,"cmd.txt");
				protectList.add(cmdFile);
				this.createCmd(replacements, cmdFile,0);
				
				TFThread thread = new TFThread(chromRunDir,this.failmax, counter, this.heartbeat, protectList, deleteEarly, this.logFile);
				this.daemon.addJob(thread);
				counter++;
				
			}
			
			//Wait for command to finish
			try {
				this.daemon.join();
				if (this.daemon.getFailed()) {
					System.exit(1);
				}
				Thread.sleep(5000);
			} catch (InterruptedException ie) {
				logFile.writeErrorMessage("[TFExomeVariant] Daemon interrupted",true);
				System.exit(1);
			}
		}
	
		//Process results
		runPostProcess(sampleList,deleteList,fullVcfList,runDirectoryList);

	}
	
	private void runFull(ArrayList<TFSampleInfo> sampleList, File mainJobDir) {
		File variantFile = new File(mainJobDir,this.study + ".raw.vcf.gz");
		HashSet<File> deleteList = new HashSet<File>();
		HashSet<File> runDirectoryList = new HashSet<File>();
		ArrayList<File> rawVcfList = new ArrayList<File>();
		
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		
		if (!variantFile.exists()) {
			this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(0).getName(), 1, this.jobs);
			this.daemon.start();
			
			runDirectoryList.add(mainJobDir);
			
		    int counter  = 1;
			ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
			ArrayList<String> bamList = new ArrayList<String>(); //bam names 
			
			for (String sampleName: samples.keySet()) {
				TFSampleInfo tfsi = samples.get(sampleName).get(0);
				
				TFFileObject tfoReduceBam = tfsi.getFileObject(TFConstants.FILE_REDUCE_BAM);
				TFFileObject tfoReduceBai = tfsi.getFileObject(TFConstants.FILE_REDUCE_BAI);
				
				File fileReduceBam = tfoReduceBam.createDestForFileObject(mainJobDir);
				File fileReduceBai = tfoReduceBai.createDestForFileObject(mainJobDir);
				
				//Add bam/bai to preserve
				protectList.add(fileReduceBam);
				protectList.add(fileReduceBai);
				deleteList.add(fileReduceBam);
				deleteList.add(fileReduceBai);
				
				//Add bam name to list
				bamList.add(tfoReduceBam.getFileName());
			}
			
			//Create bam String
			String bamString = "";
			for (String b: bamList) {
				bamString += " -I " + b;
			}
			
			//Create replacement tokens
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("STUDY", this.study);
			replacements.put("BAM_LIST", bamString);
			
			if (targetFile == null) {
				replacements.put("TARGETS","");
			} else {
				File localTarget = new File(mainJobDir,targetFile.getName());
				replacements.put("TARGETS","-L " + targetFile.getName());
				this.cpFile(targetFile, localTarget);
				protectList.add(localTarget);
			}
			replacements.putAll(this.properties);
			
			//Create cmd.txt file
			File cmdFile = new File(mainJobDir,"cmd.txt");
			protectList.add(cmdFile);
			this.createCmd(replacements, cmdFile, 0);
			
			TFThread thread = new TFThread(mainJobDir,this.failmax, counter, this.heartbeat, protectList, this.logFile);
			this.daemon.addJob(thread);
			
			//Wait for command to finish
			try {
				this.daemon.join();
				if (this.daemon.getFailed()) {
					System.exit(1);
				}
				Thread.sleep(5000);
			} catch (InterruptedException ie) {
				logFile.writeErrorMessage("[TFExomeVariant] Daemon interrupted",true);
				System.exit(1);
			}
			
			File fullVcf = new File(mainJobDir,this.study + ".raw.vcf.gz");
			File fullVcfIdx = new File(mainJobDir,this.study + ".raw.vcf.gz.tbi");
			
			rawVcfList.add(fullVcf);
			deleteList.add(fullVcf);
			deleteList.add(fullVcfIdx);
		} else {
			rawVcfList.add(variantFile);
		}
		
		runPostProcess(sampleList,deleteList,rawVcfList,runDirectoryList);
			
	}

	
	private void runByChunk(ArrayList<TFSampleInfo> sampleList, File mainJobDir) {
		//Initialize run specific lists
		HashSet<File> deleteList = new HashSet<File>();
		ArrayList<File>	fullVcfList = new ArrayList<File>();
		HashSet<File> runDirectoryList = new HashSet<File>();
		
		HashMap<String,ArrayList<TFSampleInfo>> samples = new HashMap<String,ArrayList<TFSampleInfo>>();
		
		for (TFSampleInfo si: sampleList) {
			if (!samples.containsKey(si.getSampleName())) {
				samples.put(si.getSampleName(), new ArrayList<TFSampleInfo>());
			}
			samples.get(si.getSampleName()).add(si);
		}
		
		
	
		//Count capture space
		double targetSpace = countTargetSpace();
		logFile.writeInfoMessage(String.format("[TFExomeVariant] There are %,.0f bases in the target region",targetSpace));
		
		//Count targetReads
		double readsPerBaseAll = 0;
		
		logFile.writeInfoMessage("[TFExomeVariant] Counting reads in targets");
		for (int i=0; i<sampleList.size();i++) {
			
			double readsPerSample = this.countReadsInTargets(sampleList.get(i));
			int readLength = this.getReadLength(sampleList.get(i).getFileObject(TFConstants.FILE_REDUCE_BAM).getFinalPath());
			
			double readsPerBase = readsPerSample * readLength / targetSpace;
			readsPerBaseAll += readsPerBase;
		}
		
		//Calculate spilt
		int basesPerChunk = Math.round((int)(coveragePerChunk * 10 / readsPerBaseAll));
		int chunks = 0;
		if (targetSpace % basesPerChunk == 0) {
			chunks = (int)(targetSpace / basesPerChunk);
		} else {
			chunks = (int)(targetSpace / basesPerChunk + 1);
		}
		
		logFile.writeInfoMessage(String.format("[TFExomeVariant] Combined reads per base is %,.2f. Splitting data into %,d chunks, with %,d bases per chunk. (%,.0f reads per chunk max).",
				readsPerBaseAll,chunks,basesPerChunk,this.coveragePerChunk));
		
		//Create target regions
		ArrayList<BedEntry> bedFile = this.getTargetRegions();
		int bedFileIndex = 0;
		
		ArrayList<Integer> chunksToRun = new ArrayList<Integer>();
		for (int i=0;i<chunks;i++) {
			File chunkRunDir = new File(mainJobDir,"JOB_" + String.valueOf(i) + "_variant");
			File chunkVcf = new File(chunkRunDir,this.study + ".raw.vcf.gz");
			
			if (!chunkVcf.exists()) {
				chunksToRun.add(i);
			} else {
				fullVcfList.add(chunkVcf);
			}
		}
		
		if (chunksToRun.size() > 0) {
			//Initialize daemon
			this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(0).getName(), chunksToRun.size(), this.jobs);
			this.daemon.start();
			
			int counter = 1;
			
			for (Integer chunk: chunksToRun) {
				//Check if there is space on the cluster for another job
				boolean firstPass = true;
				while(true) {
					if (firstPass && this.daemon.isMaxed()) {
						logFile.writeInfoMessage("[TFExomeVariant] Pausing splitting, waiting for servers to catch up " + String.valueOf(this.daemon.getActive()));
						firstPass = false;
					} else if (!this.daemon.isMaxed()) {
						logFile.writeInfoMessage("[TFExomeVariant] Resuming splitting " + String.valueOf(this.daemon.getActive()));
						break;
					}
					
					try {
						Thread.sleep(5000);
					} catch (InterruptedException e) {
						logFile.writeErrorMessage("[TFExomeVariant] Interupted while waiting to resume bam splitting, exiting",true);
						System.exit(1);
					}
				}
				
				//Create job directory.
				File chunkRunDir = new File(mainJobDir,"JOB_" + String.valueOf(chunk) + "_variant");
				if (chunkRunDir.exists()) {
					this.deleteFolder(chunkRunDir);
				}
				chunkRunDir.mkdir();
				runDirectoryList.add(chunkRunDir);
				
				File chunkIntervals = new File(chunkRunDir,"interval.bed");
				ArrayList<BedEntry> intervalBed = new ArrayList<BedEntry>();
				int basesInChunk = 0;
				
				//Create chunk specific lists
				ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
				ArrayList<File> deleteEarly = new ArrayList<File>();
				ArrayList<String> bamList = new ArrayList<String>(); //bam names 
				
				
				//Create interval based on chunk size
				ArrayList<String> samtoolsRegion = null;
				while (true) {
					BedEntry workingEntry = bedFile.get(bedFileIndex);
					int entrySize = workingEntry.getEnd() - workingEntry.getStart();
					
					if (entrySize + basesInChunk < basesPerChunk) {
						intervalBed.add(workingEntry);
						bedFileIndex+=1;
						basesInChunk += entrySize;
						if (bedFileIndex == (bedFile.size())) {
							samtoolsRegion = this.writeTargetRegions(intervalBed,chunkIntervals);
							break;
						}
					} else if (entrySize + basesInChunk == basesPerChunk) {
						intervalBed.add(workingEntry);
						bedFileIndex += 1;
						basesInChunk += entrySize;
						samtoolsRegion = this.writeTargetRegions(intervalBed, chunkIntervals);
						break;
					} else {
						int overhang = basesPerChunk - basesInChunk;
						int newBoundary = workingEntry.getStart() + overhang;
						BedEntry newbe = new BedEntry(workingEntry.getChrom(),workingEntry.getStart(),newBoundary);
						intervalBed.add(newbe);
						bedFile.get(bedFileIndex).setStart(newBoundary);
						samtoolsRegion = this.writeTargetRegions(intervalBed,chunkIntervals);
						break;
					}
				}
				
				//Create chunked bam files
				for (String sampleName: samples.keySet()) {
					TFSampleInfo si = samples.get(sampleName).get(0);
				
					//Create file for original bams/split bams/output files
					File bam = si.getFileObject(TFConstants.FILE_REDUCE_BAM).createDestForFileObject(mainJobDir);
					File chunkBam = new File(chunkRunDir,si.getSampleName() + "." + String.valueOf(chunk) + ".bam");
					File chunkBai = new File(chunkRunDir,si.getSampleName() + "." + String.valueOf(chunk) + ".bam.bai");
				
					//Run split
					this.logFile.writeInfoMessage("[TFExomeVariant] Splitting bam file (" + si.getSampleName() + ") by chunk (" + String.valueOf(chunk) + ")");
					this.splitByChunk(samtoolsRegion, bam, chunkBam);
					
					//Add files to the delete list or cleanup protection list
					protectList.add(chunkBam);
					protectList.add(chunkBai);
					
					deleteEarly.add(chunkBam);
					deleteEarly.add(chunkBai);
					
					
					//Add bam name to list
					bamList.add(si.getSampleName() + "." + String.valueOf(chunk)  + ".bam");

				}
				
				//Create files for chunked vcf files
				File fullVcf = new File(chunkRunDir,this.study + ".raw.vcf.gz");
				File fullVcfIdx = new File(chunkRunDir,this.study + ".raw.vcf.gz.tbi");

				//Add files to cleanup list
				deleteList.add(fullVcf);
				deleteList.add(fullVcfIdx);
				
				//Add vcf files to merge lists
				fullVcfList.add(fullVcf);
				
				//Create bam String
				String bamString = "";
				for (String b: bamList) {
					bamString += " -I " + b;
				}
				bamString = bamString.substring(1);
				
				//Create replacement tokens
				HashMap<String,String> replacements = new HashMap<String,String>();
				replacements.put("STUDY", this.study);
				replacements.put("BAM_LIST", bamString);
				if (targetFile == null) {
					replacements.put("TARGETS","");
				} else {
					//File localTarget = new File(chunkRunDir,targetFile.getName());
					replacements.put("TARGETS","-L " + chunkIntervals.getName());
					protectList.add(chunkIntervals);
				}
				replacements.putAll(this.properties);
				
				//Create cmd.txt file
				File cmdFile = new File(chunkRunDir,"cmd.txt");
				protectList.add(cmdFile);
				this.createCmd(replacements, cmdFile, 0);
				
				//Run
				TFThread thread = new TFThread(chunkRunDir,this.failmax, counter, this.heartbeat, protectList, deleteEarly, this.logFile);
				this.daemon.addJob(thread);
				counter++;
				
			}
			
			//Wait for command to finish
			try {
				this.daemon.join();
				if (this.daemon.getFailed()) {
					System.exit(1);
				}
				Thread.sleep(5000);
			} catch (InterruptedException ie) {
				logFile.writeErrorMessage("[TFExomeVariant] Daemon interrupted",true);
				System.exit(1);
			}
		}
		
		
		
		runPostProcess(sampleList,deleteList,fullVcfList,runDirectoryList);
	}
	
	
	private void runPostProcess(ArrayList<TFSampleInfo> sampleList, HashSet<File> deleteList, ArrayList<File> fullVcfList, HashSet<File> runDirectoryList) {
		//Clean/Create merge directory
		TFSampleInfo repSample = sampleList.get(0);
		
		File mergeDirectory = repSample.getFileObject(TFConstants.FILE_VCF_RAW).getWorkingDirectory();
		if (mergeDirectory.exists()) {
			this.deleteFolder(mergeDirectory);
		}
		mergeDirectory.mkdir();
		
		//Clean/Create final directory
		File finalDir = repSample.getFileObject(TFConstants.FILE_VCF_RAW).getFinalDirectory();
		if (finalDir.exists()) {
			this.deleteFolder(finalDir);
		}
		finalDir.mkdir();
		
		//Merge output files
		File rawVcfSource = repSample.getFileObject(TFConstants.FILE_VCF_RAW).getWorkingPath();
		File rawVcfSourceIdx = repSample.getFileObject(TFConstants.FILE_VCF_RAW_IDX).getWorkingPath();
		logFile.writeInfoMessage("[TFExomeVariant] Merging raw vcf files");
		if (fullVcfList.size() == 1) {
			this.cpFile(fullVcfList.get(0), rawVcfSource);
		} else {
			File tempFile = new File(mergeDirectory,"temp.vcf");
			File tempFileGz = new File(mergeDirectory,"temp.vcf.gz");
			this.mergeVcf(fullVcfList,tempFile);
			this.bgzip(tempFile);
			this.moveFile(tempFileGz, rawVcfSource);
			deleteList.addAll(fullVcfList);
		}
		this.tabix(rawVcfSource);
		
		this.finalDirectory.mkdirs();
		this.jobDirectory.mkdir();
		
		if (variantBest) {
			logFile.writeInfoMessage("[TFExomeVariant] Moving vcf files");
			
			TFFileObject tfoRawVcf = repSample.getFileObject(TFConstants.FILE_VCF_RAW);
			TFFileObject tfoRawVcfIdx = repSample.getFileObject(TFConstants.FILE_VCF_RAW_IDX);
			
			this.moveFile(tfoRawVcf.getWorkingPath(),tfoRawVcf.getFinalPath());
			this.moveFile(tfoRawVcfIdx.getWorkingPath(),tfoRawVcfIdx.getFinalPath());
		} else {
			int counter  = 1;
			ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
			protectList.add(rawVcfSource);
			protectList.add(rawVcfSourceIdx);
			
			//Add merged to runDirectoryList
			runDirectoryList.add(mergeDirectory);
			
			//Initialize daemon
			this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(1).getName(), 1, this.jobs);
			this.daemon.start();
			
			//Create replacement tokens
			HashMap<String,String> replacements = new HashMap<String,String>();
			replacements.put("STUDY", this.study);
			replacements.putAll(this.properties);
			
			//Create cmd.txt file
			File cmdFile = new File(mergeDirectory,"cmd.txt");
			protectList.add(cmdFile);
			this.createCmd(replacements, cmdFile, 1);
			
			TFThread thread = new TFThread(mergeDirectory,this.failmax, counter, this.heartbeat, protectList, this.logFile);
			this.daemon.addJob(thread);
			
			//Wait for command to finish
			try {
				this.daemon.join();
				if (this.daemon.getFailed()) {
					System.exit(1);
				}
				Thread.sleep(5000);
			} catch (InterruptedException ie) {
				logFile.writeErrorMessage("[TFExomeVariant] Daemon interrupted",true);
				System.exit(1);
			}
			
			//Initialze TFO objects
			TFFileObject tfoRawVcf = repSample.getFileObject(TFConstants.FILE_VCF_RAW);
			TFFileObject tfoRawVcfIdx = repSample.getFileObject(TFConstants.FILE_VCF_RAW_IDX);
			TFFileObject tfoFilterVcf = repSample.getFileObject(TFConstants.FILE_VCF_FILTER);
			TFFileObject tfoFilterVcfIdx = repSample.getFileObject(TFConstants.FILE_VCF_FILTER_IDX);
			TFFileObject tfoPassingVcf = repSample.getFileObject(TFConstants.FILE_VCF_PASSING);
			TFFileObject tfoPassingVcfIdx = repSample.getFileObject(TFConstants.FILE_VCF_PASSING_IDX);
			
			//Move output files
			logFile.writeInfoMessage("[TFExomeVariant] Moving vcf files");
			this.moveFile(tfoRawVcf.getWorkingPath(), tfoRawVcf.getFinalPath());
			this.moveFile(tfoRawVcfIdx.getWorkingPath(), tfoRawVcfIdx.getFinalPath());
			this.moveFile(tfoFilterVcf.getWorkingPath(), tfoFilterVcf.getFinalPath());
			this.moveFile(tfoFilterVcfIdx.getWorkingPath(), tfoFilterVcfIdx.getFinalPath());
			this.moveFile(tfoPassingVcf.getWorkingPath(), tfoPassingVcf.getFinalPath());
			this.moveFile(tfoPassingVcfIdx.getWorkingPath(), tfoPassingVcfIdx.getFinalPath());
		}
	
		this.cleanup(runDirectoryList, deleteList);
	}
	
	
	
	
	private ArrayList<String> writeTargetRegions(ArrayList<BedEntry> bedFile,File outFile) {
		ArrayList<String> regions = new ArrayList<String>();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			
			ArrayList<String> order = new ArrayList<String>();
			HashMap<String,Integer> startPos = new HashMap<String,Integer>();
			HashMap<String,Integer> endPos = new HashMap<String,Integer>();
			
			for (BedEntry be: bedFile) {
				String chrom = be.getChrom();
				int start = be.getStart();
				int end = be.getEnd();
				if (!order.contains(chrom)) {
					order.add(chrom);
					startPos.put(chrom, Integer.MAX_VALUE);
					endPos.put(chrom, Integer.MIN_VALUE);
				}
				if (start < startPos.get(chrom)) {
					startPos.put(chrom, start);
				}
				if (end > endPos.get(chrom)) {
					endPos.put(chrom, end);
				}
				bw.write(be.writeEntry());
			}
			
			bw.close();
			
			for (String c: order) {
				regions.add(" " + c + ":" + Integer.valueOf(startPos.get(c)) + "-" + Integer.valueOf(endPos.get(c)));
			}
			
			
		} catch (FileNotFoundException fnfe) {
			this.logFile.writeErrorMessage("[writeTargetRegions] Count not find target file", true);
			System.exit(1);
		} catch (IOException ioex) {
			this.logFile.writeErrorMessage("[writeTargetRegions] Error reading target file",true);
			System.exit(1);
		}
		return regions;
	}
	
	private ArrayList<BedEntry> getTargetRegions() {
		ArrayList<BedEntry> bedFile = new ArrayList<BedEntry>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.targetFile));
		
			String line = null;
			
			while((line = br.readLine()) != null) {
				String[] items = line.split("\t");
				BedEntry be = new BedEntry(items[0],Integer.parseInt(items[1]),Integer.parseInt(items[2]));
				bedFile.add(be);
			}
			
			br.close();
			
			
		}catch (FileNotFoundException fnfe) {
			this.logFile.writeErrorMessage("[getTargetRegion] Count not find target file", true);
		} catch (IOException ioex) {
			this.logFile.writeErrorMessage("[getTargetRegion] Error reading target file",true);
		}
		
		Collections.sort(bedFile,new BedEntryComparator());
		return bedFile;
	}
	
	
	private double countTargetSpace() {
		double count = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.targetFile));
		
			String line = null;
			
			while((line = br.readLine()) != null) {
				String[] items = line.split("\t");
				double chunk = Double.parseDouble(items[2]) - Double.parseDouble(items[1]);
				count += chunk;
			}
			
			br.close();
		} catch (FileNotFoundException fnfe) {
			this.logFile.writeErrorMessage("[countTargetSpace] Count not find target file", true);
		} catch (IOException ioex) {
			this.logFile.writeErrorMessage("[countTargetSpace] Error reading target file",true);
		}
		
		return count;
	}
	
	private double countReadsInTargets(TFSampleInfo si) {
//		double n = -1;
//		String fileName = si.getFile(TFConstants.FILE_REALIGN_SAMPLE_BAM).getAbsolutePath();
//		
//		try {
//			ProcessBuilder pb = new ProcessBuilder("/home/u0855942/applications/samtools/samtools","view","-c","@","10",fileName);
//			
//			Process p = pb.start();
//			
//			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
//			String line = null;
//			String result = "";
//			while((line = br.readLine()) != null) {
//				result += line;
//				
//			}
//			
//			int val = p.waitFor(); 
//			
//			n = Double.parseDouble(result);
//					
//			if (val != 0) {
//				logFile.writeErrorMessage("[countReadsInTargets] Error when counting reads in bam file\n",true);
//				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//				String line2 = null;
//				while((line2 = br2.readLine()) != null) {
//					System.out.println(line2);
//				}
//				System.exit(1);
//			}
//			
//		} catch (IOException ioex) {
//			logFile.writeErrorMessage("[countReadsInTargers] IO Exception while trying to count reads in bam: " + fileName,true);
//			System.exit(1);
//		} catch (InterruptedException ieex) {
//			logFile.writeErrorMessage("[splitByChrom] Process was interrupted while trying count reads in bam: " + fileName,true);
//			System.exit(1);
//		}
		
		return 50000000;
	}
	
	private void splitByChrom(String chrom, File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[splitByChrom] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			ProcessBuilder pb = new ProcessBuilder(properties.get("SAM_PATH_LOCAL") + "/samtools","view","-b","-h","-@","10",source.getAbsolutePath(),chrom);
			Process p = pb.start();
			
			
			
			BufferedInputStream bis = new BufferedInputStream(p.getInputStream());
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(dest));
			
			
			byte[] buffer = new byte[1024*1024*10];
			int n = -1;
			
			while((n = bis.read(buffer))!=-1) {
			  bos.write(buffer,0,n);
			}
		

			int val = p.waitFor();
			bos.close();
			bis.close();
			
			if (val != 0) {
				logFile.writeErrorMessage("[splitByChrom] Error while splitting the bam file: " + chrom + " "+ source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
			
			
			pb = new ProcessBuilder("/tomato/app/samtools/samtools","index",dest.getAbsolutePath());
			p = pb.start();
			
			val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[splitByChrom] Error while indexing the bam file: " + chrom + " "+ source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[splitByChrom] IO Exception while trying to split the bam file: " + chrom + " " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[splitByChrom] Process was interrupted while trying to split bam files: " + chrom + " " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
	
	private void splitByChunk(ArrayList<String> samRegion, File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[splitByChunk] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			String[] commandArray = {properties.get("SAM_PATH_LOCAL") + "/samtools","view","-b","-h","-@","10",source.getAbsolutePath()};
			ArrayList<String> commandList = new ArrayList<String>(Arrays.asList(commandArray));
			commandList.addAll(samRegion);
			
			
			ProcessBuilder pb = new ProcessBuilder(commandList);
			Process p = pb.start();
			
			BufferedInputStream bis = new BufferedInputStream(p.getInputStream());
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(dest));
			
			
			byte[] buffer = new byte[1024*1024*10];
			int n = -1;
			
			while((n = bis.read(buffer))!=-1) {
			  bos.write(buffer,0,n);
			}
		

			int val = p.waitFor();
			bos.close();
			bis.close();
			
			if (val != 0) {
				logFile.writeErrorMessage("[splitByChunk] Error while splitting the bam file: " + samRegion + " "+ source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
			pb = new ProcessBuilder(properties.get("SAM_PATH_LOCAL") + "/samtools","index",dest.getAbsolutePath());
			p = pb.start();
			
			val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[splitByChunk] Error while indexing the bam file: " + samRegion + " "+ source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[splitByChunk] IO Exception while trying to split the bam file: " + samRegion + " " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[splitByChunk] Process was interrupted while trying to split bam files: " + samRegion + " " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
	
	
	public ArrayList<TFSampleInfo> get1KGenomes() {
		File oneKDir = new File(this.properties.get("BACKGROUND_PATH"));
		if (oneKDir.canRead()) {
			File[] bams = oneKDir.listFiles();
			ArrayList<TFSampleInfo> oneKTFList = new ArrayList<TFSampleInfo>();
			Pattern bamPattern = Pattern.compile("(.+?)\\.reduced\\.bam$");
			for (File bam: bams) {
				String name = bam.getName();
				Matcher m = bamPattern.matcher(name);
				if (m.matches()) {
					String sampleName = m.group(1);
					TFSampleInfo nsi = new TFSampleInfo(sampleName,sampleName,"",this.logFile);
					TFFileObject tfoReduceBam = new TFFileObject(sampleName + ".reduced.bam",oneKDir,oneKDir);
					TFFileObject tfoReduceBai = new TFFileObject(sampleName + ".reduced.bai",oneKDir,oneKDir);
					
					nsi.setFileObject(TFConstants.FILE_REDUCE_BAM, tfoReduceBam);
					nsi.setFileObject(TFConstants.FILE_REDUCE_BAI, tfoReduceBai);
		
					oneKTFList.add(nsi);
				}
				
			}
			return oneKTFList;
		} else {
			return null;
		}
		
	}
	
	private int getReadLength(File bamFile){
//		SAMFileReader samReader = null;
//		
//		int readLength = 0;
//		int counter = 0;
//			
//		System.out.println(bamFile.getAbsolutePath());
//		samReader = new SAMFileReader(bamFile);
//		samReader.setValidationStringency(ValidationStringency.SILENT);
//		SAMRecordIterator it = samReader.iterator();
//
//		while (it.hasNext()) {
//			
//			try {
//				SAMRecord sam = it.next();
//				
//				int len = sam.getReadString().length();
//				if (len > readLength) {
//					readLength = len;
//				}
//				
//				counter++;
//				if (counter > 10000) {
//					break;
//				}
//			} catch (SAMFormatException e) {
//				System.out.println("Error reading sam");
//			}
//		}
//			
//		samReader.close();
	
//		return readLength;
		return 101;
	}
			
	
	
	class BedEntry {
		private String chrom;
		private int start;
		private int end;
		
		public BedEntry(String chrom, int start, int end) {
			this.chrom = chrom;
			this.start = start;
			this.end = end;
		}
		
		public String getChrom() {
			return this.chrom;
		}
		
		public int getStart() {
			return this.start;
		}
		
		public int getEnd() {
			return this.end;
		}
		
		public void setStart(int start) {
			this.start = start;
		}
		
		public String writeEntry() {
			String outline = String.format("%s\t%d\t%d\n",this.chrom,this.start,this.end);
			return outline;
		}

	}
	
	class BedEntryComparator implements Comparator<BedEntry> {
		
		private String[] chroms = {"1","2","3","4","5","6","7","8",
                "9","10","11","12","13","14","15",
                "16","17","18","19","20","21","22",
                "X","Y","MT","chr1","chr2","chr3","chr4","chr5",
                "chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX","chrY"};

		@Override
		public int compare(BedEntry be1, BedEntry be2) {

			
			int rank1 = Arrays.asList(chroms).indexOf(be1.getChrom());
			int rank2 = Arrays.asList(chroms).indexOf(be2.getChrom());
			
			if (rank1 != rank2) {
				return rank1 - rank2;
			} else {
				return be1.getStart() - be2.getStart();
			}
		}
		
	}

}
