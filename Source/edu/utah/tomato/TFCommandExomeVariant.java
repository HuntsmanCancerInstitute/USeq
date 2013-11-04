package edu.utah.tomato;

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
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;



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
	}

	
	
	
	private void runByChrom(ArrayList<TFSampleInfo> sampleList) {
		this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(0).getName(), hg19Chrom.length, this.jobs);
		this.daemon.start();
		
		int counter = 1;
		
		File runDirectory = new File(this.rootDirectory,"JOB_" + this.study + "_variant");
		runDirectory.mkdir();
		
		//Create target regions
		ArrayList<BedEntry> bedFile = this.getTargetRegions();
		
		
		//Initialize run specific lists
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<File>	fullVcfList = new ArrayList<File>();
		
		for (TFSampleInfo si: sampleList) {
			//Create output files
			File bam = new File(runDirectory,si.getSampleID() + ".bam");
			File bai = new File(runDirectory,si.getSampleID() + ".bai");
			
			deleteList.add(bam);
			deleteList.add(bai);
			
			//Create links
			this.createLink(si.getFile(TFConstants.FILE_REDUCE_BAM), bam);
			this.createLink(si.getFile(TFConstants.FILE_REDUCE_BAI), bai);
		}
		
		for (String chrom: this.hg19Chrom) {
			boolean firstPass = true;
			
			while(true) {
				if (firstPass && this.daemon.isMaxed()) {
					logFile.writeInfoMessage("Pausing splitting, waiting for servers to catch up");
					firstPass = false;
				} else if (!this.daemon.isMaxed()) {
					logFile.writeInfoMessage("Resuming splitting");
					break;
				}
				
				try {
					Thread.sleep(5000);
				} catch (InterruptedException e) {
					logFile.writeErrorMessage("Interupted while waiting to resume bam splitting, exiting",true);
					System.exit(1);
				}
			}
			
			//Create chromosome specific directory
			File chromRunDir = new File(runDirectory,chrom);
			chromRunDir.mkdir();
			
			//Create interval file
			File chromIntervals = new File(chromRunDir,"interval.bed");
			ArrayList<BedEntry> intervalBed = new ArrayList<BedEntry>();
			for (BedEntry interval: bedFile) {
				if (interval.getChrom().equals(chrom)) {
					intervalBed.add(interval);
				}
			}
			
			if (intervalBed.size() == 0) {
				logFile.writeInfoMessage("There are no targets for chromosome: " + chrom + " moving on.");
				this.daemon.decrementFinishedJobs();
				continue;
			}
			
			this.writeTargetRegions(intervalBed, chromIntervals);
			
			
			
			
			//Create chrom specific lists
			ArrayList<File> keepers = new ArrayList<File>(); //files to preserve on cleanup
			ArrayList<String> bamList = new ArrayList<String>(); //bam names 
			ArrayList<File> deleteEarly = new ArrayList<File>();
			
			for (TFSampleInfo si: sampleList) {
				//Create file for original bams/split bams/output files
				File bam = new File(runDirectory,si.getSampleID() + ".bam");
				File chromBam = new File(chromRunDir,si.getSampleID() + "." + chrom + ".bam");
				File chromBai = new File(chromRunDir,si.getSampleID() + "." + chrom + ".bam.bai");
			
				//Run split
				this.logFile.writeInfoMessage("Splitting bam file (" + si.getSampleID() + ") by chromosome (" + chrom + ")");
				this.splitByChrom(chrom, bam, chromBam);
				
				//Add bam/bai to preserve
				keepers.add(chromBam);
				keepers.add(chromBai);
				
				deleteEarly.add(chromBam);
				deleteEarly.add(chromBai);
				//deleteList.add(chromBam);
				//deleteList.add(chromBai);
				
				//Add bam name to list
				bamList.add(si.getSampleID() + "." + chrom + ".bam");

			}
			
			
			
			
			
			
			//Create file for original bams/split bams/output files
			File fullVcf = new File(chromRunDir,this.study + ".vcf");
	
			//Add files to cleanup list
			deleteList.add(fullVcf);
			
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
				//File localTarget = new File(chromRunDir,targetFile.getName());
				replacements.put("TARGETS","-L " + chromIntervals.getName());
				//this.cpFile(targetFile, localTarget);
				deleteList.add(chromIntervals);
				keepers.add(chromIntervals);
			}
			replacements.putAll(this.properties);
			
			//Create cmd.txt file
			File cmdFile = new File(chromRunDir,"cmd.txt");
			keepers.add(cmdFile);
			
			//Create cmd.txt file
			this.createCmd(replacements, cmdFile,0);
			
			//Run this shit!
			TFThread thread = new TFThread(chromRunDir,this.failmax, counter, this.heartbeat, keepers, deleteEarly, this.logFile);
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
			logFile.writeErrorMessage("Daemon interrupted",true);
			System.exit(1);
		}
		
		runVqsr(runDirectory,deleteList,fullVcfList);

	}
	
	private void runFull(ArrayList<TFSampleInfo> sampleList) {
		this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(0).getName(), 1, this.jobs);
		this.daemon.start();
		
		File runDirectory = new File(this.rootDirectory,"JOB_" + this.study + "_variant");
		runDirectory.mkdir();
		
	    int counter  = 1;
		ArrayList<File> keepers = new ArrayList<File>(); //files to preserve on cleanup
		ArrayList<String> bamList = new ArrayList<String>(); //bam names 
		ArrayList<File> deleteList = new ArrayList<File>(); //files to delete on completeion
		
		
		for (TFSampleInfo si: sampleList) {
			
			//Create file for original bams/split bams/output files
			File bam = new File(runDirectory,si.getSampleID() + ".bam");
			File bai = new File(runDirectory,si.getSampleID() + ".bai");
			
			//Create links 
			this.createLink(si.getFile(TFConstants.FILE_REDUCE_BAM), bam);
			this.createLink(si.getFile(TFConstants.FILE_REDUCE_BAI), bai);
			
		
			//Add bam/bai to preserve
			keepers.add(bam);
			keepers.add(bai);
			deleteList.add(bam);
			deleteList.add(bai);
			
			//Add bam name to list
			bamList.add(si.getSampleID() + ".bam");
			
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
			File localTarget = new File(runDirectory,targetFile.getName());
			replacements.put("TARGETS","-L " + targetFile.getName());
			this.cpFile(targetFile, localTarget);
			keepers.add(localTarget);
		}
		replacements.putAll(this.properties);
		
		//Create cmd.txt file
		File cmdFile = new File(runDirectory,"cmd.txt");
		keepers.add(cmdFile);
		
		//Create cmd.txt file
		this.createCmd(replacements, cmdFile, 0);
		
		//Run this shit!
		TFThread thread = new TFThread(runDirectory,this.failmax, counter, this.heartbeat, keepers, this.logFile);
		this.daemon.addJob(thread);
		
		//Wait for command to finish
		try {
			this.daemon.join();
			Thread.sleep(5000);
		} catch (InterruptedException ie) {
			logFile.writeErrorMessage("Daemon interrupted",true);
			System.exit(1);
		}
		
		
		ArrayList<File> rawVcfList = new ArrayList<File>();
		File fullVcf = new File(runDirectory,this.study + ".vcf");
		rawVcfList.add(fullVcf);
		
		runVqsr(runDirectory,deleteList,rawVcfList);
			
	}

	
	private void runByChunk(ArrayList<TFSampleInfo> sampleList) {
		//Initialize run specific lists
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<File>	fullVcfList = new ArrayList<File>();
		
		//Create outer variant directory
		File runDirectory = new File(this.rootDirectory,"JOB_" + this.study + "_variant");
		runDirectory.mkdir();
		
		
		
		
		if (this.targetFile == null) {
			logFile.writeInfoMessage("Calling variations across all regions");
		} else {
			logFile.writeInfoMessage("Calling variations using: " + targetFile.getAbsolutePath());
		}
			
		
		for (TFSampleInfo si: sampleList) {
			//Create output files
			File bam = new File(runDirectory,si.getSampleID() + ".bam");
			File bai = new File(runDirectory,si.getSampleID() + ".bai");
			
			deleteList.add(bam);
			deleteList.add(bai);
			
			//Create links
			this.createLink(si.getFile(TFConstants.FILE_REDUCE_BAM), bam);
			this.createLink(si.getFile(TFConstants.FILE_REDUCE_BAI), bai);
		}
		
		
		//Count capture space
		double targetSpace = countTargetSpace();
		logFile.writeInfoMessage(String.format("There are %,.0f bases in the target region",targetSpace));
		
		//Count targetReads
		double readsPerBaseAll = 0;
		
		logFile.writeInfoMessage("Counting reads in targets");
		for (int i=0; i<sampleList.size();i++) {
			
			double readsPerSample = this.countReadsInTargets(sampleList.get(i));
			int readLength = this.getReadLength(sampleList.get(i).getFile(TFConstants.FILE_REDUCE_BAM));
			
			double readsPerBase = readsPerSample * readLength / targetSpace;
			readsPerBaseAll += readsPerBase;
//			logFile.writeInfoMessage(String.format("There are %,.0f reads in the target region, with an average of %,.2f reads/base in sample %s. "
//					+ "Read length %d",
//					readsPerSample,readsPerBase,sampleList.get(i).getSampleID(),readLength));
		}
		
		//Calculate spilt
		int basesPerChunk = Math.round((int)(coveragePerChunk * 10 / readsPerBaseAll));
		int chunks = 0;
		if (targetSpace % basesPerChunk == 0) {
			chunks = (int)(targetSpace / basesPerChunk);
		} else {
			chunks = (int)(targetSpace / basesPerChunk + 1);
		}
		
		
		logFile.writeInfoMessage(String.format("Combined reads per base is %,.2f. Splitting data into %,d chunks, with %,d bases per chunk. (%,.0f reads per chunk max).",
				readsPerBaseAll,chunks,basesPerChunk,this.coveragePerChunk));
		
		//Create target regions
		ArrayList<BedEntry> bedFile = this.getTargetRegions();
		int bedFileIndex = 0;
		
		//Initialize daemon
		this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(0).getName(), chunks, this.jobs);
		this.daemon.start();
		
		int counter = 1;
		
		for (int i=0;i<chunks;i++) {
			boolean firstPass = true;
			while(true) {
				if (firstPass && this.daemon.isMaxed()) {
					logFile.writeInfoMessage("Pausing splitting, waiting for servers to catch up " + String.valueOf(this.daemon.getActive()));
					firstPass = false;
				} else if (!this.daemon.isMaxed()) {
					logFile.writeInfoMessage("Resuming splitting " + String.valueOf(this.daemon.getActive()));
					break;
				}
				
				try {
					Thread.sleep(5000);
				} catch (InterruptedException e) {
					logFile.writeErrorMessage("Interupted while waiting to resume bam splitting, exiting",true);
					System.exit(1);
				}
			}
			
			//Chunking variables
			File chunkRunDir = new File(runDirectory,String.valueOf(i));
			chunkRunDir.mkdir();
			
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
				
				if (entrySize + basesInChunk <= basesPerChunk) {
					intervalBed.add(workingEntry);
					bedFileIndex+=1;
					basesInChunk += entrySize;
					if (bedFileIndex == (bedFile.size())) {
						samtoolsRegion = this.writeTargetRegions(intervalBed,chunkIntervals);
						break;
					}
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
			for (TFSampleInfo si: sampleList) {
				//Create file for original bams/split bams/output files
				File bam = new File(runDirectory,si.getSampleID() + ".bam");
				File chunkBam = new File(chunkRunDir,si.getSampleID() + "." + String.valueOf(i) + ".bam");
				File chunkBai = new File(chunkRunDir,si.getSampleID() + "." + String.valueOf(i) + ".bam.bai");
			
				//Run split
				this.logFile.writeInfoMessage("Splitting bam file (" + si.getSampleID() + ") by chunk (" + String.valueOf(i) + ")");
				this.splitByChunk(samtoolsRegion, bam, chunkBam);
				
				//Add files to the delete list or cleanup protection list
				protectList.add(chunkBam);
				protectList.add(chunkBai);
				
				deleteEarly.add(chunkBam);
				deleteEarly.add(chunkBai);
				
				//Add files to cleanup list
				//deleteList.add(chunkBam);
				//deleteList.add(chunkBai);
				
				//Add bam name to list
				bamList.add(si.getSampleID() + "." + String.valueOf(i)  + ".bam");

			}
			
			//Create files for chunked vcf files
			File fullVcf = new File(chunkRunDir,this.study + ".vcf");
			File fullVcfIdx = new File(chunkRunDir,this.study + ".vcf.idx");
			
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
			
			//Create cmd.txt file
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
			logFile.writeErrorMessage("Daemon interrupted",true);
			System.exit(1);
		}
		
		runVqsr(runDirectory,deleteList,fullVcfList);
		
			
		//Fin
	}
	
	
	private void runVqsr(File runDirectory, ArrayList<File> deleteList, ArrayList<File> fullVcfList) {
		/* Clean up the results from the individual SNP call runs
		 */
		
		File mergeDirectory = new File(runDirectory,"merged_calls");
		mergeDirectory.mkdir();
		
		//Clean destination directory
		for(File file: mergeDirectory.listFiles()){
			if (file.exists()) {
				file.delete();
			}
		}
		
		//Make destination files
		File fullVcfSource = new File(mergeDirectory,this.study + ".vcf");
		File filterVcfSource = new File(mergeDirectory,this.study + ".filtered.vcf");
		File passingVcfSource = new File(mergeDirectory,this.study + ".passing.vcf");

		
		//Merge output files
		logFile.writeInfoMessage("Merging raw vcf files");
		if (fullVcfList.size() == 1) {
			this.moveFile(fullVcfList.get(0), fullVcfSource);
		} else {
			this.mergeVcf(fullVcfList,fullVcfSource);
		}
		
		//Delete stragglers
		for (File f: deleteList) {
			if (f.exists()) {
				this.deleteFile(f);
			}
		}
		
		int counter  = 1;
		ArrayList<File> keepers = new ArrayList<File>(); //files to preserve on cleanup
		keepers.add(fullVcfSource);
		
		//Initialize daemon
		this.daemon = new TFThreadDaemon(this.logFile, this.templateFiles.get(1).getName(), 1, this.jobs);
		this.daemon.start();
		
		//Create replacement tokens
		HashMap<String,String> replacements = new HashMap<String,String>();
		replacements.put("STUDY", this.study);
		replacements.putAll(this.properties);
		
		//Create cmd.txt file
		File cmdFile = new File(mergeDirectory,"cmd.txt");
		keepers.add(cmdFile);
		
		//Create cmd.txt file
		this.createCmd(replacements, cmdFile, 1);
		
		//Run this shit!
		TFThread thread = new TFThread(mergeDirectory,this.failmax, counter, this.heartbeat, keepers, this.logFile);
		this.daemon.addJob(thread);
		
		//Wait for command to finish
		try {
			this.daemon.join();
			Thread.sleep(5000);
		} catch (InterruptedException ie) {
			logFile.writeErrorMessage("Daemon interrupted",true);
			System.exit(1);
		}
		
		//Make destination directory
		File varDir = new File(this.rootDirectory,"Variants");
		varDir.mkdir();
		for(File file: varDir.listFiles()) file.delete();
		File jobDir = new File(varDir,"Jobs");
		
		//Make destination files
		File fullVcfDest = new File(varDir,this.study + ".raw.vcf");
		File filterVcfDest = new File(varDir,this.study + ".filterFieldSetAll.vcf");
		File passingVcfDest = new File(varDir,this.study + ".filterFieldSetPassing.vcf");
		
		File fullVcfDestGz = new File(varDir,this.study + ".raw.vcf.gz");
		File filterVcfDestGz = new File(varDir,this.study + ".filterFieldSetAll.vcf.gz");
		File passingVcfDestGz = new File(varDir,this.study + ".filterFieldSetPassing.vcf.gz");
		
		//Merge output files
		logFile.writeInfoMessage("Moving vcf files");
		this.moveFile(fullVcfSource, fullVcfDest);
		this.moveFile(filterVcfSource,filterVcfDest);
		this.moveFile(passingVcfSource, passingVcfDest);
		
		//Compress
		logFile.writeInfoMessage("Compressing vcf files");
		this.bgzip(fullVcfDest);
		this.bgzip(filterVcfDest);
		this.bgzip(passingVcfDest);
		
		//Index
		logFile.writeInfoMessage("Indexing vcf files");
		this.tabix(fullVcfDestGz);
		this.tabix(filterVcfDestGz);
		this.tabix(passingVcfDestGz);
		
		
		File existDir = new File(jobDir,mergeDirectory.getName());
		if (existDir.exists()) {
			deleteFolder(existDir);
		}
		
		this.moveFile(runDirectory, jobDir);
	}
	
	@Override
	public void run(ArrayList<TFSampleInfo> sampleList) {
		//Use the final list object.
		if (this.use1KGenomes) {
			logFile.writeInfoMessage("Attempting to load 1K genomes samples");
			ArrayList<TFSampleInfo> onekList = get1KGenomes();
			if (onekList == null) {
				logFile.writeWarningMessage("You don't have the permissions to access the 1K genomes archive, moving on withouth them");
			} else {
				sampleList.addAll(onekList);
			}
		}
		
		System.out.println(sampleList.size());
		
		if (splitType.equals("chrom")) {
			runByChrom(sampleList);
		} else if (splitType.equals("chunk")) {
			runByChunk(sampleList);
		} else if (splitType.equals("none")) {
			runFull(sampleList);
		} else {
			logFile.writeErrorMessage("Don't recognize run type: " + splitType, true);
			System.exit(1);
		}
		
		//Clean up reduced files if specified
		this.logFile.writeInfoMessage("Deleting raw alignments");
		if (this.deleteReducedBams && this.isFull) {
			for (TFSampleInfo si: sampleList) {
				if (si.getFile(TFConstants.FILE_REDUCE_BAM).exists()) {
					si.getFile(TFConstants.FILE_REDUCE_BAM).delete();
				}
				if (si.getFile(TFConstants.FILE_REDUCE_BAI).exists()) {
					si.getFile(TFConstants.FILE_REDUCE_BAI).delete();
				}
			}
		}
		

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
		
		return 20000000;
	}
	
	private void splitByChrom(String chrom, File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[splitByChrom] Expected file does not exist: " + source.getAbsolutePath(),true);
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
				logFile.writeErrorMessage("[splitByChrom] Expected file does not exist: " + source.getAbsolutePath(),true);
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
				String dir = bam.getParent();
				String name = bam.getName();
				Matcher m = bamPattern.matcher(name);
				if (m.matches()) {
					TFSampleInfo nsi = new TFSampleInfo(m.group(1),this.logFile);
					nsi.setFile(TFConstants.FILE_REDUCE_BAM, new File(dir,m.group(1) + ".reduced.bam"));
					nsi.setFile(TFConstants.FILE_REDUCE_BAI, new File(dir,m.group(1) + ".reduced.bai"));
					oneKTFList.add(nsi);
				}
				
			}
			return oneKTFList;
		} else {
			return null;
		}
		
	}
	
	public int getReadLength(File bamFile){
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
                "X","Y","MT"};

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
