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
import java.util.HashSet;


import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import edu.utah.seq.data.sam.SamAlignment;


public class TFCommandExomeVariantRaw extends TFCommand {
	
	private String splitType = "none";
	private String study = null;
	private String[] hg19Chrom = {"chr1","chr2","chr3","chr4","chr5","chr6",
			"chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
			"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
			"chrX","chrY","chrM"};
	private File targetFile;
	private int coveragePerChunk = 1000000000;
	

	
	
	
	public TFCommandExomeVariantRaw(File templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, String study, String splitType, File targetFile) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress);
		this.splitType = splitType;
		this.study = study;
		this.targetFile = targetFile;
	}

	
	
	
	private void runByChrom(ArrayList<TFSampleInfo> sampleList) {
		this.daemon = new TFThreadDaemon(this.logFile, this.commandString, hg19Chrom.length, this.jobs);
		this.daemon.start();
		
		int counter = 1;
		
		File runDirectory = new File(this.rootDirectory,"JOB_" + this.study + "_variant");
		runDirectory.mkdir();
		
		
		
		//Initialize run specific lists
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<File>	fullVcfList = new ArrayList<File>();
		ArrayList<File>	filteredVcfList = new ArrayList<File>();
		ArrayList<File>	passingVcfList = new ArrayList<File>();
		
		for (TFSampleInfo si: sampleList) {
			//Create output files
			File bam = new File(runDirectory,si.getSampleName() + ".bam");
			File bai = new File(runDirectory,si.getSampleName() + ".bai");
			
			deleteList.add(bam);
			deleteList.add(bai);
			
			//Create links
			this.createLink(si.getFile(TFConstants.FILE_BAM), bam);
			this.createLink(si.getFile(TFConstants.FILE_BAI), bai);
		}
		
		for (String chrom: this.hg19Chrom) {
			//Create chromosome specific directory
			File chromRunDir = new File(runDirectory,chrom);
			chromRunDir.mkdir();
			
			//Create chrom specific lists
			ArrayList<File> keepers = new ArrayList<File>(); //files to preserve on cleanup
			ArrayList<String> bamList = new ArrayList<String>(); //bam names 
			
			
			for (TFSampleInfo si: sampleList) {
				//Create file for original bams/split bams/output files
				File bam = new File(runDirectory,si.getSampleName() + ".bam");
				File chromBam = new File(chromRunDir,si.getSampleName() + "." + chrom + ".bam");
				File chromBai = new File(chromRunDir,si.getSampleName() + "." + chrom + ".bam.bai");
			
				//Run split
				this.logFile.writeInfoMessage("Splitting bam file (" + si.getSampleName() + ") by chromosome (" + chrom + ")");
				this.splitByChrom(chrom, bam, chromBam);
				
				//Add bam/bai to preserve
				keepers.add(chromBam);
				keepers.add(chromBai);
				
				//Add files to cleanup list
				deleteList.add(chromBam);
				deleteList.add(chromBai);
				
				//Add bam name to list
				bamList.add(si.getSampleName() + "." + chrom + ".bam");

			}
			
			//Create file for original bams/split bams/output files
			File fullVcf = new File(chromRunDir,this.study + ".vcf");
			File filterVcf = new File(chromRunDir,this.study + ".filtered.vcf");
			File passVcf = new File(chromRunDir,this.study + ".passing.vcf");
	
			//Add files to cleanup list
			deleteList.add(fullVcf);
			deleteList.add(filterVcf);
			deleteList.add(passVcf);
			
			//Add vcf files to merge lists
			fullVcfList.add(fullVcf);
			filteredVcfList.add(filterVcf);
			passingVcfList.add(passVcf);
			
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
				File localTarget = new File(chromRunDir,targetFile.getName());
				replacements.put("TARGETS","-L " + targetFile.getName());
				this.cpFile(targetFile, localTarget);
				deleteList.add(localTarget);
				keepers.add(localTarget);
			}
			
			//Create cmd.txt file
			File cmdFile = new File(chromRunDir,"cmd.txt");
			keepers.add(cmdFile);
			
			//Create cmd.txt file
			this.createCmd(replacements, cmdFile);
			
			//Run this shit!
			TFThread thread = new TFThread(chromRunDir,this.failmax, counter, this.heartbeat, keepers, this.logFile);
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
		
	    cleanupRun(fullVcfList,filteredVcfList,passingVcfList,deleteList,runDirectory);
		
	}
	
	private void runFull(ArrayList<TFSampleInfo> sampleList) {
		this.daemon = new TFThreadDaemon(this.logFile, this.commandString, 1, this.jobs);
		this.daemon.start();
		
		File runDirectory = new File(this.rootDirectory,"JOB_" + this.study + "_variant");
		runDirectory.mkdir();
		
	    int counter  = 1;
		ArrayList<File> keepers = new ArrayList<File>(); //files to preserve on cleanup
		ArrayList<String> bamList = new ArrayList<String>(); //bam names 
		ArrayList<File> deleteList = new ArrayList<File>(); //files to delete on completeion
		
		
		for (TFSampleInfo si: sampleList) {
			
			//Create file for original bams/split bams/output files
			File bam = new File(runDirectory,si.getSampleName() + ".bam");
			File bai = new File(runDirectory,si.getSampleName() + ".bai");
			
			//Create links 
			this.createLink(si.getFile(TFConstants.FILE_BAM), bam);
			this.createLink(si.getFile(TFConstants.FILE_BAI), bai);
			
		
			//Add bam/bai to preserve
			keepers.add(bam);
			keepers.add(bai);
			deleteList.add(bam);
			deleteList.add(bai);
			
			//Add bam name to list
			bamList.add(si.getSampleName() + ".bam");
			
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
		
		//Create cmd.txt file
		File cmdFile = new File(runDirectory,"cmd.txt");
		keepers.add(cmdFile);
		
		//Create cmd.txt file
		this.createCmd(replacements, cmdFile);
		
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
		
		
		File varDir = new File(this.rootDirectory,"Variants");
		varDir.mkdir();
		for(File file: varDir.listFiles()) file.delete();
		File jobDir = new File(varDir,"Jobs");
		
		
		//Make Source files
		File fullVcf = new File(runDirectory,this.study + ".vcf");
		File filterVcf = new File(runDirectory,this.study + ".filtered.vcf");
		File passVcf = new File(runDirectory,this.study + ".passing.vcf");
	
		
		//Make destination files
		File fullVcfDest = new File(varDir,this.study + ".full.vcf");
		File fullVcfDestGz = new File(varDir,this.study + ".full.vcf.gz");
		File filterVcfDest = new File(varDir,this.study + ".filtered.vcf");
		File filterVcfDestGz = new File(varDir,this.study + ".filtered.vcf.gz");
		File passingVcfDest = new File(varDir,this.study + ".passing.vcf");
		File passingVcfDestGz = new File(varDir,this.study + ".passing.vcf.gz");
		
		
		//move files
		this.moveFile(fullVcf, fullVcfDest);
		this.moveFile(filterVcf, filterVcfDest);
		this.moveFile(passVcf, passingVcfDest);
		
		//Compress
		this.bgzip(fullVcfDest);
		this.bgzip(filterVcfDest);
		this.bgzip(passingVcfDest);
		
		//Index
		this.tabix(fullVcfDestGz);
		this.tabix(filterVcfDestGz);
		this.tabix(passingVcfDestGz);
		

		for(File i : deleteList) {
			i.delete();
		}
		
		this.moveFile(runDirectory, jobDir);
			
			
			
	}

	
	private void runByChunk(ArrayList<TFSampleInfo> sampleList) {
		//Initialize run specific lists
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<File>	fullVcfList = new ArrayList<File>();
		ArrayList<File>	filteredVcfList = new ArrayList<File>();
		ArrayList<File>	passingVcfList = new ArrayList<File>();
		
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
			File bam = new File(runDirectory,si.getSampleName() + ".bam");
			File bai = new File(runDirectory,si.getSampleName() + ".bai");
			
			deleteList.add(bam);
			deleteList.add(bai);
			
			//Create links
			this.createLink(si.getFile(TFConstants.FILE_BAM), bam);
			this.createLink(si.getFile(TFConstants.FILE_BAI), bai);
		}
		
		
		//Count capture space
		int targetSpace = countTargetSpace();
		logFile.writeInfoMessage(String.format("There are %,d bases in the target region",targetSpace));
		
		//Count targetReads
		float readsPerBaseAll = 0;
		
		logFile.writeInfoMessage("Counting reads in targets");
		for (int i=0; i<sampleList.size();i++) {
			
			int readsPerSample = this.countReadsInTargets(sampleList.get(i));
			int readLength = this.getReadLength(sampleList.get(i).getFile(TFConstants.FILE_BAM));
			
			float readsPerBase = (float)readsPerSample * readLength / targetSpace;
			readsPerBaseAll += readsPerBase;
			logFile.writeInfoMessage(String.format("There are %,d reads in the target region, with an average of %,.2f reads/base in sample %s",
					readsPerSample,readsPerBase,sampleList.get(i).getSampleName()));
		}
		
		//Calculate spilt
		int basesPerChunk = Math.round((coveragePerChunk * 10 / readsPerBaseAll));
		int chunks = targetSpace / basesPerChunk + 1;
		
		logFile.writeInfoMessage(String.format("Combined reads per base is %,.2f. Splitting data into %,d chunks, with %,d bases per chunk. (%,d reads per chunk max).",
				readsPerBaseAll,chunks,basesPerChunk,this.coveragePerChunk));
		
		//Create target regions
		ArrayList<BedEntry> bedFile = this.getTargetRegions();
		int bedFileIndex = 0;
		
		//Initialize daemon
		this.daemon = new TFThreadDaemon(this.logFile, this.commandString, chunks, this.jobs);
		this.daemon.start();
		
		int counter = 1;
		
		for (int i=0;i<chunks;i++) {
			//Chunking variables
			File chunkRunDir = new File(runDirectory,String.valueOf(i));
			chunkRunDir.mkdir();
			
			File chunkIntervals = new File(chunkRunDir,"interval.bed");
			ArrayList<BedEntry> intervalBed = new ArrayList<BedEntry>();
			int basesInChunk = 0;
			
			//Create chunk specific lists
			ArrayList<File> protectList = new ArrayList<File>(); //files to preserve on cleanup
			ArrayList<String> bamList = new ArrayList<String>(); //bam names 
			
			
			//Create interval based on chunk size
			while (true) {
				BedEntry workingEntry = bedFile.get(bedFileIndex);
				int entrySize = workingEntry.getEnd() - workingEntry.getStart();
				
				if (entrySize + basesInChunk <= basesPerChunk) {
					intervalBed.add(workingEntry);
					bedFileIndex+=1;
					basesInChunk += entrySize;
					if (bedFileIndex == bedFile.size()) {
						this.writeTargetRegions(intervalBed,chunkIntervals);
						break;
					}
				} else {
					int overhang = basesPerChunk - basesInChunk;
					int newBoundary = workingEntry.getStart() + overhang;
					BedEntry newbe = new BedEntry(workingEntry.getChrom(),workingEntry.getStart(),newBoundary);
					intervalBed.add(newbe);
					workingEntry.setStart(newBoundary);
					this.writeTargetRegions(intervalBed,chunkIntervals);
					break;
				}
			}
			
			//Create chunked bam files
			for (TFSampleInfo si: sampleList) {
				//Create file for original bams/split bams/output files
				File bam = new File(runDirectory,si.getSampleName() + ".bam");
				File chunkBam = new File(chunkRunDir,si.getSampleName() + "." + String.valueOf(i) + ".bam");
				File chunkBai = new File(chunkRunDir,si.getSampleName() + "." + String.valueOf(i) + ".bam.bai");
			
				//Run split
				this.logFile.writeInfoMessage("Splitting bam file (" + si.getSampleName() + ") by chunk (" + String.valueOf(i) + ")");
				this.splitByChunk(chunkIntervals, bam, chunkBam);
				
				//Add files to the delete list or cleanup protection list
				protectList.add(chunkBam);
				protectList.add(chunkBai);
				
				//Add files to cleanup list
				deleteList.add(chunkBam);
				deleteList.add(chunkBai);
				
				//Add bam name to list
				bamList.add(si.getSampleName() + "." + String.valueOf(i)  + ".bam");

			}
			
			//Create files for chunked vcf files
			File fullVcf = new File(chunkRunDir,this.study + ".vcf");
			File filterVcf = new File(chunkRunDir,this.study + ".filtered.vcf");
			File passVcf = new File(chunkRunDir,this.study + ".passing.vcf");
			File fullVcfIdx = new File(chunkRunDir,this.study + ".vcf.idx");
			File filterVcfIdx = new File(chunkRunDir,this.study + ".filtered.vcf.idx");
			File passVcfIdx = new File(chunkRunDir,this.study + ".passing.vcf.idx");
			
			//Add files to cleanup list
			deleteList.add(fullVcf);
			deleteList.add(filterVcf);
			deleteList.add(passVcf);
			deleteList.add(fullVcfIdx);
			deleteList.add(filterVcfIdx);
			deleteList.add(passVcfIdx);
			
			//Add vcf files to merge lists
			fullVcfList.add(fullVcf);
			filteredVcfList.add(filterVcf);
			passingVcfList.add(passVcf);
			
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
			
			//Create cmd.txt file
			File cmdFile = new File(chunkRunDir,"cmd.txt");
			protectList.add(cmdFile);
			
			//Create cmd.txt file
			this.createCmd(replacements, cmdFile);
			
			//Run
			TFThread thread = new TFThread(chunkRunDir,this.failmax, counter, this.heartbeat, protectList, this.logFile);
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
		
		cleanupRun(fullVcfList,filteredVcfList,passingVcfList,deleteList,runDirectory);
		
		
		//Fin
	}
	
	@Override
	public void run(ArrayList<TFSampleInfo> sampleList) {
		
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
		

	}
	
	
	private void cleanupRun(ArrayList<File> fullVcfList, ArrayList<File> filteredVcfList, ArrayList<File> passingVcfList, 
			ArrayList<File> deleteList, File runDirectory) {
			//Make destination directory
			File varDir = new File(this.rootDirectory,"Variants");
			varDir.mkdir();
			for(File file: varDir.listFiles()) file.delete();
			File jobDir = new File(varDir,"Jobs");
			
			//Make destination files
			File fullVcfDest = new File(varDir,this.study + ".vcf");
			File filterVcfDest = new File(varDir,this.study + ".filtered.vcf");
			File passingVcfDest = new File(varDir,this.study + ".passing.vcf");
			
			File fullVcfDestGz = new File(varDir,this.study + ".vcf.gz");
			File filterVcfDestGz = new File(varDir,this.study + ".filtered.vcf.gz");
			File passingVcfDestGz = new File(varDir,this.study + ".passing.vcf.gz");
			
			
			//Merge output files
			logFile.writeInfoMessage("Merging vcf files");
			this.mergeVcf(fullVcfList, fullVcfDest);
			this.mergeVcf(filteredVcfList,filterVcfDest);
			this.mergeVcf(passingVcfList, passingVcfDest);
			
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
			
			//Delete stragglers
			for (File f: deleteList) {
				this.deleteFile(f);
			}
			
			
			this.moveFile(runDirectory, jobDir);
				
	}
	
	private void writeTargetRegions(ArrayList<BedEntry> bedFile,File outFile) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			
			for (BedEntry be: bedFile) {
				bw.write(be.writeEntry());
			}
			
			bw.close();
		} catch (FileNotFoundException fnfe) {
			this.logFile.writeErrorMessage("[writeTargetRegions] Count not find target file", true);
		} catch (IOException ioex) {
			this.logFile.writeErrorMessage("[writeTargetRegions] Error reading target file",true);
		}
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
	
	
	private int countTargetSpace() {
		int count = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.targetFile));
		
			String line = null;
			
			while((line = br.readLine()) != null) {
				String[] items = line.split("\t");
				int chunk = Integer.parseInt(items[2]) - Integer.parseInt(items[1]);
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
	
	private int countReadsInTargets(TFSampleInfo si) {
		int n = -1;
		String fileName = si.getFile(TFConstants.FILE_BAM).getAbsolutePath();
		
		try {
			ProcessBuilder pb = new ProcessBuilder("/home/u0855942/applications/samtools/samtools","view","-c","-L",
					this.targetFile.getAbsolutePath(),fileName);
			
			Process p = pb.start();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String line = null;
			String result = "";
			while((line = br.readLine()) != null) {
				result += line;
				
			}
			
			int val = p.waitFor(); 
			
			n = Integer.parseInt(result);
					
			if (val != 0) {
				logFile.writeErrorMessage("[countReadsInTargets] Error when counting reads in bam file\n",true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[countReadsInTargers] IO Exception while trying to count reads in bam: " + fileName,true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[splitByChrom] Process was interrupted while trying count reads in bam: " + fileName,true);
			System.exit(1);
		}
		
		return n;
	}
	
	private void splitByChrom(String chrom, File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[splitByChrom] Expected file does not exist: " + source.getAbsolutePath(),true);
			}
			
			ProcessBuilder pb = new ProcessBuilder("/home/u0855942/applications/samtools/samtools","view","-b","-h","-@","8",source.getAbsolutePath(),chrom);
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
	
	private void splitByChunk(File intervals, File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[splitByChrom] Expected file does not exist: " + source.getAbsolutePath(),true);
			}
			
			ProcessBuilder pb = new ProcessBuilder("/home/u0855942/applications/samtools/samtools","view","-b","-h","-@","8",source.getAbsolutePath(),"-L",intervals.getAbsolutePath());
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
				logFile.writeErrorMessage("[splitByChunk] Error while splitting the bam file: " + intervals.getAbsolutePath() + " "+ source.getAbsolutePath(),true);
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
				logFile.writeErrorMessage("[splitByChunk] Error while indexing the bam file: " + intervals.getAbsolutePath() + " "+ source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[splitByChunk] IO Exception while trying to split the bam file: " + intervals.getAbsolutePath() + " " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[splitByChunk] Process was interrupted while trying to split bam files: " + intervals.getAbsolutePath() + " " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
	
	public int getReadLength(File bamFile){
		SAMFileReader samReader = null;
		
		int readLength = 0;
		int counter = 0;
			
		samReader = new SAMFileReader(bamFile);

		SAMRecordIterator it = samReader.iterator();

		while (it.hasNext()) {
			SAMRecord sam = it.next();
			int len = sam.getReadString().length();
			if (len > readLength) {
				readLength = len;
			}
			
			counter++;
			if (counter > 10000) {
				break;
			}
		}
			
		samReader.close();
	
		return readLength;
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
		
		private String[] chroms = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                "chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                "chrX","chrY","chrM"};

		@Override
		public int compare(BedEntry be1, BedEntry be2) {
			// TODO Auto-generated method stub
			
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
