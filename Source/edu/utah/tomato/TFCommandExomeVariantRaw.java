package edu.utah.tomato;

import java.io.BufferedInputStream;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;


public class TFCommandExomeVariantRaw extends TFCommand {
	
	private boolean splitChrom = false;
	private String study = null;
	private String[] hg19Chrom = {"chr1","chr2","chr3","chr4","chr5","chr6",
			"chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
			"chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
			"chrX","chrY","chrM"};
	private File targetFile;

	
	
	
	public TFCommandExomeVariantRaw(File templateFile, File rootDirectory,
			String commandString, String commandType, TFLogger logFile,
			String email, Integer wallTime, Integer heartbeat, Integer failmax,
			Integer jobs, boolean suppress, String study, boolean splitChrom, File targetFile) {
		super(templateFile, rootDirectory, commandString, commandType, logFile, email,
				wallTime, heartbeat, failmax, jobs, suppress);
		this.splitChrom = splitChrom;
		this.study = study;
		this.targetFile = targetFile;
	}

	
	
	
	private void runByChrom(ArrayList<TFSampleInfo> sampleList, File runDirectory) {
		this.daemon = new TFThreadDaemon(this.logFile, this.commandString, hg19Chrom.length, this.jobs);
		this.daemon.start();
		
		int counter = 1;
		
		//Initialize run specific lists
		ArrayList<File> deleteList = new ArrayList<File>();
		ArrayList<File>	fullVcfList = new ArrayList<File>();
		ArrayList<File>	filteredVcfList = new ArrayList<File>();
		ArrayList<File>	passingVcfList = new ArrayList<File>();
		
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
				this.logFile.writeInfoMessage("Splitting vcf file (" + si.getSampleName() + ") by chromosome (" + chrom + ")");
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
			File fullVcfGz = new File(chromRunDir,this.study + ".vcf.gz");
			File filterVcf = new File(chromRunDir,this.study + ".filtered.vcf");
			File filterVcfGz = new File(chromRunDir,this.study + ".filtered.vcf.gz");
			File passVcf = new File(chromRunDir,this.study + ".passing.vcf");
			File passVcfGz = new File(chromRunDir,this.study + ".passing.vcf.gz");
		
			//Add files to cleanup list
			deleteList.add(fullVcf);
			deleteList.add(fullVcfGz);
			deleteList.add(filterVcf);
			deleteList.add(filterVcfGz);
			deleteList.add(passVcf);
			deleteList.add(passVcfGz);
			
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
		
		//Make destination directory
		File varDir = new File(this.rootDirectory,"variants");
		varDir.mkdir();
		
		//Make destination files
		File fullVcfDest = new File(varDir,this.study + ".full.vcf");
		File fullVcfDestGz = new File(varDir,this.study + ".full.vcf.gz");
		File filterVcfDest = new File(varDir,this.study + ".filtered.vcf");
		File filterVcfDestGz = new File(varDir,this.study + ".filtered.vcf.gz");
		File passingVcfDest = new File(varDir,this.study + ".passing.vcf");
		File passingVcfDestGz = new File(varDir,this.study + ".passing.vcf.gz");
		
		//Merge output files
		this.mergeVcf(fullVcfList, fullVcfDest);
		this.mergeVcf(filteredVcfList,filterVcfDest);
		this.mergeVcf(passingVcfList, passingVcfDest);
		
		//Compress
		this.bgzip(fullVcfDest);
		this.bgzip(filterVcfDest);
		this.bgzip(passingVcfDest);
		
		//Index
		this.tabix(fullVcfDestGz);
		this.tabix(filterVcfDestGz);
		this.tabix(passingVcfDestGz);
		
		//Delete stragglers
		for (File f: deleteList) {
			if (f.exists()) {
				f.delete();
			}
		}
	}
	
	private void runFull(ArrayList<TFSampleInfo> sampleList, File runDirectory) {
		this.daemon = new TFThreadDaemon(this.logFile, this.commandString, 1, this.jobs);
		this.daemon.start();
		
		
	    int counter  = 1;
		ArrayList<File> keepers = new ArrayList<File>(); //files to preserve on cleanup
		ArrayList<String> bamList = new ArrayList<String>(); //bam names 
		
		
			for (TFSampleInfo si: sampleList) {
				//Create file for original bams/split bams/output files
				File bam = new File(runDirectory,si.getSampleName() + ".bam");
				File bai = new File(runDirectory,si.getSampleName() + ".bai");
			
				//Add bam/bai to preserve
				keepers.add(bam);
				keepers.add(bai);
				
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
			
			//Make Source files
			File fullVcf = new File(runDirectory,this.study + ".vcf");
			File fullVcfGz = new File(runDirectory,this.study + ".vcf.gz");
			File filterVcf = new File(runDirectory,this.study + ".filtered.vcf");
			File filterVcfGz = new File(runDirectory,this.study + ".filtered.vcf.gz");
			File passVcf = new File(runDirectory,this.study + ".passing.vcf");
			File passVcfGz = new File(runDirectory,this.study + ".passing.vcf.gz");
			
			//Make destination directory
			File varDir = new File(this.rootDirectory,"variants");
			varDir.mkdir();
			
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
			
			//Delete stragglers
			fullVcfGz.delete();
			filterVcfGz.delete();
			passVcfGz.delete();
			
			
			
			
	}

	@Override
	public void run(ArrayList<TFSampleInfo> sampleList) {
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
			
			//Create links
			this.createLink(si.getFile(TFConstants.FILE_BAM), bam);
			this.createLink(si.getFile(TFConstants.FILE_BAI), bai);
		}
		
		
		if (splitChrom) {
			runByChrom(sampleList, runDirectory);
		} else {
			runFull(sampleList, runDirectory);
		}
		
		
	}
	
	private void splitByChrom(String chrom, File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[splitByChrom] Expected file does not exist: " + source.getAbsolutePath(),true);
			}
			
			ProcessBuilder pb = new ProcessBuilder("/tomato/app/samtools/samtools","view","-b","-h",source.getAbsolutePath(),chrom);
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
	
	

}
