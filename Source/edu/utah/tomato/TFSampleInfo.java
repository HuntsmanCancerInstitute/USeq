package edu.utah.tomato;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TFSampleInfo {

	private String sampleName = null;
	private String sampleID = null;
	private String puID = null;
	private boolean qual64 = false;
	private HashMap<String,File> fileList = null;
	private TFLogger logFile = null;
	private int recordCount = 0;
	private int recordTargetCount = 0;
	private int readLength = 0;
	
	public TFSampleInfo(String sampleID, TFLogger logFile) {
		this.sampleID = sampleID;
		this.logFile = logFile;
		this.fileList = new HashMap<String,File>();
		this.fileList.put(TFConstants.FILE_FASTQ1, null);
		this.fileList.put(TFConstants.FILE_FASTQ2, null);
		this.fileList.put(TFConstants.FILE_BAI, null);
		this.fileList.put(TFConstants.FILE_BAM, null);
		this.fileList.put(TFConstants.FILE_SAM, null);
		this.fileList.put(TFConstants.FILE_METRICS, null);
		this.fileList.put(TFConstants.FILE_LANE_BAM,null);
		this.fileList.put(TFConstants.FILE_LANE_BAI,null);
		this.fileList.put(TFConstants.FILE_SAMPLE_BAM, null);
		this.fileList.put(TFConstants.FILE_SAMPLE_BAI, null);
		this.fileList.put(TFConstants.FILE_SPLIT_LANE_BAM,null);
		this.fileList.put(TFConstants.FILE_SPLIT_LANE_BAI, null);
		this.fileList.put(TFConstants.FILE_REALIGN_SAMPLE_BAM,null);
		this.fileList.put(TFConstants.FILE_REALIGN_SAMPLE_BAI, null);
		this.fileList.put(TFConstants.FILE_REDUCE_BAM,null);
		this.fileList.put(TFConstants.FILE_REDUCE_BAI, null);
	
	}
	
	public static ArrayList<TFSampleInfo> identifyAndValidateSamples(File directory, String commandType, TFLogger logFile) {
		//Create containers for valid/ignored files
		ArrayList<TFSampleInfo> sampleList = new ArrayList<TFSampleInfo>();
		ArrayList<File> ignored = new ArrayList<File>();
		
		//Get patterns for the commandType
		ArrayList<StringPattern> patterns = TFSampleInfo.createPatterns(commandType, logFile);
		
		//Iterate through file in the directory
		File[] contents = directory.listFiles();
			
		for (File file: contents) {
			if (file.isDirectory()) {
				continue;
			}
			
			boolean found = false;
			
			for (StringPattern pattern: patterns) {
				Matcher m = pattern.getPattern().matcher(file.getName());
				if (m.matches()) {
					TFSampleInfo.matchFiles(file, sampleList, m, pattern, logFile);
					found = true;
				}
			}
			
			if(!found) {
				ignored.add(file);
			}
		}
		
		//Write out ignored files
		//for (File ignore: ignored) {
			//logFile.writeInfoMessage("Ignored file: " + ignore.getAbsolutePath());
		//}
		
		//Make sure each sample has the appropriate files
		for (TFSampleInfo sample: sampleList) {
			validateSample(sample,commandType,logFile);
		}
	
		return sampleList;
	
	}
	
	
	private static void matchFiles(File f, ArrayList<TFSampleInfo> sampleList, Matcher m, StringPattern p, TFLogger logFile) {
		TFSampleInfo used = null;
		for (TFSampleInfo sample: sampleList) {
			if (sample.comparePrefix(m.group(1))) {
				used = sample;
			}
		}
		if (used == null) {
			used = new TFSampleInfo(m.group(1),logFile);
			sampleList.add(used);
		}
		
		logFile.writeInfoMessage("Used file: " + used.getSampleID() + " " +  f.getAbsolutePath());
		used.setFile(p.getName(), f);
	}
	
	private static ArrayList<StringPattern> createPatterns(String commandType, TFLogger logFile) {
		ArrayList<StringPattern> patterns = new ArrayList<StringPattern>();
			
		if (commandType.equals(TFConstants.ANALYSIS_EXOME_ALIGN_BWA) || commandType.equals(TFConstants.ANALYSIS_EXOME_ALIGN_NOVO)) {
			patterns.add(new StringPattern(TFConstants.FILE_FASTQ1,Pattern.compile("(.+?)_1\\..+\\.gz$")));
			patterns.add(new StringPattern(TFConstants.FILE_FASTQ2,Pattern.compile("(.+?)_2\\..+\\.gz$")));
		} else if (commandType.equals(TFConstants.ANALYSIS_EXOME_METRICS)) {
			patterns.add(new StringPattern(TFConstants.FILE_SPLIT_LANE_BAM,Pattern.compile("(.+?)\\.split.bam$")));
			patterns.add(new StringPattern(TFConstants.FILE_SPLIT_LANE_BAI,Pattern.compile("(.+?)\\.split.bai$")));
			patterns.add(new StringPattern(TFConstants.FILE_BAM,Pattern.compile("(.+?)\\.raw\\.bam$")));
			patterns.add(new StringPattern(TFConstants.FILE_BAI,Pattern.compile("(.+?)\\.raw\\.bai$")));
		} else if (commandType.equals(TFConstants.ANALYSIS_EXOME_VARIANT_RAW) || commandType.equals(TFConstants.ANALYSIS_EXOME_VARIANT_VQSR)) {
			patterns.add(new StringPattern(TFConstants.FILE_REDUCE_BAM,Pattern.compile("(.+?)\\.reduce\\.bam$")));
			patterns.add(new StringPattern(TFConstants.FILE_REDUCE_BAI,Pattern.compile("(.+?)\\.reduce\\.bai$")));
		} else {
			logFile.writeErrorMessage("<TFSampleInfo> Command type is not recognized by sample info parser: " + commandType, true);
			System.exit(1);
		}
		return patterns;
	}
	
	
	private static void validateSample(TFSampleInfo si, String commandType, TFLogger logFile) {
		if (commandType.equals(TFConstants.ANALYSIS_EXOME_ALIGN_BWA) || commandType.equals(TFConstants.ANALYSIS_EXOME_ALIGN_NOVO)) {
		    if (si.fileList.get(TFConstants.FILE_FASTQ1) == null || si.fileList.get(TFConstants.FILE_FASTQ2) == null) {
		    	logFile.writeErrorMessage("<TFSampleInfo> Found only one paired-end file, please make sure both are in the directory.  Offending sample name: " + si.getSampleName(),false);
				System.exit(1);
		    }
		    String[] parts = si.getSampleID().split("_",2);
		    if (parts.length != 2) {
		    	logFile.writeErrorMessage("<TFSampleInfo> The Fastq file name is not properly formed. Should have an underscore separating the sample name and the flowcell name: SAMPLE_FLOWCELL_1.fastq and SAMPLE_FLOWCELL_2.fastq",false);
		    	System.exit(1);
		    }
		    si.setPuID(parts[1]);
		    si.setSampleName(parts[0]);
		    
		     
		} else if (commandType.equals(TFConstants.ANALYSIS_EXOME_METRICS)) {
			if (si.fileList.get(TFConstants.FILE_BAM) == null || si.fileList.get(TFConstants.FILE_BAI) == null || si.fileList.get(TFConstants.FILE_SPLIT_LANE_BAM) == null || 
					si.fileList.get(TFConstants.FILE_SPLIT_LANE_BAI) == null) {
				logFile.writeErrorMessage("<TFSampleInfo> Could not find bam/bai/sam set.  Please make sure that the bam/bai/sam combination is present for each prefix. Offending sample name: " + si.getSampleName(),false);
				System.exit(1);
			}
		} else if (commandType.equals(TFConstants.ANALYSIS_EXOME_VARIANT_RAW) || commandType.equals(TFConstants.ANALYSIS_EXOME_VARIANT_VQSR)) {
			if (si.fileList.get(TFConstants.FILE_REDUCE_BAM) == null || si.fileList.get(TFConstants.FILE_REDUCE_BAI) == null) {
				logFile.writeErrorMessage("<TFSampleInfo> Could not find bam/bai set.  Please make sure that the bam/bai combination is present for each prefix. Offending sample name: " + si.getSampleName(),false);
				System.exit(1);
			}
		} else {
			logFile.writeErrorMessage("<TFSampleInfo> Do not recognize analysis type: " + commandType,true);
			System.exit(1);
		}
	}
	
	public void setFile(String name, File file) {
		if (!this.fileList.containsKey(name)) {
			logFile.writeErrorMessage("<TFSampleInfo> File type does not exist in sample info object: " + name, true);
			System.exit(1);
		}	
		this.fileList.put(name, file);
	}
	
	public File getFile(String name) {
		if (!this.fileList.containsKey(name)) {
			logFile.writeErrorMessage("<TFSampleInfo> File type does not exist in sample info object: " + name, true);
			System.exit(1);
		}
		return this.fileList.get(name);
	}
		
	public void setQual64() {
		this.qual64 = true;
	}
		
	public boolean isQual64() {
		return this.qual64;
	}
		
	
	public String getSampleName() {
		return this.sampleName;
	}
		
	private boolean comparePrefix(String prefix) {
		if (this.sampleID.equals(prefix)) {
			return true;
		} else {
			return false;
		}
	}
	
	public void setRecordCount(int count) {
		this.recordCount = count;
	}
	
	public int getRecordCount() {
		return this.recordCount;
	}
	
	public void setRecordTargetCount(int count) {
		this.recordTargetCount = count;
	}
	
	public int getRecordTargetCount() {
		return this.recordTargetCount;
	}
	
	public int getReadLength() {
		return this.readLength;
	}
	
	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}
	
	public String getSampleID() {
		return sampleID;
	}

	public String getPuID() {
		return puID;
	}
	
	public void setPuID(String id) {
		this.puID = id;
	}
	
	public void setSampleName(String name) {
		this.sampleName = name;
	}

	
	

}
