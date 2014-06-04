package edu.utah.tomato.model;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.tomato.util.TFLogger;

public class TFSampleInfo {

	private String sampleName = null;
	private String sampleID = null;
	private String puID = null;
	private boolean qual64 = false;
	private HashMap<String,TFFileObject> fileList = null;
	private TFLogger logFile = null;
	private int recordCount = 0;
	private int recordTargetCount = 0;
	private int readLength = 0;
	

	public TFSampleInfo(String sampleID, String sampleName, String puID, TFLogger logFile) {
		this.sampleID = sampleID;
		this.sampleName = sampleName;
		this.puID = puID;
		this.logFile = logFile;
		this.fileList = new HashMap<String,TFFileObject>();
	}

	
	public void setFileObject(String name, TFFileObject fileObject) {
//		if (!this.fileList.containsKey(name)) {
//			logFile.writeErrorMessage("<TFSampleInfo> File type does not exist in sample info object: " + name, true);
//			System.exit(1);
//		}	
		this.fileList.put(name, fileObject);
	}
	
	public boolean termExists(String name) {
		if (this.fileList.get(name) != null) {
			return true;
		} else {
			return false;
		}
	}
	
	public boolean finalFileExists(String name) {
		if (termExists(name) && this.fileList.get(name).getFinalPath().exists()) {
			return true;
		} else {
			return false;
		}
	}
	
	public TFFileObject getFileObject(String name) {
		if (!this.fileList.containsKey(name)) {
			logFile.writeErrorMessage("[TFSampleInfo] File type does not exist in sample info object: " + name, true);
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
		
	public boolean comparePrefix(String prefix) {
		if (this.sampleID.equals("unknown")) {
			if (this.sampleName.equals(prefix)) {
				return true;
			} else {
				return false;
			}
		} else {
			if (this.sampleID.equals(prefix)) {
				return true;
			} else {
				return false;
			}
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
	
	public void removeFileObject(String name) {
		if (this.fileList.containsKey(name)) {
			this.fileList.remove(name);
		}
	}
	

	
	

}
