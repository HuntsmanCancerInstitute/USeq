package edu.utah.tomato.model;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;


	

public class TFFileObject {
	private String fileName = null;
	private File workingPath = null;
	private File finalPath = null;
	private File finalDirectory = null;
	private File workingDirectory = null;
	
	
	public TFFileObject(String fileName, File finalDirectory, File workingDirectory) {
		this.fileName = fileName;
		
		//Set up processing directory
		this.workingDirectory = workingDirectory;
		this.finalDirectory = finalDirectory;
		
		//Create paths
		this.workingPath = new File(workingDirectory, fileName);
		this.finalPath = new File(finalDirectory,fileName);
			
	}
	
//	public File createDestForFileObject(TFFileObject tffo) {
//		File linkDir = new File(tffo.workingDirectory,this.fileName);
//		return linkDir;
//	}
	
	public String getFileName() {
		return this.fileName;
	}
	
	
	public File createDestForFileObject(File directory) {
		File linkDir = new File(directory, this.fileName);
		return linkDir;
	}
	
	public File createDestForFileObject(File directory, String[] search, String[] replace) {
		String linkName = this.fileName;
		for (int i=0;i<search.length;i++) {
			linkName = linkName.replace(search[i], replace[i]);
		}
		
		File linkDir = new File(directory, linkName);
		return linkDir;
	}
	
	public boolean doesFinalExist() {
		if (finalPath != null && finalPath.exists()) {
			return true;
		} else {
			return false;
		}
	}
	
	public boolean doesWorkingExist() {
		if (workingPath != null && workingPath.exists()) {
			return true;
		} else {
			return false;
		}
	}
	
	
	public File getFinalDirectory() {
		return this.finalDirectory;
	}

	public File getFinalPath() {
		return this.finalPath;
	}
	
	public File getWorkingDirectory() {
		return this.workingDirectory;
	}
	
	public File getWorkingPath() {
		return this.workingPath;
	}
	
}
