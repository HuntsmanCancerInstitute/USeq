package util.apps;
import java.io.*;
import java.util.*;

import util.gen.IO;

public class FileMover {
	
	//fields
	private File destinationRoot;
	private File uploadDirectory;
	private int checkIntervalMinutes = 1;
	private long maximumAgeMilliseconds = 24 * 60 * 60 * 1000;
	
	//methods
	public void scan(){
		while (true){
			//get current time
			long now = System.currentTimeMillis();
			//look for directories in the uploadDirectory
			File[] directories = IO.extractOnlyDirectories(uploadDirectory);
			//for each directory
			for (File dir: directories){
				//is it older than max age?
				long lastModified = dir.lastModified();
				long diff = now - lastModified;
				if (diff > maximumAgeMilliseconds) IO.deleteDirectory(dir);
				else {
					//any files within?
					File[] files = IO.extractOnlyFiles(dir);
					if (files != null){
						//look for matching destination dir
						
					}
				}
				
			}
		}
	}
	
}
