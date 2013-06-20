package edu.utah.tomato;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class TFLogger {
	private BufferedWriter logWriter = null;
	private String logLevel = null;
	

	public TFLogger(File logDirectory, String logPrefix, String logLevel) {
		this.logLevel = logLevel;

		DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");
		DateFormat year = new SimpleDateFormat("yyyy");
		DateFormat month = new SimpleDateFormat("MM");
		DateFormat day = new SimpleDateFormat("dd");
		Date date = new Date();
		String logName = "log_" + dateFormat.format(date) + "_" + logPrefix + ".txt";
		
	    
	    
	    File logDir = new File(logDirectory,"logs/" + year.format(date) + "/" + month.format(date) + "/" + day.format(date));
	 
	    if (!logDir.exists()) {
	    	logDir.mkdirs();
	    }
	    
	    try {
	    	logWriter = new BufferedWriter(new FileWriter(new File(logDir,logName)));
	    } catch (IOException ioex) {
	    	writeErrorMessage("Could not open TomatoFarmer logging file",true);
	    }
	    
	}
	
	private void writeToLog(String message) {
		if (this.logWriter == null) {
			return;
		}
		try {
			this.logWriter.write(message + "\n");
		} catch (IOException ioex) {
			writeErrorMessage("Error writing to TomatoFarmer log file.  Exiting run, please contact core with the run log: bioinformaticscore@utah.edu\n",true);
			System.exit(1);
		}
	}
	
	public void writeInputMessage(String message) {
		String fullMessage = "INPUT: " + message + ": ";
		System.out.print(fullMessage);
		writeToLog(fullMessage);
	}
	
	public void writeSystemMessage(String message) {
		String fullMessage = "SYSTEM: " + message;
		writeToLog(fullMessage);
	}
	
	public void writeErrorMessage(String message, boolean internal) {
		
		if (internal) {
			String fullMessage = "ERROR (internal): " + message + ". Exiting run, please contact core with the run log: bioinformaticscore@utah.edu\n";
			System.out.println(fullMessage);
			writeToLog(fullMessage);
		} else {
			String fullMessage = "ERROR (user): " + message + ". Exiting run, please correct errors and re-submit. Contact the core with questions: bioinformaticscore@utah.edu\n";
			System.out.println(fullMessage);
			writeToLog(fullMessage);
		}
	}
	
	public void writeWarningMessage(String message) {
		String fullMessage = "WARNING: " + message;
		writeToLog(fullMessage);
		if (this.logLevel.equals("INFO") || this.logLevel.equals("WARNING")) {
			System.out.println(fullMessage);
		}
	}
	
	public void writeInfoMessage(String message) {
		String fullMessage = "INFO: " + message;
		writeToLog(fullMessage);
		if (this.logLevel.equals("INFO")) {
			System.out.println(fullMessage);
		}
	}
	
	public void closeLogger() {
		try {
			logWriter.close();
		} catch (IOException ioex) {
			writeErrorMessage("Failed to close TomatoFarmer logging file",true);
		}
	}

}
