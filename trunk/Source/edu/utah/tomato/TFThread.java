package edu.utah.tomato;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.Callable;
import java.util.regex.Pattern;

public class TFThread implements Callable<Boolean> {
		private File logFile = null;
		private TFLogger tfLogger = null;
		private File startFile = null;
		private File errorFile = null;
		private File stdoutFile = null;
		private File endFile = null;
		private File tmpFailFile = null;
		private static volatile Integer failCount = 0;
		private int failMax = 0;
		private int threadnumber = 0;
		private int heartbeat = 0;
		private File directory = null;
		private ArrayList<File> keepers = null;
		
		public TFThread(File directory, int failMax, int threadnumber, int heartbeat, ArrayList<File> keepers, TFLogger tfLogger) {
			this.logFile = new File(directory,"log.txt");
			this.tmpFailFile = new File(directory,"tmp.fail");
			this.errorFile = new File(directory,"stderr.txt");
			this.stdoutFile = new File(directory,"stdout.txt");
			this.startFile = new File(directory,"b");
			this.endFile = new File(directory,"a");
			this.heartbeat = heartbeat;
			this.failMax = failMax;
			this.threadnumber = threadnumber;
			this.directory = directory;
			this.keepers = keepers;
			this.tfLogger = tfLogger;
		}
		
		private void writeErrorMessage(String message, boolean internal) {
			tfLogger.writeErrorMessage("[T" + this.threadnumber + "] " + message,internal);
		}
		private void writeInfoMessage(String message) {
			tfLogger.writeInfoMessage("[T" + this.threadnumber + "] " + message);
		}
		
		private void transferExisting() {
			this.writeInfoMessage("Cleaning directory.");
			
			DateFormat dateFormat = new SimpleDateFormat("yyyyMMddHHmmss");
			Date date = new Date();
			String logName = "logs_" + dateFormat.format(date);
			
			File destDir = new File(this.directory,logName);
			
			if (this.stdoutFile.exists()) {
				if (!destDir.exists()) {
					destDir.mkdir();
				}
				if (this.stdoutFile.exists()) {
					moveFile(stdoutFile,new File(destDir,"stdout.txt"));
				}
			}
			
			if (this.errorFile.exists()) {
				if (!destDir.exists()) {
					destDir.mkdir();
				}
				if (this.errorFile.exists()) {
					moveFile(errorFile,new File(destDir,"stderr.txt"));
				}
			}
			
			if (this.logFile.exists()) {
				if (!destDir.exists()) {
					destDir.mkdir();
				}
				if (this.logFile.exists()) {
					moveFile(logFile,new File(destDir,"log.txt"));
				}
			}
			
			
			File[] fileList = this.directory.listFiles();
			for (File file: fileList) {
				if (!file.isDirectory() && !this.keepers.contains(file)) {
					file.delete();
				}
			}
			
		}
		
		private void moveFile(File source, File dest) {
			try {
				if (!source.exists()) {
					this.writeErrorMessage("Expected file does not exist: " + source.getAbsolutePath(),true);
				}
				Process p = Runtime.getRuntime().exec("mv " + source.getAbsolutePath() + " " + dest.getAbsolutePath());
				int val = p.waitFor();
				if (val != 0) {
					this.writeErrorMessage("System could not move your alignment file " + source.getAbsolutePath(),true);
					System.exit(1);
				}
			} catch (IOException ioex) {
				this.writeErrorMessage("IO Exception while trying move alignment file: " + source.getAbsolutePath(),true);
				System.exit(1);
			} catch (InterruptedException ieex) {
				this.writeErrorMessage("Process was interrupted while trying to create a link to your fastq file: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
		}
		
		
		public Boolean call()  {
			
			
			boolean analysisFailed = false;
			boolean analysisFinished = false;
			try {
				while(!analysisFinished) {
					this.transferExisting();
					try {
						this.writeInfoMessage("Submitting job in directory: "  + this.directory.getAbsolutePath());
						Process p = Runtime.getRuntime().exec("touch " + startFile.toString());
						p.waitFor();
					} catch (IOException ioex) {
						this.writeErrorMessage("Could not send tomato start signal.",true);
						analysisFailed = true;
						break;
					} catch (InterruptedException iex2) {
						this.writeErrorMessage("Start signal interrupted",true);
						analysisFailed = true;
						break;
					}
					
					boolean foundStderr = false; //If we find the stderr file
					int heartbeat = 0;
					
					while(!foundStderr) {
						boolean tomatoFailed = false;
						Thread.sleep(5000);
						heartbeat++;
						if (heartbeat % (this.heartbeat * 60 / 5) == 0) {
							this.writeInfoMessage("Heartbeat: " + Math.round((float)heartbeat / 12) + " minutes.");
						}
						
											
						if (errorFile.exists()) {
							foundStderr = true;
							
							this.writeInfoMessage("Found log file!");
							Thread.sleep(5000);
							
							if (!this.logFile.exists()) {
								this.writeInfoMessage("T" + this.threadnumber + " Could not find log file, assuming error");
								tomatoFailed = true;
							} else {
								BufferedReader br = new BufferedReader(new FileReader(logFile));
								String line = "";
								String lastLine = "";
								while ((line=br.readLine()) != null) {
									lastLine = line;
								}
								br.close();
								
								if (Pattern.matches(".+?Job failed.*", lastLine)) {
									this.writeInfoMessage("Tomato job finished with errors");
									tomatoFailed = true;
								} else if (Pattern.matches(".+?Job completed successfully.*", lastLine)) {
									this.writeInfoMessage("Tomato job finished correctly");
									analysisFinished = true;
								} else {
									this.writeErrorMessage("Don't recognize the final line of the log file, assuming failed",true);
									tomatoFailed = true;
								}
							}
						} else if (this.tmpFailFile.exists()) {
							foundStderr = true;
							this.writeInfoMessage("Found tmp.fail, resubmitting");
						}
							
							
						if (tomatoFailed) {
							failCount += 1;
							if (failCount > this.failMax) {
								this.writeInfoMessage("Too many tomato run failures, killing all jobs and exiting");
								analysisFinished = true;
								analysisFailed = true;
							} else {
								this.writeInfoMessage("Resubmitting tomato job.  Resubmission " + failCount + " of " + this.failMax);

							}
						} 
					}
				}
			} catch (InterruptedException iex) {
				this.writeInfoMessage("Thread interrupted, killing job");
				try {
					Process p = Runtime.getRuntime().exec("touch " + endFile.toString());
					p.waitFor();
				} catch (IOException ioex) {
					this.writeErrorMessage("Could not send tomato termination signal, please kill the jobs yourself",true);
				} catch (InterruptedException iex2) {
					this.writeErrorMessage("Termination signal interrupted, please kill jobs yourself",true);
				}
				analysisFailed = true;
			} catch (IOException ioex) {
				this.writeErrorMessage("Could not read the log file, please contact core",true);
				analysisFailed = true;
			}
			return analysisFailed;
		}

}
