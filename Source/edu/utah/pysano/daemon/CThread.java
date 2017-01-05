package edu.utah.pysano.daemon;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Properties;
import java.util.concurrent.Callable;
import java.util.regex.Pattern;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

import edu.utah.pysano.utils.Logger;

public class CThread implements Callable<Boolean> {
		private File logFile = null;
		private Logger tfLogger = null;
		private File startFile = null;
		private File errorFile = null;
		private File stdoutFile = null;
		private File endFile = null;
		private File tmpFailFile = null;
		private File running = null;
		private File pbsFile = null;
		private static volatile Integer failCount = 0;
		private int failMax = 0;
		private int threadnumber = 0;
		private int heartbeat = 0;
		private File directory = null;
		private ArrayList<File> keepers = null;
		private ArrayList<File> deleteEarly = null;
		private String email = "";
		
		public static void setFailCount(int fc) {
			failCount = fc;
		}
		
		
		public CThread(File directory, int failMax, int threadnumber, int heartbeat, ArrayList<File> keepers, Logger tfLogger, String email) {
			//this.logFile = new File(directory,"log.txt");
			this.tmpFailFile = new File(directory,"tmp.fail");
			this.errorFile = new File(directory,"stderr.txt");
			//this.stdoutFile = new File(directory,"stdout.txt");
			
			this.logFile = new File(directory,"stdout.txt");
			//this.errorFile = new File(directory,"stdout.txt");
			
			this.pbsFile = new File(directory,"pbs.sh");
			this.startFile = new File(directory,"b");
			this.endFile = new File(directory,"a");
			this.running = new File(directory,"r");
			this.heartbeat = heartbeat;
			this.failMax = failMax;
			this.threadnumber = threadnumber;
			this.directory = directory;
			this.keepers = keepers;
			this.tfLogger = tfLogger;
			this.keepers.add(this.running);
			this.keepers.add(this.endFile);
			this.email = email;
		}
		
		private void writeErrorMessage(String message, boolean internal) {
			tfLogger.writeErrorMessage("[T" + this.threadnumber + "] " + message,internal);
		}
		private void writeInfoMessage(String message) {
			tfLogger.writeInfoMessage("[T" + this.threadnumber + "] " + message);
		}
		
		private void writeWarningMessage(String message) {
			tfLogger.writeWarningMessage("[T" + this.threadnumber + "] " + message);
		}
		
		public void cleanDirectory() {
			this.writeInfoMessage("Cleaning directory.");
			
			File[] fileList = this.directory.listFiles();
			for (File file: fileList) {
				if (!file.isDirectory() && !this.keepers.contains(file)) {
					file.delete();
				}
			}
			if (endFile.exists()) {
				endFile.delete();
			}
		}
		
		private void transferExisting() {
			this.writeInfoMessage("Backing up logs.");
			
			DateFormat dateFormat = new SimpleDateFormat("yyyyMMddHHmmss");
			Date date = new Date();
			String logName = "logs_" + dateFormat.format(date);
			
			File destDir = new File(this.directory,logName);
			
			
			if (this.pbsFile.exists()) {
				if (!destDir.exists()) {
					destDir.mkdir();
				}
				moveFile(pbsFile,new File(destDir,"pbs_old.sh"));
			}
			
			if (this.errorFile.exists()) {
				if (!destDir.exists()) {
					destDir.mkdir();
				}
				moveFile(errorFile,new File(destDir,"stderr.txt"));
			}
			
			if (this.logFile.exists()) {
				if (!destDir.exists()) {
					destDir.mkdir();
				}
				moveFile(logFile,new File(destDir,"stdout.txt"));
			}
			
			this.cleanDirectory();
			
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
						this.writeErrorMessage("Could not send pysano start signal.",true);
						analysisFailed = true;
						break;
					} catch (InterruptedException iex2) {
						this.writeErrorMessage("Start signal interrupted",true);
						analysisFailed = true;
						break;
					}
					
					boolean foundStderr = false; //If we find the stderr file
					boolean foundR = false;
					int heartbeat = 0;
					
					while(!foundStderr) {
						boolean tomatoFailed = false;
						Thread.sleep(5000);
						heartbeat++;
						if (heartbeat % (this.heartbeat * 60 / 5) == 0) {
							this.writeInfoMessage("Heartbeat: " + Math.round((float)heartbeat / 12) + " minutes.");
						}
						
											
						if (logFile.exists() && !running.exists()) {
							foundStderr = true;
							
							this.writeInfoMessage("Found log file!");
							Thread.sleep(1000);
							
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
								
								//if (Pattern.matches(".+?Job failed.*", lastLine)) {
								if (Pattern.matches("failure", lastLine)) {
									this.writeInfoMessage("Pysano job finished with errors");
									tomatoFailed = true;
								} 
								//else if (Pattern.matches(".+?Job completed successfully.*", lastLine)) {
								else if (Pattern.matches("success", lastLine)) {
									this.writeInfoMessage("Pysano job finished correctly");
									analysisFinished = true;
								} else {
									this.writeWarningMessage("Don't recognize the final line of the log file, assuming failed and resubmitting");
									tomatoFailed = true;
								}
							}
						} else if (this.tmpFailFile.exists() && !running.exists()) {
							foundStderr = true;
							tomatoFailed = true;
							this.writeInfoMessage("Pysano job finished with errors (tmp.fail found)");
						} else if (!running.exists() && !logFile.exists() && !errorFile.exists() && !startFile.exists() && !endFile.exists() && foundR) {
							foundStderr = true;
							this.writeInfoMessage("Pysano job directory has no control or output files, resubmitting job");
						}
						
						if (running.exists()) {
							foundR = true;
						}
//						} else if (this.logFile.exists() && !startFile.exists() && !running.exists() && !errorFile.exists()) {
//							BufferedReader br = new BufferedReader(new FileReader(logFile));
//							String line = "";
//							String lastLine = "";
//							while ((line=br.readLine()) != null) {
//								lastLine = line;
//							}
//							br.close();
//							
//							if (Pattern.matches(".+?Job failed.*", lastLine)) {
//								this.writeInfoMessage("Apparent tomato failure, resubmitting.");
//								foundStderr = true;
//							}
//						}
//						} else if (deleteEarly != null && logFile.exists() && !errorFile.exists() && !tmpFailFile.exists()) {
//							BufferedReader br = new BufferedReader(new FileReader(logFile));
//							String line = "";
//							String lastLine = "";
//							int counter = 0;
//							while ((line=br.readLine()) != null) {
//								counter += 1;
//								lastLine = line;
//							}
//							br.close();
//							
//							if (counter > 4 && !Pattern.matches(".+?Job failed.*", lastLine)) {
//								this.writeInfoMessage("T" + this.threadnumber + " Upload appears to be finished, deleting temp files");
//								File[] fileList = this.directory.listFiles();
//								for (File file: fileList) {
//									if (!file.isDirectory() && !this.deleteEarly.contains(file)) {
//										file.delete();
//									}
//								}
//							}
//							
//							this.deleted = true;
//							
//						}
							
							
						if (tomatoFailed) {
							failCount += 1;
							if (failCount > this.failMax) {
								this.writeErrorMessage("Too many Pysano run failures, killing all jobs and exiting", false);
								analysisFinished = true;
								analysisFailed = true;
							} 
//							else if (this.deleted) {
//								this.writeErrorMessage("Some necessary files were deleted to save space, and the job can't be resubmitted. "
//										+ "(I know this is a dumb-sounding error, but it was done to save space, sorry.",true);
//								analysisFailed = true;
//							} 
							else {
								this.writeInfoMessage("Resubmitting Pysano job.  Resubmission " + failCount + " of " + this.failMax);
								String message = "TF detected a command failure. The command will be re-submitted to pysano, but it might be worth checking the logs to determine the cause of the error.  There are " + (failMax - failCount) + " re-submits remaining of " + failMax + ".";
								try {
									CThread.postMail(email, "TF Command Failure", message);
								} catch (MessagingException e) {
									tfLogger.writeErrorMessage("[CThread] Failed to send ending email", true);
									e.printStackTrace();
								}
								
								

							}
						} 
					}
				}
			} catch (InterruptedException iex) {
				this.writeInfoMessage("Thread interrupted, killing job");
				try {
					Process p = Runtime.getRuntime().exec("touch " + endFile.getAbsolutePath());
					p.waitFor();
				} catch (IOException ioex) {
					this.writeErrorMessage("Could not send Pysano termination signal, please kill the jobs yourself",true);
				} catch (InterruptedException iex2) {
					this.writeErrorMessage("Termination signal interrupted, please kill jobs yourself",true);
				}
				analysisFailed = true;
			} catch (IOException ioex) {
				this.writeErrorMessage("Could not read the log file, please contact core",true);
				analysisFailed = true;
			}
			if (deleteEarly != null && !analysisFailed) {
				for (File f: deleteEarly) {
					if (f.exists()) {
						f.delete();
					}
				}
			}
			
			
			return analysisFailed;
		}
		
		public static void postMail(String recipients, String subject, String message) throws MessagingException {

			//set the host smtp address
			Properties props = new Properties();
			props.put("mail.smtp.host", "hci-mail.hci.utah.edu");

			//create some properties and get the default Session
			Session session = Session.getDefaultInstance(props, null);

			//create message
			Message msg = new MimeMessage(session);

			//set the from and to address
			InternetAddress addressFrom = new InternetAddress("noreply@hci.utah.edu");
			msg.setFrom(addressFrom);
			msg.setRecipients(Message.RecipientType.TO, InternetAddress.parse(recipients, false));

			//setting the Subject and Content type
			msg.setSubject(subject);
			msg.setContent(message, "text/plain");
			Transport.send(msg);
		}

}

