package edu.utah.tomato;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class TFThreadDaemon extends Thread {
	
	private TFLogger logFile = null;
	private List<Callable<Boolean>> taskList = null; 
	private ExecutorService pool = null;
	private ExecutorCompletionService<Boolean> service;
	private String commandString = null;
	private Integer totalJobs = null;
	private boolean failed = false;

	/** TFThread Daemon creates a Thread derived class that manages individual tomato job threads.  Once the daemon is
	 * created and started, you can give the daemon new jobs to manage.  Note that the daemon is started with a set 
	 * number of jobs and will terminate once that limit is reached.  The daemon will also terminate if there are too many
	 * job errors
	 * 
	 * @param logFile TFLogger file, used for messaging
	 * @param commandString Name of command, only used for logging
	 * @param totalJobs Number of jobs the daemon will handle.  Daemon exits once this is reached.
	 * @param concurrentJobs Number of concurrent jobs daemon will handle.  this limits cluster workload
	 */
	public TFThreadDaemon(TFLogger logFile, String commandString, Integer totalJobs, Integer concurrentJobs) {
		this.logFile = logFile;
		this.commandString = commandString;
		this.totalJobs = totalJobs;
		this.taskList = new ArrayList<Callable<Boolean>>();
		this.pool = Executors.newFixedThreadPool(concurrentJobs);
		this.service = new ExecutorCompletionService<Boolean>(pool);
		
		
	}
	
	private void writeErrorMessage(String message, boolean internal) {
		logFile.writeErrorMessage("[DAEMON] " + message,internal);
	}
	private void writeInfoMessage(String message) {
		logFile.writeInfoMessage("[DAEMON] " + message);
	}
	
	
	public void addJob(Callable<Boolean> job) {
		taskList.add(job);
	}
	
	public boolean getFailed() {
		return this.failed;
	}

	
	@Override
	public void run() {
		try {
			int finishedJobs = 0;
			boolean failed = false;
			this.writeInfoMessage("Creating daemon for TomatoFarmer command: " + this.commandString);
			while(true) {
				if (taskList.size() > 0) {
					Callable<Boolean> job = taskList.remove(0);
					service.submit(job);
				}
				
				Future<Boolean> result = service.poll();
				if (result != null) {
					finishedJobs++;
					try {
						if (result.get() == true && failed == false) {
							this.failed = true;
							Thread.currentThread().interrupt();
						}
					} catch (ExecutionException ex) {
						this.writeErrorMessage("Execution exception while waiting for jobs to finish",true);
						System.out.println(ex.getMessage());
						ex.printStackTrace();
						System.exit(1);
					} 
				}
				
				if (finishedJobs == this.totalJobs) {
					Thread.currentThread().interrupt();
				}
				
				if (Thread.interrupted()) {
				    throw new InterruptedException();
				}
				
				Thread.sleep(5000);
			}
					
		} catch (InterruptedException irrex) {
			this.writeInfoMessage("Shutting down job pool");
			this.pool.shutdownNow();
			try {
				boolean retval = pool.awaitTermination(15, java.util.concurrent.TimeUnit.SECONDS);
				if (!retval) {
					this.writeErrorMessage("Time limit reached, exiting without thread termination",true);
				}
			} catch (InterruptedException iex) {
				this.writeErrorMessage("Interrupted while waiting for termination signal, exiting without thread termination", true);
			}
			this.writeInfoMessage("Pool terminated");
		}
	}


}
