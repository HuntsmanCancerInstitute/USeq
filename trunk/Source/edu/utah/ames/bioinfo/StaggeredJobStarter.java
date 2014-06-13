package edu.utah.ames.bioinfo;

import java.util.Date;
import java.util.Timer;
import java.util.TimerTask;

public class StaggeredJobStarter {

	public static void main(String[] args) {
		StaggeredJobStarter executingTask = new StaggeredJobStarter();
		executingTask.start();
	}
	
	long delay = 10*1000; //delay in ms : 10 * 1000 ms = 10 sec.
	LoopTask task = new LoopTask();
	Timer timer = new Timer("TaskName");
	
	public void start() {
		timer.cancel();
		timer = new Timer("TaskName");
		Date executionDate = new Date(); // no params = now
		timer.scheduleAtFixedRate(task, executionDate, delay);
	}
	
	private class LoopTask extends TimerTask {
		
		public void run() {
			System.out.println("This message will print every 10 seconds");
		}
	}
}
