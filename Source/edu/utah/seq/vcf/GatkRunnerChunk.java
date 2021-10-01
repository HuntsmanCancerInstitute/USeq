package edu.utah.seq.vcf;
import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.CopyOnWriteArrayList;

import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class GatkRunnerChunk implements Runnable{
	
	//fields
	private boolean complete = false;
	private boolean failed = false;
	private GatkRunner gatkRunner;
	private Bed[] regions;
	private String name;
	private String cmd;
	private File tempVcf;
	private File tempBed;
	private File tempLog;
	private File bamOut = null;
	private String oOption = " -o ";
	
	public GatkRunnerChunk(Bed[] regions, GatkRunner gatkRunner, String uniqueName) throws Exception{
		this.gatkRunner = gatkRunner;
		this.regions = regions;
		name = uniqueName;
		if (gatkRunner.isUseUpperCaseO()) oOption = " -O ";
	}
	
	public void run(){
		
		//write out bed
		tempBed = new File (gatkRunner.getSaveDirectory(), name+"_temp.bed");
		Bed.writeToFile(regions, tempBed);

		//make temp vcf
		tempVcf = new File (gatkRunner.getSaveDirectory(), name+"_temp.vcf");
		
		//make log file
		tempLog = new File (gatkRunner.getSaveDirectory(), name+"_temp.log");
		
		//output the bam?
		if (gatkRunner.bamOut) {
			bamOut = new File (gatkRunner.getSaveDirectory(), name+"_temp.bam");
			File bai = new File (gatkRunner.getSaveDirectory(), name+"_temp.bai");
			bamOut.deleteOnExit();
			bai.deleteOnExit();
		}

		//build and execute gatk call
		if (gatkRunner.useLowerCaseL) {
			String[] t = Misc.WHITESPACE.split(gatkRunner.getGatkArgs());
			t[t.length-1] = " -l "+tempBed+ oOption +tempVcf +" "+t[t.length-1];
			cmd = Misc.stringArrayToString(t, " ");
		}
		else cmd = gatkRunner.getGatkArgs()+ " -L "+tempBed+ oOption +tempVcf;
		if (bamOut != null) cmd = cmd + " -bamout "+bamOut;
		
		System.out.println(cmd);
		String[] out = IO.executeViaProcessBuilder(Misc.WHITESPACE.split(cmd), tempLog);
		
		//any errors? ignore StatusLogger errors
		for (String line: out) {
			if (line.toLowerCase().contains("error") && line.contains("StatusLogger") == false){
				Misc.printArray(out);
				System.err.println("\n\nAborting, GATK threw an error when executing:\n"+cmd+"\nErrorLine: "+line+"\n");
				failed = true;
				break;
			}
		}

		//finish
		complete = true;
	}
	
	//getters and setters
	public String getCmd() {
		return cmd;
	}
	public String getName() {
		return name;
	}
	
	public boolean isComplete() {
		return complete;
	}

	public boolean isFailed() {
		return failed;
	}

	public File getTempVcf() {
		return tempVcf;
	}

	public File getTempBed() {
		return tempBed;
	}

	public File getBamOut() {
		return bamOut;
	}

	public File getTempLog() {
		return tempLog;
	}
}
