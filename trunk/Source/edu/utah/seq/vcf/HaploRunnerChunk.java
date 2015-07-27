package edu.utah.seq.vcf;
import java.io.File;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class HaploRunnerChunk extends Thread{
	
	//fields
	private boolean complete = false;
	private boolean failed = false;
	private HaploRunner haploRunner;
	private Bed[] regions;
	private String name;
	private String cmd;
	private File tempVcf;
	private File tempBed;
	
	public HaploRunnerChunk(Bed[] regions, HaploRunner haploRunner, String uniqueName) throws Exception{
		this.haploRunner = haploRunner;
		this.regions = regions;
		name = uniqueName;
		start();
	}
	
	public void run(){
		
		//write out bed
		tempBed = new File (haploRunner.getSaveDirectory(), name+"_temp.bed");
		Bed.writeToFile(regions, tempBed);

		//make temp vcf
		tempVcf = new File (haploRunner.getSaveDirectory(), name+"_temp.vcf");

		//build and execute haplo call
		cmd = haploRunner.getHaploArgs()+" -L "+tempBed+" -o "+tempVcf;
		String[] out = IO.executeCommandLineReturnAll(Misc.WHITESPACE.split(cmd));
		System.out.println(cmd);
		
		//any errors?
		for (String line: out) {
			if (line.toLowerCase().contains("error")){
				Misc.printArray(out);
				System.err.println("\nGATK threw an error when executing:\n"+cmd);
				failed = true;
				break;
			}
		}

		//finish
		complete = true;
	}
	
	//getters and setters
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
}
