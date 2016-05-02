package edu.utah.seq.vcf;
import java.io.File;
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
	public String getName() {
		return name;
	}

	private String cmd;
	public String getCmd() {
		return cmd;
	}

	private File tempVcf;
	private File tempBed;
	
	public GatkRunnerChunk(Bed[] regions, GatkRunner gatkRunner, String uniqueName) throws Exception{
		this.gatkRunner = gatkRunner;
		this.regions = regions;
		name = uniqueName;
	}
	
	public void run(){
		//write out bed
		tempBed = new File (gatkRunner.getSaveDirectory(), name+"_temp.bed");
		Bed.writeToFile(regions, tempBed);

		//make temp vcf
		tempVcf = new File (gatkRunner.getSaveDirectory(), name+"_temp.vcf");

		//build and execute gatk call
		cmd = gatkRunner.getGatkArgs()+" -L "+tempBed+" -o "+tempVcf;
		System.out.println(cmd);
		String[] out = IO.executeViaProcessBuilder(Misc.WHITESPACE.split(cmd), false);
		
		//any errors?
		for (String line: out) {
			if (line.toLowerCase().contains("error")){
				Misc.printArray(out);
				System.err.println("\n\nGATK threw an error when executing:\n"+cmd);
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
