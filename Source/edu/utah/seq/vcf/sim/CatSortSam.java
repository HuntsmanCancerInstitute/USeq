package edu.utah.seq.vcf.sim;

import java.io.File;
import java.util.ArrayList;

import edu.utah.seq.data.sam.PicardSortSam;
import util.gen.IO;
import util.gen.Misc;

public class CatSortSam implements Runnable {

	private boolean failed = false;
	private ArrayList<File> gzSamsToMerge;
	private File finalBamFile;
	private int id;
	
	public CatSortSam (ArrayList<File> gzSamsToMerge, File finalBamFile, int id) {
		this.gzSamsToMerge = gzSamsToMerge;
		this.finalBamFile = finalBamFile;
		this.id = id;
	}
	
	public void run() {	
		try {
			IO.pl("\tGenerating "+finalBamFile.getName());
			File finalSam = new File (finalBamFile.getParentFile(), "temp_"+id+"_"+Misc.removeExtension(finalBamFile.getName())+".sam.gz");
			finalSam.deleteOnExit();
			IO.concatinateFiles(gzSamsToMerge, finalSam);
			new PicardSortSam (finalSam, finalBamFile);
		} catch (Exception e) {
			failed = true;
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem concatinating sams or sorting sams ");
		} 
	}

	public boolean isFailed() {
		return failed;
	}

}
