package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import util.gen.IO;

/**Info related to a directory of alignment files.*/
public class AlignmentDataset2{

	private boolean dirExists = false;
	private boolean complete = false;
	private boolean rna = false;
	private File cramFile = null;
	private File cramIndexFile = null;
	private File bamFile = null;
	private File bamIndexFile = null;
	private File bedFile = null;
	private File bpPileupFile = null;
	private File bpIndexFile = null;
	private ArrayList<String> info = null;

	public AlignmentDataset2 (File alignDir, ArrayList<String> info, boolean rna) throws IOException{
		this.info = info;
		this.rna = rna;
		if (alignDir.exists()) {
			dirExists = true;
			loadAlignmentFiles(alignDir);
		}
		else info.add("\t"+alignDir.getName()+" not found");
		
	}
	
	/**Looks for the COMPLETE, .cram, .bam, PassRC.bed.gz, bp.txt.gz, and their indexes
	 * @throws IOException */
	private void loadAlignmentFiles(File alignDir) throws IOException {
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(alignDir);
		if (nameFile.size() == 0) return;

		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			String dirName = "Bam";
			if (new File(alignDir, dirName).exists()== false) dirName = "Alignment";
			File[] cram = IO.extractFiles(new File(alignDir, dirName), ".cram");
			File[] cramIndex = IO.extractFiles(new File(alignDir, dirName), ".crai");
			File[] bam = IO.extractFiles(new File(alignDir, dirName), ".bam");
			File[] bamIndex = IO.extractFiles(new File(alignDir, dirName), ".bai");
			File[] bp = IO.extractFiles(new File(alignDir, dirName), "bp.txt.gz");
			File[] bpIndex = IO.extractFiles(new File(alignDir, dirName), "bp.txt.gz.tbi");
			
			File[] passingBed = IO.extractFiles(new File(alignDir, "QC"), "PassRC.bed.gz");
			if (passingBed == null || passingBed.length !=1) passingBed = IO.extractFiles(new File(alignDir, "QC"), "Pass.bed.gz");


			String problem = null;
			
			//load cram and crai, must be present
			if (cram == null || cram.length !=1) {
				problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find just one XXX.cram file in the Alignment/ in "+alignDir;
				info.add(problem);
			}
			if (cramIndex == null || cramIndex.length !=1) {
				problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find just one XXX.crai index file in the Alignment/ in "+alignDir;
				info.add(problem);	
			}
			if (problem == null) {
				cramFile = cram[0];
				cramIndexFile = cramIndex[0];
			}
			
			//load bam and bai if present, might not be available which is OK
			if (bam != null && bamIndex != null && bam.length ==1 && bamIndex.length == 1) {
				bamFile = bam[0];
				bamIndexFile = bamIndex[0];
			}
			
			//load bp.txt.gz and the tbi
			if (bp == null || bp.length !=1) {
				problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find just one bam pileup XXX.bp.txt.gz file in the Alignment/ in "+alignDir;
				info.add(problem);
			}
			if (bpIndex == null || bpIndex.length !=1) {
				problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find just one bam pileup XXX.bp.txt.gz.tbi index file in the Alignment/ in "+alignDir;
				info.add(problem);
			}
			if (problem == null) {
				bpPileupFile = bp[0];
				bpIndexFile = bpIndex[0];
			}
			
			//load passing bed for DNA
			if (rna == false) {
				if (passingBed == null || passingBed.length !=1) {
					problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find the XXX.PassRC.bed.gz file in QC/ in "+alignDir;
					info.add(problem);
				}
				else bedFile = passingBed[0];
			}
			
			if (problem == null) {
				info.add("\tCOMPLETE "+alignDir);
				complete = true;
			}
		}
	}
	
	public File[] checkNumberFiles(File dir, String extension, int requiredNumberFiles) throws IOException {
		File[] f = IO.extractFiles(dir, extension);
		if (f.length != requiredNumberFiles) {
			return null;
		}
		return f;
	}

	public boolean isDirExists() {
		return dirExists;
	}

	public boolean isComplete() {
		return complete;
	}

	public boolean isRna() {
		return rna;
	}

	public File getCramFile() {
		return cramFile;
	}

	public File getBedFile() {
		return bedFile;
	}

	public File getCramIndexFile() {
		return cramIndexFile;
	}

	public File getBpPileupFile() {
		return bpPileupFile;
	}

	public File getBpIndexFile() {
		return bpIndexFile;
	}

	public File getBamFile() {
		return bamFile;
	}

	public File getBamIndexFile() {
		return bamIndexFile;
	}
}
