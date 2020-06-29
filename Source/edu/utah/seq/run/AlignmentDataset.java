package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;

/**Info related to a directory of alignment files.*/
public class AlignmentDataset{

	private boolean dirExists = false;
	private boolean complete = false;
	private boolean rna = false;
	private File bamFile = null;
	private File bamIndexFile = null;
	private File bedFile = null;
	private File gVcfFile = null;
	private File gVcfIndexFile = null;
	private ArrayList<String> info = null;

	public AlignmentDataset (File alignDir, ArrayList<String> info, boolean rna) throws IOException{
		this.info = info;
		this.rna = rna;
		if (alignDir.exists()) {
			dirExists = true;
			loadAlignmentFiles(alignDir);
		}
		else info.add("\t"+alignDir.getName()+" not found");
		
	}
	
	/**Looks for and returns the required files in the alignment dir, returns null if not found.
	 * Looks for the COMPLETE, .bam, _Pass.bed.gz, _Haplo.g.vcf.gz, and it's index
	 * @throws IOException */
	private void loadAlignmentFiles(File alignDir) throws IOException {
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(alignDir);
		if (nameFile.size() == 0) return;

		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final bam and bed
			File[] bam = IO.extractFiles(new File(alignDir, "Bam"), ".bam");
			File[] bamIndex = IO.extractFiles(new File(alignDir, "Bam"), ".bai");
			File[] passingBed = IO.extractFiles(new File(alignDir, "QC"), "_Pass.bed.gz");
			File[] gvcf = IO.extractFiles(new File(alignDir, "Vcfs"), "_Haplo.g.vcf.gz");
			File[] gvcfIndex = IO.extractFiles(new File(alignDir, "Vcfs"), "_Haplo.g.vcf.gz.tbi");
			
			//check bam
			String problem = null;
			if (bam == null || bam.length !=1) {
				File[] crams = IO.extractFiles(new File(alignDir, "Bam"), ".cram");
				if (crams != null && crams.length >0) {
					problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find the XXX.bam file in the Bam/ in "+
						alignDir+ ". Did find XXX.cram, convert this to bam, touch COMPLETE in the Bam/ dir, and restart.";
				}
				else problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find the XXX.bam file in the Bam/ in "+alignDir;
				info.add(problem);
			}
			else {
				//check for index
				if (bamIndex == null || bamIndex.length !=1) {
					problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find the XXX.bai index file in the Bam/ in "+alignDir;
					info.add(problem);
				}
				else {
					bamFile = bam[0];
					bamIndexFile = bamIndex[0];
				}
			}
			
			if (rna==false) {
				if (passingBed == null || passingBed.length !=1) {
					problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find the XXX_Pass.bed.gz file in QC/ in "+alignDir;
					info.add(problem);
				}
				else if (gvcf == null || gvcf.length !=1) {
					problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find the XXX_Haplo.g.vcf.gz file in Vcfs/ in "+alignDir;
					info.add(problem);
				}
				else if (gvcfIndex == null || gvcfIndex.length !=1) {
					problem = "\tERROR - The alignemnt was marked COMPLETE but failed to find the XXX_Haplo.g.vcf.gz.tbi index file in Vcfs/ in "+alignDir;
					info.add(problem);
				}
			}
			
			if (problem == null) {
				info.add("\tCOMPLETE "+alignDir);
				complete = true;
				if (rna == false) {
					bedFile = passingBed[0];
					gVcfFile = gvcf[0];
					gVcfIndexFile = gvcfIndex[0];
				}
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

	public File getBamFile() {
		return bamFile;
	}

	public File getBedFile() {
		return bedFile;
	}

	public File getgVcfFile() {
		return gVcfFile;
	}

	public File getgVcfIndexFile() {
		return gVcfIndexFile;
	}

	public File getBamIndexFile() {
		return bamIndexFile;
	}
}
