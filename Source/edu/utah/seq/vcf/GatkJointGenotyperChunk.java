package edu.utah.seq.vcf;
import java.io.File;
import java.util.ArrayList;
import util.gen.IO;
import util.gen.Misc;

public class GatkJointGenotyperChunk implements Runnable{

	//fields
	private boolean complete = false;
	private boolean failed = false;
	private GatkJointGenotyper gatkRunner;
	private String threadId;
	private String[] cmd = null;
	private File tempVcf;
	private File tempBed;
	private File tempGdbLog;
	private File tempJGLog;
	private File tempChrGdb;
	private File tmpDirThread;
	private ArrayList<File> vcfs = new ArrayList<File>();
	private int numRetries = 3;


	public GatkJointGenotyperChunk(GatkJointGenotyper gatkRunner, String uniqueName) throws Exception{
		this.gatkRunner = gatkRunner;
		threadId = uniqueName;
	}

	public void run(){

		try {

			while ((tempBed = gatkRunner.getNextBed())!= null){

				long startTime = System.currentTimeMillis();
	
				String name = tempBed.getName().replace(".bed", "");

				//make temp vcf
				tempVcf = new File (gatkRunner.getTmpDirectory(), name+"_temp.vcf");

				//make log files
				tempGdbLog = new File (gatkRunner.getTmpDirectory(), name+"_temp.gdb.log");
				tempJGLog = new File (gatkRunner.getTmpDirectory(), name+"_temp.jg.log");

				//make chr_gdb
				tempChrGdb = new File (gatkRunner.getTmpDirectory(), name+"_temp.gdb");

				//make tmpDir for this thread
				tmpDirThread = new File (gatkRunner.getTmpDirectory(), name+"_temp.tmp");
				tmpDirThread.mkdirs();

				int gbRam = (int)Math.round(gatkRunner.getNumberGBPerChunk());

				//build and execute GenomicsDBImport call
				cmd = new String[]{
						gatkRunner.getGatkExecutable().getCanonicalPath(), 
						"--java-options", "-Xmx"+gbRam+"G",
						"GenomicsDBImport",
						"--genomicsdb-workspace-path", tempChrGdb.getCanonicalPath(),
						"--reference", gatkRunner.getFasta().getCanonicalPath(),
						"--sample-name-map", gatkRunner.getSampleMap().getCanonicalPath(),
						"--intervals", tempBed.getCanonicalPath(),
						"--tmp-dir", tmpDirThread.getCanonicalPath(),
						"--batch-size", "50",
						"--genomicsdb-vcf-buffer-size", "163840",
						"--merge-input-intervals", "true",
						"--bypass-feature-reader", "true"
				};
				int exitCode = 1;
				for (int i=1; i<= numRetries; i++) {
					IO.pl(name+"."+i+" Launching GenomicsDBImport "+Misc.stringArrayToString(cmd, " "));
					exitCode = IO.executeViaProcessBuilderReturnExit(cmd, tempGdbLog);
					if (exitCode == 0) break;
					else {
						IO.deleteDirectory(tempChrGdb);
						tempChrGdb = new File (gatkRunner.getTmpDirectory(), name+"_temp.gdb");
					}
				}

				//finish and calc run time
				double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
				if (exitCode == 0) IO.pl(name+" GenomicsDBImport DONE "+Math.round(diffTime)+" Min");
				else {
					IO.el(name+" GenomicsDBImport FAILED "+Math.round(diffTime)+" Min");
					complete = false;
					failed = true;
					gatkRunner.setKeepRunning(false);
					return;
				}

				//build and execute GenomicsDBImport call
				cmd = new String[]{
						gatkRunner.getGatkExecutable().getCanonicalPath(), 
						"--java-options", "-Xmx"+gbRam+"G",
						"GenotypeGVCFs",
						"--reference", gatkRunner.getFasta().getCanonicalPath(),
						"--tmp-dir", tmpDirThread.getCanonicalPath(),
						"-V", "gendb://"+tempChrGdb.getCanonicalPath(),
						"--allow-old-rms-mapping-quality-annotation-data",
						"-O", tempVcf.getCanonicalPath()
				};
				
				for (int i=1; i<= numRetries; i++) {
					IO.pl(name+"."+i+" Launching GenotypeGVCFs "+Misc.stringArrayToString(cmd, " "));
					exitCode = IO.executeViaProcessBuilderReturnExit(cmd, tempJGLog);
					if (exitCode == 0) break;
					else {
						tempVcf.delete();
						tempVcf = new File (gatkRunner.getTmpDirectory(), name+"_temp.vcf");
					}
				}

				//finish and calc run time
				diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
				if (exitCode == 0) {
					complete = true;
					IO.deleteDirectory(tempChrGdb);
					IO.deleteDirectory(tmpDirThread);
					IO.pl(name+" GenotypeGVCFs DONE "+Math.round(diffTime)+" Min");
					vcfs.add(tempVcf);
				}
				else {
					IO.el(name+" GenotypeGVCFs FAILED "+Math.round(diffTime)+" Min");
					complete = false;
					failed = true;
					gatkRunner.setKeepRunning(false);
					return;
				}
			}
		} catch (Exception e) {
			IO.el("ERROR in GatkJointGenotyperChunk "+threadId+" with "+tempBed);
			e.printStackTrace();
			complete = false;
			failed = true;
		}
	}

	//getters and setters
	public String getCmd() {
		return Misc.stringArrayToString(cmd, " ");
	}

	public boolean isComplete() {
		return complete;
	}

	public boolean isFailed() {
		return failed;
	}

	public ArrayList<File> getVcfs() {
		return vcfs;
	}
}
