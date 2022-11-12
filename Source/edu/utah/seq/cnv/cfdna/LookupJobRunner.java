package edu.utah.seq.cnv.cfdna;

public class LookupJobRunner implements Runnable {

	//fields
	private LiquidBiopsyCA lbca;
	private boolean failed = false;
	private LookupJob[] lookupJob = null;

	public LookupJobRunner(LiquidBiopsyCA lbca) {
		this.lbca = lbca;
	}

	public void run() {	
		try {
			//get next chunk of work
			while ((lookupJob = lbca.fetchNextJob())!=null){
				for (LookupJob job: lookupJob) job.doWork();
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem running LookupJobs" );
			e.printStackTrace();
		}
	}

	public boolean isFailed() {
		return failed;
	}

}