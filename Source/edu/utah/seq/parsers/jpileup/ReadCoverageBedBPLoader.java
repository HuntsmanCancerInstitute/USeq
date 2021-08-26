package edu.utah.seq.parsers.jpileup;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import util.bio.annotation.Bed;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;

public class ReadCoverageBedBPLoader implements Runnable { 

	//fields
	private boolean failed = false;
	private int minimumReadDepth = 0;
	private PrintWriter bedOut= null;
	private Bed[] regions = null;
	private File coverageFile = null;
	private Histogram histogram = null;
	private BamPileupTabixLoaderSingle loader = null;


	public ReadCoverageBedBPLoader (File bamPileup, File tempDir, int maxCoverageCalculated, int minimumReadDepth, int loaderIndex, Bed[] regions) throws IOException{
		this.minimumReadDepth = minimumReadDepth;
		this.regions = regions;

		//create writers
		coverageFile = new File(tempDir, Misc.getRandomString(10)+"_"+loaderIndex+"_temp.bed");
		coverageFile.deleteOnExit();
		bedOut = new PrintWriter (new FileWriter(coverageFile));

		//create tabix reader
		loader = new BamPileupTabixLoaderSingle(bamPileup, 0);

		//create histogram to count read depth
		histogram = new Histogram(0, maxCoverageCalculated, maxCoverageCalculated);
	}

	public void run() {	
		Bed region = null;
		try {
			//for each region
			int counter = 0;	
			for (int x=0; x< regions.length; x++) {
				region = regions[x];
				if (counter++ > 1000) {
					counter = 0;
					IO.p(".");
				}
				//for each found bpileup line
				ArrayList<BpileupLine> lines = loader.fetchBpileupRecords(region);
				for (BpileupLine line: lines) {
					BaseCount bc = line.getSamples()[0];
					int cov = (int)bc.getPassingReadCoverage();
					histogram.count(cov);
					if (cov >= minimumReadDepth) bedOut.println(line.getBed());
				}
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem processing bam pileup lines for "+region);
			e.printStackTrace();
		} finally {
			bedOut.close();
			loader.getTabixReader().close();
		}
	}

	public boolean isFailed() {
		return failed;
	}

	public File getCoverageFile() {
		return coverageFile;
	}

	public Histogram getHistogram() {
		return histogram;
	}
}
