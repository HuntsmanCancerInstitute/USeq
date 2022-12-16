package edu.utah.seq.cnv.cfdna;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Data related to one cfDNA and matched germline sample set for all of the genes with their group of unique observations.*/
public class LiquidBiopsySample {
	
	//fields
	private String sampleName;
	private File cfAlignment;
	private File germlineAlignment;
	private LookupJob[] cfJobs = null;
	private SamReader cfReader = null;
	private LookupJob[] germlineJobs = null;
	private SamReader germlineReader = null;
	private File resultsDir = null;
	private Gzipper cfBedOut = null;
	private Gzipper germlineBedOut = null;
	
	//constructor
	public LiquidBiopsySample (File parentDir, File resultsDir) throws IOException {
		sampleName = parentDir.getName();
		this.resultsDir = resultsDir;
		
		//extract and check for the 4 required files (cfAlignment, cfIndex, normAlignment, normIndex)
		File[] allFiles = IO.extractFiles(parentDir);
		if (allFiles.length !=4) throw new IOException("ERROR: failed to find just four files representing the cfDNA and germline alignment files and their indexes for sample -> "+parentDir);
		boolean normIndexFound = false;
		boolean cfIndexFound = false;
		for (File f: allFiles) {
			String name = f.getName().toLowerCase();
			//is it the germline sample?
			if (name.contains("germ")) {
				if (name.endsWith(".bam") || name.endsWith(".cram")) {
					if (germlineAlignment != null) throw new IOException("ERROR: found two 'germ' alignment files, should only be one, for sample -> "+parentDir);
					else germlineAlignment = f;
				}
				else {
					if (name.endsWith(".crai") || name.endsWith(".bai")) normIndexFound = true;
					else throw new IOException("ERROR: found an unknown file type in sample -> "+parentDir);
				}
			}
			//must be the cf sample
			else {
				if (name.endsWith(".bam") || name.endsWith(".cram")) {
					if (cfAlignment != null) throw new IOException("ERROR: found two 'cf' alignment files, should only be one, for sample -> "+parentDir);
					else cfAlignment = f;
				}
				else {
					if (name.endsWith(".crai") || name.endsWith(".bai")) cfIndexFound = true;
					else throw new IOException("ERROR: found an unknown file type in sample -> "+parentDir);
				}
			}
		}
		if (normIndexFound == false || cfIndexFound == false || cfAlignment == null || germlineAlignment == null) {
			throw new IOException("ERROR: failed to find all 4 required files (cfAlignment, cfIndex, normAlignment, normIndex) in sample -> "+parentDir);
		}
	}
	
	public LookupJob[] buildCFJobs (LinkedHashMap<String, ArrayList<Bed>> blockNameRegions, SamReaderFactory samFactory, int minMQ) throws IOException {
		//make a reader
		cfReader = samFactory.open(cfAlignment);
		String name = Misc.removeExtension(cfAlignment.getParentFile().getName()+"_"+cfAlignment.getName());
		File bed = new File(resultsDir,  name+".bed.gz");
		cfBedOut = new Gzipper(bed);
		LookupJob[] jobs = new LookupJob[blockNameRegions.size()];
		int index = 0;
		for (String blockName: blockNameRegions.keySet()) {
			jobs[index++] = new LookupJob(blockNameRegions.get(blockName), cfReader, minMQ, blockName+"_"+name, cfBedOut);
		}
		cfJobs = jobs;
		return cfJobs;
	}
	
	public LookupJob[] buildGermlineJobs (LinkedHashMap<String, ArrayList<Bed>> blockNameRegions, SamReaderFactory samFactory, int minMQ) throws IOException {
		//make a reader
		germlineReader = samFactory.open(germlineAlignment);
		String name = Misc.removeExtension(germlineAlignment.getParentFile().getName()+"_"+germlineAlignment.getName());
		File bed = new File(resultsDir,  name+".bed.gz");
		germlineBedOut = new Gzipper(bed);
		LookupJob[] jobs = new LookupJob[blockNameRegions.size()];
		int index = 0;
		for (String blockName: blockNameRegions.keySet()) {
			jobs[index++] = new LookupJob(blockNameRegions.get(blockName), germlineReader, minMQ, blockName+"_"+name, germlineBedOut);
		}
		germlineJobs = jobs;
		return germlineJobs;
	}
	
	public float[] getCfCounts(int totalNumberRegions) {
		float[] allCounts = new float[totalNumberRegions];
		int counter = 0;
		for (LookupJob lj: cfJobs) {
			int[] counts = lj.getCounts();
			for (int i=0; i< counts.length; i++) allCounts[counter++] = counts[i];
		}
		return allCounts;
	}
	
	public float[] getGermlineCounts(int totalNumberRegions) {
		float[] allCounts = new float[totalNumberRegions];
		int counter = 0;
		for (LookupJob lj: germlineJobs) {
			int[] counts = lj.getCounts();
			for (int i=0; i< counts.length; i++) allCounts[counter++] = counts[i];
		}
		return allCounts;
	}

	public String getSampleName() {
		return sampleName;
	}

	public File getCfAlignment() {
		return cfAlignment;
	}

	public File getGermlineAlignment() {
		return germlineAlignment;
	}

	public SamReader getCfReader() {
		return cfReader;
	}

	public SamReader getGermlineReader() {
		return germlineReader;
	}

	public LookupJob[] getCfJobs() {
		return cfJobs;
	}

	public LookupJob[] getGermlineJobs() {
		return germlineJobs;
	}

	public Gzipper getCfBedOut() {
		return cfBedOut;
	}

	public Gzipper getGermlineBedOut() {
		return germlineBedOut;
	}

	public void addCfCounts(PrintWriter out) {
		//Gene_Region#_Type count
		//for each lookup job
		for (LookupJob lj: cfJobs) {
			String gene = lj.getRegions().get(0).getName();
			int[] counts = lj.getCounts();
			float[] nCounts = lj.getNormalizedCounts();
			for (int i=0; i<counts.length; i++)out.println(gene+"_"+i+"S\t"+counts[i]+"\t"+nCounts[i]);
		}
	}
	public void addGermlineCounts(PrintWriter out) {
		//Gene_Region#_Type count
		//for each lookup job
		for (LookupJob lj: germlineJobs) {
			String gene = lj.getRegions().get(0).getName();
			int[] counts = lj.getCounts();
			float[] nCounts = lj.getNormalizedCounts();
			for (int i=0; i<counts.length; i++)out.println(gene+"_"+i+"G\t"+counts[i]+"\t"+nCounts[i]);
		}
	}
	
	public void addNormalizedLog2Ratio(PrintWriter out) {
		//Gene_Region#_Type count
		//for each lookup job
		for (int i=0; i< germlineJobs.length; i++) {
			float[] sCounts = germlineJobs[i].getNormalizedCounts();
			float[] gCounts = cfJobs[i].getNormalizedCounts();
			String gene = cfJobs[i].getRegions().get(0).getName();
			//for each region in the gene group
			for (int j=0; j<sCounts.length; j++) {
				double rto = sCounts[j]/ gCounts[j];
				rto = Num.log2(rto);
				out.println(gene+"\t"+j+"\t"+rto);
			}
		}
	}
	
	
	
		
}
