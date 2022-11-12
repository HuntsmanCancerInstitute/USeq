package edu.utah.seq.cnv.cfdna;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class LiquidBiopsySample {
	
	//fields
	private String sampleName;
	private File cfAlignment;
	private File normalAlignment;
	private LookupJob[] cfJobs = null;
	private SamReader cfReader = null;
	private LookupJob[] normalJobs = null;
	private SamReader normalReader = null;
	private File resultsDir = null;
	private Gzipper cfBedOut = null;
	private Gzipper normalBedOut = null;
	
	//constructor
	public LiquidBiopsySample (File parentDir, File resultsDir) throws IOException {
		sampleName = parentDir.getName();
		this.resultsDir = resultsDir;
		
		//extract and check for the 4 required files (cfAlignment, cfIndex, normAlignment, normIndex)
		File[] allFiles = IO.extractFiles(parentDir);
		if (allFiles.length !=4) throw new IOException("ERROR: failed to find just four files representing the cfDNA and normal alignment files and their indexes for sample -> "+parentDir);
		boolean normIndexFound = false;
		boolean cfIndexFound = false;
		for (File f: allFiles) {
			String name = f.getName().toLowerCase();
			//is it the normal sample?
			if (name.contains("norm")) {
				if (name.endsWith(".bam") || name.endsWith(".cram")) {
					if (normalAlignment != null) throw new IOException("ERROR: found two 'norm' alignment files, should only be one, for sample -> "+parentDir);
					else normalAlignment = f;
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
		if (normIndexFound == false || cfIndexFound == false || cfAlignment == null || normalAlignment == null) {
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
	
	public LookupJob[] buildNormalJobs (LinkedHashMap<String, ArrayList<Bed>> blockNameRegions, SamReaderFactory samFactory, int minMQ) throws IOException {
		//make a reader
		normalReader = samFactory.open(normalAlignment);
		String name = Misc.removeExtension(normalAlignment.getParentFile().getName()+"_"+normalAlignment.getName());
		File bed = new File(resultsDir,  name+".bed.gz");
		normalBedOut = new Gzipper(bed);
		LookupJob[] jobs = new LookupJob[blockNameRegions.size()];
		int index = 0;
		for (String blockName: blockNameRegions.keySet()) {
			jobs[index++] = new LookupJob(blockNameRegions.get(blockName), normalReader, minMQ, blockName+"_"+name, normalBedOut);
		}
		normalJobs = jobs;
		return normalJobs;
	}

	public String getSampleName() {
		return sampleName;
	}

	public File getCfAlignment() {
		return cfAlignment;
	}

	public File getNormalAlignment() {
		return normalAlignment;
	}

	public SamReader getCfReader() {
		return cfReader;
	}

	public SamReader getNormalReader() {
		return normalReader;
	}

	public LookupJob[] getCfJobs() {
		return cfJobs;
	}

	public LookupJob[] getNormalJobs() {
		return normalJobs;
	}

	public Gzipper getCfBedOut() {
		return cfBedOut;
	}

	public Gzipper getNormalBedOut() {
		return normalBedOut;
	}

}
