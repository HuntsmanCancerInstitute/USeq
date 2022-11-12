package edu.utah.seq.cnv.cfdna;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Num;

public class LookupJob {

	//fields
	private ArrayList<Bed> regions = null;
	private SamReader reader = null;
	private boolean failed = false;
	private int minMappingQuality = 0;
	private String jobName = null;
	private int[] counts = null;
	private Gzipper out = null;
	private int numMidInt = 0;
	private int numMidNotInt = 0;

	public LookupJob (ArrayList<Bed> regionBlock, SamReader reader, int minMappingQuality, String jobName, Gzipper out) throws FileNotFoundException, IOException {
		this.reader = reader;
		this.regions = regionBlock;
		this.minMappingQuality = minMappingQuality;
		this.jobName = jobName;
		counts = new int[regions.size()];
		this.out = out;
	}


	public void doWork() throws IOException {	

		HashSet<String> countedReadNames = new HashSet<String>();
		HashMap<String, ArrayList<SAMRecord>> readNameAlignments = new HashMap<String, ArrayList<SAMRecord>>();

		//for each region
		for (int i=0; i< counts.length; i++) {
			Bed region = regions.get(i);
			readNameAlignments.clear();

			//fetch alignments, will throw exception if the chr isn't in the index
			SAMRecordIterator it = reader.queryOverlapping(region.getChromosome(), region.getStart()-1, region.getStop()+1);

			//Make pairs
			while (it.hasNext()){
				SAMRecord sam = it.next();
				if (sam.getMappingQuality() < minMappingQuality || sam.isSecondaryOrSupplementary() || sam.getDuplicateReadFlag() || sam.getReadFailsVendorQualityCheckFlag()) continue;
				String readName = sam.getReadName();
				ArrayList<SAMRecord> reads = readNameAlignments.get(readName);
				if (reads == null) {
					reads = new ArrayList<SAMRecord>();
					readNameAlignments.put(readName, reads);
				}
				reads.add(sam);
			}
			it.close();

			//for each pair
			for (String readName: readNameAlignments.keySet()) {
				//already been counted?
				if (countedReadNames.contains(readName)) continue;

				//find middle and see if it intersects this region
				boolean intersectsRegion = false;
				ArrayList<SAMRecord> sams = readNameAlignments.get(readName);
				int size = sams.size();
				int middle = 0;
				String name = null;
				if (size == 1) {
					SAMRecord sam = sams.get(0);
					middle = calculateMiddle(sam);
					intersectsRegion = intersects(region, middle);
					name="S_";
				}
				else if (size == 2) {
					middle = calculateMiddle(sams);
					intersectsRegion = intersects(region, middle);
					name="P_";
				}
				else throw new IOException("\nERROR: seeing more than two alignments passing thresholds for "+readName+" with job "+jobName);

				if (intersectsRegion) {
					counts[i]++;
					countedReadNames.add(readName);
					out.println(sams.get(0).getReferenceName()+"\t"+middle+"\t"+(middle+1)+"\t"+name+sams.get(0).getReadName()+"\t0\t.");
					numMidInt++;
				}
				else numMidNotInt++;
			}
		}
		countedReadNames = null;
		readNameAlignments = null;
		IO.pl(jobName+"\t"+numMidInt+"\t"+numMidNotInt);
		IO.pl("\t"+Num.intArrayToString(counts, ","));
	}

	private int calculateMiddle(ArrayList<SAMRecord> sams) throws IOException {
		int smallest = Integer.MAX_VALUE;
		int biggest = -1;
		for (SAMRecord sam : sams) {
			int start = sam.getAlignmentStart();
			if (start < smallest) smallest = start;
			if (start > biggest) biggest = start;
			int end = sam.getAlignmentEnd();
			if (end < smallest) smallest = end;
			if (end > biggest) biggest = end;
		}
		double start = smallest;
		double stop = biggest;
		double middle = (stop-start)/2.0 + start;
		return (int)Math.round(middle);
	}

	private boolean intersects(Bed region, int middle) {
		if (middle < region.getStart() || middle >= region.getStop()) return false;
		return true;
	}

	private int calculateMiddle(SAMRecord sam) throws IOException {
		double start = sam.getAlignmentStart();
		double stop = sam.getAlignmentEnd();		
		//start is always less than stop
		double middle = (stop-start)/2.0 + start;
		return (int)Math.round(middle);
	}

	public boolean isFailed() {
		return failed;
	}

	public int[] getCounts() {
		return counts;
	}

	public ArrayList<Bed> getRegions() {
		return regions;
	}



}
