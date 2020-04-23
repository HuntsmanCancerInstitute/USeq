package edu.utah.seq.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import util.gen.IO;
import util.gen.Misc;

public class Graphite {
	
	private File graphiteExecutable = null;
	private File fastaIndex = null;
	private File tempDir = null;
	private SamReaderFactory samFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	public Graphite(File graphiteExecutable, File fastaIndex, File tempDir) {
		this.graphiteExecutable = graphiteExecutable;
		this.fastaIndex = fastaIndex;
		this.tempDir = tempDir;
	}
	
	public HashMap<String, float[]> annotate (File vcfFile, File bamFile) throws IOException {
			String readGroupId = findRgId(bamFile);
			
			String[] cmd = {
					graphiteExecutable.getCanonicalPath(),
					"-v", vcfFile.getCanonicalPath(),
					"-b", bamFile.getCanonicalPath(),
					"-s", readGroupId,
					"-o", tempDir.getCanonicalPath(),
					"-f", fastaIndex.getCanonicalPath()
			};
			
			String[] issues = IO.executeViaProcessBuilder(cmd, true);
			if (issues.length !=0) {
				//ignore '[W::hts_idx_load2] The index file is older than...' issues, this is just a warning
				for (String s: issues) if (s.contains("older")==false) throw new IOException("Error executing:\n"+Misc.stringArrayToString(cmd, " ")+"\nIssues: "+Misc.stringArrayToString(issues, " "));
			}
			File gVcf = new File (tempDir, Misc.removeExtension(vcfFile.getName())+".vcf");			
			
			if (gVcf.exists() == false) throw new IOException("Error, failed to find the Graphite output vcf file "+gVcf+"\n"+Misc.stringArrayToString(cmd, " "));
			HashMap<String, float[]> keyAfDp = parseGraphiteVcf(gVcf);
			gVcf.delete();
			
			return keyAfDp;
	}

	private HashMap<String, float[]> parseGraphiteVcf(File gVcf) throws IOException {
		HashMap<String, float[]> keyStats = new HashMap<String, float[]>();
		BufferedReader in = IO.fetchBufferedReader(gVcf);
		String line;
		String[] fields;
		String[] lastSample;
		String[] counts;
		while ((line = in.readLine()) != null) {
			if (line.startsWith("#")) continue;
			fields = Misc.TAB.split(line);
			//the first value is the DP_NFP, the second DP4_NFP (fRef,rRef,fAlt,rAlt)
			//:5946:3131,2738,43,34:0:0,0,0,0:0:0,0,0,0:0:0,0,0,0:0:0,0,0,0:36:10,26,0,0:.
			lastSample = Misc.COLON.split(fields[fields.length-1].substring(1));
			counts = Misc.COMMA.split(lastSample[1]);
			//check that there are four
			if (counts.length !=4) throw new IOException("Error parsing DP4_NFP from "+line+" in "+gVcf);
			float total = Float.parseFloat(lastSample[0]);
			float alt = Float.parseFloat(counts[2]) + Float.parseFloat(counts[3]);
			float af = alt/total;
			//#CHROM POS REF ALT 0 1 3 4
			String key = fields[0]+"_"+ fields[1]+"_"+fields[3]+"_"+fields[4];
			keyStats.put(key, new float[] {af, total});
		}
		in.close();
		return keyStats;
	}

	private String findRgId(File bamFile) throws IOException {
		SamReader samReader = samFactory.open(bamFile);
		SAMFileHeader header = samReader.getFileHeader();
		List<SAMReadGroupRecord> readGroups = header.getReadGroups();
		if (readGroups.size() !=1) {
			samReader.close();
			throw new IOException("Error: more than one read group in "+bamFile);
		}
		String rg = readGroups.get(0).getId();
		samReader.close();
		return rg;
	}
	
	/*
	public static void main(String[] args) throws IOException {
		if (args.length != 5) IO.pl("Graphite fasta tempDir vcfFile bamFile");
		else {
			Graphite g = new Graphite(new File(args[0]), new File(args[1]), new File(args[2]));
			g.annotate(new File(args[3]), new File(args[4]));
		}
	}*/
}
