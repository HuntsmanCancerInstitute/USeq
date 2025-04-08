package edu.utah.seq.run.tempus;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class TempusPatient {
	
	//fields
	private String testID = null;
	private TempusDataWrangler cdw = null;
	private ArrayList<String[]> objectInfo = new ArrayList<String[]>();
	private ArrayList<TempusJsonParser> jsonDatasets = new ArrayList<TempusJsonParser>();
	private ArrayList<String> dnaPaths = new ArrayList<String>();
	private ArrayList<String> rnaPaths = new ArrayList<String>();
	private ArrayList<String> vcfPaths = new ArrayList<String>();
	private File testDir = null;
	private File clinReportDir = null;
	private boolean ready = true;
	//TL-20-A73CBE_T_DSQ1_2.fastq.gz
	private static final Pattern dnaSource = Pattern.compile(".+_([TN])_DSQ.+");
	
	
	
	public TempusPatient(String testID, TempusDataWrangler cdw) {
		this.testID = testID;
		this.cdw = cdw;
	}
	
	public void addObjectLine(String[] t) {
		objectInfo.add(t);
	}

	public void parseFileLines(int minHours) throws Exception {
		if (loadTarObjecArrayLists() == false) return;
		checkDatasets();
	}
	
	private void checkDatasets() throws IOException {
		boolean oneJson = checkJson();
		boolean oneDna = dnaPaths.size() == 1;
		boolean oneRna = rnaPaths.size() == 1;
		boolean someVcfs = vcfPaths.size() > 0;
		//must have one json, and one dna fastq tar; rna is optional
		if (oneJson == false || oneDna == false || someVcfs == false) ready = false;
		else ready = true;
		
		//pull Json names
		ArrayList<String> jsonFileNames = new ArrayList<String>();
		for (TempusJsonParser tjp: jsonDatasets) jsonFileNames.add(tjp.getJsonFile().getName());
		
		//IO.pl("\tJson\t"+oneJson+"\t"+ jsonFileNames);
		//IO.pl("\tDNA\t"+oneDna+"\t"+ dnaPaths);
		//IO.pl("\tRNA\t"+oneRna+"\t"+ rnaPaths);
		//IO.pl("\tVcfs\t"+someVcfs+"\t"+ vcfPaths);
		//IO.pl("\tOK\t"+ready);
	}

	private boolean checkJson() {
		//need to have only one DNA test
		int numDNA = 0;
		for (TempusJsonParser tjp: jsonDatasets) if (tjp.isDNATest()) numDNA++;
		if (numDNA == 1) return true;
		return false;
	}

	private boolean loadTarObjecArrayLists() {
		
		for (String[] tokens : objectInfo) {
			//2018-11-16 13:47:32 6877541013 TL-18-29F99A/TL-18-29F99A/DNA/FastQ/TL-18-29F99A_TL-18-29F99A-DNA-fastq.tar.gz
			//    0          1         2                  3
			String objectName = tokens[3];
			if (objectName.contains("/DNA/")) dnaPaths.add(objectName);
			else if (objectName.contains("/RNA/"))  rnaPaths.add(objectName);
			else {
				cdw.getErrorMessages().add("Couldn't source the fastq type from "+tokens[3]);
				IO.pl("ERROR");
				return false;
			}
		}
		return true;
	}

	public boolean isReady() {
		return ready;
	}

	public void downloadDatasets() throws Exception {
		
		//download vcfs
		for (String vcf: vcfPaths) {
			String name = vcf.substring(vcf.lastIndexOf("/")+1);
			File v = new File(clinReportDir, name);
			cdw.cp(vcf, v, false);
		}
		
		//download and process fastqs
		File fastqDir = new File (testDir, "Fastq");
		fastqDir.mkdir();
		
		//must be one DNA, might contain tumor and normal	
		File tarDna = new File (fastqDir, "/TarDNA/"+fetchName(dnaPaths.get(0)));
		File unpacked = new File(tarDna.getCanonicalPath()+".unpacked");
		
		if (unpacked.exists() == false) {
			//download it
			cdw.cp(dnaPaths.get(0), tarDna, true);
			//untar, delete, replace with empty placeholder
			unPackIt(tarDna);
		}
		
		//move the fastq.gz files?
		moveDNAFastq(tarDna.getParentFile(), fastqDir);
		

		
		//any RNA?
		if (rnaPaths.size() == 1) {
			File tarRna =  new File (fastqDir, "/TarRNA/"+fetchName(rnaPaths.get(0)));
			unpacked = new File(tarRna.getCanonicalPath()+".unpacked");
			
			if (unpacked.exists() == false) {
				//download it
				cdw.cp(rnaPaths.get(0), tarRna, true);
				//untar, delete, replace with empty placeholder
				unPackIt(tarRna);
			}

			//any fastq.gz files?
			moveRNAFastq(tarRna.getParentFile(), fastqDir);
		}
	}
	
	private void moveRNAFastq(File rnaDir, File fastqDir) throws IOException {
		File[] gz = IO.extractFiles(rnaDir, "q.gz");
		if (gz.length !=0) {
			if (gz.length == 2) {
				File tumorFastq = new File(fastqDir, "TumorRNA");
				tumorFastq.mkdir();
				for (File f: gz) {
					File moved = new File(tumorFastq, f.getName());
					f.renameTo(moved);
				}
			}
			else throw new IOException("Failed to find two RNA fastq files in "+rnaDir);
		}
	}

	
	private void moveDNAFastq(File dnaDir, File fastqDir) throws IOException {
		File[] gz = IO.extractFiles(dnaDir, "q.gz");		
		if (gz.length !=0) {
			File tumorFastq = null;
			File normalFastq = null;
			int numT = 0;
			int numN = 0;
			for (File f: gz) {
				Matcher mat = dnaSource.matcher(f.getName());
				if (mat.matches()) {
					String source = mat.group(1);
					if (source.equals("T")) {
						if (tumorFastq == null) {
							tumorFastq = new File(fastqDir, "TumorDNA");
							tumorFastq.mkdir();
						}
						File moved = new File(tumorFastq, f.getName());
						f.renameTo(moved);
						numT++;
					}
					else if (source.equals("N")) {
						if (normalFastq == null) {
							normalFastq = new File(fastqDir, "NormalDNA");
							normalFastq.mkdir();
						}
						File moved = new File(normalFastq, f.getName());
						f.renameTo(moved);
						numN++;
					}
					else throw new IOException("Failed to extract the correct DNA source (T or N) from "+f);
				}
				else throw new IOException("Failed to extract the DNA source (T or N) from "+f);
			}
			//check the numbers
			if (numT!=0 && numT!=2) throw new IOException("Failed to find two TumorDNA fastq files in "+fastqDir);
			if (numN!=0 && numN!=2) throw new IOException("Failed to find two NormalDNA fastq files in "+fastqDir);
		}
		
	}

	private void unPackIt(File tarGzFile) throws IOException {
		
		//untar it
		String tarPath = tarGzFile.getCanonicalPath();
		String[] cmd = {
				"tar", "-xf", tarPath, 
				"-C", tarGzFile.getCanonicalFile().getParent()+"/"
		};		
		int exitCode = IO.executeViaProcessBuilderReturnExit(cmd);		
		if (exitCode !=0) throw new IOException("ERROR: failed to untar "+tarGzFile);

		//delete it
		if (tarGzFile.delete() == false) throw new IOException("ERROR: failed to delete "+tarGzFile);
		
		//create placeholder
		File unpacked = new File(tarPath+".unpacked");
		unpacked.createNewFile();
	}

	public String fetchName (String tarPath) {
		String[] t = Misc.FORWARD_SLASH.split(tarPath);
		return (t[t.length-1]);
	}

	public void makeJobDirsMoveJson(String coreId) throws Exception {
		
		testDir = new File (cdw.getJobsDirectory(), coreId+"/Tempus/"+testID);
		testDir.mkdirs();
		clinReportDir = new File (testDir, "/ClinicalReport/");
		clinReportDir.mkdirs();
		clinReportDir = clinReportDir.getCanonicalFile();
		if (clinReportDir.exists() == false) throw new IOException("\nError: failed to make "+clinReportDir);
		// for some reason the files won't move!  Wonder if it's a latency issue with creating the dir, nope seems to be an issue with rw disk mounts
		//copy then set to delete, the report(s) into the ClinicalReport folder	
		for (TempusJsonParser tjp : jsonDatasets) {
			File deId = tjp.getDeidentifiedJsonFile();
			File inClinRep = new File(clinReportDir, deId.getName());
			if (IO.copyViaFileChannel(deId, inClinRep) == false) throw new IOException("\nError: failed to move "+deId+" to "+inClinRep);
			else deId.deleteOnExit();
		}
		
	}


	public String getTestID() {
		return testID;
	}

	public ArrayList<TempusJsonParser> getJsonDatasets() {
		return jsonDatasets;
	}

	public void addVcfLine(String[] tokens) {
		// TODO Auto-generated method stub
		
	}

	public ArrayList<String> getVcfPaths() {
		return vcfPaths;
	}

}
