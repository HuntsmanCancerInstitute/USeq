package edu.utah.seq.run.tempus;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.TimeUnit;
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
	
	private File testDir = null;
	private File clinReportDir = null;
	private boolean ready = true;
	
	
	
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
		//must have one json, and one dna fastq tar; rna is optional
		if (oneJson == false || oneDna == false) ready = false;
		else ready = true;
		
		//pull Json names
		ArrayList<String> jsonFileNames = new ArrayList<String>();
		for (TempusJsonParser tjp: jsonDatasets) jsonFileNames.add(tjp.getJsonFile().getName());
if (jsonDatasets.size()>1) System.exit(0);
		
		IO.pl("\tJson\t"+oneJson+"\t"+ jsonFileNames);
		IO.pl("\tDNA\t"+oneDna+"\t"+ dnaPaths);
		IO.pl("\tRNA\t"+oneRna+"\t"+ rnaPaths);
		IO.pl("\tOK\t"+ready);
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
		File fastqDir = new File (testDir, "Fastq");
		fastqDir.mkdir();
		
		//must be one DNA, might contain tumor and normal		
		cdw.cp(dnaPaths.get(0), new File (fastqDir, "/TarDNA/"+fetchName(dnaPaths.get(0))), true);
			
		//any RNA?
		if (rnaPaths.size() == 1) {
			cdw.cp(rnaPaths.get(0), new File (fastqDir, "/TarRNA/"+fetchName(rnaPaths.get(0))), true);
		}
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
		//move the report(s) into the ClinicalReport folder
		for (TempusJsonParser tjp : jsonDatasets) {
			File deId = tjp.getDeidentifiedJsonFile();
			deId.renameTo(new File(clinReportDir, deId.getName()));
		}
		
	}


	public String getTestID() {
		return testID;
	}

	public ArrayList<TempusJsonParser> getJsonDatasets() {
		return jsonDatasets;
	}

}
