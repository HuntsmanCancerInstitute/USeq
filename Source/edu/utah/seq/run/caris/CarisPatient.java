package edu.utah.seq.run.caris;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.TimeUnit;
import util.gen.IO;
import util.gen.Misc;

public class CarisPatient {
	
	//fields
	private String testID = null;
	private CarisDataWrangler cdw = null;
	private boolean tooYoung = false;
	private ArrayList<String[]> objectInfo = new ArrayList<String[]>();
	private ArrayList<String> vcfNames = new ArrayList<String>();
	private ArrayList<String> xmlNames = new ArrayList<String>();
	private ArrayList<String> pdfNames = new ArrayList<String>();
	private ArrayList<String> dnaNames = new ArrayList<String>();
	private ArrayList<String> rnaNames = new ArrayList<String>();
	
	private CarisXmlParser carisXml = null;
	
	private File testDir = null;
	private File clinReportDir = null;
	private boolean ready = true;
	
	public static final SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	public CarisPatient(String testID, CarisDataWrangler cdw) {
		this.testID = testID;
		this.cdw = cdw;
	}
	
	public void addObjectLine(String[] tokens) {
		objectInfo.add(tokens);
	}

	public void parseFileLines(int minHours) throws Exception {
		if (checkAges(minHours) == false) return;
		if (loadArrayLists() == false) return;
		checkDatasets();
	}
	
	private void checkDatasets() throws IOException {
		boolean oneVcf = vcfNames.size() == 1;
		boolean oneXml = xmlNames.size() == 1;
		boolean twoDna = dnaNames.size() == 2;
		boolean oneRna = rnaNames.size() == 1;
		//must have a vcf, xml, and two dna fastq files; if rna present, must have two files
		if (oneVcf == false || oneXml == false || twoDna == false || oneRna == true) ready = false;
		//if false see if it's just an RNA submission
		if (ready == false) {
			if (oneXml && rnaNames.size()==2 && oneVcf==false && twoDna==false) ready = true;
		}
		IO.pl("\tXml\t"+oneXml+"\t"+ xmlNames);
		IO.pl("\tVcf\t"+oneVcf+"\t"+ vcfNames);
		IO.pl("\tDNA\t"+twoDna+"\t"+ dnaNames);
		IO.pl("\tRNA\t"+(oneRna==false) +"\t"+ rnaNames);
		IO.pl("\tOK\t"+ready);
	}

	private boolean loadArrayLists() {
		//2022-03-12 11:33:49 2077226966 RNA_TN22-109832_S29_L001_R2_001.fastq.gz
		//    0         1          2                         3
		for (String[] tokens : objectInfo) {
			String fileName = tokens[3];
			if (fileName.endsWith(".fastq.gz")) {
				if (fileName.startsWith("RNA_")) rnaNames.add(fileName);
				else if (fileName.startsWith("DNA_")) dnaNames.add(fileName);
				else {
					cdw.getErrorMessages().add("Couldn't source the fastq type from "+tokens[3]);
					ready = false;
					return false;
				}
			}
			else if (fileName.endsWith(".xml")) xmlNames.add(fileName);
			else if (fileName.endsWith(".vcf")) vcfNames.add(fileName);
			else if (fileName.endsWith(".pdf")) pdfNames.add(fileName);
			//stuff to ignore
			else if (fileName.endsWith(".log")) {}
			else {
				cdw.getErrorMessages().add("Unknown file type type "+tokens[3]);
				ready = false;
				return false;
			}
		}
		return true;
	}

	public boolean checkAges(int minHours) throws ParseException {
		//check ages, must all be > min hours
		//2022-03-12 11:33:49 2077226966 RNA_TN22-109832_S29_L001_R2_001.fastq.gz
		//    0         1          2                         3
		for (String[] tokens : objectInfo) {
			Date objectDate = sdf.parse(tokens[0]+" "+tokens[1]);
			long diff = System.currentTimeMillis() - objectDate.getTime();
			long diffHours = TimeUnit.HOURS.convert(diff, TimeUnit.MILLISECONDS);
			if (diffHours< minHours) {
				IO.pl("\tTooYoung\t"+tokens[3]);
				ready = false;
				tooYoung = true;
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
		//cp in vcf
		if (vcfNames.size()==1) cdw.cp(vcfNames.get(0), new File (clinReportDir, vcfNames.get(0)));
		//cp in xml, already done so skip
		//cdw.cp(xmlNames.get(0), new File (clinDir, xmlNames.get(0)));
		//cp in DNA
		if (dnaNames.size()==2) {
			cdw.cp(dnaNames.get(0), new File (fastqDir, "/TumorDNA/"+dnaNames.get(0)));
			cdw.cp(dnaNames.get(1), new File (fastqDir, "/TumorDNA/"+dnaNames.get(1)));
		}
		//any RNA?
		if (rnaNames.size() == 2) {
			//cp in RNA0
			cdw.cp(rnaNames.get(0), new File (fastqDir, "/TumorRNA/"+rnaNames.get(0)));
			//cp in RNA1
			cdw.cp(rnaNames.get(1), new File (fastqDir, "/TumorRNA/"+rnaNames.get(1)));
		}
	}
	
	public void fetchXmlAndLoad() throws Exception {
		//cp from S3 the xml to the PHI directory
		File xml = new File (cdw.getPhiDirectory(), xmlNames.get(0));
		cdw.cp(xmlNames.get(0), xml);
				
		//load patient info
		carisXml = new CarisXmlParser(xml);
	}
	public void makeJobDirsMoveXml(String coreId) throws Exception {
		//parse the date from the xml, TN21-109147_2021-02-24_11_18.xml and TN20-170109_2021-01-20_21.31.xml
		String xmlFileName = carisXml.getXmlFile().getName();
		String[] tokens = Misc.UNDERSCORE.split(xmlFileName);
		String date = tokens[1];
		if (date.startsWith("20")==false || date.contains("-")==false) throw new IOException("\nERROR: failed to parse the date from "+xmlFileName);
		
		testDir = new File (cdw.getJobsDirectory(), coreId+"/Caris/"+testID+"_"+date+"/");
		testDir.mkdirs();
		clinReportDir = new File (testDir, "/ClinicalReport/");
		clinReportDir.mkdirs();
		//move the report into the ClinicalReport folder
		File deIdXml = carisXml.getDeidentifiedXmlFile();
		deIdXml.renameTo(new File(clinReportDir, deIdXml.getName()));
	}

	public boolean isTooYoung() {
		return tooYoung;
	}

	public CarisXmlParser getCarisXml() {
		return carisXml;
	}

	public String getTestID() {
		return testID;
	}

}
