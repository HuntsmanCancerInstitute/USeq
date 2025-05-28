package edu.utah.seq.run.tempus.v3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.pmr.PMRSearch;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonCollection;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonSummary;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Order;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Report;
import util.gen.IO;
import util.gen.Misc;

/**Container for all of the json reports related to a particular tempusOrderId and it's associated 
 * tumor sample. Typically a capture panel assay run on a tumor and it's matched normal, transcriptome 
 * on the tumor, and other tests: PDL-1, HRD, MMR, etc. Also included are the processed vcf files and fastq tar 
 * archives downloaded from Tempus. This whole thing is packaged up as a TJob for TNRunner to work on.
 * We're now seeing multiple tumor samples submitted for the same patient.  Primary and several cfDNA
 * samples.  Each TJob is run independently. */
public class TempusDataWranglerTumorV3 {
	
	private TempusV3JsonCollection collection = null;
	
	private ArrayList<File> jsonFiles = new ArrayList<File>();
	//size tab partialUri, e.g. 10524	TL-25-ACMHTPFA0I/DNA/TL-25-ACMHTPFA0I_20250423.soma.pindel.vcf
	private ArrayList<String> dnaAwsInfo = new ArrayList<String>();
	private ArrayList<String> rnaAwsInfo = new ArrayList<String>();
	private ArrayList<String> vcfAwsInfo = new ArrayList<String>();
	private boolean ready = true;
	
	// AAK6HT3yad/Tempus/25tnlyzo_20250527/
	private File testDir = null;
	// AAK6HT3yad/Tempus/25tnlyzo_20250527/ClinicalReport/
	private File clinReportDir = null;
	// AAK6HT3yad/Tempus/25tnlyzo_20250527/Fastq/  NormalDNA TumorDNA TumorRNA
	private File fastqDir = null;
	private static final Pattern dnaSource = Pattern.compile(".+_([TN])_DSQ.+");
	
	// aws --profile tempus s3 cp s3://bucket/
	public void addVcfTarDownloadCmds(ArrayList<String> cmdsToExecute, String awsCmdProfBucket, boolean verbose, ArrayList<File> tarFiles) throws IOException {
		//vcfs
		for (String vcfInfo: vcfAwsInfo) {
			addDownloadCmd(vcfInfo, clinReportDir, cmdsToExecute, awsCmdProfBucket, verbose, tarFiles);
		}
		
		//don't download a tar if the fastq already exist
		//dna tars, put all in TumorDNA, some will be NormalDNA after untarring
		if (dnaAwsInfo.size()!=0) {
			File tumorDNA = new File(fastqDir, "TumorDNA");
			tumorDNA.mkdirs(); //might already exist
			for (String dnaInfo: dnaAwsInfo) {
				File[] gzFq = IO.extractFiles(tumorDNA, "q.gz"); //any already present?
				if (gzFq.length != 2) addDownloadCmd(dnaInfo, tumorDNA, cmdsToExecute, awsCmdProfBucket, verbose, tarFiles);
			}
		}
		//rna tars
		if (rnaAwsInfo.size()!=0) {
			File tumorRNA = new File(fastqDir, "TumorRNA");
			tumorRNA.mkdirs();
			for (String rnaInfo: rnaAwsInfo) {
				File[] gzFq = IO.extractFiles(tumorRNA, "q.gz"); //any already present?
				if (gzFq.length != 2) addDownloadCmd(rnaInfo, tumorRNA, cmdsToExecute, awsCmdProfBucket, verbose, tarFiles);
			}
		}
	}
	
	private void addDownloadCmd(String awsInfo, File saveDir, ArrayList<String> cmdsToExecute, String awsCmdProfBucket, boolean verbose, ArrayList<File> tarFiles) throws IOException {
		String[] sizeAwsPath = Misc.WHITESPACE.split(awsInfo);
		long size = Long.parseLong(sizeAwsPath[0]);
		File toSave = new File (saveDir, sizeAwsPath[1].substring(sizeAwsPath[1].lastIndexOf('/')));
		//does the file exist with the correct size?
		if (toSave.exists() && toSave.length()==size) {
			if (verbose) IO.pl("File exists and same size, skipping: "+toSave);
		}
		else {
			// aws --profile tempus s3 cp s3://bucket/path/to/file.gz /full/path/to/save/file.gz
			String cmd = awsCmdProfBucket+sizeAwsPath[1]+" "+toSave.getCanonicalPath();
			cmdsToExecute.add(cmd);
			if (verbose) IO.pl("Adding '"+cmd+"'");
			if (toSave.getName().endsWith("tar.gz")) tarFiles.add(toSave);
		}
	}

	public void writeDeIdentifiedJsons() throws Exception {
		ArrayList<String> accessionIds = new ArrayList<String>();
		//accessionId_testCode_date_deid_firstNamePhy_lastNamePhy_gender.json
		//ClinicalReport/TL-24-Z84GNRRF_XT.V4_2024-01-23_deid_Neeraj_Agarwal_M.json		
		for (TempusV3JsonSummary sum: collection.getJsonSummaries()) {
			TempusV3Order order = sum.getTempusOrder();
			String name = order.getAccessionId()+"_"+
					order.getTestCode()+"_"+
					sum.getTempusReport().getSignOutDateNoTime()+"_deid_"+
					order.getFirstNamePhysician()+"_"+
					order.getLastNamePhysician()+ "_"+
					sum.getTempusPatient().getGender()+".json";
			File deIdFile = new File (clinReportDir, name);
			IO.writeString(sum.getMainObject().toString(5), deIdFile);
			accessionIds.add(order.getAccessionId());
		}
		//save manifest file
		//25tnlyzo_ClinicalVars_TL-25-6KOAZ9D9H1_TL-25-Q4TTL72L1T_TL-25-T3L8TUDC61_TL-25-UYR2M4OVJR_TL-25-VA94DV6NPJ.manifest.txt

	}
	
	public void writeManifest() throws Exception {
		ArrayList<String> accessionIds = new ArrayList<String>();
				
		for (TempusV3JsonSummary sum: collection.getJsonSummaries()) accessionIds.add(sum.getTempusOrder().getAccessionId());
		
		//collect all the files
		ArrayList<String> files = new ArrayList<String>();
		for (File f: IO.extractFiles(clinReportDir, ".json")) files.add(f.getCanonicalPath());
		for (File f: IO.extractFiles(clinReportDir, ".vcf")) files.add(f.getCanonicalPath());  //these will eventually be gzipped
		for (File f: IO.fetchAllFilesRecursively(fastqDir, "q.gz")) files.add(f.getCanonicalPath());
		
		//25tnlyzo_20250527_TL-25-6KOAZ9D9H1_TL-25-Q4TTL72L1T_TL-25-T3L8TUDC61_TL-25-UYR2M4OVJR_TL-25-VA94DV6NPJ.manifest.txt
		String[] aids = Misc.stringArrayListToStringArray(accessionIds);
		Arrays.sort(aids);
		String name = testDir.getName()+"_"+Misc.stringArrayToString(aids, "_")+ ".manifest.txt";
		IO.writeString(Misc.stringArrayListToString(files, "\n"), new File(clinReportDir, name));
	}

	public TempusDataWranglerTumorV3(TempusDataWranglerV3 tdw, TempusV3JsonCollection collection) {
		this.collection = collection;
		String vcfs = tdw.getVcfKeys();
		String tars = tdw.getTarKeys();
		
		boolean rnaTestPresent = false;
		boolean dnaTestPresent = false;
		
		//pull json accessionIds and jsonPaths
		for (TempusV3JsonSummary sum: collection.getJsonSummaries()) {
			jsonFiles.add(sum.getTempusReport().getJsonFile());
			String accessionId = sum.getTempusOrder().getAccessionId();
			//search for vcf file objects and tar objects
			Pattern pat = PMRSearch.fetchSearchPattern(accessionId, false, true);
			TreeSet<String> vcfMatches = PMRSearch.searchKeysString(vcfs, pat);
			vcfAwsInfo.addAll(vcfMatches);
			//ditto for tars
			TreeSet<String> tarMatches = PMRSearch.searchKeysString(tars, pat); 
			for (String tar: tarMatches) {
				if (tar.contains("/DNA/")) dnaAwsInfo.add(tar);
				else if (tar.contains("/RNA/")) rnaAwsInfo.add(tar);
			}
			String testCode = sum.getTempusOrder().getTestCode();
			if (testCode.startsWith("RS")) rnaTestPresent = true;
			else if (testCode.startsWith("XF") || testCode.startsWith("XT") || testCode.startsWith("XE")) dnaTestPresent = true;
			//else IO.pl("Unrecog: "+testCode);
		}
		
		//check if tars/ vcfs are present given the test codes observed
		
		//Rna tar not needed for processing so don't check. Still want to download it though.
		//if (rnaTestPresent && rnaAwsInfo.size()==0) ready = false;

		if (dnaTestPresent) {
			if (dnaAwsInfo.size()==0 || vcfAwsInfo.size()==0) ready = false;
		}
	}
	
	/**Looks for 2 q.gz files.*/
	public void checkRnaFastqDir() throws IOException {
		if (rnaAwsInfo.size()==0) return; // no rna
		File rnaDir = new File(fastqDir, "TumorRNA");
		File[] gz = IO.extractFiles(rnaDir, "q.gz");
		if (gz.length != 2) throw new IOException("Failed to find two TumorRNA fastq files in "+rnaDir);
	}

	
	public void moveDNAFastq() throws IOException {
		if (dnaAwsInfo.size()==0) return;
		File[] gz = IO.extractFiles(fastqDir, "q.gz");	
		File tumorFastq = new File(fastqDir, "TumorDNA");
		tumorFastq.mkdir();
		
		File normalFastq = null;
		
		if (gz.length !=0) {
			int numT = 0;
			int numN = 0;
			for (File f: gz) {
				Matcher mat = dnaSource.matcher(f.getName());
				if (mat.matches()) {
					String source = mat.group(1);
					if (source.equals("T")) {
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
			//check the numbers, must always be 2 tumorDNA fastqs, sometimes 2 normalDNA fastqs
			if (numT!=2) throw new IOException("Failed to find two TumorDNA fastq files in "+tumorFastq);
			if (numN!=0 && numN!=2) throw new IOException("Failed to find two NormalDNA fastq files in "+normalFastq);
		}
	}
	
	public void makeJobDirectories(File jobsDirectory, String procDate) throws IOException {
		String pmrId = collection.getJsonSummaries().get(0).getTempusPatient().getPmrId();
		testDir = new File(jobsDirectory, pmrId+"/Tempus/"+collection.getTempusOrderId()+"_"+procDate);
		testDir.mkdirs();
		clinReportDir = new File(testDir,"ClinicalReport");
		clinReportDir.mkdir();
		if (clinReportDir.exists() == false) throw new IOException("\nError: failed to make "+clinReportDir);
		fastqDir = new File(testDir,"Fastq");
		fastqDir.mkdir();
		if (fastqDir.exists() == false) throw new IOException("\nError: failed to make "+fastqDir);
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("TempusOrderId:\t"+ collection.getTempusOrderId());
		sb.append("\nJsonReports: (TestCode Accession JsonFile SignOutDate)\n");
		for (TempusV3JsonSummary sum: collection.getJsonSummaries()) {
			sb.append("\t");
			sb.append(sum.getTempusOrder().getTestCode());
			sb.append("\t");
			sb.append(sum.getTempusOrder().getAccessionId());
			sb.append("\t");
			sb.append(sum.getTempusReport().getJsonFile().getName());
			sb.append("\t");
			sb.append(sum.getTempusReport().getSignoutDate());
			sb.append("\n");
		}
		sb.append("DnaInfo:\t"+dnaAwsInfo);
		sb.append("\nRnaInfo:\t"+rnaAwsInfo);
		sb.append("\nVcfInfo:\t"+vcfAwsInfo);
		return sb.toString();
	}

	public boolean isReady() {
		return ready;
	}

	public TempusV3JsonCollection getCollection() {
		return collection;
	}

	public File getClinReportDir() {
		return clinReportDir;
	}

	public void setClinReportDir(File clinReportDir) {
		this.clinReportDir = clinReportDir;
	}

	public File getFastqDir() {
		return fastqDir;
	}

	public void setFastqDir(File fastqDir) {
		this.fastqDir = fastqDir;
	}

	public File getTestDir() {
		return testDir;
	}

	public void setTestDir(File testDir) {
		this.testDir = testDir;
	}







	
	
	/*

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


	 */

}
