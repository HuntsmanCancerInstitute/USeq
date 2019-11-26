package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class TNSample {

	private String id;
	private File rootDir;
	private TNRunner tnRunner;
	private boolean forceRestart;
	private boolean restartFailedJobs;
	private boolean softRestart;
	private static final Pattern slurmJobId = Pattern.compile("slurm-(\\d+).out");
	private ArrayList<String> info = new ArrayList<String>();
	private boolean deleteBamConcordance = false;

	//Fastq
	private FastqDataset tumorDNAFastq = null;
	private FastqDataset normalDNAFastq = null;
	private FastqDataset tumorTransFastq = null;
	private boolean fastqIssue = false;

	//Alignment and passing region bed files
	private File[] tumorDNABamBedGvcf = null;
	private File[] normalDNABamBedGvcf = null;
	private File tumorRNABam = null;

	//Somatic variant Vcf calls
	private File somaticVariants = null;
	private File mergedSomaticVariants = null;

	//Joint genotyping
	private File[] germlineVcf = null;

	//Job status
	private boolean failed = false;
	private boolean running = false;
	private int numJobsLaunched = 0;


	public TNSample(File rootDir, TNRunner tnRunner) throws IOException{
		this.rootDir = rootDir;
		this.tnRunner = tnRunner;
		this.forceRestart = tnRunner.isForceRestart();
		this.restartFailedJobs = tnRunner.isRestartFailed();
		this.softRestart = tnRunner.isSoftRestart();

		//assign id
		id = rootDir.getName();
		info.add("\n"+id);

		//look for fastq folders and check they contain 2 xxx.gz files, some may be null
		boolean foundFastq = checkFastq();

		//launch available alignments
		if (foundFastq) checkAlignments();
		else {
			fastqIssue = true;
			return;
		}

		//launch t/n somatic DNA analysis?
		if (tumorDNABamBedGvcf != null && normalDNABamBedGvcf != null && tnRunner.getSomaticVarCallDocs() != null) somDNACall();
		
		//merge somatic vars with clinical test results?
		if (tnRunner.getClinicalVcfDocs()!= null && somaticVariants != null) parseMergeClinicalVars();

		//annotate somatic vcf?
		if (somaticVariants != null && tnRunner.getVarAnnoDocs() != null) annotateSomaticVcf();

		//bam concordance?
		if (tnRunner.getBamConcordanceDocs() != null) bamConcordance();

		//Annotate normal vcf
		if (tnRunner.getVarAnnoDocs() != null) annotateGermlineVcf();

		//copy ratio/ number
		if (tnRunner.getCopyRatioDocs() != null) copyRatioAnalysis();
		
		//msi
		if (tnRunner.getMsiDocs() != null) msiAnalysis();

		//print messages?
		if (failed || tnRunner.isVerbose()) {
			Misc.printArray(info);
			IO.pl("Failures found? "+failed);
			IO.pl("Running jobs? "+running);
		}
	}

	
		private void parseMergeClinicalVars() throws IOException {
			info.add("Checking clinical variant integration...");
			
			//look for json file
			File[] jsonTestResults = IO.extractFiles(new File(rootDir, "ClinicalReport"), ".json");
			if (jsonTestResults == null || jsonTestResults.length !=1) {
				info.add("\t\tJsonReport\tFAILED to find one xxx.json clinical test report file");
				failed = true;
				return;
			}
			
			//make dir, ok if it already exists
			File jobDir = new File (rootDir, "SomaticVariantCalls/"+id+"_ClinicalVars");
			jobDir.mkdirs();
			
			File[] toLink = new File[]{somaticVariants, new File(somaticVariants+".tbi"), jsonTestResults[0]};

			//any files?
			HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
			if (nameFile.size() == 0) launch(jobDir, toLink, tnRunner.getClinicalVcfDocs());

			//OK some files are present
			//COMPLETE
			else if (nameFile.containsKey("COMPLETE")){
				//find the final vcf file
				File[] vcf = IO.extractFiles(new File(jobDir, "Vcfs"), "_final.vcf.gz");
				if (vcf == null || vcf.length !=1) {
					clearAndFail(jobDir, "\tThe clinical variant parsing and merging workflow was marked COMPLETE but failed to find the xxx_final.vcf.gz file in the Vcfs/ in "+jobDir);
					return;
				}
				mergedSomaticVariants = vcf[0];
				//remove the linked files
				for (File f: toLink) new File(jobDir, f.getName()).delete();
				info.add("\t\tCOMPLETE "+jobDir);
			}
			
			else checkJob(nameFile, jobDir, toLink, tnRunner.getClinicalVcfDocs());
		}

	public void annotateGermlineVcf() throws IOException {
		//look for genotyped vcf
		File vcf = new File (rootDir, "GermlineVariantCalling/"+id+"_NormalDNA/"+id+"_NormalDNA_Hg38_JointGenotyped.vcf.gz");
		if (vcf.exists() == false) return;
		germlineVcf = new File[]{vcf, new File(rootDir, "GermlineVariantCalling/"+id+"_NormalDNA/"+id+"_NormalDNA_Hg38_JointGenotyped.vcf.gz.tbi") };

		info.add("Checking germline variant annotation...");		

		//make dir, ok if it already exists
		File jobDir = new File (rootDir.getCanonicalFile(), "GermlineVariantCalling/"+id+"_NormalDNAAnno");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			IO.writeString(tnRunner.getGermlineAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
			launch(jobDir, germlineVcf, tnRunner.getVarAnnoDocs());
		}

		//OK some files are present
		else if (nameFile.containsKey("COMPLETE")){
			//remove the linked vcfs and the config file
			new File(jobDir, vcf.getName()).delete();
			new File(jobDir, vcf.getName()+".tbi").delete();
			new File(jobDir, "annotatedVcfParser.config.txt").delete();
			info.add("\t\tCOMPLETE "+jobDir);
		}
		else {
			if (checkJob(nameFile, jobDir, germlineVcf, tnRunner.getVarAnnoDocs())){
				IO.writeString(tnRunner.getGermlineAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
			}
		}

	}

	private void bamConcordance() throws IOException {
		info.add("Checking bam concordance...");
		File jobDir = new File (rootDir, "SampleConcordance/"+id+"_BamConcordance");

		//Started a new alignment? then delete the concordance if it exists
		if (deleteBamConcordance && jobDir.exists()){
			info.add("\t\tNew alignments started.");
			cancelDeleteJobDir(IO.fetchNamesAndFiles(jobDir), jobDir, info, true);
			return;
		}

		//check to see if the alignments are present, sometimes these datasets are never coming
		boolean goT = false;
		boolean goN = false;
		boolean goTT = false;
		int numExist = 0;
		if (tumorDNAFastq.isFastqDirExists()){
			if (tumorDNABamBedGvcf != null) goT = true;
			numExist++;
		}
		else goT = true;
		if (normalDNAFastq.isFastqDirExists()){
			if (normalDNABamBedGvcf != null) {
				goN = true;
				//skip the mock normal
				if (normalDNABamBedGvcf[0].getName().equals("NA12878_NormalDNA_Hg38_final.bam")) numExist--;
			}
			numExist++;
		}
		else goN = true;
		if (tumorTransFastq.isFastqDirExists()){
			if (tumorRNABam != null) goTT = true;
			numExist++;
		}
		else goTT = true;
		if (numExist < 2) return;
		if (goT == false || goN == false || goTT == false) return;

		//make dir, ok if it already exists
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			removeBamConcordanceLinks(jobDir);
			createBamConcordanceLinks(jobDir);
			launch(jobDir, null, tnRunner.getBamConcordanceDocs());
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			removeBamConcordanceLinks(jobDir);
			info.add("\t\tCOMPLETE "+jobDir);
		}
		
		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getBamConcordanceDocs())){
				removeBamConcordanceLinks(jobDir);
				createBamConcordanceLinks(jobDir);
			}
		}
		
	}

	private void createBamConcordanceLinks(File jobDir) throws IOException{
		File canJobDir = jobDir.getCanonicalFile();
		if (tumorDNABamBedGvcf != null) {
			File index = fetchBamIndex(tumorDNABamBedGvcf[0]);
			Files.createSymbolicLink(new File(canJobDir, "tumorDNA.bam").toPath(), tumorDNABamBedGvcf[0].toPath());
			Files.createSymbolicLink(new File(canJobDir, "tumorDNA.bai").toPath(), index.toPath());
		}
		//don't include NA12878_NormalDNA_Hg38_final.bam, this is used as a mock normal 
		if (normalDNABamBedGvcf != null && normalDNABamBedGvcf[0].getName().equals("NA12878_NormalDNA_Hg38_final.bam")==false) {
			File index = fetchBamIndex(normalDNABamBedGvcf[0]);
			Files.createSymbolicLink(new File(canJobDir, "normalDNA.bam").toPath(), normalDNABamBedGvcf[0].toPath());
			Files.createSymbolicLink(new File(canJobDir, "normalDNA.bai").toPath(), index.toPath());
		}
		if (tumorRNABam != null) {
			File index = fetchBamIndex(tumorRNABam);
			Files.createSymbolicLink(new File(canJobDir, "tumorRNA.bam").toPath(), tumorRNABam.toPath());
			Files.createSymbolicLink(new File(canJobDir, "tumorRNA.bai").toPath(), index.toPath());
		}
		//pull in Foundation bams?
		ArrayList<File> toLink = fetchOtherBams();
		if (toLink != null && toLink.size() != 0){
			for (File f: toLink) Files.createSymbolicLink(new File(canJobDir, "F_"+f.getName()).toPath(), f.toPath());
		}
	}

	private ArrayList<File> fetchOtherBams() throws IOException {
		HashMap<String, File> otherPatientDatasets = tnRunner.getOtherBams();
		if (otherPatientDatasets != null){
			File matchedPatientDir = otherPatientDatasets.get(id);
			if (matchedPatientDir == null) return null;
			ArrayList<File> otherBams = IO.fetchAllFilesRecursively(matchedPatientDir, ".bam");

			//watch out for dups due to linking and failed jobs
			HashMap<String, File> uniqueFiles = new HashMap<String, File>();
			for (File ob: otherBams) uniqueFiles.put(ob.getCanonicalPath(), ob.getCanonicalFile());

			ArrayList<File> toLink = new ArrayList<File>();
			for (File ob: uniqueFiles.values()){
				//DNA bam and RNA bam
				if (ob.getName().endsWith("_final.bam") || ob.getName().endsWith("_Hg38.bam")) {
					File index = fetchBamIndex(ob);
					if (index != null){
						toLink.add(ob);
						toLink.add(index);
					}
				}
			}
			return toLink;
		}
		return null;
	}

	private void removeBamConcordanceLinks(File jobDir) throws IOException{
		IO.deleteFiles(jobDir, ".bam");
		IO.deleteFiles(jobDir, ".bai");
	}

	private void annotateSomaticVcf() throws IOException {
		info.add("Checking somatic variant annotation...");	
		
		//waiting on clinical vars?
		if (tnRunner.getClinicalVcfDocs() != null && mergedSomaticVariants == null) {
			running = true;
			return;
		}
		
		//make dir, ok if it already exists
		File jobDir = new File (rootDir, "SomaticVariantCalls/"+id+"_Anno");
		jobDir.mkdirs();
		
		File[] toLink = null;
		if (mergedSomaticVariants != null)toLink = new File[]{mergedSomaticVariants, new File(mergedSomaticVariants+".tbi")};
		else toLink = new File[]{somaticVariants, new File(somaticVariants+".tbi")};
		

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			IO.writeString(tnRunner.getSomaticAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
			launch(jobDir, toLink, tnRunner.getVarAnnoDocs());
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//remove the linked vcfs and the config file
			for (File f: toLink) new File(jobDir, f.getName()).delete();
			new File(jobDir, "annotatedVcfParser.config.txt").delete();
			info.add("\t\tCOMPLETE "+jobDir);
		}
		
		else {
			if (checkJob(nameFile, jobDir, toLink, tnRunner.getVarAnnoDocs())){
				IO.writeString(tnRunner.getSomaticAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
			}
		}
	}

	private void somDNACall() throws IOException {
		info.add("Checking somatic variant calling...");		

		//make dir, ok if it already exists
		File jobDir = new File (rootDir,"SomaticVariantCalls/"+id+"_Illumina");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			createSomaticVariantLinks(jobDir);
			launch(jobDir, null, tnRunner.getSomaticVarCallDocs());
		}
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final vcf file
			File[] vcf = IO.extractFiles(new File(jobDir, "Vcfs"), "_final.vcf.gz");
			if (vcf == null || vcf.length !=1) {
				clearAndFail(jobDir, "\tThe somatic variant calling was marked COMPLETE but failed to find the final vcf file in the Vcfs/ in "+jobDir);
				return;
			}
			somaticVariants = vcf[0];
			//remove the linked fastq
			removeSomaticLinks(jobDir);
			info.add("\t\tCOMPLETE "+jobDir);
		}
		else {
			if (checkJob(nameFile,jobDir,null, tnRunner.getSomaticVarCallDocs())){
				removeSomaticLinks(jobDir);
				createSomaticVariantLinks(jobDir);
			}
		}

	}
	
	private void msiAnalysis() throws IOException {
		
		info.add("Checking MSI status analysis...");
		
		//look for tumor and normal bam files
		if (tumorDNABamBedGvcf == null || normalDNABamBedGvcf == null) {
			info.add("\t\tMissing one or both T/N alignment files.");
			return;
		}

		//make dir, ok if it already exists
		File jobDir = new File (rootDir,"SomaticVariantCalls/"+id+"_Msi");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			createMsiLinks(jobDir);
			launch(jobDir, null, tnRunner.getMsiDocs());
		}
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final txt results file
			File[] res = IO.extractFiles(jobDir, "_Mantis.txt");
			if (res == null || res.length !=1) {
				clearAndFail(jobDir, "\tThe msi status calling was marked COMPLETE but failed to find the final xxx_Mantis.txt file in "+jobDir);
				return;
			}
			//remove the linked fastq
			removeMsiLinks(jobDir);
			info.add("\t\tCOMPLETE "+jobDir);
		}
		else {
			if (checkJob(nameFile,jobDir,null, tnRunner.getMsiDocs())){
				createMsiLinks(jobDir);
			}
		}
	}
	
	private void createMsiLinks(File jobDir) throws IOException {
		//remove any linked bam and bed files
		removeMsiLinks(jobDir);
		//soft link in the new ones
		Path tumorBam = tumorDNABamBedGvcf[0].toPath();
		Path tumorBai = fetchBamIndex(tumorDNABamBedGvcf[0]).toPath();
		Path normalBam = normalDNABamBedGvcf[0].toPath();
		Path normalBai = fetchBamIndex(normalDNABamBedGvcf[0]).toPath();
		//must have index as xxx.bam.bai, won't work as xxx.bai
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bam").toPath(), tumorBam);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bam.bai").toPath(), tumorBai);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bam").toPath(), normalBam);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bam.bai").toPath(), normalBai);
	}
	
	private void removeMsiLinks(File jobDir) throws IOException{
		File f = jobDir.getCanonicalFile();
		new File(f, "tumor.bam").delete();
		new File(f, "tumor.bam.bai").delete();
		new File(f, "normal.bam").delete();
		new File(f, "normal.bam.bai").delete();
	}

	private void copyRatioAnalysis() throws IOException {
		if (tnRunner.getCopyRatioDocs() == null) return;
		info.add("Checking copy ratio calling...");	
		
		//look for germline vcf
		if (germlineVcf == null) {
			File vcf = new File (rootDir, "GermlineVariantCalling/"+id+"_NormalDNA/"+id+"_NormalDNA_Hg38_JointGenotyped.vcf.gz");
			if (vcf.exists() == false) {
				info.add("\tNo germline vcf.");
				return;
			};
			germlineVcf = new File[]{vcf, new File(rootDir, "GermlineVariantCalling/"+id+"_NormalDNA/"+id+"_NormalDNA_Hg38_JointGenotyped.vcf.gz.tbi") };
		}

		//look for tumor normal alignments
		if (tumorDNABamBedGvcf == null || normalDNABamBedGvcf == null) return;

			

		//make dir, ok if it already exists
		File jobDir = new File (rootDir, "CopyAnalysis/"+id+"_GATKCopyRatio");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			createCopyRatioLinks(jobDir);
			launch(jobDir, null, tnRunner.getCopyRatioDocs());
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final vcf file
			File[] segs = IO.extractFiles(new File(jobDir, "Results"), ".called.seg.xls");
			if (segs == null || segs.length !=1) {
				clearAndFail(jobDir, "\tThe copy ratio calling was marked COMPLETE but failed to find the final seg file in the Results/ dir in "+jobDir);
				return;
			}
			//remove the linked fastq
			removeCopyRatioLinks(jobDir);
			info.add("\t\tCOMPLETE "+jobDir);
		}
		
		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getCopyRatioDocs())){
				createCopyRatioLinks(jobDir);
			}
		}
	}

	private void createCopyRatioLinks(File jobDir) throws IOException {
		//remove any linked files
		removeCopyRatioLinks(jobDir);

		//pull gender from the xxxInfo.json.gz or .json info files
		int gender = fetchGender(jobDir);
		
		if (gender == 0) throw new IOException("ERROR: failed to find or parse the gender/sex info from the .json file in "+rootDir);
		Path bkg = null;
		if (gender == 1) bkg = tnRunner.getMaleBkg().toPath();
		else if (gender == 2 ) bkg = tnRunner.getFemaleBkg().toPath();
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"bkgPoN.hdf5").toPath(), bkg);

		//soft link in the new ones
		Path tumorBam = tumorDNABamBedGvcf[0].toPath();
		Path tumorBai = fetchBamIndex(tumorDNABamBedGvcf[0]).toPath();
		Path normalBam = normalDNABamBedGvcf[0].toPath();
		Path normalBai = fetchBamIndex(normalDNABamBedGvcf[0]).toPath();

		Path normalVcf = germlineVcf[0].toPath();
		Path normalTbi = germlineVcf[1].toPath();
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bam").toPath(), tumorBam);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bai").toPath(), tumorBai);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bam").toPath(), normalBam);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bai").toPath(), normalBai);

		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.vcf.gz").toPath(), normalVcf);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.vcf.gz.tbi").toPath(), normalTbi);
	}

	//throw new IOException("ERROR: failed to find the xxx_Info.json.gz file in "+rootDir);
	private int fetchGender(File jobDir) throws IOException {
		int gender = 0;
	
		//is it from Avatar?
		File[] info = IO.extractFiles(rootDir, "Info.json.gz");
		if (info == null || info.length !=1) {
			//try to fetch from Tempus
			File d = new File(rootDir, "ClinicalReport");
			if (d.exists() == false) return 0;
			info = IO.extractFiles(d, ".json");
			if (info == null || info.length !=1) return 0;
		}
		String[] lines = IO.loadFile(info[0]);
		
		for (String s: lines){
			if (s.contains("Gender") || s.contains("\"sex\"")){
				if (s.contains("F")) gender = 2;
				else if (s.contains("M")) gender = 1;
			}
		}
		
		return gender;
	}


	private void removeCopyRatioLinks(File jobDir) throws IOException{
		File f = jobDir.getCanonicalFile();
		new File(f, "tumor.bam").delete();
		new File(f, "tumor.bai").delete();
		new File(f, "bkgPoN.hdf5").delete();
		new File(f, "normal.bam").delete();
		new File(f, "normal.bai").delete();
		new File(f, "normal.vcf.gz").delete();
		new File(f, "normal.vcf.gz.tbi").delete();
	}

	private void checkAlignments() throws IOException {
		info.add("Checking alignments...");
		//fastq available, align
		tumorDNABamBedGvcf = DNAAlignQC(tumorDNAFastq, tnRunner.getMinReadCoverageTumor());
		normalDNABamBedGvcf = DNAAlignQC(normalDNAFastq, tnRunner.getMinReadCoverageNormal());
		RNAAlignQC(tumorTransFastq);
	}

	/**For RNASeq Alignments and Picard CollectRNASeqMetrics*/
	private void RNAAlignQC(FastqDataset fd) throws IOException {
		if (fd.getFastqs() == null || tnRunner.getRNAAlignQCDocs() == null) return;
		info.add("\t"+fd.getName()+":");

		//make dir, ok if it already exists
		File jobDir = new File (rootDir, "Alignment/"+id+"_"+fd.getName()).getCanonicalFile();
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			createDNAAlignLinks(fd.getFastqs(), jobDir);
			launch(jobDir, null, tnRunner.getRNAAlignQCDocs());
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final bam and bed
			File bamDir = new File(jobDir, "Bam");
			File[] finalBam = IO.extractFiles(bamDir, ".bam");

			if (finalBam == null || finalBam.length !=1) {
				clearAndFail(jobDir, "\tThe alignemnt was marked COMPLETE but failed to find the bam file in "+bamDir);
			}
			else {
				//remove the linked fastq
				new File(jobDir, "1.fastq.gz").delete();
				new File(jobDir, "2.fastq.gz").delete();
				info.add("\t\tCOMPLETE "+jobDir);
				tumorRNABam = finalBam[0];
			}
		}
		
		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getRNAAlignQCDocs())){
				createDNAAlignLinks(fd.getFastqs(), jobDir);
			}
		}
	}
	
	private void clearAndFail(File jobDir, String infoLine) throws IOException{
		info.add(infoLine);
		failed = true;
		removeProgressFiles(jobDir);
		new File(jobDir, "FAILED").createNewFile();
	}


	/**For DNAs, diff read coverage for passing bed generation
	 * @throws IOException */
	private File[] DNAAlignQC(FastqDataset fd, int minReadCoverage) throws IOException {
		if (fd.getFastqs() == null || tnRunner.getDNAAlignQCDocs() == null) return null;
		info.add("\t"+fd.getName()+":");

		//make dir, ok if it already exists
		File jobDir = new File (rootDir, "Alignment/"+id+"_"+fd.getName());
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			createDNAAlignLinks(fd.getFastqs(), jobDir);
			deleteBamConcordance = true;
			IO.writeString("-c "+ minReadCoverage, new File(jobDir, "sam2USeq.config.txt"));
			launch(jobDir, null, tnRunner.getDNAAlignQCDocs());
		}

		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final bam and bed
			File[] finalBam = IO.extractFiles(new File(jobDir, "Bam"), "_final.bam");
			File[] passingBed = IO.extractFiles(new File(jobDir, "QC"), "_Pass.bed.gz");
			File[] gvcf = IO.extractFiles(new File(jobDir, "Vcfs"), "_Haplo.g.vcf.gz");
			File[] gvcfIndex = IO.extractFiles(new File(jobDir, "Vcfs"), "_Haplo.g.vcf.gz.tbi");

			//check for problems
			String problem = null;
			if (finalBam == null || finalBam.length !=1) {
				problem = "\tThe alignemnt was marked COMPLETE but failed to find the final bam file in the Bam/ in "+jobDir;
			}
			else if (passingBed == null || passingBed.length !=1) {
				problem = "\tThe alignemnt was marked COMPLETE but failed to find the passing bed file in QC/ in "+jobDir;
			}
			else if (gvcf == null || gvcf.length !=1) {
				problem = "\tThe alignemnt was marked COMPLETE but failed to find the gvcf file in Vcfs/ inside "+jobDir;
			}
			if (problem != null){
				clearAndFail(jobDir, problem);
				return null;
			}

			//remove the linked fastq and the sam2USeq.config.txt
			new File(jobDir, "1.fastq.gz").delete();
			new File(jobDir, "2.fastq.gz").delete();
			new File(jobDir, "sam2USeq.config.txt").delete();
			info.add("\t\tCOMPLETE "+jobDir);
			return new File[]{finalBam[0], passingBed[0], gvcf[0], gvcfIndex[0]};
		}
		
		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getDNAAlignQCDocs())){
				createDNAAlignLinks(fd.getFastqs(), jobDir);
				deleteBamConcordance = true;
				IO.writeString("-c "+ minReadCoverage, new File(jobDir, "sam2USeq.config.txt"));
			}
		}
		return null;
	}

	private void createDNAAlignLinks(File[] fastq, File alignDir) throws IOException {
		//remove any linked fastq files
		File f = alignDir.getCanonicalFile();
		new File(f, "1.fastq.gz").delete();
		new File(f,"2.fastq.gz").delete();

		//soft link in the new ones
		Path real1 = fastq[0].toPath();
		Path real2 = fastq[1].toPath();
		Files.createSymbolicLink(new File(alignDir.getCanonicalFile(),"1.fastq.gz").toPath(), real1);
		Files.createSymbolicLink(new File(alignDir.getCanonicalFile(),"2.fastq.gz").toPath(), real2);
	}
	


	private void createSomaticVariantLinks(File jobDir) throws IOException {
		//remove any linked bam and bed files
		removeSomaticLinks(jobDir);

		//soft link in the new ones
		Path tumorBam = tumorDNABamBedGvcf[0].toPath();
		Path tumorBai = fetchBamIndex(tumorDNABamBedGvcf[0]).toPath();
		Path tumorBed = tumorDNABamBedGvcf[1].toPath();
		Path normalBam = normalDNABamBedGvcf[0].toPath();
		Path normalBai = fetchBamIndex(normalDNABamBedGvcf[0]).toPath();
		Path normalBed = normalDNABamBedGvcf[1].toPath();
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bam").toPath(), tumorBam);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bai").toPath(), tumorBai);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bed.gz").toPath(), tumorBed);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bam").toPath(), normalBam);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bai").toPath(), normalBai);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bed.gz").toPath(), normalBed);
	}

	private void removeSomaticLinks(File jobDir) throws IOException{
		File f = jobDir.getCanonicalFile();
		new File(f, "tumor.bam").delete();
		new File(f, "tumor.bai").delete();
		new File(f, "tumor.bed.gz").delete();
		new File(f, "normal.bam").delete();
		new File(f, "normal.bai").delete();
		new File(f, "normal.bed.gz").delete();
	}

	public static boolean checkQueue(HashMap<String, File> nameFile, File jobDir, ArrayList<String> info) throws IOException{
		//find slurm script(s), pull the ID, and check it is still in the queue
		ArrayList<File> jobs = new ArrayList<File>();
		for (String s: nameFile.keySet()) {
			if (s.startsWith("slurm")) {
				jobs.add(nameFile.get(s));
			}
		}
		if (jobs.size() == 0) {
			info.add("\t\tThe job was marked as STARTED but couldn't find the slurm-xxx.out file in "+jobDir);
			return false;
		}
		String slurm = null;
		if (jobs.size() == 1) slurm = jobs.get(0).getName();
		else {
			//multiple jobs take last
			File[] toSort = new File[jobs.size()];
			jobs.toArray(toSort);
			Arrays.sort(toSort);
			slurm = toSort[toSort.length-1].getName();
		}
		//pull the job id
		Matcher mat = slurmJobId.matcher(slurm);
		if (mat.matches() == false) throw new IOException("\tFailed to parse the job id from  "+slurm);
		String jobId = mat.group(1);
		//check the queue
		info.add("\t\tChecking the slurm queue for "+jobId);
		String[] res = IO.executeViaProcessBuilder(new String[]{"squeue", "-j", jobId}, false);
		if (res.length < 2) {
			new File(jobDir,"STARTED").delete();
			new File(jobDir,"FAILED").createNewFile();
			info.add("The job was marked as STARTED but failed to find the "+slurm+" job in the queue for "+jobDir+", marking FAILED.");
			return false;
		}
		info.add("\t\tSTARTED and in queue, "+jobDir+"\n\t\t\t"+res[1]);
		return true;
	}

	public static void cancelDeleteJobDir(HashMap<String, File> nameFile, File jobDir, ArrayList<String> info, boolean deleteJobDir){
		//send a scancel
		ArrayList<String> sb = new ArrayList<String>();
		for (String s: nameFile.keySet()) {
			if (s.startsWith("slurm")) {
				Matcher mat = slurmJobId.matcher(s);
				if (mat.matches()) {
					sb.add(mat.group(1));
					nameFile.get(s).delete();
				}
			}
		}
		if (sb.size() !=0){
			info.add("\t\tCanceling slurm jobs"+sb+ " from "+jobDir);
			sb.add(0, "scancel");
			IO.executeViaProcessBuilder(Misc.stringArrayListToStringArray(sb), false);
		}
		//delete dir
		if (deleteJobDir){
			info.add("\t\tDeleting "+jobDir);
			IO.deleteDirectory(jobDir);
			if (jobDir.exists()) IO.deleteDirectoryViaCmdLine(jobDir);
			jobDir.mkdirs();
		}
	}

	public static File copyInWorkflowDocs(File[] workflowDocs, File jobDir) throws IOException {
		File shellScript = null;
		for (File f: workflowDocs) {
			File copy = new File(jobDir, f.getName());
			IO.copyViaFileChannel(f, copy);
			if (copy.getName().endsWith(".sh")) shellScript = f;
		}
		if (shellScript == null) throw new IOException("Failed to find the workflow xxx.sh file in "+workflowDocs[0].getParent());
		return shellScript;
	}

	/**Attempts to find a xxx.bai file in the same dir as the xxx.bam file.*/
	public static final Pattern bamName = Pattern.compile("(.+)\\.bam");
	public static File fetchBamIndex(File bam) throws IOException{
		String path = bam.getCanonicalPath();
		Matcher mat = bamName.matcher(path);
		File index = null;
		if (mat.matches()) index = new File(mat.group(1)+".bai");
		if (index == null || index.exists() == false) throw new IOException("Failed to find the xxx.bai index file for "+bam);
		return index;
	}

	public static void removeProgressFiles(File jobDir) {
		new File(jobDir, "FAILED").delete();
		new File(jobDir, "COMPLETE").delete();
		new File(jobDir, "STARTED").delete();
		new File(jobDir, "QUEUED").delete();
	}

	private boolean checkFastq() throws IOException {
		info.add("Checking Fastq availability...");
		File fastqDir = makeCheckFile(rootDir, "Fastq");
		tumorDNAFastq = new FastqDataset(fastqDir, "TumorDNA", info);
		normalDNAFastq = new FastqDataset(fastqDir, "NormalDNA", info);
		tumorTransFastq = new FastqDataset(fastqDir, "TumorRNA", info);
		if (tumorDNAFastq.isFastqDirExists()) info.add("\tTumorDNA\t"+ (tumorDNAFastq.getFastqs() != null));
		if (normalDNAFastq.isFastqDirExists()) info.add("\tNormalDNA\t"+ (normalDNAFastq.getFastqs() != null));
		if (tumorTransFastq.isFastqDirExists()) info.add("\tTumorRNA\t"+ (tumorTransFastq.getFastqs() != null));
		if (tumorDNAFastq.getFastqs() != null || normalDNAFastq.getFastqs() != null || tumorTransFastq.getFastqs() != null) return true;
		return false;
	}

	public static File makeCheckFile(File parentDir, String fileName) throws IOException {
		File f = new File(parentDir, fileName);
		if (f.exists() == false) return null;
		return f;
	}
	
	/**Returns whether the job was restarted.*/
	private boolean checkJob(HashMap<String, File> nameFile, File jobDir, File[] toSoftLink, File[] runDocs) throws IOException {
		//force a restart?
		if (nameFile.containsKey("FAILED")){
			if (forceRestart  || restartFailedJobs ){
				//cancel any slurm jobs and delete the directory
				TNSample.cancelDeleteJobDir(nameFile, jobDir, info, true);
				restart(jobDir, toSoftLink, runDocs);
				return true;
			}
			else if (softRestart){
				//cancel any slurm jobs and delete the directory
				TNSample.cancelDeleteJobDir(nameFile, jobDir, info, false);
				restart(jobDir, toSoftLink, runDocs);
				return true;
			}
			//FAILED but no restart
			else if (nameFile.containsKey("FAILED")){
				info.add("\t\tFAILED "+jobDir);
				failed = true;
				return false;
			}
		}
		//QUEUED
		else if (nameFile.containsKey("QUEUED")){
			if (forceRestart) {
				//cancel any slurm jobs and delete the directory
				TNSample.cancelDeleteJobDir(nameFile, jobDir, info, true);
				restart(jobDir, toSoftLink, runDocs);
				return true;
			}
			info.add("\t\tQUEUED "+jobDir);
			running = true;
		}
		//STARTED
		else if (nameFile.containsKey("STARTED")) {
			if (TNSample.checkQueue(nameFile, jobDir, info) == false) failed = true;
			else running = true;
		}
		//hmm no status files, probably something went wrong on the cluster? mark it as FAILED
		else {
			info.add("\t\tMarking as FAILED, no job status files in "+jobDir);
			new File(jobDir, "FAILED").createNewFile();
			failed = true;
		}
		return false;
	}
	
	private void restart(File jobDir, File[] toSoftLink, File[] runDocs) throws IOException{
		//launch it
		launch(jobDir, toSoftLink, runDocs);
		new File(jobDir, "RESTARTED").createNewFile();
		info.add("\t\tRESTARTED");
	}
	
	private void launch(File jobDir, File[] toLink, File[] docs) throws IOException{
		info.add("\t\tLaunching "+jobDir);
		running = true;
		if (toLink != null) IO.createSymbolicLinks(toLink, jobDir);
		//replace any launch scripts with the current
		File shellScript = TNSample.copyInWorkflowDocs(docs, jobDir);
		//clear any progress files
		TNSample.removeProgressFiles(jobDir);
		//squeue the shell script
		new File(jobDir, "QUEUED").createNewFile();
		String alignDirPath = jobDir.getCanonicalPath();
		String[] output = IO.executeViaProcessBuilder(new String[]{"sbatch", "-J", alignDirPath.replace(tnRunner.getPathToTrim(), ""), "-D", alignDirPath, shellScript.getCanonicalPath()}, false);
		for (String o: output) info.add("\t\t"+o);
		numJobsLaunched++;
	}

	public FastqDataset getTumorDNAFastq() {
		return tumorDNAFastq;
	}

	public FastqDataset getNormalDNAFastq() {
		return normalDNAFastq;
	}

	public File getTumorRNABamBed() {
		return tumorRNABam;
	}

	public File[] getNormalDNABamBedGvcf() {
		return normalDNABamBedGvcf;
	}

	public File getRootDir() {
		return rootDir;
	}

	public String getId() {
		return id;
	}

	public boolean isFailed() {
		return failed;
	}

	public boolean isRunning() {
		return running;
	}

	public int getNumJobsLaunched() {
		return numJobsLaunched;
	}

	public boolean isFastqIssue() {
		return fastqIssue;
	}
}
