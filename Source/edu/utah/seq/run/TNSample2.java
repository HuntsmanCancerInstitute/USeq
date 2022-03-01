package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class TNSample2 {

	private String id;
	private File rootDir;
	private TNRunner2 tnRunner;
	private boolean restartFailedJobs;
	private static final Pattern slurmJobId = Pattern.compile("slurm-(\\d+).out");
	private ArrayList<String> info = new ArrayList<String>();
	private boolean deleteSampleConcordance = false;

	//Fastq
	private FastqCramDataset tumorDNAFastqCram = null;
	private FastqCramDataset normalDNAFastqCram = null;
	private FastqCramDataset tumorTransFastqCram = null;
	private boolean fastqIssue = false;

	//Alignment and passing region bed files
	private AlignmentDataset2 tumorDNAAlignment = null;
	private AlignmentDataset2 normalDNAAlignment = null;
	private AlignmentDataset2 tumorRNAAlignment = null;

	//Somatic variant Vcf calls
	private File somaticVariants = null;
	private File mergedSomaticVariants = null;

	//Joint genotyping
	private File[] germlineVcf = null;

	//Job status
	private boolean failed = false;
	private boolean running = false;
	private int numJobsLaunched = 0;


	public TNSample2(File rootDir, TNRunner2 tnRunner) throws IOException{
		this.rootDir = rootDir;
		this.tnRunner = tnRunner;
		this.restartFailedJobs = tnRunner.isRestartFailed();

		//assign id
		id = rootDir.getName();
		info.add("\n"+id);

		//look for fastq folders and check they contain 2 xxx.gz files, some may be null
		checkFastq();
		
		//RNA fusion, no alignments needed
		if (tnRunner.getRNAFusionDocs() != null) rnaFusionAnalysis();

		if (checkAlignments()) {

			//launch t/n somatic DNA analysis?
			if (tnRunner.getSomaticVarCallDocs() != null) somDNACall();

			//merge somatic vars with clinical test results?
			if (tnRunner.getClinicalVcfDocs()!= null && somaticVariants != null) parseMergeClinicalVars();

			//annotate somatic vcf?
			if (somaticVariants != null && tnRunner.getVarAnnoDocs() != null) annotateSomaticVcf();

			//sample concordance?
			if (tnRunner.getSampleConcordanceDocs() != null) sampleConcordance();
			
			//haplotype calling on normal?
			if (tnRunner.getHaplotypeCallDocs() != null) haplotypeCalling();

			//Annotate normal germline vcf
			if (tnRunner.getVarAnnoDocs() != null) {
				annotateGermlineVcf("GATK");
				annotateGermlineVcf("Illumina");
			}

			//copy ratio/ number
			if (tnRunner.getCopyRatioDocs() != null) copyRatioAnalysis();

			//msi
			if (tnRunner.getMsiDocs() != null) msiAnalysis();
		}

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
				info.add("\tJsonReport\tFAILED to find one xxx.json clinical test report file");
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
				info.add("\tCOMPLETE "+jobDir);
			}
			
			else checkJob(nameFile, jobDir, toLink, tnRunner.getClinicalVcfDocs());
		}

	public void annotateGermlineVcf(String name) throws IOException {
		info.add("Checking "+name+" germline variant annotation...");
		
		//look for genotyped vcf
		File dir = new File (rootDir, "GermlineVariantCalling/"+id+"_"+ name);
		if (dir.exists() == false) return;
		File[] res = IO.extractFiles(dir, ".vcf.gz");
		if (res == null || res.length !=1) return;
		else germlineVcf = new File[] {res[0], new File(dir, res[0].getName()+".tbi")};	

		//make dir, ok if it already exists
		File jobDir = new File (rootDir.getCanonicalFile(), "GermlineVariantCalling/"+id+"_"+name+"_Anno");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			IO.writeString(tnRunner.getGermlineAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
			IO.writeString(tnRunner.getGermlineVcfCallFreq(), new File(jobDir, "vcfCallFrequency.config.txt"));
			launch(jobDir, germlineVcf, tnRunner.getVarAnnoDocs());
		}

		//OK some files are present
		else if (nameFile.containsKey("COMPLETE")){
			//remove the linked vcfs and the config file
			new File(jobDir, germlineVcf[0].getName()).delete();
			new File(jobDir, germlineVcf[1].getName()).delete();
			info.add("\tCOMPLETE "+jobDir);
		}
		else {
			if (checkJob(nameFile, jobDir, germlineVcf, tnRunner.getVarAnnoDocs())){
				IO.writeString(tnRunner.getGermlineAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
				IO.writeString(tnRunner.getGermlineVcfCallFreq(), new File(jobDir, "vcfCallFrequency.config.txt"));
			}
		}
	}

	private void sampleConcordance() throws IOException {
		info.add("Checking sample concordance...");
		File jobDir = new File (rootDir, "SampleConcordance/"+id+"_SampleConcordance");

		//Started a new alignment? then delete the concordance if it exists
		if (deleteSampleConcordance && jobDir.exists()){
			info.add("\tNew alignments started.");
			cancelDeleteJobDir(IO.fetchNamesAndFiles(jobDir), jobDir, info, true);
			return;
		}
		
		//is there a tumor DNA alignment?
		if (tumorDNAFastqCram.isCramFastqDirExists()) {
			//what's the status of the alignment
			if (tumorDNAAlignment == null || tumorDNAAlignment.isComplete() == false) {
				info.add("\tWaiting for tumor DNA alignment.");
				return;
			}
		}
		
		//is there a normal DNA alignment?
		if (normalDNAFastqCram.isCramFastqDirExists()) {
			//what's the status of the alignment
			if (normalDNAAlignment == null || normalDNAAlignment.isComplete() == false) {
				info.add("\tWaiting for normal DNA alignment.");
				return;
			}
		}
		
		//is there a tumor RNA alignment?
		if (tumorTransFastqCram.isCramFastqDirExists()) {
			//what's the status of the alignment
			if (tumorRNAAlignment == null || tumorRNAAlignment.isComplete() == false) {
				info.add("\tWaiting for tumor RNA alignment.");
				return;
			}
		}
		
		//OK all alignments are complete see which is available for concordance
		//check to see if the alignments are present, sometimes these datasets are never coming
		boolean goT = false;
		boolean goN = false;
		boolean goTT = false;
		int numExist = 0;
		
		if (tumorDNAAlignment != null && tumorDNAAlignment.isComplete()) {
			goT = true;
			numExist++;
		}
		
		if (normalDNAAlignment != null && normalDNAAlignment.isComplete() && normalDNAAlignment.getCramFile().getName().contains("NA12878")==false) {
			goN = true;
			numExist++;
		}
		
		if (tumorRNAAlignment != null && tumorRNAAlignment.isComplete()) {
			goTT = true;
			numExist++;
		}
		
		//too few?
		if (numExist < 2) return;

		//make dir, ok if it already exists
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			removeSampleConcordanceLinks(jobDir);
			createSampleConcordanceLinks(jobDir, goT, goN, goTT);
			launch(jobDir, null, tnRunner.getSampleConcordanceDocs());
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			removeSampleConcordanceLinks(jobDir);
			info.add("\tCOMPLETE "+jobDir);
		}
		
		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getSampleConcordanceDocs())){
				removeSampleConcordanceLinks(jobDir);
				createSampleConcordanceLinks(jobDir, goT, goN, goTT);
			}
		}
		
	}

	private void createSampleConcordanceLinks(File jobDir, boolean goT, boolean goN, boolean goTT) throws IOException{
		File bPFileDir = new File (jobDir.getCanonicalFile(), "BamPileupFiles");
		bPFileDir.mkdir();
		if (goT) {
			Files.createSymbolicLink(new File(bPFileDir, "tumorDNA.bp.txt.gz").toPath(), tumorDNAAlignment.getBpPileupFile().toPath());
			Files.createSymbolicLink(new File(bPFileDir, "tumorDNA.bp.txt.gz.tbi").toPath(), tumorDNAAlignment.getBpIndexFile().toPath());
		}
		//don't include NA12878_NormalDNA_Hg38_final.bam, this is used as a mock normal 
		if (goN) {
			Files.createSymbolicLink(new File(bPFileDir, "normalDNA.bp.txt.gz").toPath(), normalDNAAlignment.getBpPileupFile().toPath());
			Files.createSymbolicLink(new File(bPFileDir, "normalDNA.bp.txt.gz.tbi").toPath(), normalDNAAlignment.getBpIndexFile().toPath());
		}
		if (goTT) {
			Files.createSymbolicLink(new File(bPFileDir, "tumorRNA.bp.txt.gz").toPath(), tumorRNAAlignment.getBpPileupFile().toPath());
			Files.createSymbolicLink(new File(bPFileDir, "tumorRNA.bp.txt.gz.tbi").toPath(), tumorRNAAlignment.getBpIndexFile().toPath());
		}
		//link in gender file
		//is it from Avatar?
		File[] info = IO.extractFiles(rootDir, "Info.json.gz");
		if (info == null || info.length !=1) {
			//try to fetch from Tempus/ Caris/ Foundation
			File d = new File(rootDir, "ClinicalReport");
			if (d.exists() == false) info= null;
			info = IO.extractFiles(d, ".json");
			if (info == null || info.length !=1) info = IO.extractFiles(d, ".xml");
			if (info == null || info.length !=1) info= null;
		}
		
		if (info == null) throw new IOException( "ERROR: failed to find a json file containing gender information for SampleConcordance. See TNRunner menu.");
		
		Files.createSymbolicLink(new File(jobDir, "gender."+info[0].getName()).toPath(), info[0].toPath());
	}

	private void removeSampleConcordanceLinks(File jobDir) throws IOException{
		File bPFileDir = new File (jobDir.getCanonicalFile(), "BamPileupFiles");
		IO.deleteDirectory(bPFileDir);
		File[] toDel = IO.extractFilesStartingWith(jobDir, "gender.");
		if (toDel != null && toDel.length ==1) toDel[0].delete();
		
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
		if (mergedSomaticVariants != null) toLink = new File[]{mergedSomaticVariants, new File(mergedSomaticVariants+".tbi")};
		else toLink = new File[]{somaticVariants, new File(somaticVariants+".tbi")};
		
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			IO.writeString(tnRunner.getSomaticAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
			IO.writeString(tnRunner.getSomaticVcfCallFreq(), new File(jobDir, "vcfCallFrequency.config.txt"));
			launch(jobDir, toLink, tnRunner.getVarAnnoDocs());
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//remove the linked vcfs and the config file
			for (File f: toLink) new File(jobDir, f.getName()).delete();
			new File(jobDir, "annotatedVcfParser.config.txt").delete();
			new File(jobDir, "vcfCallFrequency.config.txt").delete();
			info.add("\tCOMPLETE "+jobDir);
		}
		
		else {
			if (checkJob(nameFile, jobDir, toLink, tnRunner.getVarAnnoDocs())){
				IO.writeString(tnRunner.getSomaticAnnotatedVcfParser(), new File(jobDir, "annotatedVcfParser.config.txt"));
			}
		}
	}

	private void somDNACall() throws IOException {
		if (tumorDNAAlignment == null) return;
		if (normalDNAAlignment == null) {
			//still waiting for normal alignment
			if (normalDNAFastqCram.isCramFastqDirExists()) return;
			//no non matched so can't run som calling
			if (tnRunner.getNonMatchedNormal() == null) return;
		}
		
		//OK good to go
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
			info.add("\tCOMPLETE "+jobDir);
		}
		else {
			if (checkJob(nameFile,jobDir,null, tnRunner.getSomaticVarCallDocs())) createSomaticVariantLinks(jobDir);
		}

	}
	
	private void haplotypeCalling() throws IOException {
		
		info.add("Checking germline haplotype calling status...");
		
		//look for normal cram files
		if (normalDNAAlignment == null) {
			info.add("\tMissing N alignment files.");
			return;
		}

		//make dir, ok if it already exists
		File jobDir = new File (rootDir,"GermlineVariantCalling/"+id+"_GVCF");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			createHaplotypeCallingLinks(jobDir);
			launch(jobDir, null, tnRunner.getHaplotypeCallDocs());
		}
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final Vcf  file, Final/Vcfs/
			File vcfDir = new File(jobDir,"Vcfs");
			File[] res = IO.extractFiles(vcfDir, ".g.vcf.gz");
			if (res == null || res.length !=1) {
				clearAndFail(jobDir, "\tThe haplotype status calling was marked COMPLETE but failed to find the final xxx.g.vcf.gz file in "+vcfDir);
				return;
			}
			
			//remove the linked files
			removeHaplotypeCallingLinks(jobDir);
			info.add("\tCOMPLETE "+jobDir);
		}
		else {
			if (checkJob(nameFile,jobDir,null, tnRunner.getHaplotypeCallDocs())) createHaplotypeCallingLinks(jobDir);
		}
	}
	
	private void createHaplotypeCallingLinks(File jobDir) throws IOException {
		//remove any linked alignment files
		removeHaplotypeCallingLinks(jobDir);
		//cram files
		Path normalCram = normalDNAAlignment.getCramFile().toPath();
		Path normalCrai = normalDNAAlignment.getCramIndexFile().toPath();
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.cram").toPath(), normalCram);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.crai").toPath(), normalCrai);
	}
	
	private void removeHaplotypeCallingLinks(File jobDir) throws IOException{
		File f = jobDir.getCanonicalFile();
		new File(f, "normal.cram").delete();
		new File(f, "normal.crai").delete();
	}
	
	private void msiAnalysis() throws IOException {
		info.add("Checking MSI status analysis...");
		
		//look for tumor and normal bam files
		if (tumorDNAAlignment == null || normalDNAAlignment == null) {
			info.add("\tMissing one or both T/N alignment files.");
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
			info.add("\tCOMPLETE "+jobDir);
		}
		else {
			if (checkJob(nameFile,jobDir,null, tnRunner.getMsiDocs())) createMsiLinks(jobDir);
		}
	}

	private void rnaFusionAnalysis() throws IOException {
		info.add("Checking RNA Fusion status...");
		//look for tumor RNA fastq
		if (tumorTransFastqCram == null || tumorTransFastqCram.isGoodToAlign() == false) {
			info.add("\tMissing one or more of the tumor RNA fastq data");
			return;
		}

		//make dir, ok if it already exists
		File jobDir = new File (rootDir,"RNAAnalysis/"+id+"_STARFusion");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);

		//any files?
		if (nameFile.size() == 0) {
			createDNAAlignLinks(tumorTransFastqCram.getCramFastqs(), jobDir);
			launch(jobDir, null, tnRunner.getRNAFusionDocs());
		}

		//OK some files are present, complete?
		else if (nameFile.containsKey("COMPLETE")){
			//find the final tsv file
			File tsv = new File(jobDir, "Spreadsheets/star-fusion.fusion_predictions.abridged.coding_effect.tsv.gz");
			if (tsv.exists() == false) {
				clearAndFail(jobDir, "\tThe STAR Fusion workflow was marked COMPLETE but failed to find Spreadsheets/star-fusion.fusion_predictions.abridged.coding_effect.tsv.gz in "+jobDir);
				return;
			}
			
			//remove the linked files
			File f = jobDir.getCanonicalFile();
			new File(f, "1.fastq.gz").delete();
			new File(f,"2.fastq.gz").delete();
			new File(f, "rawSeq.cram").delete();
			info.add("\tCOMPLETE "+jobDir);
		}

		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getRNAFusionDocs())){
				createDNAAlignLinks(tumorTransFastqCram.getCramFastqs(), jobDir);
			}
		}
	}
	
	private void createMsiLinks(File jobDir) throws IOException {
		//remove any linked alignment files
		removeMsiLinks(jobDir);
		//look for bams, might not be present
		if (tumorDNAAlignment.getBamFile() != null) {
			Path tumorBam = tumorDNAAlignment.getBamFile().toPath();
			Path tumorBai = tumorDNAAlignment.getBamIndexFile().toPath();
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bam").toPath(), tumorBam);
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bai").toPath(), tumorBai);
		}
		else {
			//cram files
			Path tumorCram = tumorDNAAlignment.getCramFile().toPath();
			Path tumorCrai = tumorDNAAlignment.getCramIndexFile().toPath();
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.cram").toPath(), tumorCram);
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.crai").toPath(), tumorCrai);
		}
		if (normalDNAAlignment.getBamFile() != null) {
			Path normalBam = normalDNAAlignment.getBamFile().toPath();
			Path normalBai = normalDNAAlignment.getBamIndexFile().toPath();
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bam").toPath(), normalBam);
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bai").toPath(), normalBai);
		}
		else {
			//cram files
			Path normalCram = normalDNAAlignment.getCramFile().toPath();
			Path normalCrai = normalDNAAlignment.getCramIndexFile().toPath();
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.cram").toPath(), normalCram);
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.crai").toPath(), normalCrai);
		}
	}
	
	private void removeMsiLinks(File jobDir) throws IOException{
		File f = jobDir.getCanonicalFile();
		new File(f, "tumor.bam").delete();
		new File(f, "tumor.bai").delete();
		new File(f, "normal.bam").delete();
		new File(f, "normal.bai").delete();
		new File(f, "tumor.cram").delete();
		new File(f, "tumor.crai").delete();
		new File(f, "normal.cram").delete();
		new File(f, "normal.crai").delete();
	}

	private void copyRatioAnalysis() throws IOException {
		if (tnRunner.getCopyRatioDocs() == null) return;
		info.add("Checking copy ratio calling...");	
		
		//look for tumor normal alignments
		if (tumorDNAAlignment == null || normalDNAAlignment == null) return;
		
		//look for germline vcf
		if (germlineVcf == null) {
			File vcf = null;
			File vcfIndex = null;
			File gvc = new File (rootDir, "GermlineVariantCalling");
			if (gvc.exists() == false) gvc = new File (rootDir, "GermlineVariantCalls");
			if (gvc.exists()) {
				File[] dirs = IO.extractOnlyDirectories(gvc);
				for (File d: dirs) {
					String name = d.getName();
					if (name.endsWith("_GATK") || name.endsWith("_Illumina")) {
						File[] vcfs = IO.extractFiles(d, "_JointGenotyped.vcf.gz");
						vcf = vcfs[0];
						vcfIndex = new File (vcf.getCanonicalPath()+".tbi");
						if (vcfIndex.exists() == false) {
							failed = true;
							info.add("\tERROR: failed to find the germline vcf index file?! "+vcfIndex);
						}
						break;
					}
				}
			}
			if (vcf == null) {
				info.add("\tWaiting for germline vcf from "+rootDir+"/GermlineVariantCalling");
				return;
			};
			germlineVcf = new File[]{vcf, vcfIndex};
		}
		
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
			info.add("\tCOMPLETE "+jobDir);
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
		Path tumorCram = tumorDNAAlignment.getCramFile().toPath();
		Path tumorCrai = tumorDNAAlignment.getCramIndexFile().toPath();
		Path normalCram = normalDNAAlignment.getCramFile().toPath();
		Path normalCrai = normalDNAAlignment.getCramIndexFile().toPath();

		Path normalVcf = germlineVcf[0].toPath();
		Path normalTbi = germlineVcf[1].toPath();
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.cram").toPath(), tumorCram);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.crai").toPath(), tumorCrai);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.cram").toPath(), normalCram);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.crai").toPath(), normalCrai);
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
		new File(f, "tumor.cram").delete();
		new File(f, "tumor.crai").delete();
		new File(f, "bkgPoN.hdf5").delete();
		new File(f, "normal.cram").delete();
		new File(f, "normal.crai").delete();
		new File(f, "normal.vcf.gz").delete();
		new File(f, "normal.vcf.gz.tbi").delete();
	}
	
	/**Attempts to load the tumor and normal DNA File arrays and the tumor RNA*/
	private boolean checkAlignments() throws IOException {
		info.add("Checking alignments...");

		File alignDir = new File (rootDir, "Alignment");
		boolean complete = true;

		//load with align dirs, if they exist check launch them
		HashSet<String> keys = new HashSet<String>();
		
		if (alignDir.exists()) {
			for (File ads: IO.extractOnlyDirectories(alignDir)) {
				if (ads.getName().endsWith("_NormalDNA")) {
					if (keys.contains("NormalDNA")) throw new IOException ("ERROR: found >1 NormalDNA alignments in "+alignDir);
					keys.add("NormalDNA");
					normalDNAAlignment = DNAAlignQC(normalDNAFastqCram, ads);
					if (normalDNAAlignment.isComplete() == false) complete = false;
				}
				else if (ads.getName().endsWith("_TumorDNA")) {
					if (keys.contains("TumorDNA")) throw new IOException ("ERROR: found >1 TumorDNA alignments in "+alignDir);
					keys.add("TumorDNA");
					tumorDNAAlignment = DNAAlignQC(tumorDNAFastqCram, ads);
					if (tumorDNAAlignment.isComplete() == false) complete = false;
				}
				else if (ads.getName().endsWith("_TumorRNA")) {
					if (keys.contains("TumorRNA")) throw new IOException ("ERROR: found >1 TumorRNA alignments in "+alignDir);
					keys.add("TumorRNA");
					tumorRNAAlignment = RNAAlignQC(tumorTransFastqCram, ads);
					if (tumorRNAAlignment.isComplete() == false) complete = false;
				}
			}
		}

		//launch any not already running
		if (tnRunner.getDNAAlignQCDocs() != null) {
			if (tumorDNAFastqCram.isGoodToAlign() && keys.contains("TumorDNA") == false) {
				tumorDNAAlignment = DNAAlignQC(tumorDNAFastqCram, null);
				complete = false;
			}
			if (normalDNAFastqCram.isGoodToAlign() && keys.contains("NormalDNA") == false) {
				normalDNAAlignment = DNAAlignQC(normalDNAFastqCram, null);
				complete = false;
			}
		}
		
		if (tnRunner.getRNAAlignQCDocs() != null && tumorTransFastqCram.isGoodToAlign() && keys.contains("TumorRNA") == false) {
			tumorRNAAlignment = RNAAlignQC(tumorTransFastqCram, null);
			complete = false;
		}
		return complete;
	}

	/**For RNASeq Alignments and Picard CollectRNASeqMetrics*/
	private AlignmentDataset2 RNAAlignQC(FastqCramDataset fd, File jobDir) throws IOException {

		//make dir if null
		if (jobDir == null) {
			jobDir = new File (rootDir, "Alignment/"+id+"_"+fd.getName()).getCanonicalFile();
			jobDir.mkdirs();
		}

		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		AlignmentDataset2 ad = new AlignmentDataset2(jobDir, info, true);
		
		//any files?
		if (nameFile.size() == 0) {
			createDNAAlignLinks(fd.getCramFastqs(), jobDir);
			launch(jobDir, null, tnRunner.getRNAAlignQCDocs());
		}

		//OK some files are present, complete?
		else if (nameFile.containsKey("COMPLETE")){
			if (ad.isComplete() == false) clearAndFail(jobDir, null);
			//remove any linked fastq files
			File f = jobDir.getCanonicalFile();
			new File(f, "1.fastq.gz").delete();
			new File(f,"2.fastq.gz").delete();
			new File(f, "rawSeq.cram").delete();
		}
		
		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getRNAAlignQCDocs())){
				createDNAAlignLinks(fd.getCramFastqs(), jobDir);
			}
		}
		return ad;
	}
	
	private void clearAndFail(File jobDir, String infoLine) throws IOException{
		if (infoLine != null) info.add(infoLine);
		failed = true;
		removeProgressFiles(jobDir);
		new File(jobDir, "FAILED").createNewFile();
	}

	
	/**For DNAs, diff read coverage for passing bed generation
	 * @throws IOException */
	private AlignmentDataset2 DNAAlignQC(FastqCramDataset fd, File jobDir) throws IOException {

		//make dir if null
		if (jobDir == null) {
			jobDir = new File (rootDir, "Alignment/"+id+"_"+fd.getName()).getCanonicalFile();
			jobDir.mkdirs();
		}
		
		AlignmentDataset2 ad = new AlignmentDataset2(jobDir, info, false);
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		
		//no files so launch new alignment
		if (nameFile.size() == 0) {
			createDNAAlignLinks(fd.getCramFastqs(), jobDir);
			deleteSampleConcordance = true;
			launch(jobDir, null, tnRunner.getDNAAlignQCDocs());
		}

		//files present, is it complete?
		else if (nameFile.containsKey("COMPLETE")) {
			//were all the files found?
			if (ad.isComplete() == false) clearAndFail(jobDir, null);
			
			//remove any linked fastq or cram files
			File f = jobDir.getCanonicalFile();
			new File(f, "1.fastq.gz").delete();
			new File(f,"2.fastq.gz").delete();
			new File(f, "rawSeq.cram").delete();
		}
		
		//files present but not complete, check job
		else {
			if (checkJob(nameFile, jobDir, null, tnRunner.getDNAAlignQCDocs())){
				createDNAAlignLinks(fd.getCramFastqs(), jobDir);
				deleteSampleConcordance = true;
			}
		}
		return ad;
	}

	private void createDNAAlignLinks(File[] cramFastq, File alignDir) throws IOException {
		File f = alignDir.getCanonicalFile();
		Path real1 = cramFastq[0].toPath();
		//fastq?
		if (cramFastq.length == 2) {
			//remove any links
			new File(f, "1.fastq.gz").delete();
			new File(f,"2.fastq.gz").delete();
			//soft link in the new ones
			Path real2 = cramFastq[1].toPath();
			Files.createSymbolicLink(new File(alignDir.getCanonicalFile(),"1.fastq.gz").toPath(), real1);
			Files.createSymbolicLink(new File(alignDir.getCanonicalFile(),"2.fastq.gz").toPath(), real2);
		}
		//must be cram
		else {
			//remove any links
			new File(f, "rawSeq.cram").delete();
			//soft link in the new ones
			Files.createSymbolicLink(new File(alignDir.getCanonicalFile(),"rawSeq.cram").toPath(), real1);
		}
	}
	
	private void createSomaticVariantLinks(File jobDir) throws IOException {
		//remove any linked bam and bed files
		removeSomaticLinks(jobDir);

		//soft link in the new ones
		Path tumorCram = tumorDNAAlignment.getCramFile().toPath();
		Path tumorCrai = tumorDNAAlignment.getCramIndexFile().toPath();
		Path tumorBed = tumorDNAAlignment.getBedFile().toPath();
		
		Path normalCram = null;
		Path normalCrai = null;
		Path normalBed = null;
		// link non matched normal?
		if (normalDNAAlignment == null) {
			normalCram = tnRunner.getNonMatchedNormal()[0].toPath();
			normalCrai = fetchBamIndex(tnRunner.getNonMatchedNormal()[0]).toPath();
			normalBed = tnRunner.getNonMatchedNormal()[1].toPath();
			//indicate this is non matched
			File f = new File (rootDir,"SomaticVariantCalls/NON_MATCHED_NORMAL_GERMLINE_CONTAMINATION");
			if (f.exists()== false) f.createNewFile();
		}
		else {
			normalCram = normalDNAAlignment.getCramFile().toPath();
			normalCrai = normalDNAAlignment.getCramIndexFile().toPath();
			normalBed = normalDNAAlignment.getBedFile().toPath();
		}
		
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.cram").toPath(), tumorCram);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.crai").toPath(), tumorCrai);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"tumor.bed.gz").toPath(), tumorBed);
		Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bed.gz").toPath(), normalBed);
		if (normalDNAAlignment == null) {
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bam").toPath(), normalCram);
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.bai").toPath(), normalCrai);
		}
		else {
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.cram").toPath(), normalCram);
			Files.createSymbolicLink(new File(jobDir.getCanonicalFile(),"normal.crai").toPath(), normalCrai);
		}
	}

	private void removeSomaticLinks(File jobDir) throws IOException{
		File f = jobDir.getCanonicalFile();
		new File(f, "tumor.cram").delete();
		new File(f, "tumor.crai").delete();
		new File(f, "tumor.bed.gz").delete();
		new File(f, "normal.cram").delete();
		new File(f, "normal.crai").delete();
		new File(f, "normal.bed.gz").delete();
		new File(f, "normal.bam").delete();
		new File(f, "normal.bai").delete();
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
			info.add("\tThe job was marked as STARTED but couldn't find the slurm-xxx.out file in "+jobDir);
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
		String[] res = IO.executeViaProcessBuilder(new String[]{"squeue", "-j", jobId}, false);
		if (res.length < 2) {
			new File(jobDir,"STARTED").delete();
			new File(jobDir,"FAILED").createNewFile();
			info.add("The job was marked as STARTED but failed to find the "+slurm+" job in the queue for "+jobDir+", marking FAILED.");
			return false;
		}
		info.add("\tSTARTED and in queue, "+jobDir+"\t"+res[1]);
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
			info.add("\tCanceling slurm jobs"+sb+ " from "+jobDir);
			sb.add(0, "scancel");
			IO.executeViaProcessBuilder(Misc.stringArrayListToStringArray(sb), false);
		}
		//delete dir
		if (deleteJobDir){
			info.add("\tDeleting "+jobDir);
			IO.deleteDirectory(jobDir);
			if (jobDir.exists()) IO.deleteDirectoryViaCmdLine(jobDir);
			jobDir.mkdirs();
		}
	}

	public static File copyInWorkflowDocs(File[] workflowDocs, File jobDir) throws IOException {
		File shellScript = null;
		for (File f: workflowDocs) {
			if (f.isDirectory()) continue;
			File copy = new File(jobDir, f.getName());
			IO.copyViaFileChannel(f, copy);
			if (copy.getName().endsWith(".sh")) shellScript = f;
		}
		if (shellScript == null) throw new IOException("Failed to find the workflow xxx.sh file in "+workflowDocs[0].getParent());
		return shellScript;
	}

	/**Attempts to find a xxx.bai or xxx.crai file in the same dir as the xxx.bam/cram file.*/
	public static File fetchBamIndex(File bam) throws IOException{
		String path = bam.getCanonicalPath();
		String tPath = path.substring(0, path.length()-4);
		File index = new File (tPath+".bai");
		if (index.exists()) return index;
		index = new File (path+".bai");
		if (index.exists()) return index;
		index = new File (tPath+".crai");
		if (index.exists()) return index;
		index = new File (path+".crai");
		if (index.exists() == false) throw new IOException("Failed to find the xxx.bai/crai index file for "+bam);
		return index;
	}

	public static void removeProgressFiles(File jobDir) {
		new File(jobDir, "FAILED").delete();
		new File(jobDir, "COMPLETE").delete();
		new File(jobDir, "STARTED").delete();
		new File(jobDir, "QUEUED").delete();
	}

	private void checkFastq() throws IOException {
		info.add("Checking fastq availability...");
		File fastqDir = makeCheckFile(rootDir, "Fastq");
		tumorDNAFastqCram = new FastqCramDataset(fastqDir, "TumorDNA", info);
		normalDNAFastqCram = new FastqCramDataset(fastqDir, "NormalDNA", info);
		tumorTransFastqCram = new FastqCramDataset(fastqDir, "TumorRNA", info);
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
			if (restartFailedJobs ){
				//cancel any slurm jobs and delete the directory
				TNSample2.cancelDeleteJobDir(nameFile, jobDir, info, true);
				restart(jobDir, toSoftLink, runDocs);
				return true;
			}
			//FAILED but no restart
			else if (nameFile.containsKey("FAILED")){
				info.add("\tFAILED "+jobDir);
				failed = true;
				return false;
			}
		}
		//QUEUED
		else if (nameFile.containsKey("QUEUED")){
			info.add("\tQUEUED "+jobDir);
			running = true;
		}
		//RUNME
		else if (nameFile.containsKey("RUNME")){
			info.add("\tRUNME "+jobDir);
			running = true;
		}
		//STARTED
		else if (nameFile.containsKey("STARTED")) {
			if (TNSample2.checkQueue(nameFile, jobDir, info) == false) failed = true;
			else running = true;
		}
		//hmm no status files, probably something went wrong on the cluster? mark it as FAILED
		else {
			info.add("\tMarking as FAILED, no job status files in "+jobDir);
			new File(jobDir, "FAILED").createNewFile();
			failed = true;
		}
		return false;
	}
	
	private void restart(File jobDir, File[] toSoftLink, File[] runDocs) throws IOException{
		//launch it
		launch(jobDir, toSoftLink, runDocs);
		new File(jobDir, "RESTARTED").createNewFile();
		info.add("\tRESTARTED");
	}
	
	private void launch(File jobDir, File[] toLink, File[] docs) throws IOException{
		if (tnRunner.isSbatch()) info.add("\tLAUNCHING "+jobDir);
		else info.add("\tSETTING UP "+jobDir);

		running = true;
		if (toLink != null) IO.createSymbolicLinks(toLink, jobDir);
		//replace any launch scripts with the current
		File shellScript = TNSample2.copyInWorkflowDocs(docs, jobDir);
		//clear any progress files
		TNSample2.removeProgressFiles(jobDir);

		if (tnRunner.isSbatch()) {
			//squeue the shell script
			new File(jobDir, "QUEUED").createNewFile();
			String alignDirPath = jobDir.getCanonicalPath();
			String[] output = IO.executeViaProcessBuilder(new String[]{"sbatch", "-J", alignDirPath.replace(tnRunner.getPathToTrim(), ""), "-D", alignDirPath, shellScript.getCanonicalPath()}, false);
			for (String o: output) info.add("\t\t"+o);
			numJobsLaunched++;
		}
		else new File(jobDir, "RUNME").createNewFile();
	}

	public FastqCramDataset getTumorDNAFastq() {
		return tumorDNAFastqCram;
	}

	public FastqCramDataset getNormalDNAFastq() {
		return normalDNAFastqCram;
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


	public AlignmentDataset2 getNormalDNADataset() {
		return normalDNAAlignment;
	}
}
