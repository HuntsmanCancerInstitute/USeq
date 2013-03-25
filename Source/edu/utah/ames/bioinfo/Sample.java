package edu.utah.ames.bioinfo;

import java.io.File;
import java.util.HashSet;

/**
 * This class holds info about a sequenced sample. Use these methods to populate a cmd.txt
 * file for Tomato alignment and to create an analysis report in GNomEx. If tabbed content
 * of the input autoalign report changes, just modify the index numbers of the fields below.
 * 
 * @author darren.ames@hci.utah.edu
 *
 */
public class Sample {
	
	//fields from autoalign report
	private String runID;
	private String lane;
	private String requestNum; 
	private String requester; 
	private String lab;
	private String sequencingApplicationCode; 
	private String singleOrPairedEnd;
	private String numberOfCycles;
	private String sampleName;
	private String organism;
	private String genome; 
	private String requesterEmail;
	private String experimentFolder;
	private String requestYear;
	private String projectName;
	private String sampleID;
	private String sequenceLaneNumber;
	private String fastqFileName;
	private String adapterSequence;
	
	//index fields
	public static final int laneIndex = 8;
	public static final int requestNumIndex = 2; 
	public static final int requesterIndex = 1;
	public static final int labIndex = 0;
	public static final int sequencingApplicationCodeIndex = 9;
	public static final int singleOrPairedEndIndex = 11;
	public static final int numberOfCyclesIndex = 10;
	public static final int organismIndex = 6;
	public static final int genomeIndex = 7;
	public static final int projectNameIndex = 3;
	public static final int sampleIDIndex = 4;
	public static final int sampleNameIndex = 5;
	public static final int requestYearIndex = 12;
	public static final int requesterEmailIndex = 13;
	public static final int sequenceLaneNumberIndex = 14;
	public static final int fastqFileNameIndex = 15;
	public static final int adapterSequenceIndex = 16;
	
	//worker fields
	private String sequencingApplication;
	private String versionedGenome;
	private String labNum;
	private String buildCode;
	private String sampleInfo;
	private String fastqFilePath;
	private String analysisNumber;
	private String analysisNumberPath;
	private String cmdFilePath;
	private String novoindex;
	private String params;
	private String indexCode;
	private File[] paramsFile;
	private File[] fastqFiles;
	private HashSet<String> emails = new HashSet<String>();
	private boolean isPairedEnd = false;
	private boolean isBisSeq = false;
	private boolean isRNASeq = false;
	private boolean isChIPSeq = false;
	private boolean isSmallRNA = false;
	private boolean isNewGenome = false;
	private boolean sequenceFilesAreNew = false;
	private boolean hasGenomeBuild = false;
	private boolean hasNovoindex = false;
	private boolean reverseStrand = false;
	
	//constructor
	public Sample(String[] dataValue) {

		lane = dataValue[laneIndex];
		requestNum = dataValue[requestNumIndex];
		requester = dataValue[requesterIndex];
		lab = dataValue[labIndex];
		sequencingApplicationCode = dataValue[sequencingApplicationCodeIndex];
		singleOrPairedEnd = dataValue[singleOrPairedEndIndex];
		numberOfCycles = dataValue[numberOfCyclesIndex];
		organism = dataValue[organismIndex];
		genome = dataValue[genomeIndex];
		sampleID = dataValue[sampleIDIndex];
		projectName = dataValue[projectNameIndex];
		requestYear = dataValue[requestYearIndex];
		requesterEmail = dataValue[requesterEmailIndex];
		sequenceLaneNumber = dataValue[sequenceLaneNumberIndex];
		fastqFileName = dataValue[fastqFileNameIndex];
		adapterSequence = dataValue[adapterSequenceIndex];
	}

	public String getFastqFilePath() {
		return fastqFilePath;
	}

	public void setFastqFilePath(String fastqFilePath) {
		this.fastqFilePath = fastqFilePath;
	}

	public String getCmdFilePath() {
		return cmdFilePath;
	}

	public void setCmdFilePath(String cmdFilePath) {
		this.cmdFilePath = cmdFilePath;
	}

	public File[] getFastqFiles() {
		return fastqFiles;
	}

	public void setFastqFiles(File[] fastqFiles) {
		this.fastqFiles = fastqFiles;
	}

	public HashSet<String> getEmails() {
		return emails;
	}

	public void setEmails(HashSet<String> emails) {
		this.emails = emails;
	}

	public boolean isBisSeq() {
		return isBisSeq;
	}

	public void setBisSeq(boolean isBisSeq) {
		this.isBisSeq = isBisSeq;
	}

	public boolean isRNASeq() {
		return isRNASeq;
	}

	public void setRNASeq(boolean isRNASeq) {
		this.isRNASeq = isRNASeq;
	}

	public boolean isChIPSeq() {
		return isChIPSeq;
	}

	public void setChIPSeq(boolean isChIPSeq) {
		this.isChIPSeq = isChIPSeq;
	}

	public boolean isSmallRNA() {
		return isSmallRNA;
	}

	public void setSmallRNA(boolean isSmallRNA) {
		this.isSmallRNA = isSmallRNA;
	}

	public boolean isSequenceFilesAreNew() {
		return sequenceFilesAreNew;
	}

	public void setSequenceFilesAreNew(boolean sequenceFilesAreNew) {
		this.sequenceFilesAreNew = sequenceFilesAreNew;
	}

	public boolean isHasGenomeBuild() {
		return hasGenomeBuild;
	}

	public void setHasGenomeBuild(boolean hasGenomeBuild) {
		this.hasGenomeBuild = hasGenomeBuild;
	}

	public boolean isHasNovoindex() {
		return hasNovoindex;
	}

	public void setHasNovoindex(boolean hasNovoindex) {
		this.hasNovoindex = hasNovoindex;
	}

	public File[] getParamsFile() {
		return paramsFile;
	}

	public void setParamsFile(File[] paramsFile) {
		this.paramsFile = paramsFile;
	}

	public String getNovoindex() {
		return novoindex;
	}

	public void setNovoindex(String novoindex) {
		this.novoindex = novoindex;
	}

	public boolean isNewGenome() {
		return isNewGenome;
	}

	public void setNewGenome(boolean isNewGenome) {
		this.isNewGenome = isNewGenome;
	}

	public String getRunID() {
		return runID;
	}

	public void setRunID(String runID) {
		this.runID = runID;
	}

	public String getLane() {
		return lane;
	}

	public void setLane(String lane) {
		this.lane = lane;
	}

	public String getRequestNumber() {
		return requestNum;
	}

	public void setRequestNumber(String requestNum) {
		this.requestNum = requestNum;
	}

	public String getRequester() {
		return requester;
	}

	public void setRequester(String requester) {
		this.requester = requester;
	}

	public String getLab() {
		return lab;
	}

	public void setLab(String lab) {
		this.lab = lab;
	}

	public String getSequencingApplicationCode() {
		return sequencingApplicationCode;
	}

	public void setSequencingApplicationCode(String sequencingApplicationCode) {
		this.sequencingApplicationCode = sequencingApplicationCode;
	}

	public String getSequencingApplication() {
		return sequencingApplication;
	}

	public void setSequencingApplication(String sequencingApplication) {
		this.sequencingApplication = sequencingApplication;
	}

	public String getSingleOrPairedEnd() {
		return singleOrPairedEnd;
	}

	public void setSingleOrPairedEnd(String singleOrPairedEnd) {
		this.singleOrPairedEnd = singleOrPairedEnd;
	}

	public String getNumberOfCycles() {
		return numberOfCycles;
	}

	public void setNumberOfCycles(String numberOfCycles) {
		this.numberOfCycles = numberOfCycles;
	}

	public String getOrganism() {
		return organism;
	}

	public void setOrganism(String organism) {
		this.organism = organism;
	}

	public String getGenome() {
		return genome;
	}

	public void setGenome(String genome) {
		this.genome = genome;
	}

	public String getExperimentFolder() {
		return experimentFolder;
	}

	public void setExperimentFolder(String experimentFolder) {
		this.experimentFolder = experimentFolder;
	}

	public String getProjectName() {
		return projectName;
	}

	public void setProjectName(String projectName) {
		this.projectName = projectName;
	}

	public String getSampleID() {
		return sampleID;
	}

	public void setSampleID(String sampleID) {
		this.sampleID = sampleID;
	}

	public String getRequestYear() {
		return requestYear;
	}

	public void setRequestYear(String requestYear) {
		this.requestYear = requestYear;
	}

	public String getSampleName() {
		return sampleName;
	}

	public void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}
	
	public String getLabNum() {
		return labNum;
	}
	
	public void setLabNum(String labNum) {
		this.labNum = labNum;
	}

	public String getSampleInfo() {
		return sampleInfo;
	}

	public void setSampleInfo(String sampleInfo) {
		this.sampleInfo = sampleInfo;
	}

	public String getBuildCode() {
		return buildCode;
	}

	public void setBuildCode(String buildCode) {
		this.buildCode = buildCode;
	}

	public String getParams() {
		return params;
	}

	public void setParams(String alignParams) {
		this.params = alignParams;
	}

	public String getIndexCode() {
		return indexCode;
	}

	public void setIndexCode(String indexCode) {
		this.indexCode = indexCode;
	}

	public String getAnalysisNumber() {
		return analysisNumber;
	}

	public void setAnalysisNumber(String analysisNumber) {
		this.analysisNumber = analysisNumber;
	}

	public String getRequesterEmail() {
		return requesterEmail;
	}

	public void setRequesterEmail(String requesterEmail) {
		this.requesterEmail = requesterEmail;
	}
	
	public String getSequenceLaneNumber() {
		return sequenceLaneNumber;
	}

	public void setSequenceLaneNumber(String sequenceLaneNumber) {
		this.sequenceLaneNumber = sequenceLaneNumber;
	}

	public String getAnalysisNumberPath() {
		return analysisNumberPath;
	}

	public void setAnalysisNumberPath(String analysisNumberPath) {
		this.analysisNumberPath = analysisNumberPath;
	}

	public String getVersionedGenome() {
		return versionedGenome;
	}

	public void setVersionedGenome(String versionedGenome) {
		this.versionedGenome = versionedGenome;
	}

	public boolean isPairedEnd() {
		return isPairedEnd;
	}

	public void setPairedEnd(boolean isPairedEnd) {
		this.isPairedEnd = isPairedEnd;
	}

	public boolean isReverseStrand() {
		return reverseStrand;
	}

	public void setReverseStrand(boolean reverseStrand) {
		this.reverseStrand = reverseStrand;
	}

	public String getFastqFileName() {
		return fastqFileName;
	}

	public void setFastqFileName(String fastqFileName) {
		this.fastqFileName = fastqFileName;
	}

	public String getAdapterSequence() {
		return adapterSequence;
	}

	public void setAdapterSequence(String adapterSequence) {
		this.adapterSequence = adapterSequence;
	}
}
