package edu.utah.seq.query;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

import util.gen.IO;
import util.gen.Misc;

public class DataSources {

	//fields related to filtering files that intersect a particular set of regions
	private TreeSet<String> availablePrimaryDataTypes = new TreeSet<String>();
	private TreeSet<String> availableSecondaryDataTypes = new TreeSet<String>();
	private TreeSet<File> availableDataFiles = new TreeSet<File>();
	private TreeSet<String> availableFileTypes = new TreeSet<String>();
	
	private String[] trimmedDataFileNames = null;
	private HashMap<File, String> dataFileDisplayName = new HashMap<File, String>();
	private String pathToTrimmedFile = null;
	private long recordsLoaded = 0;
	private long recordsSkipped = 0;
	
	//methods
	/**Loads the filter with the file, extension, parent, and grand parent dir name*/
	void addFileToFilter(File file) {
		//add to data sources
		availableDataFiles.add(file);

		//add parent dir name
		String pn = file.getParentFile().getName();
		if (pn != null && pn.trim().length()!=0)availablePrimaryDataTypes.add(pn);

		//add grand parent dir name
		String gpn = file.getParentFile().getParentFile().getName();
		if (gpn != null && gpn.trim().length()!=0) availableSecondaryDataTypes.add(gpn);

		//add extension minus .gz
		String name = file.getName();
		if (name.endsWith(".vcf.gz")) availableFileTypes.add("vcf");
		else if (name.endsWith(".bed.gz")) availableFileTypes.add("bed");
		else if (name.endsWith(".maf.txt.gz")) availableFileTypes.add("maf");

		//should never hit this
		//TODO: need to convert this to an Exception! 
		else Misc.printErrAndExit("\nERROR: unrecognized file type contact admin! "+file);
	}
	
	
	/**Returns all the data types (parent dir names).*/
	public String fetchPrimaryDataTypes(String delimiter){
		Iterator it = availablePrimaryDataTypes.iterator();
		return buildStringFromIterator(it, delimiter);
	}

	/**Returns all the data types (grand parent dir names).*/
	public String fetchSecondaryDataTypes(String delimiter){
		Iterator it = availableSecondaryDataTypes.iterator();
		return buildStringFromIterator(it, delimiter);
	}
	
	/**Returns all the data files, full path.*/
	public String fetchDataFiles(String delimiter){
		Iterator it = availableDataFiles.iterator();
		return buildStringFromIterator(it, delimiter);
	}
	
	/**Returns all the data files, full path.*/
	public String fetchDataFilesRelative(String delimiter){
		if (trimmedDataFileNames == null){
			Iterator<File> it = availableDataFiles.iterator();
			ArrayList<File> fullPath = new ArrayList<File>();
			while (it.hasNext()) fullPath.add(it.next());
			File[] f = new File[fullPath.size()];
			fullPath.toArray(f);
			trimmedDataFileNames = IO.trimCommonParentDirs(f);
			//set path with divider
			int index = f[0].toString().indexOf(trimmedDataFileNames[0]);
			pathToTrimmedFile = f[0].toString().substring(0, index);
			//create hash for easy lookup
			for (int i=0; i< f.length; i++) dataFileDisplayName.put(f[i], trimmedDataFileNames[i]);
		}
		StringBuilder sb = new StringBuilder(trimmedDataFileNames[0]);
		for (int i=1; i< trimmedDataFileNames.length; i++){
			sb.append(delimiter);
			sb.append(trimmedDataFileNames[i]);
		}
		return sb.toString();
	}
	
	/**Returns all the file types (extensions minus the .gz so bed, vcf, others as they are implemented).*/
	public String fetchFileTypes(String delimiter){
		Iterator it = availableFileTypes.iterator();
		return buildStringFromIterator(it, delimiter);
	}
	
	public String fetchSummary(){
		StringBuilder sb = new StringBuilder();
		sb.append("Available Filters and Cmds:");
		sb.append("\n\tPrintFilters\tprints current default filters");
		sb.append("\n\tFetchData = true/false to pull data from disk, default true");
		sb.append("\n\tMatchVcf = true/false to require the chr, pos, ref, and at least one alt to match, default false");
		sb.append("\n\tPrintJson = true/false to print results to stdout, default true");
		sb.append("\n\tPrimaryDataTypes = "+fetchPrimaryDataTypes(","));
		sb.append("\n\tSecondaryDataTypes = "+fetchSecondaryDataTypes(","));
		sb.append("\n\tFileTypes = "+fetchFileTypes(","));
		sb.append("\n\tDataSources = "+fetchDataFilesRelative(","));
		sb.append("\n\n(Note, for key=val filters, use key=val1,val2,... no spaces. Only data\nsources that pass all of the filters will be included in the output.)");

		return sb.toString();
	}
	
	public static String buildStringFromIterator(Iterator<Object> it, String delimiter){
		StringBuilder sb = new StringBuilder();
		sb.append(it.next().toString());
		while (it.hasNext()){
			sb.append(delimiter);
			sb.append(it.next().toString());
		}
		return sb.toString();
	}

	
	//getters and setters
	public TreeSet<String> getAvailablePrimaryDataTypes() {
		return availablePrimaryDataTypes;
	}
	public void setAvailablePrimaryDataTypes(TreeSet<String> availablePrimaryDataTypes) {
		this.availablePrimaryDataTypes = availablePrimaryDataTypes;
	}
	public TreeSet<String> getAvailableSecondaryDataTypes() {
		return availableSecondaryDataTypes;
	}
	public void setAvailableSecondaryDataTypes(TreeSet<String> availableSecondaryDataTypes) {
		this.availableSecondaryDataTypes = availableSecondaryDataTypes;
	}
	public TreeSet<File> getAvailableDataFiles() {
		return availableDataFiles;
	}
	public void setAvailableDataFiles(TreeSet<File> availableDataFiles) {
		this.availableDataFiles = availableDataFiles;
	}
	public TreeSet<String> getAvailableFileTypes() {
		return availableFileTypes;
	}
	public void setAvailableFileTypes(TreeSet<String> availableFileTypes) {
		this.availableFileTypes = availableFileTypes;
	}
	public String[] getTrimmedDataFileNames() {
		return trimmedDataFileNames;
	}
	public void setTrimmedDataFileNames(String[] trimmedDataFileNames) {
		this.trimmedDataFileNames = trimmedDataFileNames;
	}
	public HashMap<File, String> getDataFileDisplayName() {
		return dataFileDisplayName;
	}
	public void setDataFileDisplayName(HashMap<File, String> dataFileDisplayName) {
		this.dataFileDisplayName = dataFileDisplayName;
	}
	public String getPathToTrimmedFile() {
		return pathToTrimmedFile;
	}
	public void setPathToTrimmedFile(String pathToTrimmedFile) {
		this.pathToTrimmedFile = pathToTrimmedFile;
	}
	public long getRecordsLoaded() {
		return recordsLoaded;
	}
	public void setRecordsLoaded(long recordsLoaded) {
		this.recordsLoaded = recordsLoaded;
	}
	public long getRecordsSkipped() {
		return recordsSkipped;
	}
	public void setRecordsSkipped(long recordsSkipped) {
		this.recordsSkipped = recordsSkipped;
	}
	
}
