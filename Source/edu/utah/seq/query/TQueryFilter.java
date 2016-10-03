package edu.utah.seq.query;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.Iterator;
import util.gen.IO;
import util.gen.Misc;

public class TQueryFilter {
	
	/**This hash is loaded with the name of the parent data directories observed with each data file.
	 * Make this something descriptive like, foundation, h1k, or exome, etc. Users can select 1 or more to restrict 
	 * what data is returned. */
	private TreeSet<String> availablePrimaryDataTypes = new TreeSet<String>();
	private TreeSet<String> primaryDataTypesToReturn = new TreeSet<String>();
	boolean filterOnPrimaryDataTypes = false;
	
	/**This hash is loaded with the name of the grand parent data directories observed with each data file.
	 * Make this something descriptive like, foundation, h1k, or exome, etc. Users can select 1 or more to restrict 
	 * what data is returned. */
	private TreeSet<String> availableSecondaryDataTypes = new TreeSet<String>();
	private TreeSet<String> secondaryDataTypesToReturn = new TreeSet<String>();
	boolean filterOnSecondaryDataTypes = false;
	
	/**This hash is loaded with the File data sources. */
	private TreeSet<File> availableDataFiles = new TreeSet<File>();
	private TreeSet<File> dataFilesToReturn = new TreeSet<File>();
	boolean filterOnDataFiles = false;
	
	/**This hash contains the file types/ extensions found, just two so far: vcf and bed.*/
	private TreeSet<String> availableFileTypes = new TreeSet<String>();
	private TreeSet<String> fileTypesToReturn = new TreeSet<String>();
	boolean filterOnFileTypes = false;
	
	boolean matchVcf = false;
	
	private String[] trimmedDataFileNames = null;
	private HashMap<File, String> dataFileDisplayName = new HashMap<File, String>();
	private String pathToTrimmedFile = null;
	private TQuery tQuery = null;
	
	//constructors
	public TQueryFilter(TQuery tQuery){
		this.tQuery = tQuery;
	}
	
	/**Turns off all the filters and loads the *ToReturn to all available*/
	public void reset() {
		filterOnSecondaryDataTypes = false;
		filterOnPrimaryDataTypes = false;
		filterOnDataFiles = false;
		filterOnFileTypes = false;
		matchVcf = false;
		secondaryDataTypesToReturn.addAll(availableSecondaryDataTypes);
		primaryDataTypesToReturn.addAll(availablePrimaryDataTypes);
		dataFilesToReturn.addAll(availableDataFiles);
		fileTypesToReturn.addAll(availableFileTypes);
	}
	
	public String fetchSummary(){
		StringBuilder sb = new StringBuilder();
		sb.append("Available Filters and Cmds:");
		sb.append("\n(Note, for key=val filters, use key=val1,val2,... no spaces. Only data\nsources that pass all of the filters will be included in the output.)");

		sb.append("\n\nReset\tsets all filters to default");
		sb.append("\nPrintFilters\tprints current filters");
		sb.append("\nFetchData = true or false to pull data from disk");
		sb.append("\nMatchVcf = true or false to require the chr, pos, ref, and at least one alt to match");
		sb.append("\nPrintJson = true or false to print results to stdout");
		sb.append("\nPrimaryDataTypes = "+fetchPrimaryDataTypes(","));
		sb.append("\nSecondaryDataTypes = "+fetchSecondaryDataTypes(","));
		sb.append("\nFileTypes = "+fetchFileTypes(","));
		sb.append("\nDataSources = "+fetchDataFilesRelative(","));
		return sb.toString();
	}
	
	public String getCurrentFilters(String delimiter){
		StringBuilder sb = new StringBuilder();
		sb.append("fetchData="+tQuery.isFetchData());
		sb.append(delimiter);
		sb.append("printJson="+tQuery.isPrintJson());
		sb.append(delimiter);
		sb.append("matchVcf="+matchVcf);
		
		if (filterOnPrimaryDataTypes){
			sb.append(delimiter);
			sb.append("primaryDataTypes=");
			Iterator it = primaryDataTypesToReturn.iterator();
			sb.append(buildStringFromIterator(it, ","));
		}
		if (filterOnSecondaryDataTypes){
			sb.append(delimiter);
			sb.append("secondaryDataTypes=");
			Iterator it = secondaryDataTypesToReturn.iterator();
			sb.append(buildStringFromIterator(it, ","));
		}
		if (filterOnFileTypes){
			sb.append(delimiter);
			sb.append("fileTypes=");
			Iterator it = fileTypesToReturn.iterator();
			sb.append(buildStringFromIterator(it, ","));
		}
		if (filterOnDataFiles){
			sb.append(delimiter);
			sb.append("dataSources=");
			Iterator it = dataFilesToReturn.iterator();
			sb.append(buildStringFromIterator(it, ","));
		}
		return sb.toString();	
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
	
	private String buildStringFromIterator(Iterator<Object> it, String delimiter){
		StringBuilder sb = new StringBuilder();
		sb.append(it.next().toString());
		while (it.hasNext()){
			sb.append(delimiter);
			sb.append(it.next().toString());
		}
		return sb.toString();
	}
	
	public void filter(HashMap<File, ArrayList<TabixQuery>> fileTabixQueries){
		//filter?
		if (filterOnPrimaryDataTypes || filterOnSecondaryDataTypes || filterOnDataFiles|| filterOnFileTypes){
			ArrayList<File> toKeep = new ArrayList<File>();
			//for each file
			for (File f: fileTabixQueries.keySet()){
				boolean addIt = true;
				
				//filter on primary data types, parent directory name
				if (filterOnPrimaryDataTypes) addIt = primaryDataTypesToReturn.contains(f.getParentFile().getName());
				
				//filter on secondary data types, grand parent directory name
				if (addIt && filterOnSecondaryDataTypes) addIt = secondaryDataTypesToReturn.contains(f.getParentFile().getParentFile().getName());
				
				//filter on data sources
				if (addIt && filterOnDataFiles) addIt = dataFilesToReturn.contains(f);
				
				//filter on file types/ extensions
				if (addIt && filterOnFileTypes ){
					String name = f.getName();
					if (name.endsWith(".vcf.gz")) addIt = fileTypesToReturn.contains("vcf");
					else if (name.endsWith(".bed.gz")) addIt = fileTypesToReturn.contains("bed");
					else if (name.endsWith(".maf.txt.gz")) addIt = fileTypesToReturn.contains("maf");
					//should never hit this
					else Misc.printErrAndExit("\nERROR: unrecognized file type contact admin! "+f);
				}
				if (addIt) toKeep.add(f);
			}
			
		//OK toss those not in toKeep
			//System.out.println("Num Files Pre : "+fileTabixQueries.size());
			fileTabixQueries.keySet().retainAll(toKeep);
			//System.out.println("Num Files Post: "+fileTabixQueries.size());
		}
	}
	
	//getters and setters
	public TreeSet<String> getAvailablePrimaryDataTypes() {
		return availablePrimaryDataTypes;
	}
	public void setAvailablePrimaryDataTypes(TreeSet<String> availablePrimaryDataTypes) {
		this.availablePrimaryDataTypes = availablePrimaryDataTypes;
	}
	public TreeSet<String> getPrimaryDataTypesToReturn() {
		return primaryDataTypesToReturn;
	}
	public void setPrimaryDataTypesToReturn(TreeSet<String> primaryDataTypesToReturn) {
		this.primaryDataTypesToReturn = primaryDataTypesToReturn;
	}
	public boolean isFilterOnPrimaryDataTypes() {
		return filterOnPrimaryDataTypes;
	}
	public void setFilterOnPrimaryDataTypes(boolean filterOnPrimaryDataTypes) {
		this.filterOnPrimaryDataTypes = filterOnPrimaryDataTypes;
	}
	
	
	public TreeSet<String> getAvailableSecondaryDataTypes() {
		return availableSecondaryDataTypes;
	}
	public void setAvailableSecondaryDataTypes(TreeSet<String> a) {
		this.availableSecondaryDataTypes = a;
	}
	public TreeSet<String> getSecondaryDataTypesToReturn() {
		return secondaryDataTypesToReturn;
	}
	public void setSecondaryDataTypesToReturn(TreeSet<String> secondaryDataTypesToReturn) {
		this.secondaryDataTypesToReturn = secondaryDataTypesToReturn;
	}
	public boolean isFilterOnSecondaryDataTypes() {
		return filterOnSecondaryDataTypes;
	}
	public void setFilterOnSecondaryDataTypes(boolean filterOnSecondaryDataTypes) {
		this.filterOnSecondaryDataTypes = filterOnSecondaryDataTypes;
	}
	
	public TreeSet<String> getAvailableFileTypes() {
		return availableFileTypes;
	}

	public void setAvailableFileTypes(TreeSet<String> availableFileTypes) {
		this.availableFileTypes = availableFileTypes;
	}

	public TreeSet<String> getFileTypesToReturn() {
		return fileTypesToReturn;
	}

	public void setFileTypesToReturn(TreeSet<String> fileTypesToReturn) {
		this.fileTypesToReturn = fileTypesToReturn;
	}

	public boolean isFilterOnFileTypes() {
		return filterOnFileTypes;
	}

	public void setFilterOnFileTypes(boolean filterOnFileTypes) {
		this.filterOnFileTypes = filterOnFileTypes;
	}

	public TreeSet<File> getAvailableDataFiles() {
		return availableDataFiles;
	}

	public void setAvailableDataFiles(TreeSet<File> availableDataFiles) {
		this.availableDataFiles = availableDataFiles;
	}

	public TreeSet<File> getDataFilesToReturn() {
		return dataFilesToReturn;
	}

	public void setDataFilesToReturn(TreeSet<File> dataFilesToReturn) {
		this.dataFilesToReturn = dataFilesToReturn;
	}

	public boolean isFilterOnDataFiles() {
		return filterOnDataFiles;
	}

	public void setFilterOnDataFiles(boolean filterOnDataFiles) {
		this.filterOnDataFiles = filterOnDataFiles;
	}


	public String getPathToTrimmedFile() {
		return pathToTrimmedFile;
	}

	public HashMap<File, String> getDataFileDisplayName() {
		return dataFileDisplayName;
	}

	public boolean isMatchVcf() {
		return matchVcf;
	}

	public void setMatchVcf(boolean matchVcf) {
		this.matchVcf = matchVcf;
	}




	
}
