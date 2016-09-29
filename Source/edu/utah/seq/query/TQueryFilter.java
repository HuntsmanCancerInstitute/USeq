package edu.utah.seq.query;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class TQueryFilter {
	
	/**This hash is loaded with the name of the parent data directories observed with each data file.
	 * Make this something descriptive like, foundation, h1k, or exome, etc. Users can select 1 or more to restrict 
	 * what data is returned. */
	private TreeSet<String> availableDataTypes = new TreeSet<String>();
	private TreeSet<String> dataTypesToReturn = new TreeSet<String>();
	boolean filterOnDataTypes = false;
	
	/**This hash is loaded with the File data sources. */
	private TreeSet<File> availableDataFiles = new TreeSet<File>();
	private TreeSet<File> dataFilesToReturn = new TreeSet<File>();
	boolean filterOnDataFiles = false;
	
	/**This hash contains the file types/ extensions found, just two so far: vcf and bed.*/
	private TreeSet<String> availableFileTypes = new TreeSet<String>();
	private TreeSet<String> fileTypesToReturn = new TreeSet<String>();
	boolean filterOnFileTypes = false;
	public static final Pattern FILE_TYPE = Pattern.compile(".+\\.(.+)\\.gz$");
	
	private String[] trimmedDataFileNames = null;
	private String pathToTrimmedFile = null;
	private TQuery tQuery = null;
	
	//constructors
	public TQueryFilter(TQuery tQuery){
		this.tQuery = tQuery;
	}
	
	/**Turns off all the filters and loads the *ToReturn to all available*/
	public void reset() {
		filterOnDataTypes = false;
		filterOnDataFiles = false;
		filterOnFileTypes = false;
		dataTypesToReturn.addAll(availableDataTypes);
		dataFilesToReturn.addAll(availableDataFiles);
		fileTypesToReturn.addAll(availableFileTypes);
	}
	
	public String fetchSummary(){
		StringBuilder sb = new StringBuilder();
		sb.append("Available Filters and Cmds:");
		sb.append("\n(Note, for key=val filters, use key=val1,val2,... no spaces. Only data\nsources that pass all of the filters will be included in the output.)");

		sb.append("\n\n\treset\tsets all filters to default");
		sb.append("\n\tprintFilters\tprints current filters");
		sb.append("\n\tfetchData = true or false to pull data from disk");
		sb.append("\n\tprintJson = true or false to print results to stdout");
		sb.append("\n\tdataTypes = "+fetchDataTypes(","));
		sb.append("\n\tfileTypes = "+fetchFileTypes(","));
		sb.append("\n\tdataSources = "+fetchDataFilesRelative(","));
		return sb.toString();
	}
	
	public String getFilters(String delimiter){
		StringBuilder sb = new StringBuilder();
		sb.append("fetchData="+tQuery.isFetchData());
		sb.append(delimiter);
		sb.append("printJson="+tQuery.isPrintJson());
		
		if (filterOnDataTypes){
			sb.append(delimiter);
			sb.append("dataTypes=");
			Iterator it = dataTypesToReturn.iterator();
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
	public String fetchDataTypes(String delimiter){
		Iterator it = availableDataTypes.iterator();
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
		if (filterOnDataTypes || filterOnDataFiles|| filterOnFileTypes){
			Matcher mat;
			ArrayList<File> toKeep = new ArrayList<File>();
			//for each file
			for (File f: fileTabixQueries.keySet()){
				boolean addIt = true;
				
				//filter on data types, parent directory name
				if (filterOnDataTypes) addIt = dataTypesToReturn.contains(f.getParentFile().getName());
				
				//filter on data sources
				if (addIt && filterOnDataFiles) addIt = dataFilesToReturn.contains(f);
				
				//filter on file types/ extensions
				if (addIt && filterOnFileTypes ){
					mat = FILE_TYPE.matcher(f.getName());
					mat.matches();
					addIt = fileTypesToReturn.contains(mat.group(1));
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
	public TreeSet<String> getAvailableDataTypes() {
		return availableDataTypes;
	}

	public void setAvailableDataTypes(TreeSet<String> availableDataTypes) {
		this.availableDataTypes = availableDataTypes;
	}

	public TreeSet<String> getDataTypesToReturn() {
		return dataTypesToReturn;
	}

	public void setDataTypesToReturn(TreeSet<String> dataTypesToReturn) {
		this.dataTypesToReturn = dataTypesToReturn;
	}

	public boolean isFilterOnDataTypes() {
		return filterOnDataTypes;
	}

	public void setFilterOnDataTypes(boolean filterOnDataTypes) {
		this.filterOnDataTypes = filterOnDataTypes;
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




	
}
