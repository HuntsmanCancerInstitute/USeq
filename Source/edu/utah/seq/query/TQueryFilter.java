package edu.utah.seq.query;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TQueryFilter {
	
	/**This hash is loaded with the name of the parent data directories observed with each data file.
	 * Make this something descriptive like, foundation, h1k, or exome, etc. Users can select 1 or more to restrict 
	 * what data is returned. */
	private HashSet<String> availableDataTypes = new HashSet<String>();
	private HashSet<String> dataTypesToReturn = new HashSet<String>();
	boolean filterOnDataTypes = false;
	
	/**This hash is loaded with the File data sources. */
	private HashSet<File> availableDataFiles = new HashSet<File>();
	private HashSet<File> dataFilesToReturn = new HashSet<File>();
	boolean filterOnDataFiles = false;
	
	/**This hash contains the file types/ extensions found, just two so far: vcf and bed.*/
	private HashSet<String> availableFileTypes = new HashSet<String>();
	private HashSet<String> fileTypesToReturn = new HashSet<String>();
	boolean filterOnFileTypes = false;
	
	public static final Pattern FILE_TYPE = Pattern.compile(".+\\.(.+)\\.gz$");
	
	//constructors
	public TQueryFilter(){}
	
	
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
		sb.append("Available Filters:");
		sb.append("\n\tData Type: "+fetchDataTypes(", "));
		sb.append("\n\tFile Type: "+fetchFileTypes(", "));
		sb.append("\n\tData Source: "+fetchDataFiles(", "));
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
			System.out.println("Num Files Pre : "+fileTabixQueries.size());
			fileTabixQueries.keySet().retainAll(toKeep);
			System.out.println("Num Files Post: "+fileTabixQueries.size());
		}
	}
	
	//getters and setters
	public HashSet<String> getAvailableDataTypes() {
		return availableDataTypes;
	}

	public void setAvailableDataTypes(HashSet<String> availableDataTypes) {
		this.availableDataTypes = availableDataTypes;
	}

	public HashSet<String> getDataTypesToReturn() {
		return dataTypesToReturn;
	}

	public void setDataTypesToReturn(HashSet<String> dataTypesToReturn) {
		this.dataTypesToReturn = dataTypesToReturn;
	}

	public boolean isFilterOnDataTypes() {
		return filterOnDataTypes;
	}

	public void setFilterOnDataTypes(boolean filterOnDataTypes) {
		this.filterOnDataTypes = filterOnDataTypes;
	}

	public HashSet<String> getAvailableFileTypes() {
		return availableFileTypes;
	}

	public void setAvailableFileTypes(HashSet<String> availableFileTypes) {
		this.availableFileTypes = availableFileTypes;
	}

	public HashSet<String> getFileTypesToReturn() {
		return fileTypesToReturn;
	}

	public void setFileTypesToReturn(HashSet<String> fileTypesToReturn) {
		this.fileTypesToReturn = fileTypesToReturn;
	}

	public boolean isFilterOnFileTypes() {
		return filterOnFileTypes;
	}

	public void setFilterOnFileTypes(boolean filterOnFileTypes) {
		this.filterOnFileTypes = filterOnFileTypes;
	}

	public HashSet<File> getAvailableDataFiles() {
		return availableDataFiles;
	}

	public void setAvailableDataFiles(HashSet<File> availableDataFiles) {
		this.availableDataFiles = availableDataFiles;
	}

	public HashSet<File> getDataFilesToReturn() {
		return dataFilesToReturn;
	}

	public void setDataFilesToReturn(HashSet<File> dataFilesToReturn) {
		this.dataFilesToReturn = dataFilesToReturn;
	}

	public boolean isFilterOnDataFiles() {
		return filterOnDataFiles;
	}

	public void setFilterOnDataFiles(boolean filterOnDataFiles) {
		this.filterOnDataFiles = filterOnDataFiles;
	}




	
}
