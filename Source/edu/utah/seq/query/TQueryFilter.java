package edu.utah.seq.query;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.Iterator;
import util.gen.Misc;

public class TQueryFilter {
	
	/**This hash is loaded with the name of the parent data directories observed with each data file.
	 * Make this something descriptive like, foundation, h1k, or exome, etc. Users can select 1 or more to restrict 
	 * what data is returned. */
	private TreeSet<String> primaryDataTypesToReturn = new TreeSet<String>();
	boolean filterOnPrimaryDataTypes = false;
	
	/**This hash is loaded with the name of the grand parent data directories observed with each data file.
	 * Make this something descriptive like, foundation, h1k, or exome, etc. Users can select 1 or more to restrict 
	 * what data is returned. */
	private TreeSet<String> secondaryDataTypesToReturn = new TreeSet<String>();
	boolean filterOnSecondaryDataTypes = false;
	
	/**This hash is loaded with the File data sources. */
	private TreeSet<File> dataFilesToReturn = new TreeSet<File>();
	boolean filterOnDataFiles = false;
	
	/**This hash contains the file types/ extensions found, just two so far: vcf and bed.*/
	private TreeSet<String> fileTypesToReturn = new TreeSet<String>();
	boolean filterOnFileTypes = false;
	
	private boolean matchVcf = false;
	private boolean fetchData = true;
	private boolean printJson = true;
	private boolean scoresSet = false;
	private int bpPadding = 0;
	
	//constructors
	public TQueryFilter(){
	}

	public String getCurrentFilters(String delimiter){
		StringBuilder sb = new StringBuilder();
		sb.append("fetchData="+ fetchData);
		sb.append(delimiter);
		sb.append("printJson="+ printJson);
		sb.append(delimiter);
		sb.append("matchVcf="+matchVcf);
		
		if (filterOnPrimaryDataTypes){
			sb.append(delimiter);
			sb.append("primaryDataTypes=");
			Iterator it = primaryDataTypesToReturn.iterator();
			sb.append(DataSources.buildStringFromIterator(it, ","));
		}
		if (filterOnSecondaryDataTypes){
			sb.append(delimiter);
			sb.append("secondaryDataTypes=");
			Iterator it = secondaryDataTypesToReturn.iterator();
			sb.append(DataSources.buildStringFromIterator(it, ","));
		}
		if (filterOnFileTypes){
			sb.append(delimiter);
			sb.append("fileTypes=");
			Iterator it = fileTypesToReturn.iterator();
			sb.append(DataSources.buildStringFromIterator(it, ","));
		}
		if (filterOnDataFiles){
			sb.append(delimiter);
			sb.append("dataSources=");
			Iterator it = dataFilesToReturn.iterator();
			sb.append(DataSources.buildStringFromIterator(it, ","));
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
	public boolean isMatchVcf() {
		return matchVcf;
	}
	public void setMatchVcf(boolean matchVcf) {
		this.matchVcf = matchVcf;
	}

	public boolean isFetchData() {
		return fetchData;
	}

	public void setFetchData(boolean fetchData) {
		this.fetchData = fetchData;
	}

	public boolean isPrintJson() {
		return printJson;
	}

	public void setPrintJson(boolean printJson) {
		this.printJson = printJson;
	}

	public boolean isScoresSet() {
		return scoresSet;
	}

	public void setScoresSet(boolean scoresSet) {
		this.scoresSet = scoresSet;
	}

	public int getBpPadding() {
		return bpPadding;
	}

	public void setBpPadding(int bpPadding) {
		this.bpPadding = bpPadding;
	}
	
}
