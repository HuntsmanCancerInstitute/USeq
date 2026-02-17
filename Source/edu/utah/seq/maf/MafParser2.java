package edu.utah.seq.maf;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;
import util.gen.*;


/**
 * Parses variant maf files.  
 * @author davidnix*/
public class MafParser2 {

	//user defined fields
	private File mafFile;
	private HashMap<String, Integer> headerNameIndex = new HashMap<String, Integer>();


	public MafParser2(File mafFile) throws IOException{
		this.mafFile = mafFile;
		loadHeader();
	}

	public TreeSet<String> fetchUniqueColumnValues(String columnName) throws IOException{
		TreeSet<String> uni = new TreeSet<String>();
		
		BufferedReader in = IO.fetchBufferedReader(mafFile);

		//pull column index
		Integer nameIndex = headerNameIndex.get(columnName);
		if (nameIndex == null) throw new IOException("ERROR: failed to find '"+columnName+"' in the parsed header for "+mafFile);
		int index = nameIndex.intValue();
		
		//for each line in the file
		String line;
		while ((line = in.readLine()) != null){
			//skip blanks
			if (line.trim().length() == 0 || line.startsWith("#") || line.startsWith("Hugo_Symbol")) continue;
			String[] fields = Misc.TAB.split(line);
			String value = fields[index].trim();
			if (value.length()>0) uni.add(value);
		}
		in.close();
		return uni;
	}

	private void loadHeader() throws IOException{
		BufferedReader in = IO.fetchBufferedReader(mafFile);

		//for each line in the file
		String line;
		while ((line = in.readLine()) != null){
			//skip blanks
			if (line.trim().length() == 0 || line.startsWith("#")) continue;
			//header
			if (line.startsWith("Hugo_Symbol")) {
				String[] hugo = Misc.TAB.split(line);
				for (int x=0; x< hugo.length; x++) headerNameIndex.put(hugo[x], x);
				in.close();
				return;
			}
		}
		in.close();
		throw new IOException("ERROR: failed to parse the header from "+mafFile);

	}
	
	public ArrayList<String[]> fetchGeneLines(TreeSet<String> hugoSymbols) throws IOException{
		BufferedReader in = IO.fetchBufferedReader(mafFile);
		ArrayList<String[]> lines = new ArrayList<String[]>();

		//pull column index
		Integer hugoGeneSymbolIndex = headerNameIndex.get("Hugo_Symbol");
		if (hugoGeneSymbolIndex == null) throw new IOException("ERROR: failed to find 'Hugo_Symbol' in the parsed header for "+mafFile);
		int index = hugoGeneSymbolIndex.intValue();
		
		//for each line in the file
		String line;
		while ((line = in.readLine()) != null){
			//skip blanks
			if (line.trim().length() == 0 || line.startsWith("#") || line.startsWith("Hugo_Symbol")) continue;
			String[] fields = Misc.TAB.split(line);
			String geneName = fields[index];
			if (hugoSymbols.contains(geneName)) lines.add(fields);
			
		}
		in.close();
		return lines;
	}



	public HashMap<String, Integer> getHeaderNameIndex() {
		return headerNameIndex;
	}

	public TreeSet<String> fetchUniqueColumnValues(ArrayList<String[]> mafLines, String columnName) throws IOException {
		TreeSet<String> uni = new TreeSet<String>();

		//pull column index
		Integer nameIndex = headerNameIndex.get(columnName);
		if (nameIndex == null) throw new IOException("ERROR: failed to find '"+columnName+"' in the parsed header for "+mafFile);
		int index = nameIndex.intValue();
		
		//for each provided line 
		for (String[] fields : mafLines){
			String value = fields[index].trim();
			if (value.length()>0) uni.add(value);
		}
		return uni;
	}

	public ArrayList<String[]> filterLines(ArrayList<String[]> mafLines, TreeSet<String> toFilterFor, String columnName) throws IOException {
		//pull column index
		Integer nameIndex = headerNameIndex.get(columnName);
		if (nameIndex == null) throw new IOException("ERROR: failed to find '"+columnName+"' in the parsed header for "+mafFile);
		int index = nameIndex.intValue();
		ArrayList<String[]> lines = new ArrayList<String[]>();
		//for each provided line 
		for (String[] fields : mafLines){
			String value = fields[index].trim();
			if (toFilterFor.contains(value)) lines.add(fields);
		}
		return lines;
	}




}
