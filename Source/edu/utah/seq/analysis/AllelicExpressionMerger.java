package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.*;

/** 
 * @author Nix    
 * */
public class AllelicExpressionMerger {

	//user defined fields
	private File first;
	private File second;
	private File composite;
	
	//constructor
	/**Stand alone.*/
	public AllelicExpressionMerger(String[] args){

		//set fields
		processArgs(args);
		
		TreeMap <String, String> firstData = load(first);
		loadSecond(second, firstData);
		
		
		//write em out
		try {
			Gzipper out = new Gzipper(composite);
			out.println("gene\tsnp.id\talt.dp\tref.dp");
			for (String key: firstData.keySet()) out.println(firstData.get(key));
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private TreeMap<String, String> load(File f) {
		TreeMap<String, String> hash = new TreeMap<String, String>();
		try {
			BufferedReader in = IO.fetchBufferedReader(f);
			String line;
			String[] t;
			while ((line = in.readLine()) != null){
				if (line.startsWith("gene")) continue;
				hash.put(fetchKey(line), line);
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return hash;
	}
	
	private void loadSecond(File f, TreeMap<String, String> map) {
		try {
			BufferedReader in = IO.fetchBufferedReader(f);
			String line;
			int numRecords = 0;
			int numShared = 0;
			while ((line = in.readLine()) != null){
				if (line.startsWith("gene")) continue;
				String key = fetchKey(line);
				numRecords++;
				if (map.containsKey(key)){
					System.out.println(key);
					System.out.println("\t"+map.get(key));
					System.out.println("\t"+ line);
					numShared++;
				}
				else map.put(key, line);
			}
			in.close();
			System.out.println("\n"+map.size()+"\tIn 1st\t"+first.getName());
			System.out.println(numRecords+"\tIn 2nd\t"+second.getName());
			System.out.println(numShared+"\tShared\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public String fetchKey(String line){
		//gene	snp.id	alt.dp	ref.dp
		//ENSG00000117984_CTSD	11_1774666_T_C_0	7	17
		String[] t = Misc.TAB.split(line);
		String[] coor = Misc.UNDERSCORE.split(t[1]);
		StringBuilder sb = new StringBuilder(t[0]);
		for (int i=0; i< 4; i++){
			sb.append("\t");
			sb.append(coor[i]);
		}
		return sb.toString();
	}


	public static void main(String[] args) {

		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AllelicExpressionMerger(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': first = new File(args[++i]); break;
					case 's': second = new File(args[++i]); break;
					case 'm': composite = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (first == null || first.canRead() == false) Misc.printErrAndExit("\nError: can't find your first GeneiASE file?\n");
		if (second == null || second.canRead() == false) Misc.printErrAndExit("\nError: can't find your second GeneiASE file?\n");
		if (composite == null ) Misc.printErrAndExit("\nError: can't find your tabixed bed file?\n");
		

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Allelic Expression Merger :  Sept 2016                  **\n" +
				"**************************************************************************************\n" +
				"App to merge two GeneiASE tables from the AlleleicExpressionDetector and the \n"+
				"AllelicExpressionRNASeqWriter. Where geneName coor duplicates are found, writes out\n"+
				"the first's record to the merged file.\n\n"+

				"Required Arguments:\n"+
				"-f First GeneiASE table (gene snp.id alt.dp ref.dp)\n"+
				"-s Second GeneiASE table (gene snp.id alt.dp ref.dp)\n"+
				"-m Merged output table results.\n"+
				
				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/AllelicExpressionMerger -f snpTable.txt\n"+
				"-s rnaSeqTable.txt -m mergedGeneiASETable.txt\n\n" +

				"**************************************************************************************\n");

	}
}
