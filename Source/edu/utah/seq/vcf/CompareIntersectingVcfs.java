package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;

/**Compares variant lists showing which match by chr pos ref alt
 * @author Nix
 * */
public class CompareIntersectingVcfs {

	//user fields
	private File[] vcfFiles;
	private File results;
	
	//internal fields
	private TreeMap<String, ArrayList<String[]>> varRec = null;

	//constructor
	public CompareIntersectingVcfs(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		try {
			//process args
			processArgs(args);
			
			walkFiles();
			
			printSpreadsheet();
			
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing a vcf file\n");
		}


		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void printSpreadsheet() throws FileNotFoundException, IOException {
		Gzipper out = new Gzipper(results);
		
		//print header
		out.print("CHROM\tPOS\tREF\tALT\tNumInts");
		String[] fileNames = IO.fetchFileNames(vcfFiles);
		fileNames = Misc.trimCommonEnd(fileNames);
		for (String fn: fileNames) {
			out.print("\t");
			out.print(fn);
		}
		out.println();
		
		//for each var record
		for (String key: varRec.keySet()){
			ArrayList<String[]> hits = varRec.get(key);
			
			//print chr pos ref alt 
			out.print(key);
			out.print("\t");
			out.print(hits.size());
			
			//create String[] to represent each file
			String[] cells = new String[vcfFiles.length];
			
			//for each hit load the appropriate index in the String[]
			
			for (String[] hit: hits){
				int index = Integer.parseInt(hit[0]);
				StringBuilder sb = new StringBuilder();
				sb.append(hit[1]); sb.append(" ");
				sb.append(hit[2]); sb.append(" ");
				sb.append(hit[3]); 
				cells[index] = sb.toString();
			}
			
			//print the String[]
			for (String c: cells){
				out.print("\t");
				if (c != null) out.print(c);
				else out.print(".");
			}
			
			//close the line
			out.println();
		}
		out.close();
	}

	public void walkFiles() throws IOException{
		varRec = new TreeMap<String, ArrayList<String[]>>();
		
		String[] fields;
		String line;
		StringBuilder sb;
		for (int i=0; i< vcfFiles.length; i++){
			
			//walk and load records
			BufferedReader in = IO.fetchBufferedReader(vcfFiles[i]);
			while((line = in.readLine()) != null){
				line = line.trim();
				if (line.startsWith("#") || line.length() ==0) continue;
				fields = Misc.TAB.split(line);
				
				//parse key
				sb = new StringBuilder();
				sb.append(fields[0]); sb.append("\t");//chrom
				sb.append(fields[1]); sb.append("\t");//pos
				sb.append(fields[3]); sb.append("\t");//ref
				sb.append(fields[4]); //alt
				String key = sb.toString();

				//fetch or create AL
				ArrayList<String[]> al = varRec.get(key);
				if (al == null){
					al = new ArrayList<String[]>();
					varRec.put(key, al);
				}
				
				//add record, qual filter info
				String[] rec = {new Integer(i).toString(), fields[5], fields[6], fields[7]}; 
				al.add(rec);
			}
			in.close();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CompareIntersectingVcfs(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'r': results = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a directory containing vcf files to compare.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Compare Intersecting Vcfs : June 2019                     **\n" +
				"**************************************************************************************\n" +
				"Compares vcf files by creating a master list of variants and then scores each for the\n"+
				"presense of the same CHROM POS ALT REF in each vcf file.\n\n" +

				"Options:\n"+
				"-v A directory of vcf files to compare (xxx.vcf(.gz/.zip OK)).\n"+
				"-r Name of a spreadsheed results file, should end in xxx.txt.gz\n"+
				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/CompareIntersectingVcfs -v VCFs/\n" +
				"       -r comparisonVcf.txt.gz\n\n"+

		"**************************************************************************************\n");

	}
}
