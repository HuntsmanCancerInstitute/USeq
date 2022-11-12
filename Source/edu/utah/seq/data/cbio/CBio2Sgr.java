package edu.utah.seq.data.cbio;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class CBio2Sgr {

	//fields
	private File[] tsvFiles = null;

	public CBio2Sgr(String[] args){
		processArgs(args);

		for (File f: tsvFiles) parseTsv(f);
	}


	private void parseTsv(File f) {
		IO.pl("DatasetName\tNumSkippedRecords\tNumUniqueRecords\tNumTotalParsedRecords");
		IO.p(f.getName());
		String name = Misc.removeExtension(f.getName());
		File svg = new File (f.getParentFile(), name+".sgr");
		File bed = new File (f.getParentFile(), name+".bed");
		try {
			PrintWriter out = new PrintWriter( new FileWriter(svg));
			PrintWriter bedOut = new PrintWriter( new FileWriter(bed));
			CBioVar[] vars = fetchVars(f);
			for (CBioVar v: vars) {
				out.println(v.getSvg());
				bedOut.println(v.getBed());
			}
			out.close();
			bedOut.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printExit("\nFailed to parse: "+f);
		}

	}


	private CBioVar[] fetchVars(File cf) throws IOException {
		BufferedReader in = IO.fetchBufferedReader(cf);

		//find and parse the Sample ID line
		String line = null;
		boolean found = false;
		while ((line = in.readLine())!=null) {
			if (line.contains("Sample ID")) {
				found = true;
				break;
			}
		}
		if (found == false) throw new IOException("\nERROR: failed to find the header line containing 'Sample ID' in "+cf);
		String[] fields = Misc.TAB.split(line);
		HashMap<String, Integer> keyIndex = new HashMap<String, Integer>();
		for (int i=0; i< fields.length; i++) keyIndex.put(fields[i], i);
		Integer sampleIdIndex = keyIndex.get("Sample ID");
		Integer chromosomeIndex = keyIndex.get("Chromosome");
		Integer startPosIndex = keyIndex.get("Start Pos");
		Integer endPosIndex = keyIndex.get("End Pos");
		Integer refIndex = keyIndex.get("Ref");
		Integer varIndex = keyIndex.get("Var");
		
		if (sampleIdIndex==null || chromosomeIndex==null || startPosIndex==null || endPosIndex==null || refIndex==null || varIndex==null) {
			throw new IOException("\nFAILED to find one or more of the column indexes: \n'Sample ID', 'Chromosome', 'Start Pos', 'End Pos', 'Ref', 'Var' in:\n"+line);
		}

		//parse the rest
		int numTotal = 0;
		int unique = 0;
		int numSkipped = 0;
		ArrayList<CBioVar> vars = new ArrayList<CBioVar>();
		HashSet<String> keys = new HashSet<String>();
		try {
			while ((line = in.readLine())!=null) {
				numTotal++;
				String[] f = Misc.TAB.split(line);
				//no chrom or start pos
				if (f[chromosomeIndex].contains("NA") || f[startPosIndex].contains("-1") || f[endPosIndex].length()==0) {
					numSkipped++;
					continue;
				}
				CBioVar cbv = new CBioVar(f[sampleIdIndex], f[chromosomeIndex], Integer.parseInt(f[startPosIndex])-1, Integer.parseInt(f[endPosIndex]), f[refIndex], f[varIndex]);
				String key = cbv.fetchKey();
				if (keys.contains(key) == false) {
					unique++;
					vars.add(cbv);
					keys.add(key);
				}
			}
			in.close();
		} catch ( Exception e) {
			e.printStackTrace();
			IO.el("\nPROBLEM LINE: "+line+"\n");
			throw new IOException (e.getMessage());
		}
		IO.pl("\t"+ numSkipped+"\t"+unique+"\t"+ numTotal);

		CBioVar[] finalVars = new CBioVar[vars.size()];
		vars.toArray(finalVars);
		Arrays.sort(finalVars);
		return finalVars;
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CBio2Sgr(args);
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
					case 'd': tsvFiles = IO.extractFiles(args[++i], ".tsv"); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//genome version?
		if (tsvFiles == null || tsvFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any tsvFiles in your -d directory.");




	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                CBio 2 Sgr : July 2011                             **\n" +
				"**************************************************************************************\n" +
				"Converts the cBioPortal mutation table download tsv file to a center position svg file\n" +
				"for subsequent mutation density plotting. Removes duplicate sample variants.\n"+

				"Options:\n"+
				"-d Directory containing xxx.tsv files for parsing.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/CBio2Sgr -d GeneTableFiles \n\n"+

				"**************************************************************************************\n");

	}		


}
