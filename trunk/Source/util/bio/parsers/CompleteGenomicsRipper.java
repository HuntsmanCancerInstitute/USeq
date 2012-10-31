package util.bio.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class CompleteGenomicsRipper {

	//fields
	private File cgFile;

	//constructor
	public CompleteGenomicsRipper (String[] args){
		//make file object
		cgFile = new File (args[0]);

		//parse
		parseFile();
	}


	private void parseFile() {
		try {
			BufferedReader in = IO.fetchBufferedReader(cgFile);
			HashMap<String, PrintWriter> varTypesPrintWriters = new HashMap<String, PrintWriter>();
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			int printCounter = 0;
			int skipCounter =0;
			while ((line = in.readLine()) !=null){
				line = line.trim();
				if (line.length() == 0 || line.startsWith("#")) {
					System.out.println("Skipped "+line);
					skipCounter++;
					continue;
				}
				tokens = tab.split(line);
				//define the type
				String varType = tokens[6];
				//fetch PrintWriter
				PrintWriter out = null;
				//does it exist
				if (varTypesPrintWriters.containsKey(varType)) out = varTypesPrintWriters.get(varType);
				else{
					//no make new PrintWriter
					String name = cgFile.getName();
					name = Misc.removeExtension(name);
					name = name + "_" + varType +".bed";
					File f = new File (cgFile.getParentFile(), name);
					out = new PrintWriter (new FileWriter (f));
					varTypesPrintWriters.put(varType, out);
				}
				
				//print line to out chr, start, stop, name
				printCounter++;
				out.println(tokens[2]+"\t"+tokens[3]+"\t"+tokens[4]+"\t"+tokens[5]+"_"+tokens[6]);
			}

			//close writers
			for (PrintWriter pw: varTypesPrintWriters.values()) {
				System.out.println("\nClosing "+pw);
				pw.close();
			}

			//print varTypes
			System.out.println("Found the following varTypes:");
			System.out.println(varTypesPrintWriters.keySet());
			
			System.out.println("PC "+printCounter+"\tSC "+skipCounter);
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}



	}


	// main method called with this class
	public static void main (String[] arg){
		new CompleteGenomicsRipper (arg);
	}

}
