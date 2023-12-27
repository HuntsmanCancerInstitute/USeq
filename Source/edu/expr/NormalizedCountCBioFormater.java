package edu.expr;

import java.io.*;
import java.util.regex.*;

import edu.utah.hci.misc.Gzipper;

import java.util.*;
import util.gen.IO;
import util.gen.Misc;

public class NormalizedCountCBioFormater {
	
	private File ensSymIdFile;
	private HashMap<String,String> ensSymIds;
	private File normCountFile;
	private String dupKeyValDelim = ":";
	private Pattern dupDelim = Pattern.compile(dupKeyValDelim);
	
	public NormalizedCountCBioFormater (String[] args){
		processArgs(args);
		
		//load key value file skipping duplicate keys and making duplicate values unique
		IO.pl("Loading Ensembl 2 gene symbol lookup table...");
		ensSymIds = IO.loadFileIntoHashMapUniqueValues(ensSymIdFile, dupKeyValDelim);
		IO.pl("Number unique ensembl gene IDs "+ensSymIds.size());
		
		//walk the norm count file
		IO.pl("\nParsing "+normCountFile.getName());
		parseAndSaveMatches();
		
		IO.pl("\nDone!\n");
	}
	
	
	//here, test and run it!
	
	public void parseAndSaveMatches(){
			try {
				HashSet<String> geneSymbols = new HashSet<String>();
				File results = new File (normCountFile.getParentFile(), Misc.removeExtension(normCountFile.getName())+"_ForCBio.txt.gz");
				Gzipper out = new Gzipper (results);
				BufferedReader in = IO.fetchBufferedReader(normCountFile);
				
				//read in the header, this is just a dump of the column names
				String line = in.readLine();
				String[] tokens = Misc.TAB.split(line);
				int numSamples = tokens.length;
				IO.pl("\t"+numSamples+"\t#Samples");
				out.println("Hugo_Symbol\t"+line);
				
				//counters
				int totalLines = 0;
				int numMatches = 0;
				int numNoMatches = 0;
				int numBadColumns = 0;
				int numDuplicateSymbolsSkipped = 0;
				
				while ((line = in.readLine())!=null){
					totalLines++;
					tokens = Misc.TAB.split(line);
					if (tokens.length-1 != numSamples) {
						numBadColumns++;
						IO.el("\tSkipping, incorrect # columns -> "+line);
						continue;
					}
					
					//look for the ensembl id in the lookup hash
					String ensId = tokens[0];
					String geneSymbol = ensSymIds.get(ensId);
					if (geneSymbol!= null) {
						
						//does it contain a value delimiter? If so check it
						if (geneSymbol.contains(dupKeyValDelim)) {
							String[] symbolKey = dupDelim.split(geneSymbol);
							//already seen?
							if (geneSymbols.contains(symbolKey[0])) {
								IO.el("\tSkipping, duplicate symbol -> "+geneSymbol+ " line # "+totalLines);
								numDuplicateSymbolsSkipped++;
								continue;
							}
							else {
								geneSymbols.add(symbolKey[0]);
								geneSymbol = symbolKey[0];
							}
						}
						
						//save the line
						numMatches++;
						out.print(geneSymbol);
						for (int i=1; i< tokens.length; i++) {
							out.print("\t");
							out.print(tokens[i]);
						}
						out.println();
					}
					else numNoMatches++;
				}
				IO.pl("\t"+totalLines+"\t# Data lines");
				IO.pl("\t"+numMatches+"\t# Matches written to file");
				IO.pl("\t"+ numNoMatches+"\t# No Matches");
				IO.pl("\t"+ numBadColumns+"\t# Incorrect sample number column lines");
				IO.pl("\t"+ numDuplicateSymbolsSkipped+"\t# Duplicate symbols");
				
				//close the reader and writer IO
				in.close();
				out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'n': normCountFile = new File(args[i+1]); i++; break;
					case 'e': ensSymIdFile = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					e.printStackTrace();
				}
			}
		}
		//check agilentFiles
		if (ensSymIdFile == null ) Misc.printExit("\nCannot find your tab delimited ensembl geneSymbol ID lookup file?\n");
		if (normCountFile == null) Misc.printExit("\nCannot find your DESeq2 normalized count file?\n");
	}	
	
	public static void printDocs(){ 
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                          NormalizedCountCBioFormater:  August 2023               **\n" +
				"**************************************************************************************\n" +
				"Parses a DESeq2 or RSEM output file containing normalized counts, converts the ensembl\n"+
				"gene names to HUGO symbols, and saves the output ready for import into cBioPortal.\n"+
				
				"\nRequired Options:\n"+
				"-e Path a tab delimited txt file containing Ensembl IDs and Gene Symbols.\n" +
				"-n Path to a normalized count file from DESeq2, e.g. \n"+
				"   write.table(assay(vst(dds, blind=FALSE)), file = 'vst.txt',quote=FALSE, sep ='\\t')'\n" +
				
				"\nExample: java -jar pathTo/Apps/NormalizedCountCBioFormater -e ensembIDsGeneSym.txt\n"+
				"   -n vst.txt\n" +	
		"**************************************************************************************\n");		
	}
	
	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new NormalizedCountCBioFormater(args);
	}
	

}
