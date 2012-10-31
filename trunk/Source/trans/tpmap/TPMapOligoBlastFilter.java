package trans.tpmap;

import java.util.regex.*;
import java.io.*;
import java.util.*;

import util.gen.*;

/**
 * Application for remapping and filtering a text version bpmap file with BLAST.
 * 
 * This is wickedly slow :( !  Can't use -B or megablast so stuck with this junk. If you don't care about
 * mismatches then use the MummerFilter applicaiton.
 * 
 * Filters a tpmap file only printing out oligo lines that have one unique match to the sequenced genome.
 * The orientation direction is replaced with the number of 1bp matches for files appended with Bls or 
 * total number of matches (1bp + exact) for files appended with NoTrmBls.
 * 
 * To launch on Sapo cluster: java -Xmx1500M nix/util/gen/JQSub cd /home/sapo/nix/NixSourceCode/ return java -Xmx1500M nix/affy/TPMapOligoBlastFilter -d /home/sapo/nix/Affy/BLAST/whole_genome.fasta -s /home/sapo/software/seqanal/alignment/ncbi/bin/blastall -b /home/sapo/nix/Affy/TPMap/Extra/extra.tpmap
 */
public class TPMapOligoBlastFilter {
	//fields
	private File blastAll;
	private File tpmapFile; 
	private File database;
	
	public TPMapOligoBlastFilter(String[] args){
		processArgs(args);
		
		try {
			long start = System.currentTimeMillis();
			String line;
			String[] tokens;
			int counter = 0;
			File tempFastaFile = new File (tpmapFile.getCanonicalPath()+"tmp");
			BufferedReader in = new BufferedReader(new FileReader(tpmapFile));
			PrintWriter outNewBpmap = new PrintWriter(new FileWriter(tpmapFile.getCanonicalPath()+"Bls"));
			PrintWriter outNewBpmapNoTrim = new PrintWriter(new FileWriter(tpmapFile.getCanonicalPath()+"NoTrmBls"));
			PrintWriter outNoMatch = new PrintWriter(new FileWriter(tpmapFile.getCanonicalPath()+"NoMatch"));
			BlastFilterResult blastResults;
			int[] histo = new int[1001];
			int numRepeats = 0;
			String blastAllString = IO.getFullPathName(blastAll);
			String databaseString = IO.getFullPathName(database);
			HashSet exactMatchRepeats = new HashSet(1000);
			HashSet noExactMatches = new HashSet(1000);
			
			//run thru tpmap, blasting each oligo, if any repeats found set ori to be the number
			while ((line = in.readLine()) !=null) {           
				tokens = line.split("\\s+");
				PrintWriter out = new PrintWriter(new FileWriter(tempFastaFile));
				out.println(">x\n"+tokens[0]);
				out.close();
				//[numMatches, exactMatches, oneBPMismatches, startExactMatch]
				blastResults = blastForOneBPRepeats(blastAllString, databaseString, tempFastaFile);
				
				if (blastResults.getExactMatches()==1){
					//reset start position, chromosome, and put number of 1 bp mismatches in orientation
					tokens[3] = ""+blastResults.getStartExactMatch();
					tokens[1] = blastResults.getOneBPMisMatches()+"";
					tokens[2] = "chr" + blastResults.getChromosome();
					line = Misc.stringArrayToString(tokens,"\t");
					outNewBpmap.println(line);	
					outNewBpmapNoTrim.println(line);
				}
				
				//no exact match
				else if (blastResults.getExactMatches() == 0) {
					noExactMatches.add(tokens[0]);
					tokens[1] = Integer.toString(blastResults.getOneBPMisMatches());
					line = Misc.stringArrayToString(tokens,"\t");
					outNoMatch.println(line);
				}
				
				//write to file lines with 1 or more exact matches, put number of matches (1bp and exact) in orientation column
				else {
					//reset start position, chromosome, and put number of 1 bp mismatches in orientation
					tokens[3] = ""+blastResults.getStartExactMatch();
					tokens[1] = blastResults.getMatches()+"";
					tokens[2] = "chr" + blastResults.getChromosome();
					line = Misc.stringArrayToString(tokens,"\t");
					outNewBpmapNoTrim.println(line);
					//record number of exact match repeats
					exactMatchRepeats.add(tokens[0]);
				}
				
				//increment histogram counter
				numRepeats = blastResults.getMatches();
				if (numRepeats>1000) numRepeats = 1000;
				histo[numRepeats]++;
				counter++;
			}
			
			//close streams
			in.close();
			outNewBpmap.close();
			outNewBpmapNoTrim.close();
			outNoMatch.close();
			
			//print summary results
			StringBuffer sb = new StringBuffer();
			sb.append("\tNumber of tpmap lines BLASTed "+(counter)+"\n");
			sb.append("\tNumber of oligos removed that had more than one exact match "+exactMatchRepeats.size()+"\n");
			sb.append("\tNumber of oligos removed that had no exact match "+noExactMatches.size()+"\n");
			//sb.append("Histogram of matches (exact and or 1bp mismatches) where 0 is no match, 1 is a single match, 2... matches:\n");
			System.out.println(sb);
			//Histogram.printHistogram(histo);
			int elapse = (int)(System.currentTimeMillis()-start)/1000;
			String fin = "Finished "+elapse+" seconds";
			System.out.println(fin);
			
			//save file
			PrintWriter outFile = new PrintWriter(new FileWriter(new File(IO.getFullPathName(tpmapFile)+"Stats.txt")));
			outFile.println(sb);
			Histogram.printHistogram(histo,outFile);
			outFile.println(fin);
			outFile.close();
			
			//remove temp fasta file
			tempFastaFile.delete();
			
		} catch (Exception e) {e.printStackTrace();}
		
	}
	
	/**This has been tweaked to find 1bp mismatches in a 25mer to a genome.  No indels.
	 * Returns the number of matches found in the database, thus a unique oligo would be found once.
	 * Returns -1 if something went wrong.
	 * 
	 * Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start,q. stop, s. start, s. stop, e-value, bit score
	 *     0         1            2             3              4             5           6       7        8        9      10       11
	 * */
	public static BlastFilterResult blastForOneBPRepeats(String blastAll, String database, File fastaFile){
		int numMatches = 0;
		int exactMatches = 0;
		int oneBPMismatches = 0;
		int startExactMatch = 0;
		String chromosome = "";
		try {
			//make command line
			///home/sapo/software/seqanal/alignment/ncbi/bin/blastall -p blastn -i /home/sapo/nix/Affy/TPMapFiles/tempDeleteMe -d dmel3.fasta 
			//		-W 7 -F F -e .002  -q -1 -g F -m 8
			String[] commandArray ={blastAll,"-p", "blastn", "-i",fastaFile.getCanonicalPath(),"-d", database,
					"-W","12","-F","F","-e","0.002", "-q","-1","-g","F","-m","8"};
			Runtime rt = Runtime.getRuntime();
			//rt.traceInstructions(true); //for debugging
			//rt.traceMethodCalls(true); //for debugging
			Process p = rt.exec(commandArray);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			//BufferedReader data = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			String line;
			String[] tokens;
			int num;
			int length;
			while ((line = data.readLine()) != null){
				//System.out.println(line);
				tokens = line.split("\\s+");
				num = Integer.parseInt(tokens[4]);
				//is it a one base pair mismatch? watch ends
				length = 1+ Integer.parseInt(tokens[7]) - Integer.parseInt(tokens[6]);
				if (num !=0 || (length != 25)) oneBPMismatches++;
				//no its by default an exact match
				else {
					exactMatches++;
					startExactMatch = Integer.parseInt(tokens[8]);
					chromosome = tokens[1];
				}
				numMatches++;
			}
			data.close();
			//attempting to avoid the dreaded java.io.IOException: java.io.IOException: Too many open files
			//	due to firing thousands of jobs at once
			p.waitFor();
			p.getInputStream().close();
			p.getErrorStream().close();
			p.getOutputStream().close();
			p=null;
			data = null;
			rt = null;
			System.runFinalization();
			System.gc();
		}catch (Exception e){
			e.printStackTrace();
			numMatches = -1;
		}
		return new BlastFilterResult (chromosome, numMatches, exactMatches, oneBPMismatches, startExactMatch);
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		String[] celDirs = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': tpmapFile = new File(args[i+1]); i++; break;
					case 'd': database = new File(args[i+1]); i++; break;
					case 's': blastAll = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check to see if they entered required params
		if (tpmapFile==null || tpmapFile.canRead() == false){
			System.out.println("\nCannot find your tpmap file!\n");
			System.exit(0);
		}
		if (database==null){
			System.out.println("\nCannot find your database file!\n");
			System.exit(0);
		}
		if (blastAll==null || blastAll.canRead() == false){
			System.out.println("\nCannot find your blastall program!\n");
			System.exit(0);
		}
	}
	public static void printDocs(){
		System.out.println("\n"+ 
				"**************************************************************************************\n" +
				"**                         TPMapOligoBlastFilter: Feb 2005                          **\n" +
				"**************************************************************************************\n" +
				"TOBF takes an Affymetrix text bpmap file and BLASTs each oligo against a BLAST\n" +
				"database.  Use formatdb on a complete genome.fasta file to create the database.\n" +
				"BOBF writes three files, one that assigns new coordinates to the oligo and replaces the\n" +
				"ori column with the number of 1bp mismatches and exact matches.  A second file is\n" +
				"written containing oligos with only one exact match in the genome, the ori column is\n" +
				"replaced with the number of 1bp mismatches.  The third contains oligos with no exact\n" +
				"match to the genome.  The number of 1bp mismatches is assigned to the ori column.\n" +
				"This combo program is quite slow, 0.45sec\n" +
				"per oligo, so break up your tpmap file into separate files using the FileSplitter\n" +
				"program and farm the jobs out to a cluster. After processing, combine the files using\n" +
				"FileJoiner and sort using the TPMapSort program. If you don't care about 1bp mismatches\n" +
				"use MUMMERMapper, it should take about an hour. \n\n"+
				
				"Required Parameters:\n"+
				"-b Full path file text for text bpmap file.\n" +
				"-d Full path file text for BLAST database.\n" +
				"-s Full path file text for the blastall program.\n" +
				"\n" +
				"Example: java -Xmx256M -jar pathTo/T2/Apps/TPMapOligoBlastFilter -b /affy/tpmap1 -s\n" +
				"      /ncbi/bin/blastall -d /seq/dmel/whole_genome.fasta\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length !=6){
			printDocs();
			System.exit(0);
		}
		System.out.println("Launching TPMapOligoBlastFilter...");
		new TPMapOligoBlastFilter(args);
	}
	
	
}
