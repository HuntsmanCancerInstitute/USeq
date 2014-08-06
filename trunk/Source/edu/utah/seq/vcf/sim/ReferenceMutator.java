package edu.utah.seq.vcf.sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Given a table of snps (OMNI chip SNPMap.txt) and directory of reference fasta files.  Introduces the non reference base to the fasta.
 * @author Nix*/
public class ReferenceMutator {

	//fields
	private HashMap<String, File> chromFasta;
	private File saveDirectory;
	private File snpTableFile;
	private HashMap<String, Integer> snpType;
	private char[][] snpTypesSplit;
	private HashMap<String, ArrayList<SnpInfo>> chromSnps;
	private Gzipper bed;

	public ReferenceMutator (String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		makeSnpHash();

		//parse snpTable by chromosome
		System.out.println("Loading snp map table...");
		parseSnpTable();

		//for each chromosome of snps, load and modify fasta
		System.out.println("\nModifying fasta files...");
		modifyFastas();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void modifyFastas() {
		try {
			bed = new Gzipper(new File (saveDirectory, Misc.removeExtension(snpTableFile.getName()) + "_Ref2Alt.bed.gz"));
			for (String chrom : chromSnps.keySet()){
				//any fasta?
				File fasta = chromFasta.get(chrom);
				System.out.println(chrom);
				if (fasta == null) {
					System.out.println("\tSkipping, no fasta!");
				}
				else {
					MultiFastaParser mfp = new MultiFastaParser(fasta);
					Gzipper out = new Gzipper(new File (saveDirectory, chrom+".fasta.gz"));
					out.println("> "+mfp.getNames()[0]);
					String mutSeq = mutate(chrom, mfp.getSeqs()[0].toUpperCase(), chromSnps.get(chrom), out);
					out.println(mutSeq);
					out.close();
				}
			}
			bed.close();
			System.out.println();
		} catch (Exception e) {
			e.printStackTrace();
		} 
	}

	private String mutate(String chromName, String chromSeq, ArrayList<SnpInfo> snps, Gzipper out) throws IOException {
		char[] seq = chromSeq.toCharArray();
		
		int numMuts = snps.size();
		//int counter = 0;
		for (int i=0; i< numMuts; i++){
			SnpInfo si = snps.get(i);
			int position = si.position;
			int type = si.type;
			char refSeq = chromSeq.charAt(position);
			//char refSeq = seq[position];
			char[] pair = snpTypesSplit[type];
			//check base, must have an exact match for types 1 and 3
			if (type == 1 || type == 3){
				if (pair[0] == refSeq) seq[position] = pair[1];
				else if (pair[1] == refSeq) seq[position] = pair[0];
				else System.err.println("\tWARNING: did not find an exact match, skipping: position "+position+" type "+type + " ref "+refSeq+" pair "+pair[0]+" "+pair[1]+" "+si.name);
			}
			else {
				//check for exact match, if found flip to other
				if (pair[0] == refSeq) seq[position] = pair[1];
				else if (pair[1] == refSeq) seq[position] = pair[0];
				else {
					//check for opp base
					char compRefSeq = Seq.getComplement(refSeq);
					if (pair[0] == compRefSeq) seq[position] = Seq.getComplement(pair[1]);
					else if (pair[1] == compRefSeq) seq[position] = Seq.getComplement(pair[0]);
					else System.err.println("\tWARNING: did not find an appropriate match, skipping: position "+position+" type "+type + " ref "+refSeq+" pair "+pair[0]+" "+pair[1]+" "+si.name);
				}
			}
			//write out bed line of changes
			bed.println(chromName+"\t"+position+"\t"+(position+1)+"\t"+si.name+"_"+refSeq+"_"+seq[position]);

		}
		return new String(seq);
	}

	public void makeSnpHash(){
		String[] snpTypeF = {"A/G", "A/T", "A/C", "G/C", "G/T", "T/C"};
		String[] snpTypeR = {"G/A", "T/A", "C/A", "C/G", "T/G", "C/T"};
		snpTypesSplit = new char[][]{{'A','G'}, {'A','T'}, {'A','C'}, {'G','C'}, {'G','T'}, {'T','C'}};

		snpType = new HashMap<String, Integer>();
		for (int i=0; i< snpTypeF.length; i++) snpType.put(snpTypeF[i], new Integer(i));
		for (int i=0; i< snpTypeR.length; i++) snpType.put(snpTypeR[i], new Integer(i));
	}


	/*Creates a hash of chromosome: arraylist of position and type*/
	private void parseSnpTable() {
		chromSnps = new HashMap<String, ArrayList<SnpInfo>>();
		try {
			BufferedReader in = IO.fetchBufferedReader(snpTableFile);
			String[] tokens;
			String line = in.readLine();
			ArrayList<SnpInfo> al = null;
			//1483    kgp10007872     15      65325698        [A/C]   TOP     BOT
			while ((line = in.readLine()) != null){
				tokens = Misc.TAB.split(line);
				if (tokens.length != 7) Misc.printErrAndExit("\nProblem parsing snp line "+ line);
				if (tokens[2].equals("MT")) tokens[2] = "M";
				String chrom = "chr"+ tokens[2];
				al = chromSnps.get(chrom);
				if (al == null){
					al = new ArrayList<SnpInfo>();
					chromSnps.put(chrom, al);
				}
				String snp = tokens[4].substring(1, 4);
				Integer type = snpType.get(snp);
				if (type == null) System.err.println("\tWARNING: Failed to find a snp type for "+snp+" from "+line+". Skipping!");
				else {
					int position = Integer.parseInt(tokens[3]) - 1;  //they are using one based numbering!
					SnpInfo si = new SnpInfo(position, type, new String(tokens[1]));
					al.add(si);
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private class SnpInfo {
		int position;
		int type;
		String name;

		public SnpInfo (int position, int type, String name){
			this.position = position;
			this.type = type;
			this.name = name;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ReferenceMutator(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File fastaDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': saveDirectory = new File(args[++i]); break;
					case 'f': fastaDir = new File(args[++i]); break;
					case 't': snpTableFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//fetch fastas
		if (fastaDir == null || fastaDir.exists() == false) Misc.printErrAndExit("\nCan't find your fasta file directory? Aborting!\n");
		chromFasta = Seq.fetchChromosomeFastaFileHashMap(fastaDir);

		if (snpTableFile == null || snpTableFile.canRead() == false) Misc.printErrAndExit("\nCan't find table of snps? Aborting!\n");
		if (saveDirectory == null ) Misc.printErrAndExit("\nCan't find your save directory? Aborting!\n");
		saveDirectory.mkdirs();

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Reference Mutator  : Aug 2014                         **\n" +
				"**************************************************************************************\n" +
				"Takes a directory of fasta chromosome sequence files and converts the reference allele\n"+
				"to the alternate provided by a snp mapping table.\n\n" +

				"Required:\n"+
				"-f Full path to a directory containing chromosome specific fasta files. zip/gz OK.\n"+
				"-t Full path to a snp mapping table.\n"+
				"-s Full path to a directory to save the alternate fasta files.\n\n" +

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/ReferenceMutator -f /Hg19/Fastas\n"+
				"    -s /Hg19/AltFastas/ -t /Hg19/omni2.5SnpMap.txt\n\n" +

				"**************************************************************************************\n");

	}
}
