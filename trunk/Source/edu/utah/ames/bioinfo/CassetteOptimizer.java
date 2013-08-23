package edu.utah.ames.bioinfo;

/**
 * @author darren.ames@hci.utah.edu
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import util.gen.FisherExact;
import util.gen.Misc;
import util.gen.Num;


public class CassetteOptimizer {
	
	//fields
	private static String inputFastq1;
	private static String inputFastq2;
	private HashMap<String,Cartridge> cartridgeMap = new HashMap<String,Cartridge>();
	//private static int fivePrimeConsSeqErrorRate;
	//private HashMap<String, Integer> aaSeq;
	//private HashMap<String, Integer> DNAStringCounts;
	private String readSequence;
	private String MATCHCASSETTE = "[ACGT]+TGACGCGTCT([ACGT]{18}TG[CT])".toUpperCase(); //matches cassette (18bp+TC[C/T] + last 10bp of adjacent 5' invariant region
	private static int minReadCountReportAA = 10;
	private double cutOff= 0.05; //p-value cutoff

	public static void main(String[] args) throws Exception {
		
		//check for args
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		CassetteOptimizer co = new CassetteOptimizer(args);
		co.readFiles();
	}
	
	//constructor
	public CassetteOptimizer(String[] args) {
		processArgs(args);
	}
	
	public void readFiles() throws Exception {
		
		String[] files = {inputFastq1,inputFastq2};
		for (int i = 0; i < files.length; i++) {
			this.parseFastqFile(files[i],i);
		}
		//output file name
		String outName = null;
		String dirName = null;
		
		//fetch input path
		Pattern p2 = Pattern.compile("/.+/");
		Matcher m2 = p2.matcher(inputFastq1.toString());
		if (m2.find()) {
			dirName = m2.group();
		}

		//strip path from input file
		int index = inputFastq1.lastIndexOf(File.separatorChar);
		String name = inputFastq1.substring(index+1);

		//strip extension from input file
		Pattern pat = Pattern.compile(".+(?=.fastq((\\.)+)*)");
		Matcher mat = pat.matcher(name);
		if (mat.find()) {
			outName = mat.group();
		}
		//good output		
		BufferedWriter seq = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.good.txt")));
		
		FisherExact fe = null;
		
		ArrayList<Cartridge> cartridgeList = new ArrayList<Cartridge>(cartridgeMap.values());
		
		for (Cartridge c : cartridgeList) {
			if (fe == null) {

				fe = new FisherExact((int)c.getMax());
			}
			//calculate p-values
			c.calculateFisher(fe);
		}
		//sort p-values in descending order
		Collections.sort(cartridgeList, new Cartridge.CartridgeComparator());
		
		double[] pValueList = new double[cartridgeList.size()];
		for (int i = 0; i < pValueList.length; i++) {
			pValueList[i] = cartridgeList.get(i).getpValue();
			
		}
		//calculate FDR values on p-values
		Num.benjaminiHochbergCorrect(pValueList);
		for (int i = 0; i < pValueList.length; i++) {
			cartridgeList.get(i).setFDR(pValueList[i]);
			//System.out.println(pValueList[i]);
		}
		//reverse order so most significant results printed first
		Collections.reverse(cartridgeList);
		for (Cartridge c : cartridgeList) {
			String val = c.getFisher(cutOff);
			if (val != null) {
				seq.write(val);
			}
		}
		seq.close();
	}
	
	/**
	 * Parses the input fastq file
	 * @throws Exception
	 */
	public void parseFastqFile(String fastq, int sampleNum) throws Exception {
		
		//output file name
		String outName = null;
		String dirName = null;
		
		//fetch input path
		Pattern p2 = Pattern.compile("/.+/");
		Matcher m2 = p2.matcher(fastq.toString());
		if (m2.find()) {
			dirName = m2.group();
		}
				
		//strip path from input file
		int index = fastq.lastIndexOf(File.separatorChar);
		String name = fastq.substring(index+1);
		
		//strip extension from input file
		Pattern pat = Pattern.compile(".+(?=.fastq((\\.)+)*)");
		Matcher mat = pat.matcher(name);
		if (mat.find()) {
			outName = mat.group();
		}
		
		//read in gzipped fastq file
		BufferedReader br = new BufferedReader(new InputStreamReader
				(new GZIPInputStream(new FileInputStream(fastq))));
		
		//make hashmap to hold amino acid sequences of target region
		Map<String, Integer> aaSeq = new HashMap<String, Integer>();
		
		String line;
		
		int i = 0; //total line counter
		int read = 0; //total read counter
		int j = 0; //reads containing cassette (region of interest)
		int short5prime = 0; //reads with truncated 5' region 
		
		//make a big sloppy pile of buffered writers to output various things
		//Reads long enough, but missing cassette
		BufferedWriter noCassette = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.noCassette.txt")));
		noCassette.write("Length-filtered reads missing cassette\n\n");
		
		//Reads too short to contain cassette
		BufferedWriter trunc5Prime = new BufferedWriter(new FileWriter(new File(dirName, outName + "optimizer.trunc5Prime.txt")));
		trunc5Prime.write("Reads with truncated 5' region, therefore missing cassette\n\n");
		
		//polypeptides failing minimum total count per unique string
		BufferedWriter failCount = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.failCount.txt")));
		failCount.write("Amino acid sequences failing minimum read count per motif\n\n");
		
		//output file for sequence lines only
		BufferedWriter seq = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.seq.txt"))); 
		
		//polypeptide strings containing internal stop codons
		BufferedWriter stop = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.internalStop.txt")));
		stop.write("Polypeptide strings containing an internal stop\n\n)");
		
		//polypeptides containing K, R, or C residues in any of positions 1-6 (position 7 not included)
		BufferedWriter withKRC = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.withKRC.txt")));
		withKRC.write("Polypeptide strings containing K, R, or C residues in any of positions 1-6\n\n");
		
		//polypeptides NOT containing K, R, or C residues in any of positions 1-6 (positions 7 is ok)
		BufferedWriter withoutKRC = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.withoutKRC.txt")));
		withoutKRC.write("Polypeptide strings NOT containing K, R, or C residues in any of positions 1-6\n\n");
		
		while ((line = br.readLine()) != null) {
			//get total line #
			i++;
			
			//take read sequence line
			if (i%4 == 2) {
				/**
				 * //too short?
				//111 comes from 5' conserved region (92-6[barcode]=86;86+4[wiggle room]=90)) 
				//21 (invariable length of region of interest)
				//90+21=111
				if (line.length() < 111) {
					tooShort.write(String.format("%s\n", line));
					continue;
				}
			 * 
			 */
				//increment read#
				read++;
				readSequence = line;
				
				//grab sequence lines shorter than length of entire 5' conserved region
				if (line.length() < 86) {
					//write read sequences to file
					seq.write(String.format("%s\n", line));
					trunc5Prime.write(String.format("%s\n", line));
					short5prime++;
				}
				//grab all else, truncating length to 86bp
				else {
					seq.write(String.format("%s\n", line.substring(0, 86)));
				}
				
				//check for cassette and portion of adjacent 5' invariant region
				Pattern p = Pattern.compile(MATCHCASSETTE);
				Matcher m = p.matcher(readSequence);
				if (m.find()) {
					j++;
					
					//check for proper reading frame in cassette
					if (m.group(1).length() % 3 == 0) {

						String aa = this.getAminoAcid(m.group(1));
						
						if (aa.length() > 7) { //check for stop codons (skip them)
							stop.write(aa + "\n");
							continue;
						}
						//add amino acid seqs to hashmap with count
						if (cartridgeMap.containsKey(aa)) {
							cartridgeMap.get(aa).increment(sampleNum);
						}
						//set count to 1 if amino acid string not contained in hashmap
						else {
							Cartridge car = new Cartridge(aa, sampleNum);
							cartridgeMap.put(aa, car);
						}
					}
				}
				//missing region of interest
				else {
					//write to out file
					noCassette.write(String.format("%s\n", readSequence));
				}
			}
		}
		//close buffered reader
		br.close();
		
		//make buffered writer for output file
		//BufferedWriter out = new BufferedWriter(new FileWriter(new File(dirName, outName + ".optimizer.txt")));
		
		//print outfile header
		//out.write("CassetteOptimizer run of " + name + "\n\n");
		
		int w = 0;	//counter for key:value pairs passing minReadCountReportAA
		int x = 0; 	//counter for total # unique key:value pairs
		//print HashMap key:value pairs in descending order of value
		List<Entry<String, Integer>> outHash = entriesSortedByValues(aaSeq);
		for (Entry<String, Integer> entry : outHash) {
			x++;
			// only print amino acid strings >= 10 (default)
			if (entry.getValue() >= minReadCountReportAA) {
				w++;
				System.out.println(String.format("%s\t%s", entry.getKey(), entry.getValue()));
				
				//write results to file
			//	out.write(String.format("%s\t%s\n", entry.getKey(), entry.getValue()));
			}
			//amino acid string fails minReadCountReportAA (default >= 10)
			else {
				//write to out file		
			//	failCount.write(String.format("%s\t%s\n", entry.getKey(), entry.getValue()));
			}
		}
		//write filtering stats to file
	//	out.write("\nTotal reads: " + read);
	//	out.write("\nLength-filtered reads containing region of interest: " + j);
	//	out.write("\nTotal unique amino acid sequences: " + x);
	//	out.write("\nAmino acid sequences passing minimum read count: " + w);
		
		trunc5Prime.write("\nReads with a truncated 5' region, therefore missing cassette: " + short5prime);
		//print filtering stats
		System.out.println("\nTotal reads: " + read);
		System.out.println("Length-filtered reads containing region of interest: " + j);
		System.out.println("Total unique amino acid sequences: " + x);
		System.out.println("Amino acid sequences passing minimum read count: " + w);
		
		//close output files
	//	out.close();
		stop.close();
		noCassette.close();
		failCount.close();
		seq.close();
		trunc5Prime.close();
	}
	
	/**
	 * Sort by value of hashmap
	 * @param map
	 * @return
	 */
	static <K,V extends Comparable<? super V>> List<Entry<K,V>> entriesSortedByValues(Map<K,V> map) {

		List<Entry<K,V>> sortedEntries = new ArrayList<Entry<K,V>>(map.entrySet());

		Collections.sort(sortedEntries, new Comparator<Entry<K,V>>() {
			@Override
			public int compare(Entry<K,V> e1, Entry<K,V> e2) {
				return e2.getValue().compareTo(e1.getValue());
			}
		}
				);

		return sortedEntries;
	}
	
	/**
	 * Gets translated amino acid sequence from nucleotide string
	 * @param nucleotide
	 * @return
	 */
	public String getAminoAcid(String nuc) {
		String aminoAcid = "";
		for (int i = 0; i < nuc.length(); i+=3) {
			String codon = nuc.substring(i, i+3);
			aminoAcid += CassetteOptimizer.aminoAcidTable.get(codon);
		}
		//System.out.println(aminoAcid);
		return aminoAcid;
	}
	
	/**
	 * Check if the given DNA sequence contains correct neucleotide codes (A, C, G, T)
	 * @param dna
	 * @return
	 */
	public boolean valid(String nuc) {
		
		char[] validCodes = {'T', 'C', 'A', 'G'};
		//Assume that DNA is a valid string
		boolean validOrNot = true; 
		
		for (int position = 0; position < nuc.length(); position++) {
			char nextChar = nuc.charAt(position);
			if (!(nextChar == validCodes[0] || nextChar == validCodes[1] || nextChar == validCodes[2] || nextChar == validCodes[3])) {
					validOrNot = false;
					break;
			}
		}
		return (validOrNot);
	}
	
	/**
	 * HashMap containing the amino acid codon table
	 */
	@SuppressWarnings("serial")
	public static HashMap<String, String> aminoAcidTable = new HashMap<String, String>(){{
		//put in codon translations
		//alanine
		put("GCT","A");
		put("GCC","A");
		put("GCA","A");
		put("GCG","A");
		//arginine
		put("CGT","R");
		put("CGC","R");
		put("CGA","R");
		put("CGG","R");
		put("AGA","R");
		put("AGG","R");
		//asparagine
		put("AAT","N");
		put("AAC","N");
		//aspartic acid
		put("GAT","D");
		put("GAC","D");
		//cysteine
		put("TGT","C");
		put("TGC","C");
		//glutamine
		put("CAA","Q");
		put("CAG","Q");
		//glutamic acid
		put("GAA","E");
		put("GAG","E");
		//glycine
		put("GGT","G");
		put("GGC","G");
		put("GGA","G");
		put("GGG","G");
		//histidine
		put("CAT","H");
		put("CAC","H");
		//isoleucine
		put("ATT","I");
		put("ATC","I");
		put("ATA","I");
		//leucine
		put("TTA","L");
		put("TTG","L");
		put("CTT","L");
		put("CTC","L");
		put("CTA","L");
		put("CTG","L");
		//lysine
		put("AAA","K");
		put("AAG","K");
		//methionine
		put("ATG","M");
		//phenylalanine
		put("TTT","F");
		put("TTC","F");
		//proline
		put("CCT","P");
		put("CCC","P");
		put("CCA","P");
		put("CCG","P");
		//serine
		put("TCT","S");
		put("TCC","S");
		put("TCA","S");
		put("TCG","S");
		put("AGT","S");
		put("AGC","S");
		//threonine
		put("ACT","T");
		put("ACC","T");
		put("ACA","T");
		put("ACG","T");
		//tryptophan
		put("TGG","W");
		//tyrosine
		put("TAT","Y");
		put("TAC","Y");
		//valine
		put("GTT","V");
		put("GTC","V");
		put("GTA","V");
		put("GTG","V");
		//start(M)
		put("ATG","M");
		//stop
		put("TAA","<stop-ochre>");
		put("TGA","<stop-opal>");
		put("TAG","<stop-amber>");
	}};

	/**
	 * This method will process each argument and assign new variables.
	 * @param args
	 */
	public void processArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		String programArgs = Misc.stringArrayToString(args, ",");
		boolean verbose = false;
		if (verbose) System.out.println("\nArguments: " + programArgs + "\n");
		for (int i = 0; i < args.length; i++) {
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()) {
				char test = args[i].charAt(1);
				try {
					switch (test) {
					case 'c': inputFastq1 = new String(args[++i]); break; //control file
					case 't': inputFastq2 = new String(args[++i]); break; //treatment file
					case 'n': minReadCountReportAA = new Integer(args[++i]); break;
					case 'p': cutOff = new Double(args[++i]); break;//p-value cutoff
					default: Misc.printErrAndExit("\nProblem--unknown option used!" + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
		//print error and exit if no input file specified
		if (inputFastq1 == null) Misc.printErrAndExit("\nPlease provide full path to fastq file(s) to parse.\n");
	}
	
	public static void printDocs() {
		System.out.println("\n" +
				"**********************************************************************************\n" +
				"**                         CassetteOptimizer: July 2013                         **\n" +
				"**********************************************************************************\n" + 
				"" +
				"\nParameters: \n\n" +
				" -c control fastq file (.gz ok)\n" +
				" -t treatment fastq file (.gz ok)\n" +
				" -n minimum number of reads to report translated amino acid string\n" +
				" -p p-value threshold (default is 0.05)\n" +
				" Default: 10 reads\n" +

				"Usage:\n\n" +
				"java -jar pathTo/CassetteOptimizer.jar -c A.fastq.gz -t B.fastq.gz -n 15\n\n" +
				"Questions or comments? Contact: darren.ames@hci.utah.edu\n" +
				"**********************************************************************************\n");
	}
}
