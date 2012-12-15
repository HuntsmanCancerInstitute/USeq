package util.bio.seq;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import trans.anno.GenomicRegion;
import util.bio.parsers.MultiFastaParser;
import util.gen.IO;
import util.gen.Misc;

/**
 * For manipulating nucleic acid sequences.
 */
public class Seq {
	
	public static final Pattern chrNumber = Pattern.compile("chr[123456789][1234567890]?");
	public static final Pattern chrLetter = Pattern.compile("chr[XYM]T?");
	public static final HashMap<String,Integer> asci2FastQScore = asci2FastQScore();
	public static final String[] ORDERED_ASCII = new String[]{"!", "\"", "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", 
		"6", "7", "8", "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
		"O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_", "`", "a", "b", "c", "d", "e", "f", "g", 
		"h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "{", "|", "}", "~"};	
	
	/**Returns a vertebrate codon:aa HashMap, all capitals, single letter aa's.*/
	public static HashMap<String,String> fetchCodonAAHashMap(){
		String[] x = {"GCT","A","GCC","A","GCA","A","GCG","A","CGT","R","CGC","R","CGA","R","CGG","R","AGA","R","AGG","R","AAT",
				"N","AAC","N","GAT","D","GAC","D","TGT","C","TGC","C","CAA","Q","CAG","Q","GAA","E","GAG","E","GGT","G","GGC","G",
				"GGA","G","GGG","G","CAT","H","CAC","H","ATT","I","ATC","I","ATA","I","TTA","L","TTG","L","CTT","L","CTC","L",
				"CTA","L","CTG","L","AAA","K","AAG","K","ATG","M","TTT","F","TTC","F","CCT","P","CCC","P","CCA","P","CCG","P",
				"TCT","S","TCC","S","TCA","S","TCG","S","AGT","S","AGC","S","ACT","T","ACC","T","ACA","T","ACG","T","TGG","W",
				"TAT","Y","TAC","Y","GTT","V","GTC","V","GTA","V","GTG","V","AAT","B","AAC","B","GAT","B","GAC","B","CAA","Z",
				"CAG","Z","GAA","Z","GAG","Z","TAA","*","TAG","*","TGA","*"};
		return Misc.createHashMap(x);
	}
	
	/**Returns a ambigBase:BASES DNA HashMap, all capitals.*/
	public static HashMap<String,String> fetchDNAAmbigBase2BasesHashMap(){
		String[] x = {"R","GA", "Y","TC", "S","GC", "W","TA", "K","GT", "M","AC", "D","GTA", "H","TAC", "B","GTC", "V","GAC", "N","GATC", "X","GATC"};
		return Misc.createHashMap(x);
	}
	
	/**Returns null if the ambiguous base isn't found or it doesn't code for two bases.
	 * Otherwise returns the non refBase.*/
	public static String findOppositeBase(String refBase, String ambigSymbol, HashMap<String,String> dnaAmbig2BasesHashMap){
		String bases = dnaAmbig2BasesHashMap.get(ambigSymbol);
		if (bases == null || bases.length() !=2) return null;
		return bases.replaceFirst(refBase, "");
	}
	
	/**Returns a map of chromosome text and file provided a directory of xxx.fa(.gz/.zip OK) or xxx.fasta(.gz/.zip OK)*/
	public static HashMap<String, File> fetchChromosomeFastaFileHashMap (File fastaDirectory){
		File[][] tot = new File[6][];
		tot[0] = IO.extractFiles(fastaDirectory,".fasta");
		tot[1] = IO.extractFiles(fastaDirectory,".fasta.gz");
		tot[2] = IO.extractFiles(fastaDirectory,".fasta.zip");
		tot[3] = IO.extractFiles(fastaDirectory,".fa");
		tot[4] = IO.extractFiles(fastaDirectory,".fa.gz");
		tot[5] = IO.extractFiles(fastaDirectory,".fa.zip");
		File[] f = IO.collapseFileArray(tot);
		HashMap<String,File> fastaFiles = new HashMap<String,File>();
		for (int i=0; i< f.length; i++){
			String chr = f[i].getName();
			chr = chr.replace(".gz", "");
			chr = chr.replace(".zip", "");
			chr = chr.replace(".fasta", "");
			chr = chr.replace(".fa", "");
			fastaFiles.put(chr, f[i]);
		}
		return fastaFiles;
	}
	
	/**Assumes dnaSequence is all upper case. See fetchCodonAAHashMap(). Returns single letter representations for aa. Does not
	 * recognize non GATC characters, will call such codons an X. Will not return last partial codon.*/
	public static String translate(String dnaSequence, HashMap<String,String> codon2AAHashMap){
		int size = dnaSequence.length();
		if (size < 3) return "";
		StringBuilder aaSeq = new StringBuilder();
		for (int i=0; i< size; i+=3){
			int stop = i+3;
			if (stop <= size) {
				String aa = codon2AAHashMap.get(dnaSequence.substring(i, stop));
				if (aa == null) aaSeq.append("X");
				else aaSeq.append(aa);
			}
			else break;
		}
		return aaSeq.toString();
	}
	
	/**Returns a map of asci text character to it's associated Sanger fastq base quality score.
	 * ! = 0*/
	public static HashMap<String,Integer> asci2FastQScore(){	
		String[] acii = new String[]{"!", "\"", "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", 
			"6", "7", "8", "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
			"O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_", "`", "a", "b", "c", "d", "e", "f", "g", 
			"h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "{", "|", "}", "~"};	
		HashMap<String,Integer> map = new HashMap<String,Integer>();
		for (int i=0; i< acii.length; i++) map.put(acii[i], new Integer(i));
		return map;
	}
	
	/**Converts the ascii quality scores to numeric scores. Assumes Sanger fastq.*/
	public static int[] convertScores(String seqQual){
		int[] scores = new int[seqQual.length()];
		for (int i=0; i< seqQual.length(); i++){
			String sub = seqQual.substring(i, i+1);
			Integer val = asci2FastQScore.get(sub);
			if (val != null) scores[i] = val.intValue();
			else System.err.println("\nError converting seq quality character -> "+sub+" from "+seqQual);
		}
		return scores;
	}
	
	/**Attempts to extract chr1,2,3...22 or chrX,Y,M,MT from the String.
	 * Returns chromosome or null.*/
	public static String extractChromosomeName(String x){
		Matcher mat = chrNumber.matcher(x);
		if (mat.find()) return mat.group();
		mat = chrLetter.matcher(x);
		if (mat.find()) return mat.group();
		//check for controls
		if (x.indexOf("chrCtrls") !=-1) return "chrCtrls";
		return null;
	}
	
	/**Uses Seq.extractChromosomeName() to extract a chromosome text from each file.
	 * Adds this to a HashMap containing key: value, chrom text: File.
	 * Returns null if two files with the same chromosome text are found.
	 * Skips any files in which a chromosome text cannot be extracted.*/
	public static HashMap<String, File> makeChromosomeNameFileHash(File[] files){
		HashMap<String, File> hash = new HashMap<String, File>();
		for (int i=0; i< files.length; i++){
			String chromName = extractChromosomeName(files[i].getName());
			if (chromName != null) {
				if (hash.containsKey(chromName)) {
					return null;
				}
				hash.put(chromName, files[i]);
			}
		}
		return hash;
	}
	
	/**Converts a DNA sequence into a boolean[], everything not g or c is recorded as false.
	 * Returns null if a fasta is not found.*/
	public static boolean[] fetchGCContent(char[] chromosomeSequence){
		boolean[] gcContent = new boolean[chromosomeSequence.length];
		for (int i=0; i<chromosomeSequence.length; i++) {
				if (chromosomeSequence[i] == 'g' || chromosomeSequence[i] == 'c') gcContent[i] = true;
				else gcContent[i] = false;
			}
		return gcContent;
	}

	/**Writes a binary sequence, gatc, anything else is assumed to be n.
	 * @return true if sucessful, false if something bad happened.*/
	public static boolean writeBinarySequence(String seq, File file){
		try {
			FileOutputStream fos = new FileOutputStream(file);
			DataOutputStream dos = new DataOutputStream( new BufferedOutputStream (fos));
			//write text of array
			dos.writeUTF("Sequence");
			//find start indexes and the number of non gatc bases, save in ArrayList
			Pattern pat = Pattern.compile("[^gatc]");
			Matcher mat = pat.matcher(seq.toLowerCase());
			ArrayList nStartNumberPairs = new ArrayList(); 
			int numberOfNs = 1;
			int startIndex = 0;
			int index = 0;
			boolean open = false;
			while (mat.find()) {
				index = mat.start();	
				//continuing text of N's?
				if (open && index == (startIndex + numberOfNs)) numberOfNs++;
				//no, new start
				else{
					if (open){					
						//close old start; add startIndex and numberOfNs
						nStartNumberPairs.add(new Integer(startIndex));
						nStartNumberPairs.add(new Integer(numberOfNs));
					}
					else open = true;
					//open new
					startIndex = index;
					numberOfNs = 1;
				}	
			}
			//add last pair?
			if (open){				
				nStartNumberPairs.add(new Integer(startIndex));
				nStartNumberPairs.add(new Integer(numberOfNs));
			}
			int numberOfPairs = nStartNumberPairs.size();			
			//write number of N pairs
			dos.writeInt(numberOfPairs);
			//write index number pairs
			for (int x=0; x<numberOfPairs; x++){
				dos.writeInt(((Integer)nStartNumberPairs.get(x)).intValue());
			}
			//strip sequence of non gatc bases
			seq = mat.replaceAll("");
			//write size of stripped sequence, 
			int sizeStrippedSeq = seq.length();
			dos.writeInt(sizeStrippedSeq);
			//append c's to stop to bring up to a multiple of four
			int basesToAdd = 4 - sizeStrippedSeq % 4;			
			if (basesToAdd != 4){
				if (basesToAdd == 1) seq = seq + "c";
				else if (basesToAdd == 2) seq = seq + "cc";
				else seq = seq + "ccc";
			}
			int sizeFilledSequence = seq.length();
			//write size of filled in sequence
			dos.writeInt(sizeFilledSequence);
			//creat hashmap of text byte
			HashMap map = makeByte4BaseMap();
			//write 4 bp bytes
			for (int i=0; i< sizeFilledSequence; i+=4){
				String sub = seq.substring(i,i+4);
				dos.writeByte(((Byte)map.get(sub)).byteValue());
			}
			dos.close();
			fos.close();
			return true;
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return false; 
		}
	}
	
	/**Reads a binary sequence file returning gatc or n, lower case.
	 * @return null if something bad happened.*/
	public static String readBinarySequence(File file){
		try {
			FileInputStream fis = new FileInputStream(file);
			DataInputStream dis = new DataInputStream( new BufferedInputStream(fis ));
			//read array text
			if (dis.readUTF().equals("Sequence") == false) return null;
			//read number of N's startIndex:length pairs
			int numberOfPairs = dis.readInt();
			//read in startIndex, length pairs
			int[] nIndexPairs = null;
			int totalSizeNs = 0;
			if (numberOfPairs !=0){
				nIndexPairs = new int[numberOfPairs];
				for (int i=0; i< numberOfPairs; i++){
					//read start index
					nIndexPairs[i++]= dis.readInt();
					//read number of Ns
					nIndexPairs[i]= dis.readInt();
					//add to total Ns
					totalSizeNs += nIndexPairs[i];
					
				}
			}
			//read size of stripped sequence
			int sizeStrippedSeq = dis.readInt();
			//read size of filled sequence and convert to number of bytes
			int sizeFilledSequence = dis.readInt();
			int numberOfSeqBytes = sizeFilledSequence/4;
			//load sequence into StringBuffer
			StringBuffer sb = new StringBuffer(sizeFilledSequence + totalSizeNs);
			for (int i=0; i<numberOfSeqBytes; i++){
				byte b = dis.readByte();
				sb.append(all4BaseCombinations[b + 128]);
			}
			//remove filler
			sb.delete(sizeStrippedSeq, sizeFilledSequence);
			//insert n's
			if (nIndexPairs != null) {
				for (int i=0; i<numberOfPairs; i++) {
					//find start index and length
					int startIndex = nIndexPairs[i++];
					int numberNs = nIndexPairs[i];
					//build char[] of ns
					char[] ns = new char[numberNs];
					Arrays.fill(ns,'n');					
					//insert n's
					sb.insert(startIndex,ns);
				}
			}
			dis.close();
			fis.close();
			return sb.toString();
		}
		catch (Exception ioe){
			ioe.printStackTrace();
			return null;   
		} 
	}
	
	/**Makes a hash map of each 4bp combination in the all4BaseCombinations String[] and a unique Byte.*/
	public static HashMap makeByte4BaseMap (){
		HashMap map = new HashMap(256);
		for (int x=0; x<all4BaseCombinations.length; x++){
			byte b = (byte)(x -128);
			map.put(all4BaseCombinations[x], new Byte(b));
		}
		return map;
	}
	
	/**All possible 4 base combinations.*/
	public static final String[] all4BaseCombinations = {"gggg","ggga","gggt","gggc","ggag","ggaa","ggat","ggac",
		"ggtg","ggta","ggtt","ggtc","ggcg","ggca","ggct","ggcc","gagg","gaga","gagt","gagc","gaag","gaaa","gaat",
		"gaac","gatg","gata","gatt","gatc","gacg","gaca","gact","gacc","gtgg","gtga","gtgt","gtgc","gtag","gtaa",
		"gtat","gtac","gttg","gtta","gttt","gttc","gtcg","gtca","gtct","gtcc","gcgg","gcga","gcgt","gcgc","gcag",
		"gcaa","gcat","gcac","gctg","gcta","gctt","gctc","gccg","gcca","gcct","gccc","aggg","agga","aggt","aggc",
		"agag","agaa","agat","agac","agtg","agta","agtt","agtc","agcg","agca","agct","agcc","aagg","aaga","aagt",
		"aagc","aaag","aaaa","aaat","aaac","aatg","aata","aatt","aatc","aacg","aaca","aact","aacc","atgg","atga",
		"atgt","atgc","atag","ataa","atat","atac","attg","atta","attt","attc","atcg","atca","atct","atcc","acgg",
		"acga","acgt","acgc","acag","acaa","acat","acac","actg","acta","actt","actc","accg","acca","acct","accc",
		"tggg","tgga","tggt","tggc","tgag","tgaa","tgat","tgac","tgtg","tgta","tgtt","tgtc","tgcg","tgca","tgct",
		"tgcc","tagg","taga","tagt","tagc","taag","taaa","taat","taac","tatg","tata","tatt","tatc","tacg","taca",
		"tact","tacc","ttgg","ttga","ttgt","ttgc","ttag","ttaa","ttat","ttac","tttg","ttta","tttt","tttc","ttcg",
		"ttca","ttct","ttcc","tcgg","tcga","tcgt","tcgc","tcag","tcaa","tcat","tcac","tctg","tcta","tctt","tctc",
		"tccg","tcca","tcct","tccc","cggg","cgga","cggt","cggc","cgag","cgaa","cgat","cgac","cgtg","cgta","cgtt",
		"cgtc","cgcg","cgca","cgct","cgcc","cagg","caga","cagt","cagc","caag","caaa","caat","caac","catg","cata",
		"catt","catc","cacg","caca","cact","cacc","ctgg","ctga","ctgt","ctgc","ctag","ctaa","ctat","ctac","cttg",
		"ctta","cttt","cttc","ctcg","ctca","ctct","ctcc","ccgg","ccga","ccgt","ccgc","ccag","ccaa","ccat","ccac",
		"cctg","ccta","cctt","cctc","cccg","ccca","ccct","cccc"};


	
	/**Deletes any non IUP characters*/
	public static String filterDNASequence(String seq){
		return seq.replaceAll("[^ACGTBDHKMNRSVWYXacgtbdhkmnrsvwyx]","");
	}

	/**Deletes any non GATCNX characters*/
	public static String filterDNASequenceStrict(String seq){
		return seq.replaceAll("[^ACGTNXacgtnx]","");
	}
	
	/**Deletes any non IUP characters but leaves whitespaces*/
	public static String filterDNASeqLeaveWS(String seq){
		return seq.replaceAll("[^ACGTBDHKMNRSVWYXacgtbdhkmnrsvwyx\\s]","");
	}

	/**Returns a sub sequence given a relative start and stop, the bp for the first base. 
	 * Stop is inclusive. 
	 * 0 starts.*/
	public static String fetchSubSequence(int start, int stop, int bpFirstBase, String sequence){
		int begin = start-bpFirstBase;
		int end = stop-bpFirstBase;
		if (begin < 0 ) begin = 0;
		if (end < 0 ) return "";
		if (end > sequence.length()) return sequence.substring(begin);
		return sequence.substring(begin, end+1);
	}
	
	/**Returns a sub sequence given a relative start and stop, the bp for the first base. 
	 * Interbase coordinates, 0 starts, stop is excluded. */
	public static String fetchSubSequenceInterbaseCoordinates(int start, int stop, int bpFirstBase, String sequence){
		int begin = start-bpFirstBase;
		int end = stop-bpFirstBase;
		if (begin < 0 ) begin = 0;
		if (end < 0 ) return "";
		if (end > sequence.length()) return sequence.substring(begin);
		return sequence.substring(begin, end);
	}
	
	/**Generates identitiy dashes (ie"|||  | |||") between two aligned sequences*/
	public static String genDashes(String seq1, String seq2) {
		int len = seq1.length();
		StringBuffer dashes = new StringBuffer(len);
		for (int i = 0; i < len; i++) {
			Character s1 = new Character(seq1.charAt(i));
			Character s2 = new Character(seq2.charAt(i));
			int comp = s1.compareTo(s2);
			if (comp == 0) {
				dashes.append("|");
			} else {
				dashes.append(" ");
			}
		}
		return new String(dashes);
	}
	/**Given a sequence returns the number of G's,A's,T's,C's,N's as an int[5].
	 * Any non word chars are deleted, afterward any non GATC chars are counted 
	 * as Ns.  Case, space insensitive.*/
	public static int[] countBases(String sequence){
		int g =0;int a =0;int t =0;int c =0;int n =0;
		String seq = sequence.replaceAll("\\W","");
		seq = seq.toLowerCase();
		int len = seq.length();
		for (int i=0; i<len; i++){
			char test = seq.charAt(i);
			if(test=='g') g++;
			else if(test=='a') a++;
			else if(test=='t') t++;
			else if(test=='c') c++;
			else n++;
		}	
		return new int[]{g,a,t,c,n};
	}
	
	/**Calculates the fraction bases (GATCN) in the String*/
	public static double[] calculateFractionBases(String seq){
		//G's,A's,T's,C's,N's
		int[] counts = Seq.countBases(seq);
		double total = seq.length();
		double[] fractionGATC = new double[5];
		fractionGATC[0] = ((double)(counts[0]))/total;
		fractionGATC[1] = ((double)(counts[1]))/total;
		fractionGATC[2] = ((double)(counts[2]))/total;
		fractionGATC[3] = ((double)(counts[3]))/total;
		fractionGATC[4] = ((double)(counts[4]))/total;
		return fractionGATC;
	}

	/**Generates a matrix of the number of As Cc Gg Ts (top) by 1,2,3,4...positions in the
	 * motif (side) observed
	 * in all the Strings of the String[].  The String[] should contain equal length Strings,
	 * comprised entirely of GATC.*/
	public static int[][] makeFrequencyMatrix(String[] hits){
		int lenOfAHit = hits[0].length();
		int len = hits.length;
		int[][] matrix = new int[lenOfAHit][4];
		System.out.println("      \tAs\tCs\tGs\tTs");
		for (int i=0; i<lenOfAHit; i++){
			//set base counters
			int As =0;
			int Cs =0;
			int Gs =0;
			int Ts =0;
			//run through each hit looking at a particular base
			for (int j=0; j<len; j++){
				char test = hits[j].charAt(i);
				switch (test){
					case 'a':
					case 'A': As++; break;
					case 'c':
					case 'C': Cs++; break;
					case 'g':
					case 'G': Gs++; break;
					case 't':
					case 'T': Ts++; break;
					default: System.out.println("\n\nFatal error: odd base in makeFrequencyMactrix base -> "+test+"\n\n"); System.exit(0);
				}
			}
			System.out.println("Pos: "+(i+1)+"\t"+As+"\t"+Cs+"\t"+Gs+"\t"+Ts);
			int[] counts = {As, Cs, Gs, Ts};
			matrix[i] = counts;
		}
		return matrix;
	}



	/**Takes a DNA seq and reverse comps it, ambiguous symbols OK.
	 * Will warn if it finds an unrecognized base.  Works with ' GATCRYWSKMBDHVNX .- '
	 * case and space insensitive.*/	
	public static String reverseComplementDNA(String seq){
		int seqLen = seq.length();
		StringBuffer rcSeq = new StringBuffer(seqLen);
		char test;
		for (int i=seqLen-1; i>=0; i--){
			test = seq.charAt(i);
			switch (test){
				case 'A': rcSeq.append('T'); break;
				case 'C': rcSeq.append('G'); break;
				case 'G': rcSeq.append('C'); break;
				case 'T': rcSeq.append('A'); break;
				case 'c': rcSeq.append('g'); break;
				case 'a': rcSeq.append('t'); break;
				case 'g': rcSeq.append('c'); break;
				case 't': rcSeq.append('a'); break;
				case 'n':
				case 'N':
				case 'x':
				case 'X':
				case '-':
				case '.':
				case ' ':
				case 'S':
				case 's':
				case 'w':
				case 'W': rcSeq.append(test); break;
				case 'r': rcSeq.append('y'); break;
				case 'R': rcSeq.append('Y'); break;
				case 'y': rcSeq.append('r'); break;
				case 'Y': rcSeq.append('R'); break;
				case 'm': rcSeq.append('k'); break;
				case 'M': rcSeq.append('K'); break;
				case 'k': rcSeq.append('m'); break;
				case 'K': rcSeq.append('M'); break;
				case 'b': rcSeq.append('v'); break;
				case 'B': rcSeq.append('V'); break;
				case 'd': rcSeq.append('h'); break;
				case 'D': rcSeq.append('H'); break;
				case 'h': rcSeq.append('d'); break;
				case 'H': rcSeq.append('D'); break;
				case 'v': rcSeq.append('b'); break;
				case 'V': rcSeq.append('B'); break;
				default: rcSeq.append(test); System.err.println("\nWarning: odd base in revComp-> '"+test+
				"' from "+seq+" Reverse Complement possibly incorrect!\n");
			}
		}
		return rcSeq.toString();
	}
	
	/**Takes a RNA seq and reverse comps it, ambiguous symbols OK.
	 * Will warn if it finds an unrecognized base.  Works with ' GATCURYWSKMBDHVNX .- '
	 * case and space insensitive.*/	
	public static String reverseComplementRNA(String seq){
		int seqLen = seq.length();
		StringBuffer rcSeq = new StringBuffer(seqLen);
		char test;
		for (int i=seqLen-1; i>=0; i--){
			test = seq.charAt(i);
			switch (test){
				case 'A': rcSeq.append('U'); break;
				case 'C': rcSeq.append('G'); break;
				case 'G': rcSeq.append('C'); break;
				case 'U': rcSeq.append('A'); break;
				case 'T': rcSeq.append('A'); break;
				case 'c': rcSeq.append('g'); break;
				case 'a': rcSeq.append('u'); break;
				case 'g': rcSeq.append('c'); break;
				case 't': rcSeq.append('a'); break;
				case 'u': rcSeq.append('a'); break;
				case 'n':
				case 'N':
				case 'x':
				case 'X':
				case '-':
				case '.':
				case ' ':
				case 'S':
				case 's':
				case 'w':
				case 'W': rcSeq.append(test); break;
				case 'r': rcSeq.append('y'); break;
				case 'R': rcSeq.append('Y'); break;
				case 'y': rcSeq.append('r'); break;
				case 'Y': rcSeq.append('R'); break;
				case 'm': rcSeq.append('k'); break;
				case 'M': rcSeq.append('K'); break;
				case 'k': rcSeq.append('m'); break;
				case 'K': rcSeq.append('M'); break;
				case 'b': rcSeq.append('v'); break;
				case 'B': rcSeq.append('V'); break;
				case 'd': rcSeq.append('h'); break;
				case 'D': rcSeq.append('H'); break;
				case 'h': rcSeq.append('d'); break;
				case 'H': rcSeq.append('D'); break;
				case 'v': rcSeq.append('b'); break;
				case 'V': rcSeq.append('B'); break;
				default: rcSeq.append(test); System.err.println("\nWarning: odd base in revComp-> '"+test+
				"' from "+seq+" Reverse Complement possibly incorrect!\n");
			}
		}
		return rcSeq.toString();
	}
	
	/**Takes a DNA seq and complement it, ambiguous symbols OK.
	 * Will warn if it finds an unrecognized base.  Works with ' GATCRYWSKMBDHVNX .- '
	 * case and space insensitive.*/	
	public static String complementDNA(String seq){
		int seqLen = seq.length();
		StringBuffer cSeq = new StringBuffer(seqLen);
		char test;
		for (int i=0; i<seqLen; i++){
			test = seq.charAt(i);
			switch (test){
				case 'a': cSeq.append('t'); break;
				case 'A': cSeq.append('T'); break;
				case 'c': cSeq.append('g'); break;
				case 'C': cSeq.append('G'); break;
				case 'g': cSeq.append('c'); break;
				case 'G': cSeq.append('C'); break;
				case 't': cSeq.append('a'); break;
				case 'T': cSeq.append('A'); break;
				case 'n':
				case 'N':
				case 'x':
				case 'X':
				case '-':
				case '.':
				case ' ':
				case 'S':
				case 's':
				case 'w':
				case 'W': cSeq.append(test); break;
				case 'r': cSeq.append('y'); break;
				case 'R': cSeq.append('Y'); break;
				case 'y': cSeq.append('r'); break;
				case 'Y': cSeq.append('R'); break;
				case 'm': cSeq.append('k'); break;
				case 'M': cSeq.append('K'); break;
				case 'k': cSeq.append('m'); break;
				case 'K': cSeq.append('M'); break;
				case 'b': cSeq.append('v'); break;
				case 'B': cSeq.append('V'); break;
				case 'd': cSeq.append('h'); break;
				case 'D': cSeq.append('H'); break;
				case 'h': cSeq.append('d'); break;
				case 'H': cSeq.append('D'); break;
				case 'v': cSeq.append('b'); break;
				case 'V': cSeq.append('B'); break;
				default: cSeq.append(test); System.err.println("\nWarning: odd base in comp-> '"+test+
				"' Complement possibly incorrect!\n");
			}
		}
		return cSeq.toString();
	}


}
