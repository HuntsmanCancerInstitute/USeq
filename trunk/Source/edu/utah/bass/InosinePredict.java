package edu.utah.bass;

import java.io.*;
import util.gen.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.parsers.*;
import util.bio.seq.Seq;

/**@author david.nix@hci.utah.edu*/
public class InosinePredict {

	//fields
	private File matrixFile;
	private File sequenceFastaFile;
	private double[][] matrix;
	private String matrixName;
	private boolean includeOppositeStrand = true;
	private File[] pngFiles;
	private File[] epsFiles;
	private int[] widths;
	private File table;
	private File zipArchive;
	private File saveDirectory;
	private MultiFastaParser mfp;

	//constructors
	/**For integration with InosinePredictServlet.java*/
	public InosinePredict (File sequenceFile, MultiFastaParser mfp, File matrixFile, File saveDirectory, File zipArchive){
		this.sequenceFastaFile = sequenceFile;
		this.mfp = mfp;
		this.matrixFile = matrixFile;
		this.saveDirectory = saveDirectory;
		this.zipArchive = zipArchive;
		
		launch();
		
	}
	
	/**Stand alone use*/
	public InosinePredict (String[] args){
		processArgs(args);

		//load the fasta file
		mfp = new MultiFastaParser(sequenceFastaFile);
		if (mfp.isFastaFound() == false){
			Misc.printErrAndExit("\nError: your (multi) fasta file is malformed, correct and restart.\n");
		}
		
		//print some info
		System.out.println("Fasta file -> " +sequenceFastaFile);
		System.out.println("Scoring against -> " +matrixName+" matrix");
		System.out.println("\nNegative values/ gray bases denote unscorable bases due to end proximity or ambiguous bps in the score block.");
		System.out.println("Graphs are scaled to the % editing at 20% A -> I.\n");

		launch();
		
		System.out.println("Done!\n");
	}
	
	public void launch(){
		//load matrix
		boolean loaded = loadMatrix();
		if (loaded == false){
			System.err.println("\nError: your matrix is malformed see the example:\n\n"+exampleMatrix+"\n");
			return;
		}
		
		//make a print writer for the results table
		PrintWriter out = null;
		try {
			table = new File (saveDirectory, Misc.removeExtension(sequenceFastaFile.getName())+".xls");
			out = new PrintWriter( new FileWriter( table));
			out.println("#Fasta file -> " +sequenceFastaFile);
			out.println("#Scoring against -> " +matrixName);
			out.println("#Negative values denote unscorable A bases due to proximity to an end or ambiguous bps in the score block.\n");

			//for each sequence
			String[] headers = mfp.getNames();
			String[] seqs = mfp.getSeqs();
			pngFiles = new File[seqs.length];
			epsFiles = new File[seqs.length];
			widths = new int[seqs.length];
			char[] compSeq = null;
			double[] compBaseProb = null;

			for (int i=0; i< seqs.length; i++){
				//calculate probabilities
				double[] baseProb = scanSeq(seqs[i]);

				if (includeOppositeStrand){
					String revCompSeq = Seq.reverseComplementRNA(seqs[i]);
					double[] baseProbRevComp = scanSeq(revCompSeq);
					//reverse the order so 3' -> 5'
					compSeq = new StringBuffer(revCompSeq).reverse().toString().toCharArray();
					compBaseProb = new double[baseProbRevComp.length];
					int index =0;
					for (int x= baseProbRevComp.length-1; x>= 0; x--) {
						compBaseProb[index] = baseProbRevComp[x];
						index++;
					}
				}

				//print png
				String name = Misc.removeExtension(headers[i].trim().replaceAll("\\s+", "_"));
				pngFiles[i] = new File(saveDirectory, name+"_"+i+".png");
				epsFiles[i] = new File(saveDirectory, name+"_"+i+".eps");
				SequenceGraph sg = new SequenceGraph(baseProb, compBaseProb, seqs[i].toCharArray(), compSeq, pngFiles[i], epsFiles[i], matrixName);
				widths[i] = sg.getWidth();

				//print results, -1's could not be scored (too close to an end or an ambigous base found in pattern area)
				out.println(">"+headers[i]);
				char[] s = seqs[i].toCharArray();
				for (int x=0; x< s.length; x++){
					out.print(s[x]+"\t");
				}
				out.println();
				for (int x=0; x< s.length; x++){
					out.print(baseProb[x]+"\t");
				}
				out.println();
				
				if (includeOppositeStrand){
					for (int x=0; x< compSeq.length; x++){
						out.print(compSeq[x]+"\t");
					}
					out.println();
					for (int x=0; x< compBaseProb.length; x++){
						out.print(compBaseProb[x]+"\t");
					}
					out.println();
				}
				out.println();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			out.close();
		}
		
		//archive
		if (zipArchive != null){
			File[][] toCollapse = new File[4][];
			toCollapse[0] = pngFiles;
			toCollapse[1] = epsFiles;
			toCollapse[2] = new File[]{table};
			toCollapse[3] = new File[]{sequenceFastaFile};
			File[] toZip = IO.collapseFileArray(toCollapse);
			IO.zip(toZip, zipArchive);
		}
	}

	public double[] scanSeq(String seq){
		double[] baseProb = new double[seq.length()];
		int len = seq.length()-4;
		int[] seqIndexes = seqToIndex(seq.toLowerCase());

		//scan first 4 for A's
		for (int j=0; j<4; j++){
			//for each A, set to -1, cannot be scored
			if (seqIndexes[j] ==0) baseProb[j] = -1;
		}
		//scan middle
		for (int j = 4; j < len; j++){
			//for each A
			if (seqIndexes[j] ==0){
				double[] p = new double[8];
				p[0] = matrix[seqIndexes[j-4]][0] ;
				p[1] = matrix[seqIndexes[j-3]][1] ;
				p[2] = matrix[seqIndexes[j-2]][2] ;
				p[3] = matrix[seqIndexes[j-1]][3] ;
				p[4] = matrix[seqIndexes[j+1]][4] ;
				p[5] = matrix[seqIndexes[j+2]][5] ;
				p[6] = matrix[seqIndexes[j+3]][6] ;
				p[7] = matrix[seqIndexes[j+4]][7];
				baseProb[j] = calculateProb(p);
			}
		}
		//scan last
		for (int j=len; j<seq.length(); j++){
			//for each A, set to -1, cannot be scored
			if (seqIndexes[j] ==0) baseProb[j] = -1;
		}
		return baseProb;
	}

	public static double calculateProb( double[] p){
		Arrays.sort(p);
		if (p[0] == -1) return -1;
		else {
			double prod = p[0];
			for (int i=1; i< p.length; i++) prod = prod * p[i];
			return (prod *20);
		}
	}

	public static int[] seqToIndex(String seqLowerCase){
		int len = seqLowerCase.length();
		int[] indexes = new int[len];
		for (int i=0; i< len; i++){
			char c = seqLowerCase.charAt(i);
			if (c == 'a') indexes[i] = 0;
			else if (c == 'c') indexes[i] = 1;
			else if (c == 'g') indexes[i] = 2;
			else if (c == 'u') indexes[i] = 3;
			else if (c == 't') indexes[i] = 3;
			else indexes[i] = 4;
		}
		return indexes;
	}


	public static final String exampleMatrix = 
		"#Example matrix, tab delimited, order sensitive\n"+
		"hADAR1\tL4\tL3\tL2\tL1\tR1\tR2\tR3\tR4\n"+
		"A\t0.71194\t1.61822\t0.57216\t23.75355\t1.42806\t0.91581\t0.88952\t1.17913\n"+
		"C\t1.12399\t1.3486\t0.53298\t5.86694\t1.46775\t0.90437\t0.95155\t1.04024\n"+
		"G\t0.87164\t1.4734\t0.54996\t1.19222\t2.13172\t0.949\t0.93087\t1.10156\n"+
		"U\t1\t1\t1\t44.03081\t1\t1\t1\t1\n";

	/**For parsing a matrix of cooeficients.
	 * hADAR1	L4	L3	L2	L1	R1	R2	R3	R4
	 * A	0.71194	1.61822	0.57216	23.75355	1.42806	0.91581	0.88952	1.17913
	 * C	1.12399	1.3486	0.53298	5.86694	1.46775	0.90437	0.95155	1.04024
	 * G	0.87164	1.4734	0.54996	1.19222	2.13172	0.949	0.93087	1.10156
	 * U	1	1	1	44.03081	1	1	1	1*/
	public boolean loadMatrix(){
		try {
			Pattern tab = Pattern.compile("\t");
			matrix = new double[5][8];
			boolean goodMatrix = false;
			BufferedReader in = IO.fetchBufferedReader(matrixFile);
			String line;
			String[] tokens;
			String prior = null;
			while ((line= in.readLine()) != null){
				if (line.startsWith("A\t")){
					//parse name
					if (prior != null) {
						tokens = tab.split(prior);
						matrixName = tokens[0];
					}
					else matrixName = Misc.removeExtension(matrixFile.getName());
					//Parse A
					tokens = tab.split(line);
					for (int x=0; x<8; x++) matrix[0][x] = Double.parseDouble(tokens[x+1]);
					//Parse C
					line= in.readLine();
					tokens = tab.split(line);
					if (tokens[0].equals("C") == false) return false;
					for (int x=0; x<8; x++) matrix[1][x] = Double.parseDouble(tokens[x+1]);
					//Parse G
					line= in.readLine();
					tokens = tab.split(line);
					if (tokens[0].equals("G") == false) return false;
					for (int x=0; x<8; x++) matrix[2][x] = Double.parseDouble(tokens[x+1]);
					//Parse U
					line= in.readLine();
					tokens = tab.split(line);
					if (tokens[0].equals("U") == false) return false;
					for (int x=0; x<8; x++) matrix[3][x] = Double.parseDouble(tokens[x+1]);
					//set last row as -1 for bad bases
					for (int x=0; x<8; x++) matrix[4][x] = -1;
					goodMatrix = true;
				}
				prior = line;
			}
			in.close();
			return goodMatrix;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}


	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new InosinePredict(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': matrixFile = new File(args[++i]); break;
					case 'f': sequenceFastaFile = new File(args[++i]); break;
					case 'p': Misc.printExit("\n"+exampleMatrix+"\n");
					case 'o': includeOppositeStrand = false; break;
					case 's': saveDirectory = new File(args[++i]); saveDirectory.mkdirs(); break;
					case 'z': zipArchive = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		if (matrixFile == null || matrixFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your matrix file? -> "+matrixFile+"\n");
		if (sequenceFastaFile == null || sequenceFastaFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your sequence fasta file? -> "+sequenceFastaFile+"\n");
		if (saveDirectory == null) saveDirectory = sequenceFastaFile.getParentFile();
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Inosine Predict: Aug 2010                          **\n" +
				"**************************************************************************************\n" +
				"IP estimates the likelihood of ADAR RNA editing using the multiplicative 4L,4R model\n" +
				"described in Eggington et. al. 2010.\n" +

				"\nOptions:\n"+
				"-f Multi fasta file containing sequence(s) to score.\n"+
				"-m Maxtrix scoring file.\n" +
				"-p Print an example matrix.\n" +
				"-o Don't include the opposite strand.\n"+
				"-s Save directory, defaults to parent of the fasta file.\n"+
				"-z Name of a zip archive to create containing the results.\n"+

				"\n"+
				"Example: java -Xmx2G -jar pathTo/USeq/Apps/InosinePredict -m \n" +
				"    ~/ADARMatrix/hADAR1-D.matrix.txt -f ~/SeqsToScore/candidates.fasta.gz\n\n" +

		"**************************************************************************************\n");		
	}

	public File getMatrixFile() {
		return matrixFile;
	}

	public File getSequenceFastaFile() {
		return sequenceFastaFile;
	}

	public double[][] getMatrix() {
		return matrix;
	}

	public String getMatrixName() {
		return matrixName;
	}

	public boolean isIncludeOppositeStrand() {
		return includeOppositeStrand;
	}

	public File[] getPngFiles() {
		return pngFiles;
	}

	public File[] getEpsFiles() {
		return epsFiles;
	}

	public File getTable() {
		return table;
	}

	public File getZipArchive() {
		return zipArchive;
	}

	public File getSaveDirectory() {
		return saveDirectory;
	}

	public MultiFastaParser getMfp() {
		return mfp;
	}

	public static String getExamplematrix() {
		return exampleMatrix;
	}

	public int[] getWidths() {
		return widths;
	}
}
