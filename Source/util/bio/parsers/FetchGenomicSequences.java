package util.bio.parsers;

import trans.anno.*;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import trans.main.*;
import util.bio.calc.*;

import util.bio.parsers.MultiFastaParser;
import util.gen.*;
import util.bio.seq.*;
/**
 * Fetches the genomic sequences for a list of regions.
 *
 */
public class FetchGenomicSequences {
	//fields
	private File[] firstRegionsFiles = null;
	private BindingRegion[] one;
	private File genomicSequenceDirectory = null;
	private int addSubtract = 0;
	private boolean reverseComplement = false;
	private boolean fastaFormat = false;

	public FetchGenomicSequences(String[] arguments){
		processArgs(arguments);
		System.out.println("\nLaunching...");
		//for each file
		for (int x=0; x<firstRegionsFiles.length; x++){
			//parse BindingRegions s
			one = IntersectRegions.parseRegionsFile(firstRegionsFiles[x]);
			//sort by chromosome
			Arrays.sort(one);
			System.out.println("\nLoading for "+firstRegionsFiles[x].getName());
			loadSequences();
			//sort by original order
			Arrays.sort( one, new BindingRegionComparator());
			//print to file
			String fasta = ".txt";
			if (fastaFormat) fasta = ".fasta";
			File loaded = new File (firstRegionsFiles[x].getParentFile(), Misc.removeExtension(firstRegionsFiles[x].getName())+"_FGS"+fasta);
			StringBuffer sb = new StringBuffer();
			for (int i=0; i< one.length; i++){
				sb.append(one[i].getChromosome());
				sb.append("\n");
			}
			IO.writeString(sb.toString(), loaded);
		}
	}

	public void loadSequences(){
		String chrom = "";
		String chromosomeSequence = null;
		for (int i=0; i< one.length; i++){
			//different chromosome? load info
			if (one[i].getChromosome().equals(chrom) == false){
				chrom = one[i].getChromosome();
				if (i!=0) System.out.println();
				System.out.println("\tLoading "+chrom);
				File chromFastaFile = new File (genomicSequenceDirectory,chrom+".fasta");
				if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,chrom+".fa");
				if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,chrom+".fasta.gz");
				if (chromFastaFile.exists() ==  false ) chromFastaFile = new File (genomicSequenceDirectory,chrom+".fa.gz");
				if (chromFastaFile.exists() ==  false ){
					Misc.printExit("Error: Cannot find genomic sequence file for -> "+chrom +" Aborting.");
				}
				MultiFastaParser fastaParser = new MultiFastaParser();
				fastaParser.parseIt(chromFastaFile);
				chromosomeSequence = fastaParser.getSeqs()[0];
			}
			//calc ends
			int end = one[i].getEnd()+ addSubtract;
			boolean complete = true;
			if (end >= chromosomeSequence.length()){
				end = chromosomeSequence.length();
				complete = false;
			}
			int start = one[i].getStart() - addSubtract;
			if (start < 0) {
				start = 0;
				complete = false;
			}
			String seq = new String (chromosomeSequence.substring(start, end)); //needed to force a new text and not a reference to something huge
			if (reverseComplement) seq = Seq.reverseComplementDNA(seq);
			String sumLine;
			if (fastaFormat){
				sumLine = ">"+one[i].getChromosome()+":"+one[i].getStart()+"-"+one[i].getEnd()+"_"+ complete+ "\n"+ seq;
			}
			else sumLine = one[i].simpleSummaryLine()+"\t"+start+"\t"+ end+"\t"+ complete+ "\t"+ seq;
			System.out.print(".");
			one[i].setChromosome(sumLine);
		}
		System.out.println();
	}





	public static void main(String[] args){
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new FetchGenomicSequences(args);
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': firstRegionsFiles = IO.extractFiles(new File (args[++i])); break;
					case 'b': addSubtract = Integer.parseInt(args[i+1]); i++; break;
					case 's': genomicSequenceDirectory = new File (args[i+1]); i++; break;
					case 'r': reverseComplement = true; break;
					case 'a': fastaFormat = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//look for required parameters
		if (firstRegionsFiles == null || genomicSequenceDirectory == null ){
			Misc.printExit("\nPlease complete one or more of the following required parameters: -f or -c .\n");
		}

		
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           FetchGenomicSequences: Feb 2013                        **\n" +
				"**************************************************************************************\n" +
				"Given a file containing genomic coordinates, fetches and saves the sequence (column\n" +
				"output: chrom origStart origStop fetchedStart fetchedStop completeFetch seq).\n\n"+

				"-f Full path to a file or directory containing tab delimited chrom, start,\n" +
				"        stop text files.  Interbabase coordinates (zero based, stop excluded).\n"+
				"-s Full path directory text containing containing genomic fasta files. The fasta\n" +
				"        header defines the name of the sequence, not the file name. \n"+
				"-b Fetch flanking bases, defaults to 0. Will set start to zero or stop to last base if\n" +
				"        boundaries are exceeded.\n"+
				"-r Reverse complement fetched sequences, defaults to returning the + genomic strand.\n"+
				"-a Output fasta format.\n"+


				"\nExample: java -Xmx1000M -jar pathTo/T2/Apps/FetchGenomicSequences -f /data/miRNAs.txt\n" +
				"      -s /genomes/human/v35.1/ -b 5000 -r   \n\n" +

				"\n" +
		"**************************************************************************************\n");		
	}

	/**@return double[2]{total bases,median length}.*/
	public static double[] statBindingRegionArray(BindingRegion[] brs){
		//total bases
		double totalBases = 0;
		double[] lengths = new double[brs.length];
		for (int i=0; i< brs.length; i++){
			lengths[i] = brs[i].getLength();
			totalBases += lengths[i];
		}
		//calculate median length
		Arrays.sort(lengths);
		double medianLength = Num.median(lengths);
		return new double[] {totalBases, medianLength};
	}

	/**@return total number of bases in the BindingRegion[], stop-start+1.*/
	public static int numberBases(BindingRegion[] brs){
		//total bases
		int totalBases = 0;
		for (int i=0; i< brs.length; i++){
			totalBases += brs[i].getLength();
		}
		//calculate median length
		return totalBases;
	}

}
