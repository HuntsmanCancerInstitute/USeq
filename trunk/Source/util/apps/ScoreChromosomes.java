package util.apps;
import meme.*;
import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.parsers.*;
import util.bio.parsers.*;
import util.gen.*;

/**
 * Scores chromosomes for the presence of transcription factor binding sites.
 *
 */
public class ScoreChromosomes {

	//fields
	private File flySeqDir;
	private File[] chromosomeFiles;
	private File bindingSiteFile;
	private String[] bindingSites;		//aligned and trimmed
	private MotifScanner motifScanner;
	private double cutOff = -1;
	private long totalNumHits = 0;
	private long totalSeqLength = 0;
	private boolean printBarFiles = false;
	private File barFileDirectory = null;
	private String versionedGenome;
	private boolean printHits = false;

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ScoreChromosomes(args);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Score Chromosomes: Oct  2012                           **\n" +
				"**************************************************************************************\n" +
				"SC scores chromosomes for the presence of transcription factor binding sites. Use the\n" +
				"following options:\n\n" +

				"-g The full path directory text to the split genomic sequences (i.e. chr2L.fasta, \n"+
				"      chr3R.fasta...), FASTA format.\n" +
				"-t Full path file text for the FASTA file containing aligned trimmed examples of\n"+
				"      transcription factor binding sites.  A log likelihood position specific\n"+
				"      probability matrix will be generated from these sequences and used to scan the\n"+
				"      chromosomes for hits to the matrix.\n"+
				"-s Score cut off for the matrix. Defaults to the score of the lowest scoring sequence\n"+
				"      used in making the LLPSPM.\n"+
				"-p Print hits to screen, default is no.\n"+
				"-v Provide a versioned genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases, if you would like to write graph LLPSPM\n" +
				"      scores in xxx.bar format for direct viewing in IGB.\n" +
				"\n" +
				"Example: java -Xmx4000M -jar pathTo/T2/Apps/ScoreChromosomes -g /my/affy/Hg18Seqs/ -t \n" +
				"      /my/affy/fgf8.fasta -s 4.9 -v H_sapiens_Mar_2006\n\n" +

		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign any new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'g': flySeqDir = new File(args[i+1]); i++; break;
					case 't': bindingSiteFile = new File(args[i+1]); i++; break;
					case 's': cutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'v': versionedGenome = args[++i]; printBarFiles = true; break;
					case 'p': printHits = true; break;
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

		if (flySeqDir ==null || flySeqDir.isDirectory()==false){
			System.out.println("\nEnter the full path text for the directory  containing split genomic sequences " +
					"(i.e. chr2L.fasta, chr3R.fasta...).  The file prefex names should match the chromosome names " +
			"in the Intervals (ie chr2L, chr3R...)\n");
			System.exit(0);
		}
		if (bindingSiteFile==null || bindingSiteFile.exists()==false){
			System.out.println("\nEnter a full path file text for the fasta file containing aligned trimmed examples of transcription factor binding sites.\n");
			System.exit(0);
		}
	}



	public ScoreChromosomes(String[] args){
		processArgs(args);
		//set genomic chromosome seq directory and files
		chromosomeFiles = IO.extractFiles(flySeqDir, ".fasta");
		if (chromosomeFiles == null || chromosomeFiles.length == 0) chromosomeFiles = IO.extractFiles(flySeqDir, ".fasta.gz");

		//get binding sites
		MultiFastaParser fastaParser = new MultiFastaParser(bindingSiteFile);
		bindingSites = fastaParser.getSeqs();
		//convert binding sites to upper case
		for (int i=0; i< bindingSites.length; i++) bindingSites[i] = bindingSites[i].toUpperCase();

		//make a MotifScanner
		double[][] PSPM = MemeMotif.makePSPM(bindingSites);
		double[][] llPSPM = MemeMotif.makeLLPSPM(PSPM, 0.25, 0.25, 0.25, 0.25);
		motifScanner = new MotifScanner(llPSPM);

		if (cutOff == -1) {
			cutOff = motifScanner.findLowestScoringSeq(bindingSites);
			System.out.println("\nScore cut off set to lowest binding sequence -> "+cutOff);
		}
		

		//make a bar file directory?
		if (printBarFiles){
			barFileDirectory = new File (bindingSiteFile.getParent(), "ScrChrmBarFiles_"+Misc.removeExtension(bindingSiteFile.getName())+"_"+Num.formatNumber(cutOff, 3));
			barFileDirectory.mkdir();
		}

		//scan em
		scanChromosomes();

		System.out.println("\nFor a cutoff of "+cutOff);
		System.out.println("\tTotal # Hits: "+totalNumHits);
		System.out.println("\tTotal Length: "+totalSeqLength);

		double aveKb = 1000*((double)totalNumHits/(double)totalSeqLength);
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(3);
		System.out.println("\tHits per kb : "+f.format(aveKb));
		System.out.println();
	}


	public void scanChromosomes(){
		MotifHit[] hits;
		totalNumHits = 0;
		totalSeqLength = 0;
		String chromSeq;
		int seqLength;
		int numHits;
		File file;
		//make bar parser?
		BarParser bp = null;
		if (printBarFiles){
			bp = new BarParser();
			//set hashmap
			HashMap<String,String> map = new HashMap<String,String>();
			map.put(BarParser.DESCRIPTION_TAG, "ScoreChromosomes results with a LLPSP threshold of "+cutOff);
			map.put(BarParser.SOURCE_TAG, "Matrix sequence file "+bindingSiteFile.getName());
			map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			map.put(BarParser.UNIT_TAG, "Log likelihood position specific probablilities.");
			bp.setTagValues(map);
			bp.setVersionedGenome(versionedGenome);
			bp.setZipCompress(true);
		}

		//for each chromosome
		for (int i=0; i<chromosomeFiles.length; i++){
			//get chromosome
			file = chromosomeFiles[i];
			System.out.println("\nScanning chromosome: "+file.getName());
			MultiFastaParser mfp = new MultiFastaParser(file);
			chromSeq = mfp.getSeqs()[0].toUpperCase();
			seqLength = chromSeq.length();
			totalSeqLength += seqLength;
			System.out.println("\tLength: "+seqLength);		
			hits = motifScanner.scoreSequence(cutOff,chromSeq);
			//make bar files?
			if (printBarFiles && hits.length !=0){
				//load centered positions and scores
				int[] positions = new int[hits.length];
				float[] values = new float[hits.length];
				for (int x=0; x< hits.length; x++){
					positions[x] = (int)(Math.round(((double)(hits[x].getStop() - hits[x].getStart()))/2)) + hits[x].getStart();
					values[x] = new Double(hits[x].getScore()).floatValue();
					//if (values[x] > 11.1f) System.out.println(hits[x]);
				}
				bp.setBasePositions(positions);
				bp.setValues(values);
				String chr = Misc.removeExtension(file.getName());
				bp.setChromosome(chr);
				bp.setBarFile(new File (barFileDirectory, chr+".bar"));
				try {
					bp.writeSimpleBarFile();
				} catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nProblem writing bar files for chromosome "+chr);
				}
			}
			//print hits?
			if (printHits){
				String chr = Misc.removeExtension(file.getName())+"\t";
				String name = "\t"+Misc.capitalizeFirstLetter(Misc.removeExtension(bindingSiteFile.getName()))+"\t";
				//chr, start, stop, ,text, score, strand, seq
				for (int x=0; x<hits.length; x++){
					String strand = "\t+\t";
					if (hits[x].getOrientation() == 1) strand = "\t-\t";
					System.out.println(chr+hits[x].getStart()+"\t"+hits[x].getStop()+name+hits[x].getScore()+strand+hits[x].getSeq());
				}
			}
			numHits = hits.length;
			System.out.println("\t# Hits: "+numHits);
			//calculate median score
			if (numHits>0){
				double[] scores = new double[hits.length];
				for (int x=0; x< hits.length; x++) scores[x] = hits[x].getScore();
				Arrays.sort(scores);
				System.out.println("\tMedian Score: "+Num.median(scores));
			}
			totalNumHits +=numHits;

		}
	}


}
