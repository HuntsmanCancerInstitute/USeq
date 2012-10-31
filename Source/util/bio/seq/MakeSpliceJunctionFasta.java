package util.bio.seq;
import java.io.*;
import java.util.*;
import java.util.regex.*;

import util.gen.*;
import util.bio.parsers.*;
import util.bio.annotation.*;

/**Generates a fasta file for each.*/
public class MakeSpliceJunctionFasta {

	//fields
	private HashMap<String,File> fastaFiles;
	private int sequenceLengthRadius = 0;
	private HashMap<String,UCSCGeneLine[]> chromSpecificGeneLines;
	private File resultsFastaFile;
	private File resultsGeneFile;
	private boolean addExonSequence = false;
	private String[] ns;
	

	public MakeSpliceJunctionFasta(String[] args){
		//process user args and load data
		processArgs(args);

		//for each chromosome of gene models
		try{
			PrintWriter out = new PrintWriter (new FileWriter (resultsFastaFile));
			PrintWriter outGenes = new PrintWriter (new FileWriter (resultsGeneFile));
			System.out.println("Processing...");
			Iterator<String> it = chromSpecificGeneLines.keySet().iterator();
			while (it.hasNext()){
				String chr = it.next();
				System.out.println("\n\t" + chr);

				//fetch the sequence
				MultiFastaParser mfp = new MultiFastaParser(fastaFiles.get(chr));
				String seq = mfp.getSeqs()[0];

				//make hashmap to hold junction names
				HashSet<String> junctions = new HashSet<String>();

				//for each gene line
				UCSCGeneLine[] models = chromSpecificGeneLines.get(chr);
				boolean printCheck = true;
				for (int i=0; i< models.length; i++){
					ExonIntron[] exons = models[i].getExons();
					if (printCheck && exons.length != 1) {
						String name = models[i].getDisplayName();
						if (name== null) name = models[i].getName();
						else name = name+" "+models[i].getName();
						System.out.println("Check "+name +", "+exons.length + " exons, "+models[i].getStrand()+" strand");
					}
					//any junctions?
					if (exons.length == 1) continue;
					//save header to gene splice file
					outGenes.print(chr+"\t"+models[i].getName()+"\t");
					//for each pairing
					for (int j=0; j< exons.length; j++){
						//define the stop
						int end = exons[j].getEnd();
						//create left side sequence, watch for exon overrun
						int begin = end-sequenceLengthRadius;
						if (begin < 0) begin = 0;
						String leftSide = seq.substring(begin, end);
						int diffLeft = sequenceLengthRadius - exons[j].getLength();
						if (begin !=0 && diffLeft > 0){
							leftSide = ns[diffLeft]+ leftSide.substring(diffLeft);
						}

						//for all subsequent starts
						for (int k=j+1; k< exons.length; k++){
							int start = exons[k].getStart();

							//make text and see if already printed
							String name = chr+"_"+end+"_"+start;
							//save to gene file
							outGenes.print(end+"_"+start+",");
							if (junctions.contains(name)) continue;
							junctions.add(name);
							
							//make rightside
							int stop = start+sequenceLengthRadius;
							if (stop >= seq.length()) stop = seq.length()-1;
							String rightSide = seq.substring(start, stop);
							int diffRight = sequenceLengthRadius - exons[k].getLength();
							if (diffRight > 0){
								rightSide = rightSide.substring(0, exons[k].getLength()) + ns[diffRight];
							}
							
							//make seqFusion, problem with grabbing intronic
							String seqFus = leftSide + rightSide;

							//print
							out.println(">"+name);
							out.println(seqFus);
							if (printCheck){
								System.out.println("\t>"+name);
								System.out.println("\t"+leftSide +" "+ rightSide);
							}
						}
					}
					//close line
					outGenes.println();
					printCheck = false;
				}
			}
			out.close();
			outGenes.close();
		} catch (Exception e){
			e.printStackTrace();
		}


	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MakeSpliceJunctionFasta(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		File[] fastas = null;
		File ucscTableFile = null;
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fastas = IO.extractFiles(new File (args[++i])); break;
					case 's': sequenceLengthRadius = Integer.parseInt(args[++i]); break;
					case 'u': ucscTableFile = new File (args[++i]); break;
					case 'r': resultsFastaFile = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check params
		if (sequenceLengthRadius == 0) Misc.printExit("\nPlease enter a sequence length radius.\n");
		if (ucscTableFile == null || ucscTableFile.canRead() == false) Misc.printExit("\nCannot find or read your UCSC gene table?!\n");
		if (fastas == null || fastas.length ==0) Misc.printExit("\nCannot find your chromosome specific xxx.fasta files?\n");
		if (resultsFastaFile == null) Misc.printExit("\nPlease provide a full path file text for saving the spliced fasta results file.\n");
		//load UCSC table
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(ucscTableFile,0);
		chromSpecificGeneLines = reader.getChromSpecificGeneLines();
		//load fastaFiles
		fastaFiles = new HashMap<String,File>();
		Pattern chrom = Pattern.compile("(.+)\\.fa.*");
		for (int i=0; i< fastas.length; i++){
			Matcher mat = chrom.matcher(fastas[i].getName());
			if (mat.matches()) fastaFiles.put(mat.group(1), fastas[i]);
		}
		//check that all chroms are there in the fasta directory
		Iterator<String> it = chromSpecificGeneLines.keySet().iterator();
		while (it.hasNext()){
			String chr = it.next();
			if (fastaFiles.containsKey(chr) == false) Misc.printExit("\nCould not find fasta sequence for "+chr+"! Either delete the offending lines from " +
			"the UCSC gene table or put an appropriate chromosme in the fasta directory.\n");
		}
		//make geneID splices file
		resultsGeneFile = new File (resultsFastaFile.getParentFile(), Misc.removeExtension(resultsFastaFile.getName())+"_Gene2Splices.txt");
		
		//make Ns
		ns = makeNs(sequenceLengthRadius);
	}	
	
	public static String[] makeNs(int number){
		String[] ns = new String[number+1];
		StringBuilder sb = new StringBuilder();
		for (int i=0; i< ns.length; i++){
			ns[i] = sb.toString();
			sb.append("N");
		}
		return ns;
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Make Splice Junction Fasta: Nov 2010                     **\n" +
				"**************************************************************************************\n" +
				
				"\nDEPRECIATED, don't use!  See MakeTranscriptome app!\n"+
				
				"MSJF creates a multi fasta file containing sequences representing all possible linear\n" +
				"splice junctions. The header on each fasta is the chr_endPosExonA_startPosExonB. The\n" +
				"length of sequence collected from each junction is 2x the radius. A word of warning,\n" +
				"be very careful about the coordinate system used in the gene table to define the\n" +
				"start and stop of exons.  UCSC uses interbase and this is assumed in this app. Check\n" +
				"a few of the junctions to be sure correct splices were made. All junction sequences\n" +
				"are from the top/ plus strand of the genome, they are not reverse complemented. Exon\n" +
				"sequence shorter than the radius will be appended with Ns.\n\n" +

				"Options:\n"+
				"-f Fasta file directory, should contain chromosome specific xxx.fasta files.\n" +
				"-u UCSC gene table file, full path. See, http://genome.ucsc.edu/cgi-bin/hgTables\n"+
				"-s Sequence length radius.\n" +
				"-r Results fasta file, full path.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/MakeSpliceJunctionFasta -s 32 \n" +
				"      -f /Genomes/Hg18/Fastas/ -u /Anno/Hg18/ucscKnownGenes.txt -r\n" +
				"      /Genomes/Hg18/Fastas/hg18_32_splices.fasta \n\n" +

		"************************************************************************************\n");

	}
}
