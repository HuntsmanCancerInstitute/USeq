package util.bio.seq;
import java.io.*;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.annotation.*;
import util.gen.*;
import util.bio.parsers.*;

/**Takes a UCSC gene table and a directory of xxx.fasta files and replaces the exonic sequence with Ns.*/
public class MaskExonsInFastaFiles {
	
	private File ucscTableFile = null;
	private File fastaFileDirectory = null;
	private File saveDirectory = null;
	private HashMap<String,File> fastaFiles;
	private HashMap<String, UCSCGeneLine[]> chromSpecGenes;

	public  MaskExonsInFastaFiles(String[] args) {
		
		processArgs(args);
		
		//for each chromosome
		
		UCSCGeneLine gene = null;
		try {
		for (String chrom : chromSpecGenes.keySet()){
			//get gene lines
			
			UCSCGeneLine[] sub = (UCSCGeneLine[]) chromSpecGenes.get(chrom);
			
			//parse fasta and get sequence
			MultiFastaParser fastaParser = new MultiFastaParser(fastaFiles.get(chrom));
			char[] seq = fastaParser.getSeqs()[0].toCharArray();
			
			System.out.println(chrom+"\t"+seq.length);
			
			//for each line, get exons
			for (int i=0; i< sub.length; i++){
				gene = sub[i];
				ExonIntron[] exons = gene.getExons();
				//for each exon, mask seq
				for (int j=0; j< exons.length; j++){
					//for each base, note stop is excluded, interbase coordinates
					int start = exons[j].getStart();
					int stop = exons[j].getEnd();
					for (int k=start; k< stop; k++){
						seq[k] = 'N';
					}
				}
			}
			
			//write exon masked fasta
			String name = fastaParser.getNames()[0];
			File masked = new File (saveDirectory, name+".fasta");
			String fastaSeq = ">"+name+"\n"+new String(seq);
			IO.writeString(fastaSeq, masked);
		}
		} catch (Exception e){
			System.out.println("\nProblem processing:\n"+gene);
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MaskExonsInFastaFiles(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fastaFileDirectory = new File (args[++i]); break;
					case 'u': ucscTableFile = new File (args[++i]); break;
					case 's': saveDirectory = new File (args[++i]); break;
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
		if (ucscTableFile == null || ucscTableFile.canRead() == false) Misc.printExit("\nCannot find or read your UCSC gene table?!\n");
		if (fastaFileDirectory == null || fastaFileDirectory.isDirectory() == false) Misc.printExit("\nCannot find your directory of chromosome specific xxx.fasta files?\n");
		if (saveDirectory == null) Misc.printExit("\nPlease indicate a directory to save the results?\n");
		saveDirectory.mkdirs();

		//load UCSC table
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(ucscTableFile,0);
		chromSpecGenes = reader.getChromSpecificGeneLines();

		//load fastaFiles
		fastaFiles = Seq.fetchChromosomeFastaFileHashMap(fastaFileDirectory);

		//check that all chroms are there in the fasta directory
		StringBuilder missingChromosomes = new StringBuilder();
		for (String chr : chromSpecGenes.keySet()){
			if (fastaFiles.containsKey(chr) == false) {
				missingChromosomes.append(" ");
				missingChromosomes.append(chr);
			}
		}
		if (missingChromosomes.length() !=0) Misc.printExit("\nCould not find fasta sequence files for :"+missingChromosomes+"!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Mask Exons In Fasta Files: June 2011                      **\n" +
				"**************************************************************************************\n" +
				"Replaces the exonic sequence with Ns.\n\n"+

				"Options:\n"+
				"-f Fasta file directory, one per chromosome (e.g. chrX.fasta or chrX.fa, .gz/.zip OK)\n"+
				"-u UCSC RefFlat gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (geneName transcriptName chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 ENST00000329454 chr1 + \n" +
				"       16203317 16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 .\n"+
				"-s Save directory, full path.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MaskExonsInFastaFiles -f \n" +
				"      /Genomes/Hg18/Fastas/ -u /Anno/Hg18/ensemblTranscripts.txt.ucsc -s \n" +
				"      /Genomes/Hg18/MaskedFastas/\n\n" +

		"************************************************************************************\n");

	}

}
