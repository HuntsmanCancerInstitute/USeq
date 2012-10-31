package util.bio.seq;
import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.*;
import util.gen.*;
import util.bio.parsers.*;

/**Takes a bed file and a directory of xxx.fasta files and replaces the regions sequence with Ns. Interbase coordinates*/
public class MaskRegionsInFastaFiles {

	private File bedFile = null;
	private File fastaFileDirectory = null;
	private File saveDirectory = null;
	private HashMap<String,File> fastaFiles;
	HashMap<String,RegionScoreText[]> regionsHash;
	boolean maskRegions = true;

	public  MaskRegionsInFastaFiles(String[] args) {

		processArgs(args);

		//for each chromosome
		try {
			for (String chrom : regionsHash.keySet()){

				//get regions 
				RegionScoreText[] sub = regionsHash.get(chrom);

				//parse fasta and get sequence
				MultiFastaParser fastaParser = new MultiFastaParser(fastaFiles.get(chrom));
				char[] seq = fastaParser.getSeqs()[0].toCharArray();

				System.out.println(chrom+"\t"+seq.length);

				//make boolean of masked bases
				boolean[] toMask = new boolean[seq.length];
				
				
				if (maskRegions){
					Arrays.fill(toMask, false);
					//for each region flip toMask to true
					for (int i=0; i< sub.length; i++){
						int start = sub[i].getStart();
						if (start < 0) start = 0;
						int stop = sub[i].getStop();
						if (stop > toMask.length) stop = toMask.length;
						for (int k=start; k< stop; k++) toMask[k] = true;
					}
				}
				else {
					//mask non regions
					Arrays.fill(toMask, true);
					//for each region flip toMask to false
					for (int i=0; i< sub.length; i++){
						int start = sub[i].getStart();
						if (start < 0) start = 0;
						int stop = sub[i].getStop();
						if (stop > toMask.length) stop = toMask.length;
						for (int k=start; k< stop; k++) toMask[k] = false;
					}
				}
				
				//mask sequence
				for (int i=0; i< toMask.length; i++) {
					if (toMask[i]) {
						seq[i] = 'N';
					}
				}
				
				//write masked fasta
				String name = fastaParser.getNames()[0];
				File masked = new File (saveDirectory, name+".fasta");
				String fastaSeq = ">"+name+"\n"+new String(seq);
				IO.writeString(fastaSeq, masked);
			}
		} catch (Exception e){
			System.out.println();
			e.printStackTrace();
			System.exit(1);
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MaskRegionsInFastaFiles(args);
	}		

	/**This method will process each argument and assign new variables*/
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
					case 'b': bedFile = new File (args[++i]); break;
					case 's': saveDirectory = new File (args[++i]); break;
					case 'r': maskRegions = false; break;
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
		if (bedFile == null || bedFile.canRead() == false) Misc.printExit("\nCannot find or read your bed file?!\n");
		if (fastaFileDirectory == null || fastaFileDirectory.isDirectory() == false) Misc.printExit("\nCannot find your directory of chromosome specific xxx.fasta files?\n");
		if (saveDirectory == null) Misc.printExit("\nPlease indicate a directory to save the results?\n");
		saveDirectory.mkdirs();

		//load regions, no strand
		regionsHash = Bed.parseBedFile(bedFile, true);

		//load fastaFiles
		fastaFiles = Seq.fetchChromosomeFastaFileHashMap(fastaFileDirectory);

		//check that all chroms are there in the fasta directory
		StringBuilder missingChromosomes = new StringBuilder();
		for (String chr : regionsHash.keySet()){
			if (fastaFiles.containsKey(chr) == false) {
				missingChromosomes.append(" ");
				missingChromosomes.append(chr);
			}
		}
		if (missingChromosomes.length() !=0) Misc.printExit("\nAborting! Could not find fasta sequence files for :"+missingChromosomes+"!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Mask Regions In Fasta Files: Dec 2011                      **\n" +
				"**************************************************************************************\n" +
				"Replaces the region (or non region) sequence with Ns. Interbase coordinates.\n\n"+

				"Options:\n"+
				"-f Fasta file directory, one per chromosome (e.g. chrX.fasta or chrX.fa, .gz/.zip OK)\n"+
				"-b Bed file of regions to mask.\n"+
				"-s Save directory, full path.\n"+
				"-r Mask sequence not in regions, reverse mask.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MaskRegionsInFastaFiles -f \n" +
				"      /Genomes/Hg18/Fastas/ -b /Anno/Hg18/badRegions.bed -s \n" +
				"      /Genomes/Hg18/MaskedFastas/\n\n" +

		"************************************************************************************\n");

	}

}
