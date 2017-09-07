package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Converts a vcf file into a simple bed w/ w/o padding.*/
public class CorrectVCFEnds {

	private File vcfFile;
	private File bedFile;
	private File correctedVcfFile;
	public static Pattern END_SWAP = Pattern.compile("END=(\\d+)");
	public static Pattern MC = Pattern.compile(";*MC=\\w+:\\d+-\\d+");
	
	public CorrectVCFEnds (String[] args) {

		processArgs(args);
		
		//parse bed and add to a hash based on the score
		HashMap<Integer, Bed> indexBed = new HashMap<Integer, Bed>();
		Bed[] bed;
		if (bedFile.exists()) {
			bed = Bed.parseFile(bedFile, 0, 0);
			for (Bed b: bed) {
				if (b.getChromosome().startsWith("chr") == false) b.setChromosome("chr"+b.getChromosome());
				indexBed.put((int)b.getScore(), b);
			}
		}
		
		//walk vcf
		try {
			Gzipper out = new Gzipper(correctedVcfFile);
			BufferedReader in = IO.fetchBufferedReader(vcfFile);
			String[] fields;
			int index = 0;
			String line;
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.trim().length() == 0 || line.startsWith("##INFO=<ID=MC")) continue;
				if (line.startsWith("#")) out.println(line);
				else {
					//look for chr
					if (line.startsWith("chr") == false) line = "chr"+line;
					//any index present?
					Bed b = indexBed.get(index);
					if (b != null){
						fields = Misc.TAB.split(line);
						//check chromosome
						if (fields[0].equals(b.getChromosome()) == false) Misc.printErrAndExit("Error: the chromosome doesn't match for:\n\tbed: "+b.toString()+"\n\tvcf: "+line);
						//check the id
						if (fields[2].equals(b.getName()) == false) Misc.printErrAndExit("Error: the ID/Name doesn't match for:\n\tbed: "+b.toString()+"\n\tvcf: "+line);
						//attempt a match
						Matcher mat = Bed.END.matcher(line);
						if (mat.matches() == false) Misc.printErrAndExit("Error: failed to find and END=xxx for:\n\tbed: "+b.toString()+"\n\tvcf: "+line);
						//swap it
						mat = END_SWAP.matcher(line);
						String swapped = mat.replaceFirst("END="+b.getStop());
						//drop any MC=chr6:152080638-152080754
						mat = MC.matcher(swapped);
						String deleted = mat.replaceFirst("");
						out.println(deleted);
					}
					else out.println(line);
					index++;
				}
			}
			out.close();
			System.out.println("Inserted "+indexBed.size()+" end positions.");
		} catch (Exception e) {
			correctedVcfFile.deleteOnExit();
			e.printStackTrace();
			Misc.printErrAndExit("\nError converting END coordinates!\n");
		} 
	}

	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CorrectVCFEnds(args);
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
					case 'v': vcfFile = new File(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'r': correctedVcfFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

		//pull files
		if (vcfFile == null || vcfFile.exists() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file\n");
		if (correctedVcfFile == null) Misc.printExit("\nError: please provide a vcf results file to write the modified records, xxx.vcf.gz\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Correct VCF Ends: July 2017                       **\n" +
				"**************************************************************************************\n" +
				"Use to correct the END=xxx tags in a Crossmap vcf . Removes any MC tags. Adds chr.\n"+

				"\nRequired Options:\n"+
				"-v Path to the Crossmap vcf file.\n" +
				"-b Path to the VCF2Bed -> Crossmap bed file\n"+
				"-r Path to save the modifed gzipped file. \n"+

				"\nEx: java -jar USeq_XXX/Apps/CorrectVCFEnds -v b38.vcf -b b38.bed -r finalB38.vcf.gz\n\n" +
				"**************************************************************************************\n");

	}

}
