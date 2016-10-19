package edu.utah.seq.analysis.ase;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import htsjdk.samtools.*;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.*;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamLayoutForMutation;

/** Application for merging the GeneiASE results and input data into a comprehensive spreadsheet.
 * @author Nix
 * */
public class GeneiASEParser {

	//user defined fields
	private File resultsFile;
	private File dataFile;
	private File outputFile;
	private GeneiASEResult[] gr;
	private GeneiASEData[] gd;
	private TreeMap<String, GeneiASEGene> gGene;

	//constructor
	public GeneiASEParser(String[] args){
		//set fields
		processArgs(args);
		
		parse();
		
		print();

	}

	private void print() {
		try {
			PrintWriter out = new PrintWriter (new FileWriter(outputFile));
			out.println("Gene\t#Alt\t#Ref\tAltFreq\tFDR\t#Calls\tIndividual Calls (chr pos alt ref #alt #ref)");
			for (GeneiASEGene g: gGene.values()) out.println(g);
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void parse() {
		gr = GeneiASEResult.loadResults(resultsFile);
		gd = GeneiASEData.loadData(dataFile);
		
		gGene = new TreeMap<String, GeneiASEGene>();
		
		for (GeneiASEResult r: gr) gGene.put(r.getGene(), new GeneiASEGene(r));
		for (GeneiASEData d: gd){
			GeneiASEGene gene = gGene.get(d.getGene());
			if (gene == null) continue; //Misc.printErrAndExit("\nError: the data has a gene not found in the results "+d.getGene());
			gene.getData().add(d);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new GeneiASEParser(args);
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
					case 'r': resultsFile = new File(args[++i]); break;
					case 'd': dataFile = new File(args[++i]); break;
					case 'o': outputFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (resultsFile == null || resultsFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your GeneiASE results file?\n");
		if (dataFile == null || dataFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your GeneiASE input data table?\n");
		if (outputFile == null ) Misc.printErrAndExit("\nError: can't find your output file?\n");
		

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              GeneiASE Parser:  Sept 2016                         **\n" +
				"**************************************************************************************\n" +
				"Combines the GeneiASE results file with the input data file.\n\n"+

				"Required Options:\n"+
				"-r GeneiASE results output file\n"+
				"-d GeneiASE input data file\n"+
				"-o Output file for the summary spreadsheet\n"+
				
				
				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/AllelicExpressionDetector -b Bam/RPENormal/\n"+
				"-n D002-14,D005-14,D006-14,D009-14 -d GenotypingResults.txt.gz -s SNPMap_Ref2Alt_Int.txt\n"+
				"-r RPENormal -t ~/Anno/b37EnsGenes7Sept2016_Exons.bed.gz\n\n" +

				"**************************************************************************************\n");

	}
}
