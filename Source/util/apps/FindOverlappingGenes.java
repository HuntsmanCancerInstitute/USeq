package util.apps;
import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.analysis.BaseContext;
import edu.utah.seq.analysis.BisStat;

import trans.tpmap.WindowMaker;
import util.gen.*;
import util.bio.parsers.*;

public class FindOverlappingGenes {

	private File ucscGeneFile;
	private File resultsFile;

	public FindOverlappingGenes(String[] args){
		//process the user arguments
		processArgs(args);

		//load table and split by chromosome
		UCSCGeneModelTableReader genes = new UCSCGeneModelTableReader(ucscGeneFile, 0);
		HashMap<String,UCSCGeneLine[]> lines = genes.splitByChromosomeAndStrand();

		try {
			PrintWriter out = new PrintWriter (new FileWriter (resultsFile));
			out.println("#NamePlusGene\tStart\tStop\tNameMinusGene\tStart\tStop\tChr\tStartInt\tStopInt\tType\tLength");

			//for each plus stranded chromosome
			Iterator<String> it = lines.keySet().iterator();
			while (it.hasNext()){
				String chromStrand = it.next();

				//skip negative strand 
				if (chromStrand.endsWith("-")) continue;
				String chrom = chromStrand.substring(0,chromStrand.length()-1);

				//fetch stranded data
				UCSCGeneLine[] plus = lines.get(chromStrand);
				UCSCGeneLine[] minus = lines.get(chrom+"-");
				if (plus == null || minus == null) {
					System.err.println("WARNING: Missing standed annotation for "+chrom+". Skipping.");
					continue;
				}

				//look for intersection
				for (int j=0; j< plus.length; j++){
					UCSCGeneLine testPlus = plus[j];
					int start = testPlus.getTxStart();
					int stop = testPlus.getTxEnd();
					int size = stop-start;
					for (int k=0; k< minus.length; k++){
						UCSCGeneLine testMinus = minus[k];
						int[] startStop = testMinus.fetchOverlap(start, stop);
						if (startStop != null){
							//make bed line for region of intersection
							//System.out.println(chrs[i]+"\t"+startStop[0]+"\t"+startStop[1]+"\t"+testPlus.getNames("_")+":"+testMinus.getNames("_")+"\t0\t.");
							//NamePlusGene Start Stop NameMinusGene Start Stop StartInt StopInt
							String type;
							int lengthOverlap = startStop[1]-startStop[0];
							if (lengthOverlap == size || lengthOverlap == testMinus.getLength()) type = "Contained";
							else if (testMinus.getTxStart()> start && testMinus.getTxStart() <= stop) type = "Convergent";
							else type = "Divergent";
							out.println(testPlus.getNames("_")+"\t"+start+"\t"+stop+"\t"+testMinus.getNames("_")+"\t"+testMinus.getTxStart()+"\t"+testMinus.getTxEnd()+"\t"+chrom+"\t"+startStop[0]+"\t"+startStop[1]+"\t"+type+"\t"+lengthOverlap);

						}
					}
				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FindOverlappingGenes(args);
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
					case 'u': ucscGeneFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (ucscGeneFile == null || ucscGeneFile.exists() == false) Misc.printErrAndExit("\nError: cannot find your UCSC formatted gene table?\n");
		String name = ucscGeneFile.getName();
		name = name.replace(".gz", "");
		name = name.replace(".txt", "");
		name = name.replace(".ucsc", "");
		resultsFile = new File (ucscGeneFile.getParentFile(), name+"_OverlappingGenes.ucsc");

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Find Overlapping Genes: Oct 2010                       **\n" +
				"**************************************************************************************\n" +
				"Finds overlapping genes that converge, diverge, or contain one another given a UCSC\n" +
				"gene table.\n\n" +

				"Options:\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (name1 name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). NOTE:\n"+
				"       this table should contain only one composite transcript per gene (e.g. Use\n" +
				"       Ensembl genes NOT transcripts. See MergeUCSCGeneTable app.). \n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/FindOverlappingGenes -u \n" +
				"      /data/zv8EnsemblGenes.ucsc.gz\n\n" +

		"**************************************************************************************\n");

	}

}
