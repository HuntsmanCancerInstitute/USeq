package util.bio.parsers;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.apps.FindOverlappingGenes;
import util.bio.annotation.ExonIntron;
import util.gen.IO;
import util.gen.Misc;

/**Collapses transcripts into consensus gene merging those with the same text from a UCSC table.*/
public class MergeUCSCGeneTable {

	private File ucscGeneFile;
	private File resultsFile;

	/**Download ensembl gene table from UCSC and replace the ENST text column with the name2.*/
	public  MergeUCSCGeneTable(String[] args) {
		//process user arguments
		processArgs(args);
		
		//load models
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader (ucscGeneFile, 0);
		UCSCGeneLine[] lines = reader.getGeneLines();

		//if displayName is present then need to replace with name
		if (lines[0].getDisplayName() != null){
			for (int i=0; i< lines.length; i++) lines[i].setName(lines[i].getDisplayName());
		}

		HashMap<String,UCSCGeneLine> genes = new HashMap<String,UCSCGeneLine>();


		//merge
		for (int i=0; i< lines.length; i++){
			//does it exist
			if (genes.containsKey(lines[i].getName())){
				//merge exons
				UCSCGeneLine old = genes.get(lines[i].getName());
				ExonIntron[] merged = ExonIntron.merge(old.getExons(), lines[i].getExons());
				old.setExons(merged);
				//reset tx start and stop
				if (old.getTxStart() > lines[i].getTxStart()) old.setTxStart(lines[i].getTxStart());
				if (old.getTxEnd() < lines[i].getTxEnd()) old.setTxEnd(lines[i].getTxEnd());
				//reset cds start stop
				if (old.getCdsStart() > lines[i].getCdsStart()) old.setCdsStart(lines[i].getCdsStart());
				if (old.getCdsEnd() < lines[i].getCdsEnd()) old.setCdsEnd(lines[i].getCdsEnd());
			}
			else genes.put(lines[i].getName(), lines[i]);
		}
		
		System.out.println(lines.length + " transcripts collapsed to "+genes.size()+" genes.\n");

		//print
		try {
			PrintWriter out = new PrintWriter (new FileWriter (resultsFile));
			Iterator<String> it = genes.keySet().iterator();
			while (it.hasNext()){
				out.println(genes.get(it.next()));
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
		new MergeUCSCGeneTable (args);
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
		resultsFile = new File (ucscGeneFile.getParentFile(), name+"_Merged.ucsc");

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Merge UCSC Gene Table: Feb  2013                       **\n" +
				"**************************************************************************************\n" +
				"Merges transcript models that share the same gene name (in column 0). Maximizes exons,\n" +
				"minimizes introns. Assumes interbase coordinates.\n\n" +

				"Options:\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (geneName name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). \n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MergeUCSCGeneTable -u \n" +
				"      /data/zv8EnsemblGenes.ucsc.gz\n\n" +

		"**************************************************************************************\n");

	}


}
