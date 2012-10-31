package util.bio.annotation;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;

import util.bio.parsers.*;

public class ExportIntronicRegions {

	private File ucscTableFile;
	private int minSize = 60;
	private boolean subtractOne = false;	
	private UCSCGeneModelTableReader genes;

	public ExportIntronicRegions (String[] args){
		processArgs(args);
		System.out.println("\nLaunching...\n\tMinimum Size "+minSize+"\n\tSubtract one? "+subtractOne+"\n");
		//load UCSC table and split by chromosome
		int num =0;
		if (subtractOne) num = 1;
		genes = new UCSCGeneModelTableReader(ucscTableFile, num);
		genes.splitByChromosome();
		//for each chromosome, mask, then fetch segment, attempt to build introns
		buildAndPrintIntrons();
		
		System.out.println("\nDone!");

	}

	public void buildAndPrintIntrons(){
		try{
			File intronFile = new File (Misc.removeExtension(ucscTableFile.getCanonicalPath()) + "_Introns.bed");
			PrintWriter out = new PrintWriter ( new FileWriter( intronFile ));
			//for each chromosome
			HashMap chromGenes = genes.getChromSpecificGeneLines();
			Iterator it = chromGenes.keySet().iterator();
			while (it.hasNext()){
				String chrom = (String) it.next();
				System.out.println(chrom);
				UCSCGeneLine[] sub = (UCSCGeneLine[]) chromGenes.get(chrom);
				//make mask
				boolean[] mask = makeExonMask(sub);
				//fetch introns
				String[] introns = fetchIntrons(mask, sub);
				//print introns
				for (int i=0; i< introns.length; i++){
					out.println(chrom+"\t"+introns[i]);
				}
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public String[] fetchIntrons (boolean[] mask, UCSCGeneLine[] genes) {
		//use hash to only get unique introns
		LinkedHashSet introns = new LinkedHashSet((int)genes.length/2);
		for (int i=0; i< genes.length; i++){
			int start = genes[i].getTxStart();
			int length = genes[i].getTxEnd() - start +1;
			if (length <=0) Misc.printExit(genes[i].toString());
			boolean[] sub = new boolean[length];
			System.arraycopy(mask, start, sub, 0, length);
			//fetch blocks of false
			int[][] startStop = ExportIntergenicRegions.fetchFalseBlocks(sub,0,0);
			//add smallest to starts and ends to return to real coordinates and make introns
			if (startStop.length>0){
				for (int j=0; j< startStop.length; j++){
					int size = startStop[j][1]- startStop[j][0] +1;
					if (size >= minSize) introns.add((startStop[j][0]+start)+"\t"+(startStop[j][1]+start));
				}
			}
		}
		//convert hash to array
		Iterator it = introns.iterator();
		String[] in = new String[introns.size()];
		int counter = 0;
		while (it.hasNext()){
			in[counter++] = (String) it.next();
		}
		return in;
	}

	public boolean[] makeExonMask(UCSCGeneLine[] genes){
		//find largest stop
		int largest = -1;
		for (int i=0; i< genes.length; i++){
			ExonIntron[] exons = genes[i].getExons();
			int lastBp = exons[exons.length-1].getEnd();
			if (lastBp > largest) largest = lastBp;
		}
		//make boolean to hold bps
		boolean[] bps = new boolean[largest+1];
		//for each region set booleans to true
		for (int i=0; i< genes.length; i++){
			ExonIntron[] exons = genes[i].getExons();
			//for each exon
			for (int k=0; k< exons.length; k++){
				//for each base
				int end = exons[k].getEnd() + 1;
				int start = exons[k].getStart();
				for (int j= start; j< end; j++){
					bps[j] = true;
				}
			}

		}
		return bps;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ExportIntronicRegions(args);
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
					case 'g': ucscTableFile = new File (args[i+1]); i++; break;
					case 'm': minSize = Integer.parseInt(args[i+1]); i++; break;
					case 's': subtractOne = true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (ucscTableFile == null || ucscTableFile.canRead() == false) Misc.printExit("\nError: cannot find your UCSC table file!\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Export Intronic Regions    June 2007                     **\n" +
				"**************************************************************************************\n" +
				"EIR takes a UCSC Gene table and fetches the most conservative/ smallest intronic\n" +
				"regions. Base coordinates are assumed to be stop inclusive, not interbase.\n\n"+

				"Parameters:\n"+
				"-g Full path file text for the UCSC Gene table.\n"+
				"-m Minimum acceptable intron size, those smaller will be tossed, defaults to 60bp\n"+
				"-s Subtract one from the stop coordinates of your UCSC table to convert from interbase.\n\n"+

				"Example: java -Xmx1000M -jar pathTo/T2/Apps/ExportIntronicRegions -s -m 100 -g\n"+
				"                 /user/Jib/ucscPombe.txt\n\n"+
		"**************************************************************************************\n");		
	}
}

