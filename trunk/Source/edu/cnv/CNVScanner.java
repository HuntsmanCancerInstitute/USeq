package edu.cnv;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import trans.misc.*;
import edu.utah.seq.parsers.*;
import util.gen.*;
import trans.roc.*;
import trans.anno.*;
import util.bio.annotation.*;

public class CNVScanner {
	
	//fields
	private File barDirectory;
	//int numberProbes = 2;
	//float minimumScore = 0.2630344f; //1.5x
	//float maximumScore = -0.2630344f;
	private int numberProbes = 3;
	private float minimumScore = 2.8f;
	private float maximumScore = 1.2f;
	private File resultsDirectory;
	private ArrayList<CNV> cnvsAL = new ArrayList<CNV>();
	private CNV[] cnvs;
	private File ucscGeneFile;
	private int neighborhood = 10000;
	private boolean antiLog2 = false;
	
	
	public CNVScanner(String[] args){
		//parse args
		processArgs(args);
		
		//scan chromosomes
		scanChromosomes();
		
		//make coordinates
		if (cnvs.length == 0) Misc.printExit("No CNVs found!");
		System.out.println(cnvs.length+" CNVs found");
		Coordinate[] coor = new Coordinate[cnvs.length];
		for (int i=0; i< cnvs.length; i++) coor[i] = new Coordinate(cnvs[i].getGrGraph().getChromosome(), cnvs[i].getStart(), cnvs[i].getStop());
		
		//find neighboring genes
		FindNeighboringGenes fng = new FindNeighboringGenes(ucscGeneFile, coor, neighborhood, true);
		String[] effGenes = fng.fetchGeneNames();
		
		//print
		try{
			PrintWriter spreadSheet = new PrintWriter (new FileWriter (new File(resultsDirectory,barDirectory.getName()+".xls")));
			spreadSheet.println("Chrom\tStart\tStop\tLength\t#Probes\tMedianCN\tProbeBasePositions\tProbeCNs\tIntersectedGenes\tParams: Min#Probes="+numberProbes+" MinCN="+minimumScore+" MaxCN="+maximumScore+" Neighborhood="+neighborhood);
			PrintWriter bed = new PrintWriter (new FileWriter (new File(resultsDirectory,barDirectory.getName()+".bed")));
			bed.println("#Chrom\tStart\tStop\tParams: Min#Probes="+numberProbes+" MinCN="+minimumScore+" MaxCN="+maximumScore+" Neighborhood="+neighborhood);

			for (int i=0; i< cnvs.length; i++) {
				spreadSheet.println(cnvs[i]+"\t"+effGenes[i]);
				bed.println(cnvs[i].toStringBed(""+(i+1)));
			}
			spreadSheet.close();
			bed.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		
		//save CNVs
		File cnvFile = new File(resultsDirectory, barDirectory.getName()+".cnv");
		IO.saveObject(cnvFile, cnvs);
	}
	
	public void scanChromosomes(){
		File[] barFiles = fetchBarFiles(barDirectory);
		for (int i=0; i< barFiles.length; i++){
			scanChromosome(barFiles[i]);
		}
	}
	
	public void scanChromosome(File barFile){
		BarParser bp = new BarParser();
		bp.readBarFile(barFile, true);
		int[] positions = bp.getBasePositions();
		float[] values = bp.getValues();
		if (antiLog2) values = Num.antiLog(values, 2);
		String chromosome = bp.getChromosome();
		boolean in = false;
		CNV currentCNV = null;
		//scan for increased
		for (int i=0; i< positions.length; i++){
			//does value meet threshold?
			if (values[i] >= minimumScore){
				if (in == false){
					//start new
					currentCNV = new CNV(chromosome);
					in = true;
				}
				//add value
				currentCNV.grs.add(new Gr(positions[i],values[i]));
			}
			else {
				//doesn't make it
				in = false;
				//add currentCNV?
				if (currentCNV !=null && currentCNV.grs.size() >= numberProbes) cnvsAL.add(currentCNV);
				//reset current
				currentCNV = new CNV(chromosome);
			}
		}
		//close 
		in = false;
		if (currentCNV !=null && currentCNV.grs.size() >= numberProbes) cnvsAL.add(currentCNV);
		currentCNV = new CNV(chromosome);
		
		//scan for decreased
		for (int i=0; i< positions.length; i++){
			//does value meet threshold?
			if (values[i] <= maximumScore){
				if (in == false){
					//start new
					currentCNV = new CNV(chromosome);
					in = true;
				}
				//add value
				currentCNV.grs.add(new Gr(positions[i],values[i]));
			}
			else {
				//doesn't make it
				in = false;
				//add currentCNV?
				if (currentCNV !=null && currentCNV.grs.size() >= numberProbes) cnvsAL.add(currentCNV);
				//reset current
				currentCNV = new CNV(chromosome);
			}
		}
		//close
		if (currentCNV !=null && currentCNV.grs.size() >= numberProbes) cnvsAL.add(currentCNV);
		
		//convert
		cnvs = new CNV[cnvsAL.size()];
		cnvsAL.toArray(cnvs);
		//load
		for (int i=0; i< cnvs.length; i++) cnvs[i].loadGrGraph();
		//sort
		Arrays.sort(cnvs, new ComparatorCNVScore());
	}
	

	
	public static File[] fetchBarFiles(File barDirectory){
		//fetch files
		File[] files = IO.extractFiles(barDirectory, ".bar");
		if (files == null || files.length == 0) {
			files = IO.extractFiles(barDirectory, ".bar.zip");
			if (files == null || files.length == 0) {
				System.out.println("\t\tError: no xxx.bar or xxx.bar.zip files found in -> "+barDirectory+"\n");
				return null;
			}	
		}
		return files;
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CNVScanner(args);
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
					case 'b': barDirectory = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'g': ucscGeneFile = new File(args[++i]); break;
					case 'n': numberProbes = Integer.parseInt(args[++i]); break;
					case 'm': minimumScore = Float.parseFloat(args[++i]); break;
					case 'x': maximumScore = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		System.out.println("\tBar Dir:\t"+barDirectory);
		System.out.println("\tResults Dir:\t"+resultsDirectory);
		System.out.println("\tUCSC Gene Tbl:\t"+ucscGeneFile);
		System.out.println("\tMin # probes:\t"+numberProbes);
		System.out.println("\tMin Score:\t"+minimumScore);
		System.out.println("\tMax Score:\t"+maximumScore);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            CNV Scanner: January 2009                             **\n" +
				"**************************************************************************************\n" +
				"Takes a directory containing chromosome specific xxx.bar files and identifies CNV loci\n" +
				"that contain the minimum number of adjacent probes that exceed the set scores. These\n" +
				"can then be intersected using the IntersectCNVs application.\n\n" +

				"Options:\n"+
				"-r Results directory\n"+
				"-b Bar file directory\n"+
				"-n Minimum number of probes, defaults to 3\n"+
				"-m Minimum score for enriched CNVs, defaults to 2.8\n"+
				"-x Maximum score for reduced CNVs, defaults to 1.2\n"+
				"-g UCSC RefFlat gene table file, See http://genome.ucsc.edu/cgi-bin/hgTables\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/CNVAn/Apps/CNVScanner -r /Leuk/Results -b\n" +
				"      /Leuk/MIPBarData -n 4 -m 3 -x 1 -g /Anno/Hg18/refSeq.txt \n\n" +

		"**************************************************************************************\n");

	}
}
