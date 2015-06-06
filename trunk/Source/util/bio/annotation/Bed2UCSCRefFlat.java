package util.bio.annotation;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.Region;
import util.gen.*;
import util.bio.parsers.*;

public class Bed2UCSCRefFlat {

	private File ucscTableFile;	
	private File bedFile;
	private File resultsFile;
	private boolean removeUtrs = true;
	private HashMap<String, ArrayList<Coordinate>> geneName2Coordinate = new HashMap<String, ArrayList<Coordinate>>();
	private HashSet<String> addedRegions = new HashSet<String>();
	private HashMap<String,UCSCGeneLine[]> chrGenes;
	private HashMap<String,Coordinate[]> chrBed;
	private HashSet<String> allChroms;
	private HashMap<String, UCSCGeneLine> geneName2GeneLine = new HashMap<String, UCSCGeneLine>();
	
	public Bed2UCSCRefFlat (String[] args){
		//Parse args
		processArgs(args);
		
		//parse genes
		parseFiles();
		
		//for each bed chrom
		intersect();
		
		//build ucsc table
		buildGeneTable();

	}

	private void buildGeneTable() {
		System.out.println("GeneName\tCoordinates\tNumAssociatedBedRegions");
		Gzipper out;
		try {
			out = new Gzipper(resultsFile);
			//for each gene name
			for (String geneName: geneName2Coordinate.keySet()){
				//fetch gene line and bed regions that int it
				UCSCGeneLine geneLine = geneName2GeneLine.get(geneName);
				ArrayList<Coordinate> bedRegions = geneName2Coordinate.get(geneName);
				UCSCGeneLine newGene = new UCSCGeneLine (geneName, geneLine.getStrand(), bedRegions);
				out.println(newGene);
				System.out.println(geneName+"\t"+ newGene.coordinates()+ "\t"+ bedRegions.size());
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem saving results file!\n");
		} 

	}

	private void intersect() {
		//for each bed chrom
		for (String chr: chrBed.keySet()){
			Coordinate[] bed = chrBed.get(chr);
			UCSCGeneLine[] genes = chrGenes.get(chr);
			if (genes == null) Misc.printErrAndExit("No genes for chrom "+chr);
			
			//for each gene
			for (UCSCGeneLine gene: genes){
				//for each region
				for (Coordinate region: bed){
					//does it intersect?
					if (gene.intersects(region.getStart(), region.getStop())){
						//does it intersect an exon?
						if (gene.intersectsAnExon(region.getStart(), region.getStop())){
							ArrayList<Coordinate> al = geneName2Coordinate.get(gene.getName());
							if (al == null){
								al = new ArrayList<Coordinate>();
								geneName2Coordinate.put(gene.getName(), al);
							}
							al.add(region);
							geneName2GeneLine.put(gene.getName(), gene);
							//check if already added
							String regionName = region.toString();
							if (addedRegions.contains(regionName)) Misc.printErrAndExit("Region has already been associated with a gene \n"+
							regionName+"\n"+gene.getName());
							addedRegions.add(regionName);
							//
							//if (region.getStart() == 2770353) System.out.println(gene.getName());
						}
					}
				}
			}
			
			//look to see if any regions were not added
			for (Coordinate region:bed){
				String regionName = region.toString();
				if (addedRegions.contains(regionName) == false) Misc.printErrAndExit("Region failed to associate with a gene "+regionName);
			}
		}
		
	}


	private void parseFiles() {
		UCSCGeneModelTableReader tr = new UCSCGeneModelTableReader(ucscTableFile, 0);
		if (removeUtrs) tr.trimExonsOfUTRBPs();
		chrGenes = tr.getChromSpecificGeneLines();
		
		//parse bed
		Coordinate[] regions = Coordinate.parseFile(bedFile, 0, 0);
		Arrays.sort(regions);
		chrBed = Coordinate.splitByChromosome(regions);
		
		//collect all chrom
		allChroms = new HashSet<String>();
		allChroms.addAll(chrBed.keySet());
		allChroms.addAll(chrGenes.keySet());
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Bed2UCSCRefFlat(args); 
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
					case 'u': ucscTableFile = new File (args[i+1]); i++; break;
					case 'b': bedFile = new File (args[i+1]); i++; break;
					case 'r': resultsFile = new File (args[i+1]); i++; break;
					case 't': removeUtrs = false; break;
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
		if (bedFile == null || bedFile.canRead() == false) Misc.printExit("\nError: cannot find your bed file!\n");
		if (resultsFile == null) Misc.printExit("\nError: please provide a file path to save the results.\n");
		
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Bed 2 UCSC RefFlat   June 2015                        **\n" +
				"**************************************************************************************\n" +
				"Takes a bed file and a UCSC gene table, intersects them and assigns each bed region\n"+
				"to a gene, then builds a new gene table using the bed region coordinates. Note, each\n"+
				"bed region must intersect only one gene. Modify the input gene table\n"+
				"(MergeUCSCGeneTable and manually trim) based on the errors. Lastly, all bed regions\n"+
				"must be assigned to genes.\n\n"+

				"Parameters:\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (geneName name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). \n"+
				"-b Bed file of regions to intersect with the gene table.\n"+
				"-t Don't remove UTRs if present, from the gene table.\n"+
				"-r Results file.\n"+

				"\nExample: java -Xmx2G -jar pathTo/USeq/Apps/Bed2UCSCRefFlat -u refSeqJun2015.ucsc\n" +
				"      -b targetRegionsFat.bed -r targetRegionsFat.ucsc\n"+

				"**************************************************************************************\n");		
	}		



}

