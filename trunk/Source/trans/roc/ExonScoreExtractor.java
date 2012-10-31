package trans.roc;
import java.io.*;

import trans.main.Interval;
import trans.misc.GrGraph;
import trans.misc.Util;
import util.gen.*;
import java.util.*;

/**Extracts scores from xxx.bar files for particular groups of exons.*/
public class ExonScoreExtractor {

	//fields
	private ControlGene[] controlGenes;
	private HashMap chromNamesHash;
	private String currChrom = "";
	private Gr[][] currGrs = null;
	private int sizeOfOligoMinusOne = 24;
	
	//constructor
	public ExonScoreExtractor(String[] args){
		
		//fetch sorted control genes
		File cgFile = new File(args[0]);
		controlGenes = ControlGene.parseControlGeneFile(cgFile);
		ControlGene.makeExons(controlGenes, sizeOfOligoMinusOne);
		System.out.println(controlGenes.length+" Control genes parsed.");
		
		//fetch chromosome gr directories
		File chromosomeDirectory = new File (args[1]);
		chromNamesHash = IO.fetchNamesAndDirectories(chromosomeDirectory);
		System.out.println(chromNamesHash.size()+ " Chromosome directories found.");

		//process intensities
		extractScores();
		
	}
	
	public void extractScores(){
		//run through control genes extracting scores from each chromsome directory
		for (int i=0; i< controlGenes.length; i++){
			//if needed, load chromosome gr files
			loadChromosomeGrs(controlGenes[i].getChromosome());
			//extract intensities from each file for the particular control gene
			System.out.println("\n"+controlGenes[i].getName()+"\t"+controlGenes[i].getNotes());
			for (int j=0; j< currGrs.length; j++){
				ArrayList intensities = new ArrayList(500);
				//scan positions for intersection and collect their associated intensities 
				for (int k=0; k< currGrs[j].length; k++){
					int grPos = currGrs[j][k].getPosition();
					if (controlGenes[i].containsPoint(grPos)) {
						//System.out.println(currChrom+"\t"+grPos+"\t"+currGrs[j][k].getScore());
						intensities.add(new Float(currGrs[j][k].getScore()));
					}
				}
				float[] scores = Num.arrayListOfFloatToArray(intensities);
				scores = Num.antiLog(scores, 2);
				Arrays.sort(scores);
				//print results: index numOligos median scores
				System.out.println(j+"\t"+scores.length+"\t"+ Num.median(scores));
				
			}
		}
	}
	
	public void loadChromosomeGrs(String chromosome){
		if (currChrom.equals(chromosome) == false) {
			File chromDir = (File)chromNamesHash.get(chromosome);
			if (chromDir == null) Misc.printExit("\nCannot find chrom directory for "+chromosome+ "\n");
			currChrom = chromosome;
			//parse each bar file in the chromosome directory
			File[] barFiles = IO.extractFiles(chromDir, "bar");
			Arrays.sort(barFiles);
			currGrs = new Gr[barFiles.length][];
			System.out.println("\nLoading bar files for "+chromosome);
			for (int j=0; j< barFiles.length; j++){
				System.out.println("\t"+j+"\t"+barFiles[j].getName());
				currGrs[j] = Gr.parseGrGraph(Util.readSimpleGrBarFile(barFiles[j]));
			}
		}
	}
	
	
	
	
	
	
	
	
	
	public static void main (String[] args){
		System.out.println("\nLaunching...\n");
		new ExonScoreExtractor(args);
	}
	

	
	
}
