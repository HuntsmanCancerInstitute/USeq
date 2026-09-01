package edu.utah.kohli;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class CopyNumberRiskScoreCalculator {

	//User fields
	private File[] bedFilesToParse = null;
	private String[] genesToCheck = {"AR+_AR-Enh+", "COL22A1+", "MYC+", "NOTCH1+", "PIK3CA+", "PIK3CB+", "TMPRSS2-", "TP53-", "ZBTB16-", "NCOR1-", "NKX3.1-_NKX3-1-"};
	private double[] coxRegCoef = {0.775, 0.904, 0.451, 0.039, 0.199, 3.039, -0.073, 0.358, 0.829, 0.683, -0.371};
	
	//Internal fields
	private HashMap<String, HashSet<String>> dataSetCopyAltGenes = new HashMap<String, HashSet<String>>();
	private TreeSet<String> allObservedGenes = new TreeSet<String>();


	public CopyNumberRiskScoreCalculator(String[] args) {

		try {

			processArgs(args);
			
			IO.pl("Tested Genes: "+Misc.stringArrayToString(genesToCheck, ", "));
			IO.pl("Cox Reg Coef: "+Misc.doubleArrayToString(coxRegCoef, ", ")+"\n");

			parseBedFiles();

			scoreDatasets();

			IO.pl("\nAll observed MG-CNV genes\n\t"+Misc.treeSetToString(allObservedGenes, ", "));

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: running the CopyNumberRiskScoreCalculator\n");
		}
		IO.pl("Done");

	}

	private void scoreDatasets() {
		IO.pl("Dataset\tMG-CNV Risk Score\tObserved Genes");
		//for each dataset
		for (String dataset: dataSetCopyAltGenes.keySet() ) {
			HashSet<String> copyAltGenes = dataSetCopyAltGenes.get(dataset);
			ArrayList<String> observedGenes = new ArrayList<String>();
			IO.p(dataset);
			
			double riskScore = 0;
			//for each gene to check
			for (int i=0; i< genesToCheck.length; i++) {
				//check if it is a split gene like AR or AR-Enh; only score one
				String[] splitGenes = Misc.UNDERSCORE.split(genesToCheck[i]);
				boolean found = false;
				for (int x=0; x<splitGenes.length; x++) {
					String geneToCheck = splitGenes[x];
					if (copyAltGenes.contains(geneToCheck)) {
						if (found==false) riskScore += coxRegCoef[i];
						found = true;
						allObservedGenes.add(geneToCheck);
						observedGenes.add(geneToCheck);
					}
				}
			}
			
			IO.pl("\t"+riskScore+"\t"+Misc.stringArrayListToString(observedGenes, ", "));
			
		}
		
	}

	private void parseBedFiles() {
		//parse bed files
		for (File bedFile: bedFilesToParse) {
			//126006-01-001-C4D15_25KB_Hg38.called.seg.pass.bed
			//  dataset            window
			String[] underscore = Misc.UNDERSCORE.split(bedFile.getName());
			String datasetName = underscore[0];
			HashSet<String> copyGenes = dataSetCopyAltGenes.get(datasetName);
			if (copyGenes == null) {
				copyGenes = new HashSet<String>();
				dataSetCopyAltGenes.put(datasetName, copyGenes);
			}
			addCalls(bedFile, copyGenes);
		}
	}
	
	private void addCalls(File bedFile, HashSet<String> copyGenes) {
		Bed[] bedRegions = Bed.parseFile(bedFile, 0, 0);
		for (Bed b: bedRegions) {
			//numOb=64;lg2Tum=0.2725;lg2Norm=-0.0285;genes=TTTY17C,TTTY17B,TTTY17A
			String[] split = b.getName().split("genes=");
			String genes = split[1];
			for (String g: Misc.COMMA.split(genes)) {
				//exclude any antisense stuff or empty call sets '.'
				if (g.contains("-AS") == false) {
					//strand is used to represent + amp, - del
					copyGenes.add(g+b.getStrand());
				}
			}
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CopyNumberRiskScoreCalculator(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File bedDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedDir = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull bed files
		if (bedDir == null || bedDir.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a bed file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(bedDir, ".bed");
		tot[1] = IO.extractFiles(bedDir,".bed.gz");
		tot[2] = IO.extractFiles(bedDir,".vcf.zip");
		bedFilesToParse = IO.collapseFileArray(tot);
		if (bedFilesToParse == null || bedFilesToParse.length ==0 || bedFilesToParse[0].canRead() == false) {
			Misc.printExit("\nError: cannot find your xxx.bed(.zip/.gz OK) file(s)!\n");
		}

			}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Copy Number Risk Score Calculator: April 2026               **\n" +
				"**************************************************************************************\n" +
				"Merges bed files from the GATK-USeq workflow called with different window sizes. See\n"+
				"https://github.com/HuntsmanCancerInstitute/Workflows/tree/master/Hg38RunnerWorkflows2\n"+
				"/CopyAnalysis and calculates the prostate cancer 11 gene risk score using the\n"+
				"hazard ratios from Huang et.al. 2022 MDPI Cancers.\n"+

				"\nOptions:\n"+
				"-b Directory containing the xxx.called.seg.pass.bed.gz files. The file name will be\n"+
				"     split on '_' to define the dataset name (e.g.\n"+
				"     126006-01-001_25KB_Hg38.called.seg.pass.bed .gz/.zip OK). Datasets with the same\n"+
				"     name are merged.\n"+

				"\nExample: java -jar pathTo/USeq/Apps/CopyNumberRiskScoreCalculator -b PassingBeds/  \n\n"+ 


				"**************************************************************************************\n");

	}
}
