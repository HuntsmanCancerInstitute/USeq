package edu.utah.seq.cnv;

import java.io.BufferedReader;
import java.io.File;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class GatkCopyRatioFilter {

	//user defined fields
	private File[] bedFiles = null;
	private File resultsDirectory = null;
	private double minAbsLg2TumorCopyRatio = 0.585;  //1.5x
	private double maxAbsLg2NormalCopyRatio = 0.202; //1.15x
	private double minCopyRatioMeanTNRatios = 0.322; //1.25x
	private File geneNamesToKeep = null;
	
	//internal fields
	private HashSet<String> geneNames = null;
	private Gzipper impactedGenes = null;
	
	public GatkCopyRatioFilter(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);
		
		try {
			
			//load gene names to keep if present
			if (geneNamesToKeep != null) geneNames = IO.loadFileIntoHashSet(geneNamesToKeep);
			
			//start up impactedGenes.txt.gz
			impactedGenes = new Gzipper( new File(resultsDirectory, "impactedGenes.txt.gz"));
			
			//for each bed file
			IO.pl("File\tPassing\tFailing\tImpactedGenes");
			for (File b: bedFiles) parse(b);
			
			
		} catch (Exception e) {
			IO.el("\nERROR running CopyRatioAggregator\n");
			e.printStackTrace();
		} finally {
			if (impactedGenes != null) impactedGenes.closeNoException();
		}

		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}
	


	private void parse(File b) throws Exception{
		String name = Misc.removeExtension(b.getName());
		IO.p(name);
		
		BufferedReader in = IO.fetchBufferedReader(b);
		Gzipper out = new Gzipper(new File(resultsDirectory, name+".filt.bed.gz"));
		String line = null;
		String[] fields = null;
		String[] info = null;
		int numFail = 0;
		int numPass = 0;
		HashSet<String> igs = new HashSet<String>();
		
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.length()==0) continue;
			if (line.startsWith("#")) out.println(line);
			else {
				fields = Misc.TAB.split(line);
				
				//does it fail the minCopyRatioMeanTNRatios
				double rto = Double.parseDouble(fields[4]);
				if (Math.abs(rto) < minCopyRatioMeanTNRatios) {
					numFail++;
					continue;
				}
				//parse the info 
				//numOb=23;lg2Tum=-1.553;lg2Norm=0.4019;genes=CCDC102B,LOC101927430,LOC643542,DSEL,TMX3
				info = Misc.SEMI_COLON.split(fields[3]);
				double lg2Tum = Double.NaN;
				double lg2Norm = Double.NaN;
				String[] genes = null;
				for (String i: info) {
					if (i.startsWith("lg2Tum=")) lg2Tum = Double.parseDouble(i.substring(7));
					else if (i.startsWith("lg2Norm=")) lg2Norm = Double.parseDouble(i.substring(8));
					else if (i.startsWith("genes=")) genes = Misc.COMMA.split(i.substring(6));
				}
				if (Double.isNaN(lg2Tum) || Double.isNaN(lg2Norm)) throw new Exception("\nERROR: Is this a xxx.called.seg.pass.bed.gz file from the GatkCalledSegmentAnnotator app? Failed to find lg2Tum= or lg2Norm= in "+line);
				
				//does it pass the tumor and normal thresholds
				if (Math.abs(lg2Tum) < minAbsLg2TumorCopyRatio || Math.abs(lg2Norm) > maxAbsLg2NormalCopyRatio) numFail++;
				else {
					numPass++;
					out.println(line);
					//any genes to parse?
					if (genes!= null && genes.length !=0) {
						//filter for genes to keep?
						if (geneNames == null) for (String g: genes) igs.add(g);
						else {
							for (String g: genes) {
								if (geneNames.contains(g)) igs.add(g);
							}
						}
					}
				}
			}
		}
		//print results for the file
		IO.pl("\t"+ numPass+"\t"+numFail+"\t"+igs.size());
		
		//save the impacted genes
		impactedGenes.println(name+ Misc.hashSetToString(igs, "\t"));
		
		//close the IO
		out.close();
		in.close();
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new GatkCopyRatioFilter(args);
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
					case 'b': bedFiles = IO.extractFiles(new File(args[++i]), ".bed.gz"); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'g': geneNamesToKeep = new File(args[++i]); break;
					case 'm': minCopyRatioMeanTNRatios = Double.parseDouble(args[++i]); break;
					case 'c': minAbsLg2TumorCopyRatio = Double.parseDouble(args[++i]); break;
					case 'x': maxAbsLg2NormalCopyRatio = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		
		//check results dir
		if (resultsDirectory == null) Misc.printErrAndExit("\nError: provide a directory to save results.\n");
		resultsDirectory.mkdirs();
		
		if (bedFiles == null || bedFiles.length==0) Misc.printErrAndExit("\nError: please provide a directory containing xxx.bed.gz files or a single file for parsing.\n");

		printOptions();
		
	}	
	


	private void printOptions() {
		IO.pl("Settings:");
		IO.pl("\t-b Bed file #\t"+bedFiles.length);
		IO.pl("\t-r ResultsDirectory\t"+resultsDirectory);
		IO.pl("\t-g Genes to keep in impactedGenes.txt file\t"+geneNamesToKeep);
		IO.pl("\t-m MinimumAbsLg2TNRatioOfCopyRatios\t"+minCopyRatioMeanTNRatios);
		IO.pl("\t-c MinimumAbsLg2TumorCopyRatio\t"+minAbsLg2TumorCopyRatio);
		IO.pl("\t-x MaximumAbsLg2NormalCopyRatio\t"+maxAbsLg2NormalCopyRatio);
		IO.pl();
	}

	public static void printDocs(){
		
		
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Gatk Copy Ratio Filter: Nov 2023                          **\n" +
				"**************************************************************************************\n" +
				"Filters and aggregates the copy ratio calls from the GatkCalledSegmentAnnotator. Use\n"+
				"it to refine the calls in batch from the GCSA.\n"+

				"\nOptions:\n"+
				"-b Directory containing the xxx.bed.gz files from the GatkCalledSegmentAnnotator\n"+
				"-r Directory to save the filtered bed files and the impacted genes aggregate file\n"+
				"-g File containing gene names to keep in the impacted gene aggregate file. One gene\n"+
				"      per line.\n"+
				"-c Minimum absolute tumor log2 copy ratio, defaults to 0.585 (1.5x)\n"+
				"-x Maximum absolute normal log2 copy ratio, defaults to 0.202 (1.15x)\n"+
				"-m Minimum absolute log2 TN ratio of copy ratios, defaults to 0.322 (1.25x)\n"+
				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/GatkCopyRatioFilter -r FilteredResults/\n"+
				"       -b PassingBedFiles -g ~/UCSC/stndHg38RefSeqGenes.txt -m 0.585 \n\n" +

				"**************************************************************************************\n");
	}
}
