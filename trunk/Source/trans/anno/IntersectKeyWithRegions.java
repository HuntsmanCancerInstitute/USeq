package trans.anno;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;


/**Takes a key from a spike in experiment and sorted lists of genomicRegions to calculate TPR, FPR, and FDR for every threshold*/
public class IntersectKeyWithRegions {

	//fields
	private File keyFile;
	private File[] regionsFiles;
	private GenomicRegion[] key;
	private GenomicRegion[][] genomicRegions;
	private File resultsFile = null;
	private HashSet hitKeys = new HashSet();
	private StringBuffer multipleHits = new StringBuffer();
	private int maxGap = -1;
	private double numberInKey;
	private double maxFDR = 1;//0.2;
	private boolean subtractOne = false;

	//constructor
	public IntersectKeyWithRegions(String[] args) {

		processArgs(args);

		parseRegionFiles();

		numberInKey = key.length;
		System.out.println("\nNumber in key -> "+numberInKey);

		//for each region file
		for (int i=0; i< regionsFiles.length; i++){
			System.out.println("\nProcessing "+regionsFiles[i].getName());
			resultsFile = new File(regionsFiles[i].getParentFile(), Misc.removeExtension(regionsFiles[i].getName())+"IntWithKey.xls");
			calculateAndPrintStats(genomicRegions[i]);
			//clear hash
			hitKeys.clear();
			//Any multiple hits?
			if (multipleHits.length() !=0 ){
				System.out.println("\nWARNING: Multiple hits to same key were found and excluded!\n"+multipleHits);
				multipleHits = new StringBuffer();
			}
		}
	}

	public void calculateAndPrintStats(GenomicRegion[] r){	
		try {
			PrintWriter out = new PrintWriter ( new FileWriter(resultsFile));
			out.println("NumTP\tNumFP\tFPR\tTPR\tFDR\tChr\tStart\tStop\tEtc");
			//for each region, check if it intersects with the key and calculate stats
			double numTP =0;
			double total = 0;
			for (int j=0; j< r.length; j++){
				//check for intersection
				int hitStatus = intersectsKey(r[j]);
				if (hitStatus == 1) {
					numTP++;
				}
				else if (hitStatus == -1) {
					continue;
				}

				//calculate stats
				double tpr = numTP/numberInKey;
				total++;
				double numFP = total-numTP;

				double fpr = numFP/numberInKey;
				double fdr = numFP/(numFP+numTP); 
				
				//index numTP numFP fpr tpr fdr
				out.println((int)numTP+"\t"+(int)numFP+"\t"+fpr+"\t"+tpr+"\t"+fdr+"\t"+r[j]);
				
				/*
				//System.out.println("non "+total+"\t"+fpr+"\t"+tpr+"\t"+fdr);
				//print FPR TPR FDR ?
				if ( fdr !=0 && fdr <= maxFDR) {
					if (firstLine) {
						out.println(line);
						firstLine = false;
					}
					out.println((j+1)+"\t"+fpr+"\t"+tpr+"\t"+fdr);
				}
				if (firstLine) line = (j+1)+"\t"+fpr+"\t"+tpr+"\t"+fdr;
				*/
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/**Returns 1 if bpInt <= maxGap and the particular key hasn't been hit, 0 not in key, -1 hits a key that 
	 * has already been counted.*/
	public int intersectsKey (GenomicRegion genomicRegion){
		String chromosome = genomicRegion.getChromosome();
		for (int i=0; i< numberInKey; i++){
			//check chromosome
			if (key[i].getChromosome().equals(chromosome)){
				int bpInt = key[i].bpIntersectionSameChromosome(genomicRegion);
				if (bpInt <= maxGap) {
					//hash it already been hit?
					String keyCoor = key[i].toString();
					if (hitKeys.contains(keyCoor)) {
						multipleHits.append("\tkey: "+keyCoor.replaceAll("\\t"," ")+"  pick: "+genomicRegion.toString().replaceAll("\\t"," ")+"\n");
						return -1;
					}
					hitKeys.add(keyCoor);
					return 1;
				}
			}
		}
		return 0;
	}

	public void parseRegionFiles(){
		key = GenomicRegion.parseRegions(keyFile);
		if (subtractOne){
			for (GenomicRegion gr: key){
				gr.setEnd(gr.getEnd()-1);
			}
		}
		genomicRegions = new GenomicRegion[regionsFiles.length][];
		for (int i=0; i< regionsFiles.length; i++){
			genomicRegions[i] = GenomicRegion.parseRegions(regionsFiles[i]);
			if (subtractOne){
				for (GenomicRegion gr: genomicRegions[i]){
					gr.setEnd(gr.getEnd()-1);
				}
			}
		}
	}

	public static void main(String[] args){
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new IntersectKeyWithRegions(args);
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
					case 'k': keyFile = new File(args[i+1]); i++; break;
					case 'r': regionsFiles = IO.extractFiles(new File (args[i+1])); i++; break;
					case 'g': maxGap = Integer.parseInt(args[i+1]); i++; break;
					case 's': subtractOne = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//look for required parameters
		if (keyFile == null || keyFile.canRead() == false )Misc.printExit("\nCannot find your key file? "+keyFile);
		if (regionsFiles == null || regionsFiles.length == 0 || regionsFiles[0].canRead() == false){
			Misc.printExit("\nCannot find your region files? ");
		}


	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Intersect Key With Regions: July 2012                      **\n" +
				"**************************************************************************************\n" +
				"IR intersects lists of genomicRegions (chrom start stop(inclusive)) with a key, assumes the\n" +
				"lists are sorted from most confident to least confident. Multiple hits to the same key\n" +
				"region are ignored.\n\n"+

				"-k Full path file text for the key genomicRegions file, tab delimited (chr start\n" +
				"      stop(inclusive)).\n"+
				"-r Full path file text or directory containing your region files to score.\n"+
				"-g Max gap, defaults to -1. A max gap of 0 = genomicRegions must abut, negative values force\n" +
				"      overlap (ie -1= 1bp overlap, be careful not to exceed the length of the smaller\n" +
				"      region), positive values enable gaps (ie 1=1bp gap).\n"+
				"-s Subtract 1 from end coordinates.  Use for interbase.\n"+

				"\nExample: java -Xmx1500M -jar pathTo/Apps/IntersectKeyWithRegions -k /data/key.txt\n" +
				"      -r /data/HitLists/ \n" +

				"\n" +
		"**************************************************************************************\n");		
	}

}
