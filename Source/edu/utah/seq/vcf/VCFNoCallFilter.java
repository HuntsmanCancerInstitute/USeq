package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Parsing for efficient vaast runs.*/
public class VCFNoCallFilter {

	private File[] vcfFiles;
	private int startIndex = 0;
	private int stopIndex = 0;
	private int maxNumNoCalls = 0;
	private float minimumGenotypeQualityGQ = 13f;

	public VCFNoCallFilter (String[] args) {

		processArgs(args);

		System.out.println("\nName\tStartEndNames\tTotal\tPassing\tFraction");
		for (File vcf: vcfFiles) parse(vcf);
	}

	public void parse(File vcf){
		try {
			System.out.print(vcf.getName());
			VCFParser parser = new VCFParser (vcf, false, true, true);
			BufferedReader in = parser.initializeParser();
			
			String[] sampleNames = parser.getSampleNames();
			if (stopIndex==0) stopIndex = sampleNames.length;
			System.out.print("\t"+ sampleNames[startIndex]+"-"+sampleNames[stopIndex-1]);
			
			String name = Misc.removeExtension(vcf.getName());
			Gzipper passOut = new Gzipper(new File (vcf.getParentFile(), name+".passNC.vcf.gz"));
			Gzipper failOut = new Gzipper(new File (vcf.getParentFile(), name+".failNC.vcf.gz"));
			
			//write out header
			for (String h : parser.getStringComments()) {
				passOut.println(h);
				failOut.println(h);
			}
			
			VCFRecord r = null;
			int numRecords = 0;
			int numPass = 0;
			
			while ((r= parser.fetchNext(in)) != null){
				VCFSample[] samples = r.getSample();
				numRecords++;
				
				//check samples for no calls
				int numNoCalls = countNoCallsAndGQ(samples);
				if (numNoCalls > maxNumNoCalls) failOut.println(new String (r.getOriginalRecord()));
				else {
					passOut.println(new String(r.getOriginalRecord()));
					numPass++;
				}
			}
			passOut.close();
			failOut.close();
			double frac = (double)numPass/ (double)numRecords;
			System.out.println("\t"+numRecords+"\t"+numPass+"\t"+frac);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing "+vcf);
		}
	}
	
	
	private int countNoCallsAndGQ(VCFSample[] samples) {
		int numNoCalls = 0;
		for (int i=startIndex; i< stopIndex; i++){
			if (samples[i].isNoCall() || samples[i].getGenotypeQualityGQ() < minimumGenotypeQualityGQ) {
				//System.out.println("NoCall\t"+samples[i].getUnmodifiedSampleString());
				numNoCalls++;
			}
		}
		
		
		return numNoCalls;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFNoCallFilter(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'b': startIndex = Integer.parseInt(args[++i]); break;
					case 'e': stopIndex = Integer.parseInt(args[++i]); break;
					case 'm': maxNumNoCalls = Integer.parseInt(args[++i]); break;
					case 'g': minimumGenotypeQualityGQ = Float.parseFloat(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCFNoCallFilter: April 2015                          **\n" +
				"**************************************************************************************\n" +
				"Parses multi sample VCF records for too many no call or low Genotype Quality records.\n"+
				"Good for removing records where the background is poorly called.\n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-m Maximum number no call or low GQ samples to pass a record\n"+
				"-b Beginning sample index, zero based, included, defaults to 0\n"+
				"-e Ending sample index, not included, defaults to last\n"+
				"-g Minimum Genotype Quality GQ, defaults to 13\n"+

				"\nExample: java -jar pathToUSeq/Apps/VCFNoCallFilter -v /VCFFiles/ -m 5 -b 15 -e 83 \n\n" +
				"**************************************************************************************\n");

	}

}
