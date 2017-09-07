package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Fast and loose merge of vcf files with the same sample names.  Essentially hashes header to collapse.  
 * Probably won't play nice with downstream apps that can't handle mixed FORMAT and INFO records if these are different in the files to merge.
 * @author Nix*/
public class VCFMerger {
	
	private File[] vcfFiles;
	private VCFParser[] vcfParsers;
	private Gzipper out;
	private File mergedVcfFile;

	public VCFMerger(String[] args){
		try {	
			processArgs(args);

			//load and merge headers
			System.out.println("Loading vcf files:");
			vcfParsers = new VCFParser[vcfFiles.length];
			String[] sampleNames = null;
			int numberRecords = 0;
			for (int i=0; i< vcfFiles.length; i++){
				System.out.print("\t"+ vcfFiles[i].getName());
				vcfParsers[i] = new VCFParser(vcfFiles[i], true, false, false);
				//check that sample names are the same
				if (sampleNames == null) sampleNames = vcfParsers[i].getSampleNames();
				else {
					String[] newNames = vcfParsers[i].getSampleNames();
					if (sampleNames.length != newNames.length) Misc.printErrAndExit("\nError: different number of samples!\n");
					for (int j=0; j< sampleNames.length; j++){
						if (sampleNames[j].equals(newNames[j]) == false) Misc.printErrAndExit("\nError: sample names differ!\n");
					}
				}
				//increment num records
				int num = vcfParsers[i].getVcfRecords().length;
				numberRecords+= num;
				System.out.println("\t"+num);
			}
			
			System.out.println("\nMerging headers, skips those with the same ID...");
			String[] mergedHeader = VCFParser.mergeHeaders(vcfParsers, false);
			if (mergedHeader == null) Misc.printErrAndExit("\nError: hmm something is wrong when merging headers, are the #CHROM lines different?\n");
			
			//merge records
			System.out.println("Sorting vcf records...");
			VCFRecord[] mergedRecords = new VCFRecord[numberRecords];
			int index = 0;
			for (int i=0; i< vcfParsers.length; i++){
				VCFRecord[] vr = vcfParsers[i].getVcfRecords();
				for (VCFRecord v : vr) mergedRecords[index++] = v;
			}
			Arrays.sort(mergedRecords);

			//print header and records
			out = new Gzipper(mergedVcfFile);
			for (String l : mergedHeader) out.println(l);
			for (VCFRecord v: mergedRecords) out.println(v.getOriginalRecord());

			out.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}



	

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFMerger(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'o': mergedVcfFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please indicate a vcf file to filter.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");

		//final file
		if (mergedVcfFile == null) mergedVcfFile = new File (vcfFiles[0].getParentFile(), "merged.vcf.gz");
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Merger : Feb 2016                               **\n" +
				"**************************************************************************************\n" +
				"Merges VCF files with the same samples. Collapses the headers with a simple hash. Will\n"+
				"not work well with downstream apps that cannot process mixed INFO and FORMAT records.\n" +

				"\nRequired:\n"+
				"-v Full path to a vcf file (xxx.vcf(.gz/.zip OK)) or directory containing such. Note,\n"+
				"       Java often fails to parse tabix compressed vcf files.  Best to uncompress.\n\n"+
								
				"Optional:\n" +
				"-o Full path to an output vcf file, defaults to merged.vcf.gz in parent -v dir.\n" +

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFMerger -v /CancerSamples/\n\n"+

		"**************************************************************************************\n");

	}

}
