package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Fast and loose merge of vcf files with the same sample names.  Essentially hashes header to collapse.  
 * For records with the same chr pos ref alt, keeps the master and adds the idName
 * Probably won't play nice with downstream apps that can't handle mixed FORMAT and INFO records if these are different in the files to merge.
 * @author Nix*/
public class VCFConsensus {
	
	private File primaryVcf;
	private File secondaryVcf;
	private String primaryName;
	private String secondaryName;
	private VCFParser vcf1;
	private VCFParser vcf2;
	private Gzipper out;
	private File mergedVcfFile;

	public VCFConsensus(String[] args){
		try {	
			processArgs(args);

			System.out.println("Loading vcf files...");
			vcf1 = new VCFParser(primaryVcf, true, false, false);
			vcf2 = new VCFParser(secondaryVcf, true, false, false);
			
			System.out.println("Merging headers...");
			String[] mergedHeader = VCFParser.mergeHeaders(new VCFParser[]{vcf1, vcf2});
			if (mergedHeader == null) Misc.printErrAndExit("\nError: hmm something is wrong when merging headers, are the #CHROM lines different?\n");
			
			//create a hash of chromPosRefAlt
			HashMap<String, VCFRecord> secondaryRecords = new HashMap<String, VCFRecord>();
			for (VCFRecord r : vcf2.getVcfRecords()) {
				r.appendId(secondaryName);
				secondaryRecords.put(r.getChrPosRefAlt(false), r);
			}
			
			//for each primary record
			System.out.println("Combining records...");
			ArrayList<VCFRecord> toPrint = new ArrayList<VCFRecord>();
			for (VCFRecord r : vcf1.getVcfRecords()){
				r.appendId(primaryName);
				String prim = r.getChrPosRefAlt(false);
				if (secondaryRecords.containsKey(prim)){
					secondaryRecords.remove(prim);
					r.appendId(secondaryName);
				}
				toPrint.add(r);
			}
			
			//add on remaining secondaries
			toPrint.addAll(secondaryRecords.values());
			
			System.out.println("Sorting and saving vcf records...");
			VCFRecord[] mergedRecords = new VCFRecord[toPrint.size()];
			toPrint.toArray(mergedRecords);
			Arrays.sort(mergedRecords);

			//print header and records
			out = new Gzipper(mergedVcfFile);
			for (String l : mergedHeader) out.println(l);
			for (VCFRecord v: mergedRecords) {
				String[] ori = Misc.TAB.split(v.getOriginalRecord());
				ori[2] = v.getRsNumber();
				out.print(ori[0]);
				for (int i=1; i< ori.length; i++){
					out.print("\t");
					out.print(ori[i]);
				}
				out.println();
			}
			out.close();
			
			//print stats
			System.out.println(vcf1.getVcfRecords().length+"\tRecords in primary "+primaryVcf.getName());
			System.out.println(vcf2.getVcfRecords().length+"\tRecords in secondary "+secondaryVcf.getName());
			System.out.println(vcf2.getVcfRecords().length+"\tRecords in merge "+mergedVcfFile.getName());
			

		} catch (Exception e) {
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
		
		/*here
		
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
		*/
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Consensus : Nov 2016                            **\n" +
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
