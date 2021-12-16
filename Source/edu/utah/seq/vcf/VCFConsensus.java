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
	private int primaryNameIndex = 0;
	private int secondaryNameIndex = 0;
	private VCFParser primaryVcfParser;
	private VCFParser secondaryVcfParser;
	private Gzipper out;
	private File mergedVcfFile;
	private int numInCommon = 0;
	private boolean verbose = true;
	private boolean tossSampleInfo = false;
	private boolean useFirstChrom = false;

	public VCFConsensus(String[] args){
		try {	
			processArgs(args);
			doWork();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public VCFConsensus(File primaryVcf, File secondaryVcf, File mergedVcfFile, boolean tossSampleInfo, boolean verbose) throws Exception{
		this.primaryVcf = primaryVcf;
		this.secondaryVcf = secondaryVcf;
		this.mergedVcfFile = mergedVcfFile;
		this.tossSampleInfo = tossSampleInfo;
		this.verbose = verbose;
		doWork();
	}
	
	public void doWork() throws Exception {
		if (verbose) System.out.println("Loading vcf files...");
		primaryVcfParser = new VCFParser(primaryVcf, true, false, false);
		secondaryVcfParser = new VCFParser(secondaryVcf, true, false, false);
		
		if (verbose) System.out.println("Merging headers...");
		String[] mergedHeader = VCFParser.mergeHeaders(new VCFParser[]{primaryVcfParser, secondaryVcfParser}, tossSampleInfo, useFirstChrom);
		if (mergedHeader == null) Misc.printErrAndExit("\nError: hmm something is wrong when merging headers, are the #CHROM lines different?\n");
		
		//create a hash of chromPosRefAlt
		HashMap<String, VCFRecord> secondaryRecords = new HashMap<String, VCFRecord>();
		for (VCFRecord r : secondaryVcfParser.getVcfRecords()) {
			if (secondaryName!=null) r.setRsNumber(secondaryName+ "_"+ (secondaryNameIndex++));
			secondaryRecords.put(r.getChrPosRefAlt(false), r);
		}
		
		//for each primary record
		if (verbose) System.out.println("Combining records...");
		ArrayList<VCFRecord> toPrint = new ArrayList<VCFRecord>();
		for (VCFRecord r : primaryVcfParser.getVcfRecords()){
			if (primaryName!= null) r.setRsNumber(primaryName+ "_"+ (primaryNameIndex++));
			String prim = r.getChrPosRefAlt(false);
			if (secondaryRecords.containsKey(prim)){
				String sec = secondaryRecords.get(prim).getRsNumber();
				if (sec == null || sec.length()==0){
					if (secondaryName == null) Misc.printErrAndExit("\nA secondary record was found with no ID info. Please provide a secondaryName to append to the merged output, see "+secondaryRecords.get(prim).getOriginalRecord());
				}
				//does the primary have a name?
				if (r.getRsNumber() == null || r.getRsNumber().length() == 0) Misc.printErrAndExit("\nA primary record was found with no ID info. Please provide a primaryName to append to the merged output, see "+r.getOriginalRecord());
				//add the secondary id to the primary
				r.appendId(sec);
				r.appendFilter(secondaryRecords.get(prim).getFilter());
				secondaryRecords.remove(prim);
				numInCommon++;
			}
			toPrint.add(r);
		}
		
		//add on remaining secondaries
		toPrint.addAll(secondaryRecords.values());
		
		if (verbose) System.out.println("Sorting and saving vcf records...");
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
			int numToPrint = 0;
			if (tossSampleInfo) numToPrint = 8;
			else numToPrint = ori.length;
			for (int i=1; i< numToPrint; i++){
				out.print("\t");
				out.print(ori[i]);
			}
			out.println();
		}
		out.close();
		
		//print stats
		if (verbose) {
			System.out.println(primaryVcfParser.getVcfRecords().length+"\tRecords in primary "+primaryVcf.getName());
			System.out.println(secondaryVcfParser.getVcfRecords().length+"\tRecords in secondary "+secondaryVcf.getName());
			System.out.println(numInCommon+"\tRecords in common");
			System.out.println(mergedRecords.length+"\tRecords in merge "+mergedVcfFile.getName());
		}
		
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFConsensus(args);
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
					case 'p': primaryVcf = new File(args[++i]); break;
					case 's': secondaryVcf = new File(args[++i]); break;
					case 'o': mergedVcfFile = new File(args[++i]); break;
					case 'q': primaryName = args[++i]; break;
					case 't': secondaryName = args[++i]; break;
					case 'u': useFirstChrom = true; break;
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
		if (primaryVcf == null || secondaryVcf == null) Misc.printExit("\nError: please provide both a primary and secondary vcf file to merge.\n");
		
		//final file
		if (mergedVcfFile == null) Misc.printExit("\nError: please provide a file name to write the merged vcf records.\n");
				
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Consensus : Dec 2021                            **\n" +
				"**************************************************************************************\n" +
				"Merges VCF files with the approx same #CHROM line. Primary records with the same \n"+
				"chrPosRefAlt as a secondary are saved after appending the ID and FILTER, the secondary\n"+
				"is dropped. Headers are joined keeping the primary header line when the same. Run\n" +
				"iteratively with multiple VCF files you'd like to merge.  Good for combining multiple\n"+
				"variant callers run on the same sample. The ID field lists which callers found each\n"+
				"variant.\n"+

				"\nRequired:\n"+
				"-p Path to a primary vcf file (xxx.vcf(.gz/.zip OK)) to merge.\n"+
				"-s Path to a secondary vcf file (xxx.vcf(.gz/.zip OK)) to merge.\n"+
				"-o Path to an output xxx.vcf.gz file.\n" +
								
				"\nOptional:\n" +
				"-q Primary name to replace the ID column.\n" +
				"-t Secondary name to replace the ID column.\n" +
				"-u Use the primary #CHROM line and skip this header check. Not recommended!\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFConsensus -p illumina.vcf -q Strelka\n"+
				"-s stnd.indel.vcf.gz -t Scalpel -o indelCalls.vcf.gz \n\n"+

		"**************************************************************************************\n");

	}

}
