package edu.utah.seq.vcf;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.vcf.fdr.VCFFdrEstimator;
import util.gen.*;

/**Memory intensive yet fast process for collapsing records, tosses QUAL onward
 * @author Nix
 * https://bioinformatics.stackexchange.com/questions/6826/sort-vcf-by-contig-and-position-within-contig
 * */
public class VCFCollapser {

	//user fields
	private File[] vcfFiles = null;
	private File mergedVcf = null;
	private int numberRecords = 0;
	private int numberCollapsedRecords = 0;
	private HashSet<String> records = new HashSet<String>();
	private StringBuilder header = null;
	
	//constructor
	public VCFCollapser(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//parse each vcf
		IO.pl("Parsing...");
		for (File vcfFile: vcfFiles) parseVcfFile(vcfFile);
		numberCollapsedRecords = records.size();
		
		//write out merged file
		printMergedVcf();
		
		IO.pl("\n"+numberRecords+" -> "+ numberCollapsedRecords);
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void printMergedVcf() {
		Gzipper out = null;
		try {
			out = new Gzipper(mergedVcf);
			out.println("##fileformat=VCFv4.2");
			out.print(header.toString());
			out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
			for (String key: records) {
				String[] t = Misc.UNDERSCORE.split(key);
				StringBuilder sb = new StringBuilder(t[0]);
				sb.append("\t");
				sb.append(t[1]);
				sb.append("\t.\t");
				sb.append(t[2]);
				sb.append("\t");
				sb.append(t[3]);
				sb.append("\t0\t.\t.");
				out.println(sb.toString());
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		out.closeNoException();
	}

	private void parseVcfFile (File vcf) {
		IO.pl("\t"+vcf.getName());
		
		boolean addContigs = false;
		if (header == null) {
			header = new StringBuilder();
			addContigs = true;
		}
		
		BufferedReader in = null;
		String line = null;
		try {
			in = IO.fetchBufferedReader(vcf);
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
					//   0       1   2   3   4   
					numberRecords++;
					String[] fields = Misc.TAB.split(line);
					StringBuilder sb = new StringBuilder(fields[0]);
					sb.append("_");
					sb.append(fields[1]);
					sb.append("_");
					sb.append(fields[3]);
					sb.append("_");
					sb.append(fields[4]);
					records.add(sb.toString());
				}
				else if (addContigs && line.startsWith("##contig=")) {
					header.append(line);
					header.append("\n");
				}
			}
		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing "+vcf+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFCollapser(args);
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
					case 'd': forExtraction = new File(args[++i]); break;
					case 'v': mergedVcf = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		vcfFiles = VCFFdrEstimator.fetchVcfFiles(forExtraction);

	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Collapser : June 2022                           **\n" +
				"**************************************************************************************\n" +
				"Merges records with the same CHROM POS ID REF ALT and outputs an unsorted stripped\n"+
				"down vcf file. \n"+

				"\nRequired Params:\n"+
				"-d Directory containing  vcf files for merging (xxx.vcf(.gz/.zip OK))\n"+
				"-v Path to a vcf file ending in .gz to write the unsorted results. \n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/VCFCollapser -d VcfsToMerge/\n" +
				"       -v merged.vcf.gz && bcftools sort merged.vcf.gz > sorted.vcf \n"+

		"\n**************************************************************************************\n");

	}
}
