package edu.utah.seq.vcf;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.vcf.fdr.VCFFdrEstimator;
import util.gen.*;

/**Memory intensive yet fast process for appending FILTER and INFO fields to exact matching records
 * @author Nix
 * https://bioinformatics.stackexchange.com/questions/6826/sort-vcf-by-contig-and-position-within-contig
 * */
public class VCFInfoAppender {

	//user fields
	private File[] vcfFiles = null;
	private File vcfWithInfo = null;
	private File saveDirectory = null;

	private HashMap<String,String[]> keyInfo = new HashMap<String,String[]>();
	private LinkedHashSet<String> infoVcfHeader = new LinkedHashSet<String>();
	
	//constructor
	public VCFInfoAppender(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//load the vcf with filter info
		loadAnnotatedVcf();
		
		//parse each vcf
		IO.pl("Appending FILTER and INFO...");
		for (File vcfFile: vcfFiles) parseVcfFile(vcfFile);
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void loadAnnotatedVcf () {
		IO.pl("Loading the annotated vcf...");
		BufferedReader in = null;
		String line = null;
		
		try {
			in = IO.fetchBufferedReader(vcfWithInfo);
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO	FORMAT
					//   0       1   2   3   4   5      6     7   
					String[] fields = Misc.TAB.split(line);
					StringBuilder sb = new StringBuilder(fields[0]);
					sb.append("_");
					sb.append(fields[1]);
					sb.append("_");
					sb.append(fields[3]);
					sb.append("_");
					sb.append(fields[4]);
					String key = sb.toString();
					String[] fi = new String[] {fields[6], fields[7]};
					keyInfo.put(key, fi);
				}
				else if (line.startsWith("#CHROM") == false && line.startsWith("##fileformat=") == false) infoVcfHeader.add(line);
			}
		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing "+vcfWithInfo+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
	}

	private void parseVcfFile (File vcf) {
		IO.pl("\t"+vcf.getName());
		Gzipper out = null;
		BufferedReader in = null;
		String line = null;
		HashSet<String> header = new HashSet<String>();
		try {
			out = new Gzipper(new File (saveDirectory, vcf.getName()));
			in = IO.fetchBufferedReader(vcf);
			while ((line = in.readLine())!= null) {
				//a data line
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT
					//   0       1   2   3   4   5      6     7    8
					//make the key
					String[] fields = Misc.TAB.split(line);
					StringBuilder sb = new StringBuilder(fields[0]);
					sb.append("_");
					sb.append(fields[1]);
					sb.append("_");
					sb.append(fields[3]);
					sb.append("_");
					sb.append(fields[4]);
					String key = sb.toString();
					//pull the key from the info vcf
					String[] filterInfo = keyInfo.get(key);
					if (filterInfo == null) Misc.printErrAndExit("\t\tFailed to find "+key+" for "+line);
					else {
						if (filterInfo[1].endsWith(";.")) filterInfo[1] = filterInfo[1].substring(0, filterInfo[1].length()-2);
						
						sb = new StringBuilder(fields[0]);
						for (int i=1; i< 6; i++) {
							sb.append("\t");
							sb.append(fields[i]);
						}
						sb.append("\t");
						
						//FILTER
						if (filterInfo[0].equals(".") || filterInfo[0].toUpperCase().equals("PASS")) {
							//don't append filter
						}
						else {
							//append or replace
							if (fields[6].equals(".") || fields[6].toUpperCase().equals("PASS") || fields[6].length()==0) fields[6] = filterInfo[0];
							else fields[6] = filterInfo[0] + ";"+ fields[6];
						}
						sb.append(fields[6]);
						sb.append("\t");
						
						//INFO
						//append or replace				
						if (fields[7].equals(".") || fields[7].length()==0) fields[7] = filterInfo[1];
						else fields[7] = filterInfo[1] + ";"+ fields[7];
						sb.append(fields[7]);
						
						//Last
						for (int i=8; i< fields.length; i++) {
							sb.append("\t");
							sb.append(fields[i]);
						}
						out.println(sb.toString());
					}
				}
				//last header line!
				else if (line.startsWith("#CHROM")) {
					//output any new infoHeader lines
					for (String ivh: infoVcfHeader) {
						if (header.contains(ivh) == false) {
							out.println(ivh);
						}
					}
					
					// add the final CHROM line
					out.println(line);
				}
				//a comment line, print it and save to the hash
				else {
					header.add(line);
					out.println(line);
				}
			}
		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing "+vcf+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
			if (out != null) out.closeNoException();
		}
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFInfoAppender(args);
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
					case 'v': vcfWithInfo = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
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
		saveDirectory.mkdirs();
	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Info Appender : June 2022                        **\n" +
				"**************************************************************************************\n" +
				"Appends the FILTER and INFO from one vcf to the others with matching chr_pos_ref_alt.\n"+

				"\nRequired Params:\n"+
				"-d Directory containing  vcf files for appending info (xxx.vcf(.gz/.zip OK))\n"+
				"-v Path to a master vcf file containing info to append. Use the VCFCollapser to merge\n"+
				"      all the variants and then annotate it, e.g. VCFCallFrequency, VCFSpliceScanner\n"+
				"-s Save directory.\n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/VCFInfoAppender -d VcfsToModify/\n" +
				"       -v master.callFreq.vcf.gz -s VcfsAppended/ \n"+

		"\n**************************************************************************************\n");

	}
}
