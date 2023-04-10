package edu.utah.seq.vcf.anno;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;


/**Filter for annotated vcf files.  */
public class AnnotatedVcfFilter {

	private File[] vcfFiles;
	private int minimumTumorReadDepth = 0;
	private int minimumNormalReadDepth = 0;
	private double minimumTNFractionDiff = 0;
	private double minQual = 0;
	private int minimumAltReadDepth = 0;
	private double minimumTNRatio = 0;
	private double maximumNormalAltFraction = 1;
	private double minimumTumorAltFraction = 0;
	private String[] filters = null;
	private File saveDirectory = null;
	
	public AnnotatedVcfFilter (String[] args) {

		processArgs(args);
		
		System.out.println("Thresholds for Tumor and Normal:");
		System.out.println(minimumNormalReadDepth+"\tMin N alignment depth");
		System.out.println(minimumTumorReadDepth+"\tMin T alignment depth");
		System.out.println(minimumAltReadDepth+"\tMin T alt count");
		System.out.println(minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		System.out.println(minimumTNRatio+"\tMin T/N allelic fraction ratio");
		System.out.println(maximumNormalAltFraction+"\tMax N allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin T allelic fraction");
		System.out.println(minQual+"\tMin QUAL score");
		System.out.println(saveDirectory+"\tSave directory.");
		if (filters != null) System.out.println(Misc.stringArrayToString(filters, ", ") +"\t Exclusion FILTERs");
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);
			
			for (VCFRecord r: parser.getVcfRecords()){			

				//check filters?
				if (filters != null){
					for (String f: filters) {
						if (r.getFilter().contains(f)) {
							r.setFilter(VCFRecord.FAIL);
							break;
						}
						
					}
				}
				if (r.getFilter().equals(VCFRecord.FAIL)) continue;
				r.setFilter(VCFRecord.PASS);

				//check QUAL			
				if (r.getQuality()< minQual) {
					r.setFilter(VCFRecord.FAIL);
					continue;
				}

				//pull info fields
				String tumorDpString = r.getInfoObject().getInfo("T_DP");
				String normalDpString = r.getInfoObject().getInfo("N_DP");
				String tumorAfString = r.getInfoObject().getInfo("T_AF");
				String normalAfString = r.getInfoObject().getInfo("N_AF");
				if (tumorDpString==null || normalDpString==null || tumorAfString==null || normalAfString==null) {
					Misc.printErrAndExit("\nFailed to find the required INFO fields T_DP, N_DP, T_AF, N_AF in "+r.getOriginalRecord());
				}
				
				double tDp = Double.parseDouble(tumorDpString);
				double nDp = Double.parseDouble(normalDpString);
				double tAf = Double.parseDouble(tumorAfString);
				double nAf = Double.parseDouble(normalAfString);


				if (nDp < minimumNormalReadDepth || tDp < minimumTumorReadDepth){
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				
				//check alt counts
				int tumorAltCounts = (int)Math.round(tAf * tDp);
				if (tumorAltCounts < minimumAltReadDepth) {
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				
				//check allelic ratio diff
				if (minimumTNFractionDiff != 0){
					double change = tAf-nAf;
					if (change < minimumTNFractionDiff) {
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check T/N AF ratio
				if (minimumTNRatio != 0 && nAf !=0){
					double change = tAf/nAf;
					if (change < minimumTNRatio){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check normal alt fraction?
				if (maximumNormalAltFraction !=1){
					if (nAf > maximumNormalAltFraction) {
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check tumor alt fraction?
				if (minimumTumorAltFraction !=0){
					if (tAf < minimumTumorAltFraction) {
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
			}
			printRecords(parser);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void printRecords(VCFParser parser) throws Exception{
		String name = Misc.removeExtension(parser.getVcfFile().getName());
		
		Gzipper outPass = new Gzipper (new File(saveDirectory, name+".passAVF.vcf.gz"));
		Gzipper outFail = new Gzipper (new File(saveDirectory, name+".failAVF.vcf.gz"));
		
		int numFail = 0;
		int numPass = 0;

		//print header
		for (String h : parser.getVcfComments().getOriginalComments()) {
			outPass.println(h);
			outFail.println(h);
		}
		
		//print records
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(VCFRecord.FAIL)) {
				numFail++;
				outFail.println(vcf.getOriginalRecord());
			}
			else {
				numPass++;
				outPass.println(vcf.getOriginalRecord());
			}
		}
		outPass.close();
		outFail.close();
		System.out.println(numPass+"\t"+numFail);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AnnotatedVcfFilter(args);
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
					case 'v': forExtraction = new File(args[++i]).getCanonicalFile(); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'a': minimumAltReadDepth = Integer.parseInt(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
					case 'q': minQual = Double.parseDouble(args[++i]); break;
					case 'f': filters = Misc.COMMA.split(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
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
		if (saveDirectory == null) Misc.printExit("\nError: cannot find your save directory!\n");
		saveDirectory.mkdirs();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Annotated Vcf Filter: Feb 2023                          **\n" +
				"**************************************************************************************\n" +
				"Filters vcf records that have T_DP, N_DP, T_AF, and N_AF INFO field attributes.\n"+

				"\nOptions:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-s Directory to save the parsed files.\n"+
				"-t Minimum tumor allele frequency (AF), defaults to 0.\n"+
				"-n Maximum normal AF, defaults to 1.\n"+
				"-u Minimum tumor alignment depth, defaults to 0.\n"+
				"-a Minimum tumor alt count, defaults to 0.\n"+
				"-o Minimum normal alignment depth, defaults to 0.\n"+
				"-d Minimum T-N AF difference, defaults to 0.\n"+
				"-r Minimum T/N AF ratio, defaults to 0.\n"+
				"-q Minimum QUAL score, defaults to 0.\n"+
				"-f Comma delimited list of FILTER keys to fail records.\n"+
				

				"\nExample: java -jar pathToUSeq/Apps/TNScopeVCFParser -v ToParse/ -n 0.5 -u 100\n"+
				"        -o 20 -d 0.05 -r 2 -a 3 -t 50 -p -s FilteredVcfFiles/ -f\n"+
				"        multi_event_alt_allele_in_normal,str_contraction \n\n"+


				"**************************************************************************************\n");

	}

}
