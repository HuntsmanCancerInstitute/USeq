package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Simple Sentieon TNScope vcf parser.  */
public class TNScopeVCFParser {

	private File[] vcfFiles;
	private int minimumTumorReadDepth = 0;
	private int minimumNormalReadDepth = 0;
	private double minimumTNFractionDiff = 0;
	private double minQual = 0;
	private int minimumAltReadDepth = 0;
	private double minimumTNRatio = 0;
	private double maximumNormalAltFraction = 1;
	private double minimumTumorAltFraction = 0;
	private boolean excludeNonPass = false;
	private String afInfo = "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private Pattern tLod = Pattern.compile(".+;TLOD=([\\d\\.]+).+"); 
	private Pattern dp = Pattern.compile("DP=\\d+;");
	private File saveDirectory = null;
	
	public TNScopeVCFParser (String[] args) {

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
		System.out.println(excludeNonPass+"\tRemove non PASS FILTER field records.");
		System.out.println(saveDirectory+"\tSave directory.");
		
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
				VCFSample[] tumorNormal = r.getSample();
				if (tumorNormal.length !=2) Misc.printErrAndExit("\nLooks like there aren't two samples present in this VCF, aborting.\n");
				
				//check PASS?
				if (excludeNonPass){
					if (r.getFilter().toLowerCase().contains("pass") == false){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				
				r.setFilter(VCFRecord.PASS);

				//check QUAL			
				if (r.getQuality()< minQual) {
					r.setFilter(VCFRecord.FAIL);
					continue;
				}

				
				//check depth
				int normDepth = tumorNormal[1].getReadDepthDP();
				int tumDepth = tumorNormal[0].getReadDepthDP();
				


				if (normDepth < minimumNormalReadDepth || tumDepth < minimumTumorReadDepth){
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				
				//check alt counts
				int tumorAltCounts = Integer.parseInt(tumorNormal[0].getAlternateCounts());
				if (tumorAltCounts < minimumAltReadDepth) {
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				
				//get alt counts from normal
				int normAltCounts = Integer.parseInt(tumorNormal[1].getAlternateCounts());
				
				//calc AFs, don't use the AF in the sample, then set the recalc
				double normRto = (double)normAltCounts / (double)normDepth;
				double tumRto = (double)tumorAltCounts / (double)tumDepth;
				tumorNormal[1].setAlleleFreqAF(normRto);
				tumorNormal[0].setAlleleFreqAF(tumRto);
				
				//check allelic ratio diff
				if (minimumTNFractionDiff != 0){
					double change = tumRto-normRto;
					if (change < minimumTNFractionDiff) {
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check T/N AF ratio
				if (minimumTNRatio != 0 && normRto !=0){
					double change = tumRto/normRto;
					if (change < minimumTNRatio){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check normal alt fraction?
				if (maximumNormalAltFraction !=1){
					if (normRto > maximumNormalAltFraction) {
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				//check tumor alt fraction?
				if (minimumTumorAltFraction !=0){
					if (tumRto < minimumTumorAltFraction) {
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
			}
			File outFile = new File (saveDirectory, Misc.removeExtension(vcf.getName())+"_Filtered.vcf.gz");
			printRecords(parser, outFile);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void printRecords(VCFParser parser, File f) throws Exception{
		
		Gzipper out = new Gzipper (f);
		
		int numFail = 0;
		int numPass = 0;

		//print header
		boolean addInfo = true;
		for (String h : parser.getStringComments()) {
			//replace chrom line
			if (h.startsWith("#CHROM")) h= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR";
			else if (addInfo && h.startsWith("##INFO")){
				out.println(afInfo);
				out.println(dpInfo);
				addInfo = false;
			}
			out.println(h);
		}
		
		//print records
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(VCFRecord.FAIL)) numFail++;
			else {
				String[] t = Misc.TAB.split(vcf.getOriginalRecord());
				
				//add DP and AF for tumor to INFO
				String tumorAF = Num.formatNumberJustMax(vcf.getSample()[0].getAltRatio(), 4);
				t[7] = "DP=" + vcf.getSample()[0].getReadDepthDP()+ ";AF=" + tumorAF + ";"+ t[7] ;
				
				numPass++;
				t[2] = "TNScope_"+numPass;
				
				//swap tumor and normal samples
				//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tumor	normal
				//   0       1  2    3   4   5        6      7        8       9       10
				String tumor = new String(t[9]);
				String norm = new String (t[10]);
				t[9] = norm;
				t[10] = tumor;
				
				out.println(Misc.stringArrayToString(t, "\t"));
			}
		}
		out.close();
		System.out.println(numPass+"\t"+numFail);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TNScopeVCFParser(args);
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
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'a': minimumAltReadDepth = Integer.parseInt(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
					case 'l': minQual = Double.parseDouble(args[++i]); break;
					case 'p': excludeNonPass = true; break;
					case 'f': saveDirectory = new File(args[++i]); break;
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
		if (saveDirectory == null) saveDirectory = vcfFiles[0].getParentFile();
		else saveDirectory.mkdirs();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              TNScope VCF Parser: Aug 2019                        **\n" +
				"**************************************************************************************\n" +
				"Parses Sentieon TNScope VCF files. Filters for read depth, allele frequency, QUAL, alt\n"+
				"count, etc. Inserts tumor AF and DP into the INFO field. Swaps the order of the tumor\n"+
				"and normal samples to enable downstream merging.\n"+

				"\nOptions:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-f Directory to save the parsed files, defaults to the parent dir of the first vcf.\n"+
				"-t Minimum tumor allele frequency (AF), defaults to 0.\n"+
				"-n Maximum normal AF, defaults to 1.\n"+
				"-u Minimum tumor alignment depth, defaults to 0.\n"+
				"-a Minimum tumor alt count, defaults to 0.\n"+
				"-o Minimum normal alignment depth, defaults to 0.\n"+
				"-d Minimum T-N AF difference, defaults to 0.\n"+
				"-r Minimum T/N AF ratio, defaults to 0.\n"+
				"-t Minimum QUAL score, defaults to 0.\n"+
				"-p Remove non PASS filter field records.\n"+
				

				"\nExample: java -jar pathToUSeq/Apps/TNScopeVCFParser -v /VCFFiles/ -n 0.5 -u 100\n"+
				"        -o 20 -d 0.05 -r 2 -a 3 -t 50 \n\n"+


				"**************************************************************************************\n");

	}

}
