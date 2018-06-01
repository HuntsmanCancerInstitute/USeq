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

/**Simple Mutect vcf parser.  */
public class MutectVCFParser {

	private File[] vcfFiles;
	private int minimumTumorReadDepth = 0;
	private int minimumNormalReadDepth = 0;
	private double minimumTNFractionDiff = 0;
	private int minimumAltReadDepth = 0;
	private double minimumTNRatio = 0;
	private double maximumNormalAltFraction = 1;
	private double minimumTumorAltFraction = 0;
	private boolean excludeNonPass = false;
	private boolean printSpreadsheet = false;
	private String afInfo = "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private File saveDirectory = null;
	
	public MutectVCFParser (String[] args) {

		processArgs(args);
		
		System.out.println("Thresholds for Tumor and Normal:");
		System.out.println(minimumNormalReadDepth+"\tMin N alignment depth");
		System.out.println(minimumTumorReadDepth+"\tMin T alignment depth");
		System.out.println(minimumAltReadDepth+"\tMin T alt count");
		System.out.println(minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		System.out.println(minimumTNRatio+"\tMin T/N allelic fraction ratio");
		System.out.println(maximumNormalAltFraction+"\tMax N allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin T allelic fraction");
		System.out.println(excludeNonPass+"\tRemove non PASS FILTER field records.");
		System.out.println(printSpreadsheet+"\tPrint spreadsheet output.");
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
			File txt = null;
			Gzipper out = null;
			if (printSpreadsheet){
				txt = new File (saveDirectory, Misc.removeExtension(vcf.getName())+".txt.gz");
				out = new Gzipper(txt);
				out.println("#PASS\tCHROM\tPOS\tREF\tALT\tT_AF\tT_DP\tN_AF\tN_DP\tFILTER\tINFO");
			}
			
			for (VCFRecord r: parser.getVcfRecords()){			
				VCFSample[] tumNorm = r.getSample();
				if (tumNorm.length !=2) Misc.printErrAndExit("\nLooks like there aren't two samples present in this VCF, aborting.\n");
				
				boolean pass = true;
				double normRto = tumNorm[1].getAltRatio();
				double tumRto = tumNorm[0].getAltRatio();
				
				//check depth
				int normDepth = tumNorm[1].getReadDepthDP();
				int tumDepth = tumNorm[0].getReadDepthDP();
				if (normDepth < minimumNormalReadDepth || tumDepth < minimumTumorReadDepth) pass = false;
				
				//check alt counts
				int altCounts = Integer.parseInt(tumNorm[0].getAlternateCounts());
				if (altCounts < minimumAltReadDepth) pass = false;
				
				//check allelic ratio diff
				if (pass && minimumTNFractionDiff != 0){
					double change = tumRto-normRto;
					if (change < minimumTNFractionDiff) pass = false;
				}
				
				//check T/N AF ratio
				if (pass && minimumTNRatio != 0 && normRto !=0){
					double change = tumRto/normRto;
					if (change < minimumTNRatio) pass = false;
				}
				
				//check normal alt fraction?
				if (pass && maximumNormalAltFraction !=1){
					if (normRto > maximumNormalAltFraction) pass = false;
				}
				//check tumor alt fraction?
				if (pass && minimumTumorAltFraction !=0){
					if (tumRto < minimumTumorAltFraction) pass = false;
				}
				//check PASS?
				if (pass && excludeNonPass && r.getFilter().toLowerCase().contains("pass") == false){
					pass = false;
				}
				
				//build txt output
				ArrayList<String> al = new ArrayList<String>();
				al.add(pass+"");
				al.add(r.getChromosome());
				al.add((r.getPosition()+1)+"");
				al.add(r.getReference());
				al.add(Misc.stringArrayToString(r.getAlternate(), ","));
				al.add(tumRto+"");
				al.add(tumDepth+"");
				al.add(normRto+"");
				al.add(normDepth+"");
				al.add(r.getFilter());
				al.add(r.getInfoObject().getInfoString());
				String line = Misc.stringArrayListToString(al, "\t");
				if (printSpreadsheet) out.println(line);
				if (pass) r.setFilter(VCFRecord.PASS);
				else r.setFilter(VCFRecord.FAIL);

			}
			if (printSpreadsheet) out.close();
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
			//watch out for #CHROM line, switch TUMOR NORMAL to NORMAL and TUMOR to match others
			if (h.startsWith("#CHROM")){
				String[] t = Misc.TAB.split(h);
				if (t[10].equals("NORMAL") == false || t[9].equals("TUMOR") == false) Misc.printErrAndExit("Problem, #CHROM doesn't end with TUMOR and NORMAL? "+t);
				String tu = t[9];
				String no = t[10];
				t[9] = no;
				t[10] = tu;
				h = Misc.stringArrayToString(t, "\t");
			}
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
				
				//flip t[9] and t[10]
				String tu = t[9];
				String no = t[10];
				t[9] = no;
				t[10] = tu;
				numPass++;
				t[2] = "Mutect_"+numPass;
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
		new MutectVCFParser(args);
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
					case 'p': excludeNonPass = true; break;
					case 's': printSpreadsheet = true; break;
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
				"**                               Mutect VCF Parser: May 2018                        **\n" +
				"**************************************************************************************\n" +
				"Parses Mutect2 VCF files, filtering for read depth, allele frequency diff ratio, etc.\n"+
				"Inserts AF and DP into for the tumor sample into the INFO field. Changes the sample\n"+
				"order to Normal and Tumor and updates the #CHROM line.\n"+

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
				"-p Remove non PASS filter field records.\n"+
				"-s Print spreadsheet variant summary.\n"+
				

				"\nExample: java -jar pathToUSeq/Apps/MutectVCFParser -v /VCFFiles/ -t 0.05 -n 0.5 -u 100\n"+
				"        -o 20 -d 0.05 -r 2 -a 3 \n\n"+


				"**************************************************************************************\n");

	}

}
