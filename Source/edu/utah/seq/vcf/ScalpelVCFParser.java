package edu.utah.seq.vcf;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Scalpel filtering app. */
public class ScalpelVCFParser {

	private File[] vcfFiles;
	private float minimumQual = 0;
	private int minimumTumorReadDepth = 0;
	private int minimumNormalReadDepth = 0;
	private double minimumTNRatio = 0;
	private double minimumTNFractionDiff = 0;
	private double maximumNormalAltFraction = 1;
	private double minimumTumorAltFraction = 0;
	private String afFormat = "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">";
	private String afInfo = "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	
	public ScalpelVCFParser (String[] args) { 

		processArgs(args);
		
		System.out.println("Thresholds for Tumor and Normal:");
		System.out.println(minimumQual+"\tMin QUAL score");
		System.out.println(minimumTumorReadDepth+"\tMin T alignment depth");
		System.out.println(minimumNormalReadDepth+"\tMin N alignment depth");
		System.out.println(minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		System.out.println(minimumTNRatio+"\tMin T/N allelic fraction ratio");
		System.out.println(maximumNormalAltFraction+"\tMax N allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin T allelic fraction");
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}
	

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);
			
			File txt = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+".txt.gz");
			Gzipper out = new Gzipper(txt);
			out.println("#PASS\tCHROM\tPOS\tREF\tALT\tT_AF\tT_DP\tN_AF\tN_DP\tFILTER\tINFO");
			
			for (VCFRecord r: parser.getVcfRecords()){	
				VCFSample[] normTum = r.getSample();
				int normDepth = normTum[0].getReadDepthDP();
				int tumDepth = normTum[1].getReadDepthDP();	
				double normRto = normTum[0].getAltRatio();
				double tumRto = normTum[1].getAltRatio();
				boolean pass = true;
				
				//check depth
				if (normDepth < minimumNormalReadDepth || tumDepth < minimumTumorReadDepth) pass = false;
				
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
				out.println(line);
				if (pass) r.setFilter(VCFRecord.PASS);
				else r.setFilter(VCFRecord.FAIL);
			}
			
			out.close();
			File outFile = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"_Filtered.vcf.gz");
			printRecords(parser, outFile);
			

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void printRecords(VCFParser parser, File f) throws Exception{
		Gzipper out = new Gzipper (f);
		
		//write out header inserting AF
		writeHeaderWithExtraInfo(out, parser);
		
		int numFail = 0;
		int numPass = 0;
		int counter = 0;
		
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			boolean fail = vcf.getFilter().equals(VCFRecord.FAIL);
			if (fail) numFail++;
			else numPass++;
			if (fail == false){
				String orig = vcf.toString();
				//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOR......
				//   0    1   2  3   4   5     6     7      8      9     10
				String[] fields = VCFParser.TAB.split(orig);
				//reset ID
				fields[2] = "Scalpel_"+counter;
				counter++;
				//add DP and AF for tumor to INFO
				String tumorAF = formatAf (vcf.getSample()[1].getAltRatio());
				fields[7] = "DP=" + vcf.getSample()[1].getReadDepthDP()+ ";AF=" + tumorAF + ";"+ fields[7] ;
				//add af to format
				fields[8] = fields[8]+ ":AF";
				//add af to Norm and Tum
				fields[9] = fields[9]+ ":"+ formatAf (vcf.getSample()[0].getAltRatio());
				fields[10] = fields[10]+ ":"+ tumorAF;
				out.println(Misc.stringArrayToString(fields, "\t"));
			}
		}
		out.close();
		System.out.println(numPass+"\t"+numFail);
	}
	
	private String formatAf(double af){
		if (af == 0.0) return "0";
		return Num.formatNumberJustMax(af, 4);
	}
	
	private void writeHeaderWithExtraInfo(Gzipper out, VCFParser parser) throws IOException {
		for (String h : parser.getStringComments()) {
			if (h.startsWith("#CHROM")) {
				out.println(afFormat);
				out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR");
			}
			else if (h.startsWith("##INFO=<ID=DENOVO")){
				out.println(afInfo);
				out.println(dpInfo);
			}
			else if (h.startsWith("##source=scalpel")){
				String[] t = Misc.PATTERN_EQUALS.split(h);
				out.println("##source=scalpel");
				out.println("##version="+t[1]);
			}
			else out.println(h);
		}
	}


	

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ScalpelVCFParser(args);
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
					case 'm': minimumQual = Float.parseFloat(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
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
				"**                            Scalpel VCF Parser: Jan 2017                         **\n" +
				"**************************************************************************************\n" +
				"Filters Scalpel VCF INDEL files for various thresholds.  Adds tumor DP and AF values\n"+
				"to the info field.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				"-m Minimum QUAL score, defaults to 0\n"+
				"-t Minimum tumor allele frequency (AF), defaults to 0.\n"+
				"-n Maximum normal AF, defaults to 1.\n"+
				"-u Minimum tumor alignment depth, defaults to 0.\n"+
				"-o Minimum normal alignment depth, defaults to 0.\n"+
				"-d Minimum T-N AF difference, defaults to 0.\n"+
				"-r Minimum T/N AF ratio, defaults to 0.\n"+

				"\nExample: java -jar pathToUSeq/Apps/ScalpelVCFParser -v /VCFFiles/ -t 0.05 -n 0.5 -u 100\n"+
				"        -o 20 -d 0.05 -r 2\n\n"+


				"**************************************************************************************\n");

	}
}
