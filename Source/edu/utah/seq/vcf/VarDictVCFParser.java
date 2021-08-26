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

/**Simple VarDict  vcf formatter and parser. 
 */
public class VarDictVCFParser {

	private File[] vcfFiles;
	private float minimumScore = 0;
	private int minimumTumorReadDepth = 0;
	private int minimumNormalReadDepth = 0;
	private double minimumTNFractionDiff = 0;
	private double minimumTNRatio = 0;
	private double maximumNormalAltFraction = 1;
	private int minimumAltReadDepth = 0;
	private Pattern dp = Pattern.compile("DP=\\d+;*");
	private Pattern af = Pattern.compile("AF=\\d+;*");
	private boolean printSpreadsheet = false;
	private String afInfo = "##INFO=<ID=T_AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=T_DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private String nafInfo = "##INFO=<ID=N_AF,Number=1,Type=Float,Description=\"Allele Frequency for normal\">";
	private String ndpInfo = "##INFO=<ID=N_DP,Number=1,Type=Integer,Description=\"Read depth for normal\">";
	private String afFormat = "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">";
	private String gtFormat = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	private double minimumTumorAltFraction = 0;
	private boolean excludeNonSomatic = true;
	private File saveDirectory = null;

	
	public VarDictVCFParser (String[] args) { 

		processArgs(args);
		
		System.out.println("Thresholds for Tumor and Normal:");
		System.out.println(minimumNormalReadDepth+"\tMin N alignment depth");
		System.out.println(minimumTumorReadDepth+"\tMin T alignment depth");
		System.out.println(minimumAltReadDepth+"\tMin T alt count");
		System.out.println(minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		System.out.println(minimumTNRatio+"\tMin T/N allelic fraction ratio");
		System.out.println(maximumNormalAltFraction+"\tMax N allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin T allelic fraction");
		System.out.println(minimumScore+"\tMin QSI or QSS score");
		System.out.println(printSpreadsheet+"\tPrint spreadsheet output");

		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
		
		System.out.println("\nComplete!");
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
				
				VCFSample tum = r.getSample()[0];
				VCFSample norm = r.getSample()[1];
				
				//force it to calc AF
				tum.setAlleleFreqAF(-1);
				norm.setAlleleFreqAF(-1);
				
				//set AF based on counts
				double normRto = norm.getAltRatio();
				double tumRto = tum.getAltRatio();
				boolean pass = true;
				
				//check depth
				int normDepth = norm.getReadDepthDP();
				int tumDepth = tum.getReadDepthDP();
				
				ArrayList<String> failing = new ArrayList<String>();
				
				if (normDepth < minimumNormalReadDepth) {
					pass = false;
					failing.add("N_DP");
				}
				if (tumDepth < minimumTumorReadDepth) {
					pass = false;
					failing.add("T_DP");
				}
				
				//check alt counts
				int altCounts = Integer.parseInt(tum.getAlternateCounts());
				if (altCounts < minimumAltReadDepth) {
					pass = false;
					failing.add("T_Alt");
				}
				
				//check allelic ratio diff
				if (minimumTNFractionDiff != 0){
					double change = tumRto-normRto;
					if (change < minimumTNFractionDiff) {
						pass = false;
						failing.add("Diff");
					}
				}

				//check T/N AF ratio
				if (minimumTNRatio != 0 && normRto !=0){
					double change = tumRto/normRto;
					if (change < minimumTNRatio) {
						pass = false;
						failing.add("Ratio");
					}
				}
				
				//check normal alt fraction?
				if (maximumNormalAltFraction !=1){
					if (normRto > maximumNormalAltFraction) {
						pass = false;
						failing.add("N_AF");
					}
				}
				//check tumor alt fraction?
				if (minimumTumorAltFraction !=0){
					if (tumRto < minimumTumorAltFraction) {
						pass = false;
						failing.add("T_AF");
					}
				}
				//check PASS and somatic?
				if (excludeNonSomatic) {
					if (r.getUnmodifiedInfoString().contains("STATUS=StrongSomatic") == false && r.getUnmodifiedInfoString().contains("STATUS=LikelySomatic") == false || r.getFilter().equals(VCFRecord.PASS) == false) {
						pass = false;
						failing.add("STATUS");
					}
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
				else r.setFilter(Misc.stringArrayListToString(failing, ","));
			}
			if (printSpreadsheet) out.close();
			
			File pass = new File (saveDirectory, Misc.removeExtension(vcf.getName())+"_Pass.vcf.gz");
			File fail = new File (saveDirectory, Misc.removeExtension(vcf.getName())+"_Fail.vcf.gz");
			printRecords(parser, pass, fail);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void printRecords(VCFParser parser, File pass, File fail) throws Exception{
		Gzipper passOut = new Gzipper (pass);
		Gzipper failOut = new Gzipper (fail);
		
		//write out header inserting AF for pass
		writeHeaderWithExtraInfo(passOut, parser);
		for (String s: parser.getStringComments()) failOut.println(s);
		
		
		int numFail = 0;
		int numPass = 0;
		
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			String orig = vcf.toString();
			//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT  TUMOR NORMAL......
			//   0    1   2  3   4   5     6     7      8      9     10
			String[] fields = VCFParser.TAB.split(orig);
			
			if (vcf.getFilter().equals(VCFRecord.PASS)) {
				
				numPass++;
				VCFSample tum = vcf.getSample()[0];
				VCFSample norm = vcf.getSample()[1];
				
				//reset score
				fields[5] = Float.toString(vcf.getQuality());
				//reset ID
				fields[2] = "VarDict_"+numPass;
				//remove existing DP
				fields[7] = dp.matcher(fields[7]).replaceFirst("");
				//remove existing AF
				fields[7] = af.matcher(fields[7]).replaceFirst("");
				
				//modify INFO
				String normalAf = formatAf (norm.getAltRatio());
				String tumorAf = formatAf (tum.getAltRatio());
				
				//add DP and AF for tumor to INFO
				fields[7] = "T_DP=" + tum.getReadDepthDP()+ ";T_AF=" + tumorAf+ ";N_DP=" + norm.getReadDepthDP()+ ";N_AF=" + normalAf + ";"+ fields[7] ;
				
				passOut.println(Misc.stringArrayToString(fields, "\t"));
			}
			else {
				numFail++;
				fields[6] = vcf.getFilter();
				failOut.println((Misc.stringArrayToString(fields, "\t")));
			}
		}
		
		passOut.close();
		failOut.close();
		System.out.println(numPass+"\t"+numFail);
	}
	
	private String formatAf(double af){
		if (af == 0.0) return "0";
		return Num.formatNumberJustMax(af, 5);
	}

	private void writeHeaderWithExtraInfo(Gzipper out, VCFParser parser) throws IOException {
		boolean added = false;
		for (String h : parser.getStringComments()) {
			//skip existing DP and AF in info
			if (h.startsWith("##INFO=<ID=DP,") || h.startsWith("##INFO=<ID=AF,")) continue;
			if (added == false && h.startsWith("##FORMAT")) {
				out.println(afInfo);
				out.println(dpInfo);
				out.println(nafInfo);
				out.println(ndpInfo);
				out.println(gtFormat);
				out.println(afFormat);
				added = true;
			}
			out.println(h);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VarDictVCFParser(args);
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
					case 'm': minimumScore = Float.parseFloat(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'a': minimumAltReadDepth = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
					case 'p': excludeNonSomatic = false; break;
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
				"**                            VarDict VCF Parser: July 2021                         **\n" +
				"**************************************************************************************\n" +
				"Parses and filters VarDict VCF files and inserts the tumor and normal DP and AF info.\n"+
				"\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-t Minimum tumor allele frequency (AF), defaults to 0.\n"+
				"-n Maximum normal AF, defaults to 1.\n"+
				"-u Minimum tumor alignment depth, defaults to 0.\n"+
				"-a Minimum tumor alt count, defaults to 0.\n"+
				"-o Minimum normal alignment depth, defaults to 0.\n"+
				"-d Minimum T-N AF difference, defaults to 0.\n"+
				"-r Minimum T/N AF ratio, defaults to 0.\n"+
				"-p Don't remove non StrongSomatic or LikelySomatic.\n"+
				"-s Print spreadsheet variant summary.\n"+
				"-f Directory in which to save the parsed files, defaults to the parent dir of the vcfs.\n"+

				"\nExample: java -jar pathToUSeq/Apps/VarDictVCFParser -v /VCFFiles/ -t 0.03 -n 0.6 \n"+
				"-u 30 -o 10 -a 3 -d 0.03 -r 2 -e 1 \n\n"+


				"**************************************************************************************\n");

	}

}
