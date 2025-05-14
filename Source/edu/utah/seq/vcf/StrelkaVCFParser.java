package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Simple Strelka INDEL vcf formatter and parser.
 * Takes the tumor QSI or QSS score and replaces the qual. 
 */
public class StrelkaVCFParser {

	private File[] vcfFiles;
	private float minimumScore = 0;
	private int minimumTumorReadDepth = 0;
	private int minimumNormalReadDepth = 0;
	private double minimumTNFractionDiff = 0;
	private double minimumTNRatio = 0;
	private double maximumNormalAltFraction = 1;
	private int minimumAltReadDepth = 0;
	private Pattern qsiOrs = Pattern.compile(".+;QS[IS]=(\\d+);.+");
	private Pattern somEVS = Pattern.compile(".+;SomaticEVS=(\\d+).*");
	private Pattern dp = Pattern.compile("DP=\\d+;*");
	private boolean printSpreadsheet = false;
	private String afInfo = "##INFO=<ID=T_AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=T_DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private String nafInfo = "##INFO=<ID=N_AF,Number=1,Type=Float,Description=\"Allele Frequency for normal\">";
	private String ndpInfo = "##INFO=<ID=N_DP,Number=1,Type=Integer,Description=\"Read depth for normal\">";
	private String afFormat = "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">";
	private String gtFormat = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	private double minimumTumorAltFraction = 0;
	private boolean excludeNonPass = false;
	private File saveDirectory = null;
	private int stringency = 0;
	//score points
	private double[] qsiA;
	private double[] qsiB;
	private double[] qsiC;
	private double[] qssA;
	private double[] qssB;
	private double[] qssC;
	
	public StrelkaVCFParser (String[] args) { 

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
		System.out.println(stringency+"\tTiered 100X exome QSI/QSS filtering");
		System.out.println(excludeNonPass+"\tRemove non PASS FILTER field records");
		//System.out.println(replaceQUALWithEVS+"\tReplace QUAL with somaticEVS score.");
		System.out.println(printSpreadsheet+"\tPrint spreadsheet output");
		
		if (stringency != 0) setScoreTiers();
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parseLineByLine(vcf);
		}
		
		System.out.println("\nComplete!");
	}
	
	/*For indels to calc allele freq  TIR/(TAR+TIR) 
	##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2">
	##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">
	 */
	public void setAltRefCountsForIndels(VCFRecord record, VCFSample sample){
		//parse tir and tar
		String tir = sample.getFormatData("TIR");
		if (tir == null) Misc.printErrAndExit("\nError: TIR doesn't appear in this sample?! "+sample.getUnmodifiedSampleString()+"\n"+record.getOriginalRecord());
		String[] tier1_2IndelCounts = Misc.COMMA.split(tir);
		String tar = sample.getFormatData("TAR");
		if (tar == null) Misc.printErrAndExit("\nError: TAR doesn't appear in this sample?! "+sample.getUnmodifiedSampleString()+"\n"+record.getOriginalRecord());
		String[] tier1_2NonIndelCounts = Misc.COMMA.split(tar);
		
		sample.setReferenceCounts(tier1_2NonIndelCounts[0]);
		sample.setAlternateCounts(tier1_2IndelCounts[0]);
		
	}
	
	/*For snvs to calc allele freq  tier1, tier2  DP:FDP:SDP:SUBDP:AU:CU:GU:TU	1655:16:0:0:1638,1686:0,0:0,0:1,2
		##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
		##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
		##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
		##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
	 */
	public static final String[] acgtUs = {"AU", "CU", "GU", "TU"};
	public void setAltRefCountsForSnvs(VCFRecord record, VCFSample sample){
		//parse counts
		int[] acgt = new int[4];
		String ref = record.getReference();
		int nonRefCount = 0;
		int refCount = 0;
		for (int i=0; i< 4; i++){
			String u = sample.getFormatData(acgtUs[i]);
			if (u == null) Misc.printErrAndExit("\nError: "+ acgtUs[i] +" doesn't appear in this sample?! "+sample.getUnmodifiedSampleString());
			String[] tier1_2 = Misc.COMMA.split(u);
			acgt[i] = Integer.parseInt(tier1_2[0]);
			if (acgtUs[i].substring(0,1).equals(ref)) refCount = acgt[i];
			else nonRefCount+= acgt[i];
		}
		sample.setReferenceCounts(refCount+"");
		sample.setAlternateCounts(nonRefCount+"");
	}

	
	public void parseLineByLine(File vcf){
		try {

			VCFParser parser = new VCFParser (vcf, false, true, true);
			if (parser.getSampleNames()[1].equals("TUMOR") == false) Misc.printErrAndExit("Error: TUMOR doesn't appear to be the second sample in the VCF file?! "+vcf.getName());

			//make io
			BufferedReader in = parser.initializeParser();
			File outFile = new File (saveDirectory, Misc.removeExtension(vcf.getName())+"_Filtered.vcf.gz");
			Gzipper passOut = new Gzipper(outFile);
			
			Gzipper outSpreadsheet = null;
			if (printSpreadsheet){
				File txt = new File (saveDirectory, Misc.removeExtension(vcf.getName())+".txt.gz");
				outSpreadsheet = new Gzipper(txt);
				outSpreadsheet.println("#PASS\tCHROM\tPOS\tREF\tALT\tT_AF\tT_DP\tN_AF\tN_DP\tFILTER\tINFO");
			}
			
			//write out header inserting AF
			writeHeaderWithExtraInfo(passOut, parser);
			
			VCFRecord r = null;
			int numPass = 0;
			int numFail = 0;
			
			while ((r= parser.fetchNext(in)) != null){

				//look for . in ref or first alt
				if (r.getReference().equals(".") || r.getAlternate()[0].equals(".")){
					numFail++;
					continue;
				}
				
				VCFSample[] normTum = r.getSample();
				//set allele counts
				if (r.isSNP()) {
					setAltRefCountsForSnvs(r, normTum[0]);
					setAltRefCountsForSnvs(r, normTum[1]);
				}
				else {
					setAltRefCountsForIndels(r, normTum[0]);
					setAltRefCountsForIndels(r, normTum[1]);
				}
				double normRto = normTum[0].getAltRatio();
				double tumRto = normTum[1].getAltRatio();
				
				//check depth
				int normDepth = normTum[0].getReadDepthDP();
				int tumDepth = normTum[1].getReadDepthDP();
				
				if (normDepth < minimumNormalReadDepth || tumDepth < minimumTumorReadDepth) {
					numFail++;
					continue;
				}
				
				//check alt counts
				int altCounts = Integer.parseInt(normTum[1].getAlternateCounts());
				if (altCounts < minimumAltReadDepth) {
					numFail++;
					continue;
				}
				
				//check allelic ratio diff
				if (minimumTNFractionDiff != 0){
					double change = tumRto-normRto;
					if (change < minimumTNFractionDiff) {
						numFail++;
						continue;
					}
				}

				//check T/N AF ratio
				if (minimumTNRatio != 0 && normRto !=0){
					double change = tumRto/normRto;
					if (change < minimumTNRatio) {
						numFail++;
						continue;
					}
				}
				
				//check normal alt fraction?
				if (maximumNormalAltFraction !=1){
					if (normRto > maximumNormalAltFraction) {
						numFail++;
						continue;
					}
				}
				//check tumor alt fraction?
				if (minimumTumorAltFraction !=0){
					if (tumRto < minimumTumorAltFraction) {
						numFail++;
						continue;
					}
				}
				//check PASS?
				if (excludeNonPass && r.getFilter().toLowerCase().contains("pass") == false) {
					numFail++;
					continue;
				}
				
				//set QSI QSS score in QUAL
				Matcher mat = qsiOrs.matcher(r.getOriginalRecord());
				if (mat.matches() == false) Misc.printErrAndExit("QSI or QSS score doesn't appear to be present in the record? "+r.getOriginalRecord());
				float score = Float.parseFloat(mat.group(1));
				r.setQuality(score);
				if (score < minimumScore) {
					numFail++;
					continue;
				}
				
				//apply tuned exome scoring based on QSI and QSS?
				if (stringency !=0) {
					if (tunedExomeFilt(tumRto, r.isSNP(), score) == false) {
						numFail++;
						continue;
					}
				}
				
				if (printSpreadsheet) {
					//build txt output
					ArrayList<String> al = new ArrayList<String>();
					al.add("true");
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
					outSpreadsheet.println(line);
				}
				
				//save it
				numPass++;
				String orig = r.toString();
				//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOR......
				//   0    1   2  3   4   5     6     7      8      9     10
				String[] fields = VCFParser.TAB.split(orig);
				//reset score
				fields[5] = Float.toString(r.getQuality());
				//reset ID
				fields[2] = "Strelka_"+numPass;
				//add GT to format, igv is now requiring this to be first
				fields[8] = "GT:"+fields[8]+ ":AF";
				//add af to Norm and Tum
				fields[9] = "./.:"+ fields[9]+ ":"+ formatAf (r.getSample()[0].getAltRatio());
				String tumorAf = formatAf (r.getSample()[1].getAltRatio());
				fields[10] = "./.:"+ fields[10]+ ":"+ tumorAf;
				//remove existing DP
				fields[7] = dp.matcher(fields[7]).replaceFirst("");
				
				//modify INFO
				String normalAf = formatAf (r.getSample()[0].getAltRatio());
				
				//add DP and AF for tumor to INFO
				fields[7] = "T_DP=" + r.getSample()[1].getReadDepthDP()+ ";T_AF=" + tumorAf+ ";N_DP=" + r.getSample()[0].getReadDepthDP()+ ";N_AF=" + normalAf + ";"+ fields[7] ;
				
				passOut.println(Misc.stringArrayToString(fields, "\t"));

				
			}
			System.out.println(numPass+"\t"+numFail);
			passOut.close();
			if (printSpreadsheet) outSpreadsheet.close();
			in.close();

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing "+vcf);
		}
	}

	public void printRecords(VCFParser parser, File f) throws Exception{
		Gzipper out = new Gzipper (f);
		
		//write out header inserting AF
		writeHeaderWithExtraInfo(out, parser);
		
		int numFail = 0;
		int numPass = 0;
		
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(VCFRecord.FAIL)) numFail++;
			else {
				numPass++;
				String orig = vcf.toString();
				//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOR......
				//   0    1   2  3   4   5     6     7      8      9     10
				String[] fields = VCFParser.TAB.split(orig);
				//reset score
				fields[5] = Float.toString(vcf.getQuality());
				//reset ID
				fields[2] = "Strelka_"+numPass;
				//add GT to format, igv is now requiring this to be first
				fields[8] = "GT:"+fields[8]+ ":AF";
				//add af to Norm and Tum
				fields[9] = "./.:"+ fields[9]+ ":"+ formatAf (vcf.getSample()[0].getAltRatio());
				String tumorAf = formatAf (vcf.getSample()[1].getAltRatio());
				fields[10] = "./.:"+ fields[10]+ ":"+ tumorAf;
				//remove existing DP
				fields[7] = dp.matcher(fields[7]).replaceFirst("");
				
				//modify INFO
				String normalAf = formatAf (vcf.getSample()[0].getAltRatio());
				
				//add DP and AF for tumor to INFO
				fields[7] = "T_DP=" + vcf.getSample()[1].getReadDepthDP()+ ";T_AF=" + tumorAf+ ";N_DP=" + vcf.getSample()[0].getReadDepthDP()+ ";N_AF=" + normalAf + ";"+ fields[7] ;
				
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
		boolean added = false;
		for (String h : parser.getStringComments()) {
			//skip existing DP in info, this is being replaced with dp for tumor
			if (h.startsWith("##INFO=<ID=DP,")) continue;
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
		new StrelkaVCFParser(args);
	}	
	
	public void setScoreTiers(){
			if (stringency == 1){
				qsiA =  Num.calculateSlopeAndYIntercept(0.044,40,0.059,20);
				qsiB =  Num.calculateSlopeAndYIntercept(0.059,20,0.105,15);
				qsiC =  Num.calculateSlopeAndYIntercept(0.105,15,0.207,14);
				qssA =  Num.calculateSlopeAndYIntercept(0.038,34,0.056,19);
				qssB =  Num.calculateSlopeAndYIntercept(0.056,19,0.107,16);
				qssC =  Num.calculateSlopeAndYIntercept(0.107,16,0.213,16);
			}
			else if (stringency == 2){
				qsiA =  Num.calculateSlopeAndYIntercept(0.044,48,0.059,40);
				qsiB =  Num.calculateSlopeAndYIntercept(0.059,40,0.105,31);
				qsiC =  Num.calculateSlopeAndYIntercept(0.105,31,0.207,29);
				qssA =  Num.calculateSlopeAndYIntercept(0.038,50,0.056,34);
				qssB =  Num.calculateSlopeAndYIntercept(0.056,34,0.107,31);
				qssC =  Num.calculateSlopeAndYIntercept(0.107,31,0.213,31);
			}
			else if (stringency == 3){
				qsiA =  Num.calculateSlopeAndYIntercept(0.044,61,0.059,47);
				qsiB =  Num.calculateSlopeAndYIntercept(0.059,47,0.105,45);
				qsiC =  Num.calculateSlopeAndYIntercept(0.105,45,0.207,44);
				qssA =  Num.calculateSlopeAndYIntercept(0.038,67,0.056,49);
				qssB =  Num.calculateSlopeAndYIntercept(0.056,49,0.107,46);
				qssC =  Num.calculateSlopeAndYIntercept(0.107,46,0.213,46);
			}
			else Misc.printErrAndExit("\nError: stringency setting is unknown "+stringency);
		}
	
	/**Applies the score cut off from the NA12878 simulation points using interpolation.*/
	public boolean tunedExomeFilt(double af, boolean snv, double score){
		double threshold;
		if (af <=0.056){
			if (snv) threshold = Num.calculateYGivenX(qssA, af);
			else threshold = Num.calculateYGivenX(qsiA, af);
		}
		else if (af <=0.105){
			if (snv) threshold = Num.calculateYGivenX(qssB, af);
			else threshold = Num.calculateYGivenX(qsiB, af);
		}
		else {
			if (snv) threshold = Num.calculateYGivenX(qssC, af);
			else threshold = Num.calculateYGivenX(qsiC, af);
		}
		//System.out.println("AF:"+af+" Snv:"+snv+" Score:"+score+" Thresh:"+threshold + " Pass:"+(score >= threshold));
		if (score >= threshold) return true;
		return false;
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
					case 'e': stringency = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
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
				"**                           Strelka VCF Parser: April 2025                         **\n" +
				"**************************************************************************************\n" +
				"Parses Strelka VCF INDEL and SNV files, replacing the QUAl score with the QSI or QSS\n"+
				"score. Also filters for read depth, T/N alt allelic ratio and diff, ref/alt with '.',\n"+
				"and tumor and normal alt allelic ratios. Lastly, it inserts the tumor DP and AF info.\n"+
				"For somatic exome datasets sequenced at >100X unique observation read depth, try the \n"+
				"tier filtering to select for lists with approx 1%, 3-5%, and 9-15% FDR. Follow the\n"+
				"example. Not accurate with lower read depts.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-m Minimum QSI or QSS score, defaults to 0.\n"+
				"-e Apply a tuned 100X exome QSI/S stringency tier: 1 (9-15%FDR), 2 (4-6%FDR),\n"+
				"         3 (1-2%FDR), defaults to 0, no tiered filtiering (37-63%FDR).\n"+
				"-t Minimum tumor allele frequency (AF), defaults to 0.\n"+
				"-n Maximum normal AF, defaults to 1.\n"+
				"-u Minimum tumor alignment depth, defaults to 0.\n"+
				"-a Minimum tumor alt count, defaults to 0.\n"+
				"-o Minimum normal alignment depth, defaults to 0.\n"+
				"-d Minimum T-N AF difference, defaults to 0.\n"+
				"-r Minimum T/N AF ratio, defaults to 0.\n"+
				"-p Remove non PASS filter field records.\n"+
				"-s Print spreadsheet variant summary.\n"+
				"-f Directory in which to save the parsed files, defaults to the parent dir of the vcfs.\n"+

				"\nExample: java -jar pathToUSeq/Apps/StrelkaVCFParser -v /VCFFiles/ -t 0.03 -n 0.6 \n"+
				"-u 30 -o 10 -a 3 -d 0.03 -r 2 -e 1 \n\n"+


				"**************************************************************************************\n");

	}

}
