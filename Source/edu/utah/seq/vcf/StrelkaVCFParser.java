package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
	private Pattern qsiOrs = Pattern.compile(".+;QS[IS]=(\\d+);.+");
	private String afInfo = "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private String dpInfo = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private String afFormat = "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">";
	private String gtFormat = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	private double minimumTumorAltFraction = 0;
	private boolean excludeNonPass = false;
	
	public StrelkaVCFParser (String[] args) { 

		processArgs(args);
		
		System.out.println("Thresholds for Tumor and Normal:");
		System.out.println(minimumTumorReadDepth+"\tMin T alignment depth");
		System.out.println(minimumNormalReadDepth+"\tMin N alignment depth");
		System.out.println(minimumTNFractionDiff+"\tMin T-N allelic fraction diff");
		System.out.println(minimumTNRatio+"\tMin T/N allelic fraction ratio");
		System.out.println(maximumNormalAltFraction+"\tMax N allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin T allelic fraction");
		System.out.println(excludeNonPass+"\tRemove non PASS FILTER field records.");
		
		System.out.println("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
	}
	
	/*For indels to calc allele freq  TIR/(TAR+TIR) 
	##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2">
	##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">
	 */
	public void setAltRefCountsForIndels(VCFRecord record, VCFSample sample){
		//parse tir and tar
		String tir = sample.getFormatData("TIR");
		if (tir == null) Misc.printErrAndExit("\nError: TIR doesn't appear in this sample?! "+sample.getUnmodifiedSampleString());
		String[] tier1_2IndelCounts = Misc.COMMA.split(tir);
		String tar = sample.getFormatData("TAR");
		if (tar == null) Misc.printErrAndExit("\nError: TAR doesn't appear in this sample?! "+sample.getUnmodifiedSampleString());
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

	public void parse(File vcf){
		try {
			VCFParser parser = new VCFParser (vcf);
			if (parser.getSampleNames()[1].equals("TUMOR") == false) Misc.printErrAndExit("Error: TUMOR doesn't appear to be the second sample in the VCF file?! "+vcf.getName());
			
			File txt = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+".txt.gz");
			Gzipper out = new Gzipper(txt);
			out.println("#PASS\tCHROM\tPOS\tREF\tALT\tT_AF\tT_DP\tN_AF\tN_DP\tFILTER\tINFO");
			
			for (VCFRecord r: parser.getVcfRecords()){	
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
				boolean pass = true;
				
				//check depth
				int normDepth = normTum[0].getReadDepthDP();
				int tumDepth = normTum[1].getReadDepthDP();
				
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
				//check PASS?
				if (pass && excludeNonPass && r.getFilter().toLowerCase().contains("pass") == false){
					pass = false;
				}
				
				//set QSI QSS score
				Matcher mat = qsiOrs.matcher(r.getOriginalRecord());
				if (mat.matches() == false) Misc.printErrAndExit("QSI or QSS scored doesn't appear to be present in the record? "+r.getOriginalRecord());
				float score = Float.parseFloat(mat.group(1));
				r.setQuality(score);
				if (score < minimumScore) pass = false;
				
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
				fields[5] = Integer.toString((int)vcf.getQuality());
				//reset ID
				fields[2] = "Strelka_"+numPass;
				//add GT to format, igv is now requiring this to be first
				fields[8] = "GT:"+fields[8]+ ":AF";
				//add af to Norm and Tum
				fields[9] = "./.:"+ fields[9]+ ":"+ formatAf (vcf.getSample()[0].getAltRatio());
				String tumorAf = formatAf (vcf.getSample()[1].getAltRatio());
				fields[10] = "./.:"+ fields[10]+ ":"+ tumorAf;
				//add DP and AF for tumor to INFO
				fields[7] = "DP=" + vcf.getSample()[1].getReadDepthDP()+ ";AF=" + tumorAf + ";"+ fields[7] ;
				
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
			if (afInfo != null && h.startsWith("##FORMAT")) {
				out.println(afInfo);
				out.println(dpInfo);
				out.println(this.gtFormat);
				out.println(afFormat);
				afInfo = null;
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
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'd': minimumTNFractionDiff = Double.parseDouble(args[++i]); break;
					case 'r': minimumTNRatio = Double.parseDouble(args[++i]); break;
					case 'p': excludeNonPass = true; break;
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
				"**                             Strelka VCF Parser: Jan 2017                         **\n" +
				"**************************************************************************************\n" +
				"Parses Strelka VCF INDEL and SNV files, replacing the QUAl score with the QSI or QSS\n"+
				"score. Also filters for read depth, T/N alt allelic ratio and diff,\n"+
				"and tumor and normal alt allelic ratios. Lastly, it inserts the tumor DP and AF info.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-m Minimum QSI or QSS score, defaults to 0.\n"+
				"-t Minimum tumor allele frequency (AF), defaults to 0.\n"+
				"-n Maximum normal AF, defaults to 1.\n"+
				"-u Minimum tumor alignment depth, defaults to 0.\n"+
				"-o Minimum normal alignment depth, defaults to 0.\n"+
				"-d Minimum T-N AF difference, defaults to 0.\n"+
				"-r Minimum T/N AF ratio, defaults to 0.\n"+
				"-p Remove non PASS filter field records.\n"+

				"\nExample: java -jar pathToUSeq/Apps/StrelkaVCFParser -v /VCFFiles/ -m 32 -t 0.05\n"+
				"        -n 0.5 -u 100 -o 20 -d 0.05 -r 2\n\n"+


				"**************************************************************************************\n");

	}

}
