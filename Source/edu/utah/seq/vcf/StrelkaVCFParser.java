package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Simple Strelka INDEL vcf formatter and parser.
 * Takes the tumor QSI or QSS score and replaces the qual. 
 * Can also filter for minimum read depth count in Tum and Normal, as well as minimum Alt Allele fraction change. */
public class StrelkaVCFParser {

	private File[] vcfFiles;
	private float minimumScore = 0;
	private int minimumCount = 0;
	private double minimumAltTNRatio = 0;
	private double maximumNormalAltFraction = 1;
	private Pattern qsiOrs = Pattern.compile(".+;QS[IS]=(\\d+);.+");
	private String afInfo = "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">";
	private double minimumTumorAltFraction = 0;
	
	public StrelkaVCFParser (String[] args) {

		processArgs(args);
		
		System.out.println("Thresholds:");
		System.out.println(minimumScore+"\tMin QSI/ QSS score");
		System.out.println(minimumCount+"\tMin alignment depth");
		System.out.println(minimumAltTNRatio+"\tMin Allelic fraction change");
		System.out.println(maximumNormalAltFraction+"\tMax Normal alt allelic fraction");
		System.out.println(minimumTumorAltFraction+"\tMin Tumor alt allelic fraction");
		
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
			//set all to pass
			parser.setFilterFieldOnAllRecords(VCFRecord.PASS);
			
			for (VCFRecord r: parser.getVcfRecords()){	
				VCFSample[] normTum = r.getSample();
				
				//check depth
				int normDepth = normTum[0].getReadDepthDP();
				int tumDepth = normTum[1].getReadDepthDP();
				if (normDepth < minimumCount || tumDepth < minimumCount) {
					r.setFilter(VCFRecord.FAIL);
					continue;
				}
				
				//set allele counts
				if (r.isSNP()) {
					setAltRefCountsForSnvs(r, normTum[0]);
					setAltRefCountsForSnvs(r, normTum[1]);
				}
				else {
					setAltRefCountsForIndels(r, normTum[0]);
					setAltRefCountsForIndels(r, normTum[1]);
				}
				
				//check alt allelic ratio 
				if (minimumAltTNRatio != 0.0){
					double normRto = normTum[0].getAltRatio();
					if (normRto == 0) normRto = 0.001;
					double tumRto = normTum[1].getAltRatio();
					if (tumRto == 0) tumRto = 0.001;
					
					double rto = tumRto/normRto;
					
					if (rto < minimumAltTNRatio){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				
				//check normal alt fraction?
				if (maximumNormalAltFraction !=1){
					double normRto = normTum[0].getAltRatio();
					if (normRto > maximumNormalAltFraction){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				
				//check tumor alt fraction?
				if (minimumTumorAltFraction !=0){
					double tumAF = normTum[1].getAltRatio();
					if (tumAF < minimumTumorAltFraction){
						r.setFilter(VCFRecord.FAIL);
						continue;
					}
				}
				
				//set QSI QSS score
				Matcher mat = qsiOrs.matcher(r.getOriginalRecord());
				if (mat.matches() == false) Misc.printErrAndExit("QSI or QSS scored doesn't appear to be present in the record? "+r.getOriginalRecord());
				float score = Float.parseFloat(mat.group(1));
				r.setQuality(score);
				if (score < minimumScore) r.setFilter(VCFRecord.FAIL);
			}
			File out = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"_Filtered.vcf.gz");
			printRecords(parser, out);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void printRecords(VCFParser parser, File f) throws Exception{
		Gzipper out = new Gzipper (f);
		
		//write out header inserting AF
		writeHeaderWithAF(out, parser);
		
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
				//reset FILTER
				fields[6] = ".";
				//add af to format
				fields[8] = fields[8]+ ":AF";
				//add af to Norm and Tum
				fields[9] = fields[9]+ formatAf (vcf.getSample()[0].getAltRatio());
				fields[10] = fields[10]+ formatAf (vcf.getSample()[1].getAltRatio());
				
				out.println(Misc.stringArrayToString(fields, "\t"));
			}
		}
		out.close();
		System.out.println(numPass+"\t"+numFail);
	}
	
	private String formatAf(double af){
		if (af == 0.0) return ":0";
		return ":" + Num.formatNumber(af, 3);
	}

	private void writeHeaderWithAF(Gzipper out, VCFParser parser) throws IOException {
		for (String h : parser.getStringComments()) {
			if (afInfo != null && h.startsWith("##INFO")) {
				out.println(afInfo);
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
					case 'a': minimumCount = Integer.parseInt(args[++i]); break;
					case 'r': minimumAltTNRatio = Double.parseDouble(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
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
				"**                             Strelka VCF Parser: Dec 2015                         **\n" +
				"**************************************************************************************\n" +
				"Parses Strelka VCF INDEL and SNV files, replacing the QUAl score with the QSI or QSS\n"+
				"score. Also filters for minimum tumor normal read depth, T/N alt allelic ratio,\n"+
				"and tumor and normal alt allelic ratios. Lastly, it sets the FILTER field to '.' and\n"+
				"inserts the AF allele frequency.\n"+

				"\nRequired Params:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-m Minimum QSI or QSS score, defaults to 0.\n"+
				"-a Minimum alignment depth for both tumor and normal samples, defaults to 0.\n"+
				"-r Minimum alt Tum/Norm allelic ratio, defaults to 0.\n"+
				"-t Minimum tumor alt allelic fraction, defaults to 0.\n"+
				"-n Maximum normal alt allelic fraction, defaults to 1.\n"+

				"\nExample: java -jar pathToUSeq/Apps/StrelkaVCFParser -v /VCFFiles/ -m 32 -a 150\n"+
				"      -r 2 -n 0.6 -t 0.025  \n\n" +


				"**************************************************************************************\n");

	}

}
