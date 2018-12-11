package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.its.Interval;
import edu.utah.seq.its.IntervalTree;
import edu.utah.seq.useq.data.RegionScoreText;
import edu.utah.seq.useq.data.RegionScoreTextData;
import util.bio.annotation.Bed;
import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;

/**Estimates Fdr for each somatic variant based on a mockTumor vs normal contrast. 
 * @author Nix
 * */
public class VCFFdrEstimator {

	//user fields
	private File mockVcfFile;
	private File[] realTumorVcfsFiles;
	private File saveDirectory;

	private int numSnvsInMock = 0;
	private int numIndelsInMock = 0;
	private float[] mockSnvQuals = null;
	private float[] mockIndelQuals = null;
	private Histogram snvQualHist = null;
	private Histogram indelQualHist = null;
	private VCFParser testParser;

	//constructor
	public VCFFdrEstimator(String[] args) {
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//parse mock and filter
		parseMock();


		//process each test vcf
		for (int i=0; i< realTumorVcfsFiles.length; i++){
			testParser = new VCFParser(realTumorVcfsFiles[i], true, false, false);
			testParser.setRecordQUALAsScore();

			//split into snv and indel
			VCFRecord[][] snvIndel = testParser.splitVCFRecordsBySnvIndelOther();

			//sort snv and indels by score, small to large
			ComparatorVCFRecordScore scr = new ComparatorVCFRecordScore();
			Arrays.sort(snvIndel[0], scr);
			Arrays.sort(snvIndel[1], scr);

			replaceScoreWithFDR(snvIndel[0], mockSnvQuals);
			replaceScoreWithFDR(snvIndel[1], mockIndelQuals);
			
			//watch out for non snvIndel records
			if (snvIndel[2].length !=0) for (VCFRecord v: snvIndel[2]) v.setScore(-1.0f);

			//print modified records adding an FDR=xxx; to the INFO column
			printVariants(Misc.removeExtension(realTumorVcfsFiles[i].getName()));


		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void printVariants(String name) {
		Gzipper out = null;
		try {
			out = new Gzipper( new File(this.saveDirectory, name+"_VFE.vcf.gz"));
			printHeader(out, testParser.getStringComments());
			VCFRecord[] vcfRecords = testParser.getVcfRecords();
			for (VCFRecord v : vcfRecords){
				if (v.getScore() == -1.0f) out.println(v.getOriginalRecord());
				else {
					String[] fields = Misc.TAB.split(v.getOriginalRecord());
					fields[7] = "dFDR="+Num.formatNumber(v.getScore(), 3)+";"+fields[7];
					out.println(Misc.stringArrayToString(fields, "\t"));
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem processing "+name);
		} finally {
			if (out != null) out.closeNoException();
		}
		
	}

	private void printHeader(Gzipper out, String[] header) throws IOException {
		boolean add2Info = true;
		for (int i=0; i< header.length; i++){
			if (add2Info && header[i].startsWith("##INFO")){
				add2Info = false;
				out.println("##INFO=<ID=dFDR,Number=1,Type=Float,Description=\"Lowest estimated FDR observed for the given QUAL or a smaller QUAL score threshold, see USeq VCFFdrEstimator app.\">");
			}
			out.println( header[i]);
		}
		
	}

	private void replaceScoreWithFDR(VCFRecord[] vcfRecords, float[] mockQuals)  {
		double numTest = vcfRecords.length;
		double numMock = mockQuals.length;

		int mockIndex = 0;
		double numTestRemaining;
		float oldScore = vcfRecords[0].getScore();
		double numMockRemaining;
		double fdr;
		double oldFdr = numMock/numTest;
		int i=0;
		float testScore;
		int indexLastSet = 0;
		//for each VCF record
		for (; i< vcfRecords.length; i++){
			testScore =  vcfRecords[i].getScore();

			if (testScore != oldScore){
				numTestRemaining = numTest - i;
				mockIndex = numMockRemaining(mockQuals, mockIndex, testScore);
				//watch out for when there are none left!
				if (mockIndex == numMock){
					vcfRecords[i].setScore(0.0f); 
					for (; i< vcfRecords.length; i++) vcfRecords[i].setScore(0.0f);
					return;
				}
				numMockRemaining = numMock-mockIndex;
				fdr = numMockRemaining/ numTestRemaining;
				if (fdr < oldFdr) oldFdr = fdr;
				oldScore = testScore;
				indexLastSet = i;
			}
			vcfRecords[i].setScore((float)oldFdr);
		}
		//set last?
		if (indexLastSet != vcfRecords.length-1){
			i--;
			testScore =  vcfRecords[i].getScore();
			numTestRemaining = numTest - i;
			mockIndex = numMockRemaining(mockQuals, mockIndex, testScore);
			numMockRemaining = numMock-mockIndex;
			fdr = numMockRemaining/ numTestRemaining;
			if (fdr < oldFdr) oldFdr = fdr;
			vcfRecords[i].setScore((float)oldFdr);
		}
	}

	private int numMockRemaining(float[] mockQuals, int mockIndex, float score) {
		int i = mockIndex;
		for (; i < mockQuals.length; i++){
			if (mockQuals[i] >= score) break;
		}
		return i;
	}

	private void parseMock() {
		IO.pl("Mock vcf QUAL score stats:");
		ArrayList<Float> snvQuals = new ArrayList<Float>();
		ArrayList<Float> indelQuals = new ArrayList<Float>();


		VCFParser parser = new VCFParser(mockVcfFile, true, false, false);	
		//for each record
		float maxIndel = 0;
		float maxSnv = 0;  
		for (VCFRecord rec: parser.getVcfRecords()){
			float q = rec.getQuality();
			if (rec.isSNP()){
				numSnvsInMock++;
				snvQuals.add(q);
				if (q > maxSnv) maxSnv = q;
			}
			else if (rec.isDeletion() || rec.isInsertion()){
				numIndelsInMock++;
				indelQuals.add(q);
				if (q > maxIndel) maxIndel = q;
			}
		}

		mockSnvQuals = Num.arrayListOfFloatToArray(snvQuals);
		Arrays.sort(mockSnvQuals);
		mockIndelQuals = Num.arrayListOfFloatToArray(indelQuals);
		Arrays.sort(mockIndelQuals);

		snvQualHist = new Histogram(0, maxSnv*1.1, 25);
		snvQualHist.countAll(mockSnvQuals);
		indelQualHist = new Histogram(0, maxIndel*1.1, 25);
		indelQualHist.countAll(mockIndelQuals);	

		IO.pl("\t"+ numSnvsInMock+"\tNum SNVs");
		IO.pl("\t"+ numIndelsInMock+"\tNum INDELs");
		IO.pl("\nHistogram of SNV QUALs:");
		snvQualHist.printScaledHistogram();
		IO.pl("\nHistogram of INDEL QUALs:");
		indelQualHist.printScaledHistogram();
		IO.pl();
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFFdrEstimator(args);
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
					case 'm': mockVcfFile = new File(args[++i]); break;
					case 'v': forExtraction = new File(args[++i]); break;
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
		//checkfiles
		if (mockVcfFile == null || mockVcfFile.exists() == false) Misc.printErrAndExit("\nError: please provide either a mock somatic variant vcf to use to count the number of false positives.\n");

		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a vcf file or directory containing such to estimate FDRs using the mock.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		realTumorVcfsFiles = IO.collapseFileArray(tot);
		if (realTumorVcfsFiles == null || realTumorVcfsFiles.length ==0 || realTumorVcfsFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your somatic xxx.vcf(.zip/.gz) file(s)!\n");


		if (saveDirectory != null){
			saveDirectory.mkdirs();
			if (saveDirectory.isDirectory() == false || saveDirectory.exists() == false) Misc.printErrAndExit("\nCannot find or make your save directory?! "+saveDirectory);
		}
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Comparator : March 2017                         **\n" +
				"**************************************************************************************\n" +
				"Compares test vcf file(s) against a gold standard key of trusted vcf calls. Only calls\n" +
				"that fall in the common interrogated regions are compared. WARNING tabix gzipped files\n" +
				"often fail to parse correctly with java. Seeing odd error messages? Try uncompressing.\n"+
				"Be sure a score is provided in the QUAL field.\n\n" +

				"Required Options:\n"+
				"-a VCF file for the key dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-b Bed file of interrogated regions for the key dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-c VCF file for the test dataset (xxx.vcf(.gz/.zip OK)). May also provide a directory\n" +
				"       containing xxx.vcf(.gz/.zip OK) files to compare.\n"+
				"-d Bed file of interrogated regions for the test dataset (xxx.bed(.gz/.zip OK)).\n"+

				"\nOptional Options:\n"+
				"-k Use a bed file of approx key variants (chr start stop type[#alt_#ref_SNV/INS/DEL]\n"+
				"       instead of a vcf key.\n"+
				"-g Require the genotype to match, defaults to scoring a match when the alternate\n" +
				"       allele is present.\n"+
				"-f Only require the position to match, don't consider the alt base or genotype.\n"+
				"-v Use VQSLOD score as ranking statistic in place of the QUAL score.\n"+
				"-s Only compare SNPs, defaults to all.\n"+
				"-n Only compare non SNPs, defaults to all.\n"+
				"-p Provide a full path directory for saving the parsed data. Defaults to not saving.\n"+
				"-e Exclude test and key records whose FILTER field is not . or PASS. Defaults to\n" +
				"       scoring all.\n"+
				"-i Relax matches to key INDELs to include all test variants within x bps.\n"+

				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/VCFComparator -a /NIST/NA12878/key.vcf\n" +
				"       -b /NIST/NA12878/regions.bed.gz -c /EdgeBio/Exome/testHaploCaller.vcf.zip\n" +
				"       -d /EdgeBio/Exome/NimbleGenExomeV3.bed -g -v -s -e -p /CompRes/ \n\n"+

				"**************************************************************************************\n");

	}
}
