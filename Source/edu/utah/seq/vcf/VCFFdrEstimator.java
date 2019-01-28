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

/**Estimates Fdr for each somatic variant based on a backgroundTumor vs normal contrast. 
 * @author Nix
 * */
public class VCFFdrEstimator {

	//user fields
	private File bkgVcfFile;
	private File realTumorVcfsFile;
	private File fdrEstimatesVcfFile;

	private VCFRecord[][] snvIndel = null;
	private int numSnvsInBkg = 0;
	private int numIndelsInBkg = 0;
	private int numSnvsInReal = 0;
	private int numIndelsInReal = 0;
	private float[] backgroundSnvQuals = null;
	private float[] backgroundIndelQuals = null;
	private Histogram snvQualHistBkg = null;
	private Histogram indelQualHistBkg = null;
	private Histogram snvQualHistReal = null;
	private Histogram indelQualHistReal = null;
	private VCFParser testParser;

	//constructor
	public VCFFdrEstimator(String[] args) throws FileNotFoundException, IOException {
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//parse background
		parseBkg();

		parseReal();
		
		//print stats
		printStats();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void parseReal() {
		testParser = new VCFParser(realTumorVcfsFile, true, false, false);
		testParser.setRecordQUALAsScore();

		//split into snv and indel
		snvIndel = testParser.splitVCFRecordsBySnvIndelOther();

		//sort snv and indels by score, small to large
		ComparatorVCFRecordScore scr = new ComparatorVCFRecordScore();
		Arrays.sort(snvIndel[0], scr);
		Arrays.sort(snvIndel[1], scr);
		numSnvsInReal = snvIndel[0].length;
		numIndelsInReal = snvIndel[1].length;
		for (VCFRecord v: snvIndel[0]) snvQualHistReal.count(v.getScore());
		for (VCFRecord v: snvIndel[1]) indelQualHistReal.count(v.getScore());

		replaceScoreWithFDR(snvIndel[0], backgroundSnvQuals);
		replaceScoreWithFDR(snvIndel[1], backgroundIndelQuals);
		
		//watch out for non snvIndel records
		if (snvIndel[2].length !=0) for (VCFRecord v: snvIndel[2]) v.setScore(-1.0f);

		//print modified records adding an dFDR=xxx; to the INFO column
		printVariants(Misc.removeExtension(realTumorVcfsFile.getName()));
		
	}

	private void printStats() {
		IO.pl("QUAL score stats:");
		
		IO.pl("\t"+ numSnvsInReal+"\tNum Som SNVs");
		IO.pl("\t"+ numSnvsInBkg+"\tNum Bkg SNVs");
		
		IO.pl("\t"+ numIndelsInReal+"\tNum Som INDELs");
		IO.pl("\t"+ numIndelsInBkg+"\tNum Bkg INDELs");
		
		IO.pl("\nHistogram of Som SNV QUALs:");
		snvQualHistReal.printScaledHistogram();
		IO.pl("\nHistogram of Bkg SNV QUALs:");
		snvQualHistBkg.printScaledHistogram();
		
		IO.pl("\nHistogram of Som INDEL QUALs:");
		indelQualHistReal.printScaledHistogram();
		IO.pl("\nHistogram of Bkg INDEL QUALs:");
		indelQualHistBkg.printScaledHistogram();
		
		//print list FDRs
		IO.pl("\nEstimated dFDRs:");
		IO.pl("SNV Min "+snvIndel[0][numSnvsInReal-1].getScore());
		IO.pl("SNV Max "+snvIndel[0][0].getScore());
		IO.pl("INDEL Min "+snvIndel[1][numIndelsInReal-1].getScore());
		IO.pl("INDEL Max "+snvIndel[1][0].getScore());
		IO.pl();
		
	}

	private void printVariants(String name) {
		Gzipper out = null;
		try {
			out = new Gzipper( fdrEstimatesVcfFile);
			printHeader(out, testParser.getStringComments());
			VCFRecord[] vcfRecords = testParser.getVcfRecords();
			for (VCFRecord v : vcfRecords){
				if (v.getScore() == -1.0f) out.println(v.getOriginalRecord());
				else {
					String[] fields = Misc.TAB.split(v.getOriginalRecord());
					fields[7] = "dFDR="+Num.formatNumber(v.getScore(), 5)+";"+fields[7];
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

	private void replaceScoreWithFDR(VCFRecord[] vcfRecords, float[] backgroundQuals)  {
		double numTest = vcfRecords.length;
		double numBkg = backgroundQuals.length;

		int backgroundIndex = 0;
		double numTestRemaining;
		float oldScore = vcfRecords[0].getScore();
		double numBkgRemaining;
		double fdr;
		double oldFdr = numBkg/numTest;
		int i=0;
		float testScore;
		int indexLastSet = 0;

		//for each VCF record
		for (; i< vcfRecords.length; i++){
			testScore =  vcfRecords[i].getScore();

			if (testScore != oldScore){
				numTestRemaining = numTest - i;
				backgroundIndex = numBkgRemaining(backgroundQuals, backgroundIndex, testScore);
				//watch out for when there are none left!
				if (backgroundIndex == numBkg){
					vcfRecords[i].setScore(0.0f); 
					for (; i< vcfRecords.length; i++) vcfRecords[i].setScore(0.0f);
					return;
				}
				numBkgRemaining = numBkg-backgroundIndex;
				fdr = numBkgRemaining/ numTestRemaining;

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
			backgroundIndex = numBkgRemaining(backgroundQuals, backgroundIndex, testScore);
			numBkgRemaining = numBkg-backgroundIndex;
			fdr = numBkgRemaining/ numTestRemaining;
			if (fdr < oldFdr) oldFdr = fdr;
			vcfRecords[i].setScore((float)oldFdr);
		}
	}

	private int numBkgRemaining(float[] backgroundQuals, int backgroundIndex, float score) {
		int i = backgroundIndex;
		for (; i < backgroundQuals.length; i++){
			if (backgroundQuals[i] >= score) break;
		}
		return i;
	}

	private void parseBkg() {
		
		ArrayList<Float> snvQuals = new ArrayList<Float>();
		ArrayList<Float> indelQuals = new ArrayList<Float>();


		VCFParser parser = new VCFParser(bkgVcfFile, true, false, false);	
		//for each record
		float maxIndel = 0;
		float maxSnv = 0;  
		for (VCFRecord rec: parser.getVcfRecords()){
			float q = rec.getQuality();
			if (rec.isSNP()){
				numSnvsInBkg++;
				snvQuals.add(q);
				if (q > maxSnv) maxSnv = q;
			}
			else if (rec.isDeletion() || rec.isInsertion()){
				numIndelsInBkg++;
				indelQuals.add(q);
				if (q > maxIndel) maxIndel = q;
			}
		}

		backgroundSnvQuals = Num.arrayListOfFloatToArray(snvQuals);
		Arrays.sort(backgroundSnvQuals);

		backgroundIndelQuals = Num.arrayListOfFloatToArray(indelQuals);
		Arrays.sort(backgroundIndelQuals);

		snvQualHistBkg = new Histogram(0, maxSnv*1.1, 25);
		snvQualHistBkg.countAll(backgroundSnvQuals);
		indelQualHistBkg = new Histogram(0, maxIndel*1.1, 25);
		indelQualHistBkg.countAll(backgroundIndelQuals);
		
		snvQualHistReal = new Histogram(0, maxSnv*1.1, 25);
		indelQualHistReal = new Histogram(0, maxIndel*1.1, 25);

	}



	public static void main(String[] args) throws FileNotFoundException, IOException {
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
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bkgVcfFile = new File(args[++i]); break;
					case 's': realTumorVcfsFile = new File(args[++i]); break;
					case 'r': fdrEstimatesVcfFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check files
		if (bkgVcfFile == null || bkgVcfFile.exists() == false) Misc.printErrAndExit("\nError: please provide either a background variant VCF file to use in estimating FDRs.\n");
		if (realTumorVcfsFile == null || realTumorVcfsFile.exists() == false) Misc.printErrAndExit("\nError: please provide a somatic variant VCF file to use in calculating dFDRs.\n");
		if (fdrEstimatesVcfFile == null || fdrEstimatesVcfFile.getName().endsWith(".vcf.gz") == false) Misc.printErrAndExit("\nError: please provide a vcf file ending in xxx.vcf.gz to save the results.\n");


	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            VCF Fdr Estimator  :   Jan 2019                       **\n" +
				"**************************************************************************************\n" +
				"Estimates false discovery rates for each QUAL score in a somatic VCF file by counting\n"+
				"the number of records that are >= to that QUAL score in a matched background VCF. The\n"+
				"estimated FDR = #Bkg/ #Som. In cases where increasingly stringent QUAL thresholds\n"+
				"reduce the nummber of Som records but not Bkg records, the FDR increases. To control\n"+
				"for this inconsistancy, the prior FDR is assigned to the more stringent QUAL, a 'dFDR'.\n\n"+
				
				"To generate a matched background VCF file, use the SamReadDepthMatcher app to\n"+
				"subsample a high depth normal bam file to match the read depth over each exon in the\n"+
				"tumor bam file.  Run the same somatic variant calling and filtering workflow used in\n"+
				"generating the real somatic VCF file but substitute the matched depth mock tumor bam\n"+
				"for the real tumor bam. Lastly, use a low stringency set of germline variants\n"+
				"identified in the high depth normal sample to filter out any het and hom variants in\n"+
				"the bkg VCF.\n\n"+

				"Required Options:\n"+
				"-b Background VCF file (xxx.vcf(.gz/.zip OK)).\n"+
				"-s Somatic VCF file (xxx.vcf(.gz/.zip OK)).\n"+
				"-r VCF file for saving the estimated FDR results.\n"+

				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VCFFdrEstimator -b patient123Bkg.vcf.gz\n" +
				"       -v patient123Somatic.vcf.gz -r FinalVcfs/patient123SomaticFdr.vcf.gz\n\n"+

				"**************************************************************************************\n");

	}
}
