package edu.utah.seq.base;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.parsers.VCFLookUp;
import edu.utah.seq.parsers.VCFParser;
import edu.utah.seq.parsers.VCFRecord;
import edu.utah.seq.useq.data.RegionScoreText;

import util.bio.annotation.Bed;
import util.bio.seq.Seq;
import util.gen.*;

/**Compares variant lists
 * @author Nix
 * */
public class SNPComparator {

	//user fields
	private File sequencingVCFFile;
	private File arraySnpBedFile;
	private float minimumArrayScore = 0.9f;
	private double minimumVCFScore = 0.01;
	private boolean excludeFailingCalls = true;
	private String genomeVersion = "hg19";
	
	private HashMap<String, VCFLookUp> chromVCFRecords;
	private HashMap<String,RegionScoreText[]> chromRegions;
	private VCFParser vcfParser;
	private ArrayList<String> matches = new ArrayList<String>();
	private ArrayList<String> misMatches = new ArrayList<String>();
	private int totalNumberArraySnps = 0;
	private int numberArraySnpsFailingMinimumScore = 0;
	private int numberArraySnpsWithNoVcf = 0;
	private int numberVCFsFailingMinimumScore = 0;
	private int numberSNPMatches = 0;
	private int numberSNPNoMatches;
	private static final Pattern UNDERSCORE = Pattern.compile("_");
	private String url;

	public void printStats(){
		System.out.println();
		System.out.println(vcfParser.getNumberVCFRecords() + "\t# Parsed VCF records");
		System.out.println(totalNumberArraySnps+ "\t# Array SNPs");
		System.out.println(numberArraySnpsFailingMinimumScore+ "\t# Array SNPs failing minimum score "+minimumArrayScore);
		System.out.println(numberArraySnpsWithNoVcf+ "\t# Array SNPs with no VCF record");
		System.out.println(numberVCFsFailingMinimumScore+ "\t# Array SNPs with failing VCF record score "+minimumVCFScore);
		System.out.println(numberSNPMatches+ "\t# matches");
		System.out.println(numberSNPNoMatches+ "\t# mismatches");
	}
	
	//constructor
	public SNPComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}

	//methods
	
	
	private void doWork(){
		//parse vcf file
		System.out.println("Parsing VCF...");
		vcfParser = new VCFParser(sequencingVCFFile, excludeFailingCalls);
		chromVCFRecords = vcfParser.getChromVCFRecords();
		
		//parse bed file of array snp calls
		System.out.println("Parsing BED...");
		chromRegions = Bed.parseBedFile(arraySnpBedFile, true);
		
		//for each chromosome of snp calls (the truth), compare to the seq calls
		System.out.print("Comparing calls for:");
		for (String chr: chromRegions.keySet()){
			if (chromVCFRecords.containsKey(chr)){
				System.out.print(" "+chr);
				compare (chr);
			}
		}
		System.out.println();
		
		printStats();
		
		saveHits();
		
	}
	
	private void saveHits(){
		String name = Misc.removeExtension (sequencingVCFFile.getName()) +"_Int_"+ Misc.removeExtension(arraySnpBedFile.getName());
		final String header = "IGBLink\tchr\tpos\tref\tvcfAlt\tvcfGeno\tvcfScore\taSNPCall\taSNPScore";
			
		if (matches.size()!=0){
			File mFile = new File(sequencingVCFFile.getParentFile(), name+"_Matches.xls");
			matches.add(0, header);
			IO.writeArrayList(matches, mFile);
		}
		
		if (misMatches.size()!=0){
			if (name == null) name = Misc.removeExtension (sequencingVCFFile.getName()) +"_Int_"+ Misc.removeExtension(arraySnpBedFile.getName());
			File mFile = new File(sequencingVCFFile.getParentFile(), name+"_MisMatches.xls");
			misMatches.add(0, header);
			IO.writeArrayList(misMatches, mFile);
		}
		
	}

	private void compare(String chr) {
		RegionScoreText[] regionScoreTexts = chromRegions.get(chr);
		VCFLookUp vcfLookUp = chromVCFRecords.get(chr);
		totalNumberArraySnps += regionScoreTexts.length;
		
		//for each snp
		for (RegionScoreText snp: regionScoreTexts){
			//check score
			if (snp.getScore()< minimumArrayScore){
				numberArraySnpsFailingMinimumScore++;
				continue;
			}
			
			//an sequence calls? should be just one or null
			VCFRecord[] vcf = vcfLookUp.fetchVCFRecords(snp.getStart(), snp.getStop());
			if (vcf == null) {
				numberArraySnpsWithNoVcf++;
				continue;
			}
			//vcfLookUp.fetchVCFRecordsDebug(snp.getStart(), snp.getStop());
			
			//check score
			double score = Double.parseDouble(vcf[0].getSampleScore());

			if (score > minimumVCFScore){
				numberVCFsFailingMinimumScore++;
				continue;
			}
			//compare bases

			if (snpsMatch(vcf[0], snp)){
				//System.err.println("Match");
				numberSNPMatches++;
				matches.add(getHtmlLink(chr, snp.getMiddle(), snp.getText()) +"\t"+chr+"\t"+ vcf[0].toStringSimple()+"\t"+snp.getText()+"\t"+snp.getScore());
			}
			else {
				numberSNPNoMatches++;
				misMatches.add(getHtmlLink(chr, snp.getMiddle(), snp.getText()) +"\t"+chr+"\t"+ vcf[0].toStringSimple()+"\t"+snp.getText()+"\t"+snp.getScore());
				//if (snp.getText().equals("rs1556611_T_C_A_G")){
				//	System.out.println("\n"+chr+"\t"+ vcf[0].toStringSimple()+"\t"+snp.getText()+"\t"+snp.getScore()+"  "+snpsMatchDebug(vcf[0],snp));
				//}
			}
		}
		
	}
	
	public String getHtmlLink(String chr, int position, String name){
			int winStart = position - 101;
			if (winStart < 0) winStart = 0;
			int winEnd = position + 99;
			return url+ chr +"&start="+winStart+"&end="+winEnd+"\",\""+name+"\")";
	}
	
	public boolean snpsMatch (VCFRecord vcf, RegionScoreText snp){

		//get called bases at snp from seq data
		String vcfCalls = vcf.getCalledBases();

		//get snp forward calls microarray, assuming bed name value is rs#######_A1For_A2For_A1Top_A2_Top
		String[] fields = UNDERSCORE.split(snp.getText());
		
		//check forward
		String basesFor = fields[1]+fields[2];
		if (basesFor.equals(vcfCalls)) return true;
		String basesRev = fields[2]+fields[1];
		if (basesRev.equals(vcfCalls)) return true;

		//check revcomp
		basesFor = Seq.reverseComplementDNA(basesFor);
		if (basesFor.equals(vcfCalls)) return true;
		basesRev = Seq.reverseComplementDNA(basesRev);
		if (basesRev.equals(vcfCalls)) return true;

		return false;
	}
	
	public boolean snpsMatchDebug (VCFRecord vcf, RegionScoreText snp){

		//get called bases at snp from seq data
		String vcfCalls = vcf.getCalledBases();
System.out.println("\nvcfCalls"+vcfCalls);
		//get snp forward calls microarray, assuming bed name value is rs#######_A1For_A2For_A1Top_A2_Top
		String[] fields = UNDERSCORE.split(snp.getText());
		
		//check forward
		String basesFor = fields[1]+fields[2];
System.out.println("\t"+basesFor);
		if (basesFor.equals(vcfCalls)) return true;
		String basesRev = fields[2]+fields[1];
System.out.println("\t"+basesRev);
		if (basesRev.equals(vcfCalls)) return true;

		//check revcomp
		basesFor = Seq.reverseComplementDNA(basesFor);
System.out.println("\trc "+basesFor);
		if (basesFor.equals(vcfCalls)) return true;
		basesRev = Seq.reverseComplementDNA(basesRev);
System.out.println("\trc "+basesRev);
		if (basesRev.equals(vcfCalls)) return true;

		return false;
	}
	

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SNPComparator(args);
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
					case 'v': sequencingVCFFile = new File(args[++i]); break;
					case 'b': arraySnpBedFile = new File(args[++i]); break;
					case 's': minimumArrayScore = Float.parseFloat(args[++i]); break;
					case 'g': genomeVersion = args[++i];
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//checkfiles
		if (sequencingVCFFile == null || sequencingVCFFile.canRead() == false) Misc.printExit("\nError: cannot find your sequencing variant file!\n");
		if (arraySnpBedFile == null || arraySnpBedFile.canRead() == false) Misc.printExit("\nError: cannot find your array variant bed file!\n");
		
		//set hotlink
		url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             SNP Comparator : Nov 2012                            **\n" +
				"**************************************************************************************\n" +
				"Beta.\n\n" +

				"Options:\n"+
				"-v VCF file (xxx.vcf(.gz/.zip OK).\n"+
				"-b Bed file containing array snp calls.\n"+
				"-g Genome version, defaults to hg19.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/\n\n"+

		"**************************************************************************************\n");

	}
}
