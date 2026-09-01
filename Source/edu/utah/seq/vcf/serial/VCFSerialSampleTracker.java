package edu.utah.seq.vcf.serial;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.vcf.fdr.VCFFdrEstimator;
import util.gen.*;

public class VCFSerialSampleTracker {

	//user fields
	private File[] individualCallVcfs = null;
	private File[] allScoredVcfs = null;
	private File highVcf = null;
	private File highMediumVcf = null;
	private File allVcf = null;
	
	private HashMap<String, String> highVcfKeys = new HashMap<String, String>();
	private HashMap<String, String> highMediumVcfKeys = new HashMap<String, String>();
	private HashMap<String, String> allVcfKeys = new HashMap<String, String>();
	
	//constructor
	public VCFSerialSampleTracker(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		processArgs(args);

		//load var type keys
		highVcfKeys = loadVcfKeysInfo(highVcf);
		highMediumVcfKeys = loadVcfKeysInfo(highMediumVcf);
		allVcfKeys = loadVcfKeysInfo(allVcf);
		
		//load the individual sample keys
		ArrayList<HashMap<String, String>> sampleKeys = new ArrayList<HashMap<String, String>>();
		for (File f: individualCallVcfs) sampleKeys.add(loadVcfKeysInfo(f));
		
		//walk each all scored sample
		IO.pl("\nPatientID\tSampleName\tSampleNumber\tVariantID\tVariantConseq\tFoundInSample\tGene\tVAF\tDP");
		for (int i=0; i< allScoredVcfs.length; i++) {
			ArrayList<String> toPrint = new ArrayList<String>();
			//1110039_22712X28_23802X10
			String[] sn = Misc.UNDERSCORE.split(allScoredVcfs[i].getName());
			toPrint.add(sn[0]);
			toPrint.add(sn[1]+"_"+sn[2]);
			toPrint.add((i+1)+"");
			addStats(toPrint, allScoredVcfs[i], sampleKeys.get(i));
			IO.pl(Misc.stringArrayListToString(toPrint, "\t"));
		}
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	
	private static Pattern annPat = Pattern.compile(";ANN=");
	private void addStats(ArrayList<String> toPrint, File allScoredVcf, HashMap<String, String> sampleKeys) {
		//VariantID\tVariantConseq\tFoundInSample\tGene\tVAF
		String prefix = Misc.stringArrayListToString(toPrint, "\t");
		
		BufferedReader in = null;
		String line = null;
		try {
			in = IO.fetchBufferedReader(allScoredVcf);
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					toPrint.clear();
					toPrint.add(prefix);
					
					//#CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT
					//   0       1   2   3   4    5     6     7     8  
					String[] fields = Misc.TAB.split(line);
					//create the key, variantID
					StringBuilder sb = new StringBuilder(fields[0]);
					sb.append("_");
					sb.append(fields[1]);
					sb.append("_");
					sb.append(fields[3]);
					sb.append("_");
					sb.append(fields[4]);
					String key = sb.toString();
					toPrint.add(key);
					
					//Consequence?
					String con = null;
					if (highVcfKeys.containsKey(key)) con = "HIGH";
					else if (highMediumVcfKeys.containsKey(key)) con = "MODERATE";
					else if (allVcfKeys.containsKey(key)) con = "LOW";
					else con = "NA"; // should never see this.
					toPrint.add(con);
					
					//Found in sample, if so parse the gene
					String found = "F";
					String gene = "NA";
					if (sampleKeys.containsKey(key)) {
						found = "T";
						String sampleInfo = sampleKeys.get(key);
						//BKZ=21.57;BKAF=0.0173,0.017,0.017,0.016,0.0154,0.0154,0.0153,0.0148,0.0145,0.0138,0.0128,0.012,0.0067,0.0036,0.0034,0.0025,0.0023,0.0021,0.0016,0.001,0.0009,0,0,0,0,0;T_AF=0.15926;T_DP=923;N_AF=0.04294;N_DP=1211;FPV=194.168;INDEL;IDV=3;IMF=0.00325027;DP=2134;I16=41,1894,19,180,42433,1.39587e+06,2862,46216,114415,6.79773e+06,11940,716400,38661,890699,4954,123766;QS=1.90564,0.0943575;VDB=0.00443749;SGB=-63.0254;MQSB=0.806033;MQ0F=0;ANN=CAA|intron_variant|MODIFIER|ATM|ATM|transcript|NM_000051.4|protein_coding|60/62|c.8786+193_8786+194insAA||||||INFO_REALIGN_3_PRIME
						String[] ann = annPat.split(sampleInfo);
						if (ann.length==2) {
							String[] diffAnn = Misc.COMMA.split(ann[1]);
							TreeSet<String> genes = new TreeSet<String>();
							for (String a: diffAnn) {
								//ANN=CAA|intron_variant|MODIFIER|ATM|ATM|transcript|NM_000051.4|protein_coding|60/62|c.8786+193_8786+194insAA||||||INFO_REALIGN_3_PRIME
								//     0         1           2     3   4      5
								String[] p = Misc.PIPE.split(a);
								if (p.length>5) genes.add(p[3]);
							}
							if (genes.size() !=0) gene = Misc.treeSetToString(genes, ",");
				
						}
						else IO.el("WARNING, failed to find ';ann=' in "+sampleInfo+" with key "+key);
						
					}
					toPrint.add(found);
					toPrint.add(gene);
					
					//Parse AF from the all scored, AF=0.115;DP=1191;.
					String af = "NA";
					String dp = "NA";
					for (String i: Misc.SEMI_COLON.split(fields[7])) {
						if (i.startsWith("AF=")) af = i.substring(3);
						else if (i.startsWith("DP=")) dp = i.substring(3);
					}
					toPrint.add(af);
					toPrint.add(dp);
					IO.pl(Misc.stringArrayListToString(toPrint, "\t"));
				}
			}
		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing "+allScoredVcf+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
	}



	private HashMap<String, String> loadVcfKeysInfo(File vcfFile) {
		BufferedReader in = null;
		String line = null;
		HashMap<String, String> keys = new HashMap<String, String>();
		try {
			in = IO.fetchBufferedReader(vcfFile);
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL FILTER	INFO FORMAT
					//   0       1   2   3   4    5     6     7     8  
					String[] fields = Misc.TAB.split(line);
					StringBuilder sb = new StringBuilder(fields[0]);
					sb.append("_");
					sb.append(fields[1]);
					sb.append("_");
					sb.append(fields[3]);
					sb.append("_");
					sb.append(fields[4]);
					keys.put(sb.toString(), fields[7]);
				}
			}
		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing "+vcfFile+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
		return keys;
	}



	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFSerialSampleTracker(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File indiForExtraction = null;
		File allForExtraction = null;
		
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try {
					switch (test){
					case 'i': indiForExtraction = new File(args[++i]); break;
					case 'a': allForExtraction = new File(args[++i]); break;
					case 'h': highVcf = new File(args[++i]); break;
					case 'm': highMediumVcf = new File(args[++i]); break;
					case 'l': allVcf = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (indiForExtraction == null || allForExtraction == null) {
			printDocs();
			System.exit(1);
		}
		
		individualCallVcfs = VCFFdrEstimator.fetchVcfFiles(indiForExtraction);
		allScoredVcfs = VCFFdrEstimator.fetchVcfFiles(allForExtraction);
		if (individualCallVcfs.length == 0 || (individualCallVcfs.length != allScoredVcfs.length)) {
			Misc.printErrAndExit("\nERROR: the number of vcf files between the two folders differ. These should be paired.\n");
		}
		
		if (highVcf==null || highVcf.exists()==false || highMediumVcf==null || highMediumVcf.exists()==false || allVcf==null || allVcf.exists()==false) {
			Misc.printErrAndExit("\nERROR: one of your -h high, -m medium, or -l low vcfs is null or doesn't exist.\n");
		}
		
		
		//print file names
		IO.pl("Check these are correctly paired:");
		for (int i=0; i< individualCallVcfs.length; i++) {
			IO.pl(individualCallVcfs[i].getName()+" -> "+allScoredVcfs[i].getName());
		}
	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Serial Sample Tracker : May 2026                 **\n" +
				"**************************************************************************************\n" +
				"Generated a variety of statistics to assist in tracking variants over serial sampling.\n"+
				"Run this app on each patient's datasets to report:\n"+
				"  PatientID SampleName SampleNumber VariantID VariantConseq FoundInSample Gene VAF DP\n"+

				"\nRequired Params:\n"+
				"-i Directory containing individual annotated calls for each sample for one patient,\n"+
				"      (xxx.vcf(.gz/.zip OK)), e.g. from TNRunner analysis.\n"+
				"-a Directory containing all calls annotated using the VCFMpileupAnnotator for each\n"+
				"      sample (xxx.vcf(.gz/.zip OK)). Use VCFCollapser on -i to generate the common\n"+
				"      all call vcf.\n"+
				"-h Vcf file containing all of the high consequence variants, use the VCFCollapser to\n"+
				"      generate it.\n"+
				"-m Vcf file containing the high and medium consequence variants.\n"+
				"-l Vcf file containing the high, medium, low, and modifier variants.\n"+
				
				"\nExample: java -jar pathTo/USeq/Apps/VCFSerialSampleTracker -i IndiCalls_1170061\n" +
				"       -a AllCallsScored_1170061 -h high.vcf.gz -m mod.vcf.gz -l low.vcf.gz\n"+

		"\n**************************************************************************************\n");

	}
}
