package edu.utah.seq.vcf;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.vcf.fdr.VCFFdrEstimator;
import util.gen.*;

/**For SnpEff annotated vcf files, counts the number of affected genes and compares their counts to those in another group using a fisher's exact test.
 * @author Nix
 * */
public class VCFGeneFreqComparator {

	//user fields
	private File[] vcfA;
	private ParsedVcfStats[] pvsA = null;
	private HashMap<String, Integer> geneNameCountsA = null;
	
	private File[] vcfB;
	private ParsedVcfStats[] pvsB = null;
	private HashMap<String, Integer> geneNameCountsB = null;
	private Pattern geneInfo = Pattern.compile(".+;GENEINFO=([a-zA-Z_0-9\\-\\.]+):.+");
	private int maxGeneCount = 0;
	
	//constructor
	public VCFGeneFreqComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//parse vcfs collecting stats on each
		IO.pl("Loading group A (file #vars #genes)...");
		pvsA = loadVcfs(vcfA);
		geneNameCountsA = mergeCountsPrintStats(pvsA);
		
		IO.pl("\nLoading group B (file #vars #genes)...");
		pvsB = loadVcfs(vcfB);
		geneNameCountsB = mergeCountsPrintStats(pvsB);
		
		IO.pl("\nComparing groups...");
		compareGroups();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void compareGroups() {
		TreeSet<String> allGenes = new TreeSet<String>();
		allGenes.addAll(geneNameCountsA.keySet());
		allGenes.addAll(geneNameCountsB.keySet());
		double numADatasets = pvsA.length;
		double numBDatasets = pvsB.length;
		double numAVars = countVars(pvsA);
		double numBVars = countVars(pvsB);
		
		IO.pl("# A Datasets\t"+ (int) numADatasets);
		IO.pl("# B Datasets\t"+ (int) numBDatasets);
		IO.pl("# A Varaints\t"+ (int) numAVars);
		IO.pl("# B Varaints\t"+ (int) numBVars);
		
		FisherExact fe = new FisherExact((int)((2*maxGeneCount) + numAVars+ numBVars));
		
		//for each gene
		IO.pl("\nGene\t#A\t#B\tASetFreq\tAVarFreq\tBSetFreq\tBVarFreq\tSeqPVal\tVarPVal");
		for (String gene: allGenes) {
			StringBuilder sb = new StringBuilder(gene); sb.append("\t");
			Integer aCount = geneNameCountsA.get(gene);
			if (aCount == null) aCount = new Integer(0);
			Integer bCount = geneNameCountsB.get(gene);
			if (bCount == null) bCount = new Integer(0);
			sb.append(aCount); sb.append("\t");
			sb.append(bCount); sb.append("\t");
			
			double freqADatasets = 0;
			double freqAVars = 0;
			if (aCount !=0) {
				double count = aCount.doubleValue();
				freqADatasets = count/numADatasets;
				freqAVars = count/numAVars;
				sb.append(Num.formatNumber(freqADatasets, 4)); sb.append("\t");
				sb.append(Num.formatNumber(freqAVars, 4)); sb.append("\t");
			}
			else sb.append("0\t0\t");
			
			double freqBDatasets = 0;
			double freqBVars = 0;
			if (bCount !=0) {
				double count = bCount.doubleValue();
				freqBDatasets = count/numBDatasets;
				freqBVars = count/numBVars;
				sb.append( Num.formatNumber(freqBDatasets, 4)); sb.append("\t");
				sb.append( Num.formatNumber(freqBVars, 4)); sb.append("\t");
			}
			else sb.append("0\t0\t");
			
			
			double setP = fe.getTwoTailedP(aCount, (int)numADatasets, bCount, (int)numBDatasets);
			double varP = fe.getTwoTailedP(aCount, (int)numAVars, bCount, (int)numBVars);
			sb.append(Num.formatNumber(setP, 4)); sb.append("\t");
			sb.append(Num.formatNumber(varP, 4));

			IO.pl(sb.toString());
		}
	}

	private int countVars(ParsedVcfStats[] pvss) {
		int num = 0;
		for (ParsedVcfStats p: pvss) num+=p.numVars;
		return num;
	}

	private HashMap<String, Integer> mergeCountsPrintStats(ParsedVcfStats[] groupA) {
		HashMap<String, Integer> geneNameCounts = new HashMap<String, Integer>();
		double totalNumVars = 0;
		for (int i=0; i<groupA.length; i++) {
			ParsedVcfStats pvs = groupA[i];
			totalNumVars+= (double) pvs.numVars;
			IO.pl(Misc.removeExtension(pvs.vcfFile.getName())+"\t"+pvs.numVars+"\t"+pvs.geneNameCounts.size());
			for (String gn: pvs.geneNameCounts.keySet()) {
				int numHits = pvs.geneNameCounts.get(gn);
				Integer numTotalHits = geneNameCounts.get(gn);
				if (numTotalHits != null) numHits += numTotalHits;
				geneNameCounts.put(gn, new Integer(numHits));
				if (numHits > maxGeneCount) maxGeneCount = numHits;
			}
		}
		double aveVarsPerFile = totalNumVars/(double)groupA.length;
		IO.pl("\tMean # vars\t"+ Num.formatNumber(aveVarsPerFile, 2));
		return geneNameCounts;
	}

	private class ParsedVcfStats {
		File vcfFile = null;
		int numVars = 0;
		HashMap<String, Integer> geneNameCounts = null;
		
		ParsedVcfStats (File vcfFile, int numVars, HashMap<String, Integer> geneNameCounts) {
			this.vcfFile = vcfFile;
			this.numVars = numVars;
			this.geneNameCounts = geneNameCounts;
		}
	}
	
	private ParsedVcfStats[] loadVcfs(File[] vcfs) {
		ParsedVcfStats[] pvs = new ParsedVcfStats[vcfs.length];
		for (int i=0; i< vcfs.length; i++) pvs[i] = countVcfFile(vcfs[i]);
		return pvs;
	}

	private ParsedVcfStats countVcfFile (File vcf) {
		BufferedReader in = null;
		String line = null;
		HashMap<String, Integer> geneNameCounts = new HashMap<String, Integer>();
		int numVars = 0;
		try {
			in = IO.fetchBufferedReader(vcf);
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
					numVars++;
					String[] fields = Misc.TAB.split(line);
					
					//attempt to pull the clinvar gene for patho
					String gene = fetchClinvarPathoGene(fields[7]);
					if (gene == null) {
						//not found so find most freq affected
						String[] splitInfo = Misc.SEMI_COLON.split(fields[7]);
						LinkedHashMap<String, Integer> geneNameCount = fetchAffectedGenes (splitInfo);
						gene = findFirstLargest( geneNameCount);
					}
					Integer i = geneNameCounts.get(gene);
					if (i == null) i = new Integer(1);
					else i = new Integer (i.intValue()+1);
					geneNameCounts.put(gene,i);
				}
			}
			return new ParsedVcfStats(vcf, numVars, geneNameCounts);

		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing gene name from "+vcf+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
		return null;
	}
	
	private String fetchClinvarPathoGene(String info) {
		//is there a CLNSIG=Pathogenic or CLNSIG=Likely_pathogenic? Then return the associated gene
		if (info.contains("CLNSIG=Pathogenic") || info.contains("CLNSIG=Likely_pathogenic")) {
			Matcher mat = geneInfo.matcher(info);
			if (mat.matches()) return mat.group(1);
		}
		return null;
	}

	/**Walks the LHM finding the key with the max counts. For those with the same # counts, the first is returned.*/
	public static String findFirstLargest(LinkedHashMap<String, Integer> nameCount) {
		int max = 0;
		String name = null;
		for (String s: nameCount.keySet()) {
			int count = nameCount.get(s);
			if (count > max) {
				max = count;
				name = s;
			}
		}
		return name;
	}
	
	/**
	 * ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
	 * 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
	 *    0          1                2               3          4            5
	 * ANN=GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20),
	 * GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20)
	 * 
	 * */
	public static LinkedHashMap<String, Integer> fetchAffectedGenes(String[] splitInfo) throws IOException {
		//find ANN=
		String annValue = null;
		for (String i: splitInfo){
			if (i.startsWith("ANN=")) {
				annValue = i;
				break;
			}
		}
		if (annValue == null) throw new IOException("Failed to find the ANN= in "+Misc.stringArrayToString(splitInfo, " "));

		//drop ANN=
		annValue = annValue.substring(4);

		//split ANN by comma
		String[] splitAnns = Misc.COMMA.split(annValue);

		//Save observed genes, but don't double count!
		LinkedHashMap<String, Integer> geneNameHits = new LinkedHashMap<String, Integer>();
		for (String ann: splitAnns){
			//split ANN and save gene name
			String[] splitAnn = Misc.PIPE.split(ann);
			Integer i = geneNameHits.get(splitAnn[3]);
			if (i == null) i = new Integer(1);
			else i = new Integer (i.intValue()+1);
			geneNameHits.put(splitAnn[3],i);
		}
		
		return geneNameHits;
		
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFGeneFreqComparator(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtractionA = null;
		File forExtractionB = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': forExtractionA = new File(args[++i]); break;
					case 'b': forExtractionB = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		vcfA = VCFFdrEstimator.fetchVcfFiles(forExtractionA);
		vcfB = VCFFdrEstimator.fetchVcfFiles(forExtractionB);
	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            VCF Gene Freq Comparator : Aug 2020                   **\n" +
				"**************************************************************************************\n" +
				"For SnpEff annotated vcf files, counts the number of affected genes and compares their\n"+
				"counts to those in another group using a two-tailed fisher's exact test. Pre filter\n" +
				"vcfs to include relevant variants.\n"+

				"\nRequired Params:\n"+
				"-a Directory containing filtered vcf files for group A (xxx.vcf(.gz/.zip OK))\n"+
				"-a Directory containing group B vcf files \n"+

				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VCFGeneFreqComparator -a MeningiomaVcfs/\n" +
				"       -b AllCancerVcfs/ &> vgfc.results.txt \n"+

		"\n**************************************************************************************\n");

	}
}
