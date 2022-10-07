package edu.utah.seq.vcf.loh;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.jpileup.BamPileupTabixLoaderSingle;
import edu.utah.seq.parsers.jpileup.BpileupLine;
import util.gen.CombinePValues;
import util.gen.FisherExact;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**
 * @author Nix*/
public class LoH {

	//user defined fields
	private File germlineVcf = null;
	private File germlineBamPileup = null;
	private File somaticBamPileup = null;
	private double minimumQUAL = 20;
	private double minimumGQ = 20;
	private int minimumDP = 20;
	private int maxGapBetweenVars = 1000;
	private int windowBpPadding = 100;
	private File resultsDir = null;

	//window
	private float minAdjTransWindowPval = 8.239f; //.15
	private float minMeanAfWindowDiff = 0.1f;
	//het var
	private double minVarPval = 0.1;
	private double minVarAfDiff = 0.1;
	private String gqSelector = "GQ";

	//internal
	private HashMap<String, String> formatValues = new HashMap<String,String>();
	private HeterozygousVar[] hetVcfRecords = null;
	private LinkedHashMap<String,HeterozygousVar[]> chromHets = null;
	private CombinePValues combinePValues = new CombinePValues();
	private int maxSizeForFishers = 0;
	private HetWindow[] hetWindow = null;
	private String fractionLoH = null;
	private ArrayList<String> vcfHeader = new ArrayList<String>();

	//constructor
	public LoH(String[] args){
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);
			
			IO.pl("Thresholds:\n"+fetchThresholds(""));

			//fetch heterozygous germline variants
			fetchHeterozygousVars();

			//bam pileup counts from the somatic and germline datasets, could thread this!
			addBamPileupObservations();

			printAfStats();

			calculateFisherPValues();

			splitHetsByChrom();

			walkVariants();
			
			adjustPValues();
			
			scoreWindows();

			windowScoreIndividualHets();
			
			printBed();
			
			printVcf();


		} catch (Exception e) {
			e.printStackTrace();
		} finally {

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" min\n");
		}
	}
	
	private void printVcf() throws IOException {
		File vcfFile = new File(resultsDir, "loh.vcf.gz");
		Gzipper out = new Gzipper(vcfFile);
		out.println(Misc.stringArrayListToString(vcfHeader, "\n"));
		
		//for each block
		for (HetWindow w: hetWindow) if (w.isPassesThresholds()) out.print(w.toStringVcfFormat(windowBpPadding));
		
		out.close();
		
	}

	private String fetchThresholds(String prePend) {
		StringBuilder sb = new StringBuilder();
		
		//vcf
		sb.append(prePend); sb.append("minimumVcfQUAL"); sb.append("\t"); sb.append(minimumQUAL); sb.append("\n");
		sb.append(prePend); sb.append("minimumVcfGQ"); sb.append("\t"); sb.append(minimumGQ); sb.append("\n");
		sb.append(prePend); sb.append("minimumVcfDP"); sb.append("\t"); sb.append(minimumDP); sb.append("\n");
		
		//het var
		sb.append(prePend); sb.append("minVarPval"); sb.append("\t"); sb.append(minVarPval); sb.append("\n");
		sb.append(prePend); sb.append("minVarAfDiff"); sb.append("\t"); sb.append(minVarAfDiff); sb.append("\n");
		
		//window
		sb.append(prePend); sb.append("min-10Log10(AdjWindowPval)"); sb.append("\t"); sb.append(minAdjTransWindowPval); sb.append("\n");
		sb.append(prePend); sb.append("minMeanAfWindowDiff"); sb.append("\t"); sb.append(minMeanAfWindowDiff); sb.append("\n");
		sb.append(prePend); sb.append("maxGapBetweenWindowVars"); sb.append("\t"); sb.append(maxGapBetweenVars); sb.append("\n");
		sb.append(prePend); sb.append("bpPaddingForLoHWindow"); sb.append("\t"); sb.append(windowBpPadding); 
		
		return sb.toString();
	}

	private void printBed() throws IOException {
		File bedFile = new File(resultsDir, "loh.bed.gz");
		Gzipper out = new Gzipper(bedFile);
		//fraction LoH
		out.println("# FractionLoH\t"+fractionLoH);
		
		//thresholds
		out.println(fetchThresholds("# "));
		
		//for each window
		out.println("# Chr\tStart\tStop\tLoH_#Vars_MeanAFDiff_-10Log10(AdjWindowPValue)\t-10Log10(WindowPValue)");
		for (HetWindow w: hetWindow) if (w.isPassesThresholds()) out.println(w.toStringBedFormat(windowBpPadding));
		
		out.close();
	}

	private void scoreWindows() throws IOException {
		IO.pl("Scoring windows and vars...\n\nLoH Window Blocks:");
		for (HetWindow w: hetWindow) {
			if (w.getTransAdjPVal()< this.minAdjTransWindowPval || w.getMeanAfDiff() < this.minMeanAfWindowDiff) w.setPassesThresholds(false);
			else {
				w.setPassesThresholds(true);
				IO.pl(w);
			}
		}
	}

	private void windowScoreIndividualHets() throws IOException {
		//IO.pl("LoH vars:");
		//for each of the windows
		for (HetWindow w: hetWindow) {
			
			//for each of the hets, add window scores if better
			//actually now the het only belongs to one window
			for (HeterozygousVar s: w.getHetVars()) {
				if (s.getBestHetWindow() == null) s.setBestHetWindow(w);
				else {
					HetWindow pastBest = s.getBestHetWindow();
					HetWindow currBest = pastBest.compare(w);
					if (currBest!=null) s.setBestHetWindow(currBest);
				}
			}
		}
		int numLoHVars = 0;
		for (HeterozygousVar s: hetVcfRecords) {
			if (s.getBestHetWindow().isPassesThresholds()) {
				numLoHVars++;
				//IO.pl(s.toString());
			}
		}
		
		float frac = (float)numLoHVars/ (float)hetVcfRecords.length;
		fractionLoH = frac+" ("+numLoHVars+"/"+hetVcfRecords.length+")";
		IO.pl(fractionLoH+ " : Fraction heterozygous germline snvs and indels with a significant increase in their allele fraction in the tumor");
		
		
	}

	private void calculateFisherPValues() throws IOException {
		IO.pl("\nCalculating Fisher's Exact PValues...");
		FisherExact fe = new FisherExact(maxSizeForFishers);
		for (HeterozygousVar s: hetVcfRecords) {
			int[] c = s.fetchSomAltRefGermAltRefCounts();
			//fe.getRightTailedP(somAlt, somRef, germAlt, germRef), only interested in a higher somatic AF relative to germline
			double pval = fe.getRightTailedP(c[0], c[1], c[2], c[3]);
			s.setPvalue(pval);
		}
	}

	private void splitHetsByChrom() throws IOException {
		chromHets = new LinkedHashMap<String,HeterozygousVar[]>();
		ArrayList<HeterozygousVar> al = new ArrayList<HeterozygousVar>();
		String currChrom = hetVcfRecords[0].getVcfRecord()[0];

		for (int i=0; i< hetVcfRecords.length; i++){
			if (hetVcfRecords[i].getVcfRecord()[0].equals(currChrom) == false){
				HeterozygousVar[] sub = new HeterozygousVar[al.size()];
				al.toArray(sub);
				chromHets.put(currChrom, sub);
				al.clear();
				currChrom = hetVcfRecords[i].getVcfRecord()[0];
				if (chromHets.containsKey(currChrom)) throw new IOException("\nIs your germline vcf sorted? Seeing this chrom again: "+currChrom);
			}
			al.add(hetVcfRecords[i]);
		}
		//add last to hash
		HeterozygousVar[] sub = new HeterozygousVar[al.size()];
		al.toArray(sub);
		chromHets.put(currChrom, sub);

	}
	
	private void walkVariants() {
		IO.pl("Walking chromosomes...");

		ArrayList<HetWindow> hetWindowAl = new ArrayList<HetWindow>();
		ArrayList<HeterozygousVar> hetAl = new ArrayList<HeterozygousVar>();
		//for each chrom
		for (String chr: chromHets.keySet()) {
			HeterozygousVar[] hets = chromHets.get(chr);
			
			//for each het attempt to make a window of <= maxSize, might be only one het in the window
			for (int i=0; i< hets.length; i++) {
				int startBp = hets[i].getPosition();
				hetAl.clear();
				hetAl.add(hets[i]);
				
				//does that het pass thresholds?
				if (hets[i].getPvalue() <= minVarPval && hets[i].getAlleleFractionDifference() >= minVarAfDiff ) {
					

					//look at next het
					for (int j=i+1; j< hets.length; j++) {				
						//check distance, can't be too far away
						int dist = hets[j].getPosition() - startBp;
						if (dist > maxGapBetweenVars )  {
							break;
						}
						
						//check thresholds
						if (hets[j].getPvalue() > minVarPval || hets[j].getAlleleFractionDifference() < minVarAfDiff ) {	
							break;
						}
						
						//OK, add it
						hetAl.add(hets[j]);
						startBp = hets[j].getPosition();
						i=j;  
					}
				}
				
				
				HeterozygousVar[] winHets = new HeterozygousVar[hetAl.size()];
				hetAl.toArray(winHets);
				
				//calculate a combine pvalue?
				double windowPVal = -1;
				if (winHets.length == 1) windowPVal = winHets[0].getPvalue();
				else {
					double[] pvals = new double[winHets.length];
					double min = 1;
					for (int x=0; x< pvals.length; x++) {
						pvals[x] = winHets[x].getPvalue();
						if (pvals[x]< min) min = pvals[x];
					}
					//might return zero! in that case set it to the smallest
					windowPVal = combinePValues.calculateCombinePValues(pvals);
					if (windowPVal == 0.0) windowPVal = min;
					
				}
				
				HetWindow win = new HetWindow((float)Num.minus10log10(windowPVal), winHets);
				hetWindowAl.add(win);
			}
		}
		
		
		
		hetWindow = new HetWindow[hetWindowAl.size()];
		hetWindowAl.toArray(hetWindow);
		
		
	}

	private void adjustPValues() throws IOException {
		//pull pvals for adjustment
		IO.pl("BH Adjusting PValues...");
		float[] pvals = new float[hetWindow.length];
		for (int i=0; i< pvals.length; i++) pvals[i] = hetWindow[i].getTransPvalue();
		float[] adjPvals = Num.benjaminiHochbergCorrectUnsorted(pvals);
		for (int i=0; i< pvals.length; i++) {
			float ap = adjPvals[i];
			if (ap == -0.0f) ap = 0.0f;
			hetWindow[i].setAdjTransPvalue(ap);
		}
	}

	private void printAfStats() throws IOException {
		float[] germAf = new float[hetVcfRecords.length];
		float[] somAf = new float[hetVcfRecords.length];
		for (int i=0; i< hetVcfRecords.length; i++) {
			germAf[i] = (float)hetVcfRecords[i].getAlleleFractionGermline();
			somAf[i] = (float)hetVcfRecords[i].getAlleleFractionSomatic();
		}
		IO.pl("\nGermline AF:");
		Num.statFloatArray(germAf,false);
		IO.pl("Somatic AF:");
		Num.statFloatArray(somAf,false);
	}

	private void addBamPileupObservations() throws Exception {
		BamPileupTabixLoaderSingle bpg = new BamPileupTabixLoaderSingle(germlineBamPileup, 0);
		BamPileupTabixLoaderSingle bps = new BamPileupTabixLoaderSingle(somaticBamPileup, 0);
		ArrayList<HeterozygousVar> passing = new ArrayList<HeterozygousVar>();

		for (HeterozygousVar hs: hetVcfRecords) {
			ArrayList<BpileupLine> germ = bpg.fetchBpileupRecords(hs.getVcfRecord());
			hs.setGermlineBPs(germ, minimumDP);
			ArrayList<BpileupLine> som = bps.fetchBpileupRecords(hs.getVcfRecord());
			hs.setSomaticBPs(som, minimumDP);
			
			if (hs.isPassing()== false) continue;
			passing.add(hs);
			
			//check total counts
			int totalCounts = hs.getTotalRefAltCount();
			if (totalCounts > maxSizeForFishers) maxSizeForFishers = totalCounts;
		
		}

		hetVcfRecords = new HeterozygousVar[passing.size()];
		passing.toArray(hetVcfRecords);

		//close the readers
		bpg.getTabixReader().close();
		bps.getTabixReader().close();
	}
	
	private void fetchHeterozygousVars() throws Exception {
		ArrayList<HeterozygousVar> hetVcfRec = new ArrayList<HeterozygousVar>();
		int numRecords = 0;
		int numFailingRecords = 0;
		int numPassingIndel = 0;
		int numPassingHom = 0;
		int numPassingHet = 0;
		String line = null;
		BufferedReader in = IO.fetchBufferedReader(germlineVcf);

		//save and modify the header
		boolean added = false;
		//LoHBlock=start-stop,#vars,meanAfDiffSomGerm,pval,adjpval;
		String lohBlock = "##INFO=<ID=LoHBlock,Number=1,Type=String,Description=\"Loss of Heterozygosity block: startBP-stopBP, number het vars in block, "
				+ "mean allele fraction differences between the somatic and germline, Fisher's -10Log10(combine p-value) for the block, Benjamini-Hochberg "
				+ "adjusted window -10Log10(combine pvalue)\">";
		//LoHVar=afS,afG,diff,sAlt/Ref,gAlt/Ref,pval;
		String lohVar = "##INFO=<ID=LoHVar,Number=1,Type=String,Description=\"Loss of Heterozygosity Variant: allele fraction somatic, allele fraction germline, "
				+ "allele fraction difference, somatic ALT/REF counts, germline ALT/REF counts, Fisher's Exact -10log10(p-value)\">";
		while ((line = in.readLine()) != null){
			line = line.trim();
			if (line.length()==0) continue;
			vcfHeader.add(line);
			if (added == false && line.startsWith("##INFO=")) {
				added = true;
				vcfHeader.add(lohBlock);
				vcfHeader.add(lohVar);
			}
			if (line.startsWith("#CHROM")) break;
		}

		//for each data line in the file
		while ((line = in.readLine()) != null){
			numRecords++;

			//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1, Sample2.....
			//   0    1   2  3   4   5     6      7     8      9       10
			String[] fields = Misc.TAB.split(line);

			//check FILTER
			String lc = fields[6].toUpperCase();
			if (lc.contains("PASS") == false && lc.equals(".") == false) {
				numFailingRecords++;
				continue;
			}

			//check whole line QUAL
			double q = Double.parseDouble(fields[5]);
			if (q < minimumQUAL) {
				numFailingRecords++;
				continue;
			}

			//check ref and alt
			if (fields[3].equals("*") || fields[4].equals("*") || fields[4].contains(",")) {
				numFailingRecords++;
				continue;
			}

			//check first sample
			String[] format = Misc.COLON.split(fields[8]);
			Integer hetHom = checkSample(format, fields[9]);
			if (hetHom == null) {
				numFailingRecords++;
				continue;
			}
			
			if (hetHom == 1) numPassingHom++;
			
			else {
				if (fields[3].length()!=1 || fields[4].length()!=1) {
					numPassingIndel++;
					hetVcfRec.add(new HeterozygousVar(fields));
				}
				else {
					numPassingHet++;
					hetVcfRec.add(new HeterozygousVar(fields));
				}
			}
		}
		in.close();
		hetVcfRecords = new HeterozygousVar[hetVcfRec.size()];
		hetVcfRec.toArray(hetVcfRecords);

		IO.pl("\nGermline VCF Parsing Statistics:");
		IO.pl(numRecords + "\t# Records");
		IO.pl(numFailingRecords + "\t# Failing Records");
		IO.pl(numPassingHom + "\t# Passing Hom Vars and Indels");
		IO.pl(numPassingIndel + "\t# Passing Het Indels");
		IO.pl(numPassingHet + "\t# Passing Het Snvs");
	}



	/**Returns null if fails, otherwise returns 0 for het, 1 for hom */
	private Integer checkSample(String[] ids, String sample) throws Exception {

		formatValues.clear();
		String[] values = Misc.COLON.split(sample);

		//GATK often includes a final PS id that isn't populated so load from values
		for (int i=0; i< values.length; i++) formatValues.put(ids[i], values[i]);

		//check genotype
		String gt = formatValues.get("GT");
		if (gt == null) throw new Exception ("\t\tNo GT? "+sample);
		//replace any phasing info
		gt = gt.replace('|', '/');
		//no values?
		if (gt.equals("0/1") == false && gt.equals("1/1") == false  && gt.equals("1/0") == false) {
			return null;
		}

		//check GQ/ GQX
		String gqString = formatValues.get(gqSelector);
		if (gqString == null || gqString.equals(".")) return null;
		else {
			double gq = Double.parseDouble(gqString);
			if (gq < minimumGQ) return null;
		}

		if (gt.equals("0/1") || gt.equals("1/0")) return 0;
		if (gt.equals("1/1")) return 1;
		return null;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new LoH(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': germlineVcf = new File(args[++i]); break;
					case 'g': germlineBamPileup = new File(args[++i]); break;
					case 's': somaticBamPileup = new File(args[++i]); break;
					case 'o': resultsDir = new File(args[++i]); break;
					case 'm': maxGapBetweenVars = Integer.parseInt(args[++i]); break;
					case 'r': minimumDP = Integer.parseInt(args[++i]); break;
					case 'w': windowBpPadding = Integer.parseInt(args[++i]); break;
					case 'q': minimumQUAL = Double.parseDouble(args[++i]); break;
					case 'e': minimumGQ = Double.parseDouble(args[++i]); break;
					case 'x': gqSelector = "GQX"; break;
					case 'p': minAdjTransWindowPval = Float.parseFloat(args[++i]); break;
					case 'd': minMeanAfWindowDiff = Float.parseFloat(args[++i]); break;
					case 'l': minVarPval = Double.parseDouble(args[++i]); break;
					case 'f': minVarAfDiff = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (germlineVcf == null || germlineVcf.exists() == false) {
			Misc.printErrAndExit("Error: cannot find your high configence, single sample, germline variant call file "+germlineVcf);
		}
		
		if (germlineBamPileup == null || germlineBamPileup.exists() == false) {
			Misc.printErrAndExit("Error: cannot find your germline bampileup file "+germlineBamPileup);
		}
		
		if (somaticBamPileup == null || somaticBamPileup.exists() == false) {
			Misc.printErrAndExit("Error: cannot find your somatic bampileup file "+somaticBamPileup);
		}
		
		if (resultsDir == null ) {
			Misc.printErrAndExit("Error: cannot find your result directory");
		}
		resultsDir.mkdirs();
		
	}	


	/* Crazy number of preprocessing steps, need to do in a workflow
jj USeq_9.3.4.beta/Apps/VCF2Bed -v merged.vcf.gz
jj USeq_9.3.4.beta/Apps/MergeRegions -r mergedPad0bp.bed.gz -o var0bp.merged.bed
jj USeq_9.3.4.beta/Apps/MergeAdjacentRegions -b var0bp.merged.bed -r finalVars0bp.bed.gz -m 250
jj ~/USeqApps/BedTabix -v finalVars0bp.bed.gz -t ~/BioApps/HtsLib/1.15/bin/
module load samtools/1.16
samtools view -L finalVars0bp.bed.gz -@ 10 -T ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -o normal0bp.bam ../../1274232_21-0012878_NormalDNA_Hg38.cram
samtools view -L finalVars0bp.bed.gz -@ 10 -T ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -o tumor0bp.bam ../../1274232_21-0012878_TumorDNA_Hg38.cram
samtools index -b normal0bp.bam
samtools index -b tumor0bp.bam
jj ../USeq_9.3.4.beta/Apps/SamAlignmentDepthMatcher -m normal0bp.bam -s tumor0bp.bam -o tumor0bpSub.sam.gz -b finalVars0bp.bed.gz -t 10
	Or run SamReadDepthSubSampler on both tumor and normal setting it to 200?
samtools view -H tumor0bpSub.sam.gz > header.sam
samtools view tumor0bpSub.sam.gz | sort -u > body.sam
cat header.sam body.sam > tumorMatchedDept.sam
samtools sort -o tumor0bpSub.bam tumorMatchedDept.sam
rm *sam *sam.gz
samtools index -b tumor0bpSub.bam
jj ../USeq_9.3.4.beta/Apps/BamPileup -b tumor0bpSub.bam -r finalVars0bp.bed.gz -f ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -s tumor0bpSub.bp.txt.gz -t ~/BioApps/HtsLib/1.15/bin/ -q 13 -p 10
jj ../USeq_9.3.4.beta/Apps/BamPileup -b normal0bp.bam -r finalVars0bp.bed.gz -f ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -s normal0bp.bp.txt.gz -t ~/BioApps/HtsLib/1.15/bin/ -q 13 -p 10

LoH 
jj ~/USeqApps/AnnotateBedWithGenes -u /Users/u0028003/HCI/Annotations/GRCh38/hg38RefSeq9Dec2020_MergedStdChr.ucsc.gz -b ~/LoH/WithNoPadding/loh.bed -r ~/LoH/WithNoPadding/loh.anno.bed.gz -p 1000
*/


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                    LoH : Oct 2022                                **\n" +
				"**************************************************************************************\n" +
				"LoH compares heterozygous snv and indel allele counts between germline and somatic\n"+
				"sequencing datasets to identify potential loss of heterozygousity events. LoH is\n"+
				"defined here as a significant increase in the somatic allele fraction (AF) relative to\n"+
				"the matched germline. Adjacent variants passing p-value and AF difference thresholds\n"+
				"are merged and a combine p-value calculated for the window block. For best results,\n"+
				"use the USeq SamReadDepthSubSampler to limit the read depth over the variants in each\n"+
				"dataset to < 200.\n"+

				"\nRequired:\n"+
				"-v Path to a germline, single sample, variant xxx.vcf(.gz/.zip OK) file. These should\n"+
				"     be filtered, high confidence, vt normalized and decomposed short variants.\n" +
				"-g Path to the germline bampileup file, use the USeq BamPileup tool to generate the\n"+
				"     xxx.bp.txt.gz and tabix index files.\n"+
				"-s Path to the somatic bampileup file, ditto.\n"+
                "-o Path to directory to write the bed and vcf output results.\n"+

				"\nOptional:\n" +
				"-q Minimum vcf record QUAL, defaults to 20\n"+
                "-e Minimum vcf record sample GQ, defaults to 20\n"+
                "-x Use Strelka's recalibrated GQX genotype score, defaults to GQ\n"+
                "-r Minimum unique observation variant read depth for both germline and tumor samples,\n"+
                "     defaults to 20\n"+
                "-l Minimum Fisher exact p-value for including a variant in a block, defaults to 0.1\n"+
                "-f Minimum difference in AF between the somatic and germline samples for including\n"+
                "     a variant into a block, defaults to 0.1\n"+
                "-m Maximum bp gap between variants for grouping into a block, defaults to 1000\n"+
                "-p Minimum Benjamini-Hochberg adjusted window -10Log10(combine p-value), defaults to\n"+
                "     8.239 (0.15)\n"+
                "-d Minimum window mean difference in AFs, defaults to 0.1\n"+
                "-w Window BP padding for reporting significant LoH regions, defaults to 100\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/LoH -v gatkStrelkaGermline.vcf.gz \n"+
				"-g normal0bp.bp.txt.gz -s tumor0bpSub.bp.txt.gz -b loh.bed -d 0.15 -m 5000\n\n"+

				"**************************************************************************************\n");
	}

}
