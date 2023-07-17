package edu.utah.seq.vcf.loh;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.parsers.jpileup.BamPileupTabixLoaderSingle;
import edu.utah.seq.parsers.jpileup.BpileupLine;
import util.gen.FisherExact;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**
 * @author Nix*/
public class NonLoH {

	//user defined fields
	private File normalVcf = null;
	private File normalBamPileup = null;
	private File somaticBamPileup = null;
	private File[] ponBamPileups = null;
	private double minimumQUAL = 20;
	private double minimumGQ = 20;
	private int minimumDP = 50;   //keep this high for better StdDev
	private File resultsDir = null;
	private int bpPadding = 4;
	private double maxBkgAf = 0.025;
	private boolean excludeXY = true;
	
	//het var
	private double minVarPval = 0.9;
	private String gqSelector = "GQ";

	//internal
	private HashMap<String, String> formatValues = new HashMap<String,String>();
	private HeterozygousVar[] hetVcfRecords = null;
	private LinkedHashMap<String,HeterozygousVar[]> chromHets = null;
	private PassingHet[] passingHets = null;
	private int maxSizeForFishers = 0;
	private double[] tumorStats;
	private double[] normalStats;
	
	//mean, median, pseudoMedian  scalars, pse is best
	//  0     1           2  
	private int scalerToUseIndex = 2;

	//constructor
	public NonLoH(String[] args){
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);

			IO.pl("Thresholds:\n"+fetchThresholds(""));

			//fetch heterozygous snv normal variants
			fetchHeterozygousVars();

			//bam pileup counts from the somatic and normal datasets, could thread this!
			addBamPileupObservations();

			printAfStats();

			calculateFisherPValues();

			splitHetsByChrom();
			
			filterHets();
			
			normalizeTNPassingHets();
			
			normalizeBPs();
			

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		} finally {

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" min\n");
		}
	}
	
	private void normalizeBPs() throws Exception {
		IO.pl("\nCalculating normalization stats for each bam pileup file...");
		IO.pl("Name\t#Pass\t#Fail\tMean\tMedian\tPseMedian\tStdDev\tMin\tMax\t10th\t90th");
		
		//do the tumor and normal
		normalizePassingHets("Tumor", somaticBamPileup, passingHets, minimumDP/2, tumorStats);
		normalizePassingHets("Normal", normalBamPileup, passingHets, minimumDP/2, normalStats);
		
		//then do all of the normals in the PoN
		float[][] scaledCounts = new float[ponBamPileups.length][];
		for (int i=0; i< ponBamPileups.length; i++) {
			scaledCounts[i] = normalizePassingHets("PoN", ponBamPileups[i], passingHets, minimumDP/2, null);
		}
		
		//save them to file
		Gzipper out = new Gzipper (new File (this.resultsDir, "scaledCountForPoN.txt.gz"));
		for (int i=0; i< ponBamPileups.length; i++) {
			out.println("#\t"+i+"\t"+ponBamPileups[i].getCanonicalPath());
		}
		out.println("#Coordinates\t1KScaledCounts\tZScores");
		for (int i=0; i< passingHets.length; i++) {
			out.print( passingHets[i].getChr()+":"+passingHets[i].getPosition());
			double[] scs = new double[ponBamPileups.length];
			
			//output the actual counts
			for (int j=0; j< ponBamPileups.length; j++) {
				out.print("\t");
				out.print(scaledCounts[j][i]);
				scs[j] = scaledCounts[j][i];
			}
			//output the z-scores
			double[] zScores = Num.convertToZScores(scs);
			for (int j=0; j< ponBamPileups.length; j++) {
				out.print("\t");
				out.print((float)zScores[j]);
			}
			out.println();
		}
		out.close();
	}

	private float[] normalizePassingHets(String prefixName, File bamPileup, PassingHet[] hets, float minCount, double[] stats) throws Exception {
		int numSkipped = 0;
		int numPass = 0;
		float[] scaledCounts = null;
		
		//pull counts and normalize?
		if (stats == null) {
			BamPileupTabixLoaderSingle bps = new BamPileupTabixLoaderSingle(bamPileup, 0);
			ArrayList<Float> scalars = new ArrayList<Float>();
			ArrayList<Float> counts = new ArrayList<Float>();
			for (int i=0; i< hets.length; i++) {
				ArrayList<BpileupLine> lines = bps.fetchLines(hets[i].getTabixCoor());
				if (lines.size() == 1) {
					//Misc.printErrAndExit(ph.getTabixCoor()+" new lu "+ lines.get(0).getLine());
					Float count = new Float(lines.get(0).getSamples()[0].getPassingReadCoverageSnv());
					counts.add(count);
					if (count >= minCount) scalars.add(1000.0f/count);
					else numSkipped++;
				}
				else numSkipped++;
			}
			bps.getTabixReader().close();
			float[] sc = Num.arrayListOfFloatToArray(scalars);
			Arrays.sort(sc);
			//mean, median, pseudoMedian, stdDev, min, max, 10th, 90th
			stats = statFloatArray(sc);
			numPass = scalars.size();
			
			//scale the counts
			scaledCounts = new float[counts.size()];
			for (int i=0; i< scaledCounts.length; i++) {
				double scaled = stats[scalerToUseIndex] * (double)counts.get(i);
				scaledCounts[i] = (float) scaled;
			}

			
		}
		else numPass = hets.length;
		
		printNormalizedStats(prefixName+"_"+bamPileup.getName(), numPass, numSkipped, stats);
		
		//File normalized = new File(normalizedReadCovDir, "normalized_"+Misc.removeExtension(bamPileup.getName())+".txt.gz");
		//normalizeBamPileup(bamPileup, normalized, stats[2]);
		
		return scaledCounts;

	}
	
	private static void printNormalizedStats(String name, int numPass, int numSkipped, double[] stats) {
		//mean, median, pseudoMedian, stdDev, min, max, 10th, 90th
		StringBuilder sb = new StringBuilder(name);
		sb.append("\t");
		sb.append(numPass);
		sb.append("\t");
		sb.append(numSkipped);
		for (double d: stats) {
			sb.append("\t");
			sb.append(d);
		}
		IO.pl(sb);
	}

	/**Provide sorted float[]
	 * Returns mean, median, pseudoMedian, stdDev, min, max, 10th, 90th*/
	public static double[] statFloatArray(float[] sortedFloat){
		//calc mean
		double mean = Num.mean(sortedFloat); 
		//calc median, need array copy to not sort referenced float!
		double median = Num.median(sortedFloat);
		//pseudomedian
		double pseudoMedian = Num.pseudoMedian(sortedFloat);
		//standard deviation
		double stdDev = Num.standardDeviation(sortedFloat,mean);
		//calc min, max
		float[] minMax = Num.findMinMaxFloatValues(sortedFloat);
		//calc 10th
		double perc10th = Num.percentile(sortedFloat, 0.1);
		//calc 90th
		double perc90th = Num.percentile(sortedFloat,0.9);
		return new double[] {mean, median, pseudoMedian, stdDev, minMax[0], minMax[1], perc10th, perc90th};
	}

	private void normalizeBamPileup(File bpFile, File normalizedGzipFile, double scalar) throws Exception {
		IO.pl("\nNormalizing "+bpFile.getName()+" with "+scalar);
		Gzipper gzOut = new Gzipper(normalizedGzipFile);
		BufferedReader in = IO.fetchBufferedReader(bpFile);
		String line;
		ArrayList<Float> normalizedCounts = new ArrayList<Float>();
		while ((line = in.readLine())!=null) {
			if (line.startsWith("#")) gzOut.println(line);
			else {
				BpileupLine bpl = new BpileupLine(line);
				double normalizedRC = bpl.getSamples()[0].getPassingReadCoverage() * scalar;
				gzOut.println(bpl.getChr()+"\t"+bpl.getZeroPos()+"\t"+normalizedRC);
				normalizedCounts.add(new Float(normalizedRC));
			}
		}
		
		IO.pl("\tMean\tMedian\tStdDev\tMin\tMax\t10th\t90th");
		float[] sc = Num.arrayListOfFloatToArray(normalizedCounts);
		Arrays.sort(sc);
		IO.pl("\t"+ Num.statFloatArray(sc));
		in.close();
		gzOut.close();
	}

	private String fetchThresholds(String prePend) {
		StringBuilder sb = new StringBuilder();
		//vcf
		sb.append(prePend); sb.append("   minimumVcfQUAL"); sb.append("\t"); sb.append(minimumQUAL); sb.append("\n");
		sb.append(prePend); sb.append("   minimumVcfGQ"); sb.append("\t"); sb.append(minimumGQ); sb.append("\n");
		sb.append(prePend); sb.append("   minimumVcfDP"); sb.append("\t"); sb.append(minimumDP); sb.append("\n");
		//het var
		sb.append(prePend); sb.append("   minVarPval"); sb.append("\t"); sb.append(minVarPval); sb.append("\n");
		return sb.toString();
	}

	private void filterHets() throws IOException {
		int numFail = 0;
		ArrayList<PassingHet> ph = new ArrayList<PassingHet>();
		for (String chr: chromHets.keySet()) {
			if (chr.contains("X") || chr.contains("Y")) {
				if (excludeXY) continue;
			}
			for (HeterozygousVar h: chromHets.get(chr)) {
				if (h.getPvalue() < minVarPval) numFail++;
				else ph.add(new PassingHet(chr, h.getPosition(), h.getSomaticDp(), h.getGermlineDp()));
			}
		}
		passingHets = new PassingHet[ph.size()];
		ph.toArray(passingHets);
		IO.pl(ph.size()+" Passing PVal HetVars, "+numFail+" Failing HetVars");	
	}
	
	private void normalizeTNPassingHets() throws IOException {
		IO.pl("Generating normalization stats for the TN and saving a log2(normT/N) bed file...");
		//collect tumor and normal 1K scalars
		float[] tumorScalars = new float[passingHets.length];
		float[] normalScalars = new float[tumorScalars.length];
		for (int i=0; i< normalScalars.length; i++) {
			PassingHet ph = passingHets[i];
			tumorScalars[i] = (float)ph.fetchTumor1KScalar();
			normalScalars[i] = (float)ph.fetchNormal1KScalar();
		}
		
		//sort them
		Arrays.sort(tumorScalars);
		Arrays.sort(normalScalars);
		
		//mean, median, pseudoMedian, stdDev, min, max, 10th, 90th
		//  0     1           2          3     4    5     6    7
		tumorStats = statFloatArray(tumorScalars);
		normalStats = statFloatArray(normalScalars);
		
		//Output file for histo and boxplots the log2 T/N rto
		File bedFile = new File(resultsDir, "passingHetsLg2TNRto.bed.gz");
		Gzipper out = new Gzipper(bedFile);
		double ratioSum = 0;
		for (int i=0; i< normalScalars.length; i++) {
			PassingHet ph = passingHets[i];
			out.println(ph.fetchBedLine(tumorStats[scalerToUseIndex], normalStats[scalerToUseIndex]));
			double rto = ph.fetchNormalizedTNRatio(tumorStats[scalerToUseIndex], normalStats[scalerToUseIndex]);
			ratioSum+= rto;
		}
		out.close();
		IO.pl("\tLog2(aveNormalize(T/NRatios))\t"+Num.log2(ratioSum/(double)normalScalars.length));
	}
	
	private void calculateFisherPValues() throws IOException {
		IO.pl("\nCalculating Fisher's Exact PValues...");
		FisherExact fe = new FisherExact(maxSizeForFishers);
		for (HeterozygousVar s: hetVcfRecords) {
			int[] c = s.fetchSomAltRefGermAltRefCounts();
			//fe.getRightTailedP(somAlt, somRef, germAlt, germRef), only interested in a higher somatic AF relative to normal
			//double pval = fe.getRightTailedP(c[0], c[1], c[2], c[3]);
			//want pval for either up or down
			double pval = fe.getTwoTailedP(c[0], c[1], c[2], c[3]);
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
				if (chromHets.containsKey(currChrom)) throw new IOException("\nIs your normal vcf sorted? Seeing this chrom again: "+currChrom);
			}
			al.add(hetVcfRecords[i]);
		}
		//add last to hash
		HeterozygousVar[] sub = new HeterozygousVar[al.size()];
		al.toArray(sub);
		chromHets.put(currChrom, sub);
	}
	
	private void printAfStats() throws IOException {
		float[] germAf = new float[hetVcfRecords.length];
		float[] somAf = new float[hetVcfRecords.length];
		float[] germDp = new float[hetVcfRecords.length];
		float[] somDp = new float[hetVcfRecords.length];
		
		for (int i=0; i< hetVcfRecords.length; i++) {
			germAf[i] = (float)hetVcfRecords[i].getAlleleFractionGermline();
			somAf[i] = (float)hetVcfRecords[i].getAlleleFractionSomatic();
			int[] sAltRefgAltRef = hetVcfRecords[i].fetchSomAltRefGermAltRefCounts();
			somDp[i] = (float)(sAltRefgAltRef[0]+ sAltRefgAltRef[1]);
			germDp[i] = (float)(sAltRefgAltRef[2]+ sAltRefgAltRef[3]);
		}
		//sort em
		Arrays.sort(germAf);
		Arrays.sort(somAf);
		Arrays.sort(germDp);
		Arrays.sort(somDp);
		
		IO.pl("\nMean\tMedian\tStdDev\tMin\tMax\t10th\t90th");
		IO.pl("NormalAF\t"+ Num.statFloatArray(germAf));
		IO.pl("SomaticAF\t"+Num.statFloatArray(somAf));
		IO.pl("NormalDP\t"+Num.statFloatArray(germDp));
		IO.pl("SomaticDP\t"+ Num.statFloatArray(somDp));
	}
	

	private void addBamPileupObservations() throws Exception {
		IO.pl("\nChecking variant counts and background...");

		BamPileupTabixLoaderSingle bpg = new BamPileupTabixLoaderSingle(normalBamPileup, bpPadding);
		BamPileupTabixLoaderSingle bps = new BamPileupTabixLoaderSingle(somaticBamPileup, bpPadding);
		ArrayList<HeterozygousVar> passing = new ArrayList<HeterozygousVar>();
		int sizeWindow = bpPadding + bpPadding+ 1;
		int failingWindowSize = 0;
		int failingDp = 0;
		int failingBkg = 0;
		
		for (HeterozygousVar hs: hetVcfRecords) {
			char alt = hs.getVcfRecord()[4].charAt(0);

			ArrayList<BpileupLine> germ = bpg.fetchBpileupRecords(hs.getVcfRecord());
			if (germ.size() != sizeWindow) {
				failingWindowSize++;
				continue;
			}

			//check the read depth
			hs.setGermlineBPs(germ.get(bpPadding), minimumDP);
			if (hs.isPassing() == false) {
				failingDp++;
				continue;
			}
			
			//check
			boolean passBkd = checkBkg(alt,germ);
			/*if (passBkd == false) {
				IO.pl(Misc.stringArrayToString(hs.getVcfRecord(), " "));
				for (BpileupLine l: germ){
					IO.pl(l.getLine());
				}
				System.exit(0);
			}*/
			if (passBkd == false) {
				failingBkg++;
				continue;
			}
			
			ArrayList<BpileupLine> som = bps.fetchBpileupRecords(hs.getVcfRecord());
			hs.setSomaticBPs(som.get(bpPadding), minimumDP);
			if (hs.isPassing()== false) {
				failingDp++;
				continue;
			}
			
			passBkd = checkBkg(alt,som);
			if (passBkd == false) {
				failingBkg++;
				continue;
			}

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
		IO.pl("   "+failingDp+"\tFailing DP");
		IO.pl("   "+failingBkg+"\tFailing background");
		IO.pl("   "+failingWindowSize+"\tFailing bkg window size");
		IO.pl("   "+passing.size()+"\tPassing het snvs");
	}
	
	private boolean checkBkg(char alt, ArrayList<BpileupLine> germ) {
		//check before
		for (int i=0; i< bpPadding; i++) {
			double bkgAf = germ.get(i).getSamples()[0].findNonRefAltSnvBackground('x');
			if (bkgAf > maxBkgAf) return false;

		}
		//check after
		for (int i=bpPadding+1; i< germ.size(); i++) {
			double bkgAf = germ.get(i).getSamples()[0].findNonRefAltSnvBackground('x');
			if (bkgAf > maxBkgAf) return false;
		}
		//check the actual het snv
		double bkgAf = germ.get(bpPadding).getSamples()[0].findNonRefAltSnvBackground(alt);
		if (bkgAf > maxBkgAf) return false;
		return true;
	}

	private void fetchHeterozygousVars() throws Exception {
		ArrayList<HeterozygousVar> hetVcfRec = new ArrayList<HeterozygousVar>();
		int numRecords = 0;
		int numFailingRecords = 0;
		int numPassingHet = 0;
		String line = null;
		BufferedReader in = IO.fetchBufferedReader(normalVcf);
			
		while ((line = in.readLine()) != null){
			line = line.trim();
			if (line.length()==0) continue;
			if (line.startsWith("#CHROM")) break;
		}
		

		//for each data line in the file
		while ((line = in.readLine()) != null){
			numRecords++;

			//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1, Sample2.....
			//   0    1   2  3   4   5     6      7     8      9       10
			String[] fields = Misc.TAB.split(line);

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
			
			//heterozygous
			if (hetHom == 0 && fields[3].length()==1 && fields[4].length()==1) {
					numPassingHet++;
					hetVcfRec.add(new HeterozygousVar(fields));
			}
		}
		in.close();
		hetVcfRecords = new HeterozygousVar[hetVcfRec.size()];
		hetVcfRec.toArray(hetVcfRecords);

		IO.pl("Normal VCF Parsing Statistics:");
		IO.pl("   "+ numRecords + "\t# Records");
		IO.pl("   "+ numFailingRecords + "\t# Failing Records");
		IO.pl("   "+ numPassingHet + "\t# Passing Het Snvs");
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
		new NonLoH(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File ponDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': normalVcf = new File(args[++i]); break;
					case 'g': normalBamPileup = new File(args[++i]); break;
					case 's': somaticBamPileup = new File(args[++i]); break;
					case 'o': resultsDir = new File(args[++i]); break;
					case 'p': ponDir = new File(args[++i]); break;
					case 'r': minimumDP = Integer.parseInt(args[++i]); break;
					case 'q': minimumQUAL = Double.parseDouble(args[++i]); break;
					case 'e': minimumGQ = Double.parseDouble(args[++i]); break;
					case 'x': gqSelector = "GQX"; break;
					case 'l': minVarPval = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (normalVcf == null || normalVcf.exists() == false) {
			Misc.printErrAndExit("Error: cannot find your high configence, single sample, normal variant call file "+normalVcf);
		}
		
		if (normalBamPileup == null || normalBamPileup.exists() == false) {
			Misc.printErrAndExit("Error: cannot find your normal bampileup file "+normalBamPileup);
		}
		
		if (somaticBamPileup == null || somaticBamPileup.exists() == false) {
			Misc.printErrAndExit("Error: cannot find your somatic bampileup file "+somaticBamPileup);
		}
		
		if (resultsDir == null ) {
			Misc.printErrAndExit("Error: cannot find your result directory");
		}
		resultsDir.mkdirs();
		
		if (ponDir == null || ponDir.isDirectory()==false ) {
			Misc.printErrAndExit("Error: cannot find your PoN directory containing normal sample bam pileup files? "+ponDir);
		}
		ponBamPileups = IO.extractFiles(ponDir, ".gz");
		if (ponBamPileups == null || ponBamPileups.length ==0) {
			Misc.printErrAndExit("Error: failed to find any xxx.gz bam pileup files in "+ponDir);
		}
	}	


	/* Crazy number of preprocessing steps, need to do in a workflow
java -jar -Xmx100G ~/USeqApps/VCFConsensus -p *_JointGenotyped.vcf.gz -s *_NormalDNA_Hg38_JointGenotyped.vcf.gz -o gatkIllum.vcf.gz -u -c
java -jar -Xmx100G ~/USeqApps/VCF2Bed -v gatkIllum.vcf.gz
java -jar -Xmx100G ~/USeqApps/MergeRegions -r gatkIllumPad0bp.bed.gz -o gatkIllum.merged.bed
java -jar -Xmx100G ~/USeqApps/MergeAdjacentRegions -b gatkIllum.merged.bed -r finalRegions.bed.gz -m 500
java -jar -Xmx100G ~/USeqApps/BedTabix -v finalRegions.bed.gz -t ~/BioApps/HtsLib/1.15/bin/
module load samtools/1.16
samtools view -L finalRegions.bed.gz -@ 10 -T ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -o normal.bam *_NormalDNA_Hg38.cram &
samtools view -L finalRegions.bed.gz -@ 10 -T ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -o tumor.bam *_TumorDNA_Hg38.cram && echo DONE &
samtools index -b normal.bam &
samtools index -b tumor.bam &
java -jar -Xmx100G ~/USeqApps/SamReadDepthSubSampler -a normal.bam -b finalRegions.bed.gz -t 200 -p -x 1000 -q 13 &
java -jar -Xmx100G ~/USeqApps/SamReadDepthSubSampler -a tumor.bam -b finalRegions.bed.gz -t 200 -p -x 1000 -q 13 &

##jj ../USeq_9.3.4.beta/Apps/SamAlignmentDepthMatcher -m normal0bp.bam -s tumor0bp.bam -o tumor0bpSub.sam.gz -b finalVars0bp.bed.gz -t 10
##samtools view -H tumor0bpSub.sam.gz > header.sam
##samtools view tumor0bpSub.sam.gz | sort -u > body.sam
##cat header.sam body.sam > tumorMatchedDept.sam
##samtools sort -o tumor0bpSub.bam tumorMatchedDept.sam
##rm *sam *sam.gz
##samtools index -b tumor0bpSub.bam

java -jar -Xmx100G ~/USeqApps/BamPileup -b tumorRDSub200.bam -r finalRegions.bed.gz -f ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -s tumor.sub.bp.txt.gz -t ~/BioApps/HtsLib/1.15/bin/ -q 13 -p 10
java -jar -Xmx100G ~/USeqApps/BamPileup -b normalRDSub200.bam -r finalRegions.bed.gz -f ~/TNRunner/Indexes/B38IndexForBwa-0.7.17/hs38DH.fa -s normal.sub.bp.txt.gz -t ~/BioApps/HtsLib/1.15/bin/ -q 13 -p 10

java -jar -Xmx100G  ~/USeqApps/LoH -v gatkIllum.vcf.gz -g normal.sub.bp.txt.gz -s tumor.sub.bp.txt.gz -o LoH

java -jar -Xmx100G  ~/USeqApps/AnnotateBedWithGenes -u ~/TNRunner/AnnotatorData/UCSC/9Dec2020/hg38RefSeq9Dec2020_MergedStdChr.ucsc.gz -b ~/LoH/loh.bed.gz -r ~/LoH/WithNoPadding/loh.anno.bed.gz -p 1000
*/


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                   Non LoH : Nov 2022                             **\n" +
				"**************************************************************************************\n" +
				"LoH compares heterozygous snv and indel allele counts between normal and somatic\n"+
				"sequencing datasets to identify potential loss of heterozygosity events. LoH is\n"+
				"defined here as a significant increase in the somatic allele fraction (AF) relative to\n"+
				"the matched normal. Adjacent variants passing p-value and AF difference thresholds\n"+
				"are merged and a combine p-value calculated for the window block.\n"+

				"\nRequired:\n"+
				"-v Path to a normal, single sample, variant xxx.vcf(.gz/.zip OK) file. These should\n"+
				"     be filtered, high confidence, vt normalized and decomposed short variants.\n" +
				"-g Path to the normal bampileup file, use the USeq BamPileup tool to generate the\n"+
				"     xxx.bp.txt.gz and tabix index files.\n"+
				"-s Path to the somatic bampileup file, ditto.\n"+
				"-p Path to a directory containing normal sample bampileup files for the PoN\n"+
                "-o Path to directory to write the bed and vcf output result files.\n"+

				"\nOptional:\n" +
				"-c Bed file containing passing regions from the TNRunner2 CopyRatio workflow.\n"+
				"-q Minimum vcf record QUAL, defaults to 20\n"+
                "-e Minimum vcf record sample GQ, defaults to 20\n"+
                "-x Use Strelka's recalibrated GQX genotype score, defaults to GQ\n"+
                "-r Minimum unique observation variant read depth for both normal and tumor samples,\n"+
                "     defaults to 20\n"+
                "-l Minimum Fisher exact p-value for including a variant in a block, defaults to 0.1\n"+
                "-f Minimum difference in AF between the somatic and normal samples for including\n"+
                "     a variant into a block, defaults to 0.1\n"+
                "-m Maximum bp gap between variants for grouping into a block, defaults to 1000\n"+
                "-p Minimum Benjamini-Hochberg adjusted window -10Log10(combine p-value), defaults to\n"+
                "     8.239 (0.15)\n"+
                "-d Minimum window mean difference in AFs, defaults to 0.1\n"+
                "-w Window BP padding for reporting significant LoH regions, defaults to 100\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/LoH -v gatkStrelkaNormal.vcf.gz \n"+
				"-g normal.bp.txt.gz -s tumor.bp.txt.gz -o LoHResults/ -d 0.15 -m 5000 -b\n"+
				"GATKCopyRatio_Hg38.called.seg.pass.bed.gz\n\n"+

				"**************************************************************************************\n");
	}

	public int getMinimumDP() {
		return minimumDP;
	}

	public PassingHet[] getPassingHets() {
		return passingHets;
	}

}
