package edu.utah.seq.cnv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.AnnotateBedWithGenes;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class GatkCalledSegmentAnnotator {

	//user defined fields
	private File segFile;
	private File tumorCopyRatioFile;
	private File normalCopyRatioFile;
	private File tumorAlleleFreqFile;
	private File normalAlleleFreqFile;
	private File resultsDirectory;
	private double maxAbsLg2NormalCopyRatio = 0.5;
	private double minCopyRatioMeanTNRatios = 0.15;
	private double minAbsLg2TumorCopyRatio = 0.15;
	private int maxGapGeneIntersection = 1000;
	private File ucscGeneTableFile;
	
	//internal fields
	private GatkSegment[] gatkSegments;
	
	public GatkCalledSegmentAnnotator(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);
		
		try {
			parseSegments();
			
			IO.pl("Adding copy ratios...");
			addTumorCopyRatios();
			addNormalCopyRatios();
			IO.pl("Adding allele frequencies...");
			addTumorAlleleFrequencies();
			addNormalAlleleFrequencies();
			
			checkCopyRatios();
			
			IO.pl("Adding intersecting genes...");
			intersectWithGenes();
			
			saveSegResults();
			saveSpreadsheetResults();
			saveFilteredBedResults();
			

			
		} catch (IOException e) {
			IO.el("\nERROR running GatkSegmentAnnotator\n");
			e.printStackTrace();
		}

		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}
	
	private void intersectWithGenes() {
		Bed[] bufferedRegions = new Bed[gatkSegments.length];
		for (int i=0; i< gatkSegments.length; i++) {
			int start = gatkSegments[i].getStart() - maxGapGeneIntersection;
			if (start < 0) start = 0;
			bufferedRegions[i] = new Bed(gatkSegments[i].getChr(), start, gatkSegments[i].getEnd()+maxGapGeneIntersection, "", 0, '.');
			gatkSegments[i].setBedRegion(bufferedRegions[i]);
		}
		new AnnotateBedWithGenes(ucscGeneTableFile, bufferedRegions, false);
	}

	private void saveSegResults() throws IOException {
		String[] names = Misc.UNDERSCORE.split(segFile.getName());
		String passName = names[0]+"_Pass";
		String failName = names[0]+"_Fail";
		
		File pass = new File(resultsDirectory, Misc.removeExtension(segFile.getName())+".anno.seg");
		
		PrintWriter out = new PrintWriter( new FileWriter( pass));
		
		String header = "Sample\tChr\tStart\tEnd\tCR_NumPoints\tCR_LgTumorMean\tCR_LgNormalMean\tCR_LgTNRatio\tAF_TumorSize\tAF_TumorMean\tAF_NormalSize\tAF_NormalMean\tAF_TNRatio\tGATK_Call\tGATK_LgTumorGeoMean";

		out.println(header);
		
		int numPass = 0;
		int numFail = 0;
		for (GatkSegment gs : gatkSegments) {
			if (Math.abs(gs.getLgMeanTNRatios()) >= minCopyRatioMeanTNRatios && Math.abs(gs.getLogMeanNormalCopyRatios()) <= maxAbsLg2NormalCopyRatio && Math.abs(gs.getLogMeanTumorCopyRatios()) >= minAbsLg2TumorCopyRatio) {
				out.println(gs.toSeg(passName));
				numPass++;
			}
			else {
				out.println(gs.toSeg(failName));
				numFail++;
			}
		}
		out.close();

		
		IO.pl("\nNumber Passing: "+numPass+"\tFailing: "+numFail);
	}
	
	private void saveSpreadsheetResults() throws IOException {
		File pass = new File(resultsDirectory, segFile.getName()+".xls");
		PrintWriter out = new PrintWriter( new FileWriter( pass));
		
		String header = "#Chr:Start-End\tPass Thresholds\tCR Call\tNum CR Points\tGeometric Mean Tumor CR\tLg2 Mean Tumor CR\tLg2 Mean Normal CR\tLg2 Mean TN CR Ratios\tNum Tumor AF Points\tMean Tumor AF\tNum Normal AF Points\tMean Normal AF\tMean TN AFs\tGenes";

		out.println(header);

		for (GatkSegment gs : gatkSegments) {
			boolean passThres;
			if (Math.abs(gs.getLgMeanTNRatios()) >= minCopyRatioMeanTNRatios && Math.abs(gs.getLogMeanNormalCopyRatios()) <= maxAbsLg2NormalCopyRatio && Math.abs(gs.getLogMeanTumorCopyRatios()) >= minAbsLg2TumorCopyRatio) passThres = true;
			else passThres = false;
			out.println(gs.toSpreadSheet(passThres));
		}
		out.close();
	}
	
	private void saveFilteredBedResults() throws IOException {
		File pass = new File(resultsDirectory, segFile.getName()+".pass.bed");
		PrintWriter out = new PrintWriter( new FileWriter( pass));
		
		out.println("#Chr Start Stop Info Lg2(meanTumCRs/meanNormCRs) Call");
		out.println("##Info numOb = number of copy ratio observations, typically exons");
		out.println("##Info lg2Tum = log2 of mean tumor copy ratio observations");
		out.println("##Info lg2Norm = log2 of mean normal copy ratio observations");
		out.println("##Info genes = affected genes");
		out.println("##Call + for amplification, - for loss");

		for (GatkSegment gs : gatkSegments) {
			if (Math.abs(gs.getLgMeanTNRatios()) >= minCopyRatioMeanTNRatios && Math.abs(gs.getLogMeanNormalCopyRatios()) <= maxAbsLg2NormalCopyRatio && Math.abs(gs.getLogMeanTumorCopyRatios()) >= minAbsLg2TumorCopyRatio) {
				out.println(gs.toBed());
			}
		}
		out.close();
	}

	private void checkCopyRatios() {
		for (GatkSegment gs : gatkSegments) {
			gs.calculateMeans();
			if (gs.copyRatiosOK() == false){
				IO.el("Problem with number of copy ratio points? These should be identical:");
				gs.printData();
				System.exit(1);
			}
		}
	}

	private void addNormalAlleleFrequencies() throws IOException {
		TabixReader tr = new TabixReader(normalAlleleFreqFile.toString());
		String line;
		TabixReader.Iterator it;
		String q;

		//for each Segment
		for (GatkSegment seg: gatkSegments){
			try {
				q = seg.getChr()+":"+seg.getStart()+"-"+seg.getEnd();
				it = tr.query(q);
				//any hits?
				if (it != null) {
					while ((line = it.next()) != null) {
						//CONTIG	POSITION	REF_COUNT	ALT_COUNT	REF_NUCLEOTIDE	ALT_NUCLEOTIDE
						//chrY	11344928	305	166	A	T
						String[] tokens = Misc.TAB.split(line);
						float refCount = Float.parseFloat(tokens[2]);
						float altCount = Float.parseFloat(tokens[3]);
						seg.getNormalAlleleFrequencies().add(altCount/(altCount+refCount));
					}
				}
			} catch (ArrayIndexOutOfBoundsException e){}
		}
		tr.close();	
	}
	
	private void addTumorAlleleFrequencies() throws IOException {
		TabixReader tr = new TabixReader(tumorAlleleFreqFile.toString());
		String line;
		TabixReader.Iterator it;
		String q;
		
		//for each Segment
		for (GatkSegment seg: gatkSegments){
			//watch out for no hits
			try {
				q = seg.getChr()+":"+seg.getStart()+"-"+seg.getEnd();
				it = tr.query(q);
				//any hits?
				if (it != null) {
					while ((line = it.next()) != null) {
						//CONTIG	POSITION	REF_COUNT	ALT_COUNT	REF_NUCLEOTIDE	ALT_NUCLEOTIDE
						//chrY	11344928	305	166	A	T
						String[] tokens = Misc.TAB.split(line);
						float refCount = Float.parseFloat(tokens[2]);
						float altCount = Float.parseFloat(tokens[3]);
						seg.getTumorAlleleFrequencies().add(altCount/(altCount+refCount));
					}
				}
			} catch (ArrayIndexOutOfBoundsException e){}
		}
		tr.close();	
	}
	
	
	private void addTumorCopyRatios() throws IOException {
		TabixReader tr = new TabixReader(tumorCopyRatioFile.toString());
		String line;
		TabixReader.Iterator it;
		String q;

		//for each Segment
		for (GatkSegment seg: gatkSegments){
			try {
				q = seg.getChr()+":"+seg.getStart()+"-"+seg.getEnd();
				it = tr.query(q);
				//any hits?
				if (it != null) {
					while ((line = it.next()) != null) {
						//chrX	156014883	156015954	0.444810
						String[] tokens = Misc.TAB.split(line);
						seg.getTumorCopyRatios().add(Float.parseFloat(tokens[3]));
					}
				}
			} catch (ArrayIndexOutOfBoundsException e){}
		}
		tr.close();	
	}
	
	private void addNormalCopyRatios() throws IOException {
		TabixReader tr = new TabixReader(normalCopyRatioFile.toString());
		String line;
		TabixReader.Iterator it;
		String q;

		//for each Segment
		for (GatkSegment seg: gatkSegments){
			try {
				q = seg.getChr()+":"+seg.getStart()+"-"+seg.getEnd();
				it = tr.query(q);
				//any hits?
				if (it != null) {
					while ((line = it.next()) != null) {
						//chrX	156014883	156015954	0.444810
						String[] tokens = Misc.TAB.split(line);
						seg.getNormalCopyRatios().add(Float.parseFloat(tokens[3]));
					}
				}
			} catch (ArrayIndexOutOfBoundsException e){}
		}
		tr.close();	
	}


	private void parseSegments() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(segFile);
		String line;
		
		//check header
		boolean contigFound = false;
		while ((line = in.readLine()) != null){
			if (line.length() == 0 || line.startsWith("@")) continue;
			if (line.equals("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\tCALL")){
				contigFound = true;
				break;
			}
		}
		
		if (contigFound == false) throw new IOException("\nFailed to find the 'CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\tCALL' line. Is this a CallCopyRatioSegments xxx.called.seg file?");

		ArrayList<GatkSegment> segAL = new ArrayList<GatkSegment>();
		while ((line = in.readLine()) != null) segAL.add(new GatkSegment(line));
		gatkSegments = new GatkSegment[segAL.size()];
		segAL.toArray(gatkSegments);
		Arrays.sort(gatkSegments);
		
		in.close();
	}
	


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new GatkCalledSegmentAnnotator(args);
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
					case 's': segFile = new File(args[++i]); break;
					case 't': tumorCopyRatioFile = new File(args[++i]); break;
					case 'n': normalCopyRatioFile = new File(args[++i]); break;
					case 'u': tumorAlleleFreqFile = new File(args[++i]); break;
					case 'o': normalAlleleFreqFile = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'g': ucscGeneTableFile = new File(args[++i]); break;
					case 'a': maxGapGeneIntersection = Integer.parseInt(args[++i]); break;
					case 'm': minCopyRatioMeanTNRatios = Double.parseDouble(args[++i]); break;
					case 'c': minAbsLg2TumorCopyRatio = Double.parseDouble(args[++i]); break;
					case 'x': maxAbsLg2NormalCopyRatio = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		
		//check results dir
		if (resultsDirectory == null) Misc.printErrAndExit("\nError: provide a directory to save results.\n");
		resultsDirectory.mkdirs();
		
		//check files
		File[] files = new File[]{segFile, tumorCopyRatioFile, normalCopyRatioFile, tumorAlleleFreqFile, normalAlleleFreqFile, ucscGeneTableFile};
		boolean pass = true;
		for (File f: files) {
			if (f== null || f.exists() == false) {
				pass = false;
				break;
			}
		}
		if (pass == false){
			IO.pl("\n\nFAILED to find one of the following files:");
			printOptions();
			System.exit(1);
		}

		printOptions();
		
	}	
	


	private void printOptions() {
		IO.pl("Settings:");
		IO.pl("\t-s SegFile\t"+segFile);
		IO.pl("\t-t TumorCopyRatioFile\t"+tumorCopyRatioFile);
		IO.pl("\t-n NormalCopyRatioFile\t"+normalCopyRatioFile);
		IO.pl("\t-u TumorAlleleFreqFile\t"+tumorAlleleFreqFile);
		IO.pl("\t-o NormalAlleleFreqFile\t"+normalAlleleFreqFile);
		IO.pl("\t-r ResultsDirectory\t"+resultsDirectory);
		IO.pl("\t-g RefFlatGeneFile\t"+ucscGeneTableFile);
		IO.pl("\t-a MaxGapForGeneIntersection\t"+maxGapGeneIntersection);
		IO.pl("\t-m MinimumAbsLg2TNRatioOfCopyRatios\t"+minCopyRatioMeanTNRatios);
		IO.pl("\t-c MinimumAbsLg2TumorCopyRatio\t"+minAbsLg2TumorCopyRatio);
		IO.pl("\t-x MaximumAbsLg2NormalCopyRatio\t"+maxAbsLg2NormalCopyRatio);
		IO.pl();
	}

	public static void printDocs(){
		
		
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Gatk Called Segment Annotator: August 2019                  **\n" +
				"**************************************************************************************\n" +
				"Annotates GATKs CallCopyRatioSegments output with denoised copy ratio and heterozygous\n"+
				"allele frequency data from the tumor and matched normal samples. Enables filtering\n"+
				"using these values to remove copy ratio calls with high normal background. Adds \n"+
				"intersecting gene names.\n"+

				"\nRequired Options:\n"+
				"-r Results directory to save the passing and failing segments.\n"+
				"-s Called segment file from GATKs CallCopyRatioSegments app, e.g. xxx.called.seg\n"+
				"-t Tumor denoised copy ratio file, from GATKs DenoiseReadCounts app. Bgzip compress\n"+
				"      and tabix index it with https://github.com/samtools/htslib :\n"+
				"      grep -vE '(@|CONTIG)' tumor.cr.tsv > tumor.cr.txt\n"+
				"      ~/HTSLib/bgzip tumor.cr.txt\n"+
				"      ~/HTSLib/tabix -s 1 -b 2 -e 3 tumor.cr.txt.gz\n"+
				"-n Normal denoised copy ratio file, ditto.\n"+
				"-u Tumor allele frequency file, from GATKs ModelSegments app. Bgzip compress\n"+
				"      and tabix index it with https://github.com/samtools/htslib :\n"+
				"      grep -vE '(@|CONTIG)' gbm7.hets.tsv > gbm7.hets.txt\n"+
				"      ~/HTSLib/bgzip gbm7.hets.txt\n"+
				"      ~/HTSLib/tabix -s 1 -b 2 -e 2 gbm7.hets.txt.gz\n"+
				"-o Normal allele frequency file, ditto.\n"+
				"-g RefFlat UCSC gene file, run USeq's MergeUCSCGeneTable to collapse transcripts.\n"+
				
				"\nDefault Options:\n"+
				"-c Minimum absolute tumor log2 copy ratio, defaults to 0.15\n"+
				"-x Maximum absolute normal log2 copy ratio, defaults to 0.5\n"+
				"-m Minimum absolute log2 TN ratio of copy ratios, defaults to 0.15\n"+
				"-a Maximum bp gap for intersecting a segment with a gene, defaults to 1000\n"+

				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/GatkCalledSegmentAnnotator -r AnnoResults/\n"+
				"       -s gbm7.called.seg -t tumor.cr.txt.gz -n normal.cr.txt.gz -u gbm7.hets.txt.gz\n"+
				"       -o gbm7.hets.normal.txt.gz -g ~/UCSC/hg38RefSeq_Merged.refFlat.gz -a 100 \n\n" +

				"**************************************************************************************\n");
	}
}
