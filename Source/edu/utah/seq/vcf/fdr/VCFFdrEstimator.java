package edu.utah.seq.vcf.fdr;

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

/**Estimates Fdr for each somatic variant based on mock tumor vs normal contrasts. 
 * @author Nix
 * */
public class VCFFdrEstimator {

	//user fields
	private File[] mockVcfs;
	private File[] realVcfs;
	private File saveDir;
	private double minFdr = 6;
	
	private float[][] mockQuals = null;
	private float[][] realQuals = null;
	private float[] mockThresholds = null;
	private float[] mockCounts = null;

	private Histogram mockQualScoreHistogram = null;
	

	//constructor
	public VCFFdrEstimator(String[] args) throws FileNotFoundException, IOException {
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//analyze mock vcfs
		analyzeMocks();

		parseReal();
		
		//print stats
		printStats();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void parseReal() {
		File f = null;
		String line = null;
		try {
			IO.pl("\nScoring VCF files (Name #Pass #Fail), minFDR "+minFdr);
			
			//for each file
			for (File vcfFile: realVcfs) {
				f= vcfFile;
				IO.p("\t"+vcfFile.getName());
				File passFile = new File (saveDir, Misc.removeExtension(vcfFile.getName())+".pass"+minFdr+"Fdr.vcf.gz");
				Gzipper outPass = new Gzipper (passFile);
				File failFile = new File (saveDir, Misc.removeExtension(vcfFile.getName())+".fail"+minFdr+"Fdr.vcf.gz");
				Gzipper outFail = new Gzipper (failFile);
				VCFFdrDataset vcfDataset = new VCFFdrDataset(vcfFile, minFdr);
				vcfDataset.addHeader(outPass, outFail);
				vcfDataset.estimateFdrs(outPass, outFail, this);
				outPass.close();
				outFail.close();
				if (vcfDataset.getNumFail() == 0) failFile.deleteOnExit();
				if (vcfDataset.getNumPass() == 0) passFile.deleteOnExit();
				
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem processing "+f+" record "+line);
		}

		/*
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
		 */
	}

	private void printStats() {
		/*IO.pl("QUAL score stats:");
		
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
		*/
	}

	
	

	

	private void analyzeMocks() {
		//collect qual scores
		mockQuals = new float[this.mockVcfs.length][];
		float min = 1000;
		float max = 0;
		HashSet<Float> allScores = new HashSet<Float>();
		allScores.add(0.0f);
		for (int i=0; i< mockVcfs.length; i++) {
			mockQuals[i] = fetchSortedQualScores(mockVcfs[i]);
			for (float f: mockQuals[i]) {
				allScores.add(f);
				if (f<min) min = f;
				if (f>max) max = f;
			}
		}
		
		//load and print histogram
		printMockHistogram(min, max);
		
		//load real vcfs
		realQuals = new float[realVcfs.length][];
		for (int i=0; i< realVcfs.length; i++) {
			realQuals[i] = fetchSortedQualScores(realVcfs[i]);
			for (float f: realQuals[i]) allScores.add(f);
		}
		
		
		//Calculate # FPs vs Score, at every threshold, will start at zero
		calculateCountVsScoreTable(allScores);
	}
	
	



	private void calculateCountVsScoreTable(HashSet<Float> allScores) {
		File table = new File (this.saveDir, "qualVcfCountTable.txt");
		IO.pl("\nQUAL threshold VCF record count table to "+table);
		mockThresholds = new float[allScores.size()];
		mockCounts = new float[allScores.size()];
		try {
			float[] scores = Num.hashSetToFloat(allScores);
			Arrays.sort(scores);
			//print header
			PrintWriter out = new PrintWriter ( new FileWriter(table));
			out.print("QUALThreshold");
			for (int j=0; j<mockVcfs.length; j++) out.print("\t"+Misc.removeExtension(mockVcfs[j].getName()));
			out.print("\tMockMean");
			for (int j=0; j<realVcfs.length; j++) out.print("\t"+Misc.removeExtension(realVcfs[j].getName()));
			out.println();
			float numMockDatasets = (float)mockVcfs.length;
			for (int i=0; i< scores.length; i++) {
				out.print(scores[i]);
				float total = 0;
				for (int j=0; j<mockVcfs.length; j++) {
					int num = countGreaterThanEqualTo (scores[i], mockQuals[j]);
					total += num;
					out.print("\t"+num);
				}
				float mean = total/numMockDatasets;
				mockThresholds[i] = scores[i];
				mockCounts[i] = mean;
				out.print("\t"+ mean);
				
				for (int j=0; j<realVcfs.length; j++) {
					int num = countGreaterThanEqualTo (scores[i], realQuals[j]);
					out.print("\t"+num);
				}
				out.println();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem saving count QUAL table to  "+table);
		}
	}
	
	public float estimateMockCounts(float threshold){
		int index = Arrays.binarySearch(mockThresholds,threshold);
		//no exact match
		boolean exact = true;
		if (index<0) {
			exact = false;
			index = (-1* index) -1;
			if (index > mockThresholds.length-1) {
				//IO.pl("Not exact past end returning last \t"+mockCounts[mockCounts.length-1]);
				return mockCounts[mockCounts.length-1];
			}
		}
		//look for smaller indexes with same value
		float value = mockThresholds[index];
		int testIndex = index;
		while (true){
			testIndex--;
			//look for negative index
			if (testIndex< 0) {
				index = 0;
				break;
			}
			//if different value break
			if (value != mockThresholds[testIndex]) {
				index = ++testIndex;
				break;
			}
		}
		if (exact) {
			//IO.pl("Exact index\t"+index+"\t"+mockCounts[index]);
			return mockCounts[index];
		}
		else {
			//must interpolate
			int leftIndex = index -1;
			float x1 = mockThresholds[leftIndex];
			float x2 = mockThresholds[index];
			float y1 = mockCounts[leftIndex];
			float y2 = mockCounts[index];
			float interpCounts = Num.interpolateY(x1, y1, x2, y2, threshold);
			//IO.pl("Interp\t"+leftIndex+"\t"+index+ "\t"+x1+"\t"+x2+ "\t"+y1+"\t"+y2+" final "+interpCounts);
			return interpCounts;
		}
	}

	static int countGreaterThanEqualTo(float threshold, float[] sortedQuals) {
		int numLessThan = Num.countNumLessThan(threshold, sortedQuals);
		return sortedQuals.length - numLessThan;
	}

	private void printMockHistogram(float min, float max) {
		IO.pl("\nMock QUAL score histogram:");
		mockQualScoreHistogram = new Histogram(min, max, 50);
		for (int i=0; i< mockVcfs.length; i++) mockQualScoreHistogram.countAll(mockQuals[i]);
		mockQualScoreHistogram.printScaledHistogram();
	}

	private float[] fetchSortedQualScores(File vcf) {
		BufferedReader in = null;
		String line = null;
		try {
			in = IO.fetchBufferedReader(vcf);
			ArrayList<Float> quals = new ArrayList<Float>();
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
					String[] fields = Misc.TAB.split(line);
					if (fields[5].equals(".") || fields[5].length() == 0) quals.add(0.0f);
					else quals.add(Float.parseFloat(fields[5]));
				}
			}
			float[] q = Num.arrayListOfFloatToArray(quals);
			Arrays.sort(q);
			return q;

		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing QUAL score from "+vcf+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
		return null;
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
		File toParseMocks = null;
		File toParseReal = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': toParseMocks = new File(args[++i]); break;
					case 'r': toParseReal = new File(args[++i]); break;
					case 's': saveDir = new File(args[++i]); break;
					case 'f': minFdr = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check args
		if (saveDir == null) Misc.printErrAndExit("\nError: please provide a directory to save the FDR scored vcf files.\n");
		saveDir.mkdirs();
		mockVcfs = fetchVcfFiles(toParseMocks);
		realVcfs = fetchVcfFiles(toParseReal);

	}
	
	public File[] fetchVcfFiles (File forExtraction) {
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such for "+forExtraction);
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		File[] vcfFiles = IO.collapseFileArray(tot);
		Arrays.sort(vcfFiles);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s) in "+forExtraction);
		return vcfFiles;
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           VCF Fdr Estimator  :   July 2020                       **\n" +
				"**************************************************************************************\n" +
				"This app estimates false discovery rates (FDR) for each somatic VCF by counting the\n"+
				"number of records that are >= to that QUAL score in one or more mock VCF.\n"+
				"The estimated FDR = #Mock/ #Real. In cases where increasingly stringent QUAL thresholds\n"+
				"reduce the nummber of real records but not mock records, the FDR increases. To control\n"+
				"for this inconsistancy, the prior FDR is assigned to the more stringent QUAL. Use\n"+
				"the USeq VCFBkz app to generate reliable QUAL ranking scores if your somatic call does.\n"+
				"not provide them.\n"+
				
				"\nTo generate mock somatic VCF files, run 3 or more mock tumor vs normal sample sets\n"+
				"where the mock tumor samples are either uninvolved tissue or biopsies from\n"+
				"healthy volunteers. Prepare and analyze them the same way as the real tumor samples.\n\n"+

				"Options:\n"+
				"-m Directory containing one or more mock VCF files (xxx.vcf(.gz/.zip OK)).\n"+
				"-r Directory containing one or more real VCF files to score\n"+
				"-s Directory for saving the dFDR scored vcf files. \n"+
				"-f Minimum -10Log10(FDR), defaults to 6, where 6=0.25, 10=0.1, 13=0.05, 20=0.01\n"+

				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VCFFdrEstimator -m MockVcfs/ -f 10"+
				"   -r RealVcfs/ -s FDRFilteredVcfs\n\n"+

				"**************************************************************************************\n");

	}
}
