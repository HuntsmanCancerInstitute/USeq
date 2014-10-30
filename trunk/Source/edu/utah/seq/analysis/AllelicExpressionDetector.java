package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import htsjdk.samtools.*;
import util.bio.annotation.Bed;
import util.gen.*;
import edu.utah.seq.analysis.multi.PairedCondition;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamLayoutForMutation;

/** Application for identifying allelic expression based on a table of snps and bam alignments filtered for alignment bias (see ReferenceMutator and SamComparator).
 * @author Nix
 * */
public class AllelicExpressionDetector {

	//user defined fields
	private File snpBedFile;
	private File[] bamFiles;
	private String[] sampleNamesToProcess;
	private String[] treatmentSampleNames;
	private String[] controlSampleNames;
	private File fullPathToR = new File ("/usr/bin/R");
	private File snpDataFile;
	private double minimumGenCallScore = 0.2;
	private int minimumReadCoverage = 4;
	private int minimumBaseQuality = 20;
	private int minimumReplicas = 3;
	private File resultsDirectory;
	private LinkedHashMap<String, Snp> nameSnp = new LinkedHashMap<String, Snp>();
	private HashMap<String, SamReader> sampleNameSamReader = new HashMap<String, SamReader>();
	private Random random = new Random(System.currentTimeMillis());
	private int minNumSnpsForDESeq = 50;
	
	//each round of DESeq2
	private File dataTable = null;
	private int numReplicas = 0;
	private int numSnps = 0;
	private String[] deseqResults;
	
	
	//constructor
	/**Stand alone.*/
	public AllelicExpressionDetector(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load snp info via bed file
		System.out.print("Loading snp info from bed file... ");
		loadSnpInfo();
		System.out.println(nameSnp.size());

		//create sam readers
		openSamReaders();

		//walk through snp data looking for snps with sufficient coverage
		System.out.print("Reading snp data and extracting alignments");
		walkSnpData();

		//close readers
		closeSamReaders();

		//write out count table
		System.out.println("Launching DESeq for snps with: ");
		numReplicas = minimumReplicas;
		while (true){
			System.out.print("\t"+numReplicas+" replicas ");
			writeOutCountTable(numReplicas);
			if (dataTable == null) {
				System.out.println(" too few!");
				break;
			}
			executeDESeq2();
			System.out.println();
			parseDESeq2Results();
			numReplicas++;
		}
		System.out.println("Writing filtered count table for DESeq2...");
		
		//write out summary
		System.out.println("Writing summary table...");
		writeOutSummaryTable();


		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes\n");
	}


	
	private void parseDESeq2Results() {
		//name	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
		//name chr20_47731157_kgp10629917_T_C_Ref_D005-14_D006-14_D040-13_D012-14_OS_Mac_RPE
		for (int i=1; i< deseqResults.length; i++){
			String[] t = Misc.TAB.split(deseqResults[i]);
			//String[] n = Misc.  here
		}
		
	}



	private void executeDESeq2() {
		try {
			File results = new File (resultsDirectory, "resultsDESeq2_"+numReplicas+"Reps.txt");
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("#load libraries\n");
			sb.append("library(DESeq2)\n");
			sb.append("countTable = read.delim('"+dataTable.getCanonicalPath()+"', header=FALSE)\n");
			sb.append("rownames(countTable) = countTable[,1]\n");
			sb.append("countTable = countTable[,-1]\n");
			
			//add column names
			sb.append("colnames(countTable) = c('A0','B0'");
			for (int i=1; i< numReplicas; i++) sb.append(",'A"+i+"','B"+i+"'");
			sb.append(")\n");
			
			//add condition info
			sb.append("sampleInfo = data.frame(condition=as.factor(c('A','B'");
			for (int i=1; i< numReplicas; i++) sb.append(",'A','B'");
			sb.append(")))\n");

			sb.append("rownames(sampleInfo) = colnames(countTable)\n");
			sb.append("cds = DESeqDataSetFromMatrix(countData=countTable, colData=sampleInfo, design = ~condition)\n");
			sb.append("cds = DESeq(cds)\n");
			sb.append("res = results(cds, contrast = c('condition', 'A', 'B'))\n");
			sb.append("res[,5] = -10 * log10(res[,5])\n"); //transform pvals since java can't handle some of these big numbers
			sb.append("res[,6] = -10 * log10(res[,6])\n"); 
			sb.append("write.table(res, file = '"+results+"', quote=FALSE, sep='\t')\n");

			//write script to file
			File scriptFile = new File (resultsDirectory,"rScriptDESeq2_"+numReplicas+"Reps.txt");
			File rOut = new File(resultsDirectory, "rScriptDESeq2_"+numReplicas+"Reps.ROut.txt");
			IO.writeString(sb.toString(), scriptFile);

			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};

			//execute command
			IO.executeCommandLine(command);

			//check for warnings?
			deseqResults = IO.loadFile(results);
			if ((deseqResults.length-1) != numSnps) Misc.printErrAndExit("\nERROR: the number of result rows ("+(deseqResults.length-1)+") does not match the number of snps ("+numSnps+")? See temp files for details.\n");


		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error: failed to execute DESeq2.\n");
		}
	}



	private void writeOutSummaryTable() {
		try {
		
			//make gzipper
			File f = new File (resultsDirectory, "summaryTable.xls");
			Gzipper out = new Gzipper(f);
			
			//for each Snp
			out.println("Chr\tPos\tRef\tAlt\tName\tSampleNames\tRefAltCounts\tSampleGCScores");
			for (Snp s : nameSnp.values()) {
				if (s.sampleName != null) out.println(s.toString());
			}
			
			out.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	
	
	private File writeOutCountTable(int numReps) {
		try {

			//make gzipper
			dataTable = new File (resultsDirectory, "countTableForDESeq2_"+numReps+"Reps.txt.gz");
			Gzipper out = new Gzipper(dataTable);
			
			//for each Snp
			boolean flip = true;
			HashSet<String> coordinates = new HashSet<String>();

			numSnps =0;
			for (Snp s : nameSnp.values()){
				
				//minSamples
				if (s.sampleName != null) {
					int numSamples = s.sampleName.size();
					if (numSamples < numReps) continue;
					
					//already been saved?
					String coordName = s.toStringCoordinates();
					if (coordinates.contains(coordName)) continue;
					coordinates.add(coordName);
					numSnps++;
					
					//flip alt vs ref counts
					if (flip) flip = false;
					else flip = true;

					CountPair[] cps = new CountPair[numSamples];
					//for each sample
					for (int i=0; i< numSamples; i++){
						int[] refAltCounts = s.getRefAltCounts(s.gatcnCounts.get(i));
						//filter here for cases with no evidence of being hetero?
						//flip?
						if (flip) cps[i] = new CountPair(refAltCounts[1], refAltCounts[0], s.sampleName.get(i));
						else cps[i] = new CountPair(refAltCounts[0], refAltCounts[1], s.sampleName.get(i));
					}
					
					//randomize CountPairs
					Misc.randomize(cps, random);

					//print line
					StringBuilder sb = new StringBuilder(coordName);
					if (flip) sb.append("_Alt");
					else sb.append("_Ref");
					//appendSampleNames
					for (int i=0; i< numReps; i++){
						sb.append("_");
						sb.append(cps[i].sampleName);
					}
					//for (int i=0; i< maxSamples; i++){
					for (int i=0; i< numReps; i++){
						sb.append("\t");
						sb.append(cps[i].numA);
						sb.append("\t");
						sb.append(cps[i].numB);
					}
					out.println(sb);
				}
			}
			out.close();
			
			//enough snps?
			System.out.print(numSnps);
			if (numSnps < minNumSnpsForDESeq ) {
				dataTable.delete();
				dataTable = null;
			}
			return dataTable;
			
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}

	private class CountPair{
		int numA = 0;
		int numB = 0;
		String sampleName;
		private CountPair(int a, int b, String sampleName){
			numA = a;
			numB = b;
			this.sampleName = sampleName;
		}
	}

	private int findMaxNumberSamples() {
		int max=0;
		for (Snp s : nameSnp.values()){
			if (s.sampleName != null && s.sampleName.size() > max) max = s.sampleName.size();
		}
		return max;
	}

	private void walkSnpData() {
		try {
			int counter = 0;
			BufferedReader in = IO.fetchBufferedReader(snpDataFile);
			String line;
			String[] tokens;
			HashMap<String, Character> counts = new HashMap<String, Character>();

			//skip header
			while((line = in.readLine())!=null){
				if (line.startsWith("SNP Name")) break;
			}
			if (line == null) Misc.printErrAndExit("\nError: failed to skip over the data file header?\n");

			//SNPName SampleID Allele1-Forward Allele2-Forward Allele1-Top Allele2-Top GCScore
			//   0       1            2               3             4           5         6
			while((line = in.readLine())!=null){
				if (counter++ > 500000) {
					System.out.print(".");
					counter = 0;
				}
				tokens = Misc.TAB.split(line);
				if (tokens.length !=7) continue;

				//is it a het?
				if (tokens[2].equals(tokens[3])) {
					//System.out.println("Not het!");
					continue;
				}

				//pass thresholds?
				double score = Double.parseDouble(tokens[6]);
				if (score < minimumGenCallScore)  {
					//System.out.println("Failed score!");
					continue;
				}

				//present in bed file data?
				Snp snp = nameSnp.get(tokens[0]);
				if (snp == null) {
					//System.out.println("Not in bed!");
					continue;
				}

				//one of the samples to process?
				String sampleName = tokens[1];
				SamReader sr = sampleNameSamReader.get(sampleName);
				if (sr == null) {
					//System.out.println("Not a sample to test!");
					continue;
				}

				//fetch alignments
				SAMRecordIterator it = sr.queryOverlapping(snp.chromosome, snp.position, snp.position+1);
				if (it.hasNext() == false) {
					//System.out.println("No overlapping reads!");
					it.close();
					continue;
				}

				counts.clear();
				while (it.hasNext()){
					SAMRecord sam = it.next();
					SamAlignment sa = new SamAlignment(sam.getSAMString().trim(), true);
					SamLayoutForMutation layout = new SamLayoutForMutation(sa);
					int index = layout.findBaseIndex(snp.position);
					if (index == -1) {
						//System.out.println("Base Not Found?");
						//System.out.println(sa);
						//System.out.println(snp);
						continue;
					}
					//good score?
					if (layout.getQual()[index] < minimumBaseQuality) {
						//System.out.println("Bad base!");
						continue;
					}
					//a match?
					if (layout.getCall()[index] != 'M'){
						//System.out.println("Not a M!");
						continue;
					}
					//add to hash
					char base = layout.getSeq()[index];

					counts.put(sa.getName(), new Character(base));
				}
				it.close();

				if (counts.size() == 0) continue;

				//count
				int[] gatcn = countHash(counts);
				if (Num.sumIntArray(gatcn) < minimumReadCoverage) continue;
				if (snp.sampleName == null) snp.initializeALs();
				snp.sampleName.add(sampleName);
				snp.gatcnCounts.add(gatcn);
				snp.sampleGCScore.add(new Float((float)score));

				//System.out.println(snp);
				//if (counter++ > 10) System.exit(0);
			}
			System.out.println();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private int[] countHash(HashMap<String, Character> counts) {
		int[] gatcn = new int[5];
		for (Character c : counts.values()){
			if (c == 'G') gatcn[0]++;
			else if (c == 'A') gatcn[1]++;
			else if (c == 'T') gatcn[2]++;
			else if (c == 'C') gatcn[3]++;
			else if (c == 'N') gatcn[4]++;
			else Misc.printErrAndExit("\nError: unrecognized base '"+c+"' from "+counts);
		}
		return gatcn;
	}

	private void closeSamReaders() {
		try {
			for (SamReader sr: sampleNameSamReader.values()) sr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void openSamReaders() {
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

		System.out.println("Opening SAM readers...");
		for (String sampleName: sampleNamesToProcess){
			File bam = new File (bamFiles[0].getParentFile(), sampleName+".bam");
			//attempt to modify bai to avoid stupid picard warning.
			File bai = new File (bamFiles[0].getParentFile(), sampleName+".bai");
			bai.setLastModified(System.currentTimeMillis());
			bai = new File (bam.toString()+".bai");
			bai.setLastModified(System.currentTimeMillis());
			System.out.println("\t"+bam.getName());
			if (bam.canRead() == false) Misc.printErrAndExit("\nError: cannot find the corresponding sample bam file? "+bam+"\n");
			SamReader sr = factory.open(bam);
			sampleNameSamReader.put(sampleName, sr);
		}
	}

	private void loadSnpInfo() {
		Pattern underscore = Pattern.compile("_");
		Bed[] bed = Bed.parseFile(snpBedFile, 0, 0);
		for (Bed b : bed){
			//parse name: name_reference_alternate
			String[] tokens = underscore.split(b.getName());
			Snp s = new Snp (b.getChromosome(), b.getStart(), tokens[0], tokens[1], tokens[2]);
			if (nameSnp.containsKey(tokens[0])) Misc.printErrAndExit("\nError: same snp name found, aborting! "+tokens[0]+"\n");
			nameSnp.put(tokens[0], s);
			Misc.printErrAndExit(tokens[0]);
		}
	}

	private class Snp{
		String chromosome;
		int position;
		String name;
		String reference;
		String alternate;
		ArrayList<String> sampleName;
		ArrayList<int[]> gatcnCounts;
		ArrayList<Float> sampleGCScore;

		private Snp(String chromosome, int position, String name, String reference, String alternate){
			this.chromosome = chromosome;
			this.position = position;
			this.name = name;
			this.reference = reference;
			this.alternate = alternate;
		}

		private int[] getRefAltCounts(int[] gatcn){
			int ref = 0; 
			int alt = 0;
			if (reference.equals("G")) ref = gatcn[0];
			else if (reference.equals("A")) ref = gatcn[1];
			else if (reference.equals("T")) ref = gatcn[2];
			else if (reference.equals("C")) ref = gatcn[3];

			if (alternate.equals("G")) alt = gatcn[0];
			else if (alternate.equals("A")) alt = gatcn[1];
			else if (alternate.equals("T")) alt = gatcn[2];
			else if (alternate.equals("C")) alt = gatcn[3];
			return new int[]{ref,alt};
		}

		private void initializeALs(){
			sampleName = new ArrayList<String>();
			gatcnCounts = new ArrayList<int[]>();
			sampleGCScore = new ArrayList<Float>();
		}

		public String toStringCoordinates(){
			StringBuilder sb = new StringBuilder();
			sb.append(chromosome); sb.append("_");
			sb.append(position); sb.append("_");
			sb.append(name); sb.append("_");
			sb.append(reference); sb.append("_");
			sb.append(alternate); 
			
			return sb.toString();
		}

		/**Chr\tPos\tRef\tAlt\tName\tSampleNames\tRefAltCounts\tSampleGCScores*/
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append(chromosome); sb.append("\t");
			sb.append(position); sb.append("\t");
			sb.append(reference); sb.append("\t");
			sb.append(alternate); sb.append("\t");
			sb.append(name); 
			if (sampleName != null){
				sb.append("\t");
				//sampleNames
				String sampleNames = Misc.stringArrayListToString(sampleName, ",");
				sb.append(sampleNames); sb.append("\t");
				
				//ref/alt counts
				int[] refAltCounts = getRefAltCounts(gatcnCounts.get(0));
				sb.append(refAltCounts[0]); sb.append("_"); sb.append(refAltCounts[1]); 
				for (int i=1; i< gatcnCounts.size(); i++){
					sb.append(",");
					refAltCounts = getRefAltCounts(gatcnCounts.get(i));
					sb.append(refAltCounts[0]); sb.append("_"); sb.append(refAltCounts[1]);
				}
				sb.append("\t");
				
				//sample scores
				float[] s = Num.arrayListOfFloatToArray(sampleGCScore);
				String scores = Misc.floatArrayToString(s, ",");
				sb.append(scores);
			}
			return sb.toString();
		}

		public boolean passMinCoverage(){
			for (int i=0; i< gatcnCounts.size(); i++){
				if (Num.sumIntArray(gatcnCounts.get(i)) >= minimumReadCoverage) return true;
			}
			return false;
		}
	}

	public static void main(String[] args) {

		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AllelicExpressionDetector(args);
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
					case 'b': bamFiles = IO.extractFiles(new File(args[++i]), ".bam"); break;
					case 's': snpBedFile = new File(args[++i]); break;
					case 't': treatmentSampleNames = args[++i].split(","); break;
					case 'c': controlSampleNames = args[++i].split(","); break;
					case 'n': sampleNamesToProcess = args[++i].split(","); break;
					case 'd': snpDataFile = new File(args[++i]); break;
					case 'g': minimumGenCallScore = Double.parseDouble(args[++i]); break;
					case 'm': minimumReadCoverage = Integer.parseInt(args[++i]); break;
					case 'q': minimumBaseQuality = Integer.parseInt(args[++i]); break;
					case 'e': resultsDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (bamFiles == null || bamFiles[0].canRead() == false) Misc.printErrAndExit("\nError: can't find your bam alignment files?\n");
		if (snpBedFile == null || snpBedFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your SNP Map bed file generated by the ReferenceMutator?\n");
		if (snpDataFile == null || snpDataFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your SNP data file?\n");
		
		//if (treatmentSampleNames == null || controlSampleNames == null) Misc.printErrAndExit("\nError: please enter the name(s) of samples to process with -t -c\n");
		//sampleNamesToProcess = Misc.copyAndMerge(treatmentSampleNames, controlSampleNames);
		
		if (resultsDirectory == null) Misc.printErrAndExit("\nError: please enter a directory in which to save the results.\n");
		resultsDirectory.mkdirs();
		if (resultsDirectory.canWrite() == false || resultsDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: please enter a directory in which to save the results.\n");

		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}

		String errors = IO.runRCommandLookForError("library(DESeq2)", fullPathToR, resultsDirectory);
		if (errors == null || errors.length() !=0){
			Misc.printErrAndExit("\nError: Cannot find the required R library.  Did you install DESeq2?. R error message:\n\t\t"+errors+"\n\n");
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Allelic Expression Detector:  Oct 2014                      **\n" +
				"**************************************************************************************\n" +
				"Application for identifying allelic expression based on a table of snps and bam\n"+
				"alignments that have been filtered for alignment bias.  See the ReferenceMutator and\n"+
				"SamComparator apps.\n\n"+

				"Required Options:\n"+
				"-t Treatment sample names to process, comma delimited, no spaces.\n"+
				"-c Control sample names to process, comma delimited, no spaces.\n"+
				"-b Directory containing coordinate sorted bam and index files named according to their\n"+
				"      sample name.\n"+
				"-d SNP data file containing all sample snp calls.\n"+
				"-e Results directory.\n"+
				"-s SNP map bed file from the ReferenceMutator app.\n"+

				"\nDefault Options:\n"+
				"-g Minimum GenCall score, defaults to >= 0.25\n"+
				"-q Minimum alignment base quality at snp, defaults to 20\n"+
				"-m Minimum alignment read coverage, defaults to 4\n"+
				"-r Full path to R (version 3+) loaded with DESeq2, see http://www.bioconductor.org\n"+
				"       Type 'library(DESeq2) in R to see if it is installed. Defaults to '/usr/bin/R'\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/AllelicExpressionDetector -b Bam/RPENormal/\n"+
				"-n D002-14,D005-14,D006-14,D009-14 -d GenotypingResults.txt.gz -s SNPMap_Ref2Alt_Int.txt\n"+
				"-r RPENormal\n\n" +

				"**************************************************************************************\n");

	}
}
