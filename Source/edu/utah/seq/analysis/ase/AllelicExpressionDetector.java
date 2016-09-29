package edu.utah.seq.analysis.ase;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import htsjdk.samtools.*;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.*;
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
	private File fullPathToR = new File ("/usr/bin/R");
	private File snpDataFile;
	private double minimumGenCallScore = 0.25;
	private int minimumReadCoverage = 4;
	private int minimumBaseQuality = 20;
	private int minimumReplicas = 3;
	private File resultsDirectory;
	private LinkedHashMap<String, AllelicExpressionSnp> nameSnp = new LinkedHashMap<String, AllelicExpressionSnp>();
	private HashMap<String, SamReader> sampleNameSamReader = new HashMap<String, SamReader>();
	private Random random = new Random(System.currentTimeMillis());
	private int minNumSnpsForDESeq = 500;
	private boolean deleteTempFiles = true;
	private File tabixBedFile;

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

		//run iterative DESeq2
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
		//correct FDRs
		correctFDRs();

		//add gene info
		annotateWithGenes();

		//write out summary
		System.out.println("Writing summary tables...");
		writeOutSummaryTable();
		writeOutGeneiASETable();

		//save nameSnp hash
		System.out.println("Saving serialized name:Snp hash object...");
		stripAndSave();


		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes\n");
	}

	private void annotateWithGenes() {
		TabixReader tr;
		try {
			tr = new TabixReader( tabixBedFile.toString() );


			for (AllelicExpressionSnp s : nameSnp.values()) {
				if (s.getSampleName() == null) continue;
				int pos = s.getPosition() -1;
				int end = pos + 2;
				String q = s.getChromosome()+":"+pos+"-"+end;
				TabixReader.Iterator it = tr.query(q);
				if (it != null) {
					String bedLine;
					ArrayList<String> names = new ArrayList<String>();
					while ((bedLine = it.next()) != null) {
						//chr begin end name score strand
						String[] tokens = Misc.TAB.split(bedLine);
						String geneName = tokens[3];
						names.add(geneName);
					}
					if (names.size() !=0) s.setGeneNames(names);
				}
			}
			tr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}



	private void stripAndSave() {
		ArrayList<String> toRemove = new ArrayList<String>();
		for (String name: nameSnp.keySet()){
			AllelicExpressionSnp aes = nameSnp.get(name);
			if (aes.getSampleName() == null) toRemove.add(name);
		}
		for (String s: toRemove) nameSnp.remove(s);
		File ns = new File (resultsDirectory, "nameSnp.obj");
		IO.saveObject(ns, nameSnp);
		System.out.println("\t"+nameSnp.size()+" data snps saved");

	}

	/**Fixes issue with DESeq FDRs*/
	private void correctFDRs() {
		//find max FDR but not the biggie
		float maxValue = 0;
		for (AllelicExpressionSnp snp : nameSnp.values()){
			if (snp.getFdr() != Float.MAX_VALUE && snp.getFdr() > maxValue) maxValue = snp.getFdr();
		}
		//set FDR for all biggies
		for (AllelicExpressionSnp snp : nameSnp.values()){
			if (snp.getFdr() == Float.MAX_VALUE) snp.setFdr(maxValue); 
		}

	}



	private void parseDESeq2Results() {
		//name	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
		//  0       1             2           3       4       5       6
		//chr20_47731157_kgp10629917_T_C_Ref_D005-14_D006-14_D040-13_D012-14_OS_Mac_RPE
		//  0      1          2      3 4  5  6   7.....
		for (int i=1; i< deseqResults.length; i++){
			String[] t = Misc.TAB.split(deseqResults[i]);
			String[] n = Misc.UNDERSCORE.split(t[0]);
			//fetch Snp
			AllelicExpressionSnp snp = nameSnp.get(n[2]);

			//parse fdr, watch for bad values
			float currFDR = 0;
			if (t[6].equals("Inf")) currFDR = Float.MAX_VALUE;
			else if (t[6].equals("NA")) currFDR = 0f;
			else {
				currFDR = Float.parseFloat(t[6]);
			}
			//modify
			if (currFDR > snp.getFdr()) {
				snp.setFdr(currFDR);
				snp.setLog2Ratio(Float.parseFloat(t[2]));
			}
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

			if (deleteTempFiles){
				results.deleteOnExit();
				scriptFile.deleteOnExit();
				rOut.deleteOnExit();
			}

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
			out.println("Chr\tPos\tRef\tAlt\tName\tSampleNames\tRefAltCounts\tTotalRef\tTotalAlt\tSampleGCScores\tFDR\tLog2Rto\tBadHetCall?\tGeneNames");
			for (AllelicExpressionSnp s : nameSnp.values()) {
				if (s.getSampleName() != null) out.println(s.toString(false));
			}

			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void writeOutGeneiASETable() {

		try {
			//make gzipper
			Gzipper out = new Gzipper(new File (resultsDirectory, "geneiASETable.xls"));
			//header
			out.println("gene\tsnp.id\talt.dp\tref.dp");
			//use to keep the snp.id column unique even though some are shared between genes
			int snpId = 0;
			//for each gene snp
			for (AllelicExpressionSnp s : nameSnp.values()) {
				//skip those with no sample
				if (s.getSampleName() == null) continue;

				int[] refAlt = s.getTotalRefAltCounts();

				//any genes?
				ArrayList<String> geneNames = s.getGeneNames();
				if (geneNames != null){
					//make an entry for each gene name
					for (String gn : geneNames){
						out.print(gn); out.print("\t");
						out.print(s.getChrPosRefAlt()); out.print("_"); out.print(snpId++);out.print("\t");
						out.print(refAlt[1]); out.print("\t");
						out.print(refAlt[0]); out.println("\t");
					}
				}
				else {
					String noName = s.getChrPosRefAlt();
					out.print(noName); out.print("\t");
					out.print(noName); out.print("_"); out.print(snpId++); out.print("\t");
					out.print(refAlt[1]); out.print("\t");
					out.print(refAlt[0]); out.println("\t");
				}
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
			if (deleteTempFiles) dataTable.deleteOnExit();
			Gzipper out = new Gzipper(dataTable);

			//for each Snp
			boolean flip = true;
			HashSet<String> coordinates = new HashSet<String>();

			numSnps =0;
			for (AllelicExpressionSnp s : nameSnp.values()){

				//minSamples
				if (s.getSampleName() != null) {
					int numSamples = s.getSampleName().size();
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
						int[] refAltCounts = s.getRefAltCounts(s.getGatcnCounts().get(i));

						//filter here for cases with no evidence of being hetero?
						//flip?
						if (flip) cps[i] = new CountPair(refAltCounts[1], refAltCounts[0], s.getSampleName().get(i));
						else cps[i] = new CountPair(refAltCounts[0], refAltCounts[1], s.getSampleName().get(i));
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
		for (AllelicExpressionSnp s : nameSnp.values()){
			if (s.getSampleName() != null && s.getSampleName().size() > max) max = s.getSampleName().size();
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
				AllelicExpressionSnp snp = nameSnp.get(tokens[0]);
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
				SAMRecordIterator it = sr.queryOverlapping(snp.getChromosome(), snp.getPosition(), snp.getPosition()+1);
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
					int index = layout.findBaseIndex(snp.getPosition());
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
				if (snp.getSampleName() == null) snp.initializeALs();
				snp.getSampleName().add(sampleName);
				snp.getGatcnCounts().add(gatcn);
				snp.getSampleGCScore().add(new Float((float)score));

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
			AllelicExpressionSnp s = new AllelicExpressionSnp (b.getChromosome(), b.getStart(), tokens[0], tokens[1], tokens[2]);
			if (nameSnp.containsKey(tokens[0])) Misc.printErrAndExit("\nError: same snp name found, aborting! "+tokens[0]+"\n");
			nameSnp.put(tokens[0], s);
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
					case 't': tabixBedFile = new File(args[++i]); break;
					case 'n': sampleNamesToProcess = args[++i].split(","); break;
					case 'd': snpDataFile = new File(args[++i]); break;
					case 'g': minimumGenCallScore = Double.parseDouble(args[++i]); break;
					case 'm': minimumReadCoverage = Integer.parseInt(args[++i]); break;
					case 'q': minimumBaseQuality = Integer.parseInt(args[++i]); break;
					case 'p': minimumReplicas = Integer.parseInt(args[++i]); break;
					case 'e': resultsDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
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
		if (tabixBedFile == null || tabixBedFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your tabixed bed file?\n");

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
				"**                      Allelic Expression Detector:  Sept 2016                     **\n" +
				"**************************************************************************************\n" +
				"Application for identifying allelic expression based on a table of snps and bam\n"+
				"alignments that have been filtered for alignment bias.  See the ReferenceMutator and\n"+
				"SamComparator apps. Uses DESeq2 to identify differential expression between alleles.\n\n"+

				"Required Options:\n"+
				//"-t Treatment sample names to process, comma delimited, no spaces.\n"+
				//"-c Control sample names to process, comma delimited, no spaces.\n"+
				"-n Sample names to process, comma delimited, no spaces.\n"+
				"-b Directory containing coordinate sorted bam and index files named according to their\n"+
				"      sample name.\n"+
				"-d SNP data file containing all sample snp calls.\n"+
				"-e Results directory.\n"+
				"-s SNP map bed file from the ReferenceMutator app.\n"+
				"-t Tabix gz indexed bed file of exons where the name column is the gene name, see\n"+
				"      ExportExons and https://github.com/samtools/htslib\n"+

				"\nDefault Options:\n"+
				"-g Minimum GenCall score, defaults to 0.25\n"+
				"-q Minimum alignment base quality at snp, defaults to 20\n"+
				"-m Minimum alignment read coverage, defaults to 4\n"+
				"-p Minimum number replicas with heterozygous snp to score, defaults to 3\n"+
				"-r Full path to R (version 3+) loaded with DESeq2, see http://www.bioconductor.org\n"+
				"       Type 'library(DESeq2) in R to see if it is installed. Defaults to '/usr/bin/R'\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/AllelicExpressionDetector -b Bam/RPENormal/\n"+
				"-n D002-14,D005-14,D006-14,D009-14 -d GenotypingResults.txt.gz -s SNPMap_Ref2Alt_Int.txt\n"+
				"-r RPENormal -t ~/Anno/b37EnsGenes7Sept2016_Exons.bed.gz\n\n" +

				"**************************************************************************************\n");

	}
}
