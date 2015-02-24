package edu.utah.seq.cnv;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import htsjdk.samtools.*;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.analysis.DefinedRegionDifferentialSeq;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.analysis.multi.GeneCount;
import edu.utah.seq.analysis.multi.Replica;
import edu.utah.seq.data.HeatMapMaker;
import edu.utah.seq.data.HeatMapMakerPosNeg;
import edu.utah.seq.data.Info;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.apps.Bar2USeq;
import edu.utah.seq.useq.data.PositionScore;


/** App to detect CNV's in case and tumor-normal datasets.  Wrapps Alun Thomas's algorithm.
 * @author Nix
 * */
public class PoisRegCNV {

	//user defined fields
	private File saveDirectory;
	private File bamDirectory;
	private String[] caseSampleFileNames;
	private String[] controlSampleFileNames;
	private File fullPathToR = new File ("/usr/bin/R");
	private File refSeqFile;
	private int minimumCounts = 20;
	private int maxAlignmentsDepth = 50000;
	private File alunRScript = null;
	private boolean deleteTempFiles = false;
	private String genomeVersion;
	private float minimumLog2Ratio = 1;
	private float minimumAdjPVal = 0.01f;
	private int numberChunks = 20;

	//internal fields
	private File graphDirectory;
	private UCSCGeneLine[] genes;
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private HashMap<String, UCSCGeneLine> name2Gene;
	public static final Pattern BAD_NAME = Pattern.compile("(.+)/[12]$");
	private HashSet<String> flaggedGeneNames = new HashSet<String>();
	private LinkedHashSet<String> geneNamesWithMinimumCounts = new LinkedHashSet<String>();
	private String[] geneNamesToAnalyze = null;
	private Replica[] cases = null;
	private Replica[] controls = null;

	//container for cases and their exons that pass the minimumCounts threshold
	private GeneExonSample[] ges = null;
	private File countTableFile;
	private int numExonsProcessed = 0;
	private int numberPassingExons = 0;

	//from R
	private File residualsFile = null;
	private File obsExpLg2RtosFile = null;
	private File sigLevelFile = null;
	private float significanceThreshold;

	//for loading data
	private int workingFragmentNameIndexPlus = 1;
	private int workingFragmentNameIndexMinus = -1;
	private HashMap<String, Integer> workingFragNameIndex = new HashMap<String, Integer>(10000);

	//constructors
	/**Stand alone.*/
	public PoisRegCNV(String[] args){	
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);
		//launch
		run();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}


	public void run(){
		//run in non X11 mode
		System.setProperty("java.awt.headless", "true");

		//load gene models
		System.out.println("Loading regions/ gene models...");
		loadGeneModels();

		//load samples 
		loadSamples();

		//write out count table of genes passing minimum counts
		if (controls != null) writeCountTableMatrix();
		else {
			System.out.println("Writing count table and executing cnv analysis in R...");
			writeCaseCountTable();
			
			//execute Alun's RScript
			executeCaseAnalysisScript();
			
			//load results int0 GeneExonSample[]
			System.out.println("Loading results...");
			loadCaseAnalysisResults();
			
			//write out graph files for each sample
			System.out.println("Saving graphs and spreadsheets...");
			exportCaseSampleGraphs();
			
			//write out spreadsheet
			exportGeneExonSpreadSheets();
			
			System.out.println("\n"+numExonsProcessed+" Processed exons, "+numberPassingExons+" passed thresholds ("+minimumLog2Ratio+" Lg2Rto, "+minimumAdjPVal+" AdjPVal)\n");
		}
		
		//convert graphs to useq format
		new Bar2USeq(graphDirectory, true);
	}

	private void exportGeneExonSpreadSheets() {
		try {
		//collect all samples by gene
		String currGeneName = ges[0].getGene().getDisplayName();
		ArrayList<GeneExonSample> al = new ArrayList<GeneExonSample>();
		al.add(ges[0]);
		
		Gzipper outAll = new Gzipper(new File(saveDirectory, "resultsAll.xls.gz"));
		Gzipper outPass = new Gzipper(new File(saveDirectory, "resultsPass.xls.gz"));
		
		//print headers
		String gn = "Gene Name";
		if (ges[0].getGene().getName() != null) gn = gn+"\tDescription";
		StringBuilder sb = new StringBuilder(gn+"\tExon Coordinates\tExon Index\tIGB Link\tIGV Link\tPassing Samples");
		for (int i=0; i< cases.length; i++){
			String name = cases[i].getNameNumber();
			sb.append("\t"+name+" Lg2(Ob/Ex)\t"+name+" Res\t"+name+" Counts");
		}
		sb.append("\tLog2Rto "+minimumLog2Ratio);
		sb.append("\tResidual "+minimumAdjPVal+" sig threshold +/- "+significanceThreshold);
		
		outAll.println(sb);
		outPass.println(sb);
		
		for (int i=1; i< ges.length; i++){
			String nextName = ges[i].getGene().getDisplayNameThenName();
			if (nextName.equals(currGeneName) == false){
				printGeneSamples(outAll, outPass, al);
				currGeneName = nextName;
				al.clear();
			}
			al.add(ges[i]);
		}
		
		outAll.close();
		outPass.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nError printing spreadsheets.\n");
		}
	}


	private void printGeneSamples(Gzipper outAll, Gzipper outPass, ArrayList<GeneExonSample> al) throws IOException {
		//print geneName, exon coordinates, exon index, igb, igv, passsingSampleNames, sample1_Lg2Rt, sample1_Res, sample1_counts, sample2_.....
		GeneExonSample g = al.get(0);
		int exonIndex = g.getExonIndex();
		
		StringBuilder sb = new StringBuilder (fetchGeneInfo(g));
		//any samples that pass?
		String passingSamples = exonPassThresholds(exonIndex, al);
		if (passingSamples != null) sb.append (passingSamples+"\t");
		else sb.append("\t");
		
		//add first sample
		sb.append(g.getDataString());
		
		//for each subsequent sample
		for (int i=1; i< al.size(); i++){
			GeneExonSample t = al.get(i);
			
			//different exon index?
			if (t.getExonIndex() != exonIndex){
				//print out old
				outAll.println(sb);
				if (passingSamples != null) outPass.println(sb);
				
				//reset for new
				exonIndex = t.getExonIndex();
				sb = new StringBuilder (fetchGeneInfo(t));
				passingSamples = exonPassThresholds(exonIndex, al);
				if (passingSamples != null) sb.append (passingSamples);
			}
			
			sb.append("\t");
			sb.append(t.getDataString());
		}
		
		//print last
		outAll.println(sb);
		if (passingSamples != null) outPass.println(sb);
	}
	
	/**Returns null if no sample exon passes filters, otherwise comma delimited sample name list of those that pass.*/
	public String exonPassThresholds(int exonIndex, ArrayList<GeneExonSample> al){
		ArrayList<String> sampPass = new ArrayList<String>();
		for (int i=0; i<al.size(); i++){
			GeneExonSample g = al.get(i);
			//correct index?
			if (g.getExonIndex() == exonIndex){
				if (passesThresholds(g)) sampPass.add(cases[g.getSampleIndex()].getNameNumber());
			}
		}
		if (sampPass.size() !=0) {
			numberPassingExons++;
			return Misc.stringArrayListToString(sampPass, ",");
		}
		return null;
	}
	
	public String fetchGeneInfo(GeneExonSample g){
		UCSCGeneLine gene = g.getGene();
		ExonIntron[] exons = gene.getExons();
		String chr = gene.getChrom();
		int exonIndex = g.getExonIndex();
		//name
		StringBuilder sb = new StringBuilder (g.getGene().getDisplayName()); sb.append("\t");
		//description
		if (g.getGene().getName() != null) sb.append(g.getGene().getName()); sb.append("\t");
		//exon coordinates
		sb.append(chr);
		sb.append(":");
		sb.append(exons[exonIndex].getStartStopString()); sb.append("\t");
		//exon index
		sb.append(exonIndex); sb.append("\t");
		//igb link
		sb.append(fetchIGBLink(chr, exons[exonIndex])); sb.append("\t");
		//igv link
		sb.append(fetchIGVLink(chr, exons[exonIndex])); sb.append("\t");
		
		return sb.toString();
	}

	public static String fetchIGVLink(String chr, ExonIntron exon){
		int start = exon.getStart()-5000;
		if (start < 0) start = 0;
		int end = exon.getEnd()+5000;
		return "=HYPERLINK(\"http://localhost:60151/goto?locus="+chr+":"+start+ "-" + end+"\",\"IGV\")";
	}
	
	public String fetchIGBLink(String chr, ExonIntron exon){
		int start = exon.getStart()-5000;
		if (start < 0) start = 0;
		int end = exon.getEnd()+5000;
		return "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid="+chr+"&start="+start+"&end="+end+ "\",\"IGB\")";
	}

	private boolean passesThresholds(GeneExonSample g) {
		if (Math.abs(g.getObsExpLgRto()) < minimumLog2Ratio || Math.abs(g.getResidual()) < significanceThreshold ) return false;
		return true;
	}


	private void exportCaseSampleGraphs() {
		//split ges by sample
		GeneExonSample[][] sampleGES = new GeneExonSample[cases.length][numExonsProcessed];
		for (int i=0; i< ges.length; i++){
			sampleGES[ges[i].getSampleIndex()][ges[i].getGlobalExonIndex()] = ges[i];
		}
		//make dirs
		graphDirectory = new File (saveDirectory, "DataTracks");
		graphDirectory.mkdir();
		File resDir = new File (graphDirectory, "Residual");
		resDir.mkdir();
		File rtoDir = new File (graphDirectory, "ObsExpLog2Rto");
		rtoDir.mkdir();
		
		//for each write files
		for (int i=0; i< sampleGES.length; i++){
			//make sample specific dirs
			File resSampDir = new File (resDir, "res_"+cases[i].getNameNumber());
			resSampDir.mkdir();
			File rtoSampDir = new File (rtoDir, "eo_"+cases[i].getNameNumber());
			rtoSampDir.mkdir();
			saveGraphs(cases[i].getNameNumber(), sampleGES[i], resSampDir, rtoSampDir);
		}
	}


	private void saveGraphs(String sampleName, GeneExonSample[] ges, File resDir, File rtoDir) {
		//walk through data and make a SmoothingWindow[] for each chrom
		String currChr = ges[0].getGene().getChrom();
		ArrayList<SmoothingWindow> smAL = new ArrayList<SmoothingWindow>();
		for (int i=0; i< ges.length; i++){
			//make an sm
			ExonIntron exon = ges[i].getGene().getExons()[ges[i].getExonIndex()];
			SmoothingWindow sm = new SmoothingWindow(exon.getStart(), exon.getEnd(), new float[]{ges[i].getResidual(), ges[i].getObsExpLgRto()});
			
			//diff chrom?
			if (ges[i].getGene().getChrom().equals(currChr) == false) {
				//diff chrom, write out old
				Info info = new Info(sampleName, genomeVersion, currChr, ".", 0, null);
				SmoothingWindow[] sms = new SmoothingWindow[smAL.size()];
				smAL.toArray(sms);
				saveStairStepGraph(0, sms, info, resDir, "#FF0000", true); //red
				saveStairStepGraph(1, sms, info, rtoDir, "#FFFF00", true); //yellow
				
				//reset and set new
				currChr = ges[i].getGene().getChrom();
				smAL.clear();
			}
			
			smAL.add(sm);
		}
		//process last!
		Info info = new Info(sampleName, genomeVersion, currChr, ".", 0, null);
		SmoothingWindow[] sms = new SmoothingWindow[smAL.size()];
		smAL.toArray(sms);
		saveStairStepGraph(0, sms, info, resDir, "#FF0000", true); //red
		saveStairStepGraph(1, sms, info, rtoDir, "#FFFF00", true); //yellow
	}
	
	/**Saves bar heatmap/ stairstep graph files*/
	public void saveStairStepGraph (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color, boolean posNeg){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//save in info
		info.setNotes(map);
		//get heatmap positions and values
		PointData pd;		
		if (posNeg){
			HeatMapMakerPosNeg hm = new HeatMapMakerPosNeg(scoreIndex, 0, 0);
			pd = hm.makeHeatMapPositionValues(sm);
		}
		else {
			HeatMapMaker hm = new HeatMapMaker(scoreIndex,0);
			pd = hm.makeHeatMapPositionValues(sm, false);
		}
		pd.setInfo(info);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
	}



	private void loadCaseAnalysisResults() {
		//load data
		float[] residuals = Num.loadFloats(residualsFile);
		float[] obsExpRtos = Num.loadFloats(obsExpLg2RtosFile);
		float[] sigLevelArray = Num.loadFloats(sigLevelFile);
		//check
		if (residuals == null || residuals.length != ges.length || obsExpRtos == null || obsExpRtos.length != ges.length || sigLevelArray == null || sigLevelArray.length !=1){
			Misc.printErrAndExit("\nError: cannot load appropriate data from R results, check logs.\n");
		}
		//load em
		for (int i=0; i< ges.length; i++){
			ges[i].setResidual(residuals[i]);
			ges[i].setObsExpLgRto(obsExpRtos[i]);
		}
		significanceThreshold = sigLevelArray[0];
		
		//clean up
		if (deleteTempFiles){
			residualsFile.deleteOnExit();
			obsExpLg2RtosFile.deleteOnExit();
			sigLevelFile.deleteOnExit();
		}
	}


	private void executeCaseAnalysisScript() {
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("source('"+ alunRScript+ "')\n");
			sb.append("x = readdata('"+countTableFile+"',3,3)\n");
			sb.append("x = fitmodel(x)\n");
			sb.append("exportCNVData(x, '"+saveDirectory +"', testsig="+minimumAdjPVal+" )\n");

			//write script to file
			File scriptFile = new File (saveDirectory,"cnv_RScript.R");
			File rOut = new File(saveDirectory, "cnv_RScript.Rout");
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

			//look for results files
			residualsFile = new File (saveDirectory, "residuals.txt");
			obsExpLg2RtosFile = new File (saveDirectory, "obsExpLg2Rtos.txt");
			sigLevelFile = new File (saveDirectory, "sigLevel.txt");
			if (residualsFile.exists() == false || obsExpLg2RtosFile.exists() == false || sigLevelFile.exists() == false){
				Misc.printErrAndExit("\nError: cannot find the R results files? Check the cnv_RScript.Rout log file for errors.\n");
			}

			//cleanup
			if (deleteTempFiles) {
				countTableFile.deleteOnExit();
				rOut.deleteOnExit();
				scriptFile.deleteOnExit();
			}

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error: failed to execute case only cnv analysis in R.\n");
		}


	}


	private void loadSamples(){

		//load cases
		System.out.println("\nLoading case samples...");
		cases = new Replica[caseSampleFileNames.length];
		String[] shortNames = new String[cases.length];
		for (int i=0; i< caseSampleFileNames.length; i++){
			File bam = new File(bamDirectory, caseSampleFileNames[i]);
			String name = Misc.removeExtension(caseSampleFileNames[i]);
			shortNames[i] = name;
			cases[i] = new Replica(name, bam);
			System.out.print("\t"+caseSampleFileNames[i]);
			loadReplica(cases[i]);
			System.out.println("\t"+cases[i].getTotalCounts());
		}
		shortNames = Misc.trimCommon(shortNames);
		for (int i=0; i< cases.length; i++) cases[i].setNameNumber(shortNames[i]);

		//load controls
		if (controlSampleFileNames != null){
			System.out.println("Loading control samples...");
			shortNames = new String[controls.length];
			controls = new Replica[this.controlSampleFileNames.length];
			for (int i=0; i< controlSampleFileNames.length; i++){
				File bam = new File(bamDirectory, controlSampleFileNames[i]);
				String name = Misc.removeExtension(controlSampleFileNames[i]);
				shortNames[i] = name;
				controls[i] = new Replica(name, bam);
				System.out.print("\t"+controlSampleFileNames[i]);
				loadReplica(controls[i]);
				System.out.println("\t"+controls[i].getTotalCounts());
			}
			shortNames = Misc.trimCommon(shortNames);
			for (int i=0; i< cases.length; i++) controls[i].setNameNumber(shortNames[i]);
		}

		//any genes with too many reads that were excluded?
		if (flaggedGeneNames.size() !=0) {
			String[] badGeneNames = Misc.hashSetToStringArray(flaggedGeneNames);			
			for (Replica replica: cases) replica.removeFlaggedGenes(badGeneNames);
			if (controlSampleFileNames != null) for (Replica replica: controls) replica.removeFlaggedGenes(badGeneNames);
			//remove them from geneNamesWithMinimumCounts
			for (String baddie: badGeneNames) geneNamesWithMinimumCounts.remove(baddie);
		}

		System.out.println();
		System.out.println(geneNamesWithMinimumCounts.size()+" genes, with >= "+minimumCounts+" and < "+maxAlignmentsDepth+" counts, will be examined for copy number alterations.\n ");
		System.out.println("Skipped genes: "+flaggedGeneNames);
		geneNamesToAnalyze = Misc.hashSetToStringArray(geneNamesWithMinimumCounts);
	}

	public void loadBlocks(SAMRecord sam, int chrStartBp, ArrayList<Integer>[] bpNames, boolean[] badBases){

		NameInteger nameIndex = null;
		boolean addIt;
		ArrayList<int[]> blocks = DefinedRegionDifferentialSeq.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
		//add name to each bp
		for (int[] b : blocks){
			int start = b[0] - chrStartBp;
			int stop = b[1] - chrStartBp;
			//need to watch for out of bounds issues, sometimes the length of the chromosome is incorrect in the bam header.
			if (stop > badBases.length) stop = badBases.length;
			for (int i=start; i < stop; i++){
				//bad base?
				if (badBases[i]) continue;
				addIt = true;

				//never seen before?
				if (bpNames[i] == null) bpNames[i] = new ArrayList<Integer>();
				//old so check size
				else if (bpNames[i].size() == maxAlignmentsDepth) {
					badBases[i] = true;
					addIt = false;
					bpNames[i] = null;	
					if (nameIndex !=null) workingFragNameIndex.remove(nameIndex.name);
				}

				// add it?
				if (addIt) {
					if (nameIndex == null) nameIndex = fetchFragmentNameIndex(sam);
					bpNames[i].add(nameIndex.index);
				}
			}
		}
	}


	/**Fetches an old or makes a new Integer to represent the sam read name (e.g. fragment name)*/
	public NameInteger fetchFragmentNameIndex(SAMRecord sam){
		String samReadName = sam.getReadName();
		Integer index;
		Matcher mat = BAD_NAME.matcher(samReadName);
		if (mat.matches()) samReadName = mat.group(1);
		if (workingFragNameIndex.containsKey(samReadName)) index = workingFragNameIndex.get(samReadName);
		else {
			if (sam.getReadNegativeStrandFlag()) {
				index = new Integer(workingFragmentNameIndexMinus--);
			} else {
				index = new Integer(workingFragmentNameIndexPlus++);
			}
			workingFragNameIndex.put(samReadName, index);
		}
		return new NameInteger (samReadName, index);
	}

	private class NameInteger {
		String name;
		Integer index;

		public NameInteger(String name, Integer index){
			this.name = name;
			this.index = index;
		}
	}

	public void loadReplica(Replica replica){
		//the idea here is to take a sorted bam file and add the reads to an ArrayList, basically every base should have a ArrayList of reads that overlap that base
		//one can then count the number of overlapping fragments for an exon or a collection of exons by hashing the ArrayList
		//assumes each read (first or second) has the same name

		//make reader
		SAMFileReader reader = new SAMFileReader(replica.getBamFile());	
		reader.setValidationStringency(ValidationStringency.SILENT);


		//fetch chromName: length for all chroms
		HashMap<String, Integer> chromLength = new HashMap<String, Integer>();
		List<SAMSequenceRecord> seqs =  reader.getFileHeader().getSequenceDictionary().getSequences();
		for (SAMSequenceRecord sr: seqs) chromLength.put(sr.getSequenceName(), sr.getSequenceLength()+1000);

		SAMRecordIterator iterator = reader.iterator();

		HashSet<String> priorChroms = new HashSet<String>();
		String chrom = null;

		//make first chromosome
		SAMRecord sam = null;
		int chrStartBp = -1;
		ArrayList<Integer>[] bpNames = null;
		boolean[] badBases = null;
		while (iterator.hasNext()){
			sam = iterator.next();

			//unaligned? 
			if (sam.getReadUnmappedFlag()) continue;
			chrom = sam.getReferenceName();

			if (chromGenes.containsKey(chrom) == false) {
				continue;
			}		
			priorChroms.add(chrom);
			chrStartBp = sam.getUnclippedStart()-1;
			badBases = new boolean[chromLength.get(chrom) - chrStartBp];
			bpNames = new ArrayList[badBases.length];
			break;
		}
		//any data found?
		if (bpNames == null) {
			reader.close();
			return;
		} 

		//reset working fields
		workingFragNameIndex.clear();
		workingFragmentNameIndexPlus = 1;
		workingFragmentNameIndexMinus = -1;

		//load first alignment into bpNames
		loadBlocks (sam, chrStartBp, bpNames, badBases);

		//for each record
		while (iterator.hasNext()){
			sam = iterator.next();

			//unaligned? 
			if (sam.getReadUnmappedFlag()) continue;

			//same chrom?
			if (sam.getReferenceName().equals(chrom)){
				loadBlocks (sam, chrStartBp, bpNames, badBases);
			}
			else {
				//different chrom so time to scan
				loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);

				//check that new chrom from SAM is something interrogated by their gene list
				boolean reset = false;
				if (chromGenes.containsKey(sam.getReferenceName()) == false) {
					//advance until it does
					while (iterator.hasNext()){
						sam = iterator.next();
						if (chromGenes.containsKey(sam.getReferenceName())) {
							reset = true;
							break;
						}
					}
				}
				else reset = true;

				//reset
				if (reset){
					//reset working fields
					workingFragNameIndex.clear();
					workingFragmentNameIndexPlus = 1;
					workingFragmentNameIndexMinus = -1;

					chrom = sam.getReferenceName();
					if (priorChroms.contains(chrom)) Misc.printErrAndExit("\nError: your sam file isn't sorted by chromosome! Aborting.\n");
					priorChroms.add(chrom);
					chrStartBp = sam.getUnclippedStart()-1;
					badBases = new boolean[chromLength.get(chrom) - chrStartBp];
					bpNames = new ArrayList[badBases.length];

					//load
					loadBlocks (sam, chrStartBp, bpNames, badBases);
				}

			}

		}

		if (chromGenes.containsKey(chrom)) loadGeneCounts(replica, bpNames, chrStartBp, chrom, badBases);


		reader.close();
		bpNames = null;
		iterator = null;
		seqs = null;
		chromLength = null;
	}

	private void loadGeneCounts(Replica replica, ArrayList<Integer>[] bpNames, int chrStartBp, String chromosome, boolean[] badBases){
		HashMap<String, GeneCount> geneCounts = replica.getGeneCounts();
		HashSet<Integer> allReads = new HashSet<Integer>();
		HashSet<Integer> exonReads = new HashSet<Integer>();
		int lengthBpNames = bpNames.length -1;

		//for each gene in the chromosome 
		UCSCGeneLine[] chrGenes = chromGenes.get(chromosome);

		for (int i=0; i< chrGenes.length; i++){
			//flagged gene
			if (chrGenes[i].isFlagged()) continue;

			String geneName = chrGenes[i].getDisplayNameThenName();					
			allReads.clear();

			//get exons
			ExonIntron[] exons = chrGenes[i].getExons();
			int[] exonCounts = new int[exons.length];

			//for each exon
			exonLoop:
				for (int x=0; x< exons.length; x++){
					exonReads.clear();
					int start = exons[x].getStart() - chrStartBp;
					if (start < 0) start = 0;
					int end = exons[x].getEnd() - chrStartBp;
					if (end > lengthBpNames) end = lengthBpNames;
					//for each base in the exon, see if there is a read
					for (int y=start; y< end; y++){
						//bad base? if so then flag entire gene
						if (badBases[y]) {
							chrGenes[i].setFlagged(true);
							flaggedGeneNames.add(geneName);
							allReads.clear();
							break exonLoop;
						}
						if (bpNames[y] != null) exonReads.addAll(bpNames[y]);
					}
					exonCounts[x] = exonReads.size();
					allReads.addAll(exonReads);
				}
			int numCounts = allReads.size();
			if (numCounts !=0){
				GeneCount tcg = new GeneCount(numCounts, exonCounts);
				geneCounts.put(geneName, tcg);
				replica.setTotalCounts(replica.getTotalCounts() + numCounts);
				if (numCounts >= minimumCounts) geneNamesWithMinimumCounts.add(geneName);
			}
		}	
		//clean up
		allReads = null;
		exonReads = null;
	}

	/**Just for cases*/
	public void buildCaseGeneExonSamples(){
		try {
			//start counter for all exons to match Alun's code, starts with 1
			int globalExonIndex = 1;
			ArrayList<GeneExonSample> gesAL = new ArrayList<GeneExonSample>();

			//for genes with minimal counts 
			for (String geneName: geneNamesToAnalyze){
				UCSCGeneLine gene = name2Gene.get(geneName);
				int numExons = gene.getExons().length;

				//collect exon counts int[sample][exon]
				int[][] caseExonCounts = fetchExonCounts(geneName, numExons, cases);

				//for each exon
				for (int i=0; i< numExons; i++){
					//count total for this exon across all the cases, skip if too few
					int total = 0;
					for (int j=0; j<caseExonCounts.length; j++) total += caseExonCounts[j][i];
					if (total < minimumCounts) continue;
					numExonsProcessed++;
					
					//for each sample
					for (int j=0; j<caseExonCounts.length; j++){
						GeneExonSample g = new GeneExonSample(gene, globalExonIndex-1, (short)i, (short)j, caseExonCounts[j][i]);
						gesAL.add(g);
					}
					globalExonIndex++;
				}
			}
			//save GES
			ges = new GeneExonSample[gesAL.size()];
			gesAL.toArray(ges);
		}
		catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing out count table.");
		}
	}
	
	public void chunkGeneExonSamples(){
		//group by global exon index
		ArrayList<GeneExonSample>[] byGlobal = new ArrayList[numExonsProcessed];
		//for ()
		for (int i=0; i< ges.length; i++){
			
		}
	}


	/**Just for cases*/
	public void writeCaseCountTable(){
		countTableFile = new File(saveDirectory, "countTable.txt");

		try {
			PrintWriter out = new PrintWriter( new FileWriter(countTableFile));

			//start counter for all exons to match Alun's code, starts with 1
			int globalExonIndex = 1;
			ArrayList<GeneExonSample> gesAL = new ArrayList<GeneExonSample>();

			//for genes with minimal counts 
			for (String geneName: geneNamesToAnalyze){
				UCSCGeneLine gene = name2Gene.get(geneName);
				int numExons = gene.getExons().length;

				//collect exon counts int[sample][exon]
				int[][] caseExonCounts = fetchExonCounts(geneName, numExons, cases);

				//for each exon
				for (int i=0; i< numExons; i++){

					//count total for this exon across all the cases, skip if too few
					int total = 0;
					for (int j=0; j<caseExonCounts.length; j++) total += caseExonCounts[j][i];
					if (total < minimumCounts) continue;
					numExonsProcessed++;
					
					//for each sample
					StringBuilder sb = new StringBuilder(); 
					for (int j=0; j<caseExonCounts.length; j++){
						//add exon index, 1 based, global
						sb.append(globalExonIndex);  
						sb.append("\t"); 
						//add sample index, 1 based
						sb.append(j+1);
						sb.append("\t");
						//add counts
						sb.append(caseExonCounts[j][i]);
						sb.append("\n");
						//save GES
						GeneExonSample g = new GeneExonSample(gene, globalExonIndex-1, (short)i, (short)j, caseExonCounts[j][i]);
						gesAL.add(g);
					}
					out.print(sb);
					globalExonIndex++;
				}
			}
			out.close();
			//save GES
			ges = new GeneExonSample[gesAL.size()];
			gesAL.toArray(ges);
		}
		catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing out count table.");
		}
	}

	
	public void writeCountTableMatrix(){
		countTableFile = new File(saveDirectory, "countTable.txt");

		try {
			//write matrix of name, t,t,t...c,c,c,c...d,d,d, to file for genes with observations
			PrintWriter out = new PrintWriter( new FileWriter(countTableFile));

			//print header
			out.print("#GeneName\tExonIndex\tCoordinates");
			//for each sample
			for (int i=0; i<caseSampleFileNames.length; i++){
				out.print("\t");
				out.print(Misc.removeExtension(caseSampleFileNames[i]));
				if (controls != null){
					out.print("\t");
					out.print(Misc.removeExtension(controlSampleFileNames[i]));
				}
			}
			out.println();

			//print counts for genes with minimal counts in any replica
			for (String geneName: geneNamesToAnalyze){
				UCSCGeneLine gene = name2Gene.get(geneName);
				ExonIntron[] exons = gene.getExons();
				int numExons = gene.getExons().length;
				String chr= gene.getChrom();

				//collect exon counts int[sample][exon
				int[][] caseExonCounts = fetchExonCounts(geneName, numExons, cases);
				int[][] controlExonCounts = null;
				if (controls != null) controlExonCounts = fetchExonCounts(geneName, numExons, controls);

				//for each exon
				for (int i=0; i< numExons; i++){
					//print geneName and exon index
					StringBuilder sb = new StringBuilder(geneName+"\t"+i+"\t"+chr+":"+exons[i].getStartStopString());
					int total = 0;
					//for each sample
					for (int j=0; j<caseExonCounts.length; j++){
						sb.append("\t");
						sb.append(caseExonCounts[j][i]);
						total += caseExonCounts[j][i];
						if (controls != null){
							sb.append("\t");
							sb.append(controlExonCounts[j][i]);
							total += controlExonCounts[j][i];
						}
					}
					if (total >= minimumCounts) out.println(sb);
				}
			}
			out.close();
		}
		catch (Exception e){
			System.err.println("\nProblem writing out count table.");
			e.printStackTrace();
		}
	}

	private int[][] fetchExonCounts(String geneName, int numExons, Replica[] samples) {
		int[][] exonCounts = new int[samples.length][];
		int[] zero = new int[numExons];
		for (int i=0; i< samples.length; i++){
			GeneCount gc = samples[i].getGeneCounts().get(geneName);
			if (gc == null) exonCounts[i] = zero;
			else exonCounts[i] = gc.getExonCounts();
		}
		return exonCounts;
	}


	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(refSeqFile, 0);
		genes = reader.getGeneLines();
		if (genes == null || genes.length == 0) Misc.printExit("\nProblem loading your USCS gene model table or bed file? No genes/ regions?\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check gene name is unique
		if (reader.uniqueGeneNames() == false) Misc.printExit("\nDuplicate gene names were found in your gene / bed file, these must be unique.\n");
		chromGenes = reader.getChromSpecificGeneLines();
		name2Gene = new HashMap<String, UCSCGeneLine>();
		for (UCSCGeneLine gl: genes) {
			String name = gl.getDisplayNameThenName();
			name2Gene.put(name, gl);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PoisRegCNV(args);
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
					case 's': saveDirectory = new File(args[++i]); break;
					case 'b': bamDirectory = new File(args[++i]); break;
					case 'c': caseSampleFileNames = Misc.COMMA.split(args[++i]); break;
					case 'n': controlSampleFileNames = Misc.COMMA.split(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'a': alunRScript = new File(args[++i]); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'g': genomeVersion = args[++i]; break;
					case 'e': minimumCounts = Integer.parseInt(args[++i]); break;
					case 'x': maxAlignmentsDepth = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}

		//check genomeVersion
		if (genomeVersion == null) Misc.printErrAndExit("\nError: please provide a genome version, (e.g. H_sapiens_Mar_2006)\n");

		//look for bam files
		if (bamDirectory == null) Misc.printErrAndExit("\nError: cannot find your directory of bam files?\n");
		File[] bamFiles = IO.extractFiles(bamDirectory, ".bam");
		if (bamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any bam files?\n");
		OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, false);

		//check for case names
		if (caseSampleFileNames == null ) {
			File[] bams = IO.fetchFilesRecursively(bamDirectory, ".bam");
			if (bams == null || bams.length == 0) Misc.printErrAndExit("\nError: cannot find any xxx.bam files to process?!\n");
			caseSampleFileNames = new String[bams.length];
			for (int i=0; i< bams.length; i++) caseSampleFileNames[i] = bams[i].getName();
		}
		else {
			for (int i=0; i< caseSampleFileNames.length; i++){
				File bam = new File(bamDirectory, caseSampleFileNames[i]);
				if (bam.exists() == false) Misc.printErrAndExit("Error: cannot find the "+caseSampleFileNames[i]+" in the bam directory?!\n");
			}
		}

		//check controls
		if (controlSampleFileNames != null){
			//same number?
			if (controlSampleFileNames.length != caseSampleFileNames.length) Misc.printErrAndExit("\nError: the number of case and control sample file names differ?\n");
			for (int i=0; i< controlSampleFileNames.length; i++){
				File bam = new File(bamDirectory, controlSampleFileNames[i]);
				if (bam.exists() == false) Misc.printErrAndExit("Error: cannot find the "+controlSampleFileNames[i]+" in the bam directory?!\n");
			}
		}

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdirs();

		//check for R and required libraries, don't need it if they just want the first and last 1/3 count table
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}

		if (alunRScript == null || alunRScript.canRead()== false) {
			Misc.printErrAndExit("\nError: Cannot find or read Alun's RScript? -> "+alunRScript+"\n");
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Pois Reg CNV:   Feb 2015                            **\n" +
				"**************************************************************************************\n" +
				"Uses poisson regression and Pearson residuals to identify exons whose counts differ\n"+
				"significantly from the fitted value base on all the exon sample counts. This app wraps\n"+
				"an algorithm developed by Alun Thomas.  Data tracks are generated for the residuals\n"+
				"and log2(observed/ fitted counts) as well as detailed spreadsheets.\n"+

				"\nRequired Options:\n"+
				"-s Save directory.\n"+
				"-b BAM file directory, sorted and indexed by coordinate.\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. Tab delimited, see RefSeq Genes\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 . NOTE:\n" +
				"       this table should contain only ONE composite transcript per gene (e.g. use\n" +
				"       Ensembl genes NOT transcripts). Use the MergeUCSCGeneTable app to collapse\n" +
				"       transcripts. See http://useq.sourceforge.net/usageRNASeq.html for details.\n"+
				"-a Alun Thomas R script file.\n"+
				"-g Genome Version  (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +

				"\nDefault Options:\n"+
				"-c BAM file names for cases to process, comma delimited, no spaces, defaults to all.\n"+
				"-r Full path to R, defaults to /usr/bin/R\n"+
				"-e Minimum all sample exon count for inclusion in analysis, defaults to 20.\n"+
				"-x Max per sample exon alignment depth, defaults to 50000. Exons containing higher\n"+
				"       densities are ignored.\n"+
				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/PoisRegCNV -s PRCnvResults/ \n" +
				"      -b BamFiles -u hg19EnsGenes.ucsc.gz -a RBambedSource.R  -g H_sapiens_Feb_2009\n\n" +

				"**************************************************************************************\n");

	}
}
