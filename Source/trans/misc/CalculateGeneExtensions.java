package trans.misc;
import java.util.*;

import util.bio.wrappers.*;
import trans.roc.*;
import java.io.*;
import java.util.regex.*;
import util.bio.parsers.*;
import util.gen.*;

/**Finds isolated genes that are expressed and using Huber's segmentation algorithm, slices it into three and prints the extension.
 * Note, this is not stranded.  If you want to keep it stranded then parse the UCSC Gene Table and use appropriate xxx.bar files.*/
public class CalculateGeneExtensions {
	//fields
	private File ucscGeneTableFileAll;
	private File ucscGeneTableFileSelect;
	private String selectName;
	private File barDirectory;
	private int extension = 250; //1000 bp off either stop to look for adjacent genes
	private int extensionToSegment = 500; // 750 bp off either stop to attempt a segmentation
	private int bpShrink = 0; //shrink size of either stop of gene to enable shorter transcription
	private float threshold = 3; //log2
	private int bpPositionOffSetBar = -30;
	private int bpPositionOffSetRegion = -60;
	private int oligoLengthMinusOne = 59;
	private String version = "S_pombe_Apr_2007";
	private HashMap genesAll;
	private HashMap genesSelect;
	private File rApp;
	private File tempDataFile;
	private SegDat[] segDats;


	public CalculateGeneExtensions( String[] args){
		processArgs(args);

		//parse gene table, split by chromosome
		System.out.println("\nParsing gene tables...");
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader (ucscGeneTableFileAll, 0);
		reader.splitByChromosome();
		genesAll = reader.getChromSpecificGeneLines();
		reader = new UCSCGeneModelTableReader (ucscGeneTableFileSelect, 0);
		reader.splitByChromosome();
		genesSelect = reader.getChromSpecificGeneLines();

		//find isolated genes
		System.out.println("Finding isolated genes, +/- "+extension);
		findIsolatedGenes();
		
		//print isolated
		Iterator it = genesSelect.keySet().iterator();
		while (it.hasNext()){
			//get genes
			String chrom = (String)it.next();
			UCSCGeneLine[] chrGenes = (UCSCGeneLine[])genesSelect.get(chrom);
			//for each gene 
			for (int i=0; i< chrGenes.length; i++){
				System.out.println(chrGenes[i]);
			}
		}		
System.exit(0);

		//find expressed genes
		System.out.println("Finding expressed genes, log2 median threshold "+threshold);
		findExpressedGenes();

		//segment expressed genes
		System.out.println("Segmenting expressed genes into 3 parts");
		segmentGenes();
		System.out.println("\tNumber segmented genes, "+segDats.length);

		//flag for those where the middle segment is < lenght of gene
		flagSegmentedGenes();

		//estimate extensions
		estimateExtensions();
	}

	public void estimateExtensions(){
		try {
			System.out.println("Bp extension estimations per gene:<p><p>\n\tName\tChr\tOri\tStart\tStop\tGeneLength\tLog2(MedianExprLvl)\t5'\t3'\tURL<p>");
			File sgr = new File (ucscGeneTableFileSelect.getParentFile(),selectName+"_extensions.sgr");
			PrintWriter out = new PrintWriter (new FileWriter (sgr));
			ArrayList fivePrimeEx = new ArrayList();
			ArrayList threePrimeEx = new ArrayList();
			for (int i=0; i< segDats.length; i++){
				//print sgr
				String chrom = segDats[i].geneLine.getChrom();
				String five = null;
				String three = null;
				if (segDats[i].goodSegmentation){
					//print gene boundaries
					out.println(chrom+"\t"+segDats[i].geneLine.getTxStart()+"\t"+10);
					out.println(chrom+"\t"+segDats[i].geneLine.getTxEnd()+"\t"+10);
					//print extensions
					out.println(chrom+"\t"+segDats[i].getFirstSegBpPosConservative()+"\t"+5);
					out.println(chrom+"\t"+segDats[i].getLastSegBpPosConservative()+"\t"+5);
					//calc extensions
					int left = segDats[i].geneLine.getTxStart() - segDats[i].getFirstSegBpPosConservative();
					int right = segDats[i].getLastSegBpPosConservative() - segDats[i].geneLine.getTxEnd();
					//orientation?
					if (segDats[i].geneLine.getStrand().startsWith("+")){
						fivePrimeEx.add(new Float(left));
						threePrimeEx.add(new Float(right));
						five = left+"";
						three = right+"";
					}
					else {
						fivePrimeEx.add(new Float(right));
						threePrimeEx.add(new Float(left));
						five = right+"";
						three = left+"";
					}
				}
				//print summary line
				String name = segDats[i].geneLine.getName();
				int length = segDats[i].geneLine.getTxEnd()-segDats[i].geneLine.getTxStart() +1;
				int subEnds = (10000-length)/2;
				float median = segDats[i].geneLine.getScores()[0];
				String url = "<A HREF=\"http://localhost:7085/UnibrowControl?version="+version+"&seqid="+chrom+"&start="+
				(segDats[i].geneLine.getTxStart()-subEnds)+"&stop="+(segDats[i].geneLine.getTxEnd()+subEnds)+"\">"+name+"</A><BR>";
				System.out.println("\t"+name+"\t"+chrom+"\t"+segDats[i].geneLine.getStrand()+"\t"+segDats[i].geneLine.getTxStart()+"\t"+segDats[i].geneLine.getTxEnd()+"\t"+length+"\t"+median+"\t"+five+"\t"+three+"\t"+url);
			}
			
			//calculate summary stats on each stop
			float[] ext = Num.arrayListOfFloatToArray(fivePrimeEx);
			System.out.println("\nSummary stats on 5' stop extensions");
			Num.statFloatArray(ext, false);
			
			ext = Num.arrayListOfFloatToArray(threePrimeEx);
			System.out.println("\nSummary stats on 3' stop extensions");
			Num.statFloatArray(ext, false);
			
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void flagSegmentedGenes(){
		System.out.print("\tFlagging:");
		int num = 0;
		for (int i=0; i< segDats.length; i++){
			if (segDats[i].getFirstSegBpPosConservative()<= segDats[i].geneLine.getTxStart() + oligoLengthMinusOne && 
					segDats[i].getLastSegBpPosConservative() >= segDats[i].geneLine.getTxEnd() - oligoLengthMinusOne) {
				segDats[i].goodSegmentation = true;
			}
			else {
				System.out.print(" "+ segDats[i].geneLine.getName());
				num++;
			}
		}
		System.out.println("\n\t"+num+" Flagged");
		
	}

	public void segmentGenes(){
		//make SegDats
		makeSegDats();
		//write table of intensities
		writeTableOfIntensities();
		//segment
		segment();
		//cleanup
		tempDataFile.delete();
	}

	public void segment(){
		Segment seg = new Segment (rApp, tempDataFile, tempDataFile.getParentFile());
		try {
			//segment data file
			seg.segment();
			//instantiate writer for slices
			File sgr = new File (ucscGeneTableFileSelect.getParentFile(),selectName+"_segments.sgr");
			PrintWriter out = new PrintWriter (new FileWriter (sgr));
			//assign segment indexes to each SegDat
			int[][] segs = seg.getSegmentIndexes();
			for (int i=0; i< segDats.length; i++) {
				segDats[i].segmentIndexes = segs[i];
				out.println(segDats[i].geneLine.getChrom()+"\t"+segDats[i].getFirstSegBPPosition()+"\t10");
				out.println(segDats[i].geneLine.getChrom()+"\t"+segDats[i].getSecondSegBPPosition()+"\t10");
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void writeTableOfIntensities(){
		//make tempDataFile
		String randomWord = Passwords.createRandowWord(8);
		tempDataFile = new File (ucscGeneTableFileSelect.getParentFile(), "tempFile_SegmentColumns_"+randomWord+".xls");
		try {
			PrintWriter out = new PrintWriter (new FileWriter(tempDataFile));
			//print number of values
			StringBuffer sb = new StringBuffer();
			int maxLength = segDats[0].dataSlice.length;
			sb.append(maxLength);
			for (int i=1; i< segDats.length; i++){
				sb.append("\t");
				int length = segDats[i].dataSlice.length;
				if (length > maxLength) maxLength = length;
				sb.append(length);
			}
			out.println(sb);
			//for each row
			for (int i=0; i< maxLength; i++){
				sb = new StringBuffer();
				//for each gene/ column
				//first gene
				if (segDats[0].dataSlice.length > i) sb.append (segDats[0].dataSlice[i].getScore());
				//all remaining genes
				for (int j=1; j< segDats.length; j++){
					sb.append("\t");
					if (segDats[j].dataSlice.length > i) sb.append (segDats[j].dataSlice[i].getScore());
				}
				out.println(sb);
			}
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}	
	}

	/**Segments genes: applies median value across entire exon and intron structure,
	 * collects intensities surrounding gene,
	 * runs Hubers segmentation script on parsed gene.*/
	public void makeSegDats(){
		try {
			File preSgr = new File (ucscGeneTableFileSelect.getParentFile(),selectName+"_exprGenes.sgr");
			PrintWriter outPreSgr = new PrintWriter (new FileWriter (preSgr));
			File postSgr = new File (ucscGeneTableFileSelect.getParentFile(),selectName+"_exprGeneMedian.sgr");
			PrintWriter outPostSgr = new PrintWriter (new FileWriter (postSgr));
			//instantiate Bar parser
			ScoreParsedBars barParser = new ScoreParsedBars (barDirectory, bpPositionOffSetBar, bpPositionOffSetRegion);
			//make array to hold SegDats
			ArrayList segDatsAL = new ArrayList();
			//for each chromosome in select genes
			Iterator it = genesSelect.keySet().iterator();
			while (it.hasNext()){
				//get genes
				String chrom = (String)it.next();
				UCSCGeneLine[] chrGenes = (UCSCGeneLine[])genesSelect.get(chrom);
				//for each gene 
				for (int i=0; i< chrGenes.length; i++){
					//fetch Grs
					ArrayList al = barParser.fetchRegionGrs(chrom, chrGenes[i].getTxStart() - extensionToSegment, chrGenes[i].getTxEnd() + extensionToSegment);
					Gr[] grs = new Gr[al.size()];
					al.toArray(grs);
					//for each Gr, check to see if it is entirely contained in the gene, if so assign the median value
					int posStart = -1;
					int posEnd = -1;
					float median = chrGenes[i].getScores()[0];
					//these gr positions have been modified by the bpPositionOffSet, shifted left to start of oligo
					//shrink size too
					int start = chrGenes[i].getTxStart()+ bpShrink;
					int end = chrGenes[i].getTxEnd() + bpPositionOffSetRegion - bpShrink;
					if (end < start) {
						System.out.println("\tSkipping "+chrGenes[i].getName() + " End shrinking eliminated gene.");
					}
					for (int x =0; x< grs.length; x++){
						outPreSgr.println(chrom +"\t"+grs[x]);
						//clip it to median
						if (grs[x].getScore() > median) grs[x].setScore(median);
						//is it contained? Set to median
						if (grs[x].getPosition() >= start && grs[x].getPosition() <= end) {
							//set score to median
							grs[x].setScore(median);						
							//set start
							if (posStart == -1) posStart = x;
						}
						//not contained, has start been set? if so then check stop
						else if (posStart != -1 && posEnd == -1){
							//start has been set, stop has not, thus set stop to prior Gr position
							posEnd = x - 1;
						}
						outPostSgr.println(chrom +"\t"+grs[x]);
					}

					//save data
					segDatsAL.add(new SegDat(chrGenes[i], grs, posStart, posEnd));
				}
			}
			segDats = new SegDat[segDatsAL.size()];
			segDatsAL.toArray(segDats);

			outPreSgr.close();
			outPostSgr.close();
		} catch (Exception e){
			e.printStackTrace();
		}

	}

	/**Container for holding data related to a gene.*/
	private class SegDat{
		//fields
		UCSCGeneLine geneLine;
		Gr[] dataSlice;				//these gr positions have been modified by the bpPositionOffSet, shifted left to start of oligo
		int startGeneOligoIndex;
		int stopGeneOligoIndex;
		int[] segmentIndexes;
		boolean goodSegmentation = false;

		//constructor
		public SegDat(UCSCGeneLine geneLine, Gr[] dataSlice, int startGeneOligoIndex, int stopGeneOligoIndex ){
			this.geneLine = geneLine;
			this.dataSlice = dataSlice;
			this.startGeneOligoIndex = startGeneOligoIndex;
			this.stopGeneOligoIndex = stopGeneOligoIndex;
		}

		public int getFirstSegBpPosConservative (){
			return dataSlice[segmentIndexes[0]].getPosition();
		}

		public int getLastSegBpPosConservative () {
			return dataSlice[segmentIndexes[1]-1].getPosition() + oligoLengthMinusOne;
		}

		public int getFirstSegBPPosition(){
			int index = segmentIndexes[0];
			return halfDistance (dataSlice[index -1], dataSlice[index]);
		}

		public int getSecondSegBPPosition(){
			int index = segmentIndexes[1];
			return halfDistance (dataSlice[index -1], dataSlice[index]);
		}

		public int halfDistance (Gr first, Gr second){
			float size = second.getPosition() - first.getPosition();
			size = size/2.0f;
			return Math.round(size) + first.getPosition();
		}
	}

	public void segmentGrArray (Gr[] grs){
		//extract scores
		float[] scores = new float[grs.length];
		for (int i=0; i< scores.length; i++) scores[i] = grs[i].getScore();
		//write scores to file
		Num.writeToFile(scores, tempDataFile);

	}

	/**Find genes with no neighbors, replaces the UCSCGeneLine[] array in the genes HashMap.*/
	public void findExpressedGenes(){
		HashMap exp = new HashMap();
		//instantiate Bar parser
		ScoreParsedBars barParser = new ScoreParsedBars (barDirectory, bpPositionOffSetBar, bpPositionOffSetRegion);
		//for each chromosome in select
		Iterator it = genesSelect.keySet().iterator();
		while (it.hasNext()){
			//get genes
			String chrom = (String)it.next();
			UCSCGeneLine[] chrGenes = (UCSCGeneLine[])genesSelect.get(chrom);
			System.out.println("\tTotal "+chrom+"\t"+ chrGenes.length);
			//calculate medians
			barParser.calculateMedians(chrGenes);
			ArrayList expressedGenes = new ArrayList();
			//for each gene look for intersection after expanding
			for (int i=0; i< chrGenes.length; i++){
				if (chrGenes[i].getScores()[0] >= threshold ) {
					expressedGenes.add(chrGenes[i]);
					//System.out.println(chrGenes[i].simpleToString()+ chrGenes[i].getScore());
				}
			}
			//convert ArrayList to array and replace in hash
			chrGenes = new UCSCGeneLine[expressedGenes.size()];
			expressedGenes.toArray(chrGenes);
			System.out.println("\tExpressed "+chrom+"\t"+ chrGenes.length);
			exp.put(chrom, chrGenes);
		}
		genesSelect = exp;
	}


	/**Find genes with no neighbors, replaces the UCSCGeneLine[] array in the genes HashMap.*/
	public void findIsolatedGenes(){
		HashMap isolated = new HashMap();
		//for each chromosome
		Iterator it = genesSelect.keySet().iterator();
		while (it.hasNext()){
			String chrom = (String)it.next();
			UCSCGeneLine[] chrGenesAll = (UCSCGeneLine[])genesAll.get(chrom);
			UCSCGeneLine[] chrGenesSelect = (UCSCGeneLine[])genesSelect.get(chrom);
			System.out.println("\tTotal "+chrom+"\t"+ chrGenesSelect.length);
			ArrayList isolatedGenes = new ArrayList();
			//for each gene look for intersection after expanding
			for (int i=0; i< chrGenesSelect.length; i++){
				//expand size of gene
				int start = chrGenesSelect[i].getTxStart() - extension;
				int end = chrGenesSelect[i].getTxEnd() + extension;
				//scan all, allow one hit
				int hits = 0;
				for (int j =0; j< chrGenesAll.length; j++){
					if (intersects (chrGenesAll[j], start, end)) {
						hits++; 
						if (hits >1) break;
					}
				}
				//save?
				if (hits < 2){
					isolatedGenes.add(chrGenesSelect[i]);
				}
			}
			//convert ArrayList to array and replace in hash
			chrGenesSelect = new UCSCGeneLine[isolatedGenes.size()];
			isolatedGenes.toArray(chrGenesSelect);
			System.out.println("\tIsolated "+chrom+"\t"+ chrGenesSelect.length);
			isolated.put(chrom, chrGenesSelect);
		}
		genesSelect = isolated;
		//clean up
		genesAll = null;
		System.gc();
	}

	/**Assumes coordinates are inclusive.*/
	public boolean intersects (UCSCGeneLine first, int start, int end){
		if (end < first.getTxStart() || start > first.getTxEnd()) return false;
		return true;
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': ucscGeneTableFileAll = new File(args[i+1]);  i++; break;
					case 'g': ucscGeneTableFileSelect = new File(args[i+1]);  i++; break;
					case 'b': barDirectory = new File(args[i+1]);  i++; break;
					case 'r': rApp = new File(args[i+1]);  i++; break;
					case 's': threshold = Float.parseFloat(args[i+1]); i++; break;
					case 'e': extension = Integer.parseInt(args[i+1]); i++; break;
					case 'f': extensionToSegment = Integer.parseInt(args[i+1]); i++; break;
					case 'x': bpPositionOffSetBar = Integer.parseInt(args[i+1]); i++; break;
					case 'y': bpPositionOffSetRegion = Integer.parseInt(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//parse text
		selectName = Misc.removeExtension(ucscGeneTableFileSelect.getName());

	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Convert Nimblegen NDF 2 TPMap: Jan 2007                  **\n" +
				"**************************************************************************************\n" +
				"Converts a Nimblegen NDF text file to a tpmap.\n\n" +

				"-f Full path file text or directory containing text xxx.ndf file(s).\n\n" +

				"Example: java -Xmx1500M -jar pathTo/T2/Apps/ConvertNimblegenNDF2TPMap -f /badData/\n\n" +

		"**************************************************************************************\n");		
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new CalculateGeneExtensions(args);
	}
}
