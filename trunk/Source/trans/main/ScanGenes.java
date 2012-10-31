package trans.main;
import java.io.*;

import util.bio.parsers.*;
import trans.misc.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;
import util.bio.annotation.*;

/**
 * Given a list of genes and their exons extracts the associated intensities. 
 */
public class ScanGenes {

	//fields
	private File resultsFile;
	private File[] treatmentDirectories = null;
	private File[] controlDirectories = null;
	private File[] chromosomeOligoPositions = null;
	private float[][] treatmentIntensities;
	private float[][] controlIntensities;
	private String chromosome = null;  
	private int[] positions;
	private HashMap geneModels = null;
	private File geneFile = null;
	private int endOffSet = 23;
	private int numToSubFromEnds = 1;
	private ArrayList scoredGenes = new ArrayList();
	private boolean printIntensities = false;

	public ScanGenes(String[] args){
		//process args
		processArgs(args);

		//load gene models
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(geneFile, numToSubFromEnds);
		reader.splitByChromosome();
		geneModels = reader.getChromSpecificGeneLines();

		//for each chromosome	
		System.out.println();
		for (int i=0; i<chromosomeOligoPositions.length; i++){

			chromosome = chromosomeOligoPositions[i].getName();
			System.out.println("Testing chromosome: "+chromosome); 

			//get positions
			positions = (int[])IO.fetchObject(chromosomeOligoPositions[i]);

			//get genes for a particular chromosome, these are sorted by position, some overlap and may share same start
			UCSCGeneLine[] genes = (UCSCGeneLine[])geneModels.get(chromosome);
			if (genes == null) continue;

			//get treatment intensities float[replicas][values]
			treatmentIntensities = ScanChromosomes.fetchIntensities (treatmentDirectories, chromosome, false);
			if (treatmentIntensities == null) Misc.printExit("\nProblem with fetching treatment intensities. Any missing or extra files?\n");
			if (ScanChromosomes.checkSizes(treatmentIntensities) == false) Misc.printExit("\nYour treatment intensity arrays are of different lengths?!\n");

			//get control intensities float[replicas][values]
			controlIntensities = ScanChromosomes.fetchIntensities (controlDirectories, chromosome, false);
			if (controlIntensities == null) Misc.printExit("\nProblem with fetching control intensities. Any missing or extra files?\n");
			if (ScanChromosomes.checkSizes(controlIntensities) == false) Misc.printExit("\nYour control intensity arrays are of different lengths?!\n");
			if (treatmentIntensities[0].length != controlIntensities[0].length) Misc.printExit("\nYour treatment and control intensity arrays are of different lengths?!\n");


			//for each gene
			for (int j=0; j<genes.length; j++){

				//get Exons
				ExonIntron[] exons = genes[j].getExons();
				if (exons == null){
					System.out.println("No exons skipping\t"+genes[j]);
					scoredGenes.add(new ScoredGene(genes[j]));
					continue;
				}
				//convert to start stop indexes, some exons may be too small or not covered
				int[][] startStops = convertToIndexPositions(exons);
				if (startStops == null) {
					scoredGenes.add(new ScoredGene(genes[j]));
					continue;	//any found?
				}

				//for each startStop collect associates scores
				//calculate size
				int numberOligos = Util.countNumberOfOligoPositions(startStops);
				if (numberOligos <=0 ){
					System.out.println("No probes skipping\t"+genes[j]);
					scoredGenes.add(new ScoredGene(genes[j]));
					continue;
				}

				float[][] subT = new float[treatmentIntensities.length][numberOligos];
				float[][] subC = new float[controlIntensities.length][numberOligos];
				int counter = 0;
				for (int x=0; x< startStops.length; x++){
					int[] startStop = startStops[x];
					//for each oligo position
					for (int k=startStop[0]; k<= startStop[1]; k++){
						//collect treatment replicas
						for (int l=0; l< treatmentIntensities.length; l++) {
							subT[l][counter] = treatmentIntensities[l][k];
						}
						//collect control replicas
						for (int l=0; l< controlIntensities.length; l++) {
							subC[l][counter] = controlIntensities[l][k];
						}
						counter++;
					}
				}
				ScoredGene sg = new ScoredGene (genes[j], subT, subC);
				float[] ratios = Num.ratio(subT, subC);
				Arrays.sort(ratios);
				double median = Math.log(Num.median(ratios))/Num.log2;
				sg.setMedian(median);
				double mean = Num.mean(ratios);
				double cv = Num.standardDeviation(ratios, mean) / mean;
				sg.setMisc(cv);
				scoredGenes.add(sg);
			}

		}

		//print results
		ScoredGene[] sgs = new ScoredGene[scoredGenes.size()];
		scoredGenes.toArray(sgs);
		Arrays.sort(sgs);
		try{
			PrintWriter out = new PrintWriter ( new FileWriter (resultsFile));
			out.print("Name\tChromosome\tStrand\tTxStart\tTxStop\tLog2(Median)\tCoefVar\tTreatmentMedians\tControlMedians");
			if (printIntensities) out.print("\tTreatments\tControls");
			boolean notesPresent = false;
			if (sgs[0].getGene().getNotes()!=null) {
				out.print("\tNotes");
				notesPresent = true;
			}
			out.println();
			for (int i=0; i< sgs.length; i++){
				UCSCGeneLine gene = sgs[i].getGene();
				StringBuffer sb = new StringBuffer(gene.getName());
				sb.append("\t");
				sb.append(gene.getChrom());
				sb.append("\t");
				sb.append(gene.getStrand());
				sb.append("\t");
				sb.append(gene.getTxStart());
				sb.append("\t");
				sb.append(gene.getTxEnd());
				sb.append("\t");
				if (sgs[i].getTreatmentValues() != null){
					sb.append(sgs[i].getMedian());
					sb.append("\t");
					sb.append(sgs[i].getMisc());
					sb.append("\t");
					sb.append(Num.medianFloatArray(sgs[i].getTreatmentValues()));
					sb.append("\t");
					sb.append(Num.medianFloatArray(sgs[i].getControlValues()));
					if (printIntensities) {
						float[] treat = Num.collapseFloatArray(sgs[i].getTreatmentValues());
						float[] control = Num.collapseFloatArray(sgs[i].getControlValues());
						String t = Misc.floatArrayToString(treat, ",");
						String c = Misc.floatArrayToString(control, ",");
						sb.append("\t");
						sb.append(t);
						sb.append("\t");
						sb.append(c);
					}
					if (notesPresent){
						sb.append("\t");
						sb.append(gene.getNotes());
					}
				}
				out.println(sb);
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}

		//print median standard deviation
	}

	/**Converts an array of ExonIntron to an int[exon number][start index in bp positions array, stop index] */
	public int[][] convertToIndexPositions(ExonIntron[] exons){
		ArrayList al = new ArrayList();
		for (int x=0; x< exons.length; x++){
			int startIndex = Num.findClosestStartIndex(positions, exons[x].getStart());
			int endIndex = Num.findClosestEndIndex(positions, exons[x].getEnd()- endOffSet);
			int diff = endIndex-startIndex;
			if (diff < 0) continue;
			al.add(new int[]{startIndex, endIndex});
		}
		if (al.size()==0) return null;
		int[][] ss = new int[al.size()][2];
		for (int x=0; x< ss.length; x++){
			ss[x] = (int[])al.get(x);
		}
		return ss;
	}

	public static void main(String[] args) {
		if (args.length<4){
			printDocs();
			System.exit(0);
		}	
		new ScanGenes(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		String treatment = null;
		String control = null;
		File oligoPositions = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatment = args[i+1]; i++; break;
					case 'c': control = args[i+1]; i++; break;
					case 'r': resultsFile = new File(args[i+1]); i++; break;
					case 'o': oligoPositions = new File(args[i+1]); i++; break;
					case 'g': geneFile = new File(args[i+1]); i++; break;
					case 's': numToSubFromEnds = Integer.parseInt(args[i+1]); i++; break;
					case 'z': endOffSet = Integer.parseInt(args[i+1]) - 1; i++; break;
					case 'p': printIntensities = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//look for required parameters
		if (treatment == null || oligoPositions == null || geneFile == null || control == null){
			Misc.printExit("\nPlease complete one or more of the following required parameters: -g, -t, -c, -o, or -r .\n");
		}

		//parse treatments
		treatmentDirectories = IO.extractFiles(treatment);
		if (treatmentDirectories == null) Misc.printExit("\nProblem parsing treatment directories -> "+treatment);

		//parse control
		controlDirectories = IO.extractFiles(control);
		if (controlDirectories == null) Misc.printExit("\nProblem parsing control directories -> "+control);

		//find chromosome oligo positions
		chromosomeOligoPositions = IO.extractFiles(oligoPositions);
		if (chromosomeOligoPositions == null) Misc.printExit("\nProblem parsing oligo positions -> "+oligoPositions);

	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Scan Genes: Sept 2007                               **\n" +
				"**************************************************************************************\n" +
				"SG parses a UCSC table format gene file, see http://genome.ucsc.edu/cgi-bin/hgTables\n" +
				"(tab delimited: #text chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts\n" +
				"exonEnds) to identify the exons associated with each gene.  These are then used to\n" +
				"fetch the underlying intensity values. A median is then taken of the aveT/aveC ratio\n" +
				"scores for each gene.\n\n" +

				"-o The 'OligoPositions' directory, full path, generated by CelProcessor.\n" +
				"-r The full path file text to use in saving the results.\n" +
				"-t Treatment chip set directories, full path, comma delimited, no spaces.\n" +
				"-c Control chip set directories, full path, comma delimited, no spaces.\n" +
				"-g The full path file text for the UCSC table format gene file.\n" +
				"-s Number to subtract from ends, defaults to 1.  Used to convert UCSC interbase\n"+
				"       numbering to stop inclusive numbering. Set to 0 if already stop inclusive.\n"+
				"-z Size of oligo.  Defaults to 25.\n"+
				"-p Print associated treatment and control intensities.\n"+
				"\n"+

				"Example: java -Xmx1500M -jar pathTo/T2/Apps/ScanGenes -o /affy/OligoPositionsHWG14 -r\n" +
				"      /genes.xls -t /Cels/T/B10_ChrNorm,/Cels/T/B11_ChrNorm -c /Cels/C/B10_ChrNorm,\n" +
				"      /Cels/C/B11_ChrNorm,/Cels/C/B12_ChrNorm -g /ucscHg17RefSeq.txt -p\n\n" +

		"**************************************************************************************\n");		
	}




}
