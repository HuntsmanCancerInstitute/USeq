package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import edu.utah.seq.parsers.jpileup.BamPileupMerger;


/**
 * @author david.nix@hci.utah.edu 
 **/
public class StarFusion2Bed {

	//user defined fields
	private File starFusionFile;
	private File bedOutput;
	private int junctionPadding = 500;
	private File tabix;
	private File bgzip;
	private int numBeds = 0;

	//constructors
	public StarFusion2Bed(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			int numFusions = numBeds/2;
			System.out.println("\nDone! "+Math.round(diffTime)+" sec to parse "+numFusions+" STAR fusions\n");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void doWork() {
		BufferedReader in = null;
		PrintWriter out = null;
		try {
			IO.pl("Parsing...");
			//make IO
			in = IO.fetchBufferedReader(starFusionFile);
			File noGzipBed = new File (Misc.removeExtension(bedOutput.toString())+".bed");
			out = new PrintWriter (new FileWriter(noGzipBed));
			
			
			String[] tokens;
			String line;
			ArrayList<Bed> bedLines = new ArrayList<Bed>();
			while ((line = in.readLine()) != null) {
				if (line.startsWith("#")) {
					out.println(line);
					out.println("#Chr\tStart-"+junctionPadding+"\tStop+"+junctionPadding+"\tSTARFusionLine\tFFPM\tGeneStrand");
				}
				else {
					line = line.trim();
					if (line.length() !=0) {
						//     0            1                   2               3        4        5          6             7              8              9               10             11 
						//#FusionName	JunctionReadCount	SpanningFragCount	est_J	est_S	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots	CDS_LEFT_ID	CDS_LEFT_RANGE	CDS_RIGHT_ID	CDS_RIGHT_RANGE	PROT_FUSION_TYPE	FUSION_MODEL	FUSION_CDS	FUSION_TRANSL	PFAM_LEFT	PFAM_RIGHT
						//IRF4--FOXO1	47	24	47	10.07	ONLY_REF_SPLICE	IRF4^ENSG00000137265.15	chr6:397252:+	FOXO1^ENSG00000150907.10	chr13:40560860:-	YES_LDAS	1.1413	GT	1.9329	AG	1.7819	["INTERCHROMOSOMAL[chr6--chr13]"]	ENST00000493114.1	1-637	ENST00000379561.6	631-1968	FRAMESHIFT	chr6|+|[0]393153-393368[2]|[0]394821-395007[0]|[1]395847-395935[2]|[0]397108-397252[0]<==>chr13|-|[2]40559523-40560860[0]	atgaacctggagggcggcggccgaggcggagagttcggcatgagcgcggtgagctgcggcaacgggaagctccgccagtggctgatcgaccagatcgacagcggcaagtaccccgggctggtgtgggagaacgaggagaagagcatcttccgcatcccctggaagcacgcgggcaagcaggactacaaccgcgaggaggacgccgcgctcttcaaggcttgggcactgtttaaaggaaagttccgagaaggcatcgacaagccggaccctcccacctggaagacgcgcctgcggtgcgctttgaacaagagcaatgactttgaggaactggttgagcggagccagctggacatctcagacccgtacaaagtgtacaggattgttcctgagggagccaaaaaaggagccaagcagctcaccctggaggacccgcagatgtccatgagccacccctacaccatgacaacgccttacccttcgctcccagcccagcaggttcacaactacatgatgccacccctcgaccgaagctggagggactacgtcccggatcagccacacccggaaatcccgtaccaatgtcccatgacgtttggaccccgcggccaccactggcaaggcccagcttgtgaaaatgAATTCAATTCGTCATAATCTGTCCCTACACAGCAAGTTCATTCGTGTGCAGAATGAAGGAACTGGAAAAAGTTCTTGGTGGATGCTCAATCCAGAGGGTGGCAAGAGCGGGAAATCTCCTAGGAGAAGAGCTGCATCCATGGACAACAACAGTAAATTTGCTAAGAGCCGAAGCCGAGCTGCCAAGAAGAAAGCATCTCTCCAGTCTGGCCAGGAGGGTGCTGGGGACAGCCCTGGATCACAGTTTTCCAAATGGCCTGCAAGCCCTGGCTCTCACAGCAATGATGACTTTGATAACTGGAGTACATTTCGCCCTCGAACTAGCTCAAATGCTAGTACTATTAGTGGGAGACTCTCACCCATTATGACCGAACAGGATGATCTTGGAGAAGGGGATGTGCATTCTATGGTGTACCCGCCATCTGCCGCAAAGATGGCCTCTACTTTACCCAGTCTGTCTGAGATAAGCAATCCCGAAAACATGGAAAATCTTTTGGATAATCTCAACCTTCTCTCATCACCAACATCATTAACTGTTTCGACCCAGTCCTCACCTGGCACCATGATGCAGCAGACGCCGTGCTACTCGTTTGCGCCACCAAACACCAGTTTGAATTCACCCAGCCCAAACTACCAAAAATATACATATGGCCAATCCAGCATGAGCCCTTTGCCCCAGATGCCTATACAAACACTTCAGGACAATAAGTCGAGTTATGGAGGTATGAGTCAGTATAACTGTGCGCCTGGACTCTTGAAGGAGTTGCTGACTTCTGACTCTCCTCCCCATAATGACATTATGACACCAGTTGATCCTGGGGTAGCCCAGCCCAACAGCCGGGTTCTGGGCCAGAACGTCATGATGGGCCCTAATTCGGTCATGTCAACCTATGGCAGCCAGGCATCTCATAACAAAATGATGAATCCCAGCTCCCATACCCACCCTGGACATGCTCAGCAGACATCTGCAGTTAACGGGCGTCCCCTGCCCCACACGGTAAGCACCATGCCCCACACCTCGGGTATGAACCGCCTGACCCAAGTGAAGACACCTGTACAAGTGCCTCTGCCCCACCCCATGCAGATGAGTGCCCTGGGGGGCTACTCCTCCGTGAGCAGCTGCAATGGCTATGGCAGAATGGGCCTTCTCCACCAGGAGAAGCTCCCAAGTGACTTGGATGGCATGTTCATTGAGCGCTTAGACTGTGACATGGAATCCATCATTCGGAATGACCTCATGGATGGAGATACATTGGATTTTAACTTTGACAATGTGTTGCCCAACCAAAGCTTCCCACACAGTGTCAAGACAACGACACATAGCTGGGTGTCAGGCTGA	MNLEGGGRGGEFGMSAVSCGNGKLRQWLIDQIDSGKYPGLVWENEEKSIFRIPWKHAGKQDYNREEDAALFKAWALFKGKFREGIDKPDPPTWKTRLRCALNKSNDFEELVERSQLDISDPYKVYRIVPEGAKKGAKQLTLEDPQMSMSHPYTMTTPYPSLPAQQVHNYMMPPLDRSWRDYVPDQPHPEIPYQCPMTFGPRGHHWQGPACENEFNSS*SVPTQQVHSCAE*RNWKKFLVDAQSRGWQEREIS*EKSCIHGQQQ*IC*EPKPSCQEESISPVWPGGCWGQPWITVFQMACKPWLSQQ**L**LEYISPSN*LKC*YY*WETLTHYDRTG*SWRRGCAFYGVPAICRKDGLYFTQSV*DKQSRKHGKSFG*SQPSLITNIINCFDPVLTWHHDAADAVLLVCATKHQFEFTQPKLPKIYIWPIQHEPFAPDAYTNTSGQ*VELWRYESV*LCAWTLEGVADF*LSSP**HYDTS*SWGSPAQQPGSGPERHDGP*FGHVNLWQPGIS*QNDESQLPYPPWTCSADICS*RASPAPHGKHHAPHLGYEPPDPSEDTCTSASAPPHADECPGGLLLREQLQWLWQNGPSPPGEAPK*LGWHVH*ALRL*HGIHHSE*PHGWRYIGF*L*QCVAQPKLPTQCQDNDT*LGVRL	IRF|23-128|2.4e-46^IRF-3|250-406|2.7e-63^MH2|299-372|1.8e-06	FOXO-TAD-PARTIAL|~631-633|8.2e-23
						tokens = Misc.TAB.split(line);
						String noTabsLine = Misc.TAB.matcher(line).replaceAll(" ");
						double score = Double.parseDouble(tokens[11]);
						Bed chromStartStopLeft = parseBedLine(tokens[7], noTabsLine, score);
						Bed chromStartStopRight = parseBedLine(tokens[9], noTabsLine, score);
						bedLines.add(chromStartStopLeft);
						bedLines.add(chromStartStopRight);
					}
				}
			}
			//sort
			IO.pl("Sorting...");
			numBeds = bedLines.size();
			Bed[] toSort = new Bed[numBeds];
			bedLines.toArray(toSort);
			Arrays.sort(toSort);
			
			//write them to file
			for (Bed b: toSort) out.println(b);
		
			//close the IO
			in.close();
			out.close();
			
			//tabix index
			BamPileupMerger.compressAndIndex(bgzip, tabix, noGzipBed, true);
			
			
		} catch (IOException e) {
			e.printStackTrace();
			bedOutput.deleteOnExit();
			IO.closeNoException(in);
			IO.closeNoException(out);
			Misc.printErrAndExit("\nError parsing "+starFusionFile);
		} 
	}

	private Bed parseBedLine(String coor, String noTabsLine, double score) {
		//chr6:397252:+
		String[] t = Misc.COLON.split(coor);
		int breakPoint = Integer.parseInt(t[1]);
		int start = breakPoint- junctionPadding;
		if (start < 0) start = 0;
		int end = breakPoint + junctionPadding;
		return new Bed(t[0], start, end, noTabsLine, score, t[2].charAt(0));
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new StarFusion2Bed(args);
	}		

	/**This method will process each argument and assign new varibles
	 * @throws IOException 
	 * @throws FileNotFoundException */
	public void processArgs(String[] args) throws FileNotFoundException, IOException{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		File tabixBinDirectory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': starFusionFile = new File(args[++i]).getCanonicalFile(); break;
					case 'b': bedOutput = new File(args[++i]); break;
					case 'p': junctionPadding = Integer.parseInt(args[++i]); break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (starFusionFile == null || starFusionFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your "
				+ "starFusion file, e.g. -s star-fusion.fusion_predictions.abridged.coding_effect.tsv \n");

		if (bedOutput == null) Misc.printErrAndExit("\nError: please provide a bed results output file, -b xxx.bed.gz\n");
		
		//pull tabix and bgzip
		if (tabixBinDirectory == null) Misc.printExit("\nError: please point to your HTSlib directory containing the tabix and bgzip executables (e.g. -t ~/BioApps/HTSlib/1.10.2/bin/ )\n");
		bgzip = new File (tabixBinDirectory, "bgzip");
		tabix = new File (tabixBinDirectory, "tabix");
		//look for bgzip and tabix executables
		if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute the bgzip '"+bgzip+"' or tabix executables '"+tabix+"'\n");

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          STAR Fusion 2 Bed: June 2020                      **\n" +
				"**************************************************************************************\n" +
				"Parses STAR-Fusion output tsv files into bed format. Two bed entries are made per\n"+
				"fusion to represent the left and right breakpoints +/- junction padding. The entire\n"+
				"fusion line is placed in the bed name column, and the FFPM in the score column. See\n" +
				"https://github.com/STAR-Fusion/STAR-Fusion/wiki#Outputs for details. The final bed is\n"+
				"bgzip compressed and tabix indexed.\n"+

				"\nRequired Options:\n"+
				"-s Path to one of the xxx.tsv STAR-Fusion output files, recommend using\n" +
				"      star-fusion.fusion_predictions.abridged.coding_effect.tsv \n"+
				"-b Path to the xxx.bed.gz output results file\n"+
				"-t Path to the directory containing the compiled bgzip and tabix executables. See\n" +
				"     https://www.htslib.org\n"+

				"\nDefault Options:\n"+
				"-p BP junction padding, defaults to 500\n"+

				"\nExample: java -Xmx20G -jar pathToUSeq/Apps/StarFusion2Bed -b P13445.sf.bed.gz\n" +
				"     -s P13445_SF/star-fusion.fusion_predictions.abridged.coding_effect.tsv  \n\n" +

				"**************************************************************************************\n");

	}


}
