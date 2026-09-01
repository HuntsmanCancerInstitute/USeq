package edu.utah.seq.vcf.serial;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.vcf.fdr.VCFFdrEstimator;
import util.gen.*;

public class VCFSerialSampleInterpreter {

	//user fields
	private File inputSpreadsheet = null;
	private boolean verbose = false;
	private double minimalVaf = 0.0;
	private boolean justVariantsInA = false;
	private File saveDirectory = null;
	private double minimumFracChange = 0.05;
	private File pathToRScript = null;
	
	//Internal
	TreeMap<String, Patient> patients = new TreeMap<String, Patient>();

	
	//constructor
	public VCFSerialSampleInterpreter(String[] args) throws Exception{
		//start clock
		long startTime = System.currentTimeMillis();
		processArgs(args);

		loadSpreadsheet();
		
		comparePatientSamples();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	

	
	
	private void comparePatientSamples() {
		IO.pl("\nComparing patient sample sets, min VAF "+minimalVaf+", use only variants in A "+justVariantsInA+" ...");
		if (verbose == false) IO.pl("\nPatientID\tSampleA\tSampleB\t#PosPairsA-B\t#NegPairsA-B\tPVal\tMeanVafA\tMeanVafB\tmVafB/mVafA\tmVafB/mVafA-1");;
		//for each patient
		for (Patient p: patients.values()) {
			if (verbose) IO.pl("\nPatient\t"+p.getId());
			//p.compareSamples(verbose);
			p.compareSamples2(verbose);
		}
	}





	private void loadSpreadsheet() throws Exception {
		IO.pl("Parsing spreadsheet...");
		String[] lines = IO.loadFileIntoStringArray(inputSpreadsheet);
		int numDataLines = 0;
		HashSet<String> allLines = new HashSet<String>();
		int duplicates = 0;
		for (String l : lines) {
			if (l.startsWith("Patient")) continue;
			if (allLines.contains(l)) {
				duplicates++;
				continue;
			}
			allLines.add(l);
			
			//PatientID SampleName SampleNumber VariantID VariantConseq FoundInSample Gene VAF DP
			//    0          1           2           3         4              5         6   7   8
			String[] f = Misc.TAB.split(l);
			if (f.length != 9) Misc.printErrAndExit("\nFAILED to find 9 fields in "+l);
			
			// fetch or make patient
			Patient p = patients.get(f[0]);
			if (p==null) {
				p= new Patient(f[0], minimalVaf, justVariantsInA, saveDirectory, minimumFracChange, pathToRScript);
				patients.put(f[0], p);
			}
			
			// fetch or make sample
			Integer sampleNumber = Integer.parseInt(f[2]);
			Sample s = p.getSamples().get(sampleNumber);
			if (s == null) {
				s = new Sample(sampleNumber, f[1]);
				p.getSamples().put(sampleNumber, s);
			}
			
			// add a variant
			s.getVariants().put(f[3], new Variant(f));
			numDataLines++;
		}
		IO.pl("\tUnique data lines\t"+numDataLines);
		IO.pl("\tDuplicate data lines\t"+duplicates);
	}

	public static void main(String[] args) throws Exception {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFSerialSampleInterpreter(args);
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
				try {
					switch (test){
					case 'i': inputSpreadsheet = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': pathToRScript = new File(args[++i]); break;
					case 'm': minimalVaf = Double.parseDouble(args[++i]); break;
					case 'd': minimumFracChange = Double.parseDouble(args[++i]); break;
					case 'v': verbose = true; break;
					case 'j': justVariantsInA = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (inputSpreadsheet == null || inputSpreadsheet == null || saveDirectory == null || pathToRScript == null) {
			printDocs();
			System.exit(1);
		}
		saveDirectory.mkdirs();
	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          VCF Serial Sample Interpreter : July 2026               **\n" +
				"**************************************************************************************\n" +
				"Takes the output of the VSS Tracker to calculate comparisons between samples within\n"+
				"each patient. A paired Wilcoxon test is performed on the VAFs from a unique set of\n"+
				"variants for each comparison. Each variant in the set must have been found in at least\n"+
				"one of the original call files and be >= the minimum VAF.\n"+

				"\nRequired Params:\n"+
				"-i Input spreadsheet txt file from the VCFSerialSampleTracker.\n"+
				"-s Save directory for R plots.\n"+
				"-m Minimum VAF for founder variants.\n"+
				"-d Minimum fraction change for inclusion in stats, defaults to 0.05\n"+
				"-r Path to the twoSampleComparatorPlots.R script\n"+
				"-j Just use founder variants in sample A.\n"+
				"-v Verbose output.\n"+
				
				"\nExample: java -jar pathTo/USeq/Apps/VCFSerialSampleInterpreter -i input.txt\n" +
				"       -m 0.01 -s TxtOutputForR -r ~/VariantAnalysis/twoSampleComparatorPlots.R\n"+

		"\n**************************************************************************************\n");

	}
}
