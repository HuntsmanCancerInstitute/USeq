package edu.utah.seq.vcf.fdr;

import java.io.BufferedReader;
import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Takes scored fdr values and applies a tuned thresholding based on regression.
 * @author Nix*/
public class VCFFdrFilter {
	
	//user defined fields
	private File forExtraction = null;
	private File[] vcfFiles;
	private File saveDir;
	private String afInfoName = "T_AF";
	//private double[] afs;
	//private double[] thresholds;
	private double[][] snvBrackets;
	private double[][] indelBrackets;
	private Pattern AF = null;
	private boolean error = false;
	private String cmd = null;
	
	//5% FDR
	double[] snvAfs = null;
	double[] snvThresholds = null;
	double[] indelAfs = null;
	double[] indelThresholds = null;

	//working
	private File vcfFile;
	private int numRecords;
	private int numPassing;
	
	//constructor
	public VCFFdrFilter(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		//make brackets
		snvBrackets = new double[snvThresholds.length-1][2];
		for (int i=0; i< snvBrackets.length; i++) snvBrackets[i] =  Num.calculateSlopeAndYIntercept(snvAfs[i], snvThresholds[i],    snvAfs[i+1], snvThresholds[i+1]);
		indelBrackets = new double[indelThresholds.length-1][2];
		for (int i=0; i< indelBrackets.length; i++) indelBrackets[i] =  Num.calculateSlopeAndYIntercept(indelAfs[i], indelThresholds[i],    indelAfs[i+1], indelThresholds[i+1]);

		//for each vcf file
		IO.pl("\nParsing vcf files...");
		IO.pl("Name\tRecords\tPassed\tFailed");
		for (int i=0; i< vcfFiles.length; i++){
			//set values and clear past data
			vcfFile = vcfFiles[i];
			numRecords = 0;
			numPassing = 0;
			parse();
			IO.pl(vcfFile.getName()+"\t"+numRecords+"\t"+numPassing+"\t"+(numRecords-numPassing));
			if (error) System.exit(1);
		}
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		IO.pl("\nDone! "+Math.round(diffTime)+" min\n");
	}

	


	private void parse() {
		Gzipper passOut = null;
		Gzipper failOut = null;
		BufferedReader vcfIn = null;
		try {
			//make IO
			String name = Misc.removeExtension(vcfFile.getName());
			passOut = new Gzipper (new File (saveDir, name+".pass.vcf.gz"));
			failOut = new Gzipper (new File (saveDir, name+".fail.vcf.gz"));
			vcfIn = IO.fetchBufferedReader(vcfFile);
			
			//parse the file
			String line;
			while ((line = vcfIn.readLine())!= null) {
				if (line.startsWith("#")) {
					if (line.startsWith("#CHROM")) line = "##" +cmd+"\n"+line;
					passOut.println(line);
					failOut.println(line);
				}
				else {
					if (score(line)) passOut.println(line);
					else failOut.println(line);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			error = true;
		} finally {
			passOut.closeNoException();
			failOut.closeNoException();
			IO.closeNoException(vcfIn);
		}
	}




	private boolean score(String record) throws Exception {
		numRecords++;
		
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
		String[] fields = Misc.TAB.split(record);		

		//fetch AF
		Matcher mat = AF.matcher(fields[7]);
		if (mat.find() == false) throw new Exception("\nERROR: Failed to parse AF= number from the INFO field in this variant -> "+record+" aborting.\n");
		
		//IO.pl("\nSeven "+fields[7]);
		double af = Double.parseDouble(mat.group(1));
		double qual = Double.parseDouble(fields[5]);
		
		//snv or indel
		if (fields[3].length()!=1 || fields[4].length()!=1) return threshold(af, qual, indelAfs, indelThresholds, indelBrackets);
		return threshold(af, qual, snvAfs, snvThresholds, snvBrackets);
		
	}
	
	/**Applies the score cut off from the NA12878 simulation points using interpolation.*/
	public boolean threshold(double af, double score, double[] afs, double[] thresholds, double[][] brackets){
		
		//calculate threshold to pass given the variant allele freq 
		double threshold = -1;
		//before first bracket?
		if (af <=afs[0])threshold = thresholds[0];
		//inside one of the brackets?
		if (threshold == -1) {
			for (int i=0; i< brackets.length; i++) {
				if (af>afs[i] && af<= afs[i+1]) {
					threshold = Num.calculateYGivenX(brackets[i], af);
					break;
				}
			}
		}
		//greater than last bracket?
		if (threshold == -1) threshold = thresholds[thresholds.length-1];
		
		//System.out.println("\tAF:"+af+" Score:"+score+" Thresh:"+threshold + " Pass:"+(score >= threshold));
		if (score >= threshold) {
			numPassing++;
			return true;
		}
		return false;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFFdrFilter(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		cmd = IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ");
		IO.pl("\n"+cmd+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 's': saveDir = new File(args[++i]); break;
					case 'a': snvAfs = Num.stringArrayToDouble(args[++i], ","); break;
					case 'q': snvThresholds = Num.stringArrayToDouble(args[++i], ","); break;
					case 'b': indelAfs = Num.stringArrayToDouble(args[++i], ","); break;
					case 'r': indelThresholds = Num.stringArrayToDouble(args[++i], ","); break;
					case 'n': afInfoName = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		//check AFs and Thresholds
		if (snvAfs == null || snvThresholds == null || indelAfs == null || indelThresholds == null) Misc.printErrAndExit("\nError: provide a set of matched allele frequencies and QUAL thresholds for filtering.\n");
		if (snvAfs.length != snvThresholds.length || indelAfs.length != indelThresholds.length) Misc.printErrAndExit("\nError: the number of allele frequencies and QUAL thresholds  must be the same.\n");
		
		//save dir?
		if (saveDir == null) saveDir = vcfFiles[0].getParentFile();
		else saveDir.mkdirs();
		if (saveDir.exists()==false || saveDir.canWrite()==false) Misc.printExit("\nError: failed to find or create a writeable save directory? Aborting.\n");
		
		//pattern
		AF = Pattern.compile( afInfoName+"=([\\d\\.]+)");
		printSettings();
		
	}	
	
	public void printSettings(){
		IO.pl("Settings:");
		IO.pl(" -v Vcf file  "+ forExtraction);
		IO.pl(" -s Save dir  "+ IO.getCanonicalPath(saveDir));
		IO.pl(" -q SNV QUAL thresholds    "+ Num.doubleArrayToString(snvThresholds, ","));
		IO.pl(" -a SNV Allele frequencies "+ Num.doubleArrayToString(snvAfs, ","));
		IO.pl(" -r INDEL QUAL thresholds    "+ Num.doubleArrayToString(indelThresholds, ","));
		IO.pl(" -b INDEL Allele frequencies "+ Num.doubleArrayToString(indelAfs, ","));
		IO.pl(" -n Tumor AF INFO name "+ afInfoName);
		
	} 
	
	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                               VCF Fdr Filter : Jan 2020                          **\n" +
				"**************************************************************************************\n" +
				"Takes scored fdr values and applies a tuned thresholding based on regression to filter\n"+
				"vcf files to reach a list with a targeted fdr. Each threshold is matched to a given \n"+
				"variant AF and variant type. For example, filtering snv vcf records with tumor AFs of\n"+
				"0.05 with a QUAL score of 8.64 produces lists with a FDR of 5%\n"+

				"\nRequired:\n"+
				"-v Path to a vcf file xxx.vcf(.gz/.zip OK) or directory containing such.\n" +
				"-s Path to a directory to save the pass and fail vcf records.\n"+
				"-q SNV QUAl thresholds that generate lists of variants with a specific FDR. Comma\n"+
				"      delimited, no spaces.\n"+
				"-a SNV Allele frequencies.\n"+
				"-r INDEL QUAl thresholds.\n"+
				"-b INDEL Allele frequencies.\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/VCFFdrFilter -v Calls/ -s FDR5Calls -q \n"+
				"      8.75,8.63,8.62,8.64,8.85 -r 16.5,15.9,15.3,25.7,45.2 -a\n" + 
				"      0.005,0.0075,0.01,0.05,0.1 -b 0.005,0.0075,0.01,0.05,0.1 \n\n"+

		        "**************************************************************************************\n");
	}

}
