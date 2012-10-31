package util.bio.seq;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import util.bio.parsers.MultiFastaParser;

/**
 * Converts fasta file(s) into compressed binary boolean[]s where GATC is true all (NXgatcnx) others false.
 * Use the BinaryBooleanReaderWriter app to load these compressed files.
 */
public class ConvertFasta2RealSeqBoolean {
	//fields
	private File directory;
	private File saveDirectory;
	private Pattern nonReal = Pattern.compile("[^GATC]");
	private int totalNumRealBases = 0;
	public static final String NAME_NUM_REAL_BASES_FILE = "numberRealBases.txt";
	
	public ConvertFasta2RealSeqBoolean(String[] args){
		
		processArgs(args);
		System.out.println("Converting...");
		
		boolean txt = false;
		File[] files = IO.extractFiles(directory, "binarySeq");
		
		if (files == null || files.length == 0) {
			files = IO.extractFiles(directory, "fasta");
			txt = true;
		}
		
		if (files == null || files.length == 0) Misc.printExit("\nNo files found!.");
		
		//make 
		
		for (int i=0; i< files.length; i++){
			System.out.println("\t"+files[i]);
			String chromosomeSequence;
			if (txt) {
				MultiFastaParser fastaParser = new MultiFastaParser();
				fastaParser.parseIt(files[i]);
				if (fastaParser.isFastaFound()== false) Misc.printExit("\nError: fasta not found! Aborting.\n");
				chromosomeSequence = fastaParser.getSeqs()[0];
			}
			else chromosomeSequence = Seq.readBinarySequence(files[i]);
			
			Matcher mat = nonReal.matcher(chromosomeSequence);
			chromosomeSequence = mat.replaceAll("N");
			
			boolean[] realSeq = new boolean[chromosomeSequence.length()];
			int numReal = 0;
			for (int j=0; j< chromosomeSequence.length(); j++){
				if (chromosomeSequence.charAt(j) != 'N') {
					realSeq[j] = true;
					numReal++;
				}
			}
			totalNumRealBases+= numReal;
			//print stats
			System.out.println("\tLength\t"+chromosomeSequence.length());
			double fraction = (double)numReal / (double)chromosomeSequence.length();
			String percent = Num.formatPercentOneFraction(fraction);
			System.out.println("\tNum Real\t"+numReal+"\t("+percent+")");
			
			//save gc file
			File gcFile = new File(saveDirectory, Misc.removeExtension(files[i].getName())+".realSeq");
			IO.saveObject(gcFile, realSeq);
			IO.zipAndDelete(gcFile);
			realSeq = null;
			chromosomeSequence = null;
		}
		//save number of realBases for downstream applications
		File b = new File (saveDirectory, NAME_NUM_REAL_BASES_FILE);
		IO.writeString(totalNumRealBases+"", b);
		
	}

	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': directory = new File(args[i+1]); i++; break;
					case 's': saveDirectory = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		//check to see if they entered required params
		if (directory==null || directory.isDirectory() == false) Misc.printExit("\nCannot find your genome directory of fasta files!\n");
		if (saveDirectory == null) Misc.printExit("\nPlease enter a directory that should be used in saving the results.\n");
		if (saveDirectory.exists()) System.out.println("\nSave directory exists, potentially overwritting results....");
		else saveDirectory.mkdir();
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Convert Fasta 2 Real Seq Boolean: Dec 2007                   **\n" +
				"**************************************************************************************\n" +
				"Converts fasta file(s) into a compressed binary boolean[]s. GATC bases are true, all\n" +
				"others (gatcNXnx...) are false, case sensitive. Will also work with xxx.binarySeq\n" +
				"files. An info file is also saved containing the total number of real seqs.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path text for the directory containing xxx.fasta files.\n" +
				"-s Save directory, full path. \n"+
	
				"\n" +
				"Example: java util/apps/ConvertFasta2RealSeqBoolean -f /affy/SplitFiles/\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new ConvertFasta2RealSeqBoolean(args);
	}

}

