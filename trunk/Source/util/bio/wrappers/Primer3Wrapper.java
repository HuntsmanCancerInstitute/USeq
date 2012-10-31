package util.bio.wrappers;

import java.io.*;

import util.gen.*;

import java.util.*;
import java.util.regex.*;

/**This app wraps the command line primer3 app to make it useful for picking qPCR primers.*/
public class Primer3Wrapper {

	//fields
	private File primer3App = new File ("/nfs/transcriptome/software/noarch/T2/64Bit_Primer3_1.0.0/src/primer3_core");
	private File misPrimingLibrary = new File ("/nfs/transcriptome/software/noarch/T2/64Bit_Primer3_1.0.0/cat_humrep_and_simple.cgi.txt");
	private File seqFile = null;
	private boolean pickLarge = true;
	private Primer3Pick[] picks = null;
	
	
	//options
	private String standard = 
		"PRIMER_MIN_SIZE=20\n"+
		"PRIMER_MIN_TM=55\n"+
		"PRIMER_MAX_TM=68\n"+
		"PRIMER_MIN_GC=40\n"+
		"PRIMER_MAX_GC=70\n"+
		"PRIMER_NUM_RETURN=1\n"+
		"PRIMER_TASK=pick_pcr_primers\n"+
		"PRIMER_MAX_DIFF_TM=7.5\n";
	
	private String small=
		"PRIMER_OPT_SIZE=22\n"+
		"PRIMER_MAX_SIZE=24\n"+
		"PRIMER_PRODUCT_SIZE_RANGE=45-80\n";
		
	private String large=
		"PRIMER_OPT_SIZE=24\n"+
		"PRIMER_MAX_SIZE=30\n"+
		"PRIMER_PRODUCT_SIZE_RANGE=80-150\n"+
		"PRIMER_GC_CLAMP=2\n";
	
	
	/**Uses a shell script file to execute the primer3 app due to a problem feeding to STDIN from java.*/
	public Primer3Wrapper(String[] args){
		processArgs(args);
		//parse picks
		picks = Primer3Pick.parseFile(seqFile);
		//write temp file for primer3
		String formattedSeqs = buildExecuteDoc(picks);
		String random = Passwords.createRandowWord(10);
		File formattedSeqFile = new File (seqFile.getParentFile(), "P3WTempFile_"+ random +"Seq.txt");
		IO.writeString(formattedSeqs, formattedSeqFile);
		//make shell script
		String script = IO.getFullPathName(primer3App) + " < " + IO.getFullPathName(formattedSeqFile);
		//execute shell script
		String[] results = IO.executeShellScript(script, seqFile.getParentFile());
		if (results == null){
			Misc.printExit("\nError: problem launching shell file. Try "+script+" in a xxx.sh file");
		}
		else formattedSeqFile.delete();
		//parse results
		parsePrimer3Results(results);
		//save results
		File out = new File (IO.getFullPathName(seqFile)+"_P3Res.xls");
		StringBuffer sb = new StringBuffer(Primer3Pick.header);
		sb.append("\n");
		for (int i=0; i< picks.length; i++){
			sb.append(picks[i].getResultSummaryLine());
			sb.append("\n");
		}
		IO.writeString(sb.toString(), out);
		System.out.println("Done!");

	}
		
	/**Takes the output of the primer3 app and breaks by record loading appriopriate Primer3Pick.*/
	public void parsePrimer3Results(String[] results){
		//read in results, one record at a time, until hit a +
		HashMap map = new HashMap();
		for (int i=0; i< results.length; i++){
			//stop of record?, load pick
			if (results[i].startsWith("=")){
				Object obj = map.get("PRIMER_SEQUENCE_ID");
				if (obj == null) Misc.printExit("\nError: problem parsing primer3 results \n"+map);
				int index = Integer.parseInt((String)obj);
				picks[index].loadResults(map);
				map.clear();
			}
			else {
				String[] keyValue = results[i].split("=");
				map.put(keyValue[0], keyValue[1]);
			}
		}
	}
	
	public String buildExecuteDoc(Primer3Pick[] pp){
		StringBuffer sb = new StringBuffer();
		String primeLib = "\nPRIMER_MISPRIMING_LIBRARY="+IO.getFullPathName(misPrimingLibrary)+"\n";
		for (int i=0; i< pp.length; i++){
			sb.append("PRIMER_SEQUENCE_ID=");
			sb.append(i);
			sb.append("\nSEQUENCE=");
			sb.append(pp[i].getSequence());
			sb.append(primeLib);
			sb.append(standard);
			if (pickLarge) sb.append(large);
			else sb.append(small);
			sb.append("=\n");
		}
		return sb.toString();
	}
	
	
	public static void main(String[] args){
		//look for file
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		new Primer3Wrapper(args);
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
					case 'p': primer3App = new File(args[i+1]); i++; break;
					case 'm': misPrimingLibrary = new File(args[i+1]); i++; break;
					case 'f': seqFile = new File(args[i+1]); i++; break;
					case 's': pickLarge = false; break;
					case 'h': printDocs(); System.exit(0);
					default : Misc.printExit("\nError: unknown option! -> " + mat.group()+"\n");
					}
				}
				catch (Exception e){
					Misc.printExit("\nError: something doesn't look right with this parameter request: -"+test+"\n");
				}
			}
		}
		//check for files
		if (seqFile == null || seqFile.exists()== false) Misc.printExit("\nError: please enter a sequence file from which to pick primers.");
		if (misPrimingLibrary == null || (misPrimingLibrary != null && misPrimingLibrary.exists()== false)) Misc.printExit("\nError: cannot find your mis priming library file "+misPrimingLibrary);
		if (primer3App == null || (primer3App != null && primer3App.exists()== false)) Misc.printExit("\nError: cannot find your primer3_core application "+primer3App);
		
		//print options
		System.out.println("\nLaunching with the following options:");
		System.out.println("\tSequence File\t"+seqFile);
		System.out.println("\tMis priming library\t"+misPrimingLibrary);
		System.out.print("\tPrimer3 options for ");
		if (pickLarge) System.out.println("STANDARD size products\n\n"+large+standard);
		else System.out.println("SMALL size products\n\n"+small+standard);
		

	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Primer3 Wrapper: Dec  2006                            **\n" +
				"**************************************************************************************\n" +
				"Wrapper for the primer3 application. Extracts sequence, formats for primer3, executes\n" +
				"and parses the output to a spreadsheet. See http://frodo.wi.mit.edu/primer3/\n\n"+

				"-f Full path file text for your sequence file, tab delimited, sequence in 1st column.\n"+
				"-s Pick small product sizes (45-80bp), defaults to standard (80-150bp)\n"+
				"-p Full path file text for the primer3_core application. Defaults to\n" +
				"     /nfs/transcriptome/software/noarch/T2/64Bit_Primer3_1.0.0/src/primer3_core\n"+
				"-m Full path file text for the mispriming library. Defaults to\n"+
				"     /nfs/transcriptome/software/noarch/T2/64Bit_Primer3_1.0.0/\n" +
				"     cat_humrep_and_simple.cgi.txt\n\n"+
				
				"Example: java -jar pathTo/T2/Apps/Primer3Wrapper -f /home/dnix/seqForQPCR.txt -p\n" +
				"    /nfs/transcriptome/software/noarch/T2/64Bit_Primer3_1.0.0/src/primer3_core\n" +
				"    -m /nfs/transcriptome/software/noarch/T2/64Bit_Primer3_1.0.0/\n" +
				"    cat_humrep_and_simple.cgi.txt -s \n"+
				
		"**************************************************************************************\n");
	}

}
