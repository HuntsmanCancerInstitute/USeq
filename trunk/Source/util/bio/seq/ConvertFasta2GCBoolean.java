package util.bio.seq;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.*;

/**
 * Converts fasta file(s) into serialized boolean[]s where every base g or c is true all others false.
 */
public class ConvertFasta2GCBoolean {
	//fields
	private File directory;
	
	public ConvertFasta2GCBoolean(String[] args){
		
		processArgs(args);
		System.out.println("Converting...");
		
		boolean txt = false;
		File[] files = IO.extractFiles(directory, "binarySeq");
		
		if (files == null || files.length == 0) {
			files = IO.extractFiles(directory);
			txt = true;
		}
		
		for (int i=0; i< files.length; i++){
			System.out.println("\t"+files[i]);
			char[] chromosomeSequence;
			if (txt) {
				MultiFastaParser fastaParser = new MultiFastaParser();
				fastaParser.parseIt(files[i]);
				if (fastaParser.isFastaFound()== false) Misc.printExit("\nError: fasta not found! Aborting.\n");
				chromosomeSequence = fastaParser.getSeqs()[0].toLowerCase().toCharArray();
			}
			else chromosomeSequence = Seq.readBinarySequence(files[i]).toLowerCase().toCharArray();
			
			boolean[] gcs = Seq.fetchGCContent(chromosomeSequence);
			if (gcs == null) Misc.printExit("\nError: fasta not found! Aborting.\n");
			File gcFile = new File(files[i].getParentFile(), Misc.removeExtension(files[i].getName())+".gc");
			if (gcFile.exists()) Misc.printExit("\nError, aborting, "+gcFile +" already exists?!");
			IO.saveObject(gcFile, gcs);
			gcs = null;
		}
		
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
		if (directory==null){
			System.out.println("\nCannot find your directory!\n");
			System.exit(0);
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Convert Fasta 2 GC Boolean: Aug 2008                       **\n" +
				"**************************************************************************************\n" +
				"Converts fasta file(s) into serialized boolean[]s where every base g or c is true all\n" +
				"others false. Will also work with xxx.binarySeq files.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path text for the xxx.fasta file or directory containing such.\n" +
	
				"\n" +
				"Example: java -Xmx2000M -jar pathTo/Apps/ConvertFasta2GCBoolean -f /affy/Fastas/\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new ConvertFasta2GCBoolean(args);
	}

}

