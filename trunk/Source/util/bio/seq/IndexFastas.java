package util.bio.seq;
import java.io.*;
import java.util.regex.*;
import util.gen.*;
import util.bio.parsers.*;

public class IndexFastas {
	
	//fields
	File[] fastas = null;
	int numberOfBases = 1000000;
	File indexDirectory;
	MultiFastaParser parser;
	
	public void index(){
		//for each Fasta
		for (int i=0; i< fastas.length; i++){
			System.out.println("Processing "+fastas[i].getName());
			//attempt to parse
			parser = new MultiFastaParser(fastas[i]);
			if (parser.isFastaFound() == false) {
				System.out.println("Fasta seq not found, skipping -> "+fastas[i]);
				continue;
			}
			//for each parsed seq, index it
			int numSeqs = parser.getNumReads();
			for (int j=0; j< numSeqs; j++){
				//check text
				String name = parser.getNames()[j];
				if (name.startsWith("chr") == false){
					System.out.println("WARNING, the fasta seq text doesn't begin with 'chr', check ->  "+fastas[i]+" "+name);
					continue;
				}
				//index it
				splitIndex(name, parser.getSeqs()[j]);
			}
		}
	}
	
	public void splitIndex (String name, String seq){
		int start = 0;
		boolean go = true;
		while (go) {
			int end = start + numberOfBases;
			if (end >= seq.length()) {
				end = seq.length();
				go = false;
			}
			String subSeq = seq.substring(start, end);
			File binarySeq = new File (indexDirectory, name+"_"+start+"-"+(end-1));
			if (binarySeq.exists()) {
				System.out.println("WARNING, "+binarySeq+" already exists, skipping!");
				return;
			} 
			Seq.writeBinarySequence(subSeq, binarySeq);
			start = end;
		}
	}
	
	public IndexFastas(String[] args){
		processArgs(args);
		index();
		System.out.println("\nDone!");
	}
	
	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}	
		new IndexFastas(args);
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
					case 'f': fastas = IO.extractFiles(args[i+1], "fasta"); i++; break;
					case 'n': numberOfBases = Integer.parseInt(args[i+1]); i++; break;
					case 'i': indexDirectory = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//make index directory?
		if (indexDirectory == null) {
			indexDirectory = new File (fastas[0].getParentFile(), "IndexedSequences");
			indexDirectory.mkdir();
		}
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Index Fastas: April 2007                              **\n" +
				"**************************************************************************************\n" +
				"IF is used to break apart long sequences into small binary files for rapid access.\n\n"+
				
				"-f Full path file text to a xxx.fasta sequence or directory containing such.\n"+
				"-n Number of bases permitted in each indexed file, defaults to 1,000,000.\n"+
				"-i Alternative save directory, defaults to fasta directory.\n\n"+
				
				"\nExample: java -Xmx1500M -jar pathTo/T2/Apps/IndexFastas -f /genomes/hg17/\n\n" +
				
		"**************************************************************************************\n");		
	}
}
