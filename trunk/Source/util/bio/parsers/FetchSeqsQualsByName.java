package util.bio.parsers;
import java.io.*;
import util.gen.*;
import java.util.*;

public class FetchSeqsQualsByName {
	

	File[] seqNameFiles;
	File multiFasta;
	File multiQual;
	
	
	public FetchSeqsQualsByName( String[] args){
		//pull files with seq names, fasta headers
		seqNameFiles = IO.extractFiles(new File (args[0]));
		//pull fasta file and quality file
		multiFasta = new File (args[1]);
		multiQual = new File (args[2]);
		//parse fasta
		System.out.println("\nParsing fastas sequences...");
		parseFasta();
		//parse quality file
		System.out.println("\nParsing quality scores...");
		parseQuality();

	}

	public void parseFasta(){
		MultiFastaParser mfp = new MultiFastaParser (multiFasta);
		HashMap nameSeqs = mfp.getNamesSeqs();
		for (int i=0; i< seqNameFiles.length; i++){
			File subFasta = new File (seqNameFiles[i]+".fasta");
			try {
				PrintWriter out = new PrintWriter( new FileWriter (subFasta));
				String[] hitNames = IO.loadFileIntoStringArray(seqNameFiles[i]);
				for (int j=0; j< hitNames.length; j++){
					hitNames[j] = hitNames[j].trim();
					if (nameSeqs.containsKey(hitNames[j])){
						out.println(">"+hitNames[j]);
						String seq = (String)nameSeqs.get(hitNames[j]);
						out.println(seq);
					}
					else System.out.println("\t"+ hitNames[j]+"\t Not Found!");
				}
			out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		mfp = null;
		nameSeqs = null;
		System.gc();
	}
	
	public void parseQuality(){
		QualityFileParser mfp = new QualityFileParser (multiQual.toString());
		HashMap nameQuals = mfp.getQualNames();
		for (int i=0; i< seqNameFiles.length; i++){
			File subQual = new File (seqNameFiles[i]+".qual");
			try {
				PrintWriter out = new PrintWriter( new FileWriter (subQual));
				String[] hitNames = IO.loadFileIntoStringArray(seqNameFiles[i]);
				for (int j=0; j< hitNames.length; j++){
					hitNames[j] = hitNames[j].trim();
					if (nameQuals.containsKey(hitNames[j])){
						out.println(">"+hitNames[j]);
						String seq = (String)nameQuals.get(hitNames[j]);
						out.println(seq);
					}
					else System.out.println("\t"+ hitNames[j]+"\t Not Found!");
				}
			out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		mfp = null;
		nameQuals = null;
		System.gc();
	}
	
	public static void main(String[] args) {
		if (args.length == 0) Misc.printExit("\nEnter a full path file text for a fasta header text file/directory, fasta file, and quality file\n");
		new FetchSeqsQualsByName( args );

	}

}
