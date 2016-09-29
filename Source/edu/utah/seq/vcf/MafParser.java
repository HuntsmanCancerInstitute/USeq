package edu.utah.seq.vcf;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.annotation.Bed;
import util.gen.*;


/**
 * Parses variant maf files. Fixes header, sorts, indexes
 * @author davidnix*/
public class MafParser {

	//user defined fields
	private File[] mafFiles;
	private File outputDir;
	private File bgzip;
	private File tabix;


	public MafParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		System.out.println("Processing maf files...");
		for (int i=0; i< mafFiles.length; i++) {
			File parsed = parseFile(mafFiles[i]);
			if (tabix != null) tabix(parsed);
		}
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	
	
	private void tabix(File parsed) {
		//force compress it
		String[] cmd = new String[]{
				bgzip.toString(),
				"-f",
				parsed.toString()
		};
		String[] messages = IO.executeViaProcessBuilder(cmd, false);
		if (messages == null || messages.length !=0) Misc.printErrAndExit("\nERROR: trying to execute bgzip compression -> "+
		Misc.stringArrayToString(cmd, " ")+"\nMessage: "+Misc.stringArrayToString(messages,  "\n"));
		
		//force compress it
		cmd = new String[]{
				tabix.toString(),
				"-f", "-s", "5", "-b", "6", "-e", "7",
				parsed.toString()+".gz"
		};
		messages = IO.executeViaProcessBuilder(cmd, false);
		if (messages == null || messages.length !=0) Misc.printErrAndExit("\nERROR: trying to execute tabix -> "+
		Misc.stringArrayToString(cmd, " ")+"\nMessage: "+Misc.stringArrayToString(messages,  "\n"));
	}



	public File parseFile(File mpfFile){
		try {
			BufferedReader in = IO.fetchBufferedReader(mpfFile);
			StringBuilder header = new StringBuilder();
			ArrayList<Bed> bedAL = new ArrayList<Bed>();
			
			//for each line in the file
			String line;
			while ((line = in.readLine()) != null){
				if (line.length() == 0) continue;
				if (line.startsWith("#")) {
					header.append(line);
					header.append("\n");
				}
				else if (line.startsWith("Hugo_Symbol")) {
					header.append("#");
					header.append(line);
					header.append("\n");
					header.append(line);
					header.append("\n");
				}
				else {
					String[] t = Misc.TAB.split(line);
					String chr = t[4];
					int start = Integer.parseInt(t[5]);
					int stop = Integer.parseInt(t[6]);
					bedAL.add(new Bed(chr, start, stop, line, 0, '.'));
				}
			}
			in.close();
			System.out.println(mpfFile.getName()+"\t"+bedAL.size());
			
			//sort
			Bed[] bed = new Bed[bedAL.size()];
			bedAL.toArray(bed);
			Arrays.sort(bed);
			
			//print
			String name = mpfFile.getName();
			if (name.endsWith(".gz")) name = name.substring(0, name.length()-3);
			File fixed = new File(outputDir, name);
			PrintWriter out = new PrintWriter( new FileWriter( fixed ));
			out.print(header.toString());
			for (Bed b: bed) out.println(b.getName());
			out.close();
			
			return fixed;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MafParser(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		File tabixBinDirectory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': forExtraction = new File(args[++i]); break;
					case 'o': outputDir = new File(args[++i]); break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a maf file or directory containing such to parse.\n");
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".maf");
		tot[1] = IO.extractFiles(forExtraction,".maf.gz");
		tot[2] = IO.extractFiles(forExtraction,".maf.txt");
		tot[3] = IO.extractFiles(forExtraction,".maf.txt.gz");
		mafFiles = IO.collapseFileArray(tot);
		if (mafFiles == null || mafFiles.length ==0 || mafFiles[0].canRead() == false) Misc.printExit("\nError: cannot find any xxx.maf(.txt/.gz OK) file(s)!\n");

		//tabix
		if (tabixBinDirectory != null) {
			bgzip = new File (tabixBinDirectory, "bgzip");
			tabix = new File (tabixBinDirectory, "tabix");
			//look for bgzip and tabix executables
			if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+bgzip+" "+tabix);
		}

		
		//outputDir
		if (outputDir == null) Misc.printExit("\nError: please provide an output directory for writing the fixed files.\n");
		outputDir.mkdirs();
		
		
	}	
	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Maf Parser: Sept 2016                             **\n" +
				"**************************************************************************************\n" +
				"Parses and manipulates variant maf files. Sorts and adds a # header line. Provide a\n"+
				"path to the tabix executables (https://github.com/samtools/htslib) for indexing and\n"+
				"IGV browsing.\n\n"+

				"Options:\n"+
				"-m Path to a maf file (.maf/.txt/.gz OK) or directory containing such.\n"+
				"-o Output directory, will overwrite.\n" +
				"-t To tabix index the output, provide a path to the dir containing bgzip and tabix\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MafParser -m MafTCGAFiles/ -o Sorted/ \n"+
				"              -t ~/BioApps/HTSlib/1.3/bin/ \n\n" +

		        "**************************************************************************************\n");
	}

}
