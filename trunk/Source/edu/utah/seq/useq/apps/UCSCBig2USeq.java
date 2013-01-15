package edu.utah.seq.useq.apps;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;
import edu.utah.seq.useq.USeqUtilities;

/**Class to convert UCSC xxx.bb or xxx.bw archives to xxx.useq. Must know the genome build version.*/
public class UCSCBig2USeq extends Thread{

	private File[] ucscArchives;
	private File bigWigToWig;
	private File bigBedToBed;
	private boolean verbose = true;
	private boolean forceConversion = false;
	private File convertedUSeqArchive = null;
	private File ucscArchiveToConvert;
	private File tempTxtFile = null;
	private String versionedGenome = null;
	private Pattern tab = Pattern.compile("\\t");
	private Pattern bedLine = Pattern.compile("^\\w+\\t\\d+\\t\\d+.*$");
	private Pattern badName = Pattern.compile("[.\\s-]");
	private boolean deleteTempFiles = true;

	//bed typing
	private boolean stranded = false;
	private boolean scored = false;
	private boolean named = false;
	private boolean bed12 = false;


	//constructors
	//stand alone
	public UCSCBig2USeq (String[] args){
		try {
			processArgs(args);
			if (verbose) System.out.println("Processing:");
			
			//remove those that already exist?
			if (forceConversion == false) {
				ucscArchives = UCSCBig2USeq.removeExistingConvertedBigFiles(ucscArchives);
				if (ucscArchives.length == 0) {
					if (verbose) System.out.println("\tNo unconverted xxx.bb/bw archives were found.  Use the -f option to overwrite.\n");
					System.exit(0);
				}
			}
			
			//for each archive
			for (File u : ucscArchives){
				ucscArchiveToConvert = u;
				if (verbose) System.out.println("\t"+ucscArchiveToConvert);
				convert();
				if (deleteTempFiles) tempTxtFile.delete();
			}

			if (verbose) System.out.println("\nDone!\n");
		} catch (Exception e){
			e.printStackTrace();
			System.exit(1);
		}

	}


	//methods
	/**Returns converted useq archive or null if something bad happened.*/
	public File convert () throws Exception{

		//bw
		if (ucscArchiveToConvert.getName().endsWith(".bw")){
			//convert big wig to wig
			tempTxtFile = new File(Misc.removeExtension(ucscArchiveToConvert.getCanonicalPath())+".wig");
			String[] cmd = {bigWigToWig.getCanonicalPath(), ucscArchiveToConvert.getCanonicalPath(), tempTxtFile.getCanonicalPath()};
			executeUCSCCommand(cmd);
			//convert wig to useq
			Wig2USeq w2u = new Wig2USeq(tempTxtFile, versionedGenome, verbose, true);
			convertedUSeqArchive = w2u.getUseqArchive();
		}

		//bb
		else if (ucscArchiveToConvert.getName().endsWith(".bb")){
			//convert big bed to bed
			tempTxtFile = new File(Misc.removeExtension(ucscArchiveToConvert.getCanonicalPath())+".bed");
			String[] cmd = {bigBedToBed.getCanonicalPath(), ucscArchiveToConvert.getCanonicalPath(), tempTxtFile.getCanonicalPath()};
			executeUCSCCommand(cmd);

			//convert bed to useq
			//what kind of bed file is it?
			if (bedType(tempTxtFile) == false) throw new IOException("\nFailed to type the bed file, see -> "+ucscArchiveToConvert);
			
			//prepend chr?
			boolean prependChr = false;
			if (Text2USeq.chromosomesStartWithChr(tempTxtFile, 0) == false) prependChr = true;
			
			if (verbose){
				System.out.println(stranded+"\tstranded");
				System.out.println(scored+"\tscored");
				System.out.println(named+"\tnamed");
				System.out.println(bed12+"\tbed12");
				System.out.println(prependChr+"\tprepend 'chr'");
			}
			
			Text2USeq t2u = new Text2USeq(verbose);
			t2u.setVersionedGenome(versionedGenome);
			t2u.setChromosomeColumnIndex(0);
			t2u.setBeginningColumnIndex(1);
			t2u.setEndingColumnIndex(2);
			t2u.setConvertM(true);
			t2u.setPrependChr(prependChr);
			
			
			if (bed12){
				//-t 3,6,7,8,9,10,11 -v 4 -s 5
				t2u.setTextColumnIndexs(new int[] {3,6,7,8,9,10,11});
				t2u.setScoreColumnIndex(4);
				t2u.setStrandColumnIndex(5);
			}
			else {
				if (stranded) t2u.setStrandColumnIndex(5);
				if (scored) t2u.setScoreColumnIndex(4);
				if (named) t2u.setTextColumnIndexs(new int[]{3});
			}
			t2u.convert(tempTxtFile);
			
		}

		else throw new IOException("Unsupported UCSC big file format (e.g. xxx.bb or xxx.bw) -> "+ucscArchiveToConvert.getName());

		return convertedUSeqArchive;

	}

	public boolean bedType(File bedFile){
		stranded = false;
		scored = false;
		named = false;
		bed12 = false;

		boolean bedLineFound = false;

		try {

			//read through first 10000 lines
			//chrom start stop name score strand ....
			//  0     1     2   3     4     5
			BufferedReader in = IO.fetchBufferedReader(bedFile);
			String[] tokens;
			String line;
			int counter = 0;
			while ((line=in.readLine())!= null){
				if (counter++ > 10000) break;
				//bed line?
				if (bedLine.matcher(line).matches()==false) continue;
				bedLineFound = true;
				//split on tab
				tokens = tab.split(line);
				//bed12?
				if (tokens.length == 12){
					bed12 = true;
					break;
				}
				if (tokens.length > 3){
					//named?
					if (badName.matcher(tokens[3]).matches() == false) named = true;
					if (tokens.length > 4){
						//scored?
						if (tokens[4].equals("0") == false) scored = true;
						if (tokens.length > 5){
							//stranded?
							if (tokens[5].equals(".") == false) stranded = true;
						}
					}
				}

			}

			in.close();

		} catch (Exception e){
			System.err.println("\nProblem parsing bed type.\n");
			e.printStackTrace();
			bedLineFound = false;
		}

		return bedLineFound;
	}

	public void run(){
		try {
			convert();
		} catch (Exception e) {
			// TODO Auto-generated catch block

			e.printStackTrace();
		} finally {
			if (deleteTempFiles) tempTxtFile.delete();
		}
	}


	private void executeUCSCCommand(String[] command) throws IOException{
		/*if (verbose) {
			System.out.println("\nUnix Command:");
			for (String c : command) System.out.println(c);
			System.out.println();
		}*/
		//execute ucsc converter, nothing should come back for wigToBigWig and sort
		String[] results = USeqUtilities.executeCommandLineReturnAll(command);
		if (results.length !=0){
			//scan to see if just bedToBigBed normal output
			boolean ok = true;
			StringBuilder sb = new StringBuilder("Error message:");
			for (String c : results) {
				sb.append("\n");
				sb.append(c);
				if (c.contains("millis") == false) ok = false;
			}
			if (ok != true) {
				throw new IOException (sb.toString());
			}
		}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new UCSCBig2USeq(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		File ucscDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': forExtraction = new File(args[++i]); break;
					case 'd': ucscDir = new File (args[++i]); break;
					case 'v': versionedGenome = args[++i]; break;
					case 'e': verbose = false; break;
					case 'f': forceConversion = true; break;
					case 'h': printDocs(); System.exit(0); break;					
					default: USeqUtilities.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					USeqUtilities.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+USeqUtilities.stringArrayToString(args, " ")+"\n");

		//versioned genome?
		if (versionedGenome == null) USeqUtilities.printExit("\nError: you must supply a genome version. Goto http://genome.ucsc.edu/cgi-" +
		"bin/hgGateway load your organism to find the associated genome version (e.g. H_sapiens_Mar_2006, H_sapiens_Feb_2009).\n");

		//make files
		if (ucscDir == null || ucscDir.isDirectory() == false) USeqUtilities.printExit("\nCannot find your directory containing the UCSC wig2BigWig and bed2BigBed apps -> "+ucscDir);
		bigWigToWig = new File( ucscDir, "bigWigToWig");
		bigBedToBed = new File( ucscDir, "bigBedToBed");

		//check files
		//if (bigWigToBedGraph.canExecute() == false) USeqUtilities.printExit("\nCannot find or execute -> "+bigWigToBedGraph+"\n");
		if (bigWigToWig.canExecute() == false) USeqUtilities.printExit("\nCannot find or execute -> "+bigWigToWig+"\n");
		if (bigBedToBed.canExecute() == false) USeqUtilities.printExit("\nCannot find or execute -> "+bigBedToBed+"\n");

		//pull ucsc files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printExit("\nError: please enter a bw or bb file or directory containing such to convert!\n");
		File[][] tot = new File[2][];
		tot[0] = USeqUtilities.fetchFilesRecursively(forExtraction, ".bw");
		tot[1] = USeqUtilities.fetchFilesRecursively(forExtraction, ".bb");

		ucscArchives = IO.collapseFileArray(tot);
		if (ucscArchives == null || ucscArchives.length == 0 || ucscArchives[0].canRead() == false) Misc.printExit("\nError: cannot find or read any xxx.bb or xxx.bw file(s)!\n");

	}	
	
	/**Removes and big files that have a converted xxx.useq file in the same directory.*/
	public static File[] removeExistingConvertedBigFiles (File[] bigFiles){
		ArrayList<File> toReturn = new ArrayList<File>();
		for (File big: bigFiles){
			//remove .bb or .bw extension
			String name = big.getName().substring(0, big.getName().length()-3);
			//remove _Minus or _Plus if then exist
			if (name.endsWith("_Minus")) name = name.substring(0, name.length()-6);
			else if (name.endsWith("_Plus")) name = name.substring(0, name.length()-5);
			File f = new File (big.getParentFile(), name +".useq");
			if (f.exists()) continue;
			toReturn.add(big);
		}
		File[] f = new File[toReturn.size()];
		toReturn.toArray(f);
		return f;
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              UCSC Big 2 USeq: Jan 2013                           **\n" +
				"**************************************************************************************\n" +
				"Converts UCSC bigWig (xxx.bw) or bigBed (xxx.bb) archives to xxx.useq archives.\n" +

				"\nOptions:\n"+
				"-b Full path file/directory containing xxx.bw and xxx.bb files. Recurses through sub \n" +
				"       if a directory is given.\n" +
				"-d Full path directory containing the UCSC bigWigToBedGraph, bigWigToWig, and \n" +
				"       bigBedToBed apps, download from http://hgdownload.cse.ucsc.edu/admin/exe/ and\n" +
				"       make executable (e.g. chmod 755 /MyApps/UCSC/*).\n"+
				"-v Genome version (e.g. H_sapiens_Mar_2006), get from UCSC Browser,\n" +
				"      http://genome.ucsc.edu/FAQ/FAQreleases or IGB \n" +
				"      http://bioviz.org/igb/releases/current/igb-large.jnlp\n"+
				"-f Force conversion of xxx.bw or xxx.bb overwriting any existing xxx.useq archives.\n"+
				"       Defaults to skipping those already converted.\n"+
				"-e Only print error messages.\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/USeq2UCSCBig -v M_musculus_Jul_2007 -u\n" +
				"      /AnalysisResults/USeqDataArchives/ -d /MyApps/UCSC/\n\n" +

		"**************************************************************************************\n");

	}
}
