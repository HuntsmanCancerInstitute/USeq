package edu.utah.seq.useq.apps;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;
import edu.utah.seq.useq.USeqArchive;
import edu.utah.seq.useq.USeqUtilities;

/**Class to convert xxx.useq archives to UCSC xxx.bb or xxx.bw archives.*/
public class USeq2UCSCBig extends Thread{

	private File[] useqArchives;
	private File workingUSeqArchiveFile;
	private USeqArchive workingUSeqArchive;
	private File ucscWig2BigWig;
	private File ucscBed2BigBed;
	private boolean verbose = true;
	private int lengthExtender = 10000;
	private File chromLengths;
	private File convertedFile;
	private File tempFile;
	private File tempFileSorted;
	private boolean deleteTempFiles = false;
	private boolean forceConversion = false;

	//constructors
	//stand alone
	public USeq2UCSCBig (String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		if (verbose) System.out.println("Processing:");
		
		//remove those that already exist?
		if (forceConversion == false) {
			useqArchives = USeq2UCSCBig.removeExistingConvertedUSeqArchives (useqArchives);
			if (useqArchives.length == 0) {
				if (verbose) System.out.println("\tNo unconverted xxx.useq archives were found.  Use the -f option to overwrite.\n");
				System.exit(0);
			}
		}
		
		//for each archive
		for (File u : useqArchives){
			workingUSeqArchiveFile = u;
			if (verbose) System.out.println("\t"+workingUSeqArchiveFile);
			convert();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

	}

	//for GenoPub/ GNomEx integration using threads, create the Object, call the fetchConvertedFileNames(), then start the thread
	public USeq2UCSCBig (File ucscWig2BigWig, File ucscBed2BigBed, File useq){
		this.ucscWig2BigWig = ucscWig2BigWig;
		this.ucscBed2BigBed = ucscBed2BigBed;
		workingUSeqArchiveFile = useq;
		verbose = false;
	}

	//methods
	/**Returns converted File(s) or null if something bad happened.*/
	public ArrayList<File> convert (){
		try {
			//create archive
			workingUSeqArchive = new USeqArchive(workingUSeqArchiveFile);

			//write out chrom name: max lengths
			writeChromLengths();
			//convert graph data
			if (workingUSeqArchive.getArchiveInfo().isGraphData()) {
				return convertGraphData();
			}

			//convert region data
			else {
				return convertRegionData();
			}

		} catch (Exception e){
			e.printStackTrace();
		}
		return null;

	}

	public void run(){
		convert();
	}

	/**Fetches File objects for the converted files if they were converted.  Does not convert.*/
	public ArrayList<File> fetchConvertedFileNames(){
		//create archive
		try {
			workingUSeqArchive = new USeqArchive(workingUSeqArchiveFile);
		} catch (Exception e) {
			e.printStackTrace();
		}

		//convert graph data
		if (workingUSeqArchive.getArchiveInfo().isGraphData()) return fetchConvertedGraphNames();

		//convert region data
		else return fetchConvertRegionNames();
	}

	private ArrayList<File> fetchConvertRegionNames(){
		String name = workingUSeqArchiveFile.getName().replace(USeqUtilities.USEQ_EXTENSION_WITH_PERIOD, "");
		File convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + ".bb");
		ArrayList<File> al = new ArrayList<File>();
		al.add(convertedFile);
		return al;
	}
	
	public static boolean dataPresent(File data){
		try {
			BufferedReader in = IO.fetchBufferedReader(data);
			String line;
			while ((line = in.readLine()) !=null){
				if (line.startsWith("#") == false) return true;
			}
			in.close();
		} catch (IOException e) {
			System.err.println("Failed to make a buffered reader on data file: "+data);
			e.printStackTrace();
		}
		return false;
	}

	private ArrayList<File> convertRegionData() throws Exception{
		String name = workingUSeqArchiveFile.getName().replace(USeqUtilities.USEQ_EXTENSION_WITH_PERIOD, "");
		tempFile = new File (workingUSeqArchiveFile.getCanonicalPath() + ".bed");

		USeq2Text.print2TextFile(workingUSeqArchiveFile, tempFile, true, true);
		
		//any data written to tempFile?  this is a catch for bad bed12 files
		if (dataPresent(tempFile) == false){
			deleteAllFiles();
			if (verbose) System.err.println("\t\tNo data written to bed file, skipping!");
			return null;
		}

		//sort it using unix command
		tempFileSorted = new File (workingUSeqArchiveFile.getCanonicalPath() + ".sorted.bed");
		//String shell = "sort -k1,1 -k2,2n "+tempFile.getCanonicalPath()+" > "+tempFileSorted.getCanonicalPath();
		//String[] output = USeqUtilities.executeShellScript(shell, workingUSeqArchiveFile.getParentFile());
		String[] cmd = {"sort", "-k1,1", "-k2,2n", "-o", tempFileSorted.getCanonicalPath(), tempFile.getCanonicalPath(), };
		String[] output = IO.executeCommandLineReturnAll(cmd);
		if (output == null || output.length !=0 || tempFileSorted.exists() == false){
			throw new IOException("Error sorting bed file: "+tempFileSorted+"\nError: "+Misc.stringArrayToString(output, "\n"));
		}
		
		//convert to binary
		convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + ".bb");
		String[] command = new String[]{
				ucscBed2BigBed.getCanonicalPath(),
				tempFileSorted.getCanonicalPath(),
				chromLengths.getCanonicalPath(),
				convertedFile.getCanonicalPath()
		};
		executeUCSCCommand(command);
		deleteTempFiles();
		//return
		ArrayList<File> al = new ArrayList<File>();
		al.add(convertedFile);
		return al;
	}

	private ArrayList<File> fetchConvertedGraphNames(){
		String name = workingUSeqArchiveFile.getName().replace(USeqUtilities.USEQ_EXTENSION_WITH_PERIOD, "");
		ArrayList<File> convertedFiles = new ArrayList<File>();
		//is it stranded
		boolean stranded = workingUSeqArchive.isStranded();
		if (stranded){
			File convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + "_Plus.bw");
			convertedFiles.add(convertedFile);
			convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + "_Minus.bw");
			convertedFiles.add(convertedFile);
		}
		else {
			File convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + ".bw");
			convertedFiles.add(convertedFile);
		}
		return convertedFiles;
	}


	private ArrayList<File> convertGraphData() throws Exception{	
		
		USeq2Text useq2Text = new USeq2Text();
		useq2Text.setPrintWigFormat(true);
		String name = workingUSeqArchiveFile.getName().replace(USeqUtilities.USEQ_EXTENSION_WITH_PERIOD, "");
		ArrayList<File> convertedFiles = new ArrayList<File>();
		tempFile = new File (workingUSeqArchiveFile.getCanonicalPath() + ".wig");		

		//is it stranded
		boolean stranded = workingUSeqArchive.isStranded();
		if (stranded){

			//convert for plus?
			if (workingUSeqArchive.isPlusStrandedPresent()){				
				useq2Text.print2WigFile(workingUSeqArchiveFile, tempFile, "+");

				//convert text to binary, wigToBigWig in.wig chrom.sizes out.bw
				convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + "_Plus.bw");
				String[] command = new String[]{
						ucscWig2BigWig.getCanonicalPath(),
						tempFile.getCanonicalPath(),
						chromLengths.getCanonicalPath(),
						convertedFile.getCanonicalPath()
				};

				//execute it
				executeUCSCCommand(command);
				convertedFiles.add(convertedFile);

			}


			//convert minus?
			if (workingUSeqArchive.isMinusStrandedPresent()){
				useq2Text.print2WigFile(workingUSeqArchiveFile, tempFile, "-");

				//convert text to binary, wigToBigWig in.wig chrom.sizes out.bw
				convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + "_Minus.bw");
				String[] command = new String[]{
						ucscWig2BigWig.getCanonicalPath(),
						tempFile.getCanonicalPath(),
						chromLengths.getCanonicalPath(),
						convertedFile.getCanonicalPath()
				};
				//execute and save it
				executeUCSCCommand(command);
				convertedFiles.add(convertedFile);
			}

		}
		else {
			//convert all		
			
			useq2Text.print2WigFile(workingUSeqArchiveFile, tempFile, null);
			//convert text to binary, wigToBigWig in.wig chrom.sizes out.bw
			convertedFile = new File (workingUSeqArchiveFile.getParentFile(), name + ".bw");
			String[] command = new String[]{
					ucscWig2BigWig.getCanonicalPath(),
					tempFile.getCanonicalPath(),
					chromLengths.getCanonicalPath(),
					convertedFile.getCanonicalPath()
			};
			//execute it
			executeUCSCCommand(command);

			convertedFiles.add(convertedFile);
		}
		//cleanup
		deleteTempFiles();

		return convertedFiles;
	}

	private void executeUCSCCommand(String[] command) throws Exception{
		
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
				deleteAllFiles();
				throw new Exception (sb.toString());
			}
		}
	}

	private void writeChromLengths() throws IOException{

		HashMap<String,Integer> nameBase = workingUSeqArchive.fetchChromosomesAndLastBase();
		chromLengths = new File (workingUSeqArchive.getZipFile()+".chromLengths");		
		PrintWriter out = new PrintWriter( new FileWriter (chromLengths));
		for (String name: nameBase.keySet()){
			int length = nameBase.get(name)+ lengthExtender;
			out.print(name);
			out.print("\t");
			out.println(length);
		}
		out.close();
	}

	public void deleteAllFiles(){
		deleteTempFiles();
		if (convertedFile != null) convertedFile.delete();
	}
	public void deleteTempFiles(){	
		if (deleteTempFiles){
			if (tempFile!= null) tempFile.delete();
			if (chromLengths!= null) chromLengths.delete();
			if (tempFileSorted!= null) tempFileSorted.delete();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new USeq2UCSCBig(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File ucscDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'u': useqArchives = USeqUtilities.fetchFilesRecursively(new File(args[++i]), USeqUtilities.USEQ_EXTENSION_WITH_PERIOD); break;
					case 'd': ucscDir = new File (args[++i]); break;
					case 'f': forceConversion = true; break;
					case 'e': verbose = false; break;
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
		
		//make files
		if (ucscDir == null || ucscDir.isDirectory() == false) USeqUtilities.printExit("\nCannot find your directory containing the UCSC wig2BigWig and bed2BigBed apps -> "+ucscDir);
		ucscWig2BigWig = new File( ucscDir, "wigToBigWig");
		ucscBed2BigBed = new File( ucscDir, "bedToBigBed");

		//check files
		if (useqArchives == null || useqArchives.length == 0) USeqUtilities.printExit("\nCannot find any xxx."+USeqUtilities.USEQ_EXTENSION_NO_PERIOD+" USeq archives?\n");
		if (ucscWig2BigWig.canExecute() == false) USeqUtilities.printExit("\nCannot find or execute -> "+ucscWig2BigWig+"\n");
		if (ucscBed2BigBed.canExecute() == false) USeqUtilities.printExit("\nCannot find or execute -> "+ucscBed2BigBed+"\n");
		

	}	

	public static File[] removeExistingConvertedUSeqArchives (File[] useqFiles){
		ArrayList<File> toReturn = new ArrayList<File>();
		for (File useq: useqFiles){
			String name = useq.getName().substring(0, useq.getName().length()-5);
			File f = new File (useq.getParentFile(), name +".bw");
			if (f.exists()) continue;
			f = new File (useq.getParentFile(), name +"_Minus.bw");
			if (f.exists()) continue;
			f = new File (useq.getParentFile(), name +"_Plus.bw");
			if (f.exists()) continue;
			f = new File (useq.getParentFile(), name +".bb");
			if (f.exists()) continue;
			f = new File (useq.getParentFile(), name +"_Minus.bb");
			if (f.exists()) continue;
			f = new File (useq.getParentFile(), name +"_Plus.bb");
			if (f.exists()) continue;
			toReturn.add(useq);
		}
		File[] f = new File[toReturn.size()];
		toReturn.toArray(f);
		return f;
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              USeq 2 UCSC Big: Jan 2013                           **\n" +
				"**************************************************************************************\n" +
				"Converts USeq archives to UCSC bigWig (xxx.bw) or bigBed (xxx.bb) archives based on\n" +
				"the data type. WARNING: bigBed format conversion will clip any associated scores to\n" +
				"between 0-1000. \n" +

				"\nOptions:\n"+
				"-u Full path file/directory containing xxx.useq files. Recurses through sub \n" +
				"       if a directory is given.\n" +
				"-d Full path directory containing the UCSC wigToBigWig and bedToBigBed apps, download\n" +
				"       from http://hgdownload.cse.ucsc.edu/admin/exe/ and make executable with chmod.\n"+
				"-f Force conversion of xxx.useq to xxx.bw or xxx.bb overwriting any UCSC big files.\n"+
				"       Defaults to skipping those already converted.\n"+
				"-e Only print error messages.\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/USeq2UCSCBig -u\n" +
				"      /AnalysisResults/USeqDataArchives/ -d /Apps/UCSC/\n\n" +

		"**************************************************************************************\n");

	}
}
