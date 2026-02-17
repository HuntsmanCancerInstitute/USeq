package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import util.gen.*;

/**Converts bed into a gzipped tabix bed file.  Bit redundant but useful for calling from GNomEx.
 * @author david.nix@hci.utah.edu 
 **/
public class BedTabix{
	//fields
	File[] bedFiles;
	File bgzip;
	File tabix;
	boolean overWriteExisting = false;
	boolean deleteNonGzippedBED = true;
	boolean verbose = true;

	//stand alone constructor
	public BedTabix(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);

			convert(bedFiles);

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: failed to tabix bed files.");
		}
	}
	//GNomEx constructor
	public BedTabix (File bgzip, File tabix) throws IOException{
		this.bgzip = bgzip;
		this.tabix = tabix;
		verbose = false;
		if (bgzip.canExecute() == false || tabix.canExecute() == false) throw new IOException ("\nCannot find or execute bgzip or tabix see -> "+bgzip+" "+tabix);
	}

	public void convert(File[] filesToConvert) throws IOException{
		
		//for each bed file
		if (verbose) System.out.println("Converting:");
		for (File bed: filesToConvert){
			if (verbose) System.out.println("\t"+bed);
			File parentDir = bed.getParentFile();
			
			//check to see if it exists?
			if (overWriteExisting == false){
				String name = bed.getName();
				boolean compExists = true;
				//if the current bed isn't gzipped, look for a gzipped file
				if (name.endsWith(".gz") == false) {
					//does the gzipped version exist?
					File existingComp = new File (parentDir, name+".gz");
					if (existingComp.exists() == false) compExists = false;
					name = name+".gz.tbi";
				}
				else name = name + ".tbi";
				
				//does index already exist?
				File existingIndex = new File (parentDir, name);
				if (existingIndex.exists() && compExists) {
					if (verbose) System.out.println("\t\tIndexed files exist, skipping.");
					if (deleteNonGzippedBED && bed.getName().endsWith(".gz") == false) bed.delete();
					continue;
				}
			}
			
			//copy the file, uncompressing if needed
			File copy = new File (parentDir, Misc.removeExtension(bed.getName())+"_tempCon.bed");
			if (bed.getName().endsWith(".zip") || bed.getName().endsWith(".gz")){
				copy = IO.uncompress(bed, copy);
				if (copy == null) throw new IOException("\nFailed to uncompress "+bed);
			}
			else {
				if (IO.copyViaReadLine(bed, copy) == false) {
					copy.delete();
					throw new IOException("\nFailed to copy file "+bed);
				}
			}

			//force compress with bgzip, this will replace the copy
			String[] cmd = { bgzip.getCanonicalPath(), "-f", copy.getCanonicalPath()};
			String[] output = IO.executeCommandLineReturnAll(cmd);
			copy.delete();
			File compBED = new File (parentDir, copy.getName()+".gz");
			if (output == null || output.length != 0 || compBED.exists() == false){
				compBED.delete();
				throw new IOException("\nFailed to bgzip compress bed file "+bed+" Error: "+Misc.stringArrayToString(output, "\n"));
			}

			//tabix
			//must use -0 --sequence 1 --begin 2 --end 3, -p bed doesn't work with java tabix!!!!
			cmd = new String[]{ tabix.getCanonicalPath(), "-0", "--sequence", "1", "--begin", "2", "--end", "3", compBED.getCanonicalPath() };
			output = IO.executeCommandLineReturnAll(cmd);
			File indexBED = new File (parentDir, compBED.getName()+".tbi");
			if (output == null || output.length != 0 || indexBED.exists() == false){
				compBED.delete();
				indexBED.delete();
				throw new IOException("\nFailed to tabix index bed file "+bed+" Error: "+Misc.stringArrayToString(output, "\n"));
			}

			//all looks good so rename
			File newBED = new File (parentDir, compBED.getName().replace("_tempCon.bed.gz", ".bed.gz"));
			if (newBED.exists()) newBED.delete();
			if (compBED.renameTo(newBED) == false) {
				compBED.delete();
				indexBED.delete();
				newBED.delete();
				throw new IOException("\nFailed to rename compressed BED file "+newBED+"\n");
			}
			File newIndex = new File (parentDir, newBED.getName()+".tbi");
			if (newIndex.exists()) newIndex.delete();
			if (indexBED.renameTo(newIndex) == false) {
				compBED.delete();
				indexBED.delete();
				newBED.delete();
				newIndex.delete();
				throw new IOException("\nFailed to rename BED index file "+newIndex+"\n");
			}
			
			//delete uncompressed?
			if (deleteNonGzippedBED && bed.getName().endsWith(".gz") == false) bed.delete();

		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BedTabix(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		File tabixBinDirectory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
					case 'f': overWriteExisting = true; break;
					case 'd': deleteNonGzippedBED = false; break;
					case 'e': verbose = false; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

		//pull bed files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a bed file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.fetchFilesRecursively(forExtraction, ".bed");
		tot[1] = IO.fetchFilesRecursively(forExtraction,".bed.gz");
		tot[2] = IO.fetchFilesRecursively(forExtraction,".bed.zip");
		bedFiles = IO.collapseFileArray(tot);
		if (bedFiles == null || bedFiles.length ==0 || bedFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bed(.zip/.gz OK) file(s)!\n");

		//pull executables
		if (tabixBinDirectory == null) Misc.printExit("\nError: please point to your tabix executable directory (e.g. /Users/madonna/Desktop/BED/tabix-0.2.6/ )\n");
		bgzip = new File (tabixBinDirectory, "bgzip");
		tabix = new File (tabixBinDirectory, "tabix");
		//look for bgzip and tabix executables
		if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+bgzip+" "+tabix);


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                BedTabix: Nov 2025                                **\n" +
				"**************************************************************************************\n" +
				"Converts bed files to a SAMTools compressed bed tabix format. Recursive.\n" +

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.bed(.gz/.zip OK) file(s). Recursive!\n" +
				"-t Full path tabix directory containing the compiled bgzip and tabix executables. See\n" +
				"      http://sourceforge.net/projects/samtools/files/tabix/\n"+
				"-f Force overwriting of existing indexed bed files, defaults to skipping.\n"+
				"-d Do not delete non gzipped bed files after successful indexing, defaults to deleting.\n"+
				"-e Only print error messages.\n"+

				"\nExample: java -jar pathToUSeq/Apps/BedTabix -v /VarScan2/BEDFiles/\n" +
				"     -t /Samtools/Tabix/tabix-0.2.6/ \n\n" +

		"**************************************************************************************\n");

	}



}
