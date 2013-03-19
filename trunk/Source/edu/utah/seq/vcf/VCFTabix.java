package edu.utah.seq.vcf;

import java.io.*;
import java.util.regex.*;
import util.gen.*;

/**Converts vcf into a gzipped tabix vcf file.  Bit redundant but useful for calling from GNomEx.
 * @author david.nix@hci.utah.edu 
 **/
public class VCFTabix{
	//fields
	File[] vcfFiles2Convert;
	File bgzip;
	File tabix;
	boolean overWriteExisting = false;
	boolean deleteNonGzippedVCF = true;
	boolean verbose = true;

	//stand alone constructor
	public VCFTabix(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);

			convert(vcfFiles2Convert);

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	//GNomEx constructor
	public VCFTabix (File bgzip, File tabix) throws IOException{
		this.bgzip = bgzip;
		this.tabix = tabix;
		verbose = false;
		if (bgzip.canExecute() == false || tabix.canExecute() == false) throw new IOException ("\nCannot find or execute bgzip or tabix see -> "+bgzip+" "+tabix);
	}

	public void convert(File[] filesToConvert) throws IOException{
		
		//for each vcf file
		if (verbose) System.out.println("Converting:");
		for (File vcf: filesToConvert){
			if (verbose) System.out.println("\t"+vcf);
			File parentDir = vcf.getParentFile();
			
			//check to see if it exists?
			if (overWriteExisting == false){
				String name = vcf.getName();
				boolean compExists = true;
				//if the current vcf isn't gzipped, look for a gzipped file
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
					if (deleteNonGzippedVCF && vcf.getName().endsWith(".gz") == false) vcf.delete();
					continue;
				}
			}
			
			//copy the file, uncompressing if needed
			File copy = new File (parentDir, Misc.removeExtension(vcf.getName())+"_tempCon.vcf");
			if (vcf.getName().endsWith(".zip") || vcf.getName().endsWith(".gz")){
				copy = IO.uncompress(vcf, copy);
				if (copy == null) throw new IOException("\nFailed to uncompress "+vcf);
			}
			else {
				if (IO.copyViaReadLine(vcf, copy) == false) {
					copy.delete();
					throw new IOException("\nFailed to copy file "+vcf);
				}
			}

			//force compress with bgzip, this will replace the copy
			String[] cmd = { bgzip.getCanonicalPath(), "-f", copy.getCanonicalPath()};
			String[] output = IO.executeCommandLineReturnAll(cmd);
			copy.delete();
			File compVCF = new File (parentDir, copy.getName()+".gz");
			if (output == null || output.length != 0 || compVCF.exists() == false){
				compVCF.delete();
				throw new IOException("\nFailed to bgzip compress vcf file "+vcf+" Error: "+Misc.stringArrayToString(output, "\n"));
			}

			//tabix
			cmd = new String[]{ tabix.getCanonicalPath(), "-f", "-p", "vcf", compVCF.getCanonicalPath() };
			output = IO.executeCommandLineReturnAll(cmd);
			File indexVCF = new File (parentDir, compVCF.getName()+".tbi");
			if (output == null || output.length != 0 || indexVCF.exists() == false){
				compVCF.delete();
				indexVCF.delete();
				throw new IOException("\nFailed to tabix index vcf file "+vcf+" Error: "+Misc.stringArrayToString(output, "\n"));
			}

			//all looks good so rename
			File newVCF = new File (parentDir, compVCF.getName().replace("_tempCon.vcf.gz", ".vcf.gz"));
			if (newVCF.exists()) newVCF.delete();
			if (compVCF.renameTo(newVCF) == false) {
				compVCF.delete();
				indexVCF.delete();
				newVCF.delete();
				throw new IOException("\nFailed to rename compressed VCF file "+newVCF+"\n");
			}
			File newIndex = new File (parentDir, newVCF.getName()+".tbi");
			if (newIndex.exists()) newIndex.delete();
			if (indexVCF.renameTo(newIndex) == false) {
				compVCF.delete();
				indexVCF.delete();
				newVCF.delete();
				newIndex.delete();
				throw new IOException("\nFailed to rename VCF index file "+newIndex+"\n");
			}
			
			//delete uncompressed?
			if (deleteNonGzippedVCF && vcf.getName().endsWith(".gz") == false) vcf.delete();

		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFTabix(args);
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
					case 'd': deleteNonGzippedVCF = false; break;
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

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.fetchFilesRecursively(forExtraction, ".vcf");
		tot[1] = IO.fetchFilesRecursively(forExtraction,".vcf.gz");
		tot[2] = IO.fetchFilesRecursively(forExtraction,".vcf.zip");
		vcfFiles2Convert = IO.collapseFileArray(tot);
		if (vcfFiles2Convert == null || vcfFiles2Convert.length ==0 || vcfFiles2Convert[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		//pull executables
		if (tabixBinDirectory == null) Misc.printExit("\nError: please point to your tabix executable directory (e.g. /Users/madonna/Desktop/VCF/tabix-0.2.6/ )\n");
		bgzip = new File (tabixBinDirectory, "bgzip");
		tabix = new File (tabixBinDirectory, "tabix");
		//look for bgzip and tabix executables
		if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+bgzip+" "+tabix);


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                VCFTabix: Jan 2013                                **\n" +
				"**************************************************************************************\n" +
				"Converts vcf files to a SAMTools compressed vcf tabix format. Recursive.\n" +

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s). Recursive!\n" +
				"-t Full path tabix directory containing the compiled bgzip and tabix executables. See\n" +
				"      http://sourceforge.net/projects/samtools/files/tabix/\n"+
				"-f Force overwriting of existing indexed vcf files, defaults to skipping.\n"+
				"-d Do not delete non gzipped vcf files after successful indexing, defaults to deleting.\n"+
				"-e Only print error messages.\n"+

				"\nExample: java -jar pathToUSeq/Apps/VCFTabix -v /VarScan2/VCFFiles/\n" +
				"     -t /Samtools/Tabix/tabix-0.2.6/ \n\n" +

		"**************************************************************************************\n");

	}



}
