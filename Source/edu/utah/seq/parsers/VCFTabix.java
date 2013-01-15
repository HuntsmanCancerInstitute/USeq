package edu.utah.seq.parsers;

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
			//copy the file, uncompressing if needed
			File copy = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"_tempCon.vcf");
			if (vcf.getName().endsWith(".zip") || vcf.getName().endsWith(".gz")){
				copy = IO.uncompress(vcf, copy);
				if (copy == null) throw new IOException("\nFailed to uncompress "+vcf);
			}
			else {
				if (IO.copy(vcf, copy) == false) throw new IOException("\nFailed to copy file "+vcf);
			}

			//force compress with bgzip, this will replace the copy
			String[] cmd = { bgzip.getCanonicalPath(), "-f", copy.getCanonicalPath()};
			String[] output = IO.executeCommandLineReturnAll(cmd);
			copy.delete();
			File compVCF = new File (copy.getParentFile(), copy.getName()+".gz");
			if (output == null || output.length != 0 || compVCF.exists() == false){
				compVCF.delete();
				throw new IOException("\nFailed to bgzip compress vcf file "+vcf+" Error: "+Misc.stringArrayToString(output, "\n"));
			}

			//tabix
			cmd = new String[]{ tabix.getCanonicalPath(), "-f", "-p", "vcf", compVCF.getCanonicalPath() };
			output = IO.executeCommandLineReturnAll(cmd);
			File indexVCF = new File (copy.getParentFile(), compVCF.getName()+".tbi");
			if (output == null || output.length != 0 || indexVCF.exists() == false){
				compVCF.delete();
				indexVCF.delete();
				throw new IOException("\nFailed to tabix index vcf file "+vcf+" Error: "+Misc.stringArrayToString(output, "\n"));
			}

			//all looks good so rename
			File newVCF = new File (compVCF.getParentFile(), compVCF.getName().replace("_tempCon.vcf.gz", ".vcf.gz"));
			if (newVCF.exists()) newVCF.delete();
			if (compVCF.renameTo(newVCF) == false) {
				compVCF.delete();
				indexVCF.delete();
				newVCF.delete();
				throw new IOException("\nFailed to rename compressed VCF file "+newVCF+"\n");
			}
			File newIndex = new File (newVCF.getParentFile(), newVCF.getName()+".tbi");
			if (newIndex.exists()) newIndex.delete();
			if (indexVCF.renameTo(newIndex) == false) {
				compVCF.delete();
				indexVCF.delete();
				newVCF.delete();
				newIndex.delete();
				throw new IOException("\nFailed to rename VCF index file "+newIndex+"\n");
			}

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
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
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
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//pull vcf files
		File[][] tot = new File[3][];
		tot[0] = IO.fetchFilesRecursively(forExtraction, ".vcf");
		tot[1] = IO.fetchFilesRecursively(forExtraction,".vcf.gz");
		tot[2] = IO.fetchFilesRecursively(forExtraction,".vcf.zip");
		vcfFiles2Convert = IO.collapseFileArray(tot);
		if (vcfFiles2Convert == null || vcfFiles2Convert.length ==0 || vcfFiles2Convert[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");

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
				"Converts vcf files to a SAMTools compressed vcf tabix format. \n" +

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s). Recursive!\n" +
				"-t Full path tabix directory containing the compiled bgzip and tabix executables. See\n" +
				"      http://sourceforge.net/projects/samtools/files/tabix/\n"+

				"\nExample: java -jar pathToUSeq/Apps/VCFTabix -v /VarScan2/VCFFiles/\n" +
				"     -t /Samtools/Tabix/tabix-0.2.6/ \n\n" +

		"**************************************************************************************\n");

	}



}
