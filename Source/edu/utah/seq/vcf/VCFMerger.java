package edu.utah.seq.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class VCFMerger {

	private File[] vcfFiles;
	private VCFParser[] vcfParsers;
	private Gzipper out;
	private int fileIndexOfConsensus = 1;

	public VCFMerger(String[] args){
		try {	
			processArgs(args);

			//load each file
			vcfParsers = new VCFParser[vcfFiles.length];
			int num = 0;
			for (int i=0; i< vcfFiles.length; i++){
				System.out.println(i+ "\t"+ vcfFiles[i].getName());
				vcfParsers[i] = new VCFParser(vcfFiles[i], true, true, true);
				vcfParsers[i].appendChr();
				vcfParsers[i].setFilterFieldPeriodToTextOnAllRecords(VCFRecord.PASS);
				vcfParsers[i].filterVCFRecords(VCFRecord.PASS);
				num += vcfParsers[i].getVcfRecords().length;
				//set score to parser index
				VCFRecord[] v = vcfParsers[i].getVcfRecords();
				for (int j=0; j< v.length; j++) v[j].setScore(i);
			}


			out = new Gzipper(new File(vcfFiles[0].getParentFile(), "merged.vcf.gz"));


			//create array
			VCFRecord[] all = new VCFRecord[num];
			int index = 0;
			for (int i=0; i< vcfFiles.length; i++){
				VCFRecord[] v = vcfParsers[i].getVcfRecords();
				for (int j=0; j< v.length; j++) all[index++] = v[j];
			}

			//sort by chrom and position
			Arrays.sort(all);

			//filter
			String oldChrom = all[0].getChromosome();
			int oldPos = all[0].getPosition();
			ArrayList<VCFRecord> records = new ArrayList<VCFRecord>();
			records.add(all[0]);
			VCFRecord oldRecord = all[0];
			for (int i=1; i< all.length; i++){
				String currChrom = all[i].getChromosome();
				int currPos = all[i].getPosition();
				//diff chrom or position or alternate allele?
				if (currChrom.equals(oldChrom) == false || currPos != oldPos  || all[i].matchesAlternateAlleleGenotype(oldRecord, false) == false){
					processRecords(records);
				}
				records.add(all[i]);
				oldChrom = currChrom;
				oldPos = currPos;
				oldRecord = all[i];
			}
			//print last
			processRecords(records);

			out.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private void processRecords(ArrayList<VCFRecord> records) throws IOException {
		int num = records.size();

		//only one?
		if (num == 1) out.println(records.get(0));

		else {
			//is a consensus present?
			boolean consensusFound = false;
			for (int i=0; i< num; i++){
				if (records.get(i).getScore() == fileIndexOfConsensus){
					out.println(records.get(i));
					consensusFound = true;
					break;
				}
			}
			//no consensus print first
			if (consensusFound == false) out.println(records.get(0));
		}
		records.clear();

	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFMerger(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please indicate a vcf file to filter.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Comparator : March 2013                          **\n" +
				"**************************************************************************************\n" +
				"Not for distribution.\n\n" +

				"Required Options:\n"+
				"-v VCF files \n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/VCFComparator\n\n"+

		"**************************************************************************************\n");

	}

}
