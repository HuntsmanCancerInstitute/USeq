package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Parses a two sample tumor normal vcf file to generate ADTEx B Allele Frequency tables.*/
public class AdetexBafVCFParser {

	private File[] vcfFiles;
	private int minCoverage = 20;
	private int tumorSampleIndex = 1;
	private int normalSampleIndex = 0;

	public AdetexBafVCFParser (String[] args) {

		processArgs(args);

		System.out.println("VcfFile\tTotalRecords\tPassingRecords\tFraction");
		for (File vcf: vcfFiles){
			parse(vcf);
		}
	}

	public void parse(File vcf){
		try {
			System.out.print(vcf.getName());
			VCFParser parser = new VCFParser (vcf, false, true, true);
			BufferedReader in = parser.initializeParser();
			
			String name = Misc.removeExtension(vcf.getName());
			Gzipper out = new Gzipper(new File (vcf.getParentFile(), name+".baf.gz"));
			out.println("chrom\tSNP_name\tSNP_loc\tcontrol_BAF\ttumor_BAF\tcontrol_doc\ttumor_doc");
			VCFRecord r = null;
			int numRecords = 0;
			int numPrinted = 0;
			while ((r= parser.fetchNext(in)) != null){
				VCFSample[] samples = r.getSample();
				numRecords++;
				//do both have adequate coverage?
				int depthT = samples[tumorSampleIndex].getReadDepthDP();
				if (depthT < minCoverage) continue;
				int depthN = samples[normalSampleIndex].getReadDepthDP();
				if (depthN < minCoverage) continue;
				double depthAltT = Double.parseDouble(samples[tumorSampleIndex].getAlternateCounts());
				double depthAltN = Double.parseDouble(samples[normalSampleIndex].getAlternateCounts());
				double afT = depthAltT/(double)depthT;
				double afN = depthAltN/(double)depthN;
				/*
				chrom - chromosome name (same format as in BED or BAM file)
				SNP_loc - location of the SNP, but what coord system? 0 or 1 based?
				control_BAF - B allele frequency (BAF) at each SNP in control sample
				tumor_BAF - B allele frequency (BAF) at each SNP in tumor sample
				control_doc - Total read count at each SNP in control sample
				tumor_doc - Total read count at each SNP in tumor sample
				*/
				out.println(r.getChromosome()+"\tVar"+numRecords+"\t"+(r.getPosition()+1)+"\t"+afN+"\t"+afT+"\t"+depthN+"\t"+depthT);
				numPrinted++;
			}
			out.close();
			double frac = (double)numPrinted/ (double)numRecords;
			System.out.println("\t"+numRecords+"\t"+numPrinted+"\t"+frac);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing "+vcf);
		}
	}
	
	public void printRecordsSSC(VCFParser parser, File f, String passFail) throws Exception{
		PrintWriter out = new PrintWriter (new FileWriter (f));
		//write out header
		for (String h : parser.getStringComments()) out.println(h);
		VCFRecord[] records = parser.getVcfRecords();
		for (VCFRecord vcf: records){
			if (vcf.getFilter().equals(passFail) == false) continue;
			String orig = vcf.toString();
			String[] fields = VCFParser.TAB.split(orig);
			//reset score
			fields[5] = Integer.toString((int)vcf.getQuality());
			out.println(Misc.stringArrayToString(fields, "\t"));
		}
		out.close();
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AdetexBafVCFParser(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              ADTex VCFParser: Jan 2015                           **\n" +
				"**************************************************************************************\n" +
				"Parses VCF files containing Tumor and Normal samples to generate a ADTex B allele\n"+
				"frequency table.  VarScan2 outputs all variants in T or N regardless of somatic state.\n"+

				"\nRequired Options:\n"+
				"-v Full path file or directory containing xxx.vcf(.gz/.zip OK) file(s).\n" +
				"-m Minimum read coverage, defaults to 20\n"+
				"-n Normal sample index, defaults to 0\n"+
				"-t Tumor sample index, defaults to 1\n"+

				"\nExample: java -jar pathToUSeq/Apps/VarScanVCFParser -v /VarScan2/VCFFiles/\n\n" +
				"**************************************************************************************\n");

	}

}
