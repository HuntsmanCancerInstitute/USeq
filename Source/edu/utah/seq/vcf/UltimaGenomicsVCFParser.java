package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**UG parser. */
public class UltimaGenomicsVCFParser {

	private File[] vcfFiles;
	private double minimumTumorReadDepth = 10;
	private double minimumNormalReadDepth = 10;
	private double maximumNormalAltFraction = 0.1;
	private double minimumTumorAltFraction = 0.05;
	private double minimumAltReadDepth = 3;
	private boolean excludeNonPass = false;
	private boolean excludeNonChr = false;
	private File saveDirectory = null;
	
	private static String afInfo = "##INFO=<ID=T_AF,Number=1,Type=Float,Description=\"Allele Frequency for tumor\">";
	private static String dpInfo = "##INFO=<ID=T_DP,Number=1,Type=Integer,Description=\"Read depth for tumor\">";
	private static String nafInfo = "##INFO=<ID=N_AF,Number=1,Type=Float,Description=\"Allele Frequency for normal\">";
	private static String ndpInfo = "##INFO=<ID=N_DP,Number=1,Type=Integer,Description=\"Read depth for normal\">";
	private static String format = "GT:AD:BG_AD:BG_DP:BG_SB:BG_VAF:DP:GQ:SB:VAF:PL";

	public UltimaGenomicsVCFParser (String[] args) {

		processArgs(args);
		
		IO.pl("\nName\tPassing\tFailing");
		for (File vcf: vcfFiles){
			System.out.print(vcf.getName()+"\t");
			parse(vcf);
		}
		
		IO.pl();
	}

	public void parse(File vcf) {
		String line = null;
		try {
			//counters
			int numFail = 0;
			int numPass = 0;

			//IO
			BufferedReader in = IO.fetchBufferedReader(vcf);
			String name = Misc.removeExtension(vcf.getName());
			Gzipper out = new Gzipper (new File (saveDirectory, name+".filt.vcf.gz"));

			//for each line in the file
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.length()==0) continue;
				
				//header? just print out
				if (line.startsWith("#")) {
					if (line.startsWith("#CHROM")) {
						out.println(afInfo);
						out.println(dpInfo);
						out.println(nafInfo);
						out.println(ndpInfo);
					}
					out.println(line);
					continue;
				}

				//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT	TumorSampleName
				//   0    1   2  3   4   5     6      7     8     9   
				String[] tokens = Misc.TAB.split(line);

				if (excludeNonPass && tokens[6].toUpperCase().equals("PASS") == false) {
					numFail++;
					continue;
				}
				
				if (excludeNonChr && tokens[0].startsWith("chr") == false) {
					numFail++;
					continue;
				}
				
				//parse tumor and normal counts
				double[] tnRefAlt = parseCounts(tokens[8], tokens[9]);
				
				double tDP = tnRefAlt[0]+ tnRefAlt[1];
				if (tDP < minimumTumorReadDepth) {
					numFail++;
					continue;
				}
				double nDP = tnRefAlt[2]+ tnRefAlt[3];
				if (nDP < minimumNormalReadDepth) {
					numFail++;
					continue;
				}
				if (tnRefAlt[1] < minimumAltReadDepth) {
					numFail++;
					continue;
				}
				double tAF = tnRefAlt[1]/ tDP;
				if (tAF < minimumTumorAltFraction) {
					numFail++;
					continue;
				}
				double nAF = tnRefAlt[3]/ nDP;
				if (nAF > maximumNormalAltFraction) {
					numFail++;
					continue;
				}
				numPass++;
				
				//append T_AF T_DP N_AF N_DP
				String info = "T_AF="+Num.formatNumber(tAF, 5)+
						";N_AF="+Num.formatNumber(nAF, 5)+
						";T_DP="+(int)tDP+
						";N_DP="+(int)nDP;

				tokens[7] = info+";"+ tokens[7];
				
				//add id
				tokens[2] = "UG_"+ numPass;
				out.println(Misc.stringArrayToString(tokens, "\t"));
				
				//upper case ref and alt
				tokens[3] = tokens[3].toUpperCase();
				tokens[4] = tokens[4].toUpperCase();
			}
			in.close();
			out.close();
			IO.pl(numPass+"\t"+numFail);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: parsing UG vcf "+vcf+" \n"+line);
		} 
	}



	private double[] parseCounts(String formatLine, String values ) throws IOException {
		
		// check format
		if (formatLine.equals(format) == false) throw new IOException("ERROR: the format "+formatLine +" doesn't match "+format);

		// GT:AD:BG_AD:BG_DP:BG_SB:BG_VAF:DP:GQ:SB:VAF:PL
		//  0  1   2     3     4     5     6  7  8  9  10
		String[] t = Misc.COLON.split(values);
		
		String[] tumorRefAlt = Misc.COMMA.split(t[1]);
		double tR = Double.parseDouble(tumorRefAlt[0]);
		double tA = 0;
		//first is ref, alts follow
		for (int i=1; i< tumorRefAlt.length; i++) tA+= Double.parseDouble(tumorRefAlt[i]);
		
		String[] normalRefAlt = Misc.COMMA.split(t[2]);
		double nR = Double.parseDouble(normalRefAlt[0]);
		double nA = 0;
		for (int i=1; i< normalRefAlt.length; i++) nA+= Double.parseDouble(normalRefAlt[i]);
		
		return new double[] {tR, tA, nR, nA} ;
		
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new UltimaGenomicsVCFParser(args);
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
					case 's': saveDirectory = new File(args[++i]); break;
					case 't': minimumTumorAltFraction = Double.parseDouble(args[++i]); break;
					case 'n': maximumNormalAltFraction = Double.parseDouble(args[++i]); break;
					case 'u': minimumTumorReadDepth = Integer.parseInt(args[++i]); break;
					case 'a': minimumAltReadDepth = Integer.parseInt(args[++i]); break;
					case 'o': minimumNormalReadDepth = Integer.parseInt(args[++i]); break;
					case 'p': excludeNonPass = true; break;
					case 'c': excludeNonChr = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");

		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printErrAndExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");
		if (saveDirectory == null ) Misc.printErrAndExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");
		saveDirectory.mkdirs();
		
		IO.pl("\tminimumTumorAltFraction:\t"+ minimumTumorAltFraction);
		IO.pl("\tmaximumNormalAltFraction:\t"+ maximumNormalAltFraction);
		IO.pl("\tminimumTumorReadDepth:\t"+ minimumTumorReadDepth);
		IO.pl("\tminimumAltReadDepth:\t"+ minimumAltReadDepth);
		IO.pl("\tminimumNormalReadDepth:\t"+ minimumNormalReadDepth);
		IO.pl("\texcludeNonPass:\t"+ excludeNonPass);
		IO.pl("\texcludeNonChr:\t"+ excludeNonChr);
	}	

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Ultima Genomics VCF Parser: Jan 2026                   **\n" +
				"**************************************************************************************\n" +
				"Parses Ultima Genomics somatic VCF files, filtering for read depth and allele\n"+
				"frequency. Inserts T_AF, T_DP, N_AF, N_DP into the INFO field for integrated parsing\n"+
				"with the AnnotatedVcfParser. Upper cases REF and ALT.\n"+

				"\nRequired Options:\n"+
				"  -v Path to a file or directory containing xxx.vcf(.gz/.zip OK) file(s)\n" +
				"  -s Path to a directory to save the filtered vcfs\n"+
				"  -t Minimum tumor allele frequency, defaults to 0.05\n"+
				"  -u Minimum tumor read depth, defaults to 10\n"+
				"  -u Minimum tumor alt read depth, defaults to 3\n"+
				"  -n Maximum normal allele frequency, defaults to 0.1\n"+
				"  -o Minimum normal alignment depth, defaults to 10\n"+
				"  -p Remove non PASS filter field records\n"+
				"  -c Remove non chrXXX chromosome records\n"+
				

				"\nExample: java -jar pathToUSeq/Apps/UltimaGenomicsVCFParser -v RawVcfs/ -s FiltVcfs\n"+
				"     -p -t 0.1 -n 0.5 -u 50 -o 20 -u 3 -c \n\n"+

				"**************************************************************************************\n");

	}

}
