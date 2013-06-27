package edu.utah.seq.vcf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class VCFReporter {
	//Globals
	private File vcfInFile;
	private File outFile;
	private String pathToTabix = "/tomato/app/tabix/";
	private boolean compressOutput = false;
	private ArrayList<String> columnsToUse;
	private ArrayList<String> columnsToSkip;
	private String reportStyle = VCFInfo.UNMODIFIED;
	private boolean genKey = false;
	private HashSet<String> toIgnore = new HashSet<String>();
	private String reportFormat = "TAB";

	public VCFReporter(String[] args) {
		if (args.length == 0) {
			this.printDocs();
			System.exit(0);
		}
		
		long startTime = System.currentTimeMillis();
		
		//Parse the command line arguments
		this.processArgs(args);
		

		//Parse the vcf file for info
		VCFParser parsedVcf = new VCFParser(this.vcfInFile,false,false, true);
		
		//Grab the info part of the header
		HashMap<String,String> infoLines = parsedVcf.getVcfComments().getInfo();
		
		
		//Setup command and check for bad columns
		if (columnsToSkip != null) {
			//Array that holds entries that actually exist
			ArrayList<String> exists = new ArrayList<String>();
			
			//Check each specified column to see if it exists.
			for (String column: this.columnsToSkip) {
				if (!infoLines.containsKey(column)) {
					System.out.println("********************************** WARNING ****************************************");
					System.out.println("The VCF file does not contain the column : " + column + ", skipping" );
				} else {
					exists.add(column);
				}
			}
			
			
			//Replace original list with existing
			this.columnsToSkip = exists;
			this.columnsToUse = VCFInfo.buildToAddFromToSkip(this.columnsToSkip, parsedVcf.getStringComments(this.columnsToSkip));
		} else if (columnsToUse != null) {
			//Array that holds entries that actually exist
			ArrayList<String> exists = new ArrayList<String>();
			
			//Check each specified column to see if it exists
			for (String column: this.columnsToUse) {
				if (!infoLines.containsKey(column)) {
					System.out.println("********************************** WARNING ****************************************");
					System.out.println("The VCF file does not contain the column : " + column + ", skipping");
				} else {
					exists.add(column);
				}
			}
			this.columnsToUse = exists;
			
		} else {
			this.columnsToUse = parsedVcf.getVcfComments().getInfoOrder();
			this.reportStyle = VCFInfo.UNMODIFIED;
		}
		
		//Count vcf records
		ArrayList<File> tempFiles = new ArrayList<File>();
		int recordCount = VCFUtilities.countReads(vcfInFile);
		int chunks = recordCount / VCFUtilities.readsToChunk + 1;
		
		if (chunks > 1) {
			System.out.println("File too big to report all in one go, splitting into " + chunks + " chunks");
		}
		
		//make tempDirectory
		File tempDir = new File("tempVCF");
		if (tempDir.exists()) {
			IO.deleteDirectory(tempDir);
		}
		tempDir.mkdir();
		
		for (int i=0;i<chunks;i++) {
			System.out.println("Working on file chunk: " + (i + 1));
			//Create temporary vcf files
			
			VCFParser parser =  new VCFParser(vcfInFile, true, true, true, i,VCFUtilities.readsToChunk);
			
			if (reportFormat.equals("TAB")) {
				File tempText = new File(tempDir,"tempTxt_" + i + ".txt");
				tempFiles.add(tempText);
				writeTabDelimited(parser,tempText,i);
				
			} else if (reportFormat.equals("VCF")) {
				File tempVcf = new File(tempDir,"tempVcf_" + i + ".vcf");
				tempFiles.add(tempVcf);
				writeVCF(parser,tempVcf);
			} else {
				System.out.println("Application does not recognize the report format: " + reportFormat + ", please contact the developers");
				System.exit(1);
			}
			
			parser = null;
		}
	
		if (reportFormat.equals("TAB")) {
			VCFUtilities.catFiles(tempFiles, this.outFile);
		}
		
		if (reportFormat.equals("VCF")) {
			VCFUtilities.mergeVcf(tempFiles, this.outFile);
			if (compressOutput) {
				VCFUtilities.createTabix(this.outFile, this.pathToTabix);
			}
		}
		
		IO.deleteDirectory(tempDir);
			
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
				
	}

	public static void main(String[] args) {
		new VCFReporter(args);

	}
	
	private void writeVCF(VCFParser vcfFile, File tempOutFile) {
		vcfFile.printRecords(tempOutFile, true, this.columnsToUse, this.reportStyle);
		
	}
	
	private void writeTabDelimited(VCFParser vcfFile, File tempOutFile, int fileCount) {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(tempOutFile));
			
			//Generate the header for the table and write to file
			StringBuffer header = new StringBuffer("Chrom\tStart\tEnd\tReference\tAlt\tQual");
			for (String column: columnsToUse) {
				header.append("\t" + column);
			}
				
			for (String sample: vcfFile.getVcfComments().getSampleList()) {
				header.append("\t" + sample);
			}
			
			//Only write the header for the first file
			if (fileCount == 0) {
				bw.write(header.toString() + "\n");
			}
			
			//Write key if exists
			if (genKey) {
				BufferedWriter bwKey = new BufferedWriter(new FileWriter(outFile.getName().substring(0,outFile.getName().length()-3) + "key.txt"));
				int keyIndex = 0;
				bwKey.write("Annotation Report Column Descriptions for: " + vcfInFile.getName() + "\n\n\n");
				bwKey.write("Style: " + this.reportStyle + "\n\n");
				
				bwKey.write(String.valueOf(++keyIndex) + ". Chrom: Chromosome containing variant\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Start: Variant start coordinate (1-based)\n");
				bwKey.write(String.valueOf(++keyIndex) + ". End: Variant end coordinate\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Ref: Reference base\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Alt: Alternate base(s)\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Qual: Phred-scaled quality score for variant.  -10*log10(probability alt call is wrong) \n");
				
				int columnIndex = -1;
				for (String desc: vcfFile.getVcfComments().getInfoDesc(this.columnsToUse)) {
					columnIndex += 1;
					keyIndex += 1;
					if (desc == null) {
						continue;
					}
					bwKey.write(String.valueOf(keyIndex) + ". " + this.columnsToUse.get(columnIndex)+ ": " + desc + "\n");
					
				}
				bwKey.write(String.valueOf(++keyIndex) + ". Genotype/Coverage data:  The remaining columns contain genotype and coverage data for the samples in the multi-sample VCF file. "
						+ "The format of the genotype/coverage string is GT:genotype:coverage, so in the case of: 'GT:0/1:4,5', 0/1 would be the genotype and 4,5 would be the coverage. "
						+ "A zero in the genotype field signifies the reference base and 1 through N represent alternate alleles.  0/1 would be heterozygous for first alternate allele, 0/0 homozygous for the reference, "
						+ "and 1/1 homozygous for the first alternate. A genotype of 0/2 would mean the sample was heterozygous for the second alternate allele.  The coverage field contains the "
						+ "number of reads representing each allele, separated by commas.  The first value is the number of observed reference bases, the second value is the number of observed "
						+ "first alternate bases.  There there was more than one alternate allele observed at the position, there will be more than two coverage values.");
				
				bwKey.close();
			}
			
			//Write data to files
			VCFRecord[] records = vcfFile.getVcfRecords();
			for (VCFRecord record: records) {
				if (!toIgnore.contains(record.getInfoObject().getInfo("VarType",VCFInfo.UNMODIFIED))) {
					bw.write(record.getSpreadsheetOutput(this.columnsToUse, this.reportStyle, vcfFile.getVcfComments()));
				}
			}
			
			bw.close();
			
			
		} catch (IOException ioex) {
			System.out.println("Error writing the output for file: " + vcfInFile.getName());
			System.exit(1);
		}
	}
	
	private void processArgs(String[] args){
		//Allowed options
		HashMap <String,String> allowedStyles = new HashMap<String,String>() {{
			put("UNMODIFIED",VCFInfo.UNMODIFIED);
			put("SHORT",VCFInfo.CLEAN);
		}};
		
		HashSet <String> allowedFormats = new HashSet<String>() {{
			add("TAB");
			add("VCF");
		}};
		
		//Local variables
		File inputFile = null;
		String outputName = null;
		String desiredColumns = null;
		String unwantedColumns = null;
		String style = null;
		boolean damaging = false;
		boolean standard = false;
		
		//Start parsing command line
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': inputFile = new File(args[++i]); break;
					case 'o': outputName = args[++i]; break;
					case 'd': desiredColumns = args[++i]; break;
					case 'u': unwantedColumns = args[++i]; break;
					case 'r': style = args[++i]; break;
					case 'k': this.genKey = true; break;
					case 'x': damaging = true; break;
					case 'a': standard = true; break;
					case 'h': printDocs(); System.exit(0);
					case 'p': this.pathToTabix = args[++i]; break;
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				} catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		
		if (damaging) {
			toIgnore.add("");
			toIgnore.add("synonymous_SNV");
			toIgnore.add("unknown");
		}
		
		
		if (standard) {
			this.columnsToUse = new ArrayList<String>();
			this.columnsToUse.add("EnsemblRegion");
			this.columnsToUse.add("EnsemblName");
			this.columnsToUse.add("VarType");
			this.columnsToUse.add("VarDesc");
			this.columnsToUse.add("RefSeq");
			this.columnsToUse.add("DBSNP");
			this.columnsToUse.add("ONEK");
			this.columnsToUse.add("COSMIC");
			this.columnsToUse.add("ESP");
			this.columnsToUse.add("SIFT");
			this.columnsToUse.add("PP2_HVAR");
			this.columnsToUse.add("PP2_HVAR_P");
			this.columnsToUse.add("PP2_HDIV");
			this.columnsToUse.add("PP2_HDIV_P");
			this.columnsToUse.add("LRT");
			this.columnsToUse.add("LRT_P");
			this.columnsToUse.add("MT");
			this.columnsToUse.add("MT_P");
			this.columnsToUse.add("MA");
			this.columnsToUse.add("MA_P");
			this.columnsToUse.add("FATHMM");
			this.columnsToUse.add("GERP");
			this.columnsToUse.add("PHYLOP");
			this.columnsToUse.add("SIPHY");
			this.columnsToUse.add("SEGDUP");
			this.columnsToUse.add("GWAS");
			this.columnsToUse.add("OMIM");
			this.columnsToUse.add("V_LRS");
			this.columnsToUse.add("V_RANK");
			this.columnsToUse.add("V_FLAG");
			this.columnsToUse.add("ACMG");
			this.columnsToUse.add("NIST");
			if (desiredColumns != null) {
				System.out.println("********************************** WARNING ****************************************");
				System.out.println("Since standard settings were selected, skipping option -d\n");
			}
			if (unwantedColumns != null) {
				System.out.println("********************************** WARNING ****************************************");
				System.out.println("Since standard settings were selcted, skipping option -u\n");
			}
			if (style != null) {
				System.out.println("********************************** WARNING ****************************************");
				System.out.println("Since standard settings were selected, skipping option -r\n");
			}
			this.reportStyle = VCFInfo.CLEAN;
		} else if (desiredColumns != null) {
			this.columnsToUse = new ArrayList<String>();
			String[] ata = desiredColumns.split(",");
			for (String a: ata) {
				String cleaned = a.trim();
				this.columnsToUse.add(cleaned);
			}
			
			if (unwantedColumns != null) {
				System.out.println("********************************** WARNING ****************************************");
				System.out.println("User selected both desired columns and unwanted columns, using desired columns only\n");
			}
		} else if (unwantedColumns != null) {
			this.columnsToSkip = new ArrayList<String>();
			String[] ata = unwantedColumns.split(",");
			for (String a: ata) {
				String cleaned = a.trim();
				this.columnsToSkip.add(cleaned);
			}
		} 
		
		if (style != null) {
			if (!allowedStyles.containsKey(style.toUpperCase())) {
				System.out.println("**************************************** WARNING *****************************************************");
				System.out.println("Selected reporting style not recognized: " + style + ".  Falling back on unmodified\n");
			} else {
				this.reportStyle = allowedStyles.get(style.toUpperCase());
			}
		}
		 
		if (inputFile == null) {
			System.out.println("Input file was not specified, exiting");
			System.exit(1);
		} else if (!inputFile.exists()) {
			System.out.println("Input file does not exist, exiting");
			System.exit(1);
		} else if (Pattern.matches(".+?.vcf",inputFile.getName())) {
			this.vcfInFile = inputFile;
		} else if (Pattern.matches(".+?.vcf.gz",inputFile.getName())) {
			this.vcfInFile = VCFUtilities.unzipTabix(inputFile,this.pathToTabix);
		} else {
			System.out.println("Input file does not appear to be a xxx.vcf/xxx.vcf.gz file");
			System.exit(1);
		}
		
		
		
		if (outputName == null) {
			System.out.println("Output file was no specified, exiting");
			System.exit(1);
		} else if (Pattern.matches(".+?.vcf",outputName)) {
			this.compressOutput = false;
			this.outFile = new File(outputName);
			this.reportFormat = "VCF";
		} else if (Pattern.matches(".+?.vcf.gz",outputName)) {
			File vcfOutComp = new File(outputName);
			if (vcfOutComp.exists()) {
				System.out.println("Tabix won't overwrite an existing file, rename the output or delete exisiting file");
				System.exit(1);
			}
			this.outFile = new File(outputName.substring(0,outputName.length()-3));
			this.compressOutput = true;
			this.reportFormat = "VCF";
		} else if (Pattern.matches(".+?.txt", outputName)) {
			this.outFile = new File(outputName);
			this.reportFormat = "TAB";
		} else {
			System.out.println("Output file does not appear to be a xxx.vcf/xxx.vcf.gz/xxx.txt file");
			System.exit(1);
		}

	}
	
	
	private void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            VCF Reporter: April 2013                              **\n" +
				"**************************************************************************************\n" +
				"This application takes a VCF file as input and returns either a modified VCF file or a\n" +
				"tab-delimited text file containing user-specified and optionally formatted INFO \n" +
				"fields.  The modified VCF file is useful if you want to view annotations in IGV.  The\n" +
				"standard set of INFO fields is quite large and can't fit in the IGV window. The tab-\n" +
				"delimited text file allows the annotations to be viewing in Excel for easier sorting \n" +
				"and filtering. If the number of VCF records is greater than 500,000, the reporting \n" +
				"will be done in chunks.  The chunks are merged and compressed automatically at the end\n" +
				"of the application.\n\n" +

				"Required:\n"+
				"-v VCF file. Full path to a multi sample vcf file (xxx.vcf(.gz/.zip OK)).\n"+
				"-o Output file.  Full path to the output file. If xxx.txt is specified, output will \n" +
				"      be a tab-delimited text file.  If xxx.vcf is specified, output will be an \n" +
				"      uncompressed vcf file.  If xxx.vcf.gz is specifed, output will be a tabix \n" + 
				"      compressed and indexed vcf file.\n" +
				"\nOptional:\n"+
				"-d Desired Columns.  A comma-separated list of INFO-field names that will be reported\n" +
				"      in the output vcf.\n" +
				"-u Unwanted Columns. A comma-separated list of INFO-field names that will not be \n" +
				"      reported in the output vcf.\n" +
				"-r Reporting Style. INFO field styles.  Only two styles are currently supported, \n" + 
				"      'unmodified' and 'short'. Unmodified is used by default.  Short truncates some\n" +
				"      of the longer fields, which help visibility in IGV.\n" +
				"-a Annotations only.  Report standard annotations using the 'short' reporting style.\n" +
				"      Skip INFO fields reported by GATK to reduce clutter.  The skipped fields are \n" +
				"      used by GATK to determine variation quality and might not be useful to the \n" +
				"      general user.\n" +
				"-x Damaging only.  Only report nonsynonymous, frameshift or splicing variants.\n" +
				"-p Path to tabix directory.  Set this variable if the application is not run on \n" +
				"      moab/alta\n" +
				"\nTab-delimited only options:\n" +
				"-k Generate key.  Text document that lists descriptions of each column in the output\n" +
				"      table.\n" +
				"\n\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/VCFReporter -v 9908R.vcf \n" + 
				"      -d SIFT,LRT,MT,MT_P -r short -o 9908.ann.txt \n\n" +
				"**************************************************************************************\n");

	}

}
