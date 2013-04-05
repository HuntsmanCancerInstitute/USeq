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
		


		//Parse the vcf file
		VCFParser parsedVcf = new VCFParser(this.vcfInFile,true,true);
		
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
		
		if (reportFormat.equals("TAB")) {
			writeTabDelimited(parsedVcf);
		} else if (reportFormat.equals("VCF")) {
			writeVCF(parsedVcf);
		} else {
			System.out.println("Application does not recognize the report format: " + reportFormat + ", please contact the developers");
			System.exit(1);
		}
	
		
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
				
	}

	public static void main(String[] args) {
		new VCFReporter(args);

	}
	
	private void writeVCF(VCFParser vcfFile) {
		vcfFile.printRecords(this.outFile, true, this.columnsToUse, this.reportStyle);
		if (compressOutput) {
			VCFUtilities.createTabix(this.outFile, this.pathToTabix);
		}
	}
	
	private void writeTabDelimited(VCFParser vcfFile) {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			
			//Generate the header for the table and write to file
			StringBuffer header = new StringBuffer("Chrom\tStart\tEnd\tReference\tAlt\tQual");
			for (String column: columnsToUse) {
				header.append("\t" + column);
			}
				
			for (String sample: vcfFile.getVcfComments().getSampleList()) {
				header.append("\t" + sample);
			}
			
			bw.write(header.toString() + "\n");
			
			
			//Write key if exists
			if (genKey) {
				BufferedWriter bwKey = new BufferedWriter(new FileWriter(outFile.getName() + "_key"));
				int keyIndex = 0;
				bwKey.write("Annotation Report Column Descriptions for: " + vcfInFile.getName() + "\n\n\n");
				bwKey.write("Style: " + this.reportStyle + "\n\n");
				
				bwKey.write(String.valueOf(++keyIndex) + ". Chrom: Chromsome containing variant\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Start: Variant start coordinate (1-based)\n");
				bwKey.write(String.valueOf(++keyIndex) + ". End: Variant end coordinate\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Ref: Reference base\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Alt: Alternate base(s)\n");
				bwKey.write(String.valueOf(++keyIndex) + ". Qual: Phred-scaled quality score for variant.  -10*log10(probability alt call is wrong) \n");
				
				int columnIndex = -1;
				for (String desc: vcfFile.getVcfComments().getInfoDesc(this.columnsToUse)) {
					columnIndex += 1;
					if (desc == null) {
						continue;
					}
					bwKey.write(String.valueOf(keyIndex) + ". " + this.columnsToUse.get(columnIndex)+ ": " + desc + "\n");
					keyIndex += 1;
				}
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
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
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
			this.columnsToUse.add("PP2");
			this.columnsToUse.add("MT");
			this.columnsToUse.add("LRT");
			this.columnsToUse.add("PHYLOP");
			this.columnsToUse.add("SEGDUP");
			this.columnsToUse.add("GWAS");
			this.columnsToUse.add("OMIM");
			this.columnsToUse.add("V_LRS");
			this.columnsToUse.add("V_RANK");
			this.columnsToUse.add("V_FLAG");
			this.columnsToUse.add("ACMG");
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
				"**                          VCF Reporter: April 2013                                **\n" +
				"**************************************************************************************\n" +
				"This application takes a vcf file as input and returns either a modified VCF file or \n" +
				"a tab-delimited text file containing user-specified and formatted INFO fields.  The \n" +
				"modified VCF file is useful if you want to view annotations in IGV.  The standard set of \n" +
				"INFO fields is quite large and can't fit in the IGV window. The tab-delimited text file \n" +
				"makes the annotations much easier to view sort and filter.  Notice that some options \n" +
				"are only meaningful for the tab-delimited text file\n\n\n" +

				"Required:\n"+
				"-v VCF file. Full path to a multi sample vcf file (xxx.vcf(.gz/.zip OK)).\n"+
				"-o Output file.  Full path to the output file. If xxx.txt specified, output will be a tab-delimited\n" +
				"      spreadsheet.  If xxx.vcf specified, output will be an uncompressed vcf file.  If \n" +
				"      xxx.vcf.gz specifed, output will be a tabix compressed and indexed vcf file.\n" +
				"\nOptional:\n"+
				"-d Desired Columns.  A comma-separated list of Info-field names that will be reported \n" +
				"      in the output vcf.\n" +
				"-u Unwanted Columns. A comma-separated list of Info-field names that will not be reported \n" +
				"      in the output vcf\n" +
				"-r Reporting Style. Info field styles.  Only two styles are currently supported, unmodified \n" +
				"      and short. Unmodified is used by default.  Short truncates some of the longer fields\n" +
				"-a Annotations only.  Report standard annotations in CLEAN format.  Skip info fields reported \n" + 
				"      by GATK to reduce clutter in IGV.\n" +
				"-x Damaging only.  Only report nonsynonymous, frameshift or splicing variants\n" +
				"\nTab-delimited only options\n" +
				"-k Generate key.  Text document that lists descriptions of each column in the output table\n" +
				"\n\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/VCFReporter -v 9908R.vcf -d SIFT,PP2,LRT,MT \n" +
				"      -r CLEAN -o 9908.ann.txt \n" +
		        "**************************************************************************************\n");

	}

}
