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

public class VCFToTable {
	private File[] vcfFiles;
	private ArrayList<String> columnsToUse = new ArrayList<String>();
	private ArrayList<String> columnsToSkip = new ArrayList<String>();
	private int reportStyle = VCFInfo.UNMODIFIED;
	private boolean genKey = false;
	private HashSet<String> toIgnore = new HashSet<String>();
	
	

	public VCFToTable(String[] args) {
		if (args.length == 0) {
			this.printDocs();
			System.exit(0);
		}
		
		long startTime = System.currentTimeMillis();
		
		//Parse the command line arguments
		this.processArgs(args);
		
		//Read in the VCF file
		for (File f: vcfFiles) {
			//Parse the vcf file
			VCFParser vcfFile = new VCFParser(f,true,true);
			
			//Grab the info part of the header
			HashMap<String,String> infoLines = vcfFile.getVcfComments().getInfo();
			ArrayList<String> infoOrder = vcfFile.getVcfComments().getInfoOrder();
			
			//Setup command and check for bad columns
			if (columnsToSkip.size() != 0) {
				for (String column: this.columnsToSkip) {
					if (!infoLines.containsKey(column)) {
						System.out.println("********************************** WARNING ****************************************");
						System.out.println("The VCF file does not contain the column : " + column );
					}
				}
				this.columnsToUse = VCFInfo.buildToAddFromToSkip(this.columnsToSkip, vcfFile.getStringComments());
			} else if (columnsToUse.size() != 0) {
				for (String column: this.columnsToUse) {
					for (String col: this.columnsToUse) {
						if (!infoLines.containsKey(col)) {
							System.out.println("********************************** WARNING ****************************************");
							System.out.println("The VCF file does not contain the column : " + column );
						}
					}
				}
			} else {
				this.columnsToUse = vcfFile.getVcfComments().getInfoOrder();
				this.reportStyle = VCFInfo.UNMODIFIED;
			}
			
			
			try {
				//Open up the output file
				String fullPathName = Misc.removeExtension(f.getCanonicalPath());
				String annPath = fullPathName + "_ann.txt";
				String keyPath = fullPathName + "_key.txt";
				
				BufferedWriter bw = new BufferedWriter(new FileWriter(annPath));
				
				//Generate the header for the table and write to file
				StringBuffer header = new StringBuffer("Chrom\tStart\tEnd\tReference\tAlt\tQual");
				for (String info: infoOrder) {
					if (this.columnsToUse.contains(info)) {
						header.append("\t" +info);
					}
				}
				for (String sample: vcfFile.getVcfComments().getSampleList()) {
					header.append("\t" + sample);
				}
				
				bw.write(header.toString() + "\n");
				
				
				//Write key if exists
				if (genKey) {
					BufferedWriter bwKey = new BufferedWriter(new FileWriter(keyPath));
					int keyIndex = 0;
					bwKey.write(String.valueOf(++keyIndex) + ". Chrom: Chromsome containing variant\n");
					bwKey.write(String.valueOf(++keyIndex) + ". Start: Variant start coordinate (1-based)\n");
					bwKey.write(String.valueOf(++keyIndex) + ". End: Variant end coordinate\n");
					bwKey.write(String.valueOf(++keyIndex) + ". Ref: Reference base\n");
					bwKey.write(String.valueOf(++keyIndex) + ". Alt: Alternate base(s)\n");
					bwKey.write(String.valueOf(++keyIndex) + ". Qual: Phred-scaled quality score for variant.  -10*log10(probability alt call is wrong \n");
					
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
						bw.write(record.getSpreadsheetOutput(this.columnsToUse, this.reportStyle));
					}
				}
				
				bw.close();
				
				
			} catch (IOException ioex) {
				System.out.println("Error writing the output for file: " + f.getName());
				System.exit(1);
			}
			
			
		}
		
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
		
	}
	
	
	public static void main(String[] args) {
		new VCFToTable(args);

	}

	
	private void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          VCF To Excel Spreadsheet: March 2013                     **\n" +
				"**************************************************************************************\n" +
				"This application takes a VCF file and generates and excel spreadsheet containing the \n" +
				"data. Users can specify what info fields they want to include/exclude in the output \n" +
				"or specify nothing to get the full list.  There is also an option to generate a table \n" +
				"key along with the output table\n\n\n" +

				"Required:\n"+
				"-v VCF file(s). Full path to a multi sample vcf file or directory containing such\n" +
				"      (xxx.vcf(.gz/.zip OK)).\n\n"+
				"Optional:\n"+
				"-d Desired Columns.  A comma-separated list of Info-field names that will be reported \n" +
				"      in the output vcf.\n" +
				"-u Unwanted Columns. A comma-separated list of Info-field names that will not be reported \n" +
				"      in the output vcf\n" +
				"-s Reporting Style. Info field styles.  Only two styles are currently supported, unmodified \n" +
				"      and short. Unmodified is used by default.  Short truncates some of the longer fields\n" +
				"-k Generate key.  Text document that lists descriptions of each column in the output table\n" +
				"-x Damaging only.  Only report nonsynonymous, frameshift or splicing variants\n" +
				"-a Annotations only.  Skip info fields that are reported by GATK.\n" +
				"\n\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/VCFInfoEdit -v 9908R.vcf -a SIFT,PP2,LRT,MT \n" +
				"      -c short \n" +
		        "**************************************************************************************\n");

	}
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		
		File forExtraction = null;
		String desiredColumns = null;
		String unwantedColumns = null;
		String reporting = null;
		boolean damaging = false;
		boolean standard = false;
		
		HashMap <String,Integer> allowedStyles = new HashMap<String,Integer>() {{
			put("UNMODIFIED",VCFInfo.UNMODIFIED);
			put("SHORT",VCFInfo.SHORT);
		}};
	
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'd': desiredColumns = args[++i]; break;
					case 'u': unwantedColumns = args[++i]; break;
					case 'c': reporting = args[++i]; break;
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
		
		if (desiredColumns != null && unwantedColumns != null) {
			System.out.println("********************************** WARNING ****************************************");
			System.out.println("User selected both desired columns and unwanted columns, using desired columns only");
		}
		
		if (standard) {
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
		} else if (desiredColumns != null) {
			String[] ata = desiredColumns.split(",");
			for (String a: ata) {
				String cleaned = a.trim();
				this.columnsToUse.add(cleaned);
			}
		} else if (unwantedColumns != null) {
			String[] ata = unwantedColumns.split(",");
			for (String a: ata) {
				String cleaned = a.trim();
				this.columnsToSkip.add(cleaned);
			}
		} 
		
		
		if (reporting != null) {
			if (!allowedStyles.containsKey(reporting.toUpperCase())) {
				System.out.println("**************************************** WARNING *****************************************************");
				System.out.println("Selected reporting style not recognized: " + reporting + ".  Falling back on unmodified");
			} else {
				this.reportStyle = allowedStyles.get(reporting.toUpperCase());
			}
		}
		 
		//pull vcf files (from David)
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please indicate a vcf file to filter.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");

	}

}
