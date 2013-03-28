package edu.utah.seq.vcf;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class VCFInfoEdit {
	private File[] vcfFiles;
	private ArrayList<String> columnsToUse = new ArrayList<String>();
	private ArrayList<String> columnsToSkip = new ArrayList<String>();
	private int reportStyle = VCFInfo.UNMODIFIED;
	

	public VCFInfoEdit(String[] args) {
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
			
			if (columnsToSkip.size() != 0) {
				for (String column: this.columnsToSkip) {
					if (!infoLines.containsKey(column)) {
						System.out.println("********************************** WARNING ****************************************");
						System.out.println("The VCF file does not contain the column : " + column );
					}
				}
				this.columnsToUse = VCFInfo.buildToAddFromToSkip(this.columnsToSkip, vcfFile.getStringComments());
				vcfFile.printRecords(VCFRecord.PASS, true, this.columnsToUse, this.reportStyle);
			} else if (columnsToUse.size() != 0) {
				for (String column: this.columnsToUse) {
					for (String col: this.columnsToUse) {
						if (!infoLines.containsKey(col)) {
							System.out.println("********************************** WARNING ****************************************");
							System.out.println("The VCF file does not contain the column : " + column );
						}
					}
				}
				vcfFile.printRecords(VCFRecord.PASS, true, this.columnsToUse, this.reportStyle);
			} else {
				this.columnsToUse = vcfFile.getVcfComments().getInfoOrder();
				this.reportStyle = VCFInfo.UNMODIFIED;
				vcfFile.printRecords(VCFRecord.PASS, true, this.columnsToUse, this.reportStyle);
			}
			
			
		}
		
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
		
	}
	
	
	
	
	

	public static void main(String[] args) {
		new VCFInfoEdit(args);

	}

	
	private void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          VCF file info field editor: March 2013                   **\n" +
				"**************************************************************************************\n" +
				"Allows the user to select a subset of info fields to report.  This reduces file size \n" + 
				"and allows more useful fields to be visible in IGV.  Running with no columns selcted \n" +
				"should create an identical copy of the original VCF file. This can be used to check \n" +
				"that the VCF classes are all functional\n\n\n" +

				"Required:\n"+
				"-v VCF file(s). Full path to a multi sample vcf file or directory containing such\n" +
				"      (xxx.vcf(.gz/.zip OK)).\n\n"+
				"Optional:\n"+
				"-a Desired Columns.  A comma-separated list of Info-field names that will be reported \n" +
				"      in the output vcf.\n" +
				"-b Unwanted Columns. A comma-separated list of Info-field names that will not be reported \n" +
				"      in the output vcf\n" +
				"-c Reporting Style. Info field styles.  Only two styles are currently supported, unmodified \n" +
				"      and short. Unmodified is used by default.  Short truncates some of the longer fields\n" +
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
					case 'a': desiredColumns = args[++i]; break;
					case 'b': unwantedColumns = args[++i]; break;
					case 'c': reporting = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (desiredColumns != null && unwantedColumns != null) {
			System.out.println("********************************** WARNING ****************************************");
			System.out.println("User selected both desired columns and unwanted columns, using desired columns only");
		}
		
		//Set up annotations to test, either full set, or user-specified subset.
		if (desiredColumns != null) {
			String[] ata = desiredColumns.split(",");
			for (String a: ata) {
				String cleaned = a.trim();
				this.columnsToUse.add(a);
			}
		} else if (unwantedColumns != null) {
			String[] ata = unwantedColumns.split(",");
			for (String a: ata) {
				String cleaned = a.trim();
				this.columnsToSkip.add(a);
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
