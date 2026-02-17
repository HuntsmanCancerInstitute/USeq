package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Scans the ID column of vcf records and creates a INFO column element if any of the centers are found.
 * Foundation_12;Strelka_427    Foundation
 * Strelka_51;Tempus_4          Tempus
 * Strelka_269 					AVATAR
 * @author Nix*/
public class VCFCenterFormatter {
	
	//user defined fields
	private File forExtraction = null;
	private File[] vcfFiles;
	private File saveDir;
	private String info = "##INFO=<ID=CIV,Number=.,Type=String,Description=“Center identifying variant”>";
	private String cmd;
	private String[] centers  = new String[] {"Foundation", "Tempus", "HCI", "Caris", "Ambry", "Invitae"};
	private String[] toLookForLowerCase  = new String[] {"foundation", "tempus", "strelka", "caris", "ambry", "invitae"};
	
	//constructor
	public VCFCenterFormatter(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		//for each vcf file
		IO.pl("\nParsing vcf files...");
		for (int i=0; i< vcfFiles.length; i++){
			IO.pl("\t"+ vcfFiles[i]);
			parse (vcfFiles[i]);
		}
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		IO.pl("\nDone! "+Math.round(diffTime)+" min\n");
	}

	


	private void parse(File vcfFile) {
		Gzipper out = null;
		BufferedReader vcfIn = null;
		try {
			//make IO
			String name = Misc.removeExtension(vcfFile.getName());
			out = new Gzipper (new File (saveDir, name+".civ.vcf.gz"));
			vcfIn = IO.fetchBufferedReader(vcfFile);
			
			//parse the file
			boolean addInfo = true;
			String line;
			ArrayList<String> infoToAdd = new ArrayList<String>();
			while ((line = vcfIn.readLine())!= null) {
				//comment line?
				if (line.startsWith("#")) {
					if (line.startsWith("#CHROM")) line = "##" +cmd+"\n"+line;
					else if (addInfo == true && line.startsWith("##INFO")) {
						line = info+"\n"+line;
						addInfo = false;
					}
					
				}
				//must be a data line
				else {
					infoToAdd.clear();
					//#CHROM	POS	     ID	      REF	ALT	QUAL	FILTER	INFO	                        FORMAT	carcinoid5	carcinoid6	carcinoid7	colon1	colon2	colon3	colon4
					//chr1	  714427  rs12028261	G	A	57.36	PASS	AC=4;AF=1.00;AN=4;DB;DP=2;Dels=0.00;FS=0.000;HaplotypeScore=0.0000;MLEAC=4;MLEAF=1.00;MQ=70.00;MQ0=0;QD=28.68;VQSLOD=22.08;culprit=HaplotypeScore	GT:AD:DP:GQ:PL	./.	1/1:0,1:1:3:39,3,0	./.	./.	./.	./.	1/1:0,1:1:3:41,3,0
					//  0		1        2			3	4	5		6		7
					//pull ID
					String[] t = Misc.TAB.split(line);
					String lcIds = t[2].toLowerCase();
					//for each to look for
					for (int i=0; i< toLookForLowerCase.length; i++) {
						if (lcIds.contains(toLookForLowerCase[i])) infoToAdd.add(centers[i]);
					}
					//any found?
					if (infoToAdd.size() !=0) {
						String combineCenters = Misc.stringArrayListToString(infoToAdd, ",");
						t[7] = t[7]+ ";CIV="+ combineCenters;
						line = Misc.stringArrayToString(t, "\t");
					}
					//nope, don't modify line
				}
				out.println(line);
			}
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("ERROR parsing "+vcfFile);
		} finally {
			out.closeNoException();
			IO.closeNoException(vcfIn);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFCenterFormatter(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		cmd = IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ");
		IO.pl("\n"+cmd+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 's': saveDir = new File(args[++i]); break;
					case 'c': centers = Misc.COMMA.split(args[++i]); break;
					case 'i': toLookForLowerCase = Misc.COMMA.split(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		
		//save dir?
		if (saveDir == null) saveDir = vcfFiles[0].getParentFile();
		else saveDir.mkdirs();
		if (saveDir.exists()==false || saveDir.canWrite()==false) Misc.printExit("\nError: failed to find or create a writeable save directory? Aborting.\n");
		
		//check centers and replacements
		printSettings();
		
		//check arrays
		if (centers.length != toLookForLowerCase.length) Misc.printErrAndExit("\nError: the -c and -i array lengths differ.\n");
		
	}	
	
	public void printSettings(){
		IO.pl("Settings:");
		IO.pl(" -v Vcf file  "+ forExtraction);
		IO.pl(" -s Save dir  "+ IO.getCanonicalPath(saveDir));
		IO.pl(" -i IDs to look for  "+ Misc.stringArrayToString(toLookForLowerCase, ","));
		IO.pl(" -c Center names  "+ Misc.stringArrayToString(centers, ","));
	} 
	
	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                            VCF Center Formatter : Dec 2025                       **\n" +
				"**************************************************************************************\n" +
				"Looks for key words in each VCF record's ID field, if found appends the matching\n"+
				"center designation to the info field. The -i and -c arrays must be equal.\n"+

				"\nOptions:\n"+
				"-v Path to a vcf file xxx.vcf(.gz/.zip OK) or directory containing such.\n" +
				"-s Path to a directory to save the modified vcf files.\n"+
				"-i Comma delimited list, no spaces, lower case, of IDs to look for, defaults to \n"+
				"      foundation,tempus,strelka,caris,ambry,invitae\n"+
				"-c Comma delimited list, no spaces, of Center names to append to the CIV= info field\n"+
				"      defaults to Foundation,Tempus,HCI,Caris,Ambry,Invitae\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/VCFCenterFormatter -v Calls/ -s CCalls\n\n"+

		        "**************************************************************************************\n");
	}

}
