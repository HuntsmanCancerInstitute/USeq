package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class VCFAnnovar {
	//Hard coded paths
	private String pathToAnnovarDir = new String("/home/u0855942/annovar_version");
	private String pathToAnnovar = new String(pathToAnnovarDir + "/annotate_variation.pl");
	private String pathToRespository = new String(pathToAnnovarDir + "/humandb");
	
	//User specified options
	private File[] vcfFiles;
	private String dbSnpFile;
	
	//Shared variables
	private HashMap<String,AnnovarCommand> commandMap = new HashMap<String,AnnovarCommand>() {{
		put("REFSEQ",new AnnovarCommand("RefSeq Annotations"));
		put("ENSEMBL",new AnnovarCommand("Ensembl Annotations"));
		put("SEGDUP",new AnnovarCommand("Segdup Annotations"));
		put("GWAS",new AnnovarCommand("GWAS Annotations"));
		put("DBSNP",new AnnovarCommand("dbSNP Annotations"));
		put("PHYLOP",new AnnovarCommand("Phylop Annotations"));
		put("COSMIC",new AnnovarCommand("COSMIC Annotations"));
		put("ESP",new AnnovarCommand("ESP Annotations"));
		put("SIFT",new AnnovarCommand("Sift Annotations"));
		put("PP2",new AnnovarCommand("Polyphen2 Annotations"));
		put("ONEK",new AnnovarCommand("1K Genome Annotations"));
		put("MT",new AnnovarCommand("MutationTaster Annotations"));
		put("LRT",new AnnovarCommand("LRT Annotations"));
		put("TFBS",new AnnovarCommand("TFBS Annotations"));
		put("DGV",new AnnovarCommand("DGV Annotations"));
		put("OMIM",new AnnovarCommand("OMIM Annotations"));
	}};
	
	private HashMap<String,String> ethnicityMap = new HashMap<String,String>() {{
		put("ALL","1000g2012apr_all");
		put("EUR","1000g2012apr_eur");
		put("AFR","1000g2012apr_afr");
		put("AMR","1000g2012apr_amr");
		put("ASN","1000g2012apr_asn");
	}};
	
	private String usedEthnicity;
	private String inputname;
	private ArrayList<String> annsToRun = new ArrayList<String>();
	
	
	public VCFAnnovar(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		
		long startTime = System.currentTimeMillis();

		processArgs(args);

		//for each file
		for (int i=0; i< vcfFiles.length; i++){
			VCFParser parsedVCF = new VCFParser(vcfFiles[i], true, true);
	
			String[] fileParts = vcfFiles[i].getName().split("\\.(?=[^\\.]+$)");
	
			inputname = fileParts[0] + ".annovar.input.txt";

			//Setup annovar commands
			this.setupCommands();
			
			//Write annovar input file
			this.writeAnnovarInput(parsedVCF);
			
			//Run annovarDriver
			for (String command: this.annsToRun) {
				commandMap.get(command).runCommand(parsedVCF);
			}
			
		
			//Write VCF file
			parsedVCF.printRecords(VCFRecord.PASS,true);
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");	
	}
	
	
	
	
	private void writeAnnovarInput(VCFParser parser) {
		File annovarInput = new File(inputname);
		//annovarInput.deleteOnExit();
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(annovarInput));
			
			for(VCFRecord vr: parser.getVcfRecords()) {
				String chrom = vr.getChromosome().replace("chr", "");
				String endPosition = String.valueOf(this.getAnnovarEndPosition(vr));
				String altAllele = vr.getAlternate().split(",")[0];
				bw.write(chrom + "\t" + String.valueOf(vr.getPosition() + 1) + "\t" + endPosition + "\t" + vr.getReference() + "\t" + altAllele + "\n");
			}
			
			bw.flush();
			bw.close();
			
		} catch(IOException ioex) {
			System.out.println("Error writing to file: " + annovarInput.getName());
			System.exit(1);
			
		}
	}
	
	private Integer maxInt(Integer val1, Integer val2) {
		if (val1 >= val2) {
			return val1;
		} else {
			return val2;
		}
	}
	
	private Integer getAnnovarEndPosition(VCFRecord vr) {
		Integer endPos = vr.getPosition() + 1 + this.maxInt(vr.getReference().length()-1,0);
		return endPos;
	}
	
    private void setupCommands() {
    	String cmdName = null;
    	
    	//Annovar gene
    	cmdName = "ENSEMBL";
    	ProcessBuilder pbEnsembl = new ProcessBuilder(this.pathToAnnovar,"--geneanno","--buildver","hg19","-dbtype","ensgene",this.inputname,this.pathToRespository);
    	String il1Ensembl = new String("##INFO=<ID=VarType,Number=1,Type=String,Description=\"Exonic variant location (Annovar)\">");
    	String il2Ensembl = new String("##INFO=<ID=VarDesc,Number=1,Type=String,Description=\"Exonic variant description (Annovar)\">");
    	String il3Ensembl = new String("##INFO=<ID=EnsemblRegion,Number=1,Type=String,Description=\"Location of variation (Annovar)\">");
    	String il4Ensembl = new String("##INFO=<ID=EnsemblName,Number=1,Type=String,Description=\"Closest Ensembl gene (Annovar)\">");
    	OutputParser op1Ensembl = new OutputParser(new int[]{0,1},new String[]{"EnsemblRegion","EnsemblName"},new String[]{il3Ensembl,il4Ensembl},null,"STANDARD","variant_function");
    	OutputParser op2Ensembl = new OutputParser(new int[]{1,2},new String[]{"VarType","VarDesc"},new String[]{il1Ensembl,il2Ensembl},3,"UNSORTED","exonic_variant_function");
    	commandMap.get(cmdName).addCommand(pbEnsembl);
    	commandMap.get(cmdName).addOutputParser(op1Ensembl);
    	commandMap.get(cmdName).addOutputParser(op2Ensembl);
    	
    	//Annovar refseq
    	cmdName = "REFSEQ";
    	ProcessBuilder pbRefseq = new ProcessBuilder(this.pathToAnnovar,"--geneanno","--buildver","hg19","-dbtype","gene",this.inputname,this.pathToRespository);
    	String il1Refseq = new String("##INFO=<ID=RefSeq,Number=1,Type=String,Description=\"Closest refseq gene (Annovar)\">");
    	OutputParser op1Refseq = new OutputParser(new int[]{1},new String[]{"Refseq"},new String[]{il1Refseq},null,"STANDARD","variant_function");
    	OutputParser op2Refseq = new OutputParser("exonic_variant_function");
    	commandMap.get(cmdName).addCommand(pbRefseq);
    	commandMap.get(cmdName).addOutputParser(op1Refseq);
    	commandMap.get(cmdName).addOutputParser(op2Refseq);
    	
    	//Annovar TFBS
    	cmdName = "TFBS";
    	ProcessBuilder pbTfbs = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","tfbs",this.inputname,this.pathToRespository);
    	String il1Tfbs = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"TFBS identifier (Annovar)\">");
    	OutputParser op1Tfbs = new OutputParser(cmdName,il1Tfbs,"hg19_tfbsConsSites");
    	commandMap.get(cmdName).addCommand(pbTfbs);
    	commandMap.get(cmdName).addOutputParser(op1Tfbs);
  
    	//Annovar Segdup
    	cmdName = "SEGDUP";
    	ProcessBuilder pbSegdup = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","segdup",this.inputname ,this.pathToRespository);
    	String il1Segdup = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"Segdup identifier (Annovar)\">");
    	OutputParser op1Segdup = new OutputParser(cmdName,il1Segdup,"hg19_genomicSuperDups");
    	commandMap.get(cmdName).addCommand(pbSegdup);
    	commandMap.get(cmdName).addOutputParser(op1Segdup);
    	
    	//Annovar DGV
    	cmdName = "DGV";
    	ProcessBuilder pbDgv = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","dgv",this.inputname,this.pathToRespository);
    	String il1Dgv = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"DGV identifier (Annovar)\">");
    	OutputParser op1Dgv= new OutputParser(cmdName,il1Dgv,"hg19_dgv");
    	commandMap.get(cmdName).addCommand(pbDgv);
    	commandMap.get(cmdName).addOutputParser(op1Dgv);
    	
    	//Annovar DBSNP
    	cmdName = "DBSNP";
    	ProcessBuilder pbSnp = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",this.dbSnpFile,this.inputname,this.pathToRespository);
    	String il1Snp = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"dbSNP identifer, " + this.dbSnpFile + " (Annovar) \">");
    	OutputParser op1Snp = new OutputParser(this.dbSnpFile,il1Snp,"hg19_" + this.dbSnpFile + "_dropped");
    	OutputParser op2Snp = new OutputParser("hg19_" + this.dbSnpFile + "_filtered");
    	commandMap.get(cmdName).addCommand(pbSnp);
    	commandMap.get(cmdName).addOutputParser(op1Snp);
    	commandMap.get(cmdName).addOutputParser(op2Snp);
    	
    	//Annovar SIFT
    	cmdName = "SIFT";
    	ProcessBuilder pbSift = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype","avsift",this.inputname,this.pathToRespository);
    	String il1Sift = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=Float,Description=\"SIFT score (Annovar)\">");
    	OutputParser op1Sift = new OutputParser(cmdName,il1Sift,"hg19_avsift_dropped");
    	OutputParser op2Sift = new OutputParser("hg19_avsift_filtered");
    	commandMap.get(cmdName).addCommand(pbSift);
    	commandMap.get(cmdName).addOutputParser(op1Sift);
    	commandMap.get(cmdName).addOutputParser(op2Sift);
    	
    	//Annovar Polyphen
    	cmdName = "PP2";
    	ProcessBuilder pbPoly = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype","ljb_pp2",this.inputname,this.pathToRespository);
    	String il1Poly = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=Float,Description=\"Polyphen2 score (Annovar)\">");
    	OutputParser op1Poly = new OutputParser(cmdName,il1Poly,"hg19_ljb_pp2_dropped");
    	OutputParser op2Poly = new OutputParser("hg19_ljb_pp2_filtered");
    	commandMap.get(cmdName).addCommand(pbPoly);
    	commandMap.get(cmdName).addOutputParser(op1Poly);
    	commandMap.get(cmdName).addOutputParser(op2Poly);
    	
    	//Annovar Mutation Taster
    	cmdName = "MT";
    	ProcessBuilder pbMt = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype","ljb_mt",this.inputname,this.pathToRespository);
    	String il1Mt = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=Float,Description=\"Mutation taster score (Annovar)\">");
    	OutputParser op1Mt = new OutputParser(cmdName,il1Mt,"hg19_ljb_mt_dropped");
    	OutputParser op2Mt = new OutputParser("hg19_ljb_mt_filtered");
    	commandMap.get(cmdName).addCommand(pbMt);
    	commandMap.get(cmdName).addOutputParser(op1Mt);
    	commandMap.get(cmdName).addOutputParser(op2Mt);
    	
    	//Annovar LRT
    	cmdName = "LRT";
    	ProcessBuilder pbLrt = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype","ljb_lrt",this.inputname,this.pathToRespository);
    	String il1Lrt = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=Float,Description=\"LRT score (Annovar)\">");
    	OutputParser op1Lrt = new OutputParser(cmdName,il1Lrt,"hg19_ljb_lrt_dropped");
    	OutputParser op2Lrt = new OutputParser("hg19_ljb_lrt_filtered");
    	commandMap.get(cmdName).addCommand(pbLrt);
    	commandMap.get(cmdName).addOutputParser(op1Lrt);
    	commandMap.get(cmdName).addOutputParser(op2Lrt);
    	
    	//Annovar Phylop
    	cmdName = "PHYLOP";
    	ProcessBuilder pbPhylop = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype","ljb_phylop",this.inputname,this.pathToRespository);
    	String il1Phylop = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=Float,Description=\"Phylop score (Annovar)\">");
    	OutputParser op1Phylop = new OutputParser(cmdName,il1Phylop,"hg19_ljb_phylop_dropped");
    	OutputParser op2Phylop = new OutputParser("hg19_ljb_phylop_filtered");
    	commandMap.get(cmdName).addCommand(pbPhylop);
    	commandMap.get(cmdName).addOutputParser(op1Phylop);
    	commandMap.get(cmdName).addOutputParser(op2Phylop);
    	

    	//Annovar 1K Genomes
    	cmdName = "ONEK";
    	ProcessBuilder pbOneK = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",this.ethnicityMap.get(this.usedEthnicity),this.inputname,this.pathToRespository);
    	String il1OneK = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=Float,Description=\"1000 Genomes observation Frequency " + this.ethnicityMap.get(this.usedEthnicity) + " (Annovar) \">");
    	OutputParser op1OneK = new OutputParser(cmdName,il1OneK,"hg19_" + this.usedEthnicity + ".sites.2012_04_dropped");
    	OutputParser op2OneK = new OutputParser("hg19_" + this.usedEthnicity + ".sites.2012_04_filtered");
    	commandMap.get(cmdName).addCommand(pbOneK);
    	commandMap.get(cmdName).addOutputParser(op1OneK);
    	commandMap.get(cmdName).addOutputParser(op2OneK);
    	
    	//Annovar COSMIC
    	cmdName = "COSMIC";
    	ProcessBuilder pbCosmic = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype","cosmic63",this.inputname,this.pathToRespository);
    	String il1Cosmic = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"COSMIC annotation (Annovar)\">");
    	OutputParser op1Cosmic = new OutputParser(cmdName,il1Cosmic,"hg19_cosmic63_dropped");
    	OutputParser op2Cosmic = new OutputParser("hg19_cosmic63_filtered");
    	commandMap.get(cmdName).addCommand(pbCosmic);
    	commandMap.get(cmdName).addOutputParser(op1Cosmic);
    	commandMap.get(cmdName).addOutputParser(op2Cosmic);
    	
    	//Annovar ESP
    	cmdName = "ESP";
    	ProcessBuilder pbEsp = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype","esp6500_all",this.inputname,this.pathToRespository);
    	String il1Esp = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"ESP annotation (Annovar)\">");
    	OutputParser op1Esp = new OutputParser(cmdName,il1Esp,"hg19_esp6500_all_dropped");
    	OutputParser op2Esp = new OutputParser("hg19_esp6500_all_filtered");
    	commandMap.get(cmdName).addCommand(pbEsp);
    	commandMap.get(cmdName).addOutputParser(op1Esp);
    	commandMap.get(cmdName).addOutputParser(op2Esp);
    	
    	//Annovar GWAS
    	cmdName = "GWAS";
    	ProcessBuilder pbGwas = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","gwascatalog",this.inputname,this.pathToRespository);
    	String il1Gwas = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"GWAS catalog annotation (Annovar)\">");
    	OutputParser op1Gwas = new OutputParser(cmdName,il1Gwas,"hg19_gwasCatalog");
    	commandMap.get(cmdName).addCommand(pbGwas);
    	commandMap.get(cmdName).addOutputParser(op1Gwas);
    	
    	//Annovar OMIM
    	cmdName = "OMIM";
    	ProcessBuilder pbOmim = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-gff3attrib","-dbtype","gff3","-gff3dbfile","hg19_omim.gff3",this.inputname,this.pathToRespository);
    	String il1Omim = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"OMIM annotations (UofU)\">");
    	OutputParser op1Omim = new OutputParser(cmdName,il1Omim,"hg19_gff3");
    	commandMap.get(cmdName).addCommand(pbOmim);
    	commandMap.get(cmdName).addOutputParser(op1Omim);
 
    }
	

	public static void main(String[] args) {
		new VCFAnnovar(args);
	}
	
	private void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Multi Sample VCF Annotation: March 2013                 **\n" +
				"**************************************************************************************\n" +
				"Adds annotations from several different sources to VCF file INFO line.  Only hg19 \n " +
				"is supported at this time.\n\n\n" +

				"Required:\n"+
				"-v VCF file(s). Full path to a multi sample vcf file or directory containing such\n" +
				"      (xxx.vcf(.gz/.zip OK)).\n\n"+
				"Optional:\n"+
				"-d dbSNP database.  By default, this applications uses dbSNP 137 for annotation. Use \n" +
				"      this option along with the database identifier to use a different version, ie snp129 \n" +
				"      The dbSNP database must be in the annovar data directory\n" +
				"-e Ethnicity.  By default, the 1K frequency is calculated across all ethnicities.  \n" +
				"      If you want to restrict it to one of EUR, AFR, ASN or AMR, use this option \n" +
				"      followed by the ethnicity identifier. \n" +
				"-a Annotations to add.  By default, this application uses all available annovar \n" +
				"      annotations.  Use a comma-separated list of keys to specify a custom set. Avaiable\n" +
				"      annotations with (keys): refSeq gene (GENE), transcription factor binding site \n " +
				"      (TFBS), segmental duplicaton (SEGDUP), database of genomic variants (DGV), sift \n" +
				"      score (SIFT), polyphen score (POLY), 1K genomes frequency (TENK) and dbSNP rs \n" +
				"      (DBSNP)\n"+
				"-p Path to annovar directory.  Full path to the annovar directory, \n" +
				"\n\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/VCFAnnovar -v 9908R.vcf                 \n" +
		        "**************************************************************************************\n");

	}
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		
		File forExtraction = null;
		String annotationsToAdd = null;
		String ethnicity = null;
	
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'd': dbSnpFile = args[++i]; break;
					case 'a': annotationsToAdd = args[++i]; break;
					case 'e': ethnicity = args[++i]; break;
					case 'p': this.pathToAnnovarDir = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//Set up default annotation files
		if (dbSnpFile == null) {
			dbSnpFile = "snp137NonFlagged";
		}
		
		
		//Set up 1K ethnicities
		if (ethnicity == null) {
			this.usedEthnicity = "ALL";
			System.out.println("Using ALL 1K samples for annotation");
		} else {
			if (this.ethnicityMap.containsKey(ethnicity)) {
				this.usedEthnicity = ethnicity;
				System.out.println("Using " + ethnicity + " 1K samples for annotation");
				
			} else {
				System.out.println("Don't recognize the following ethnicity: " + ethnicity + ", exiting");
				System.exit(1);
			}
		}
		
		
		//Set up annotations to test, either full set, or user-specified subset.
		if (annotationsToAdd != null) {
			String[] ata = annotationsToAdd.split(",");
			for (String a: ata) {
				String cleaned = a.trim();
				if (!commandMap.keySet().contains(cleaned)) {
					System.out.println("The application does not recognize the annotation value: " + cleaned + ", skipping");
					continue;
				}
				this.annsToRun.add(cleaned);
			}
		} else {
			//Set up default command, command map keys aren't sorted, so manually set it up
			String[] ctr = new String[]{"ENSEMBL","REFSEQ","DBSNP","ONEK","COSMIC","ESP","SIFT","PP2","MT","LRT","PHYLOP","SEGDUP","GWAS","OMIM"};
			//String[] ctr = new String[]{"ENSEMBL","REFSEQ"};
			this.annsToRun.addAll(Arrays.asList(ctr));
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
	
	
	private class AnnovarCommand {
		private ProcessBuilder command = null;
		private ArrayList<OutputParser> processList = new ArrayList<OutputParser>();
		private String commandName = null;
		
		public AnnovarCommand(String name) {
			this.commandName = name;
		}
		
		public void addCommand(ProcessBuilder command) {
			this.command = command;
		}
		
		public void addOutputParser(OutputParser op) {
			processList.add(op);
		}
		
		public void runCommand(VCFParser parsedVCF) {
			if (command == null || processList.size() == 0) {
				System.out.println("Annovar command is malformed, skipping to next command: " + this.commandName);
				
			} else {
				long startTime = System.currentTimeMillis();
				System.out.println("\n\n***************************************");
				System.out.println("* Starting: " + this.commandName );
				System.out.println("***************************************\n\n");
				try {
					
					command.redirectErrorStream(true);
					Process p = command.start();
					
					try {
						BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
						String line = null;
						while ((line = br.readLine()) != null) {
							System.out.println(line);
						}
						
						int retVal = p.waitFor();
						if (retVal != 0) {
							System.out.println("Annovar command failed, moving to next command: " + this.commandName);
						
						} else {
							for (OutputParser op: processList) {
								op.parseOutput(parsedVCF);
							}
						}
					} catch (IOException ioex ) {
						System.out.println("Failed to read annovar output/error stream : " + ioex.getMessage());
						ioex.printStackTrace();
						System.exit(1);
					} catch (InterruptedException irex) {
						System.out.println("Annovar was interruped before it was finished: " + irex.getMessage());
						irex.printStackTrace();
						System.exit(1);
					} 
					
					
				} catch (IOException ioex) {
					System.out.println("IO error during annovar command execution: " + ioex.getMessage());
					ioex.printStackTrace();
					System.exit(1);
				}
				double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
				System.out.println("\nCommand Done! "+Math.round(diffTime)+" seconds\n");	
			}
		} 
		
		
	}

	private class OutputParser {
		private String[] availableMethodNames = new String[] {"STANDARD","UNSORTED"};
		private HashSet<String> avaiableMethods = new HashSet<String>(Arrays.asList(availableMethodNames));
		
		private File annovarRunOutput;
		private int[] columns;
		private String[] ids;
		private Integer start;
		private String method;
		private String[] infoLine;
		private VCFParser parsedVCF;
		
		/** The argument-less version of this command simply marks the output file for deletion.
		 */
		public OutputParser(String extension) {
			annovarRunOutput = new File(inputname + "." + extension);
			annovarRunOutput.deleteOnExit();
	
		}
		
		/** Create the annovar OutputParse object. This layout is by far the most common.  The annovar output file
		 * has a subset of the original data and the data is contained in just the second column of the output file
		 * @param id		column id
		 * @param infoLine  text describing the new info file to put in the VCF file
		 * @param extension extension of the annovar file
		 */
		public OutputParser(String id, String infoLine, String extension) {
			annovarRunOutput = new File(inputname + "." + extension);
			annovarRunOutput.deleteOnExit();
			
			this.infoLine = new String[] {infoLine};
			this.start = 2;
			this.columns = new int[] {1};
			this.method = "UNSORTED";
			this.ids = new String[] {id};
			
		}
		
		/** Create an annovar output OutputParse object.  This construtor can be used for the more compilcated output files
		 * most notably the refseq or ensembl annotations.
		 * 
		 * @param columns    desired column numbers from the output file.  0-based.
		 * @param ids        column ids.
		 * @param infoLine   text describing the new info field to put in the VCF file
		 * @param start      column number of the 'chromosome' column in the annovar output file.  Used in matching
		 * @param method     match method, options are 'delete', 'unsorted' or 'standard'
		 * @param extension  extension of the annovar output file
		 */
		public OutputParser(int[] columns, String[] ids, String[] infoLine, Integer start, String method, String extension) {
			//Create input/output files
			annovarRunOutput = new File(inputname + "." + extension);
			annovarRunOutput.deleteOnExit();
			
			this.infoLine = infoLine;
			this.start = start;
			this.columns = columns;
			
			if (!this.avaiableMethods.contains(method)) {
				System.out.println("Don't recognize the parsing method, exiting: " + method);
				System.exit(1);
			}
			
			this.method = method;
			
			this.ids = ids;
		}
		
		public void parseOutput(VCFParser parsedVCF) {
			this.parsedVCF = parsedVCF;
			if (this.method == "STANDARD") {
				this.insertInfo();
				this.standardMatch();
			} else if (this.method == "UNSORTED") {
				this.insertInfo();
				this.unsortedMatch();
			}
		}

		
		private void insertInfo() {
			//Add updated info lines to comments section
			int insertIndex=0;
			boolean found = false;
			Pattern p = Pattern.compile("##INFO.*");
		    ArrayList<String> comments = new ArrayList<String>(Arrays.asList(this.parsedVCF.getComments()));
			for (int i=0; i < comments.size(); i++) {
				insertIndex = i;
				Matcher m = p.matcher(comments.get(i));
				if (m.matches()) {
					found = true;
				} else if (found) {
					break;
				}
				
			}
			
			for (String info: this.infoLine) {
				comments.add(insertIndex,info);
				insertIndex++;
			}
			
			this.parsedVCF.setComments(comments.toArray(new String[comments.size()]));
		}
		
		/** The standard match method is for output files that have one line of data for each line of the 
		 * input file. This requires no fancy matching.
		 */
		private void standardMatch() {
			BufferedReader roBR = null;
			
			try {
				roBR = new BufferedReader(new FileReader(annovarRunOutput)); 
				
				for (VCFRecord record: this.parsedVCF.getVcfRecords()) {
					String roLine = roBR.readLine();
					String[] roItems = roLine.split("\t");
					
					for (int i=0; i<this.ids.length; i++) {
						record.getInfoObject().addInfo(this.ids[i], roItems[this.columns[i]].replace(' ','_'));
					}
					
				}
				
				roBR.close();
							
			} catch (IOException ioex) {
				System.out.println("Error reading files: " + ioex.getMessage());
				ioex.printStackTrace();
				System.exit(1);
			}
		}
		
		
		
		/** The unsorted match method matches annovar output files that have sporadic output.  Due to chromosome vs lexographic sorting
		 * issues between VCFs and annovar, I find it easier to just assume the output is unsorted, even though it really is.  As of now, 
		 * the output files aren't big enough to cause real memory issues, so this method should be fine.
		 */
		private void unsortedMatch() {
			HashMap<String,String[]> outputHash = new HashMap<String,String[]>();
			BufferedReader roBR = null;
			try {
				//read in the output file and place the results in a hash
				roBR = new BufferedReader(new FileReader(this.annovarRunOutput));
				
				String roLine = null;
				while((roLine = roBR.readLine()) != null) {
					String[] roItems = roLine.split("\t");
					String index = String.format("%s:%s-%s",roItems[this.start],roItems[this.start+1],roItems[this.start+2]);
					outputHash.put(index,roItems);
				}
				
				roBR.close();
				
				//Add the results to the VCF file
				for (VCFRecord vr: this.parsedVCF.getVcfRecords()) {
					String endPosition = String.valueOf(getAnnovarEndPosition(vr));
					String chrom = vr.getChromosome().replace("chr", "");
					String index = String.format("%s:%s-%s",chrom,vr.getPosition()+1,endPosition);
					if (outputHash.containsKey(index)) {
						for (int i=0; i<this.ids.length; i++) {
							vr.getInfoObject().addInfo(this.ids[i], outputHash.get(index)[this.columns[i]].replace(' ','_'));
						}
					}
				}
				
			} catch (IOException ioex) {
				System.out.println("Error reading files: " + ioex.getMessage());
				ioex.printStackTrace();
				System.exit(1);
			}
			
		}
		
	}


}
