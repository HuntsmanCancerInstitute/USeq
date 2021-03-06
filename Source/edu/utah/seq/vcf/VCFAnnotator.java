package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;



import util.gen.IO;
import util.gen.Misc;

public class VCFAnnotator {
	//Hard coded paths
	
	//private String pathToAnnovarDir = new String("/home/u0855942/annovar_version");
	private String pathToAnnovarDir;
	private String pathToAnnovar;
	private String pathToRepository;
	
	//User specified options
	private File vcfFile;
	private File vaastFile;
	private File vcfOutFile;
	private String pathToTabix = "/tomato/app/tabix/";
	
	//private String dbSnpFile;
	//private String espFile = "esp6500_all";
	//private String cosmicFile= "cosmic63";
	private boolean compressOutput = false;
	private Long genotypesPerChunk = (long)10000000;
	
	private HashMap<String,String> annovarCallHash = new HashMap<String,String>();
	
	
	//Shared variables
	private HashMap<String,AnnovarCommand> commandMap = new HashMap<String,AnnovarCommand>() {{
		put("REFSEQ",new AnnovarCommand("RefSeq Annotations"));
		put("ENSEMBL",new AnnovarCommand("Ensembl Annotations"));
		put("SEGDUP",new AnnovarCommand("Segdup Annotations"));
		put("GWAS",new AnnovarCommand("GWAS Annotations"));
		put("DBSNP",new AnnovarCommand("dbSNP Annotations"));
		put("COSMIC",new AnnovarCommand("COSMIC Annotations"));
		put("ESP",new AnnovarCommand("ESP Annotations"));
		put("ONEK",new AnnovarCommand("1K Genome Annotations"));
		put("SCORES", new AnnovarCommand("Variant scoring Annotation"));
		put("TFBS",new AnnovarCommand("TFBS Annotations"));
		put("DGV",new AnnovarCommand("DGV Annotations"));
		put("OMIM",new AnnovarCommand("OMIM Annotations"));
		//put("V_FLAG",new AnnovarCommand("V_FLAG Annotations"));
		//put("ACMG",new AnnovarCommand("ACMG Annotations"));
		put("NIST",new AnnovarCommand("NIST Annotations"));
		put("CLINVAR",new AnnovarCommand("CLINVAR annotations"));
		put("NCI60",new AnnovarCommand("NCI60 annotations"));
	}};
	
	
	private String usedEthnicity;
	private String inputname;
	private ArrayList<String> annsToRun = new ArrayList<String>();
	
	
	public VCFAnnotator(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		
		//Count vcf records
		System.out.println("Counting VCF samples and records");
		ArrayList<File> tempVcfFiles = new ArrayList<File>();
		int sampleCount = VCFUtilities.countSamples(vcfFile);
		int recordCount = VCFUtilities.countReads(vcfFile);
		//int chunks = recordCount / VCFUtilities.readsToChunk + 1;
		int chunks = (int)((recordCount * sampleCount) / this.genotypesPerChunk + 1);
		int recordsPerChunk = recordCount / chunks;
		
		System.out.println(String.format("Found %d samples and %d records in the vcf file.",sampleCount,recordCount));
		
		if (chunks > 1) {
			System.out.println("File too big to annotate all in one go, splitting into " + chunks + " chunks");
		}
		
		//make tempDirectory
		File tempDir = new File("tempVCF");
		if (tempDir.exists()) {
			IO.deleteDirectory(tempDir);
		}
		tempDir.mkdir();
		
		//Create temporary file names
		String[] fileParts = vcfFile.getName().split("\\.(?=[^\\.]+$)");
		inputname = fileParts[0] + ".annovar.input.txt";
		
		//Setup annovar commands
		this.setupCommands();
		
		for (int i=0;i<chunks;i++) {
			System.out.println("Working on file chunk: " + (i + 1));
			//Create temporary vcf files
			File tempVcf = new File(tempDir,"tempVcf_" + i + ".vcf");
			tempVcfFiles.add(tempVcf);
			
			
			//Create a parsed file per-chunk
			VCFParser parsedVCF = new VCFParser(vcfFile, true, true, true, i, recordsPerChunk);

			
			//Write annovar input file
			this.writeAnnovarInput(parsedVCF);
			
			//Run annovarDriver
			for (String command: this.annsToRun) {
				commandMap.get(command).runCommand(parsedVCF);
			}
			
			//Add VAAST output
			if (vaastFile != null) {
				this.addVaastOutput(parsedVCF);
			}
			
			
			
			//Write VCF file
			parsedVCF.printRecords(tempVcf,true);
			
			parsedVCF = null;
			
		}
		
		//Merge vcf file
		VCFUtilities.mergeVcf(tempVcfFiles, this.vcfOutFile);
		
		//delete temp files
		IO.deleteDirectory(tempDir);
		
		//DeleteLog 
		File logFile = new File(inputname + ".log");
		if (logFile.exists()){ 
			logFile.delete();
		}
		
		//Compress if needed
		if (this.compressOutput) {
			VCFUtilities.createTabix(this.vcfOutFile, this.pathToTabix);
		}
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");	
	}
	
	
	private void addVaastOutput(VCFParser parser) {
		//Parse vaast results
		System.out.println("\n\n***************************************");
		System.out.println("* Parsing VAAST output");
		System.out.println("***************************************\n\n");
		HashMap<String,String[]> vaastHash = new HashMap<String,String[]>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(vaastFile));
			
			ArrayList<String> batch = new ArrayList<String>();
			Pattern snpPat = Pattern.compile("T[UR]:\\s+(.+?)\\s+(\\d+?)@(.+?)\\s+.+");
			Pattern rankPat = Pattern.compile("RANK:(\\d+).*");
			Pattern indelPat = Pattern.compile("T:\\s+(.+?)\\s+(ins|del)-(.+?)-(\\d+?)-\\d+\\s+.+");
			
			
			String temp = null;
			while ((temp = br.readLine()) != null) {
				Matcher snpMatch = snpPat.matcher(temp);
				Matcher rankMatch = rankPat.matcher(temp);
				Matcher indelMatch = indelPat.matcher(temp);
				
				if (snpMatch.matches()) {
					String key = snpMatch.group(3) + ":" + snpMatch.group(2);
					key = key.replace("chr", "");
					
					if (!vaastHash.containsKey(key)) {
						vaastHash.put(key, new String[]{snpMatch.group(1),""});
						batch.add(key);
					}
				} else if (indelMatch.matches()) {
					String key = indelMatch.group(3) + ":" + indelMatch.group(4);
					key = key.replace("chr","");
					
					if (!vaastHash.containsKey(key)) {
						vaastHash.put(key, new String[]{indelMatch.group(1),""});
						batch.add(key);
					}
				} else if (rankMatch.matches()) {
					for (String key: batch) {
						vaastHash.get(key)[1] = rankMatch.group(1);
					}
					batch.clear();
				} 
			}
			br.close();
		} catch (IOException ioex) {
			System.out.println("Could not open vaast output");
			System.exit(1);
		}
		
		//Add VAAST results to VCF
		for (VCFRecord record: parser.getVcfRecords()) {
			String key = record.getChromosome() + ":" + String.valueOf(record.getPosition() + 1 + this.maxInt(record.getReference().length()-1,0));
			key = key.replace("chr","");
			if (vaastHash.containsKey(key)) {
				String[] vaastResult = vaastHash.get(key);
				record.getInfoObject().addInfo("V_LRS", vaastResult[0]);
				record.getInfoObject().addInfo("V_RANK",vaastResult[1]);
			}
		}
		
		//Add Info lines to VCF
		parser.getVcfComments().addInfo("##INFO=<ID=V_LRS,Number=1,Type=String,Description=\"VAAST Likelihood ratio score.  The higher the score "
				+ "the more likely the variant is disease-causing.  Scores are only listed for exonic variants.\">");
		parser.getVcfComments().addInfo("##INFO=<ID=V_RANK,Number=1,Type=String,Description=\"VAAST gene rank.  The lower the rank, the more "
				+ "more likely the gene is causal.  Rank will only be listed next to exonic variants\">");
	}
	
	private void writeAnnovarInput(VCFParser parser) {
		File annovarInput = new File(inputname);
		annovarInput.deleteOnExit();
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(annovarInput));
			
			for(VCFRecord vr: parser.getVcfRecords()) {
				String chrom = vr.getChromosome().replace("chr", "");
				//String chrom = vr.getChromosome();
				String ref = vr.getReference();
				String alt = vr.getAlternate()[0];
				int position = vr.getPosition() + 1;
				String endPosition = String.valueOf(this.getAnnovarEndPosition(vr));
				if (ref.length() == 1 && alt.length() > 1) {
					if (!ref.equals(alt.substring(0, 1))) {
						System.out.println("First base of InDel doesn't match, skipping.  The reference allele should match the first base of the alternate.");
						System.out.println("Reference: " + ref);
						System.out.println("Alternate: " + alt);
						continue;
					} else {
						ref = "-";
						alt = alt.substring(1);
					}
				} else if (ref.length() > 1 && alt.length() == 1) {
					if (!alt.equals(ref.substring(0,1))) {
						System.out.println("First base of InDel doesn't match, skipping. The alternate allele should match the first base of the reference");
						System.out.println("Reference: " + ref);
						System.out.println("Alternate: " + alt);
						continue;
					} else {
						ref = ref.substring(1);
						alt = "-";
						position += 1;
						
					}
					
				}
				
				bw.write(chrom + "\t" + String.valueOf(position) + "\t" + endPosition + "\t" + ref + "\t" + alt + "\n");
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
    	String annovarCall = null;
    	
    	//Annovar gene
    	cmdName = "ENSEMBL";
    	ProcessBuilder pbEnsembl = new ProcessBuilder(this.pathToAnnovar,"--geneanno","--buildver","hg19","-dbtype","ensgene","-splicing_threshold","25","-hgvs",this.inputname,this.pathToRepository);
    	String il1Ensembl = new String("##INFO=<ID=VarType,Number=1,Type=String,Description=\"Exonic variant effect.  If the variant is exonic, this column lists the "
    			+ "functional consequences of the change.  Possible values and (precedence): frameshift insertion (1), frameshift deletion (2), frameshift block "
    			+ "substitution (3), stopgain (4), stoploss (5), nonframeshift insertion (6), nonframeshift deletion (7), nonframeshift block substitution (8), "
    			+ "nonsynonymous SNV (9), synonymous SNV (10), unknown (11)\">");
    	String il2Ensembl = new String("##INFO=<ID=VarDesc,Number=1,Type=String,Description=\"Exonic variant description.  This column lists the gene name, "
    			+ "transcript identifier and sequence change of the corresponding transcript.  If the variant affects multiple transcripts, all will be listed. "
    			+ "If the annotationStyle of the VCF/Report is CLEAN, only the first transcript is reported.  This reduces the size of the field for easier viewing "
    			+ "in IGV\">");
    	String il3Ensembl = new String("##INFO=<ID=EnsemblRegion,Number=1,Type=String,Description=\"Variant location, using ENSEMBL gene definitions.  "
    			+ "Possible values and (precedence): exonic (1), splicing (1), ncRNA (2), UTR5 (3), UTR3 (3), intronic (4), upstream (5), downstream (5), "
    			+ "intergenic (6).\">");
    	String il4Ensembl = new String("##INFO=<ID=EnsemblName,Number=1,Type=String,Description=\"Closest ENSEMBL gene.  If the variant intersects multiple genes, "
    			+ "both will be listed.  If the variant doesn't intersect any gene, the closest upstream and downstream genes are listed with distances.\">");
    	OutputParser op1Ensembl = new OutputParser(new int[]{0,1},new String[]{"EnsemblRegion","EnsemblName"},new String[]{il3Ensembl,il4Ensembl},2,"UNSORTED","variant_function");
    	OutputParser op2Ensembl = new OutputParser(new int[]{1,2},new String[]{"VarType","VarDesc"},new String[]{il1Ensembl,il2Ensembl},3,"UNSORTED","exonic_variant_function");
    	commandMap.get(cmdName).addCommand(pbEnsembl);
    	commandMap.get(cmdName).addOutputParser(op1Ensembl);
    	commandMap.get(cmdName).addOutputParser(op2Ensembl);
    	
    	//Annovar refseq
    	cmdName = "REFSEQ";
    	ProcessBuilder pbRefseq = new ProcessBuilder(this.pathToAnnovar,"--geneanno","--buildver","hg19","-dbtype","gene",this.inputname,this.pathToRepository);
    	String il1Refseq = new String("##INFO=<ID=RefSeq,Number=1,Type=String,Description=\"Closest Refseq gene.  If the variant intersects multiple genes, "
    			+ "both will be listed.  If the variant doesn't intersect any gene, the closest upstream and downstream genes are listed with distances.\">");
    	OutputParser op1Refseq = new OutputParser(new int[]{1},new String[]{"RefSeq"},new String[]{il1Refseq},2,"UNSORTED","variant_function");
    	OutputParser op2Refseq = new OutputParser("exonic_variant_function");
    	commandMap.get(cmdName).addCommand(pbRefseq);
    	commandMap.get(cmdName).addOutputParser(op1Refseq);
    	commandMap.get(cmdName).addOutputParser(op2Refseq);
    	
    	//Annovar TFBS
    	cmdName = "TFBS";
    	ProcessBuilder pbTfbs = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","tfbs",this.inputname,this.pathToRepository);
    	String il1Tfbs = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"Transcrition factor binding sites (TFBS).  If the variant intersects "
    			+ "with a predicted TFBS, the name of the TFBS and the normalized prediction score are listed.  The annotations are the standard TFBS sites used by "
    			+ "annovar.\">");
    	OutputParser op1Tfbs = new OutputParser(cmdName,il1Tfbs,"hg19_tfbsConsSites");
    	commandMap.get(cmdName).addCommand(pbTfbs);
    	commandMap.get(cmdName).addOutputParser(op1Tfbs);
  
    	//Annovar Segdup
    	cmdName = "SEGDUP";
    	ProcessBuilder pbSegdup = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","segdup",this.inputname ,this.pathToRepository);
    	String il1Segdup = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"Segmental duplications (SEGDUP). If the variant intersects with a "
    			+ "known segmental duplication, the identifier and score are listed.  Variants that map to known segmental duplications are often alignment errors and "
    			+ "could be false positives.\">");
    	OutputParser op1Segdup = new OutputParser(cmdName,il1Segdup,"hg19_genomicSuperDups");
    	commandMap.get(cmdName).addCommand(pbSegdup);
    	commandMap.get(cmdName).addOutputParser(op1Segdup);
    	
    	//Annovar DGV
    	cmdName = "DGV";
    	ProcessBuilder pbDgv = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","dgv",this.inputname,this.pathToRepository);
    	String il1Dgv = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"Database of genomic variants (DGV).  If the variant intersects with a "
    			+ "previously reported structural variant listed in DGV, the DGV identifier is listed.  In the case of insertions and deletions, only 1bp of overlap "
    			+ "is required for a successful intersection.\">");
    	OutputParser op1Dgv= new OutputParser(cmdName,il1Dgv,"hg19_dgv");
    	commandMap.get(cmdName).addCommand(pbDgv);
    	commandMap.get(cmdName).addOutputParser(op1Dgv);
    	
    	//Annovar DBSNP
    	cmdName = "DBSNP";
    	
    	annovarCall = getAnnovarCall(cmdName);
    	
    	ProcessBuilder pbSnp = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype", annovarCall,this.inputname,this.pathToRepository);
    	String il1Snp = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"dbSNP identifier, using database: " + annovarCall + ".  If the database "
    			+ "is listed as NonFlagged, dbSNP variants with < 1% minor allele frequency or map only once to the reference assembly are not used in the annotation. \">");
    	OutputParser op1Snp = new OutputParser(cmdName,il1Snp,"hg19_" + annovarCall + "_dropped");
    	OutputParser op2Snp = new OutputParser("hg19_" + annovarCall + "_filtered");
    	commandMap.get(cmdName).addCommand(pbSnp);
    	commandMap.get(cmdName).addOutputParser(op1Snp);
    	commandMap.get(cmdName).addOutputParser(op2Snp);
    	
    	
    	//Annovar mutation scores
    	cmdName = "SCORES";
    	
    	annovarCall = this.annovarCallHash.get(cmdName);
    	
    	ProcessBuilder pbScore = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",annovarCall,this.inputname,this.pathToRepository,"-otherinfo");
    	String il1Sift = new String("##INFO=<ID=SIFT,Number=1,Type=Float,Description=\"SIFT score.  Standard SIFT scores below 0.05 are considered damaging. "
    			+ "If the annovarStyle of the VCF/Report is CLEAN, the SIFT score is reported as 1-SIFT_SCORE, making higher scores more damaging.  This is done "
    			+ "to make the SIFT scores compatible with Polyphen2, MutationTaster and other functional effect predictors.\">");
    	String il1SiftP = new String("##INFO=<ID=SIFT_P,Number=1,Type=String,Description=\"SIFT prediction. Prediction can be (D) daamging or (N) not "
    			+ "damaging.\">");
    	String il1PolyHvar = new String("##INFO=<ID=PP2_HVAR,Number=1,Type=Float,Description=\"Polyphen2 HVAR score. Scores range from 0-1.  Scores larger than "
    			+ "0.909 are considered 'probably damaging', scores between 0.447 and 0.908 are considered 'possibly damaging', scores less than 0.446 are considered "
    			+ "benign.  The PolyPhen HVAR score should be used for diagnositics of Mendelian disease\">");
    	String il1PolyHvarP = new String("##INFO=<ID=PP2_HVAR_P,Number=1,Type=String,Description=\"Polyphen2 HVAR Prediction. Predictions can be 'D' damaging, "
    			+ "'P' possibly damaging or 'B' benign.\">");
    	String il1PolyHidv = new String("##INFO=<ID=PP2_HDIV,Number=1,Type=Float,Description=\"Polyphen2 HDIV score. Scores range from 0-1.  Scores larger than "
    			+ "0.957 are considered 'probably damaging', scores between 0.453 and 0.956 are considered 'possibly damaging', scores less than 0.452 are considered "
    			+ "benign.  The PholyPhen HIDV score should be used for evaluating rare alleles at loci potentially involved in complex phenotypes\">");
    	String il1PolyHidvP = new String("##INFO=<ID=PP2_HDIV_P,Number=1,Type=String,Description=\"Polyphen2 HDIV Prediction.  Predictions can be 'D' damaging, "
    			+ "'P' possibly damaging or 'B' benign\">");
    	String il1Lrt = new String("##INFO=<ID=LRT,Number=1,Type=Float,Description=\"LRT score.  LRT scores is a likelihood ratio test of codon "
    			+ "constraint.  The scores range from 0-1 with higher scores signifying more damaging changes.\">");
    	String il1LrtP = new String("##INFO=<ID=LRT_P,Number=1,Type=String,Description=\"LRT  prediction.  Predictions can be 'D' damaging, 'N' not damaging"
    			+ " or 'U' unknown .\">");
    	String il1Mt = new String("##INFO=<ID=MT,Number=1,Type=Float,Description=\"MutationTaster score.  Higher scores are more damaging.\">");
    	String il1MtP = new String("##INFO=<ID=MT_P,Number=1,Type=String,Description=\"MutationTaster prediction.  Predictions can be 'A' disease_causing_automatic, "
    			+ "'D' disease_causing, 'N' polymorphism or 'P' polymorphism_automatic.\">");
    	String il1Ma = new String("##INFO=<ID=MA,Number=1,Type=Float,Description=\"MutationAssessor score.  Higher scores are more damaging.\">");
    	String il1MaP = new String("##INFO=<ID=MA_P,Number=1,Type=String,Description=\"MutationAssessor prediction.  There are two possible predictions: "
    			+ "predicted functional (high, medium), predicted non-functional (low, neutral)\">");
    	String il1Fathmm = new String("##INFO=<ID=FATHMM,Number=1,Type=Float,Description=\"Fathmm score.  If a score is smaller than -1.5 the corresponding "
    			+ " NS is predicted as D(AMAGING), otherwise it is predicted as T(OLERATED)\">");
    	String il1Gerp = new String("##INFO=<ID=GERP,Number=1,Type=Float,Description=\"GERP++ score.  Generally the higher the score, the more conserved the site.\">");
    	String il1Phylop = new String("##INFO=<ID=PHYLOP,Number=1,Type=Float,Description=\"Phylop score. The PhyloP score is based on multiple alignments of 46 genomes. "
    			+ "The larger the score, the more conserved the site\">");
    	String il1Siphy = new String("##INFO=<ID=SIPHY,Number=1,Type=Float,Description=\"Siphy score. The SiPhy score is based on multiple alignments of 29 mammalian genomes. "
    			+ "The larger the score, the more conserved the site\">");
    	String il1Cadd = new String("##INFO=<ID=CADD,Number=1,Type=Float,Description=\"CADD (Combined Annotation Dependent Depletion) is a score that is based on SVM on multiple other scores. Standard format.\">");
    	String il1CaddP = new String("##INFO=<ID=CADD_P,Number=1,Type=Float,Description=\"CADD (Combined Annotation Dependent Depletion) is a score that is based on SVM on multiple other scores. Phred style score.\">");
    	
    	int[] colLocs = new int[]{0,1,4,5,2,3,6,7,8,9,10,11,22,24,19,20};
    	String[] colNames = new String[]{"SIFT","SIFT_P","PP2_HVAR","PP2_HVAR_P","PP2_HDIV","PP2_HDIV_P","LRT","LRT_P","MT","MT_P","MA","MA_P","PHYLOP","SIPHY","CADD","CADD_P"};
    	String[] colDesc = new String[]{il1Sift,il1SiftP,il1PolyHvar,il1PolyHvarP,il1PolyHidv,il1PolyHidvP,il1Lrt,il1LrtP,il1Mt,il1MtP,il1Mt,il1MtP,il1Ma,il1MaP,il1Phylop,il1Siphy, il1Cadd, il1CaddP};
    	
    	OutputParser opScore1 = new OutputParser(colLocs,colNames,colDesc,2,"hg19_" + annovarCall + "_dropped",1);
    	OutputParser opScore2 = new OutputParser("hg19_" + annovarCall + "_filtered");
    	
    	commandMap.get(cmdName).addCommand(pbScore);
    	commandMap.get(cmdName).addOutputParser(opScore1);
    	commandMap.get(cmdName).addOutputParser(opScore2);
    	

    	//Annovar 1K Genomes
    	cmdName = "ONEK";
    	
    	annovarCall = getAnnovarCall(this.usedEthnicity);
    	String suffix = getAnnovarCall("1K_SUFFIX");
    	
    	ProcessBuilder pbOneK = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",annovarCall,this.inputname,this.pathToRepository);
    	String il1OneK = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=Float,Description=\"1000 Genomes observation Frequency, using database: " + annovarCall +"\">");
    	OutputParser op1OneK = new OutputParser(cmdName,il1OneK,"hg19_" + this.usedEthnicity + ".sites." + suffix + "_dropped");
    	OutputParser op2OneK = new OutputParser("hg19_" + this.usedEthnicity + ".sites." + suffix + "_filtered");
    	commandMap.get(cmdName).addCommand(pbOneK);
    	commandMap.get(cmdName).addOutputParser(op1OneK);
    	commandMap.get(cmdName).addOutputParser(op2OneK);
    	
    	//Annovar COSMIC
    	cmdName = "COSMIC";
    	
    	annovarCall = getAnnovarCall(cmdName);
    	
    	ProcessBuilder pbCosmic = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",annovarCall,this.inputname,this.pathToRepository);
    	String il1Cosmic = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"Catalogue of Somatic Mutations in Cancer (COSMIC) annotations, using "
    			+ "database: " + annovarCall + ".  If the variant intersects with a COSMIC annotation, the variant identifier and number of observations are listed.\">");
    	OutputParser op1Cosmic = new OutputParser(cmdName,il1Cosmic,"hg19_" + annovarCall +"_dropped");
    	OutputParser op2Cosmic = new OutputParser("hg19_" + annovarCall +"_filtered");
    	commandMap.get(cmdName).addCommand(pbCosmic);
    	commandMap.get(cmdName).addOutputParser(op1Cosmic);
    	commandMap.get(cmdName).addOutputParser(op2Cosmic);
    	
    	//Annovar ESP
    	cmdName = "ESP";
    	
    	annovarCall = getAnnovarCall(cmdName);
    	
    	ProcessBuilder pbEsp = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",annovarCall,this.inputname,this.pathToRepository);
    	String il1Esp = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"NHLBI exome sequencing project (ESP) annotations using: " + annovarCall + ". "
    			+ "Variation frequency in the 6500 exomes sequenced in the project.  The samples are from healthy individuals as well as subjects with different diseases.\">");
    	OutputParser op1Esp = new OutputParser(cmdName,il1Esp,"hg19_" + annovarCall + "_dropped");
    	OutputParser op2Esp = new OutputParser("hg19_" + annovarCall + "_filtered");
    	commandMap.get(cmdName).addCommand(pbEsp);
    	commandMap.get(cmdName).addOutputParser(op1Esp);
    	commandMap.get(cmdName).addOutputParser(op2Esp);
    	
    	//Annovar ESP
    	cmdName = "CLINVAR";
    	
    	annovarCall = getAnnovarCall(cmdName);
    	
    	ProcessBuilder pbClinvar = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",annovarCall,this.inputname,this.pathToRepository);
    	String il1Clinvar = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"CLINVAR database with Variant Clinical Significance (unknown, untested, non-pathogenic, "
    			+ "probable-non-pathogenic, probable-pathogenic, pathogenic, drug-response, histocompatibility, other) and Variant disease name.\">");
    	OutputParser op1Clinvar = new OutputParser(cmdName,il1Clinvar,"hg19_" + annovarCall + "_dropped");
    	OutputParser op2Clinvar = new OutputParser("hg19_" + annovarCall + "_filtered");
    	commandMap.get(cmdName).addCommand(pbClinvar);
    	commandMap.get(cmdName).addOutputParser(op1Clinvar);
    	commandMap.get(cmdName).addOutputParser(op2Clinvar);
    	
    	//Annovar ESP
    	cmdName = "NCI60";
    	
    	annovarCall = getAnnovarCall(cmdName);
    	
    	ProcessBuilder pbNci60 = new ProcessBuilder(this.pathToAnnovar,"--filter","--buildver","hg19","-dbtype",annovarCall,this.inputname,this.pathToRepository);
    	String il1Nci60 = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"NCI-60 human tumor cell line panel exome sequencing allele frequency data.\">");
    	OutputParser op1Nci60 = new OutputParser(cmdName,il1Nci60,"hg19_" + annovarCall + "_dropped");
    	OutputParser op2Nci60 = new OutputParser("hg19_" + annovarCall + "_filtered");
    	commandMap.get(cmdName).addCommand(pbNci60);
    	commandMap.get(cmdName).addOutputParser(op1Nci60);
    	commandMap.get(cmdName).addOutputParser(op2Nci60);
    	
    	//Annovar GWAS
    	cmdName = "GWAS";
    	ProcessBuilder pbGwas = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-dbtype","gwascatalog",this.inputname,this.pathToRepository);
    	String il1Gwas = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"GWAS catalog associations.  If the variant intersects with a GWAS catalog listing, "
    			+ "the score and associated disease are listed.  The GWAS calalog stores SNP-trait associations with p-values < 1.0x10^5. This annotation is not comprehensive.\">");
    	OutputParser op1Gwas = new OutputParser(cmdName,il1Gwas,"hg19_gwasCatalog");
    	commandMap.get(cmdName).addCommand(pbGwas);
    	commandMap.get(cmdName).addOutputParser(op1Gwas);
    	
    	//Annovar OMIM
    	cmdName = "OMIM";
    	ProcessBuilder pbOmim = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-gff3attrib","-dbtype","gff3","-gff3dbfile","hg19_omim.gff3",this.inputname,this.pathToRepository);
    	String il1Omim = new String("##INFO=<ID=" + cmdName + ",Number=1,Type=String,Description=\"Online Mendelian Inheritance in Man (OMIM) annotations.  If the variant intersects "
    			+ "a OMIM gene, the gene identifier and associated diseases are listed.\">");
    	OutputParser op1Omim = new OutputParser(cmdName,il1Omim,"hg19_gff3");
    	commandMap.get(cmdName).addCommand(pbOmim);
    	commandMap.get(cmdName).addOutputParser(op1Omim);
    	
//    	//Annovar ACMG
//    	cmdName = "ACMG";
//    	ProcessBuilder pbAcmg = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-gff3attrib","-dbtype","gff3","-gff3dbfile","hg19_acmg.gff3",this.inputname,this.pathToRepository);
//    	String il1Acmg = new String("##INFO=<ID=" + cmdName + ",Number=0,Type=Flag,Description=\"American College of Medical Genetics and Genomics (ACMG) annotations.  ACMG recommends reporting variants "
//    			+ "with this flag, since they are contained within clinically relevant genes.\">");
//    	OutputParser op1Acmg = new OutputParser(cmdName,il1Acmg,"hg19_gff3",true);
//    	commandMap.get(cmdName).addCommand(pbAcmg);
//    	commandMap.get(cmdName).addOutputParser(op1Acmg);
//    	
//    	//Annovar VAAST
//    	cmdName = "V_FLAG";
//    	ProcessBuilder pbVflag = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-gff3attrib","-dbtype","gff3","-gff3dbfile","hg19_vflag.gff3",this.inputname,this.pathToRepository);
//    	String il1Vflag = new String("##INFO=<ID=" + cmdName + ",Number=0,Type=Flag,Description=\"Incendentalome annotations. Variants marked with this flag are within genes that "
//    			+ "commonly come up as positive in VAAST runs.\">");
//    	OutputParser op1Vflag = new OutputParser(cmdName,il1Vflag,"hg19_gff3",true);
//    	commandMap.get(cmdName).addCommand(pbVflag);
//    	commandMap.get(cmdName).addOutputParser(op1Vflag);
    	
    	//Annovar NIST
    	cmdName = "NIST";
    	ProcessBuilder pbNist = new ProcessBuilder(this.pathToAnnovar,"--regionanno","--buildver","hg19","-gff3attrib","-dbtype","gff3","-gff3dbfile","hg19_nist.gff3",this.inputname,this.pathToRepository);
    	String il1Nist = new String("##INFO=<ID=" + cmdName + ",Number=0,Type=Flag,Description=\"Regions that can be resolved with high certainty in the female HapMap individual NA12878 according to NIST.\">");
    	OutputParser op1Nist = new OutputParser(cmdName,il1Nist,"hg19_gff3",true);
    	commandMap.get(cmdName).addCommand(pbNist);
    	commandMap.get(cmdName).addOutputParser(op1Nist);
 
    }
    
    private String getAnnovarCall(String cmdName) {
    	if (!this.annovarCallHash.containsKey(cmdName)) {
    		System.out.println(String.format("This software expects data for %s in the key.txt file in: %s.",cmdName, this.pathToRepository));
    		System.exit(1);
    	}
    	return this.annovarCallHash.get(cmdName);
    }
	

	public static void main(String[] args) {
		new VCFAnnotator(args);
	}
	
	private void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Annotator : November 2014                       **\n" +
				"**************************************************************************************\n" +
				"VCFAnnotator adds user-specifed annotations to the VCF file INFO field.  Only hg19 is \n" +
				"supported at this time.  If your VCF file has more than 500,000 records, it will be \n" +
				"split into smaller VCF files that are annotated separately.  Once annotation is \n" +
				"complete, the individual annotated files are merged and compressed.  This application \n" +
				"uses a lot of memory when running large VCF files, so use 20gb of memory when starting\n" +
				"java.\n\n" +
		
				"Required:\n" +
				"-v VCF file. Path to a multi-sample vcf file, compressed ok (XXX.vcf/XXX.vcf.gz).\n"+
				"-o Output VCF file.  Path to the annotated vcf file, can be specifed as XXX.vcf or \n" + 
				"   XXX.vcf.gz. If XXX.vcf.gz, the file will be compressed and indexed using tabix.\n" +
				"-p Path to annovar directory.\n\n" +
				
				"Optional:\n"+
				"-d dbSNP database.  By default, this application uses dbSNP 137 for annotation. Use \n" +
				"      this option along with dbSNP database identifier to use a different version, \n" +
				"      i.e. snp129. The annovar-formatted dbSNP database must be in the annovar data \n" +
				"      directory for this option to work.\n" +
				"-e Ethnicity.  By default, the 1K frequency is calculated across all ethnicities.  \n" +
				"      If you want to restrict it to one of African (AFR), Admixed American (AMR) \n" + 
				"      East Asian (EAS), European (EUR), South Asian (SAS), use this option \n" +
				"      followed by the ethnicity identifier.\n" +
				"-a Annotations to add.  By default, this application uses a subset of annovar  \n" +
				"      annotations.  Use a comma-separated list of keys to specify a custom set. \n" + 
				"      Available annotations with (keys): ensembl gene annotations (ENSEMBL), refSeq \n" +
				"      gene names (REFSEQ), transcription factor binding sites (TFBS), segmental \n" +
				"      duplicatons (SEGDUP), database of genomic variants (DGV), variant scores \n" + 
				"      (SCORES), GWAS catalog annotations (GWAS), clinvar variations (CLINVAR), \n" +
				"      dbsnp annotations (DBSNP), nci60 annotations (NCI60), 1K genomes annotations \n" + 
				"      (ONEK), COSMIC annotations (COSMIC), ESP annotations (ESP), OMIM genes and \n" + 
				"      diseases (OMIM), and NIST callable ragions (NIST).  The SCORES option \n " + 
				"      includes SIFT, PolyPhen2, MutationTaster, MutationAssessor, LRT, PhyloP and \n" +
				"      SiPhy. The ENSEMBL option includes the columns EnsemblRegion, EnsemblName, \n" +
				"      VarType and VarDesc. \n" +
				"-r Path to annovar repository.  By default, VCFAnnotator uses the humandb directory in \n" +
				"      the annovar program directory.  Use this if you want to specify a different \n " +
				"      directory\n " +
				"-n VAAST output.  If a VAAST output file is specified, the VCF file is annotated with\n" +
				"      the VAAST variation score and gene rank.\n" +
				"-t Path to tabix directory.\n" +
				"-c Number of genotypes to process per chunk.  The total number of genotypes is the \n" +
				"      number of records in the VCF * the number of samples in the VCF.  By default \n" +
				"      the file is split into 10M genotypes.  Can take 20GB or more to run.\n" +
				"\n\n"+

				"Example: java -Xmx20G -jar pathTo/USeq/Apps/VCFAnnotator -v 9908R.vcf \n" + 
				"      -o 9908_ann.vcf.gz \n\n" +
		        "**************************************************************************************\n");

	}
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		
		File inputFile = null;
		String annotationsToAdd = null;
		String ethnicity = null;
		String vaastName = null;
		String outputFile = null;
	
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': inputFile = new File(args[++i]); break;
					case 'd': this.annovarCallHash.put("DBSNP", args[++i]); break;
					case 'a': annotationsToAdd = args[++i]; break;
					case 'e': ethnicity = args[++i]; break;
					case 'p': this.pathToAnnovarDir = args[++i]; break;
					case 'r': this.pathToRepository = args[++i]; break;
					case 't': this.pathToTabix = args[++i]; break;
					case 'n': vaastName = args[++i]; break;
					case 'o': outputFile = args[++i]; break;
					case 'c': this.genotypesPerChunk = Long.parseLong(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (pathToAnnovarDir == null) {
			System.out.println("The path to annovar (-p) must be specified. \n");
			System.exit(1);
		}
	
		if (vaastName != null) {
			vaastFile = new File(vaastName);
			if (!vaastFile.exists()) {
				System.out.println("VAAST output file could not be found, skipping");
			}
		}
		
		//Make sure files are present and accounted for.
		pathToAnnovar = new String(pathToAnnovarDir + "/annotate_variation.pl");
		if (!(new File(pathToAnnovar).exists())) {
			System.out.println(String.format("Path to annovar does not exist: %s,  exiting",pathToAnnovar));
			System.exit(1);
		}
		pathToRepository = new String(pathToAnnovarDir + "/humandb");
		File repoFile = new File(pathToRepository);
		if (!repoFile.exists()) {
			System.out.println(String.format("Can't find the annovar data repository: %s, exiting",pathToRepository));
			System.exit(1);
		}
		File indexFile = new File(repoFile,"key.txt");
		if (!indexFile.exists()) {
			System.out.println(String.format("Can't find the annovar key file ( %s ) in the data directory: %s, exiting",indexFile.getAbsolutePath(),repoFile.getAbsolutePath()));
			System.exit(1);
		}
		
		//Read in key file
		try {
			BufferedReader br = new BufferedReader(new FileReader(indexFile));
			String temp = "";
			while((temp = br.readLine()) != null) {
				String[] items = temp.split("=");
				if (!this.annovarCallHash.containsKey(items[0])) {
					this.annovarCallHash.put(items[0], items[1]);
				}	
			}
			br.close();
		} catch (IOException ioex) {
			System.out.println(String.format("Could not read key file: %s ", ioex.getMessage()));
			System.exit(1);
		}
		

		//Set up 1K ethnicities
		if (ethnicity == null) {
			this.usedEthnicity = "ALL";
			System.out.println("Using ALL 1K samples for annotation");
		} else { 
			if (this.annovarCallHash.containsKey(ethnicity)) {
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
			String[] ctr = new String[]{"ENSEMBL","REFSEQ","DBSNP","ONEK","ESP","COSMIC","NCI60","SCORES","CLINVAR","SEGDUP","GWAS","OMIM","NIST"};
			//String[] ctr = new String[]{"ENSEMBL","REFSEQ"};
			this.annsToRun.addAll(Arrays.asList(ctr));
		}
		
		if (inputFile == null) {
			System.out.println("Input file was not specified, exiting");
			System.exit(1);
		} else if (!inputFile.exists()) {
			System.out.println("Input file does not exist, exiting");
			System.exit(1);
		} else if (Pattern.matches(".+?.vcf",inputFile.getName())) {
			vcfFile = inputFile;
		} else if (Pattern.matches(".+?.vcf.gz",inputFile.getName())) {
			System.out.println("Decompressing VCF file");
			vcfFile = VCFUtilities.unzipTabix(inputFile,this.pathToTabix);
		} else {
			System.out.println("Input file does not appear to be a XXX.vcf/XXX.vcf.gz file");
			System.exit(1);
		}
		
		
		
		if (outputFile == null) {
			System.out.println("Output file was no specified, exiting");
			System.exit(1);
		} else if (Pattern.matches(".+?.vcf",outputFile)) {
			this.compressOutput = false;
			this.vcfOutFile = new File(outputFile);
		} else if (Pattern.matches(".+?.vcf.gz",outputFile)) {
			File vcfOutComp = new File(outputFile);
			if (vcfOutComp.exists()) {
				vcfOutComp.delete();
				//System.out.println("Tabix won't overwrite an existing file, rename the output or delete exisiting file");
				//System.exit(1);
			}
			this.vcfOutFile = new File(outputFile.substring(0,outputFile.length()-3));
			this.compressOutput = true;
		} else {
			System.out.println("Output file does not appear to be a XXX.vcf/XXX.vcf.gz file");
			System.exit(1);
		}
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
		private String[] availableMethodNames = new String[] {"STANDARD","UNSORTED","FLAG"};
		private HashSet<String> avaiableMethods = new HashSet<String>(Arrays.asList(availableMethodNames));
		
		private File annovarRunOutput;
		private int[] columns;
		private String[] ids;
		private Integer start;
		private String method;
		private String[] infoLine;
		private VCFParser parsedVCF;
		private int column;
		
		/** The argumentless version of this command simply marks the output file for deletion.
		 */
		public OutputParser(String extension) {
			annovarRunOutput = new File(inputname + "." + extension);
			annovarRunOutput.deleteOnExit();
	
		}
		
		/** Create the annovar OutputParse object. This layout is by far the most common.  The annovar output file
		 * has a subset of the original data and the data is contained in just the second column of the output file.
		 * @param id		column id
		 * @param infoLine  text describing the new info file to put in the VCF file
		 * @param extension extension of the annovar file
		 */
		public OutputParser(String id, String infoLine, String extension) {
			this(id,infoLine,extension,false);
			
		}
		
		/** Create the annovar OutputParse object. This layout is by far the most common.  The annovar output file
		 * has a subset of the original data and the data is contained in just the second column of the output file.
		 * If just a yes/no for presence or absence needed, set flag = true
		 * @param id		column id
		 * @param infoLine  text describing the new info file to put in the VCF file
		 * @param extension extension of the annovar file
		 * @param flag	    set to true if you just want a flag VCF field
		 */
		public OutputParser(String id, String infoLine, String extension, boolean flag) {
			annovarRunOutput = new File(inputname + "." + extension);
			annovarRunOutput.deleteOnExit();
			
			this.infoLine = new String[] {infoLine};
			this.start = 2;
			this.columns = new int[] {1};
			if (flag) {
				this.method = "FLAG";
			} else {
				this.method = "UNSORTED";
			}
			
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
		
		
		/** Yet another constructor.  This is used if the output column contains comma-delimited data.
		 
		 * 
		 * @param columns    desired column numbers from the output file.  0-based.
		 * @param ids        column ids.
		 * @param infoLine   text describing the new info field to put in the VCF file
		 * @param start      column number of the 'chromosome' column in the annovar output file.  Used in matching
		 * @param method     match method, options are 'delete', 'unsorted' or 'standard'
		 * @param extension  extension of the annovar output file
		 */
		public OutputParser(int[] columns, String[] ids, String[] infoLine, Integer start, String extension, int column) {
			//Create input/output files
			annovarRunOutput = new File(inputname + "." + extension);
			annovarRunOutput.deleteOnExit();
			
			this.infoLine = infoLine;
			this.start = start;
			this.columns = columns;
			this.column = column;
			this.ids = ids;
			
			this.method = "MULTIVALUE";
		}
		
		public void parseOutput(VCFParser parsedVCF) {
			this.parsedVCF = parsedVCF;
			if (this.method == "STANDARD") {
				this.insertInfo();
				this.standardMatch();
			} else if (this.method == "UNSORTED") {
				this.insertInfo();
				this.unsortedMatch(false,false);
			} else if (this.method == "FLAG") {
				this.insertInfo();
				this.unsortedMatch(true,false);
			} else if (this.method == "MULTIVALUE") {
				this.insertInfo();
				this.unsortedMatch(false, true);
			}
		}

		
		private void insertInfo() {
			for (String info: this.infoLine) {
				this.parsedVCF.getVcfComments().addInfo(info);
			}
			
			
			
		
		}
		
		
		/** The standard match method is for output files that have one line of data for each line of the 
		 * input file. This requires no fancy matching.
		 */
		private void standardMatch() {
			BufferedReader roBR = null;
			
			try {
				roBR = new BufferedReader(new FileReader(annovarRunOutput));
				
				//Make sure the line counts match.
				int count=0;
				String temp = null;
			    while ((temp = roBR.readLine()) != null) {
			    	count +=1;
			    } 
			    
			    if (count != this.parsedVCF.getVcfRecords().length) {
			    	System.out.println("The number of annovar annotations doesn't equal the number of VCF records."
			    			+ "  This tends to be indicative of malformed annovar input, so you will need the fix "
			    			+ "the errors mention in the annovar invalid_input file.\n");
			    }
			    roBR.close();
			    
			    //parse the output
			    roBR = new BufferedReader(new FileReader(annovarRunOutput));
				
				for (VCFRecord record: this.parsedVCF.getVcfRecords()) {
					String roLine = roBR.readLine();
					String[] roItems = roLine.split("\t");
					
					for (int i=0; i<this.ids.length; i++) {
						
						record.getInfoObject().addInfo(this.ids[i], (roItems[this.columns[i]].replace(' ','_')).replace(';',','));
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
		private void unsortedMatch(boolean flag, boolean multicolumn) {
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
					String chrom = vr.getChromosome().replace("chr", "");
					//String chrom = vr.getChromosome();
					String ref = vr.getReference();
					String alt = vr.getAlternate()[0];
					int position = vr.getPosition() + 1;
					String endPosition = String.valueOf(getAnnovarEndPosition(vr));
					if (ref.length() == 1 && alt.length() > 1) {
						if (!ref.equals(alt.substring(0, 1))) {
							System.out.println("First base of InDel doesn't match, skipping");
						} else {
							ref = "-";
							alt = alt.substring(1);
						}
					} else if (ref.length() > 1 && alt.length() == 1) {
						if (!alt.equals(ref.substring(0,1))) {
							System.out.println("First base of InDel doesn't match, skipping");
						} else {
							ref = ref.substring(1);
							alt = "-";
							position += 1;
							
						}
						
					}
					
					String index = String.format("%s:%s-%s",chrom,position,endPosition);
					
					if (outputHash.containsKey(index)) {
						for (int i=0; i<this.ids.length; i++) {
							if (flag) {
								vr.getInfoObject().addInfo(this.ids[i],"true");
							} else if (multicolumn) {
								String[] values = outputHash.get(index)[this.column].split(",");
								String value = values[this.columns[i]];
								if (!value.equals(".")) {
									vr.getInfoObject().addInfo(this.ids[i], values[this.columns[i]]);
								}
							} else {
								vr.getInfoObject().addInfo(this.ids[i], (outputHash.get(index)[this.columns[i]].replace(' ','_')).replace(';',','));
							}
							
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
