package edu.utah.seq.run.invitae;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.opencsv.CSVReader;
import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

/**Takes the two germline Invitae batch export files containing 'patients' and 'variants', uses the PHI to find/create a PMR ID, vcf, and json results
 Run sed 's/Hg38End/END/g' RQ200422_AR182965_Hg38.vcf > fixed.vcf after doing the crossmap to convert those to END
 */
public class InvitaeBulkCsvParser {

	//user defined fields
	private File jobsDirectory = null;
	private File patientsFile = null;
	private File variantsFile = null;
	private File smmDirectory = null;
	private File smmRegistryDirectory = null;
	private File indexedFasta = null;
	private File hg38Genes = null;

	//internal
	private LinkedHashMap<String, Integer> patientHeaderIndex = null;
	private LinkedHashMap<String, Integer> variantHeaderIndex = null;
	private HashMap<String, InvitaeReport> reports = new HashMap<String, InvitaeReport>();
	private LinkedHashMap<String, String> phiPmrId = null;
	private IndexedFastaSequenceFile fasta = null;
	
	//InvitaeGeneName: InvitaeGeneName Chrom Start Stop HUGOName Strand
	private HashMap<String, Bed> iGenCoorHg38 = null;
	
	//variant indexes, leave as package visibility
	int varInterpretation = -1;
	int varGene = -1;
	int varMethod = -1;
	int varTranscript = -1;
	int varCDot = -1;
	int varPDot = -1;
	int varZygosity = -1;
	int varChromosome = -1;
	int varPos = -1;
	int varRef = -1;
	int varAlt = -1;
	
	//for the vcf header
	static String altDel = "##ALT=<ID=DEL,Description=\"Deletion - decrease in copies relative to reference\">\n";
	static String altDup = "##ALT=<ID=DUP,Description=\"Duplication - increase in copies relative to reference\">\n";
	static String infoGene = "##INFO=<ID=iGene,Number=1,Type=String,Description=\"Invitae gene name\">\n";
	static String infoTranscript = "##INFO=<ID=iTrans,Number=1,Type=String,Description=\"Invitae gene transcript\">\n";
	static String infoCDot = "##INFO=<ID=iCDot,Number=1,Type=String,Description=\"Invitae cDot\">\n";
	static String infoPDot = "##INFO=<ID=iPDot,Number=1,Type=String,Description=\"Invitae pDot\">\n";
	static String infoZygosity = "##INFO=<ID=iZyg,Number=1,Type=String,Description=\"Invitae variant zygosity\">\n";
	static String infoInterpretation = "##INFO=<ID=iInterp,Number=1,Type=String,Description=\"Invitae variant effect\">\n";
	static String infoMethod = "##INFO=<ID=iMeth,Number=1,Type=String,Description=\"Invitae variant detection method\">\n";
	//must avoid END all caps since crossmap will try to convert it
	static String infoEnd = "##INFO=<ID=Hg38End,Number=1,Type=Integer,Description=\"Hg38 end position of the DUP or DEL affected gene\">\n";
	private String vcfHeader = null;
	private String bedHeader = null;
	
	public InvitaeBulkCsvParser (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);

			parsePatientsFile();
			
			parseVariantsFile();
			setVariantIndexes();
			parseVariantsInReports();
			
			collectUniquePatients();
			
			findFetchPmrIds();
			
			assignPmrsToReportsAddJsonVcf();
			
			//delete tmp dir with PHI
			IO.deleteDirectory(smmDirectory);

			//finish and calc run time
			fasta.close();
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			IO.el("\nERROR running the Invitae Data Wrangler, aborting");
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void assignPmrsToReportsAddJsonVcf() throws IOException {
		Integer indexOrderingPhysician = patientHeaderIndex.get("ordering_clinician");
		if (indexOrderingPhysician == null) throw new IOException("Failed to find the 'ordering_clinician' in "+patientHeaderIndex);
		Integer indexReportDate = patientHeaderIndex.get("report_date");
		if (indexReportDate == null) throw new IOException("Failed to find the 'report_date' in "+patientHeaderIndex);
		String[] headerNames = Misc.setToStringArray(patientHeaderIndex.keySet());
		
		for (InvitaeReport rep: reports.values()) {
			rep.setPmrId(phiPmrId.get(rep.getSubjectMatchMakerLine()));
			rep.addClinicalReport(jobsDirectory, indexOrderingPhysician, indexReportDate, headerNames);
			rep.addDataFiles(vcfHeader, bedHeader);
		}
	}

	private void findFetchPmrIds() throws IOException {
		
		//write out the search file
		File queryFile = new File (smmDirectory, "invitaePatientQueries_PHI.txt");
		File queryResDir = new File (smmDirectory, "SMMQueryResults_PHI");
		queryResDir.mkdirs();
		PrintWriter out = new PrintWriter (new FileWriter( queryFile));
		for (String line: phiPmrId.keySet()) out.println(line);
		out.close();

		//run the SMM
		String[] args = {
				"-r", smmRegistryDirectory.getCanonicalPath(),
				"-q", queryFile.getCanonicalPath(),
				"-o", queryResDir.getCanonicalPath(),
				"-v", //turn off verbosity
				"-a", //adding queries not found to registry
				"-c", //making name matching case-insensitive for EDW
				"-u",  //fill in info not found to registry entries from queries, e.g. missing mrns, otherIds
				"-t 4" //four threads, doesn't seem to limit it though?
		};
		SubjectMatchMaker smm = new SubjectMatchMaker(args);
		Subject[] queries = smm.getQuerySubjects();
		//the order of the output matches the order of the input
		if (queries.length != phiPmrId.size()) throw new IOException("ERROR: the SMM query output size "+queries.length+" doesn't match the input "+phiPmrId.size());

		Iterator<String> it = phiPmrId.keySet().iterator();
		for (int i=0; i< queries.length; i++) {
			String input = it.next();
			//IO.pl(input);
			phiPmrId.put(input, queries[i].getCoreIdNewOrMatch());
			//IO.pl(queries[i].toString()+" new? "+queries[i].isCoreIdCreated()+" topFound? "+queries[i].isTopMatchFound()+" coreId: "+queries[i].getCoreIdNewOrMatch());
		}
		//IO.pl(phiPmrId);
	}




	private void collectUniquePatients() throws IOException {
		IO.pl("Collecting unique patient info...");
		//what PHI to collect
		Integer lastName = patientHeaderIndex.get("patient_last_name");
		Integer firstName = patientHeaderIndex.get("patient_first_name");
		Integer dob = patientHeaderIndex.get("dob");
		Integer gender = patientHeaderIndex.get("gender");
		Integer mrn = patientHeaderIndex.get("mrn");
		if (lastName==null || firstName==null || dob==null || gender==null || mrn==null) throw new IOException("ERROR: failed to find one or more of : patient_last_name, "
				+ "patient_first_name, dob, gender, mrn from "+patientHeaderIndex);
		
		//for each report, fetch PHI for SMM
		phiPmrId = new LinkedHashMap<String, String>();
		for (InvitaeReport ir: reports.values()) {
			ir.setSubjectMatchMakerLine(lastName, firstName, dob, gender, mrn);
			phiPmrId.put(ir.getSubjectMatchMakerLine(),null);
		}
	}

	private void parseVariantsFile() throws Exception {

		CSVReader reader = new CSVReader(new FileReader(variantsFile));
		String[]cells = null;
		int reportIdIndex = -1;
		int numVars = 0;
		while ((cells = reader.readNext()) != null) {

			if (cells[0].trim().equals("requisition_id")) {
				IO.p("\tParsing variant header line...\n\t");
				variantHeaderIndex = new LinkedHashMap<String, Integer>();
				for (int i=0; i< cells.length; i++) variantHeaderIndex.put(cells[i], i);
				reportIdIndex = variantHeaderIndex.get("report_id");
			}
			else {
				//ok this is a data line, has it been seen before
				String reportId = cells[reportIdIndex];
				IO.p(" "+reportId);
				InvitaeReport ir = reports.get(reportId);
				if (ir == null) throw new Exception("\nERROR: failed to find a report for the variant : "+reportId);
				ir.getVariantInfo().add(cells);
				numVars++;
			}
		}
		reader.close();
		IO.pl("\n\t\t"+numVars+"\tVariants");
	}
	
	private void parseVariantsInReports() throws IOException {
		IO.pl("\tParsing variants...");
		for (InvitaeReport ir: reports.values()) {
			ir.parseVariants(this);
		}
	}

	private void setVariantIndexes() throws IOException {
		IO.pl("\n\tParsing variant indexes...");
		varInterpretation = fetchVariantIndex("interpretation");
		varGene = fetchVariantIndex("gene");
		varMethod = fetchVariantIndex("method");
		varTranscript = fetchVariantIndex("transcript");
		varCDot = fetchVariantIndex("cdot");
		varPDot = fetchVariantIndex("pdot");
		varZygosity = fetchVariantIndex("zygosity");
		varChromosome = fetchVariantIndex("chromosome");
		varPos = fetchVariantIndex("vcf_pos");
		varRef = fetchVariantIndex("vcf_ref");
		varAlt = fetchVariantIndex("vcf_alt");
	}

	private int fetchVariantIndex(String key) throws IOException {
		Integer i = variantHeaderIndex.get(key);
		if (i==null) throw new IOException("ERROR: failed to find "+i+" in the Variant Header "+variantHeaderIndex.keySet());
		return i.intValue();
	}

	private void parsePatientsFile() throws Exception {

		CSVReader reader = new CSVReader(new FileReader(patientsFile));
		String[]cells = null;
		int reportIdIndex = -1;
		int requisitionIdIndex = -1;
		int testedGenesIdIndex = -1;
		while ((cells = reader.readNext()) != null) {

			if (cells[0].trim().equals("account_name")) {
				IO.pl("\tParsing patient header line...");
				patientHeaderIndex = new LinkedHashMap<String, Integer>();
				for (int i=0; i< cells.length; i++) patientHeaderIndex.put(cells[i], i);
				reportIdIndex = patientHeaderIndex.get("report_id");
				requisitionIdIndex = patientHeaderIndex.get("requisition_id");
				testedGenesIdIndex = patientHeaderIndex.get("genes_analyzed");
				IO.p("\t");
			}
			else {
				//ok this is a data line, has it been seen before
				String reportId = cells[reportIdIndex];
				String requisitionId = cells[requisitionIdIndex];
				IO.p(" "+reportId+"_"+requisitionId);
				InvitaeReport ir = reports.get(reportId);
				if (ir != null) throw new Exception("\nERROR: found duplicate report id in the patient file: "+reportId);
				ir = new InvitaeReport(reportId, requisitionId, cells[testedGenesIdIndex], iGenCoorHg38, cells);
				reports.put(reportId, ir);
			}
		}
		reader.close();
		IO.pl("\n\t\t"+reports.size()+"\tReports\n");

	}





	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new InvitaeBulkCsvParser(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{

		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-zA-Z]");
		File tmpDir = null;
		for (int i = 0; i<args.length; i++){
			Matcher mat = pat.matcher(args[i]);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': patientsFile = new File(args[++i]); break;
					case 'v': variantsFile = new File(args[++i]); break;
					case 'j': jobsDirectory = new File(args[++i]); break;
					case 'f': indexedFasta = new File( args[++i]); break;
					case 'g': hg38Genes = new File( args[++i]); break;
					case 't': tmpDir = new File(args[++i]); break;
					case 'c': smmRegistryDirectory = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//check for files
		if (patientsFile == null || patientsFile.exists()== false) Misc.printErrAndExit("Error: cannot find your patients file? "+patientsFile);
		if (variantsFile == null || variantsFile.exists()== false) Misc.printErrAndExit("Error: cannot find your variants file? "+variantsFile);
		if (indexedFasta == null || indexedFasta.exists()== false) Misc.printErrAndExit("Error: cannot find your hg19 indexed fasta file? "+indexedFasta);
		if (hg38Genes == null || hg38Genes.exists()== false) Misc.printErrAndExit("Error: cannot find your hg38 interrogated gene coordinates file? "+hg38Genes);

		//check jobDirectory
		if (jobsDirectory !=null) jobsDirectory.mkdirs();
		if (jobsDirectory == null || jobsDirectory.exists() == false) Misc.printErrAndExit("Error: cannot find or make your jobs directory "+jobsDirectory);

		//check smm registry directory
		if (smmRegistryDirectory == null || smmRegistryDirectory.exists() == false || smmRegistryDirectory.isDirectory()==false) Misc.printErrAndExit("Error: cannot find your Subject Match Maker registry directory, "+smmRegistryDirectory);

		//check tmpDir for the subject match maker
		if (tmpDir == null) Misc.printErrAndExit("Error: cannot find your temp subject ID directory "+tmpDir);
		tmpDir.mkdirs();
		if (tmpDir.exists() == false) Misc.printErrAndExit("Error: cannot make your temp subject ID directory "+tmpDir);
		smmDirectory = new File (tmpDir, System.currentTimeMillis()+"_DeleteMe");
		smmDirectory.mkdir();
		
		//load the interrogated gene coordinates hg38 for dup and del end coordinates and the bed file of covered genes for each assay.
		loadInvitaeGeneCoor();
		
		//create fasta fetcher
		fasta = new IndexedFastaSequenceFile(indexedFasta);
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);
		
		//create headers
		buildVcfHeader();
		buildBedHeader();

	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Invitae Bulk Csv Parser : Oct 2023                     **\n" +
				"**************************************************************************************\n" +
				"The IBCP parses the two Invitae bulk export patient and variant spreadsheets into\n"+
				"Patient Molecular Repository (PMR) compatible vcf, bed, and json files. It uses the\n"+
				"patient PHI to fetch or make PMR Ids using the SubjectMatchMaker\n"+
				"app and assembles PMR directories. For each patient's test, it writes out a bed file\n"+
				"of tested genes, a json file of info from the patient cvs, and a vcf file of variants\n"+
				"including <DUP> <DEL> copy number. Run vt normalization and CrossMap to convert the\n"+
				"hg19 vcf coordinates to hg38.\n\n"+

				"Options:\n"+
				"-p Invitae patients csv file\n"+
				"-v Invitae variants csv file\n" +
				"-f Hg19 indexed fasta\n"+
				"-t Directory to place temp files with PHI for Subject ID matching\n"+
				"-c Directory containing the SubjectMatchMaker 'currentRegistry_' file. This will be\n"+
				"     updated with new patient IDs.\n"+
				"-g Hg38 all Invitae interrogated gene coordinate file.\n\n"+

				"Example: java -jar -Xmx2G pathToUSeq/Apps/InvitaeBulkCsvParser -t TempPHIDelme \n"+
				"     -g invitaeGeneRegions17Oct2023.txt -p Huntsman_patients_20230815.csv -v\n"+
				"     Huntsman_variants_20230815.csv -f ~/Indexes/hg19.fasta -c ~/PMRRepo/\n"+

				"\n**************************************************************************************\n");
	}

	public File getJobsDirectory() {
		return jobsDirectory;
	}
	
	private void loadInvitaeGeneCoor(){
		String[] lines = IO.loadFile(hg38Genes);
		iGenCoorHg38 = new HashMap<String, Bed>(lines.length);
		for (String l: lines) {
			l = l.trim();
			if (l.startsWith("#") || l.length()==0) continue;
			String[] b = Misc.TAB.split(l);
			String geneNameForMatching = Misc.WHITESPACE.matcher(b[0]).replaceAll("_");
			
			//make name concatinate for bed
			StringBuilder sb = new StringBuilder();
			sb.append(Misc.WHITESPACE.matcher(b[0]).replaceAll("_")); 
			sb.append(":"); 
			sb.append(b[4]);
			
			//(String chromosome, int start, int stop, String name, double score, char strand)
			
			Bed bed = new Bed(b[1], Integer.parseInt(b[2]), Integer.parseInt(b[3]), sb.toString(), 0.0, b[5].charAt(0));
			
			
			iGenCoorHg38.put(geneNameForMatching.toUpperCase(), bed);
		}
	}

	public void buildBedHeader() {
		StringBuilder sb = new StringBuilder();
		//add in standard header
		sb.append("# filePathInvitaeVariants="+variantsFile+"\n");
		sb.append("# filePathInvitaePatients="+variantsFile+"\n");
		sb.append("# parseDate="+Misc.getDateNoSpaces()+"\n");
		sb.append("# genomeBuild=hg38\n");
		sb.append("# chr\tStart\tStop\tInvitaeGeneName:HUGOGeneName\t0\tStrand");
		bedHeader = sb.toString();
}
	
	public void buildVcfHeader() {
			StringBuilder sb = new StringBuilder();
			//add in standard header
			sb.append("##fileformat=VCFv4.2\n");
			sb.append("##filePathInvitaeVariants="+variantsFile+"\n");
			sb.append("##filePathInvitaePatients="+patientsFile+"\n");
			sb.append("##parseDate="+Misc.getDateNoSpaces()+"\n");
			//add alt and info
			sb.append(altDel);
			sb.append(altDup);
			sb.append(infoGene);
			sb.append(infoTranscript);
			sb.append(infoCDot);
			sb.append(infoPDot);
			sb.append(infoZygosity);
			sb.append(infoInterpretation);
			sb.append(infoMethod);
			sb.append(infoEnd);
			//chrom line
			sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
			vcfHeader = sb.toString();
	}
	
	public HashMap<String, Bed> getiGenCoor() {
		return iGenCoorHg38;
	}

	public IndexedFastaSequenceFile getFasta() {
		return fasta;
	}

	public String getVcfHeader() {
		return vcfHeader;
	}

	public String getBedHeader() {
		return bedHeader;
	}



}
