package edu.utah.seq.run.ambry;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.json.JSONArray;
import org.json.JSONObject;

import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class AmbryDataWrangler {

	//user fields
	private File xlsxDirectory = null;
	private File resultsDirectory = null;
	private File fastaReference = null;
	private File transcriptReference = null;
	private File jannovarApp = null;
	private File picard = null;
	private File chainFile = null;
	private File geneRegionFile = null;
	private File ambryPmrFile = null;
	private File smmRegistryDirectory = null;
	private boolean verbose = true;
	
	//internal
	private AmbryXlsxParser[] parsers = null;
	private HashSet<String> transcriptCDots = new HashSet<String>();
	private HashMap<String, String> cDotVcf = new HashMap<String,String>();
	private HashMap<String, AmbryPatient> allPatients = new HashMap<String, AmbryPatient>();
	private HashMap<String, Bed> geneNameBed = null;
	private HashSet<String> ambryTestIds = null;
	private File individualPatientResults = null;
	private int numPriorResults = 0;
	
	//ordered files for the SubjectMatchMaker and Ambry Data Wrangler
	private ArrayList<String> patientSMMInfo = new ArrayList<String>();
	private ArrayList<File> patientSaveDir = new ArrayList<File>();
	
    static final SimpleDateFormat dateFormatter = new SimpleDateFormat("dd-MMM-yyyy", Locale.ENGLISH);
    Calendar calendar = Calendar.getInstance();
	
	
	public AmbryDataWrangler(String[] args) {
		
		try {
			processArgs(args);
			
			parseBulkSpreadsheets();
			
			convertCDotAnnotations();
			
			groupPatients();
			
			assignVcfsCoordinates();
			
			findAmbryCoorNoVcf();
			
			//crossMap2Hg38();
			picardLiftoverHg38();
			
			parseGeneBedFile();
			
			loadAmbryPmrFile();
			
			writeOutTmpResults();
			
			if (smmRegistryDirectory != null && patientSaveDir.size()!=0) {
				runSubjectMatchMaker();
				// takes a very long time so skip and do in snakemake 
				// IO.deleteDirectoryViaCmdLine(individualPatientResults);
			}
			
			IO.pl("\nAll complete!");
			
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("Ambry Bulk Processing failed!");
		}
	}
	
	private void runSubjectMatchMaker() throws IOException {
		IO.pl("\nRunning the SubjectMatchMaker with patient PHI to find or create their PMR IDs...");
		//make dir to hold phi
		File smmDir = new File (resultsDirectory, "SubjectMatchMakerDelme_PHI");
		if (smmDir.exists()==false && smmDir.mkdir() == false) throw new IOException("Error: failed to make the SMM dir "+smmDir);

		//write out the info file for smm
		File queryFile = new File (smmDir, "ambryPatientQueries_PHI.txt");
		IO.writeArrayList(patientSMMInfo, queryFile);
		File queryResDir = new File (smmDir, "SMMQueryResults_PHI");
		queryResDir.mkdir();

		//run the SMM
		String[] args = {
				"-r", smmRegistryDirectory.getCanonicalPath(),
				"-q", queryFile.getCanonicalPath(),
				"-o", queryResDir.getCanonicalPath(),
				"-v", //turn off verbosity
				"-t", "4", //limit threads to 4
				"-a", //adding queries not found to registry
				"-c", //making name matching case-insensitive for EDW
				"-u"  //fill in info not found to registry entries from queries, e.g. missing mrns, otherIds
		};
		SubjectMatchMaker smm = new SubjectMatchMaker(args);
		Subject[] queries = smm.getQuerySubjects();

		if (queries.length != patientSaveDir.size()) throw new IOException("Error: the size of the SMM queries does not match the saved patients!");
		File yJobs = new File (resultsDirectory, "YJobs");
		yJobs.mkdir();
		
		for (int i=0; i< queries.length; i++) {
			String coreId = queries[i].getCoreIdNewOrMatch();
			//AA2mF6Vy/Ambry/A032049_SL419345_SL419548_SL420681/ClinicalReport/
			File patientDir = new File (yJobs, coreId+"/Ambry/"+patientSaveDir.get(i).getName()+"/ClinicalReport/");
			patientDir.mkdirs();		
			IO.copyDirectoryRecursive(patientSaveDir.get(i), patientDir, null);
		}
		
		//delete the query dir
		IO.deleteDirectoryViaCmdLine(smmDir);
		
	}

	private void parseGeneBedFile() {
		geneNameBed = new HashMap<String, Bed>();
		Bed[] regions = Bed.parseFile(geneRegionFile, 0, 0);
		for (Bed b: regions) {
			String[] geneNames = Misc.SEMI_COLON.split(b.getName());
			for (String gn: geneNames) geneNameBed.put(gn, b);
		}
		
	}

	private void writeOutTmpResults() throws Exception {
		IO.pl("\nMerge results and write out results...");
		
		//Collect dates to sort for most recent date
		
		for (AmbryPatient ap: allPatients.values()) {
			TreeMap<Date, AmbryResult> dateResult = new TreeMap<Date, AmbryResult>();
			//sort results by date earliest old to newest most recent
			for (AmbryResult ar: ap.getAmbryResults()) {
				String signDate = ar.getReportSignedDate();
				if (signDate == null) throw new IOException("Error: Failed to find a report signed date for "+ar.getPatientKey());
				Date date = dateFormatter.parse(signDate);
				dateResult.put(date, ar);
			}
			
			HashMap<String, AmbryGeneResult> hg38VcfKeyAGR = new HashMap<String, AmbryGeneResult>();
			HashSet<String> testedGenes = new HashSet<String>();
			
			//load the results starting older to most recent, the most recent will replace an older for dups
			for (AmbryResult ar: dateResult.values()) {
				if (ar.isTestedButNoGeneResults()) addGenes(testedGenes, ar);
				else {
					boolean variantsFound = false;
					for (AmbryGeneResult ag: ar.getGeneResultsWithParsibleResults()) {
						//any coordinates?
						String[] chrPosRefAlt = ag.getHg38Coordinates();
						if (chrPosRefAlt !=null) {
							hg38VcfKeyAGR.put(Misc.stringArrayToString(chrPosRefAlt, "_"), ag);
							variantsFound = true;
						}
					}
					if (variantsFound) addGenes(testedGenes,ar);
				}
			}
			
			//any results to write out?
			if (testedGenes.size()!=0) {
				Date lastDate = dateResult.lastKey();
				AmbryResult lastAR = dateResult.get(lastDate);
				writeOutVcfAndBed(ap, testedGenes, hg38VcfKeyAGR, lastDate, lastAR);
			}
		}
		
		IO.pl("\tNumber new results for annotating:\t"+patientSaveDir.size());
		IO.pl("\tNumber old results skipped:\t"+numPriorResults);
		
	}

	private void writeOutVcfAndBed(AmbryPatient ap, HashSet<String> testedGenes, HashMap<String, AmbryGeneResult> hg38VcfKeyAGR, Date date, AmbryResult lastAR) throws Exception {
		

		//Tempus 	TL-23-8WGYPWFK_XT.V4_2023-07-25_deid_Sonam_Puri_F.json
		//Invitae	RQ140796_AR128388_2017-12-4_Marjan_Champine.json
		//Caris		TN23-212942_20231005105551_deid_Umang_Swami.xml
		//Avatar	A025977_SL392098_SL392194_SL400389_IDTv1_GI_M.json

		//AccessionNumber1_AccessionNumber2_AccessionNumber3_LatestDate_PhysicianFirstName_LastName.json
		//AccessionNumber1_AccessionNumber2_AccessionNumber3_LatestDate.bed or .vcf
		//24-642550_24-644978_2024-25-12_Sarah_Colonna.json and
		//24-642550_24-644978_2024-25-12.bed or .vcf
		
		//Placing in YJobs/PmrID/Ambry/AccessionNumber1_AccessionNumber2_AccessionNumber3_LatestDate/ClinicalReport
		// YJobs/AA2mF6Vz/Ambry/24-642550_24-644978_2024-25-12/ClinicalReport/...

		//Pull accessionNumbers
		TreeSet<String> sortedANs = new TreeSet<String>();
		for (AmbryResult ar: ap.getAmbryResults()) sortedANs.add(ar.getAccessionNumber());
		String ans = Misc.treeSetToString(sortedANs, "_");

		//Add most recent date
		StringBuilder sb = new StringBuilder(ans);
		sb.append("_");
		calendar.setTime(date);
		sb.append(calendar.get(Calendar.YEAR)); sb.append("-");
		sb.append((1+calendar.get(Calendar.MONTH))); sb.append("-");
		sb.append(calendar.get(Calendar.DAY_OF_MONTH)); 
		
		String fileNameNoPhys = sb.toString();
		
		//Is this already in the PMR

		if (ambryTestIds.contains(fileNameNoPhys)) {
			//IO.pl("\tFound "+fileNameNoPhys+" in the PMR, skipping");
			numPriorResults++;
			return;
		}
		//else IO.pl("\tWriting bed and vcf files for "+fileNameNoPhys);
		
		//make tmp dir for writing files
		File tmpPatientSaveDir = new File(individualPatientResults, fileNameNoPhys);
		if (tmpPatientSaveDir.exists()==false && tmpPatientSaveDir.mkdirs() == false) throw new IOException("Error: failed to make the temp patient dir "+tmpPatientSaveDir);
		
		//save for SMM and data wrangler
		patientSMMInfo.add(lastAR.getSubjectMatchMakerInfo());
		patientSaveDir.add(tmpPatientSaveDir);
		
		//Add physicians names
		// Colonna, Sarah, MD
		String cleanedPhysician = lastAR.getPhysician().replaceAll("'", "");
		cleanedPhysician = cleanedPhysician.replaceAll(" ", "");
		String[] lastFirstNames = Misc.COMMA.split(cleanedPhysician);
		if (lastFirstNames.length < 3) throw new IOException("Error: physician name doesn't have 3 comma delimited fields: "+lastAR.getPhysician());
		sb.append("_");
		sb.append(lastFirstNames[1]);
		sb.append("_");
		sb.append(lastFirstNames[0]);

		String fileNameWithPhys = sb.toString();

		//genes to export?
		File bedFile = new File(tmpPatientSaveDir, fileNameNoPhys+ ".bed.gz");
		Bed[] toSort = new Bed[testedGenes.size()];
		int i = 0;
		for (String geneName: testedGenes) {
			Bed b = geneNameBed.get(geneName);
			if (b==null) throw new IOException("\nError: failed to find a bed record for "+geneName);
			toSort[i++] = b;
		}
		Arrays.sort(toSort);
		Bed.writeToGzippedFile(toSort, bedFile);

		//any vcfs to save?
		if (hg38VcfKeyAGR.size()!=0) {
			ArrayList<String> vcfInfo = new ArrayList<String>();
			vcfInfo.add(AmbryGeneResult.fetchVcfHeader());

			int counter = 1;
			for (AmbryGeneResult agr: hg38VcfKeyAGR.values()) {
				String vcf = agr.getVcfLine(new Integer (counter++).toString());
				vcfInfo.add(vcf);
			}
			File vcfFile = new File(tmpPatientSaveDir, fileNameNoPhys+ ".vcf.gz");
			// this is unsorted!
			IO.writeArrayListGzipped(vcfInfo, vcfFile);
		}
		
		//write json
		File jsonFile = new File(tmpPatientSaveDir, fileNameWithPhys+ ".json");
		writeJson(ap, jsonFile);
	}



	private void writeJson(AmbryPatient ap, File fileName) {
		JSONArray ja = new JSONArray();
		for (AmbryResult r: ap.getAmbryResults()) {
			JSONObject report = new JSONObject();
			ja.put(report);
			add ("specimenCollectionDate", r.getSpecimenCollectionDate(), report);
			add ("familialVariantPreviouslyIdentified", r.getFamilialVariantPreviouslyIdentified(), report);
			add ("orderType", r.getOrderType(), report);
			add ("accessionNumber", r.getAccessionNumber(), report);
			add ("orderNumber", r.getOrderNumber(), report);
			add ("reportSignedDate", r.getReportSignedDate(), report);
			add ("reportComment", r.getReportComment(), report);
			add ("genesTested", r.getGenesTested(), report);
			add ("geneticTestPanelOrdered", r.getGeneticTestPanelOrdered(), report);
			add ("organizationName", r.getOrganizationName(), report);
			add ("physician", r.getPhysician(), report);
			add ("additionalRecipient", r.getAdditionalRecipient(), report);
		}
		//if (verbose) IO.pl(ja.toString(5));
		IO.writeString(ja.toString(5), fileName);
	}

	private void add(String name, String value, JSONObject jo) {
		if (value!=null && value.length()!=0) {
			value = value.replaceAll("\"", "'");
			jo.put(name, value);
		}
	}

	private void addGenes(HashSet<String> testedGenes, AmbryResult ar) {
		// CHEK2: c.428A>G (p.H143R) or AIP,ALK,APC,ATM,ATRIP,AXIN2,BAP1,BARD1,BMPR1A,BRCA1,BRCA2,BRIP1,CDC73,CDH1,CDK4,CDKN1B,CDKN2A,CEBPA,CFTR,CHEK2,CPA1,CTNNA1,CTRC,DDX41,DICER1,EGFR,EGLN1,EPCAM,ETV6,FH,FLCN,GATA2,GREM1,HOXB13,KIF1B,KIT,LZTR1,MAX,MBD4,MEN1,MET,MITF,MLH1,MLH3,MSH2,MSH3,MSH6,MUTYH,NF1,NF2,NTHL1,PALB2,PALLD,PDGFRA,PHOX2B,PMS2,POLD1,POLE,POT1,PRKAR1A,PRSS1,PTCH1,PTEN,RAD51B,RAD51C,RAD51D,RB1,RET,RNF43,RPS20,RUNX1,SDHA,SDHAF2,SDHB,SDHC,SDHD,SMAD4,SMARCA4,SMARCB1,SMARCE1,SPINK1,STK11,SUFU,TERT,TMEM127,TP53,TSC1,TSC2,VHL,WT1
		String gt = ar.getGenesTested();
		if (gt == null) return;
		if (gt.contains(":")) {
			String[] f = Misc.COLON.split(gt);
			testedGenes.add(f[0].trim());
		}
		else {
			String[] genes = Misc.COMMA.split(gt);
			for (String g: genes) testedGenes.add(g.trim());
		}
		
	}
	
	private void picardLiftoverHg38() throws IOException {
		IO.pl("\nPicard liftover Ambry genomic coordinates to Hg38...");

		/*
		jj /Users/u0028003/HCI/ClinicalGenomics/Ambry/RefFilesForAmbryBulkParser/picard.jar 
		LiftoverVcf 
		-I /Users/u0028003/HCI/ClinicalGenomics/Ambry/27Aug2025/BulkParsingResults/vcfToCrossMapHg19_TmpDelme.vcf  
		-O /Users/u0028003/HCI/ClinicalGenomics/Ambry/27Aug2025/BulkParsingResults/vcfToCrossMapHg38_LiftoverPass.vcf  
		--REJECT /Users/u0028003/HCI/ClinicalGenomics/Ambry/27Aug2025/BulkParsingResults/vcfToCrossMapHg38_LiftoverFail.vcf  
		-C /Users/u0028003/HCI/ClinicalGenomics/Ambry/RefFilesForAmbryBulkParser/hg19ToHg38.over.chain.gz  
		-R /Users/u0028003/HCI/ClinicalGenomics/Ambry/RefFilesForAmbryBulkParser/hs38DH.fa
		*/
		
		//create vcf from vcf coordinates from Ambry
		File vcfHg19 = createHg19Vcf();
		File vcfHg38 = new File (resultsDirectory, "vcfCMappedHg38_TmpDelme.vcf");
		File vcfHg38Reject = new File (resultsDirectory, "vcfMappedHg38_TmpDelme.reject.vcf");
		
		String[] cmd = {"java","-jar", picard.getCanonicalPath(), 
				"LiftoverVcf", 
				"-C", chainFile.getCanonicalPath(),
				"-I", vcfHg19.getCanonicalPath(), 
				"-R", fastaReference.getCanonicalPath(), 
				"-O", vcfHg38.getCanonicalPath(),
				"--REJECT", vcfHg38Reject.getCanonicalPath(),
				"--USE_JDK_DEFLATER", "--USE_JDK_INFLATER"
		};
		
		String[] runLog = IO.executeViaProcessBuilder(cmd, false);
		
		//check log
		for (String s: runLog) {
			if (s.contains("Exception")) {
				for (String x: runLog) IO.el("\t"+x);
				throw new IOException ("Error executing: "+Misc.stringArrayToString(cmd, " "));
			}
		}
		
		loadMappedResults(vcfHg38);
		
		//clean up
		vcfHg38.delete();
		vcfHg19.delete();
		vcfHg38Reject.delete();
	}

	private void crossMap2Hg38() throws IOException {
		IO.pl("\nCrossMapping Ambry genomic coordinates to Hg38...");
		//CrossMap vcf {chainFile} {input} {hg38Index} {output}
		
		//create vcf from vcf coordinates from Ambry
		File vcfHg19 = createHg19Vcf();
		
		
		
		
		File vcfHg38 = new File (resultsDirectory, "vcfCrossMappedHg38_TmpDelme.vcf");
		File vcfHg38UnMap = new File (resultsDirectory, "vcfCrossMappedHg38_TmpDelme.vcf.unmap");
		
		String[] cmd = { "crossMap.getCanonicalPath()", "vcf", 
				chainFile.getCanonicalPath(), vcfHg19.getCanonicalPath(), fastaReference.getCanonicalPath(), 
				vcfHg38.getCanonicalPath()
		};
		
		String[] runLog = IO.executeViaProcessBuilder(cmd, false);
		
		//check log
		for (String s: runLog) {
			if (s.contains("Traceback") || s.contains("Error")) {
				for (String x: runLog) IO.el("\t"+x);
				throw new IOException ("Error executing: "+Misc.stringArrayToString(cmd, " "));
			}
		}
		
		Misc.printErrAndExit("hg19 vcf written to "+vcfHg19+ "\nhg38 "+vcfHg38);
		
		loadMappedResults(vcfHg38);
		
		//clean up
		vcfHg38.delete();
		vcfHg19.delete();
		vcfHg38UnMap.delete();
	}

	private void loadMappedResults(File vcfHg38) throws IOException {
		//load hash with conversions
		String[] lines = IO.loadFile(vcfHg38);
		HashMap<String,String[]> vcfKeyHg38Coor = new HashMap<String,String[]>();
		for (String l: lines) {
			if (l.startsWith("chr")) {
				//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
				String[] fields = Misc.TAB.split(l);
				String[] coor = new String[4];
				coor[0] = fields[0];
				coor[1] = fields[1];
				coor[2] = fields[3];
				coor[3] = fields[4];
				vcfKeyHg38Coor.put(fields[2], coor);
			}
		}
		
		//add to results
		for (AmbryPatient ap: allPatients.values()) {
			for (AmbryResult ar: ap.getAmbryResults()) {
				for (AmbryGeneResult ag: ar.getGeneResultsWithParsibleResults()) {
					String vcfHg19Key = ag.getVcfHg19Key();
					if (vcfHg19Key !=null) {
						String[] hg38Coor = vcfKeyHg38Coor.get(vcfHg19Key);
						//OK if this is null
						ag.setVcfHg38ChromPosRefAlt(hg38Coor);
						//compare to Jannovar
						if (hg38Coor !=null && ag.getCDotHg38Vcf() !=null) {
							String cDot = Misc.stringArrayToString(ag.getCDotHg38Vcf(),"_");
							String liftover = Misc.stringArrayToString(hg38Coor, "_");
							if (cDot.equals(liftover)== false) throw new IOException("Error: CDot Hg38 "+cDot+" does not agree with remapped Hg38 "+liftover);
						}
					}
				}
			}
		}
	}

	private File createHg19Vcf() {

		TreeSet<String> hg19VcfKeys = new TreeSet<String>();
		for (AmbryPatient ap: allPatients.values()) {
			for (AmbryResult ar: ap.getAmbryResults()) {
				for (AmbryGeneResult ag: ar.getGeneResultsWithParsibleResults()) {
					String vcfHg19Key = ag.getVcfHg19Key();
					if (vcfHg19Key !=null) hg19VcfKeys.add(vcfHg19Key);
				}
			}
		}
		if (verbose) IO.pl("\tVcfs to Liftover: "+hg19VcfKeys.size());
		if (hg19VcfKeys.size()==0) return null;
		
		StringBuilder sb = new StringBuilder();
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		for (String k: hg19VcfKeys) {
			String[] f = Misc.UNDERSCORE.split(k);
			Double chrNum = Double.parseDouble(f[0]);
			Double pos = Double.parseDouble(f[1]);
			sb.append("chr");
			sb.append(chrNum.intValue());
			sb.append("\t");
			sb.append(pos.intValue());
			sb.append("\t");
			sb.append(k);
			sb.append("\t");
			sb.append(f[2]);
			sb.append("\t");
			sb.append(f[3]);
			sb.append("\t.\t.\t.\n");
		}
		File vcf = new File (resultsDirectory, "vcfToReMapHg19_TmpDelme.vcf");
		IO.writeString(sb.toString(), vcf);
		return vcf;
	}

	private void findAmbryCoorNoVcf() {  
		IO.pl("\nFinding results to export...");

		IO.pl("Num patients "+allPatients.size());
		int numGeneResultsWithCoord = 0;
		int numResAmbCoorNoMapped = 0;
		for (AmbryPatient ap: allPatients.values()) {
			//IO.pl("\nPatient:\t"+ap.getPatientKey());
			for (AmbryResult ar: ap.getAmbryResults()) {
				//IO.pl("\tTest:\t"+ ar.getAccessionNumber()+"\t"+ar.getGeneticTestPanelOrdered());
				for (AmbryGeneResult ag: ar.getGeneResultsWithParsibleResults()) {
					//IO.pl("\t\tResult:\t"+ag.getTranscriptCDot()+"\tAmbryChrom:\t"+ag.getChromosome()+"\tMapped:\t"+ag.getCDotHg38Vcf());
					if (ag.getChromosome()!=null || ag.getCDotHg38Vcf()!=null) numGeneResultsWithCoord++;
					if (ag.getChromosome()!=null && ag.getCDotHg38Vcf()==null) numResAmbCoorNoMapped++;
				}
			}
		}
		IO.pl("Num results w/ genomic coor: "+numGeneResultsWithCoord);
		IO.pl("Num coor results w/ no mapped: "+numResAmbCoorNoMapped);
	}

	private void assignVcfsCoordinates() {
		IO.pl("\nAssigning vcf coordinates to cDot hgvs...");

		for (AmbryPatient ap: allPatients.values()) {
			for (AmbryResult ar: ap.getAmbryResults()) {
				for (AmbryGeneResult ag: ar.getGeneResultsWithParsibleResults()) {
					String geneResultTranscriptCDot = ag.getTranscriptCDot();
					if (geneResultTranscriptCDot!=null && cDotVcf.containsKey(geneResultTranscriptCDot)) {
						//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
						String[] fields = Misc.TAB.split(cDotVcf.get(geneResultTranscriptCDot));
						String[] coor = new String[4];
						coor[0] = fields[0];
						coor[1] = fields[1];
						coor[2] = fields[3];
						coor[3] = fields[4];
						ag.setCDotHg38Vcf(coor);
					}
				}
			}
		}
	}
	
	private void groupPatients() {
		IO.pl("\nGrouping patients...");
				
		for (int i=0; i< parsers.length; i++) {
			HashMap<String, AmbryPatient> patients = parsers[i].getPatientKeyResults();
			IO.pl("\tNumber patients "+patients.size());
			for (String key: patients.keySet()) {
				AmbryPatient am =patients.get(key);
				ArrayList<AmbryResult> ar = am.getAmbryResults();
				//already exists so just add in the results
				if (allPatients.containsKey(key)) {
					ArrayList<AmbryResult> allArs = allPatients.get(key).getAmbryResults();
					allArs.addAll(ar);
				}
				//new 
				else allPatients.put(key, am);
			}
		}
		IO.pl("\t\tTotal patients "+allPatients.size());
		
	}

	private void convertCDotAnnotations() throws IOException {
		IO.pl("\nConverting "+transcriptCDots.size()+" cDot annotations to genomic coordinates...");
		//write out the toConvert.txt file
		File toCon = new File (resultsDirectory, "cDotsToConvert_TmpDelme.txt");
		IO.writeHashSet(transcriptCDots, toCon);
		
		//Build and execute the jannovar command
		File vcf = new File (resultsDirectory, "convertedCDots_TmpDelme.vcf");
		
		String[] cmd = {"java", "-jar", "-Xmx2G", jannovarApp.getCanonicalPath(),
				"hgvs-to-vcf", "-r", fastaReference.getCanonicalPath(), "-d",
				transcriptReference.getCanonicalPath(), "-i",
				toCon.getCanonicalPath(),
				"-o", vcf.getCanonicalPath()};

		String[] runLog = IO.executeViaProcessBuilder(cmd, false);
		//check log
		for (String s: runLog) {
			if (s.contains("Exception")) {
				for (String x: runLog) IO.el("\t"+x);
				throw new IOException ("Error executing: "+Misc.stringArrayToString(cmd, " "));
			}
		}
		if (verbose) for (String s: runLog) IO.pl("\t"+s);
		
		loadJannovarResults(toCon, vcf);
		
		//clean up
		toCon.delete();
		vcf.delete();
	}

	private void loadJannovarResults(File toCon, File vcf) throws IOException {
		//load results
		String[] vcfLines = IO.loadFile(vcf);
		String[] input = IO.loadFile(toCon);
		
		int chromIndex = -1;
		for (int i=0; i< vcfLines.length; i++) {
			if (vcfLines[i].startsWith("#CHROM")) {
				chromIndex = i+1;
			}
		}
		if (chromIndex == -1) throw new IOException ("\nERROR: #CHROM line not found in Jannovar output file "+vcf);
		
		for (int i=0; i<input.length; i++) {
			int index = i+ chromIndex;
			if (vcfLines[index].contains("PARSE_ERROR")== false) {
				if (verbose) IO.pl(input[i]+ " -> "+vcfLines[index]);
				cDotVcf.put(input[i], vcfLines[index]);
			}
		}
		IO.pl("\tNumber converted "+cDotVcf.size());
	}

	private void parseBulkSpreadsheets() throws Exception {
		IO.pl("Parsing Ambry bulk spreadsheets...");
		File[] xlsx = IO.extractFiles(xlsxDirectory, ".xlsx");
		
		ArrayList<AmbryXlsxParser> parsersAL = new ArrayList<AmbryXlsxParser>();
		for (int i=0; i< xlsx.length; i++) {
			if (xlsx[i].getName().startsWith("~")) continue;
			AmbryXlsxParser p = new AmbryXlsxParser(xlsx[i], verbose);
			transcriptCDots.addAll(p.getTranscriptCDots());
			parsersAL.add(p);
		}
		
		parsers = new AmbryXlsxParser[parsersAL.size()];
		parsersAL.toArray(parsers);
	}
	
	public void loadAmbryPmrFile() throws IOException {
		ambryTestIds = new HashSet<String>();
		Pattern p = Pattern.compile(".+/Ambry/(.+)/ClinicalReport/.+json");
		BufferedReader in = IO.fetchBufferedReader(ambryPmrFile);
		String line = null;
		Matcher m = null;
		while ((line = in.readLine())!=null) {
			m = p.matcher(line);
			if (m.matches()) {
				ambryTestIds.add(m.group(1));
			}
		}
		in.close();
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AmbryDataWrangler(args);
	}	
	

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " "));
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'x': xlsxDirectory = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'j': jannovarApp = new File(args[++i]); break;
					case 'f': fastaReference = new File(args[++i]); break;
					case 't': transcriptReference = new File(args[++i]); break;
					case 'p': picard = new File(args[++i]); break;
					case 'c': chainFile = new File(args[++i]); break;
					case 'g': geneRegionFile = new File(args[++i]); break;
					case 'a': ambryPmrFile = new File(args[++i]); break;
					case 's': smmRegistryDirectory = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (xlsxDirectory == null || xlsxDirectory.isDirectory()==false) Misc.printErrAndExit("\nERROR: failed to find your Xlsx directory containing Ambry bulk exports.\n");
		if (resultsDirectory == null) Misc.printErrAndExit("\nERROR: failed to find your Xlsx directory containing Ambry bulk exports.\n");
		resultsDirectory.mkdirs();
		if (resultsDirectory.exists()==false || resultsDirectory.isDirectory()==false) Misc.printErrAndExit("\nERROR: failed to find or make your results directory.\n");
		individualPatientResults = new File(resultsDirectory,"IndividualPatientResults");
		individualPatientResults.mkdir();
		if (jannovarApp == null || jannovarApp.exists()==false) Misc.printErrAndExit("\nERROR: please provide a path to your jannovar-cli-xxx.jar file.\n");
		if (fastaReference == null || fastaReference.exists()==false) Misc.printErrAndExit("\nERROR: please provide a path to your fasta reference file with fai index and dict files in the parent directory.\n");
		if (transcriptReference == null || transcriptReference.exists()==false) Misc.printErrAndExit("\nERROR: please provide a path to your jannovar serialized transcript reference file.\n");
		if (picard == null || picard.exists()==false) Misc.printErrAndExit("\nERROR: please provide a path to your picard.jar executable.\n");
		if (chainFile == null || chainFile.exists()==false) Misc.printErrAndExit("\nERROR: please provide a path to your hg19ToHg38.over.chain.gz file.\n");
		if (geneRegionFile == null || geneRegionFile.exists()==false) Misc.printErrAndExit("\nERROR: please provide a path to your bed file containing gene coordinates to use in creating an interrogated gene file for each patient.\n");
		if (ambryPmrFile == null || ambryPmrFile.exists()==false) Misc.printErrAndExit("\nERROR: please provide a path to file containing Ambry files currently present in the PMR.\n");
		

	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                        Ambry Data Wrangler : October 2025                        **\n" +
				"**************************************************************************************\n" +
				"App for parsing xlsx bulk results spreadsheets to vcf, bed, and json results. Uses\n"+
				"Jannovar to convert Ambry cDots to vcf coordinates. Uses Liftover to convert the hg19\n"+
				"Ambry chrom pos ref alt to hg38. Patient results with cDots that cannot be converted\n"+
				"and don't contain Ambry vcf info are skipped. Patient reports with no results are\n"+
				"saved. Use the patientInfoForSMM_PHI.txt file to map patient PHI to PMR IDs. Use the  \n"+
				"patientFileNames.txt to link the PMR IDs to the files written to the -r Dir.\n"+

				"\nOptions:\n"+
				"-x Directory containing xlsx files from Ambry\n" +
				"-r Directory to write out the results.\n"+
				"-j Path to your jannovar-cli jar file, see https://github.com/charite/jannovar\n"+
				"-t Path to your jannovar serialized transcript file.\n"+
				"-f Path to your matching fasta reference file with fai index and dict files in the\n"+
				"      parent directory.\n"+
				"-p Path to the picard.jar executable.\n"+
				"-c Path to the hg19ToHg38.over.chain.gz liftover chain file.\n"+
				"-g Path to the gene region bed file for all interrogated Ambry genes.\n"+
				"-a Path to a file containing the 'Ambry' files in the PMR, e.g. grep Ambry \n"+
				"      hcibioinfo-patient-molecular-repo.index.txt > ambryPmr.txt\n"+
				"-s (Optional) Path to a directory containing the currentRegistry_xxx file for the PMR\n"+
				"      and SubjectMatchMaker. This will move the parsed dirs into a YJobs/PmrID/...\n"+
				"      structure\n"+
				
                "\nExample: java -jar -Xmx2G ~/USeqApps/AmbryDataWrangler -x AmbryExpoFiles -r Results\n"+
                "   -j ~/Apps/jannovar-cli-0.36.jar -f ~/Refs/hs38DH.fa -a ~/Refs/ambryPmr.txt -t \n"+
                "   ~/Refs/refseq_curated_109_hg38.ser -p ~/Apps/picard.jar -c \n"+
                "   ~/Refs/hg19ToHg38.over.chain.gz -g ~/Refs/hg38AmbryGeneRegionsAug2025.bed \n"+
                "   -s ~/PmrRegisry/"+
                "\n"+
				"**************************************************************************************\n");
	}
}
