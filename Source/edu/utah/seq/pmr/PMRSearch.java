package edu.utah.seq.pmr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.config.Configurator;

import edu.utah.hci.misc.Util;
import edu.utah.seq.vcf.json.TempusJson2Vcf;
import edu.utah.seq.vcf.json.TempusJsonSummary;
import edu.utah.seq.vcf.json.TempusSpecimen;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Json2Vcf;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonSummary;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Specimen;
import edu.utah.seq.vcf.xml.caris.CarisXmlVcfParser;
import util.gen.IO;
import util.gen.Misc;

/**Loads the HCI Patient Molecular Repository*/
public class PMRSearch {
	
	//fields
	private File clinicalReportDir = null;
	private long maximumDaysOld = 0;
	private String awsPatientDirUri = "s3://hcibioinfo-patient-molecular-repo/Patients/"; //must end with /
	private String awsPath = "aws";
	private String profile = "default";
	private HashMap<String, String> envPropToAdd = new HashMap <String, String>();
	private boolean verbose = false;
	private int numberThreads = 10;
	private File awsRepoList = null;
	private HashMap<String, Patient> idPatient = new HashMap<String, Patient>();
	private HashMap<String, Dataset> idDataset = new HashMap<String, Dataset>();
	private TreeMap<String, HashSet<String>> sourceDatasets = null; //Avatar,Tempus,Caris,etc : TestIds
	private TreeMap<String, TreeSet<String>> sourceCombinations = null; //Patients with combinations of tests
	private ArrayList<String> cmdsForParallel = new ArrayList<String>();
	private HashMap<String, String> icd10CodeDesc = null;
	private HashMap<String, String> icdOCodeMorphology = null;
	private HashMap<String, String> icdOCodeTopology = null;
	private String inputNonInteractive = null;
	
	//names
	private String clinicalReportDirName = "ClinicalReport/";
	private String avatarSourceName = "Avatar";
	private String carisSourceName = "Caris";
	private String tempusSourceName = "Tempus";

	//for searching
	private TreeMap<String, TreeSet<String>> diagnosis_DatasetName = null;
	private String diagnosisSearchString = null;
	private TreeMap<String, TreeSet<String>> disGrpPhy_DatasetName = null;
	private String disGrpPhySearchString = null;
	private TreeMap<String, TreeSet<String>> specimenSite_DatasetName = null;
	private String specimenSiteSearchString = null;
	private TreeMap<String, TreeSet<String>> specimenId_DatasetName = null;
	private String specimenIdSearchString = null;
	private TreeMap<String, TreeSet<String>> pmrId_DatasetName = null;
	private String datasetIdSearchString = null;
	private TreeMap<String, TreeSet<String>> sex_DatasetName = null;
	private String sexSearchString = null;
	private TreeSet<String> datasetKeysFound = null;
	
	//misc
    private DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
    private LocalDate todaysDate = LocalDate.now();

	public PMRSearch (String[] args) {

		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);

			downloadAwsFileList();

			loadPatientFiles();
			
			fetchReports();

			//Avatar
			loadAvatarClinicalJsons();
			summarizeAvatarClinicalInfo();

			//Tempus
			loadTempusClinicalJsons();

			//Caris
			loadCarisClinicalXmls();  
			summarizeCarisClinicalInfo();
			
			//interactive?
			if (inputNonInteractive != null) {
				IO.pl("NonInt "+inputNonInteractive);
				String[] toRun = Misc.SEMI_COLON.split(inputNonInteractive);
				for (String i: toRun) {
					IO.pl("\tRunning "+i);
					runSearch(i);
				}
			}
			else runInteractiveSearch();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
		} catch (Exception e) {
			IO.el("\nERROR running the SearchPMR...");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	private void fetchReports() throws Exception {
		IO.pl("Looking for new reports to download...");
		fetchAvatarClinicalJsons();
		fetchTempusClinicalJsons();
		fetchCarisClinicalXmls();
		//anything to download?
		if (cmdsForParallel.size()>0) {
			IO.pl("\tDownloading "+ cmdsForParallel.size() +" new reports...");
			IO.executeViaParallel(cmdsForParallel, numberThreads);
		}
	}

	private void runInteractiveSearch() throws IOException {
		Scanner s = new Scanner(System.in);
		while (true) {  
			IO.p("\nEnter a query, hit return for the help menu, or type 'exit'");
			if (datasetKeysFound != null && datasetKeysFound.size() > 0) {
				IO.pl(". ("+datasetKeysFound.size()+" datasets selected. Use -x to clear.)");
			}
			else IO.pl();
			Util.p("> ");
			String input = s.nextLine().trim();
			String lcInput = input.toLowerCase();

			//empty? print menu
			if (input.length() == 0) UserSearch.printInteractiveMenu();
			else if (lcInput.contains("exit")) {
				s.close();
				return;
			}
			else runSearch(input);
		}
	}
	
	private void runSearch (String input) throws IOException {
		UserSearch us = new UserSearch(input);
		
		//clear prior?
		if (us.isClearPriorResults()) {
			IO.pl("Prior selected datasets cleared.");
			datasetKeysFound = null;
		}
		
		//print a data table?
		if (us.getPrintDataTable() != null) {
			printDataTable(us.getPrintDataTable());
		}
		//search
		else if (us.isGoodToSearch()) {
			if (datasetKeysFound != null) IO.pl("Selecting from prior search of "+datasetKeysFound.size()+" datasets...");
			TreeSet<String> foundSearchKeys = searchForDatasets(us);
			if (datasetKeysFound == null) datasetKeysFound = foundSearchKeys;
			else if (us.isAddMatches()) datasetKeysFound.addAll(foundSearchKeys);
			else datasetKeysFound.retainAll(foundSearchKeys);
			IO.pl("\n"+datasetKeysFound.size()+" datasets selected");
		}
		
		//print dataset files?
		if (us.isPrintMatchedDatasetURIs()) printURIs();
		
		//print dataset info?
		if (us.isPrintDatasetInfo()) printDatasetInfo();
		
		//print dataset info?
		if (us.isPrintDatasetNames()) printDatasetNames();

	}
	
	private void printDatasetInfo() throws IOException {
		IO.pl("\nDataset Info, also download and examine the ClinicalReport/xxx.json/.xml files:");
		for (String datasetId: datasetKeysFound) {
			Dataset d = idDataset.get(datasetId);
			IO.p(d.toString(datasetId));
		}
	}
	
	private void printDatasetNames() {
		IO.pl("\nDataset Names:");
		for (String datasetId: datasetKeysFound) IO.pl(datasetId);
	}
	
	private void printURIs() {
		IO.pl("\nFile URIs for each matched dataset:");
		for (String datasetId: datasetKeysFound) {
			IO.pl(datasetId);
			Dataset d = idDataset.get(datasetId);
			for (String p: d.getPartialPaths()) IO.pl("\t"+awsPatientDirUri+datasetId+"/"+p);
		}
	}

	private void printDataTable(String whatCategory) {
		TreeMap<String, TreeSet<String>> catKeys_DatasetName = null;
		
		if (whatCategory.equals("Diagnosis")) {
			if (diagnosis_DatasetName == null) makeDiagnosisMap();
			catKeys_DatasetName = diagnosis_DatasetName;
		}
		else if (whatCategory.equals("PhysicianDiseaseGroups")) {
			if (disGrpPhy_DatasetName == null) makeDiseaseGroupPhysician();
			catKeys_DatasetName = disGrpPhy_DatasetName;
		}
		else if (whatCategory.equals("SpecimenSites")) {
			if (specimenSite_DatasetName == null) makeSpecimenSiteMap();
			catKeys_DatasetName = specimenSite_DatasetName;
		}
		else if (whatCategory.equals("SpecimenIds")) {
			if (specimenId_DatasetName == null) makeSpecimenIdMap();
			catKeys_DatasetName = specimenId_DatasetName;
		}
		else if (whatCategory.equals("PmrIds")) {
			if (pmrId_DatasetName == null) makePmrIdsMap();
			catKeys_DatasetName = pmrId_DatasetName;
		}
		else if (whatCategory.equals("Sex")) {
			if (sex_DatasetName == null) makeSexMap();
			catKeys_DatasetName = sex_DatasetName;
		}
		else {
			System.out.flush();
			IO.el("\nDid not recognized '"+whatCategory+"' to print? Choose from Diagnosis, PhysicianDiseaseGroups, SpecimenSites, SpecimenIds, Sex, or PmrIds");
			System.err.flush();
		}
		printTreeMapTreeSet(catKeys_DatasetName);
	}

	/**Searches for datasets that match particular strings
	 * Categories: Diagnosis, PhysicianDiseaseGroups, SpecimenSites, SpecimenIds, Sex, PmrIds
	 * One or more search terms
	 * CaseInsensitive
	 * RequireAll to match
	 * PartialSearchTerm match*/
	private TreeSet<String> searchForDatasets(UserSearch us) {
		String keySearchString = null;
		TreeMap<String, TreeSet<String>> catKeys_DatasetName = null;
		String whatCategory = us.getDataTableToSearch();
		
		if (whatCategory.equals("Diagnosis")) {
			if (diagnosis_DatasetName == null) makeDiagnosisMap();
			catKeys_DatasetName = diagnosis_DatasetName;
			keySearchString = diagnosisSearchString;
		}
		else if (whatCategory.equals("PhysicianDiseaseGroups")) {
			if (disGrpPhy_DatasetName == null) makeDiseaseGroupPhysician();
			catKeys_DatasetName = disGrpPhy_DatasetName;
			keySearchString = disGrpPhySearchString;
		}
		else if (whatCategory.equals("SpecimenSites")) {
			if (specimenSite_DatasetName == null) makeSpecimenSiteMap();
			catKeys_DatasetName = specimenSite_DatasetName;
			keySearchString = specimenSiteSearchString;
		}
		else if (whatCategory.equals("SpecimenIds")) {
			if (specimenId_DatasetName == null) makeSpecimenIdMap();
			catKeys_DatasetName = specimenId_DatasetName;
			keySearchString = specimenIdSearchString;
		}
		else if (whatCategory.equals("PmrIds")) {
			if (pmrId_DatasetName == null) makePmrIdsMap();
			catKeys_DatasetName = pmrId_DatasetName;
			keySearchString = datasetIdSearchString;
		}
		else if (whatCategory.equals("Sex")) {
			if (sex_DatasetName == null) makeSexMap();
			catKeys_DatasetName = sex_DatasetName;
			keySearchString = sexSearchString;
		}
		else {
			System.out.flush();
			IO.el("\nDid not recognized '"+whatCategory+"'? Choose from Diagnosis, PhysicianDiseaseGroups, SpecimenSites, SpecimenIds, Sex, or PmrIds");
			System.err.flush();
			return null;
		}
		IO.pl("\nSearching data table: '"+whatCategory+"' for search term(s): '"+Misc.stringArrayToString(us.getSearchTerms(), "','")+"' caseInsensitve? "+us.isCaseInsensitive()+
				", requireAll? "+us.isAllTermsMustMatch()+", partialMatch? "+(us.isAllTermsMustMatch()==false));

		TreeSet<String> datasets = null;
		for (String toSearch: us.getSearchTerms() ) {
			Pattern pat = fetchSearchPattern(toSearch, us.isCaseInsensitive(), us.isPartialMatchesAllowed());
			TreeSet<String> datasetKeyMatch = searchKeysString(keySearchString, pat);
			if (datasets == null) datasets = datasetKeyMatch;
			else if (us.isAllTermsMustMatch()) datasets.retainAll(datasetKeyMatch);
			else datasets.addAll(datasetKeyMatch);
		}

		//IO.pl("Hits to '"+whatCategory+"' with: '"+Misc.stringArrayToString(us.getSearchTerms(), "','")+"'");
		TreeSet<String> finalPmrIds = new TreeSet<String>();
		for (String key: datasets) {
			if (us.isVerbose()) IO.pl("\t"+key+" : "+catKeys_DatasetName.get(key));
			else IO.pl("\t"+key+" : "+catKeys_DatasetName.get(key).size());
			finalPmrIds.addAll(catKeys_DatasetName.get(key));
		}
		//more than one search term?
		if (us.getSearchTerms().length > 1) {
			if (us.isVerbose()) IO.pl("Total matches for '"+whatCategory+"' "+finalPmrIds.size()+" "+finalPmrIds);
			else IO.pl("Total matches to '"+whatCategory+"' "+finalPmrIds.size());
		}
		return finalPmrIds;
	}


	/*Diagnosis search fields Avatar, Caris, Tempus
		A='Disease Type', e.g. 'BRE - Breast Cancer'
		A='Histology/Behavior', e.g. '98753 Chronic myelogenous leukemia, BCR/ABL positive'
		C=diagnosis, e.g. 'Carcinoma with neuroendocrine differentiation'
		C=pathologicDiagnosis, e.g. '"G-E junction", biopsy: Ulcerated invasive poorly differentiated large cell carcinoma (favor adenosquamous carcinoma) involving the lamina propria and submucosa.'
		T=diagnosis, e.g. 'Adenocarcinoma with mucinous differentiation'
	 */
	private void makeDiagnosisMap() {
		IO.pl("\nMaking Diagnosis summary map...");
		diagnosis_DatasetName = new TreeMap<String, TreeSet<String>>();
		String[] avaToSearch = {"Disease Type", "Histology/Behavior"};
		String[] carisToSearch = {"diagnosis", "pathologicDiagnosis"};
		
		TreeSet<String> icdTxt = new TreeSet<String>();

		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				String id = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId();
				//Avatar?
				if (d.getAvatarClinicalInfo()!=null) {
					AvatarClinicalInfo aci = d.getAvatarClinicalInfo();
					for (String keyToPull: avaToSearch) {
						//Any Tumor DNA?
						if (aci.getTumorDna() != null) {
							addDataset(aci.getTumorDna().get(keyToPull), id, diagnosis_DatasetName);
						}
						//Any Tumor RNA?
						if (aci.getTumorRna()!= null) {
							addDataset(aci.getTumorRna().get(keyToPull), id, diagnosis_DatasetName);
						}
					}
				}
				//Caris?
				else if (d.getCarisClinicalInfo()!=null) {
					LinkedHashMap<String, String> cci = d.getCarisClinicalInfo();
					for (String keyToPull: carisToSearch) {
						addDataset(cci.get(keyToPull), id, diagnosis_DatasetName);
					}
				}
				//Tempus Pre V3?
				else if (d.getTempusJsonReportInfoPreV3()!=null) {
					TempusJsonSummary tjs = d.getTempusJsonReportInfoPreV3();
					addDataset(tjs.getTempusPatient().getDiagnosis(), id, diagnosis_DatasetName);
				}
				//Tempus V3?
				else if (d.getTempusJsonReportInfoV3()!=null) {
					icdTxt.clear();
					TempusV3JsonSummary tjs = d.getTempusJsonReportInfoV3();
					for (TempusV3Specimen s : tjs.getTempusSpecimens()) {
						//path diagnosis
						if (s.getOriginPathLabDiagnosis() != null) icdTxt.add(s.getOriginPathLabDiagnosis());
						//icd 10 diagnosis
						if (s.getTempusIcd10Code() != null && icd10CodeDesc !=null) {
							for (String code: Misc.COMMA.split(s.getTempusIcd10Code())) {
								String decoded = icd10CodeDesc.get(code);
								if (decoded !=null) {
									icdTxt.add(decoded);
									s.getTempusIcd10Txt().add(decoded);
								}
								else IO.el("\tFailed to find ICD 10 Code '"+code+"' Add it!");
							}
						}	
						//icd 0 morphology
						if (s.getTempusIcdOCodeMorphology() != null && icdOCodeMorphology !=null) {
							for (String code: Misc.COMMA.split(s.getTempusIcdOCodeMorphology())) {
								String decoded = icdOCodeMorphology.get(code);
								if (decoded !=null) {
									icdTxt.add(decoded);
									s.getTempusIcdOTxtMorphology().add(decoded);
								}
								else IO.el("\tFailed to find ICD Morphology Code '"+code+"' Add it!");
							}
						}	
					}
					//add concatinate
					if (icdTxt.size() != 0) addDataset(Misc.treeSetToString(icdTxt, "; "), id, diagnosis_DatasetName);
				}
			}
		}
		diagnosisSearchString = fetchKeysSearchString (diagnosis_DatasetName.keySet());
	}

	/*PmrIds 
		All=coreId, e.g. HDD4xq3sTP 
		A=ORIENAvatarKey, e.g. A038806
		A=hciPatientId, e.g. 1102924
		A=avatarDatasetId, e.g. A033031_SL422942_SL426281_SL423826
		C=labReportID,	e.g. TN20-149935
		T=tempusReportID, e.g. TL-22-HU55Z66A
	 */
	private void makePmrIdsMap() {
		IO.pl("\nMaking Pmr ID summary map...");
		pmrId_DatasetName = new TreeMap<String, TreeSet<String>>();

		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				String id = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId();

				addDataset(p.getCoreId(), id, pmrId_DatasetName);
				addDataset(d.getDatasetId(), id, pmrId_DatasetName); //for Avatar, Tempus, Caris from the file path

				//Avatar?
				if (d.getAvatarClinicalInfo()!=null) {
					AvatarClinicalInfo aci = d.getAvatarClinicalInfo();
					addDataset(aci.getPatient().get("hciPatientId"), id, pmrId_DatasetName);
					
					//Only need to set this for one of the following, they are all the same
					//Any Tumor DNA?
					if (aci.getTumorDna() != null) {
						addDataset(aci.getTumorDna().get("ORIENAvatarKey"), id, pmrId_DatasetName);
					}
					//Any Tumor RNA?
					else if (aci.getTumorRna()!= null) {
						addDataset(aci.getTumorRna().get("ORIENAvatarKey"), id, pmrId_DatasetName);
					}
					//any Normal DNA samples
					else if (aci.getNormalDna()!= null) {
						for (HashMap<String, String> n: aci.getNormalDna()) addDataset(n.get("ORIENAvatarKey"), id, pmrId_DatasetName);
					}
					
				}
			}
		}
		datasetIdSearchString = fetchKeysSearchString (pmrId_DatasetName.keySet());
	}
	
	/*Sex 
		M, F, or NA
	 */
	private void makeSexMap() {
		IO.pl("\nMaking Sex summary map...");
		sex_DatasetName = new TreeMap<String, TreeSet<String>>();

		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				String id = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId();
				//Avatar?
				if (d.getAvatarClinicalInfo()!=null) {
					AvatarClinicalInfo aci = d.getAvatarClinicalInfo();
					addDataset(aci.getPatient().get("sex"), id, sex_DatasetName);
				}
				//Caris?
				else if (d.getCarisClinicalInfo()!=null) {
					LinkedHashMap<String, String> cci = d.getCarisClinicalInfo();
					String sex = cci.get("gender");
					if (sex != null) {
						sex = sex.toUpperCase();
						if (sex.startsWith("F")) addDataset("F", id, sex_DatasetName);
						else if (sex.startsWith("M")) addDataset("M", id, sex_DatasetName);
						else addDataset("NA", id, sex_DatasetName);
					}
				}
				//Tempus pre V3?
				else if (d.getTempusJsonReportInfoV3()!=null) {
					TempusV3JsonSummary tjs = d.getTempusJsonReportInfoV3();
					String sex = tjs.getTempusPatient().getSex();
					if (sex != null) {
						if (sex.startsWith("F")) addDataset("F", id, sex_DatasetName);
						else if (sex.startsWith("M")) addDataset("M", id, sex_DatasetName);
						else addDataset("NA", id, sex_DatasetName);
					}
				}
			}
		}
		sexSearchString = fetchKeysSearchString (sex_DatasetName.keySet());
	}


	
	/*DiseaseGroup/Physician
		A='Disease Type', e.g. 'END - Thyroid Cancer'
		C=physicianName, e.g. 'Howard Colman'
		T=physician, e.g. 'Conan Kinsey'
	 */
	private void makeDiseaseGroupPhysician() {
		IO.pl("\nMaking Disease Group and Physician name summary map...");
		disGrpPhy_DatasetName = new TreeMap<String, TreeSet<String>>();
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				String id = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId();
				//Avatar?
				if (d.getAvatarClinicalInfo()!=null) {
					AvatarClinicalInfo aci = d.getAvatarClinicalInfo();
					//Any Tumor DNA?
					if (aci.getTumorDna() != null) {
						addDataset(aci.getTumorDna().get("Disease Type"), id, disGrpPhy_DatasetName);
					}
					//Any Tumor RNA?
					if (aci.getTumorRna()!= null) {
						addDataset(aci.getTumorRna().get("Disease Type"), id, disGrpPhy_DatasetName);
					}
				}
				//Caris?
				else if (d.getCarisClinicalInfo()!=null) {
					LinkedHashMap<String, String> cci = d.getCarisClinicalInfo();
					addDataset(cci.get("physicianName"), id, disGrpPhy_DatasetName);
				}
				//Tempus pre v3
				else if (d.getTempusJsonReportInfoPreV3()!=null) {
					TempusJsonSummary tjs = d.getTempusJsonReportInfoPreV3();
					addDataset(tjs.getTempusOrder().getPhysician(), id, disGrpPhy_DatasetName);
				}
				//Tempus v3
				else if (d.getTempusJsonReportInfoV3()!=null) {
					TempusV3JsonSummary tjs = d.getTempusJsonReportInfoV3();
					addDataset(tjs.getTempusOrder().getPhysician(), id, disGrpPhy_DatasetName);
				}
			}
		}
		disGrpPhySearchString = fetchKeysSearchString (disGrpPhy_DatasetName.keySet());
	}
	
	/*SpecimenSite
		A=SpecimenSiteOfCollection, e.g. 'Anterior wall of bladder'
		A=SpecimenSiteOfOrigin, e.g. 'Base of tongue, NOS'
		A=SpecimenSiteOfOriginRollUp, e.g. 'Anterior wall of bladder'
		C=specimenSite, e.g. 'Bile duct, NOS'
		T=sampleSite, e.g. 'Abdomen'
	 */
	private void makeSpecimenSiteMap() {
		IO.pl("\nMaking Specimen Site summary map...");
		specimenSite_DatasetName = new TreeMap<String, TreeSet<String>>();
		String[] avaToSearch = {"SpecimenSiteOfCollection", "SpecimenSiteOfOrigin", "SpecimenSiteOfOriginRollUp"};
		TreeSet<String> icdTxt = new TreeSet<String>();
		
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				String id = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId();
				//Avatar?
				if (d.getAvatarClinicalInfo()!=null) {
					AvatarClinicalInfo aci = d.getAvatarClinicalInfo();
					for (String keyToPull: avaToSearch) {
						//Any Tumor DNA?
						if (aci.getTumorDna() != null) {
							addDataset(aci.getTumorDna().get(keyToPull), id, specimenSite_DatasetName);
						}
						//Any Tumor RNA?
						if (aci.getTumorRna()!= null) {
							addDataset(aci.getTumorRna().get(keyToPull), id, specimenSite_DatasetName);
						}
					}
				}
				//Caris?
				else if (d.getCarisClinicalInfo()!=null) {
					LinkedHashMap<String, String> cci = d.getCarisClinicalInfo();
					addDataset(cci.get("specimenSite"), id, specimenSite_DatasetName);
				}
				//Tempus pre v3
				else if (d.getTempusJsonReportInfoPreV3()!=null) {
					TempusJsonSummary tjs = d.getTempusJsonReportInfoPreV3();
					for (TempusSpecimen ts: tjs.getTempusSpecimens()){
						addDataset(ts.getSampleSite(), id, specimenSite_DatasetName);
					}
				}
				//Tempus V3?
				else if (d.getTempusJsonReportInfoV3()!=null) {
					icdTxt.clear();
					TempusV3JsonSummary tjs = d.getTempusJsonReportInfoV3();
					for (TempusV3Specimen s : tjs.getTempusSpecimens()) {
						//primary sample site
						if (s.getPrimarySampleSite() != null) icdTxt.add(s.getPrimarySampleSite());
						
						//icd 0 topology
						if (s.getTempusIcdOCodeTopography() != null && icdOCodeTopology !=null) {
							for (String code: Misc.COMMA.split(s.getTempusIcdOCodeTopography())) {
								String decoded = icdOCodeTopology.get(code);
								if (decoded !=null) {
									icdTxt.add(decoded);
									s.getTempusIcdOTxtTopography().add(decoded);
								}
								else IO.el("\tFailed to find ICD Topology Code '"+code+"' Add it!");
							}
						}	
					}
					//add concatinate
					if (icdTxt.size() != 0) addDataset(Misc.treeSetToString(icdTxt, "; "), id, specimenSite_DatasetName);
				}
			}
		}
		specimenSiteSearchString = fetchKeysSearchString (specimenSite_DatasetName.keySet());
	}

	/*SpecimenId
	A=ORIENSpecimenID, e.g. 11-00101.4b
	C=specimenID, e.g. 20AA-IS-2562-A7
	T=blockId, e.g. 1A
    T=caseId, e.g. ARUP - Pathology FA-21-43328

	 */
	private void makeSpecimenIdMap() {
		IO.pl("\nMaking Specimen ID summary map...");
		specimenId_DatasetName = new TreeMap<String, TreeSet<String>>();

		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				String id = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId();
				//Avatar?
				if (d.getAvatarClinicalInfo()!=null) {
					AvatarClinicalInfo aci = d.getAvatarClinicalInfo();
					//Any Normal DNA samples?
					if (aci.getTumorDna() != null) {
						for (HashMap<String, String> n: aci.getNormalDna()) {
							addDataset(n.get("ORIENSpecimenID"), id, specimenId_DatasetName);
						}
					}
					//Any Tumor DNA?
					if (aci.getTumorDna() != null) {
						addDataset(aci.getTumorDna().get("ORIENSpecimenID"), id, specimenId_DatasetName);
					}
					//Any Tumor RNA?
					if (aci.getTumorRna()!= null) {
						addDataset(aci.getTumorRna().get("ORIENSpecimenID"), id, specimenId_DatasetName);
					}
				}
				//Caris?
				else if (d.getCarisClinicalInfo()!=null) {
					LinkedHashMap<String, String> cci = d.getCarisClinicalInfo();
					addDataset(cci.get("specimenID"), id, specimenId_DatasetName);
				}
				//Tempus pre v3
				else if (d.getTempusJsonReportInfoPreV3()!=null) {
					TempusJsonSummary tjs = d.getTempusJsonReportInfoPreV3();
					for (TempusSpecimen ts: tjs.getTempusSpecimens())addDataset(ts.getCaseId(), id, specimenId_DatasetName);
				}
				//Tempus v3
				else if (d.getTempusJsonReportInfoV3()!=null) {
					TempusV3JsonSummary tjs = d.getTempusJsonReportInfoV3();
					for (TempusV3Specimen ts: tjs.getTempusSpecimens())addDataset(ts.getCaseId(), id, specimenId_DatasetName);
				}
			}
		}
		specimenIdSearchString = fetchKeysSearchString (specimenId_DatasetName.keySet());
	}

	private void addDataset(String dtValue, String id, TreeMap<String, TreeSet<String>> map) {
		if (dtValue !=null && dtValue.length()!=0 && dtValue.toLowerCase().equals("null")==false) {
			TreeSet<String> datasets = map.get(dtValue);
			if (datasets == null) {
				datasets = new TreeSet<String>();
				map.put(dtValue, datasets);
			}
			datasets.add(id);
		}
	}

	private void summarizeAvatarClinicalInfo() {
		IO.pl("\nTabulating Avatar Clinical Info...");
		TreeMap<String, TreeMap<String, Integer>> patient = new TreeMap<String, TreeMap<String, Integer>>();
		TreeMap<String, TreeMap<String, Integer>> root = new TreeMap<String, TreeMap<String, Integer>>();
		TreeMap<String, TreeMap<String, Integer>> tumorDna = new TreeMap<String, TreeMap<String, Integer>>();
		TreeMap<String, TreeMap<String, Integer>> tumorRna = new TreeMap<String, TreeMap<String, Integer>>();
		TreeMap<String, TreeMap<String, Integer>> normalDna = new TreeMap<String, TreeMap<String, Integer>>();

		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				if (d.getAvatarClinicalInfo()!=null) {
					AvatarClinicalInfo aci = d.getAvatarClinicalInfo();
					loadTreeMap(aci.getPatient(), patient);
					loadTreeMap(aci.getRoot(), root);
					loadTreeMap(aci.getTumorDna(), tumorDna);
					loadTreeMap(aci.getTumorRna(), tumorRna);
					for (HashMap<String, String> n: aci.getNormalDna()) {
						loadTreeMap(n, normalDna);
					}
				}
			}
		}
		IO.pl("Avatar Patients:");
		printHashMapTreeMap(patient);	
		IO.pl("\nAvatar Root:");
		printHashMapTreeMap(root);
		IO.pl("\nAvatar NormalDna:");
		printHashMapTreeMap(normalDna);
		IO.pl("\nAvatar TumorDna:");
		printHashMapTreeMap(tumorDna);
		IO.pl("\nAvatar TumorRna:");
		printHashMapTreeMap(tumorRna);	
	}

	public void printHashMapTreeMap(TreeMap<String, TreeMap<String, Integer>> patient) {
		for (String key: patient.keySet()) {
			TreeMap<String, Integer> vals = patient.get(key);
			if (verbose) IO.pl("\t\""+key+"\" : "+vals);
			else {
				StringBuilder sb = new StringBuilder();
				Iterator<String> keys = vals.keySet().iterator();
				for (int i=0; i< 10; i++) {
					if (keys.hasNext()) {
						String k = keys.next();
						Integer count = vals.get(k);
						sb.append("\"");
						sb.append(k);
						sb.append("\"=");
						sb.append(count);
						sb.append(", ");
						if (sb.length()> 120) break;
					}
					else break;
				}
				if (keys.hasNext()) sb.append("...");
				IO.pl("\t\""+key+"\" : "+sb);
			}
		}
	}

	private void printTreeMapTreeSet(TreeMap<String, TreeSet<String>> trees) {
		for (String key: trees.keySet()) {
			TreeSet<String> values = trees.get(key);
			IO.p("\t"+key+" ("+values.size()+")"+" : ");
			Iterator<String> it =values.iterator();
			IO.p(it.next());
			while(it.hasNext()) IO.p(", "+it.next());
			IO.pl();
		}
	}

	/**Does a partial match to the keysToSearch (delimited by ` characters, e.g. `key1`key2`key3`key22`key4`) returning the full key for each partial match*/
	public static TreeSet<String> searchKeysString(String keysToSearch, Pattern pat){
		Matcher mat = pat.matcher(keysToSearch);
		TreeSet<String> keysFound = new TreeSet<String>();
		while (mat.find()) {
			String hit = mat.group();
			keysFound.add(hit.substring(1,hit.length()-1));
		}
		return keysFound;
	}
	public static String fetchKeysSearchString(Set<String> keySet) {
		StringBuilder sb = new StringBuilder("`");
		for (String key: keySet) {
			sb.append(key);
			sb.append("`");
		}
		return sb.toString();
	}
	public static Pattern fetchSearchPattern (String toFind, boolean caseInsensitive, boolean partialMatch) {
		Pattern pat = null;
		if (partialMatch) {
			if (caseInsensitive) pat = Pattern.compile("`[^`]*" +toFind+ "[^`]*`", Pattern.CASE_INSENSITIVE);
			else pat = Pattern.compile("`[^`]*" +toFind+ "[^`]*`");
		}
		else {
			if (caseInsensitive) pat = Pattern.compile("`" +toFind+ "`", Pattern.CASE_INSENSITIVE);
			else pat = Pattern.compile("`" +toFind+ "`");
		}
		return pat;
	}

	private void loadTreeMap(HashMap<String, String> aci, TreeMap<String,TreeMap<String, Integer>> sum) {
		//for each key from a particular dataset
		//e.g. SpecimenDerivativeSource : Blood
		for (String aciKey: aci.keySet()) {
			//see if it is in the big list
			TreeMap<String, Integer> sumKey = sum.get(aciKey);
			if (sumKey == null) {
				sumKey = new TreeMap<String, Integer>();
				sum.put(aciKey, sumKey);
			}
			Integer count = sumKey.get(aci.get(aciKey));
			if (count == null) count = new Integer(1);
			else count = new Integer(count.intValue()+1);
			sumKey.put(aci.get(aciKey), count);
		}
	}

	private void loadAvatarClinicalJsons() {
		IO.pl("Loading Avatar Clinical Json files...");
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				if (d.getSource().equals(avatarSourceName)) {
					if (d.getClinicalInfoFiles() == null) {
						//IO.el("WARNING! "+p.getCoreId()+" "+d.getDatasetId()+" Avatar dataset is missing a json!");
					}
					else if (d.getClinicalInfoFiles().size() !=1) IO.el("WARNING! "+p.getCoreId()+" "+d.getDatasetId()+" Avatar dataset has more than one json!");
					else {
						AvatarClinicalInfo aci = loadAvatarJson(d.getClinicalInfoFiles().get(0));
						d.setAvatarClinicalInfo(aci);
					}
				}
			}
		}
	}

	private void loadTempusClinicalJsons() throws Exception {		
		ArrayList<File> toParseOri = new ArrayList<File>();
		ArrayList<File> toParseV3 = new ArrayList<File>();
		ArrayList<Dataset> toSetOri = new ArrayList<Dataset>();
		ArrayList<Dataset> toSetV3 = new ArrayList<Dataset>();
		
		//for each patient
		for (Patient p: idPatient.values()) {
			//look if it has a tempus dataset
			for (Dataset d: p.getIdDataSets().values()) {
				
				//is it Tempus?
				if (d.getSource().equals(tempusSourceName)) {
					
					//any clinical json reports?
					if (d.getClinicalInfoFiles() == null) IO.el("WARNING! "+p.getCoreId()+" "+d.getDatasetId()+" Tempus dataset is missing a json!");

					//yes some reports
					else {
						//not sure this is going to work
						for (File f: d.getClinicalInfoFiles()) {
							String version = fetchTempusSchemaVersion(f);
							if (version == null) throw new IOException("ERROR: failed to parse the schema version from "+f);
							if (version.startsWith("3")) {
								toParseV3.add(f);
								toSetV3.add(d);
							}
							else {
								toParseOri.add(f);
								toSetOri.add(d);
							}
						}
					}
				}
			}
		}
		//parse all of the patient tempus reports
		File[] filesOri = new File[toParseOri.size()];
		toParseOri.toArray(filesOri);
		
		File[] filesV3 = new File[toParseV3.size()];
		toParseV3.toArray(filesV3);

		//parse the json report files, this also prints out an aggregate summary
		IO.pl("Parsing Tempus pre V3...");
		TempusJson2Vcf tj = new TempusJson2Vcf(filesOri);

		//set the summaries in each of the patients tempus datasets
		ArrayList<TempusJsonSummary> summariesOri = tj.getSummaries();
		if (toSetOri.size() != summariesOri.size()) throw new Exception("\nMismatch in Dataset number and TempusJsonSummaryOri number");
		for (int i=0; i< toSetOri.size(); i++) toSetOri.get(i).setTempusJsonReportInfoPreV3(summariesOri.get(i));
		
		IO.pl("Parsing Tempus V3...");
		TempusV3Json2Vcf tj3 = new TempusV3Json2Vcf(filesV3);
		//set the summaries in each of the patients tempus datasets
		ArrayList<TempusV3JsonSummary> summariesV3 = tj3.getSummaries();
		if (toSetV3.size() != summariesV3.size()) throw new Exception("\nMismatch in Dataset number and TempusJsonSummaryV3 number");
		for (int i=0; i< toSetV3.size(); i++) {
			toSetV3.get(i).setTempusJsonReportInfoV3(summariesV3.get(i));
		}
	}
	
	public static final Pattern schemaCleaner = Pattern.compile("[\\sv,\"]");
	public static String fetchTempusSchemaVersion(File json) throws IOException {
		//"schemaVersion": "v3.3.0",
		//"schemaVersion": "v2.0.0",
		//"schemaVersion": "1.4.2",
        BufferedReader in = IO.fetchBufferedReader(json);
        String line = null;
        String version = null;	
        while ((line = in.readLine())!=null) {
        	if (line.contains("schemaVersion")) {
        		String[] split = Misc.COLON.split(line);
        		if (split.length !=2) throw new IOException("ERROR: failed to find two fields split by : with schema, "+line+" from "+json);
        		version = schemaCleaner.matcher(split[1]).replaceAll("");
        		break;
        	}
        }
        in.close();
        return version;
	}

	private void fetchTempusClinicalJsons() throws Exception {
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				//IO.pl(d.getSource()+" -> "+d.getDatasetId());
				if (d.getSource().equals(tempusSourceName)) {
					for (String partPath: d.getPartialPaths()) {
						//IO.pl(partPath);
						// there are often 2 jsons for each Tempus (one DNA, one RNA)
						if (partPath.contains(clinicalReportDirName) && partPath.endsWith("json") && partPath.contains("_RS.")==false) {
							String relativePath = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId()+"/"+partPath;
							File j = new File (clinicalReportDir, relativePath);
							if (j.exists() == false) download(awsPatientDirUri+relativePath, j);
							d.setClinicalInfoFile(j);
						}
					}
				}
			}
		}
	}

	private void fetchCarisClinicalXmls() throws Exception {
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				//IO.pl(d.getSource()+" -> "+d.getDatasetId());
				if (d.getSource().equals(carisSourceName)) {
					for (String partPath: d.getPartialPaths()) {
						//IO.pl(partPath);
						// there are often 2 jsons for each Tempus (one DNA, one RNA)
						if (partPath.contains(clinicalReportDirName) && partPath.endsWith("xml") && partPath.contains("_RS.")==false) {
							String relativePath = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId()+"/"+partPath;
							File j = new File (clinicalReportDir, relativePath);
							if (j.exists() == false) download(awsPatientDirUri+relativePath, j);
							d.setClinicalInfoFile(j);
						}
					}
				}
			}
		}
	}

	private void loadCarisClinicalXmls() throws Exception {
		IO.pl("Loading Caris clinical xml files...");
		System.out.flush();
		ArrayList<File> toParse = new ArrayList<File>();
		ArrayList<Dataset> toSet = new ArrayList<Dataset>();
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				if (d.getSource().equals(carisSourceName)) {
					if (d.getClinicalInfoFiles() == null) IO.el("\tWARNING! "+p.getCoreId()+" "+d.getDatasetId()+" Caris dataset is missing an xml, skipping!");
					else if (d.getClinicalInfoFiles().size() !=1) IO.el("\tWARNING! "+p.getCoreId()+" "+d.getDatasetId()+" Caris dataset has more than one standard test xml!");
					else {
						toParse.add(d.getClinicalInfoFiles().get(0));
						toSet.add(d);
					}
				}
			}
		}
		System.err.flush();
		File[] files = new File[toParse.size()];
		toParse.toArray(files);

		//parse the json report files, this also prints out an aggregate summary
		CarisXmlVcfParser tj = new CarisXmlVcfParser(files, icd10CodeDesc, icdOCodeMorphology, icdOCodeTopology);
		
		LinkedHashMap<String, String>[] attributes = tj.getAllReportAttributes();

		//set the clinical info in each of the datasets
		for (int i=0; i< toSet.size(); i++) toSet.get(i).setCarisXmlReportInfo(attributes[i]);
	}

	private void summarizeCarisClinicalInfo() {
		IO.pl("\nTabulating Caris Clinical Info...");
		TreeMap<String, TreeMap<String, Integer>> summary = new TreeMap<String, TreeMap<String, Integer>>();

		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				if (d.getCarisClinicalInfo()!=null) {
					LinkedHashMap<String, String> aci = d.getCarisClinicalInfo();
					loadTreeMap(aci, summary);
				}
			}
		}
		printHashMapTreeMap(summary);		
	}

	private AvatarClinicalInfo loadAvatarJson(File file) {
		return new AvatarClinicalInfo(file);

	}

	private void fetchAvatarClinicalJsons() throws Exception {
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				//IO.pl(d.getSource()+" -> "+d.getDatasetId());
				if (d.getSource().equals(avatarSourceName)) {
					for (String partPath: d.getPartialPaths()) {
						//IO.pl(partPath);
						if (partPath.contains(clinicalReportDirName) && partPath.endsWith("json")) {
							String relativePath = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId()+"/"+partPath;
							File j = new File (clinicalReportDir, relativePath);
							if (j.exists() == false) download(awsPatientDirUri+relativePath, j);

							d.setClinicalInfoFile(j);
						}
					}
				}
			}
		}
	}

	
	
	private void download(String s3Uri, File j) throws Exception {
		if (verbose) IO.pl("\tDownloading "+s3Uri+ " -> "+j);
		String cmd = awsPath + " s3 cp "+s3Uri+ " "+ j.getCanonicalPath()+ " --profile "+profile;
		cmdsForParallel.add(cmd);
	}

	private void loadPatientFiles() throws IOException {
		IO.pl("Loading Patients...");
		BufferedReader in = IO.fetchBufferedReader(awsRepoList);
		String line = null;
		String[] fields = null;
		String[] keys = null;
		while ((line = in.readLine())!=null) {
			//2022-05-31   09:04:43   2938    Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
			//  date         time     size       key
			//    0            1        2         3
			fields = Misc.TAB.split(line);
			if (fields.length != 4) throw new IOException("Failed to find 4 fields in "+line);
			
			//date filter?
			if (maximumDaysOld != 0) {
				LocalDate inputDate = LocalDate.parse(fields[0], formatter);
		        long daysOld = ChronoUnit.DAYS.between(inputDate, todaysDate);
		        if (daysOld > maximumDaysOld) continue;
			}
			

			//Patients   AA2mF6Vy   Avatar   A032049_SL419345_SL419548_SL420681     ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
			//   0          1         2                 3                                  4
			keys = Misc.FORWARD_SLASH.split(fields[3]);
			if (keys.length < 5) {
				if (line.contains("blacklisted")) continue;
				else throw new IOException("Failed to find more then 4 dirs in "+line);
			}

			Patient p = fetchPatient(keys[1]);

			p.addDataFile(keys);


		}
		in.close();

		//some stats
		loadSourceDatasetIds();
		loadSourceCombinations();
		IO.pl("\nSummary Statistics:");
		IO.pl("\t"+idPatient.size()+"\tPatients\n");
		if (verbose) for (String source: sourceDatasets.keySet()) IO.pl("\t"+sourceDatasets.get(source).size()+"\t"+source+"\t"+sourceDatasets.get(source));
		else for (String source: sourceDatasets.keySet()) IO.pl("\t"+sourceDatasets.get(source).size()+"\t"+source);
		IO.pl("\n\tCombinations:");
		if (verbose) for (String sourceCombo: sourceCombinations.keySet()) IO.pl("\t"+sourceCombinations.get(sourceCombo).size()+"\t"+sourceCombo+"\t"+sourceCombinations.get(sourceCombo));
		else for (String sourceCombo: sourceCombinations.keySet()) IO.pl("\t"+sourceCombinations.get(sourceCombo).size()+"\t"+sourceCombo);
	}

	private void loadSourceDatasetIds() {
		if (sourceDatasets != null) return;
		IO.pl("Loading Source Datasets...");
		sourceDatasets = new TreeMap<String, HashSet<String>>();
		idDataset = new HashMap<String, Dataset>();
		
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				String id = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId();
				idDataset.put(id, d);
				String source = d.getSource();
				HashSet<String> datasetIds = sourceDatasets.get(source);
				if (datasetIds == null) {
					datasetIds = new HashSet<String>();
					sourceDatasets.put(source, datasetIds);
				}
				datasetIds.add(d.getDatasetId());
			}
		}
	}

	private void loadSourceCombinations() {
		if (sourceCombinations != null) return;
		IO.pl("Loading Source Combinations...");

		sourceCombinations = new TreeMap<String, TreeSet<String>>();
		for (Patient p: idPatient.values()) {
			TreeSet<String> sources = new TreeSet<String>();
			for (Dataset d: p.getIdDataSets().values()) sources.add(d.getSource());
			String sourceKey = Misc.treeSetToString(sources, ",");

			TreeSet<String> patientIds = sourceCombinations.get(sourceKey);
			if (patientIds == null) {
				patientIds = new TreeSet<String>();
				sourceCombinations.put(sourceKey, patientIds);
			}
			patientIds.add(p.getCoreId());
		}
	}

	private Patient fetchPatient(String id) {
		Patient p = idPatient.get(id);
		if (p == null) {
			p = new Patient(id);
			idPatient.put(id, p);
		}
		return p;
	}

	/*Looks for the repo file list file, downloads it if not found and saves it.*/
	private void downloadAwsFileList() throws IOException {
		awsRepoList = new File (clinicalReportDir, "awsRepoList.txt");
		if (awsRepoList.exists() == false) fetchFilesInRepo();
		else IO.pl("\nWARNING! Using existing AWS file list " +awsRepoList+"! Delete it and restart to pull the latest.\n");
	}

	/**Uses ProcessBuilder to execute a cmd
	 * Returns exit code, 0=OK, >0 a problem
	 * @throws IOException */
	public int executeReturnExitCode(String[] command) throws Exception{
		if (verbose) IO.pl ("Executing: "+Misc.stringArrayToString(command, " "));
		ProcessBuilder pb = new ProcessBuilder(command);
		//add enviro props?
		pb.environment().putAll(envPropToAdd);
		pb.redirectErrorStream(true);
		Process proc = pb.start();
		return proc.waitFor();
	}

	/**Uses ProcessBuilder to execute a cmd, combines standard error and standard out into one and returns their output.
	 * @throws IOException */
	public String[] executeViaProcessBuilder(String[] command, boolean printStandardOut, File workingDirectory) throws IOException{
		if (verbose) IO.pl ("Executing: '"+Misc.stringArrayToString(command, " ")+"'");
		ArrayList<String> al = new ArrayList<String>();
		ProcessBuilder pb = new ProcessBuilder(command);
		//add enviro props?
		pb.environment().putAll(envPropToAdd);
		if (workingDirectory !=null) pb.directory(workingDirectory);

		pb.redirectErrorStream(true);
		Process proc = pb.start();
		BufferedReader data = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String line;
		while ((line = data.readLine()) != null){			
			al.add(line);
			if (printStandardOut) IO.pl(line);
		}
		data.close();
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}

	private void fetchFilesInRepo() throws IOException{
		IO.pl("\nFetching file list from AWS for '"+awsPatientDirUri+ "'...");
		String[] cmd = {awsPath, "s3", "ls", "--recursive", awsPatientDirUri, "--profile", profile};
		String[] res = executeViaProcessBuilder(cmd, false, null);
		PrintWriter out = new PrintWriter(new FileWriter(awsRepoList));
		for (String s: res) {
			s= s.trim();
			//empty
			if (s.length() == 0) continue;
			//folder?
			if (s.startsWith("PRE ")) continue;
			String[] fields = Misc.WHITESPACE.split(s);
			//2023-03-01   09:04:43   2938    Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
			//   date         time     size       key
			out.println(Misc.stringArrayToString(fields, "\t"));
		}
		out.close();
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PMRSearch(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException {

		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		File icd10File = null;
		File icdMorphologyFile = null;
		File icdTopologyFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': clinicalReportDir = new File(args[++i]).getCanonicalFile(); break;
					case 'a': awsPatientDirUri = args[++i]; break;
					case 'p': profile =args[++i]; break;
					case 'i': icd10File = new File(args[++i]); break;
					case 'm': icdMorphologyFile = new File(args[++i]); break;
					case 't': icdTopologyFile = new File(args[++i]); break;
					case 's': inputNonInteractive = args[++i]; break;
					case 'x': maximumDaysOld = Long.parseLong(args[++i]); break;
					case 'v': verbose = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//icd files?
		if (icd10File != null) {
			IO.pl("Loading ICD 10 codes for Diagnosis...");
			icd10CodeDesc = IO.loadFileIntoHash(icd10File, 0, 1);
			if (icd10CodeDesc == null) throw new IOException("ERROR: failed to pars ICD Diagnosis file "+icd10File);
		}
		if (icdMorphologyFile != null) {
			IO.pl("Loading ICD morphology codes for Diagnosis...");
			icdOCodeMorphology = IO.loadFileIntoHash(icdMorphologyFile, 0, 1);
			if (icdOCodeMorphology == null) throw new IOException("ERROR: failed to pars ICD morphology file "+icdMorphologyFile);
		}
		if (icdTopologyFile != null) {
			IO.pl("Loading ICD topology codes for SpecimenSites...");
			icdOCodeTopology = IO.loadFileIntoHash(icdTopologyFile, 0, 1);
			if (icdOCodeTopology == null) throw new IOException("ERROR: failed to pars ICD topology file "+icdTopologyFile);
		}

		//root patient dirs?
		if (clinicalReportDir == null ) {
			Misc.printErrAndExit("Error: failed to find your directory for saving ClinicalReport json and xml files, see -c ? "+clinicalReportDir);
		}
		clinicalReportDir.mkdirs();

		//check credentials
		File awsCredentialsDir = new File (System.getProperty("user.home")+"/.aws");
		if (awsCredentialsDir.exists() == false)  Misc.printErrAndExit("Error: failed to find the aws credentials directory "+awsCredentialsDir);
		File credentialsFile = new File (awsCredentialsDir, "credentials").getCanonicalFile();
		if (awsCredentialsDir.exists() == false)  Misc.printErrAndExit("Error: failed to find your aws credentials file "+credentialsFile);


		//set env var for the aws cli 
		envPropToAdd.put("AWS_SHARED_CREDENTIALS_FILE", credentialsFile.getCanonicalPath());

		//Only needed in Eclipse
		/*if (true) {
			envPropToAdd.put("PATH", "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin");
			awsPath="/usr/local/bin/aws";
		}*/


	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                       Patient Molecular Repo Search : Oct 2025                   **\n" +
				"**************************************************************************************\n" +
				"Interactive searching of the clinical and sample attribute information in the json/xml\n"+
				"reports in the HCI PMR /ClinicalReport/ folders to identify datasets for analysis.\n"+
				"Press return after loading to see the search menu. Assumes the AWS CLI:\n"+
				"https://docs.aws.amazon.com/cli/ is installed in your path and you have read access\n"+
				"to the PMR. For searching for particular mutations, fusions, or cnvs, use GQuery:\n"+
				"https://github.com/HuntsmanCancerInstitute/GQuery\n"+

				"\nOptions:\n"+
				"-d  Directory to save the PHI redacted clinical xml and json reports.\n"+
				"-a  S3 URI containing the patient molecular repo, defaults to\n"+
				"      s3://hcibioinfo-patient-molecular-repo/Patients/ \n"+
				"-p  AWS credential profile, defaults to 'default'\n"+
				"-i  File containing ICD 10 diagnosis codes and their descriptors. Tab delimited. For \n"+
				"       converting Tempus diagnosis codes to searchable text.\n"+
				"-m  File containing ICD 0 morphology codes. Ditto.\n"+
				"-t  File containing ICD 0 topology codes. Ditto.\n"+
				"-s  Run non interactive search using this input then exit. Surround with \"  \".\n"+
				"       Separate multiple searches with ;\n"+
				"-x  Maximum days old to parse, defaults to all.\n"+
				"-v  Verbose output.\n"+
				
				"\nSearch Options, enter at the interactive prompt or using the -s 'xxxx' option above:\n"+
				"-d  Data table to search, choose from: Diagnosis, PhysicianDiseaseGroups, SpecimenSites,\n"+
				"       SpecimenIds, Sex, PmrIds\n"+
				"-s  Search terms, comma delimited, no spaces. Surround phrases with single quotes.\n"+
				"-c  Terms are case sensitive.\n"+
				"-e  Terms are exact, no partial matching.\n"+
				"-m  All terms must match.\n"+
				"-a  Add matches to the prior result set.\n"+
				"-x  New search, clear any prior results.\n"+
				"-p  Print the contents of the named data table, see -d\n"+
				"-f  Print the available file AWS URIs for the matched datasets.\n"+
				"-i  Print the clinical and test details for the matched datasets.\n"+
				"-n  Print the matched dataset names.\n"+

				"\nExample: java -jar pathToUSeq/Apps/PMRSearch -d ~/PMRFiles/ -x 365 -i \n"+
				"   ~/ICD/ICD-10_Diagnosis.txt -m ~/ICD/ICD_Morphology.txt -t ~/ICD/ICD_Topology.txt\n"+

				"**************************************************************************************\n");
	}
}
