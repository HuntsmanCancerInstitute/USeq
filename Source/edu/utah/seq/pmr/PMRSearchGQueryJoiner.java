package edu.utah.seq.pmr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.json.JSONArray;
import org.json.JSONObject;

import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import util.gen.IO;
import util.gen.Misc;

import org.apache.poi.ss.usermodel.*;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**Merges the results from PMRSearches and GQueries into a spreadsheet report.*/
public class PMRSearchGQueryJoiner {
	
	//User provided fields
	private File pmrSearchResultDir = null;
	private File gqueryResultFile = null;
	private File patientRegistryDirectory = null;
	private File outputDirectory = null;
	private File priorReportedPmrIds = null;
	
	//Internal fields
    private String pmrResultsStartLine = "Dataset Info, also download";
	private ArrayList<String> pmrDataLines = new ArrayList<String>();

	private HashMap<String, ArrayList<String>> pmrIdPmrDatasets = new HashMap<String, ArrayList<String>>();
	private HashMap<String, ArrayList<String>> pmrIdGQueryDatasets = new HashMap<String, ArrayList<String>>();
	private HashMap<String, ArrayList<String>> searchIdResults = new HashMap<String, ArrayList<String>>();
	private HashMap<String, TreeSet<String>> sourceVcfs = new HashMap<String, TreeSet<String>>();
	private HashSet<String> pmrIdsInCommon = new HashSet<String>();
	private HashMap<String,String[]> pmrIdsPhi = null;

	public PMRSearchGQueryJoiner (String[] args) {
		long startTime = System.currentTimeMillis();
		try {
			
			processArgs(args);
			
			loadPMRSearchResults();
			if (pmrDataLines.size()==0) {
				IO.pl("\tNo results from PMRSearch to parse, exiting!");
				return;
			}

			parseGQueryJsonFile();
			if (pmrIdGQueryDatasets.size()==0) {
				IO.pl("\tNo results from GQuery to parse, exiting!");
				return;
			}
			
			parsePMRSearchResultsAll();
				
			intersectResults();
			
			loadPhi();
			
			printResultsToXlsxFile();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
			
		} catch (Exception e) {
			IO.el("\nERROR running the PMRSearchGQueryJoiner...");
			e.printStackTrace();
			System.exit(1);
		}
	}
	

	private void loadPhi() throws IOException {
		if (pmrIdsInCommon.size()==0 || patientRegistryDirectory ==null) return;
		
		IO.pl("\nFetching PHI for matching patients...");
		
		//make IO for writing out queries
		File queryFile = new File (outputDirectory, "tmpDelme_PmrIdsToFetch.txt");
		File queryResDir = new File (outputDirectory, "TmpDelme_SMMQueryResults_PHI");
		queryResDir.mkdirs();
		
		PrintWriter out = new PrintWriter (new FileWriter( queryFile));
		for (String s : pmrIdsInCommon) out.println(s);
		out.close();
		
		//run the SMM
		String[] args = {
				"-r", patientRegistryDirectory.getCanonicalPath(),
				"-q", queryFile.getCanonicalPath(),
				"-o", queryResDir.getCanonicalPath(),
				"-v"
		};
		
		SubjectMatchMaker smm = new SubjectMatchMaker(args);
		
		//in same order as that which was printed out above from workingPatients
		Subject[] queries = smm.getQuerySubjects();
		
		pmrIdsPhi = new HashMap<String, String[]>();
		for (Subject s: queries) {
			String[] phi = new String[] {s.getCoreId(), s.getLastName(), s.getFirstName(), s.getMrn(), s.getDobMonth()+"/"+s.getDobDay()+"/"+s.getDobYear()};
			pmrIdsPhi.put(s.getCoreId(), phi);
		}
		//cleanup
		queryFile.delete();
		IO.deleteDirectory(queryResDir);
	}

	private void parseGQueryJsonFile() throws IOException {
		IO.pl("Parsing GQuery results...");

		String jString = IO.loadFile(gqueryResultFile, " ", true);
		JSONObject mainJsonObject = new JSONObject(jString);
		JSONArray queryResults = mainJsonObject.getJSONArray("queryResults");
		Pattern gqVcf = Pattern.compile(".+/([^_]+)_.+vcf.gz");
		Matcher mat = null;

		//for each query, typically a genomic region, e.g. chrX:123456-1234567
		int len = queryResults.length();
		for (int i=0; i< len; i++) {
			//hits - sources
			JSONArray hits = queryResults.getJSONObject(i).getJSONArray("hits");
			int numHits = hits.length();
			for (int j=0; j< numHits; j++) {
				JSONObject hit = hits.getJSONObject(j);

				// Source 
				// Data/Hg38/Somatic/Tempus/Vcf/smW7PN4NYk_25gcijsr_20250914_Anno_Hg38.anno.vcf.gz
				String sourceName = hit.getString("source");
				//parse patient PMR ID
				mat = gqVcf.matcher(sourceName);
				if (mat.matches()) {
					String pmrId = mat.group(1);
					ArrayList<String> al = pmrIdGQueryDatasets.get(pmrId);
					if (al == null) {
						al = new ArrayList<String>();
						pmrIdGQueryDatasets.put(pmrId, al);
					}
					al.add(sourceName);
				}	
				else throw new IOException("No gquery sources match for "+sourceName);

				//fetch or make the vars container
				TreeSet<String> allVars = sourceVcfs.get(sourceName);
				if (allVars == null) {
					allVars = new TreeSet<String>();
					sourceVcfs.put(sourceName, allVars);
				}

				// Add Vcf records
				// ["chr12\t25245350\tStrelka_22;Tempus_0;12583\tC\tA\t124.77\tPASS\tBKZ=124.77;BKAF=0.0015,0.001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;T_DP=1011;T_AF=0.0406;N_DP=350;N_AF=0;SOMATIC;QSS=131;TQSS=1;NT=ref;QSS_NT=131;TQSS_NT=1;SGT=CC->AC;MQ=60.00;MQ0=0;ReadPosRankSum=-0.55;SNVSB=0.00;SomaticEVS=10.23;EG=KRAS;CL=potentiallyActionable;FE=Missense_variant_(exon_2)_-_GOF;AF=0.04;ANN=A|missense_variant|MODERATE|KRAS|KRAS|transcript|NM_033360.4|protein_coding|2/6|c.35G>T|p.Gly12Val|225/5430|35/570|12/189||,A|missense_variant|MODERATE|KRAS|KRAS|transcript|NM_004985.5|protein_coding|2/5|c.35G>T|p.Gly12Val|225/5306|35/567|12/188||;ALLELEID=27622;CLNDISDB=Human_Phenotype_Ontology:HP:0010815,MedGen:C3854181|Human_Phenotype_Ontology:HP:0010817,MONDO:MONDO:0008097,MedGen:C4552097,OMIM:163200,Orphanet:2612|MONDO:MONDO:0006279,MedGen:C1708781|Human_Phenotype_Ontology:HP:0002408,MONDO:MONDO:0007154,MedGen:C0917804,OMIM:108010,Orphanet:46724|Human_Phenotype_Ontology:HP:0030358,MONDO:MONDO:0005233,MeSH:D002289,MedGen:C0007131|MONDO:MONDO:0005192,MeSH:C562463,MedGen:C0235974,Orphanet:1333,Orphanet:217074|Human_Phenotype_Ontology:HP:0012209,MONDO:MONDO:0011908,MedGen:C0349639,OMIM:607785,Orphanet:86834|MONDO:MONDO:0021060,MedGen:C5555857,Orphanet:536391|MedGen:C3661900|Human_Phenotype_Ontology:HP:0005506,Human_Phenotype_Ontology:HP:0005544,MONDO:MONDO:0011996,MeSH:D015464,MedGen:C0279543,OMIM:608232,Orphanet:521;CLNDN=Nevus_sebaceous|Linear_nevus_sebaceous_syndrome|Lung_sarcomatoid_carcinoma|Cerebral_arteriovenous_malformation|Non-small_cell_lung_carcinoma|Carcinoma_of_pancreas|Juvenile_myelomonocytic_leukemia|RASopathy|not_provided|Chronic_myelogenous_leukemia,_BCR-ABL1_positive;CLNHGVS=NC_000012.12:g.25245350C>A;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=COSMIC:520|ClinGen:CA122540|OMIM:190070.0006|OMIM:190070.0026|UniProtKB:P01116#VAR_006840;GENEINFO=KRAS:3845;MC=SO:0001583|missense_variant;ONC=Oncogenic;ONCDISDB=Human_Phenotype_Ontology:HP:0002664,Human_Phenotype_Ontology:HP:0003008,Human_Phenotype_Ontology:HP:0006741,MONDO:MONDO:0005070,MeSH:D009369,MedGen:C0027651;ONCDN=Neoplasm;ONCREVSTAT=criteria_provided,_single_submitter;ORIGIN=3;RS=121913529;MUTATION_EFFECT=Gain-of-function;ONCOGENIC=Oncogenic;HIGHEST_LEVEL=LEVEL_3B"],
				JSONArray vcfs = hit.getJSONArray("data");
				int numVcfs = vcfs.length();
				for (int v=0; v< numVcfs; v++) allVars.add(vcfs.get(i).toString());

			}
		}
	}
	

	private void intersectResults() {
		IO.pl("Intersecting results...");
		for (String g: pmrIdGQueryDatasets.keySet()) {
			if (pmrIdPmrDatasets.containsKey(g)) pmrIdsInCommon.add(g);
		}
		IO.pl("\t"+pmrIdsInCommon.size()+"\tPatients found in both");
		
		HashSet<String> pmrIds2Skip = new HashSet<String>();
		if (priorReportedPmrIds != null) {
			IO.pl("Removing those that were already reported...");
			pmrIds2Skip = IO.loadFileIntoHashSet(priorReportedPmrIds);
			pmrIdsInCommon.removeAll(pmrIds2Skip);
			IO.pl("\t"+pmrIdsInCommon.size()+"\tPatients remaining");
		}
		pmrIds2Skip.addAll(pmrIdsInCommon);
		IO.writeHashSet(pmrIds2Skip, new File (outputDirectory, "combinePmrIdsReported.txt"));
		
	}
	
	private void printResultsToXlsxFile() throws IOException {
		IO.pl("Printing xlsx spreadsheet results...");
		
		Workbook workbook = new XSSFWorkbook();
        Font boldFont = workbook.createFont();
        boldFont.setBoldweight(boldFont.BOLDWEIGHT_BOLD);
        CellStyle boldStyle = workbook.createCellStyle();
        boldStyle.setFont(boldFont);
        
        Sheet sheet = workbook.createSheet("Joint Query Results");
		
        //add header
        String[] header = null;
		if (patientRegistryDirectory != null) {
			header = new String[] {"PmrID", "LastName", "FirstName","MRN","DoB","Use the toggles on the left to view clinical and variant info"};
		}
		else header = new String[] {"PmrID"};
		int counter = 0;
		Row row = sheet.createRow(counter++);
		for (int i=0; i< header.length; i++) {
			Cell cell = row.createCell(i);
			cell.setCellValue(header[i]);
			cell.setCellStyle(boldStyle);
		}
		// and freeze it
		sheet.createFreezePane(0, 1);
		
		//create an empty row
		sheet.createRow(counter++);
        
		//for each PmrId
		for (String i: pmrIdsInCommon) {
			Row rowD = sheet.createRow(counter++);
			
			//fetch PHI?
			String[] phi = null;
			if (patientRegistryDirectory == null) phi = new String[] {"Patient:", i};
			else phi = pmrIdsPhi.get(i);
			for (int x=0; x< phi.length; x++) {
				Cell c = rowD.createCell(x);
				c.setCellValue(phi[x]);
				c.setCellStyle(boldStyle);
			}
				
			int startGroup = counter;
			
			ArrayList<String> searchResultIds = pmrIdPmrDatasets.get(i);
			for (int r = 0; r< searchResultIds.size(); r++) {
				String sri = searchResultIds.get(r);
				
				//clinicals
				Row rowC = sheet.createRow(counter++);
				rowC.createCell(0).setCellValue("Clinical:"); rowC.createCell(1).setCellValue(sri);
				ArrayList<String> searchLines = searchIdResults.get(sri);
				ArrayList<String[]> results = null;
				if (sri.contains("Avatar")) results = filterSearchLinesAvatarSplit(searchLines);
				else if (sri.contains("Tempus")) results = filterSearchLinesTempusSplit(searchLines);
				else if (sri.contains("Caris")) results = filterSearchLinesCarisSplit(searchLines);
				for (String[] s: results) {
					Row rowS = sheet.createRow(counter++);
					rowS.createCell(1).setCellValue(s[0]);
					rowS.createCell(2).setCellValue(s[1]);
				}
			
				//variants
				for (String gd: pmrIdGQueryDatasets.get(i)) {
					Row rowV = sheet.createRow(counter++);
					rowV.createCell(0).setCellValue("Variant:"); rowV.createCell(1).setCellValue(gd);
					TreeSet<String> vars = sourceVcfs.get(gd);
					for (String vcf: vars) {
						String[] splitVcf = Misc.TAB.split(vcf);
						//out.println("\t"+vcf);
						Row rowSV = sheet.createRow(counter++);
						for (int a = 0; a< splitVcf.length; a++) rowSV.createCell(a+1).setCellValue(splitVcf[a]);
					}
				}
			}
			//create an empty row
			sheet.createRow(counter++);
			
			int stopGroup = counter-1;
			sheet.groupRow(startGroup, stopGroup);
		}
		
		//save it
		File f = null;
		if (patientRegistryDirectory!=null) f= new File(outputDirectory, "results_PHI.xlsx");
		else f= new File(outputDirectory, "results.xlsx");
		FileOutputStream fileOut = new FileOutputStream(f);
        workbook.write(fileOut);
        fileOut.close();
        
	}


	private static final String[] carisSearchTerms = {"File(s) :", "labReportID :", "physicianName :", "orderedDate :", "specimenCollectionDate :", "specimenSite :", "primarySite :", "lineage :", "subLineage :", "diagnosis :", "icd_code :"}; 	
	private ArrayList<String> filterSearchLinesCaris(ArrayList<String> searchLines) {
		ArrayList<String> toReturn = new ArrayList<String>();
		//for each search term
		for (String st: carisSearchTerms) {
			//for each data line
			for (String dl: searchLines) {
				if (dl.startsWith(st)) {
					String tabbed = dl.replaceFirst(" : ", ":\t");
					toReturn.add(tabbed);
					break;
				}
			}
		}
		return toReturn;
	}
	private ArrayList<String[]> filterSearchLinesCarisSplit(ArrayList<String> searchLines) throws IOException {
		ArrayList<String[]> toReturn = new ArrayList<String[]>();
		//for each search term
		for (String st: carisSearchTerms) {
			//for each data line
			for (String dl: searchLines) {
				if (dl.startsWith(st)) {
					String[] splitLine = dl.split(" : ",2);
					if (splitLine.length != 2) throw new IOException("Failed to split in two: "+dl);
					toReturn.add(splitLine);
					break;
				}
			}
		}
		return toReturn;
	}

	private static final String[] tempusSearchTerms = {"File(s) :", "tempusOrderId :", "accessionId :", "physician :", "signout_date :", "signoutDate :", "tumorCollectionDate :", "diagnosis :", "specimines :", "icdOTxtMorphology :", "icdOTxtTopography :", "icd10Txt :"}; 	
	private ArrayList<String> filterSearchLinesTempus(ArrayList<String> searchLines) {
		ArrayList<String> toReturn = new ArrayList<String>();
		//for each search term
		for (String st: tempusSearchTerms) {
			//for each data line
			for (String dl: searchLines) {
				if (dl.startsWith(st)) {
					String tabbed = dl.replaceFirst(" : ", ":\t");
					toReturn.add(tabbed);
					break;
				}
			}
		}
		return toReturn;
	}
	private ArrayList<String[]> filterSearchLinesTempusSplit(ArrayList<String> searchLines) throws IOException {
		ArrayList<String[]> toReturn = new ArrayList<String[]>();
		//for each search term
		for (String st: tempusSearchTerms) {
			//for each data line
			for (String dl: searchLines) {
				if (dl.startsWith(st)) {
					String[] splitLine = dl.split(" : ",2);
					if (splitLine.length != 2) throw new IOException("Failed to split in two: "+dl);
					toReturn.add(splitLine);
					break;
					
					
				}
			}
		}
		return toReturn;
	}

	private static final String[] avatarSearchTerms = {"File(s) :", "Disease Type :", "Primary/Met :", "Histology/Behavior :", "SpecimenSiteOfCollection :", "SpecimenSiteOfOrigin :", "SpecimenSiteOfOriginRollUp :"}; 	
	private ArrayList<String> filterSearchLinesAvatar(ArrayList<String> searchLines) {
		LinkedHashSet<String> filtered = new LinkedHashSet<String>();
		//for each search term
		for (String st: avatarSearchTerms) {
			//for each data line
			for (String dl: searchLines) {
				//germline?
				if (dl.contains("germline") || dl.contains("Blood")) continue;
				if (dl.startsWith(st)) {
					if (filtered.contains(dl)==false) {
						String tabbed = dl.replaceFirst(" : ", ":\t");
						filtered.add(tabbed);
					}
					break;
				}
			}
		}
		ArrayList<String> toReturn = new ArrayList<String>();
		toReturn.addAll(filtered);
		return toReturn;
	}
	
	private ArrayList<String[]> filterSearchLinesAvatarSplit(ArrayList<String> searchLines) throws IOException {
		LinkedHashMap<String, String[]> filtered = new LinkedHashMap<String, String[]>();
		//for each search term
		for (String st: avatarSearchTerms) {
			//for each data line
			for (String dl: searchLines) {
				//germline?
				if (dl.contains("germline") || dl.contains("Blood")) continue;
				if (dl.startsWith(st)) {
					if (filtered.containsKey(dl)==false) {
						String[] splitLine = dl.split(" : ",2);
						if (splitLine.length != 2) throw new IOException("Failed to split in two: "+dl);
						filtered.put(dl, splitLine);
					}
					break;
				}
			}
		}
		ArrayList<String[]> toReturn = new ArrayList<String[]>();
		toReturn.addAll(filtered.values());
		return toReturn;
	}

	private void parsePMRSearchResultsAll() throws IOException {
		IO.pl("Parsing PMRSearch results...");
		
		// PMR
		//lines with no beginning space
		// zJN3RY4JCn/Avatar/A051837_FT-SA335954_FT-SA336024D_FT-SA336024R
		// xRp2qt5MTT/Tempus/25xtaaop_20250424
		// DWe7RA6vXv/Caris/TN24-287588_20250101
		Pattern noSpace = Pattern.compile("^\\w+.+");
		Pattern pmr = Pattern.compile("(.+)/.+/(.+)");
		Matcher mat = null;
		ArrayList<String> dataLines = null;
		String fullSampleName = null;
		for (String p: pmrDataLines) {
			mat = noSpace.matcher(p);
			
			//a new dataset?
			if (mat.matches()) {
				//save old?
				if (fullSampleName != null)searchIdResults.put(fullSampleName, dataLines);
				dataLines = new ArrayList<String>();
				fullSampleName = p;
				mat = pmr.matcher(p);
				if (mat.matches()) {
					//IO.pl("\t"+mat.group(1)+"\t"+mat.group(2));
					ArrayList<String> al = pmrIdPmrDatasets.get(mat.group(1));
					if (al == null) {
						al = new ArrayList<String>();
						pmrIdPmrDatasets.put(mat.group(1), al);
					}
					al.add(p);
				}
				else throw new IOException("No pmr match for "+p);
			}	
			//nope still in old
			else dataLines.add(p.trim());
		}
		//add last
		searchIdResults.put(fullSampleName, dataLines);
		//IO.pl(pmrIdPmrDatasets);
		/*
		for (String key: searchIdResults.keySet()) {
			IO.pl("Search: "+ key);
			for (String val: searchIdResults.get(key)) {
				IO.pl("\tVal: "+val);
			}
		}*/
	}






	private void loadPMRSearchResults() throws IOException {
		IO.pl("Loading PMRSearch results...");
		String line = null;
		File[] resultsFiles = IO.extractFiles(pmrSearchResultDir, ".txt");
		for (File f: resultsFiles) {
			if (f.getName().startsWith(".")) continue;
			BufferedReader in = IO.fetchBufferedReader(f);
			//find data start line
			boolean startFound = false;
			while ((line = in.readLine()) !=null) {
				if (line.startsWith(pmrResultsStartLine)) {
					startFound = true;
					break;
				}
			}
			if (startFound == false) throw new IOException ("ERROR: Failed to find a line starting with '"+pmrResultsStartLine+"' in  "+f+"\nDid you specify -i when searching?");
			//read in any data lines, might be zero
			while ((line = in.readLine()) !=null) {
				//look for a blank line and end of results
				int size = line.trim().length();
				if (size==0) break;
				pmrDataLines.add(line);
			}
			in.close();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PMRSearchGQueryJoiner(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException {

		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': pmrSearchResultDir = new File(args[++i]); break;
					case 'g': gqueryResultFile = new File(args[++i]); break;
					case 'r': patientRegistryDirectory = new File(args[++i]); break;
					case 'x': priorReportedPmrIds = new File(args[++i]); break;
					case 'o': outputDirectory = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (pmrSearchResultDir == null || pmrSearchResultDir.exists()== false) throw new IOException("\nERROR: failed to find your directory containing PMRSearch -i -n results txt files.");
		if (gqueryResultFile == null || gqueryResultFile.exists()== false) throw new IOException("\nERROR: failed to find your GQuery -s json results file.");
		if (patientRegistryDirectory != null && patientRegistryDirectory.exists()== false) throw new IOException("\nERROR: failed to find the direcory containing a currentRegistry_xxx_PHI.txt file with patient info.");
		if (outputDirectory == null) throw new IOException("\nERROR: faile to find or make your output directory.");
		if (priorReportedPmrIds != null && priorReportedPmrIds.exists()== false) throw new IOException("\nERROR: failed to find your file of PMR IDs to skip, e.g. those already reported.");
		outputDirectory.mkdirs();
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                         PMRSearch GQuery Joiner : Nov 2025                       **\n" +
				"**************************************************************************************\n" +
				"Merges results from the PMRSearch and GQuery tools into a summary spreadsheet.\n"+

				"\nRequired Options:\n"+
				"-p  File or directory containing PMRSearch results xxx.txt with output from the \n"+
				"    -i option.\n"+
				"-g  GQuery json results file.\n"+
				"-o  Output directory for writing the results.\n"+
				"-r  (Optional) Directory containing the currentRegistry_xxx_PHI.txt file for pulling\n"+
				"       patient info.\n"+
				"-x  (Optional) File containing PMR IDs to skip, e.g. prior runs of this tool.\n"+
				
				"\nExample: java -jar pathToUSeq/Apps/PMRSearchGQueryJoiner -p PmrSearchResults -o \n"+
				"   JointResults -g gqueryKrasResults.json -r ~/PHI/Registry/ -x priorPmrs.txt\n"+

				"\n**************************************************************************************\n");
	}
}
