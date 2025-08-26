package edu.utah.seq.run.ambry;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import util.gen.IO;


public class AmbryXlsxParser {

	//fields
	private HashMap<String, AmbryPatient> patientKeyResults = new HashMap<String, AmbryPatient>();
	private HashSet<String> transcriptCDots = new HashSet<String>();
	
	//header info
	private LinkedHashMap<String, Integer> headerIndexName = null;
	
	public static final String mrn = "mrn";
	public static final String patientFirstName = "patientfirstname";
	public static final String patientLastName = "patientlastname";
	public static final String dateOfBirth = "dateofbirth";
	public static final String gender = "gender";
	public static final String specimenCollectionDate = "specimencollectiondate";
	public static final String familialVariantPreviouslyIdentified = "familialvariantpreviouslyidentified";
	
	public static final String orderType = "ordertype";
	public static final String accessionNumber = "accessionnumber";
	public static final String orderNumber = "ordernumber";
	public static final String reportSignedDate = "reportsigneddate";
	public static final String genesTested = "genestested";
	public static final String geneticTestPanelOrdered = "genetictestpanelordered";
	public static final String reportComment = "reportcomment";
	public static final String organizationName = "organizationname";
	public static final String physician = "physician";
	public static final String additionalRecipient = "additionalrecipient";
	
	public static final String geneResult1 = "generesult1";
	public static final String nucleotide_id1 = "nucleotide_id1";
	public static final String zygostiy1 = "zygostiy1";
	public static final String c_variant1 = "c_variant1";
	public static final String p_variant1 = "p_variant1";
	public static final String method1 = "method1";
	public static final String categoryDescription1 = "categorydescription1";
	public static final String chromosome1 = "chromosome1";
	public static final String genomicStart1 = "genomicstart1";
	public static final String ref1 = "ref1";
	public static final String alt1 = "alt1";
	
	public static final String geneResult2 = "generesult2";
	public static final String nucleotide_id2 = "nucleotide_id2";
	public static final String zygostiy2 = "zygostiy2";
	public static final String c_variant2 = "c_variant2";
	public static final String p_variant2 = "p_variant2";
	public static final String method2 = "method2";
	public static final String categoryDescription2 = "categorydescription2";
	public static final String chromosome2 = "chromosome2";
	public static final String genomicStart2 = "genomicstart2";
	public static final String ref2 = "ref2";
	public static final String alt2 = "alt2";
	
	public static final String geneResult3 = "generesult3";
	public static final String nucleotide_id3 = "nucleotide_id3";
	public static final String zygostiy3 = "zygostiy3";
	public static final String c_variant3 = "c_variant3";
	public static final String p_variant3 = "p_variant3";
	public static final String method3 = "method3";
	public static final String categoryDescription3 = "categorydescription3";
	public static final String chromosome3 = "chromosome3";
	public static final String genomicStart3 = "genomicstart3";
	public static final String ref3 = "ref3";
	public static final String alt3 = "alt3";
	
	public static final String geneResult4 = "generesult4";
	public static final String nucleotide_id4 = "nucleotide_id4";
	public static final String zygostiy4 = "zygostiy4";
	public static final String c_variant4 = "c_variant4";
	public static final String p_variant4 = "p_variant4";
	public static final String method4 = "method4";
	public static final String categoryDescription4 = "categorydescription4";
	public static final String chromosome4 = "chromosome4";
	public static final String genomicStart4 = "genomicstart4";
	public static final String ref4 = "ref4";
	public static final String alt4 = "alt4";
	
	public int mrnIndex = -1;
	public int patientFirstNameIndex = -1;
	public int patientLastNameIndex = -1;
	public int dateOfBirthIndex = -1;
	public int genderIndex = -1;
	public int specimenCollectionDateIndex = -1;
	public int familialVariantPreviouslyIdentifiedIndex = -1;
	public int orderTypeIndex = -1;
	public int accessionNumberIndex = -1;
	public int orderNumberIndex = -1;
	public int reportSignedDateIndex = -1;
	public int genesTestedIndex = -1;
	public int geneticTestPanelOrderedIndex = -1;
	public int reportCommentIndex = -1;
	public int organizationNameIndex = -1;
	public int physicianIndex = -1;
	public int additionalRecipientIndex = -1;
	public int geneResult1Index = -1;
	public int nucleotide_id1Index = -1;
	public int zygostiy1Index = -1;
	public int c_variant1Index = -1;
	public int p_variant1Index = -1;
	public int method1Index = -1;
	public int categoryDescription1Index = -1;
	public int chromosome1Index = -1;
	public int genomicStart1Index = -1;
	public int ref1Index = -1;
	public int alt1Index = -1;
	public int geneResult2Index = -1;
	public int nucleotide_id2Index = -1;
	public int zygostiy2Index = -1;
	public int c_variant2Index = -1;
	public int p_variant2Index = -1;
	public int method2Index = -1;
	public int categoryDescription2Index = -1;
	public int chromosome2Index = -1;
	public int genomicStart2Index = -1;
	public int ref2Index = -1;
	public int alt2Index = -1;
	public int geneResult3Index = -1;
	public int nucleotide_id3Index = -1;
	public int zygostiy3Index = -1;
	public int c_variant3Index = -1;
	public int p_variant3Index = -1;
	public int method3Index = -1;
	public int categoryDescription3Index = -1;
	public int chromosome3Index = -1;
	public int genomicStart3Index = -1;
	public int ref3Index = -1;
	public int alt3Index = -1;
	public int geneResult4Index = -1;
	public int nucleotide_id4Index = -1;
	public int zygostiy4Index = -1;
	public int c_variant4Index = -1;
	public int p_variant4Index = -1;
	public int method4Index = -1;
	public int categoryDescription4Index = -1;
	public int chromosome4Index = -1;
	public int genomicStart4Index = -1;
	public int ref4Index = -1;
	public int alt4Index = -1;


	public AmbryXlsxParser(File xlsx, boolean debug) throws Exception {
		if (debug) IO.pl("\nParsing\t'"+xlsx.getName()+"'");
		parseIt(xlsx);
		
		if (debug) {
			int numResults = 0;
			Iterator<AmbryPatient> it = patientKeyResults.values().iterator();
			while (it.hasNext()) {
				AmbryPatient p = it.next(); 
				ArrayList<AmbryResult> results = p.getAmbryResults();
				numResults += results.size();
				//if (results.size()>1) IO.pl(p.getPatientKey()+"\t"+p.getAmbryResults().size());
			}
			IO.pl("\tNum patients with parsible results: "+patientKeyResults.size());
			IO.pl("\tNum results: "+numResults);
		}
		
		collectCDots();
		
		
		
		//split patient tests by AccessionNumber
			//sort by panel
				//sort by date keeping latest (this will throw out originals and older reissues)
			//merge results for everything under the Accession Number
			//keep latest date
			//only keep if they have coordinate sortable results
				//make AN_RSD folder
				//make json, vcf, bed: AN_RSD-year-date-month_PhysicianFN_LN.xxx; vcf might be empty
		
		
		
		//look at invitae parsing?
				//invitaeBulk.README.sh
					// InvitaeBulkCsvParser
		
				//Invitae/RQ581802_AR558582/ClinicalReport/RQ581802_AR558582_2018-25-12_Sarah_Colonna.json
				//Invitae/RQ588978_AR559986/ClinicalReport/RQ588978_AR559986_2018-28-12_Sarah_Colonna.json
					// vcf and bed for each
	
		
		//merge results?
			//yes, use the PmrVcfBedMerger see the /scratch/general/pe-nfs1/u0028003/Invitae/invitaeAutoProcessing.sh
		
		//Ambry
		// AccessionNumber and OrderNumber often mirror each other but not always
		// AN appears to be the whole order, where as ON refers to ?
		// Sometimes two panels are run with the same AN and ON and ReportSignedDate
		// Reclassified reports have the same AN and ON but diff RSD
		
		// How tell when different?  Same AN, Diff RSD dates
		
		// AN_RSD_PI.json - remove any originals that have been reissued
		
		//Physician, split on , space
		//Last, First, degrees
		// Colonna, Sarah, MD
		// Sanchez, Alejandro, MD
		// Maese, Luke, MD, DO
		
		
	}
	
	private void collectCDots() {
		Iterator<AmbryPatient> it = patientKeyResults.values().iterator();
		while (it.hasNext()) {
			AmbryPatient p = it.next(); 
			ArrayList<AmbryResult> results = p.getAmbryResults();
			for (AmbryResult ar: results) {
				for (AmbryGeneResult g: ar.getGeneResultsWithParsibleResults()) {
					if (g.iscDotPass()) transcriptCDots.add(g.getTranscriptCDot());
				}
			}
		}
	}
	
	private void parseIt(File inputFile) throws Exception {
		//Open up xlsx file
		Workbook wb = WorkbookFactory.create(inputFile);	

		//Find appropriate sheet 
		Sheet sheet = wb.getSheet("Germline");
		if (sheet == null) throw new IOException("Could not find sheet 'Germline' in "+inputFile+" ?");

		//Iterate through rows
		int numRows = sheet.getPhysicalNumberOfRows();
		for (int r = 0; r< numRows; r++) {
			Row row = sheet.getRow(r);
			if (headerIndexName == null) parseHeader(row);
			else parseTestResult(row);
		} 
	}
	
	private void parseHeader(Row row) throws IOException {
		headerIndexName = new LinkedHashMap<String, Integer>();
		ArrayList<String> duplicateHeaders = new ArrayList<String>();
		int numCells = row.getLastCellNum()+1;
		for (int c=0;c < numCells; c++) {
			Cell cell = row.getCell(c);
			if (cell != null) {
				String value = cell.toString().trim().toLowerCase().replaceAll(" ", "");
				if (headerIndexName.containsKey(value)) duplicateHeaders.add(cell.toString());
				else headerIndexName.put(value, c);
			}
		}
		
		if (duplicateHeaders.size()!=0) throw new IOException("Duplicate header value(s)! Correct and restart "+duplicateHeaders);
		//IO.pl("Header: "+headerIndexName);
		
		setColumnIndexes();
	
	}

	private void setColumnIndexes() {
		mrnIndex = headerIndexName.get( mrn );
		patientFirstNameIndex = headerIndexName.get( patientFirstName );
		patientLastNameIndex = headerIndexName.get( patientLastName );
		dateOfBirthIndex = headerIndexName.get( dateOfBirth );
		genderIndex = headerIndexName.get( gender );
		specimenCollectionDateIndex = headerIndexName.get( specimenCollectionDate );
		familialVariantPreviouslyIdentifiedIndex = headerIndexName.get( familialVariantPreviouslyIdentified );
		orderTypeIndex = headerIndexName.get( orderType );
		accessionNumberIndex = headerIndexName.get( accessionNumber );
		orderNumberIndex = headerIndexName.get( orderNumber );
		reportSignedDateIndex = headerIndexName.get( reportSignedDate );
		genesTestedIndex = headerIndexName.get( genesTested );
		geneticTestPanelOrderedIndex = headerIndexName.get( geneticTestPanelOrdered );
		reportCommentIndex = headerIndexName.get( reportComment );
		organizationNameIndex = headerIndexName.get( organizationName );
		physicianIndex = headerIndexName.get( physician );
		additionalRecipientIndex = headerIndexName.get( additionalRecipient );
		geneResult1Index = headerIndexName.get( geneResult1 );
		nucleotide_id1Index = headerIndexName.get( nucleotide_id1 );
		zygostiy1Index = headerIndexName.get( zygostiy1 );
		c_variant1Index = headerIndexName.get( c_variant1 );
		p_variant1Index = headerIndexName.get( p_variant1 );
		method1Index = headerIndexName.get( method1 );
		categoryDescription1Index = headerIndexName.get( categoryDescription1 );
		if (headerIndexName.containsKey(chromosome1)) chromosome1Index = headerIndexName.get( chromosome1 );
		if (headerIndexName.containsKey(genomicStart1)) genomicStart1Index = headerIndexName.get( genomicStart1 );
		if (headerIndexName.containsKey(ref1)) ref1Index = headerIndexName.get( ref1 );
		if (headerIndexName.containsKey(alt1)) alt1Index = headerIndexName.get( alt1 );
		geneResult2Index = headerIndexName.get( geneResult2 );
		nucleotide_id2Index = headerIndexName.get( nucleotide_id2 );
		zygostiy2Index = headerIndexName.get( zygostiy2 );
		c_variant2Index = headerIndexName.get( c_variant2 );
		p_variant2Index = headerIndexName.get( p_variant2 );
		method2Index = headerIndexName.get( method2 );
		categoryDescription2Index = headerIndexName.get( categoryDescription2 );
		if (headerIndexName.containsKey(chromosome2)) chromosome2Index = headerIndexName.get( chromosome2 );
		if (headerIndexName.containsKey(genomicStart2)) genomicStart2Index = headerIndexName.get( genomicStart2 );
		if (headerIndexName.containsKey(ref2)) ref2Index = headerIndexName.get( ref2 );
		if (headerIndexName.containsKey(alt2)) alt2Index = headerIndexName.get( alt2 );
		geneResult3Index = headerIndexName.get( geneResult3 );
		nucleotide_id3Index = headerIndexName.get( nucleotide_id3 );
		zygostiy3Index = headerIndexName.get( zygostiy3 );
		c_variant3Index = headerIndexName.get( c_variant3 );
		p_variant3Index = headerIndexName.get( p_variant3 );
		method3Index = headerIndexName.get( method3 );
		categoryDescription3Index = headerIndexName.get( categoryDescription3 );
		if (headerIndexName.containsKey(chromosome3)) chromosome3Index = headerIndexName.get( chromosome3 );
		if (headerIndexName.containsKey(genomicStart3)) genomicStart3Index = headerIndexName.get( genomicStart3 );
		if (headerIndexName.containsKey(ref3)) ref3Index = headerIndexName.get( ref3 );
		if (headerIndexName.containsKey(alt3)) alt3Index = headerIndexName.get( alt3 );
		geneResult4Index = headerIndexName.get( geneResult4 );
		nucleotide_id4Index = headerIndexName.get( nucleotide_id4 );
		zygostiy4Index = headerIndexName.get( zygostiy4 );
		c_variant4Index = headerIndexName.get( c_variant4 );
		p_variant4Index = headerIndexName.get( p_variant4 );
		method4Index = headerIndexName.get( method4 );
		categoryDescription4Index = headerIndexName.get( categoryDescription4 );
		if (headerIndexName.containsKey(chromosome4)) chromosome4Index = headerIndexName.get( chromosome4 );
		if (headerIndexName.containsKey(genomicStart4)) genomicStart4Index = headerIndexName.get( genomicStart4 );
		if (headerIndexName.containsKey(ref4)) ref4Index = headerIndexName.get( ref4 );
		if (headerIndexName.containsKey(alt4)) alt4Index = headerIndexName.get( alt4 );
	}

	private void parseTestResult(Row row) throws IOException {		
		//Create patient key and fetch or make the patient
		AmbryResult result = new AmbryResult(row, this);
		if (result.getPatientKey() == null) return;
		
		//should it be parsed?
		if (result.isTestedButNoGeneResults() || result.getGeneResultsWithParsibleResults().size()!=0) {
			String patientKey = result.getPatientKey();
			AmbryPatient p = patientKeyResults.get(patientKey);
			if (p==null) {
				p = new AmbryPatient(patientKey);
				patientKeyResults.put(patientKey, p);
			}
			//add in result
			p.getAmbryResults().add(result);
		}
	}

	public static void main(String[] args) throws Exception {
		new AmbryXlsxParser(new File ("/Users/u0028003/HCI/ClinicalGenomics/Ambry/WithVCFInfo/AmbryImportFiles/HCI_OrdersResults01.01.25-03.31.25.xlsx"), true);
	}

	public HashSet<String> getTranscriptCDots() {
		return transcriptCDots;
	}

	public HashMap<String, AmbryPatient> getPatientKeyResults() {
		return patientKeyResults;
	}

}
