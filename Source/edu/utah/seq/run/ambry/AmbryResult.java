package edu.utah.seq.run.ambry;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;

import util.gen.Misc;

public class AmbryResult {
	
	private Row xlsxRow = null;
	private AmbryXlsxParser parser = null;

	private String patientKey;
	
	//populated with construction
	private String patientFirstName;	//'Kathryn' or with nickname 'Lefujio "Edrian"'
	private String patientLastName;		// Smith
	private String dateOfBirth;			// 18-Nov-1945
	private String gender;				// M or F
	private String mrn;					// 2.1999697E7, bad parsing from excel
	private String specimenCollectionDate;	// 18-Nov-1945
	private String familialVariantPreviouslyIdentified;	// empty?
	
	private String orderType;			// Panel Test or Specific Site Test, contains spaces
	private String accessionNumber;		// 24-654072
	private String orderNumber;			// 2899643
	private String reportSignedDate;	// 18-Nov-1945
	private String reportComment;		// A copy of the clinical report documenting the familial alteration (designated as 3450del4) was provided for review; however, a positive control was not available for concurrent testing. 
	
	private String genesTested;			// 'CHEK2: c.428A>G (p.H143R)' or 'AIP,ALK,APC,ATM,ATRIP,AXIN2,BAP1,BARD1,BMPR1A,BRCA1,BRCA2,BRIP1,CDC73,CDH1,CDK4,CDKN1B,CDKN2A,CEBPA,CFTR,CHEK2,CPA1,CTNNA1,CTRC,DDX41,DICER1,EGFR,EGLN1,EPCAM,ETV6,FH,FLCN,GATA2,GREM1,HOXB13,KIF1B,KIT,LZTR1,MAX,MBD4,MEN1,MET,MITF,MLH1,MLH3,MSH2,MSH3,MSH6,MUTYH,NF1,NF2,NTHL1,PALB2,PALLD,PDGFRA,PHOX2B,PMS2,POLD1,POLE,POT1,PRKAR1A,PRSS1,PTCH1,PTEN,RAD51B,RAD51C,RAD51D,RB1,RET,RNF43,RPS20,RUNX1,SDHA,SDHAF2,SDHB,SDHC,SDHD,SMAD4,SMARCA4,SMARCB1,SMARCE1,SPINK1,STK11,SUFU,TERT,TMEM127,TP53,TSC1,TSC2,VHL,WT1'
	private String geneticTestPanelOrdered;	// CancerNext-Expanded +RNAinsight or CHEK2 specific site analysis

	private String organizationName;	// Huntsman Cancer Institute (06074)
	private String physician;			// Colonna, Sarah, MD
	private String additionalRecipient;	// Counseling Assistant, Genetic BSPauley, Kristen MS, CGC
	
	private ArrayList<AmbryGeneResult> geneResultsWithParsibleResults = new ArrayList<AmbryGeneResult>();
	private boolean testedButNoGeneResults = false;
	
	public AmbryResult (Row xlsxRow, AmbryXlsxParser parser) throws IOException {
		this.xlsxRow = xlsxRow;
		this.parser = parser;
		if (parsePatientKey(xlsxRow)) setResults();
	}
	

	public String getSubjectMatchMakerInfo() throws IOException, ParseException {
		//  lastName firstName dobMonth(1-12) dobDay(1-31) dobYear(1900-2050) gender(M|F) mrn
		StringBuilder sb = new StringBuilder();
		sb.append(patientLastName);	sb.append("\t");
		
		String[] trimmedFN = Misc.FORWARD_PARENTHESIS.split(patientFirstName);
		sb.append(trimmedFN[0].trim()); sb.append("\t");
		
		Date date = AmbryDataWrangler.dateFormatter.parse(dateOfBirth);
		Calendar calendar = Calendar.getInstance();
		calendar.setTime(date);
		
		sb.append(1+calendar.get((Calendar.MONTH))); sb.append("\t");
		sb.append(calendar.get(Calendar.DAY_OF_MONTH)); sb.append("\t");
		sb.append(calendar.get(Calendar.YEAR)); sb.append("\t");
		
		sb.append(gender);
		
		if (mrn != null) {
			Double dObject = Double.parseDouble(mrn);
			sb.append("\t");
			sb.append(dObject.intValue());	
		}
		
		return sb.toString();
	}
	

	private boolean parsePatientKey(Row row) throws IOException {
		//some results there is no mrn so create a key based on first_last_dob_gender_mrn, these seem to be consistent with all Ambry reports
		mrn = getCellString(row, parser.mrnIndex);
		patientFirstName = getCellString(row, parser.patientFirstNameIndex);
		patientLastName = getCellString(row, parser.patientLastNameIndex);
		dateOfBirth = getCellString(row, parser.dateOfBirthIndex);
		gender = getCellString(row, parser.genderIndex);
		
		if (patientFirstName == null || patientLastName == null || dateOfBirth == null || gender == null) return false;
		
		if (gender.equals("M")==false && gender.equals("F")==false) throw new IOException("Error: failed to parse M or F from gender, see accession "+accessionNumber);
		
		StringBuilder keySb = new StringBuilder();
		keySb.append(patientFirstName); keySb.append("_");
		keySb.append(patientLastName); keySb.append("_");
		keySb.append(dateOfBirth); keySb.append("_");
		keySb.append(gender); keySb.append("_");
		keySb.append(mrn);
		patientKey = keySb.toString();
		return true;
	}
	
	public static String getCellString(Row row, int index) {
		if (index < 0) return null;
		Cell cell = row.getCell(index);
		if (cell == null) return null;
		String s = cell.toString().trim();
		if (s.length()==0) return null;
		return s;
	}
	
	public void setResults() {
		specimenCollectionDate = getCellString(xlsxRow, parser.specimenCollectionDateIndex );
		familialVariantPreviouslyIdentified = getCellString(xlsxRow, parser.familialVariantPreviouslyIdentifiedIndex );
		orderType = getCellString(xlsxRow, parser.orderTypeIndex );
		accessionNumber = getCellString(xlsxRow, parser.accessionNumberIndex );
		orderNumber = getCellString(xlsxRow, parser.orderNumberIndex );
		reportSignedDate = getCellString(xlsxRow, parser.reportSignedDateIndex );
		genesTested = getCellString(xlsxRow, parser.genesTestedIndex );
		geneticTestPanelOrdered = getCellString(xlsxRow, parser.geneticTestPanelOrderedIndex );
		reportComment = getCellString(xlsxRow, parser.reportCommentIndex );
		organizationName = getCellString(xlsxRow, parser.organizationNameIndex );
		physician = getCellString(xlsxRow, parser.physicianIndex );
		additionalRecipient = getCellString(xlsxRow, parser.additionalRecipientIndex );
		
		String geneResult1 = getCellString(xlsxRow, parser.geneResult1Index );
		//if this is null then nothing will follow
		if (geneResult1 == null) {
			testedButNoGeneResults = true;
			return;
		}
		else testedButNoGeneResults = false;
		
		AmbryGeneResult agr = new AmbryGeneResult(
				this,
				geneResult1,
				getCellString(xlsxRow, parser.nucleotide_id1Index),
				getCellString(xlsxRow, parser.zygostiy1Index),
				getCellString(xlsxRow, parser.c_variant1Index),
				getCellString(xlsxRow, parser.p_variant1Index),
				getCellString(xlsxRow, parser.method1Index),
				getCellString(xlsxRow, parser.categoryDescription1Index),
				getCellString(xlsxRow, parser.chromosome1Index),
				getCellString(xlsxRow, parser.genomicStart1Index),
				getCellString(xlsxRow, parser.ref1Index),
				getCellString(xlsxRow, parser.alt1Index));
		if (agr.iscDotPass() || agr.isChromPass()) geneResultsWithParsibleResults.add(agr);
		
		String geneResult2 = getCellString(xlsxRow, parser.geneResult2Index );
		if (geneResult2 == null) return;
		agr = new AmbryGeneResult(
				this,
				geneResult2,
				getCellString(xlsxRow, parser.nucleotide_id2Index),
				getCellString(xlsxRow, parser.zygostiy2Index),
				getCellString(xlsxRow, parser.c_variant2Index),
				getCellString(xlsxRow, parser.p_variant2Index),
				getCellString(xlsxRow, parser.method2Index),
				getCellString(xlsxRow, parser.categoryDescription2Index),
				getCellString(xlsxRow, parser.chromosome2Index),
				getCellString(xlsxRow, parser.genomicStart2Index),
				getCellString(xlsxRow, parser.ref2Index),
				getCellString(xlsxRow, parser.alt2Index));
		if (agr.iscDotPass() || agr.isChromPass()) geneResultsWithParsibleResults.add(agr);
		
		String geneResult3 = getCellString(xlsxRow, parser.geneResult3Index );
		if (geneResult3 == null) return;
		agr = new AmbryGeneResult(
				this,
				geneResult3,
				getCellString(xlsxRow, parser.nucleotide_id3Index),
				getCellString(xlsxRow, parser.zygostiy3Index),
				getCellString(xlsxRow, parser.c_variant3Index),
				getCellString(xlsxRow, parser.p_variant3Index),
				getCellString(xlsxRow, parser.method3Index),
				getCellString(xlsxRow, parser.categoryDescription3Index),
				getCellString(xlsxRow, parser.chromosome3Index),
				getCellString(xlsxRow, parser.genomicStart3Index),
				getCellString(xlsxRow, parser.ref3Index),
				getCellString(xlsxRow, parser.alt3Index));
		if (agr.iscDotPass() || agr.isChromPass()) geneResultsWithParsibleResults.add(agr);
		
		String geneResult4 = getCellString(xlsxRow, parser.geneResult4Index );
		if (geneResult4 == null) return;
		agr = new AmbryGeneResult(
				this,
				geneResult4,
				getCellString(xlsxRow, parser.nucleotide_id4Index),
				getCellString(xlsxRow, parser.zygostiy4Index),
				getCellString(xlsxRow, parser.c_variant4Index),
				getCellString(xlsxRow, parser.p_variant4Index),
				getCellString(xlsxRow, parser.method4Index),
				getCellString(xlsxRow, parser.categoryDescription4Index),
				getCellString(xlsxRow, parser.chromosome4Index),
				getCellString(xlsxRow, parser.genomicStart4Index),
				getCellString(xlsxRow, parser.ref4Index),
				getCellString(xlsxRow, parser.alt4Index));
		if (agr.iscDotPass() || agr.isChromPass()) geneResultsWithParsibleResults.add(agr);
	}

	public String getPatientKey() {
		return patientKey;
	}

	public Row getXlsxRow() {
		return xlsxRow;
	}

	public AmbryXlsxParser getParser() {
		return parser;
	}

	public String getPatientFirstName() {
		return patientFirstName;
	}

	public String getPatientLastName() {
		return patientLastName;
	}

	public String getDateOfBirth() {
		return dateOfBirth;
	}

	public String getGender() {
		return gender;
	}

	public String getMrn() {
		return mrn;
	}

	public String getSpecimenCollectionDate() {
		return specimenCollectionDate;
	}

	public String getFamilialVariantPreviouslyIdentified() {
		return familialVariantPreviouslyIdentified;
	}

	public String getOrderType() {
		return orderType;
	}

	public String getAccessionNumber() {
		return accessionNumber;
	}

	public String getOrderNumber() {
		return orderNumber;
	}

	public String getReportSignedDate() {
		return reportSignedDate;
	}

	public String getReportComment() {
		return reportComment;
	}

	public String getGenesTested() {
		return genesTested;
	}

	public String getGeneticTestPanelOrdered() {
		return geneticTestPanelOrdered;
	}

	public String getOrganizationName() {
		return organizationName;
	}

	public String getPhysician() {
		return physician;
	}

	public String getAdditionalRecipient() {
		return additionalRecipient;
	}


	public boolean isTestedButNoGeneResults() {
		return testedButNoGeneResults;
	}


	public ArrayList<AmbryGeneResult> getGeneResultsWithParsibleResults() {
		return geneResultsWithParsibleResults;
	}
	
}
