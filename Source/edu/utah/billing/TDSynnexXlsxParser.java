package edu.utah.billing;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.TreeMap;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import edu.utah.hci.bioinfo.smm.Util;
import util.gen.IO;
import util.gen.Misc;

public class TDSynnexXlsxParser {

	private TreeMap<String, ArrayList<Float>> accountExpenses = new TreeMap<String, ArrayList<Float>>();
	private TreeMap<String, Float> awsAccountNumberTotals = new TreeMap<String, Float>();

	//column names and their found indexes
	private String accountName = "Account";
	private int accountIndex = -1;
	private String expenseName1 = "`"; //really weird!
	private String expenseName2 = "Price"; //really weird!
	private int expenseIndex = -1;
	private int numParsedLines = 0;

	public TDSynnexXlsxParser(File dir, boolean debug) throws IOException {

		//pull all of the xlsx files, there will be several months worth
		File[] xlsxFiles = IO.extractFiles(dir, ".xlsx");
		if (xlsxFiles == null || xlsxFiles.length ==0) throw new IOException("ERROR: failed to find xlsx TDSynnex files in "+dir);

		//parse each
		for (File xlsx: xlsxFiles) parseIt(xlsx);

		//sum all of the expenses for each account
		float awsTotal = 0.0f;
		for (String s: accountExpenses.keySet()) {
			float total = 0;
			for (Float f: accountExpenses.get(s)) total += f;
			//trim name
			if (total >= 0.01f) awsAccountNumberTotals.put(s, total);
			awsTotal += total;
		}

		//print results
		if (debug) {
			IO.pl("AWS account and total expenses:");
			for (String s: awsAccountNumberTotals.keySet()) {
				IO.pl(s+"\t"+awsAccountNumberTotals.get(s));
			}
			IO.pl("\nTotal AWS: "+awsTotal);
		}
	}

	public static void main(String[] args) throws IOException {
		TDSynnexXlsxParser p = new TDSynnexXlsxParser(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2024/3-4_CBI_Feb-Mar_2024/TDSynnex/"), true);

	}

	private void parseIt(File inputFile) {
		//watchout for working and partial saved excel files
		if (inputFile.getName().startsWith("~$")) return;

		//reset
		accountIndex = -1;
		expenseIndex = -1;
		numParsedLines = 0;


		try {
			//Open up xlsx file
			Workbook wb = WorkbookFactory.create(inputFile);	

			//Find appropriate sheet 
			Sheet sheet = wb.getSheet("Detail");
			if (sheet == null) {
				sheet = wb.getSheet("Sheet1");
				if (sheet == null) throw new IOException("Could not find the 'Detail' or 'Sheet1' sheet in "+inputFile+" ?");
			}

			//Iterate through rows
			int numRows = sheet.getPhysicalNumberOfRows();
			for (int r = 0; r< numRows; r++) {
				Row row = sheet.getRow(r);
				if (row != null) parseRow(row);
			} 

		} catch (Exception e) {
			Util.el("\nERROR parsing "+inputFile);
			e.printStackTrace();
			System.exit(1);
		} 

		//check that data lines were parsed
		if (numParsedLines == 0) Misc.printErrAndExit("\nFailed to parse any data lines from "+inputFile);


	}

	private void parseRow(Row row) {

		//convert to String[]
		int numCells = row.getPhysicalNumberOfCells()+1;
		String[] cellStrings = new String[numCells];
		for (int c=0;c < numCells; c++) {
			Cell cell = row.getCell(c);
			if (cell != null) cellStrings[c] = cell.toString().trim(); 
		}

		//header row
		if (accountIndex == -1) {
			if (cellStrings[0].contains("Billing")) {
				for (int i=1; i< numCells; i++) {
					if (cellStrings[i]!= null) {
						if (cellStrings[i].equals(accountName)) accountIndex = i;
						else if (cellStrings[i].equals(expenseName2) || cellStrings[i].equals(expenseName1)) expenseIndex = i;
					}
				}
			}
			//check that both were found
			if (accountIndex == -1 || expenseIndex == -1) Misc.printErrAndExit("ERROR: failed to find the account or expense indexes in the TDSynnex xlsx sheet.");
		}

		// data line
		else {
			numParsedLines++;
			//watchout for E formatting in excel sheet for the account #
			Double lAccName = Double.parseDouble(cellStrings[accountIndex]);
			String corrAccountName = BigDecimal.valueOf(lAccName).toPlainString();
			ArrayList<Float> expenses = accountExpenses.get(corrAccountName);
			if (expenses == null) {
				expenses = new ArrayList<Float>();
				accountExpenses.put(corrAccountName, expenses);
			}
			Float cost = Float.parseFloat(cellStrings[expenseIndex]);		
			if (cost > 0.0f) expenses.add(cost);
		}
	}

	public TreeMap<String, Float> getAwsAccountNumberTotals() {
		return awsAccountNumberTotals;
	}

}
