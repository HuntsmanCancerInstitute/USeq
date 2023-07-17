package edu.utah.billing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;

import util.gen.IO;
import util.gen.Misc;

public class CarahsoftXlsxParser {

	private TreeMap<String, ArrayList<Float>> accountExpenses = new TreeMap<String, ArrayList<Float>>();
	private TreeMap<String, Float> awsAccountNumberTotals = new TreeMap<String, Float>();

	private int finalCostIndex = -1;
	private int numParsedLines = 0;

	public CarahsoftXlsxParser(File dir, boolean debug) {
		File[] xlsxFiles = IO.extractFiles(dir, ".xlsx");

		for (File xlsx: xlsxFiles) parseIt(xlsx);

		float awsTotal = 0.0f;
		for (String s: accountExpenses.keySet()) {
			float total = 0;
			for (Float f: accountExpenses.get(s)) total += f;
			//trim name
			
			if (total >= 0.01f) {
				Long acc = Long.parseLong(s);
				awsAccountNumberTotals.put(acc.toString(), total);
			}
			awsTotal += total;
		}

		if (debug) {
			IO.pl("AWS Account and expenses:");
			for (String s: awsAccountNumberTotals.keySet()) {
				IO.pl(s+"\t"+awsAccountNumberTotals.get(s));
			}
			IO.pl("\nTotal AWS: "+awsTotal);
		}
	}

	public static void main(String[] args) {
		CarahsoftXlsxParser p = new CarahsoftXlsxParser(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2023/6_BSR_June_2023/CarahsoftMay2023"), true);

	}

	private void parseIt(File inputFile) {
		//watchout for working and partial saved excel files
		if (inputFile.getName().startsWith("~$")) return;

		//reset
		finalCostIndex = -1;
		numParsedLines = 0;

		try {
			//Open up xlsx file
			Workbook wb = WorkbookFactory.create(inputFile);	

			//Find appropriate sheet 

			Sheet sheet = wb.getSheetAt(0);
			if (sheet == null) throw new IOException("Could not find sheet in "+inputFile+" ?");

			//Iterate through rows
			int numRows = sheet.getPhysicalNumberOfRows();
			for (int r = 0; r< numRows; r++) {
				Row row = sheet.getRow(r);
				if (row != null) parseRow(row);
			} 

		} catch (Exception e) {
			//keep this silent since all of the Carahsoft xlsx files are broken!
			e.printStackTrace();
		} 

		//check that data lines were parsed
		if (numParsedLines == 0) Misc.printErrAndExit("\nFailed to parse any data lines from "+inputFile+" .\nDid you open, allow Excel to fix, then save overwritting the origin? These come broken from Carahsoft.");
	}

	private void parseRow(Row row) {

		//convert to String[]
		int numCells = row.getLastCellNum()+1;
		String[] cellStrings = new String[numCells];
		for (int c=0;c < numCells; c++) {
			Cell cell = row.getCell(c);
			if (cell != null) cellStrings[c] = cell.toString(); 
		}

		//header row
		if (finalCostIndex == -1) {
			if (cellStrings[0].equals("Linked Account")) {
				for (int i=1; i< numCells; i++) {
					if (cellStrings[i]!= null && cellStrings[i].equals("Cost")) finalCostIndex = i;
				}
			}
		}
		// data line
		else {
			numParsedLines++;
			ArrayList<Float> expenses = accountExpenses.get(cellStrings[0]);
			if (expenses == null) {
				expenses = new ArrayList<Float>();
				accountExpenses.put(cellStrings[0], expenses);
			}
			Float cost = Float.parseFloat(cellStrings[finalCostIndex]);
			if (cost > 0.0f) expenses.add(cost);
		}


	}

	public TreeMap<String, Float> getAwsAccountNumberTotals() {
		return awsAccountNumberTotals;
	}

}
