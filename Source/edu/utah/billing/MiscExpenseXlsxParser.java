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

public class MiscExpenseXlsxParser {

	private TreeMap<String, ArrayList<MiscExpense>> groupNameExpense = new TreeMap<String, ArrayList<MiscExpense>>();
	private File inputFile = null;

	public MiscExpenseXlsxParser(File xlsx, boolean debug) {
		inputFile = xlsx;
		parseIt();

		if (debug) {
			for (String s: groupNameExpense.keySet()) {
				ArrayList<MiscExpense> al = groupNameExpense.get(s);
				for (MiscExpense e: al) {
					IO.pl(s+"\t"+e.getCost()+"\t"+e.getDescription());
				}
			}
		}
	}

	public static void main(String[] args) {
		new MiscExpenseXlsxParser(
				new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2024/3-4_CBI_Feb-Mar_2024/OneTimeExpenses/miscComputingExpensesFeb2024.xlsx"), 
				true);

	}

	private void parseIt() {
		try {

			//Open up xlsx file
			Workbook wb = WorkbookFactory.create(inputFile);	

			//Find appropriate sheet 
			Sheet sheet = wb.getSheetAt(0);
			if (sheet == null) throw new IOException("Could not find a sheet in "+inputFile+" ?");

			//Iterate through rows
			int numRows = sheet.getPhysicalNumberOfRows()+1;
			for (int r = 0; r< numRows; r++) {
				Row row = sheet.getRow(r);
				if (row != null) parseExpense(r, row);
			} 
		} catch (Exception e) {
			System.out.println("The expense xlsx file is not in the correct format, exiting");
			e.printStackTrace();
			System.exit(1);
		} 
	}

	private void parseExpense(int rowNumber, Row row) throws IOException {
		int numCells = row.getPhysicalNumberOfCells();

		//skip blanks
		if (numCells == 0) return;
		//error if not three
		if (numCells != 3) throw new IOException("FAILED to find 3 cells in row "+rowNumber+" in "+inputFile);

		Cell groupNameCell = row.getCell(0);
		String groupName = groupNameCell.toString().trim();
		//skip header line
		if (groupName.startsWith("INVESTIGATOR") || groupName.length()==0) return;

		double expense = row.getCell(1).getNumericCellValue();
		String description = row.getCell(2).toString().trim();
		ArrayList<MiscExpense> al = groupNameExpense.get(groupName);
		if (al == null) {
			al = new ArrayList<MiscExpense>();
			groupNameExpense.put(groupName, al);
		}
		al.add(new MiscExpense(expense, description));

	}

	public TreeMap<String, ArrayList<MiscExpense>> getGroupNameExpense() {
		return groupNameExpense;
	}
}
