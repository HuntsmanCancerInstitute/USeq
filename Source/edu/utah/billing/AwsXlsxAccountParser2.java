package edu.utah.billing;

import java.io.File;
import java.io.IOException;
import java.util.TreeMap;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import util.gen.IO;


public class AwsXlsxAccountParser2 {

	private TreeMap<String, String> awsAccountGroupName = new TreeMap<String, String>();

	public AwsXlsxAccountParser2(File xlsx, boolean debug) {
		parseIt(xlsx);
		
		if (debug) {
			for (String s: awsAccountGroupName.keySet()) {
				IO.pl(s+"\t"+awsAccountGroupName.get(s));
			}
		}
	}

	public static void main(String[] args) {
		//AwsXlsxAccountParser p = new AwsXlsxAccountParser(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2023/6_BSR_June_2023/awsAccounts.xlsx"), true);

	}

	private void parseIt(File inputFile) {
		try {

			//Open up xlsx file
			Workbook wb = WorkbookFactory.create(inputFile);	

			//Find appropriate sheet 
			Sheet sheet = wb.getSheetAt(0);
			if (sheet == null) throw new IOException("Could not find a sheet in "+inputFile+" ?");

			//Iterate through rows
			int numRows = sheet.getPhysicalNumberOfRows();
			for (int r = 0; r< numRows; r++) {
				Row row = sheet.getRow(r);
				if (row != null) addAccount(row);
			} 
		} catch (Exception e) {
			System.out.println("Aws Accounts xlsx file is not in the correct format, exiting");
			e.printStackTrace();
			System.exit(1);
		} 
	}

	private void addAccount(Row row) {
		int numCells = row.getLastCellNum()+1;
		if (numCells < 2) return;
		
		Cell groupNameCell = row.getCell(0);
		if (groupNameCell == null) return;
		String groupName = groupNameCell.toString().trim();
		if (groupName.startsWith("INVESTIGATOR") || groupName.length()==0) return;
		
		for (int c=1;c < numCells; c++) {
			Cell cell = row.getCell(c);
			if (cell != null) {
				String accountNumber = cell.toString().trim();
				if (accountNumber.length()!=0) awsAccountGroupName.put(accountNumber, groupName);
			}
		} 
	}

	public TreeMap<String, String> getAwsAccountGroupName() {
		return awsAccountGroupName;
	}


}
