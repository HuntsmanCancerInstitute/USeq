package edu.utah.billing;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.TreeMap;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import util.gen.IO;


public class AliasXlsxParser {

	private TreeMap<String, HashSet<String>> groupNameAliases = new TreeMap<String, HashSet<String>>();

	public AliasXlsxParser(File xlsx, boolean debug) {
		parseIt(xlsx);
		
		if (debug) {
			for (String s: groupNameAliases.keySet()) {
				IO.pl(s+"\t"+groupNameAliases.get(s));
			}
		}
	}

	public static void main(String[] args) {
		AliasXlsxParser p = new AliasXlsxParser(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2023/6_BSR_June_2023/accountAliases.xlsx"), true);

	}

	private void parseIt(File inputFile) {
		try {

			//Open up xlsx file
			Workbook wb = WorkbookFactory.create(inputFile);	

			//Find appropriate sheet 
			Sheet sheet = wb.getSheet("Aliases");
			if (sheet == null) throw new IOException("Could not find sheet 'Aliases' in "+inputFile+" ?");

			//Iterate through rows
			int numRows = sheet.getPhysicalNumberOfRows();
			for (int r = 0; r< numRows; r++) {
				Row row = sheet.getRow(r);
				if (row != null) {
					addAliases(row);
				}
			} 
		} catch (Exception e) {
			System.out.println("Aliases xlsx file is not in the correct format, exiting");
			e.printStackTrace();
			System.exit(1);
		} 
	}

	private void addAliases(Row row) {
		int numCells = row.getLastCellNum()+1;
		
		HashSet<String> toKeep = new HashSet<String>();
		for (int c=0;c < numCells; c++) {
			Cell cell = row.getCell(c);
			if (cell != null) {
				String value = cell.toString().trim();
				//header line?
				if (value.startsWith("#")) return;
				if (value.length()!=0) toKeep.add(value);
			}
		} 
		//Save em
		for (String s: toKeep) groupNameAliases.put(s, toKeep);
	}

	public TreeMap<String, HashSet<String>> getGroupNameAliases() {
		return groupNameAliases;
	}
}
