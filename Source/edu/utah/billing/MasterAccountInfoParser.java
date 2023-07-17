package edu.utah.billing;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.TreeMap;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import edu.utah.seq.analysis.multi.PairedCondition;
import util.gen.IO;
import util.gen.Misc;


public class MasterAccountInfoParser {
	
	private File masterAccountFile = null;
	
	//all groupNames are in this hash
	private HashMap<String, HashSet<String>> groupNameAliases = new HashMap<String, HashSet<String>>();
	
	//the groupName here is the first in the file
	private TreeMap<String, String> groupNameCancerStatus = new TreeMap<String, String>();
	private HashMap<String, Float> groupNameHoursBilled = new HashMap<String, Float>();

	public MasterAccountInfoParser(File xlsx, boolean debug) {
		this.masterAccountFile = xlsx;
		
		parseIt();
		
		if (debug) {
			for (String s: groupNameCancerStatus.keySet()) {
				IO.pl(groupNameCancerStatus.get(s)+ "\t"+ groupNameHoursBilled.get(s)+"\t"+groupNameAliases.get(s));
			}
		}
	}

	public static void main(String[] args) {
		MasterAccountInfoParser p = new MasterAccountInfoParser(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2023/6_BSR_June_2023/masterAccountInfo.xlsx"), true);
		p.saveUpdatedInfo();
	}
	
	/**Prints a txt formatted spreadsheet that can be inspected and then used to replace the old master xlsx file.*/
	public void saveUpdatedInfo(){
		File updatedMasterAccountFile = null;
		try {
			String originalName = Misc.removeExtension(masterAccountFile.getName());
			updatedMasterAccountFile = new File(masterAccountFile.getParentFile(), originalName+"_Updated.xls");
			IO.pl("\nSaving updated master account info, review in Excel and replace the original, "+updatedMasterAccountFile.getName());
			PrintWriter out = new PrintWriter( new FileWriter(updatedMasterAccountFile));
			//header
			out.println("Cancer Status [Cancer|Non-Cancer]\tTotalHoursBilled\tGroupAliasName1\tGroupAliasName2\tEtc");
			//for each group
			for (String gn: groupNameCancerStatus.keySet()) {
				StringBuilder sb = new StringBuilder(groupNameCancerStatus.get(gn));
				sb.append("\t");
				sb.append(groupNameHoursBilled.get(gn).toString());
				//aliases
				for (String a: groupNameAliases.get(gn)) {
					sb.append("\t");
					sb.append(a);
				}
				out.println(sb.toString());
			}
			out.close();
			
		} catch (Exception e){
			System.err.println("\nProblem saving updated "+masterAccountFile);
			if (updatedMasterAccountFile != null) updatedMasterAccountFile.delete();
			e.printStackTrace();
			System.exit(1);
		}
	}


	private void parseIt() {
		try {

			//Open up xlsx file
			Workbook wb = WorkbookFactory.create(masterAccountFile);	

			//Find appropriate sheet 
			Sheet sheet = wb.getSheetAt(0);
			if (sheet == null) throw new IOException("Could not find a sheet in "+masterAccountFile+" ?");

			//Iterate through rows
			int numRows = sheet.getPhysicalNumberOfRows();
			for (int r = 0; r< numRows; r++) {
				Row row = sheet.getRow(r);
				if (row != null) parseRow(row);
			} 
		} catch (Exception e) {
			System.out.println("MasterAccountInfo xlsx parsing failed, exiting");
			e.printStackTrace();
			System.exit(1);
		} 
	}

	private void parseRow(Row row) {
		int numCells = row.getLastCellNum()+1;
		LinkedHashSet<String> toKeep = new LinkedHashSet<String>();
		
		//row 0 is Cancer Status [Cancer|Non-Cancer]
		String cancerStatus = row.getCell(0).toString();
		if (cancerStatus.startsWith("Cancer Status")) return; //skip header
		
		//row 1 is Total Hours Billed, for <=70 hrs its $70/hr for more its $100
		Float totalHours = Float.parseFloat(row.getCell(1).toString());
		
		for (int c=2;c < numCells; c++) {
			Cell cell = row.getCell(c);
			if (cell != null) {
				String value = cell.toString().trim();
				if (value.length()!=0) toKeep.add(value);
			}
		} 
		//Save em
		String firstGroupName = null;
		for (String s: toKeep) {
			groupNameAliases.put(s, toKeep);
			if (firstGroupName == null) firstGroupName = s;
		}
		groupNameCancerStatus.put(firstGroupName, cancerStatus);
		groupNameHoursBilled.put(firstGroupName, totalHours);
	}

	public HashMap<String, HashSet<String>> getGroupNameAliases() {
		return groupNameAliases;
	}

	public File getMasterAccountFile() {
		return masterAccountFile;
	}

	public TreeMap<String, String> getGroupNameCancerStatus() {
		return groupNameCancerStatus;
	}

	public HashMap<String, Float> getGroupNameHoursBilled() {
		return groupNameHoursBilled;
	}
}
