package edu.utah.billing;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import util.gen.IO;
import util.gen.Misc;


public class MasterAccountInfoParser2 {

	private File masterAccountFile = null;
	private ArrayList<BillingGroup> billingGroups = new ArrayList<BillingGroup>();
	private HashMap<String, BillingGroup> aliasesBillingGroups = new HashMap<String, BillingGroup>();

	public MasterAccountInfoParser2(File xlsx, boolean debug) throws IOException {
		masterAccountFile = xlsx;
		
		parseIt();
		
		//load all of the aliases
		for (BillingGroup bg: billingGroups) {
			for (String a: bg.getAliases()) {
				if (aliasesBillingGroups.containsKey(a)) throw new IOException("ERROR: "+a+" was found on multiple lines! Fix "+masterAccountFile);
				aliasesBillingGroups.put(a, bg);
			}
		}
		
		if (debug) {
			for (BillingGroup bg: billingGroups) {
				IO.pl("Found:\t"+ bg);
			}
		}
	}

	public static void main(String[] args) throws IOException {
		MasterAccountInfoParser2 p = new MasterAccountInfoParser2(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2024/3-4_CBI_Feb-Mar_2024/masterAccountInfo.xlsx"), true);
		p.saveUpdatedInfo();
	}
	
	/**Prints a txt formatted spreadsheet that can be inspected and then used to replace the old master xlsx file.*/
	public void saveUpdatedInfo(){
		File updatedMasterAccountFile = null;
		try {
			String originalName = Misc.removeExtension(masterAccountFile.getName());
			updatedMasterAccountFile = new File(masterAccountFile.getParentFile(), originalName+"_Updated.xls");
			IO.pl("\nSaving updated master account info, review in Excel and use it for the next billing cycle, "+updatedMasterAccountFile.getName());
			PrintWriter out = new PrintWriter( new FileWriter(updatedMasterAccountFile));
			//header
			out.println("Cancer Status [Cancer|Non-Cancer]\tTotalHoursBilled\tGroupAliasName1\tGroupAliasName2\tEtc");
			//for each group
			for (BillingGroup bg: billingGroups) out.println(bg);
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
			int numRows = sheet.getPhysicalNumberOfRows()+1;
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
		if (numCells == 1) return;
		
		
		//row 0 is Cancer Status [Cancer|Non-Cancer]
		String cancerStatus = row.getCell(0).toString();
		if (cancerStatus.startsWith("Cancer Status")) return; //skip header
		boolean isCancerMember = true;
		if (cancerStatus.toLowerCase().contains("non")) isCancerMember = false;
		
		//row 1 is Total Hours Billed, for <=70 hrs its $70/hr for more its $100
		Float totalHours = Float.parseFloat(row.getCell(1).toString());
		
		LinkedHashSet<String> aliases = new LinkedHashSet<String>();
		for (int c=2;c < numCells; c++) {
			Cell cell = row.getCell(c);
			if (cell != null) {
				String value = cell.toString().trim();
				if (value.length()!=0) aliases.add(value);
			}
		} 
		
		billingGroups.add( new BillingGroup(isCancerMember, totalHours, aliases) );
	}

	public ArrayList<BillingGroup> getBillingGroups() {
		return billingGroups;
	}

	public HashMap<String, BillingGroup> getAliasesBillingGroups() {
		return aliasesBillingGroups;
	}

	
}
