package edu.utah.billing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;

import util.gen.IO;
import util.gen.Misc;

public class WafXlsxParser {
	private HashMap<String, Integer> headerKeyIndex = new HashMap<String, Integer>();
	private ArrayList<String> header = new ArrayList<String>();
	private String headerTabbed = null;
	private TreeMap<String, ArrayList<String[]>> groupNameWafLines = new TreeMap<String, ArrayList<String[]>>();
	private boolean debug = false;


	public WafXlsxParser(File xlsx, boolean debug) {
		this.debug = debug;
		
		parseIt(xlsx);
		headerTabbed = Misc.stringArrayListToString(header, "\t");

		if (debug) {
			IO.pl("Header:\t"+ header);
			IO.pl("\nHeaderIndexes:\t"+ headerKeyIndex);
			IO.pl("\nPIs: "+groupNameWafLines.size());
			for (String pi: groupNameWafLines.keySet()) {
				IO.pl(pi);
				for (String[] v: groupNameWafLines.get(pi)) {
					IO.pl("\t"+ Misc.stringArrayToString(v, "\t"));
				}
			}
		}

	}

	public static void main(String[] args) {
		//hourly
		new WafXlsxParser(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2023/6_BSR_June_2023/WAF/Bioinformatics WAF Tracking Schedule - FY23.xlsx"), true);

		//cloud
		new WafXlsxParser(new File ("/Users/u0028003/HCI/CoreAdmin/Billing/AllBillingReports/2023/6_BSR_June_2023/WAF/Bioinformatics SB Cloud WAF Tracking Schedule - FY23.xlsx"), true);
	}

	private void parseIt(File inputFile) {
		try {

			//Open up xlsx file
			Workbook wb = WorkbookFactory.create(inputFile);	

			//Find appropriate sheet 
			Sheet sheet = wb.getSheet("Account Table");
			if (sheet == null) throw new IOException("Could not find sheet 'Account Table' in "+inputFile+" ?");

			//Iterate through rows
			int numRows = sheet.getPhysicalNumberOfRows();
			for (int r = 0; r< numRows; r++) {
				Row row = sheet.getRow(r);
				if (row == null) {
					//if (debug) IO.pl();
				}
				else {
					int numCells = row.getLastCellNum()+1;
					boolean inHeader = false;
					String[] cellsToSave = new String[numCells+1];
					for (int c=0;c < numCells; c++) {
						//Get cell
						Cell cell = sheet.getRow(r).getCell(c);
						if (cell != null) {
							//trim it and replace any returns
							String cellStringValue = cell.toString().trim();
							cellStringValue = Misc.RETURN.matcher(cellStringValue).replaceAll("");
							
							//if (debug) IO.p(cellStringValue+"\t");
							if (cellStringValue.contains("INVESTIGATOR") || inHeader == true) {
								if (cellStringValue.contains("NOTES")) {
									inHeader = false;
									c = numCells;
								}
								else inHeader = true;
								//if (debug) IO.pl("Header found! Adding "+cellStringValue+" to "+c);
								headerKeyIndex.put(cellStringValue, c);
								header.add(cellStringValue);
							}
							if (cellStringValue.contains("EXPIRED")) {
								//if (debug) IO.pl("No more active WAFS");
								return;
							}
							cellsToSave[c] = cellStringValue;
						}
						else {
							//if (debug) IO.p("\t");
							if (inHeader) header.add(" ");
						}
					} 
					//if (debug) IO.pl();
					//save?
					if (cellsToSave[0]!=null && cellsToSave[0].length()!=0 && cellsToSave[0].equals("INVESTIGATOR")==false) {
						ArrayList<String[]> al = groupNameWafLines.get(cellsToSave[0]);
						if (al == null) {
							al = new ArrayList<String[]>();
							groupNameWafLines.put(cellsToSave[0], al);
						}
						al.add(cellsToSave);
					}

				}
			} 
		} catch (Exception e) {
			System.out.println("Xlsx file is not in the correct format, exiting -> "+inputFile);
			e.printStackTrace();
			System.exit(1);
		} 
	}

	public TreeMap<String, ArrayList<String[]>> getGroupNameWafLines() {
		return groupNameWafLines;
	}

	public ArrayList<String> getHeader() {
		return header;
	}

	public String getHeaderTabbed() {
		return headerTabbed;
	}
}
