package edu.utah.ames.bioinfo;

import java.io.File;
import java.io.FileOutputStream;

import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

public class MakeXlsxWorkbook {

	private Workbook workbook = null;
	private String saveDirectory = null;
	private String outputName = null;
	
	public void printStatSpreadSheed() {
		try {
			//create workbook
			workbook = new XSSFWorkbook();
			FileOutputStream fileOut = new FileOutputStream(new File(saveDirectory, outputName + ".xlsx"));
			
			//create a worksheet
			Sheet sheet0 = workbook.createSheet("Contamination stats");
			sheet0.createFreezePane( 0, 1, 0, 1 );
			
			//create row, one for each item plus header
			//Row[] rows = new Row[]
		}
		catch (Exception e) {
			System.err.println("\nProblem printing the spreadsheet!");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
}
