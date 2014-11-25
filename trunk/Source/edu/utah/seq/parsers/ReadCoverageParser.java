package edu.utah.seq.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.charts.AxisCrosses;
import org.apache.poi.ss.usermodel.charts.AxisPosition;
import org.apache.poi.ss.usermodel.charts.ChartAxis;
import org.apache.poi.ss.usermodel.charts.ChartDataSource;
import org.apache.poi.ss.usermodel.charts.DataSources;
import org.apache.poi.ss.usermodel.charts.LegendPosition;
import org.apache.poi.ss.usermodel.charts.LineChartData;
import org.apache.poi.ss.usermodel.charts.ValueAxis;
import org.apache.poi.ss.util.CellRangeAddress;
import org.apache.poi.xssf.usermodel.XSSFCell;
import org.apache.poi.xssf.usermodel.XSSFChart;
import org.apache.poi.xssf.usermodel.XSSFClientAnchor;
import org.apache.poi.xssf.usermodel.XSSFDrawing;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.apache.poi.xssf.usermodel.charts.XSSFChartLegend;
import org.openxmlformats.schemas.drawingml.x2006.chart.CTTitle;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class ReadCoverageParser {

	//fields
	private File[] dataFiles;
	private File saveFile;
	private int coverageThreshold = 10;
	private double minimumFraction = 0.95;
	private String[] names;
	private String[] trimmedNames;
	private float[][] matrix;
	private ArrayList<String> failedNames = new ArrayList<String>();
	private String headerWord = "FractionObservedWithGivenOrMoreCoverage";

	//constructors
	public ReadCoverageParser(String[] args){

		long startTime = System.currentTimeMillis();
		processArgs(args);

		loadNames();

		loadMatrix();

		flagDatasets();
		System.out.println("\nFailed samples: ");
		if (failedNames.size()!=0)Misc.printArray(failedNames);
		else System.out.println("\tAll passed");

		saveMatrix(saveFile);
		
		makeExcelDoc();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" Sec\n");
	}




	private void saveMatrix(File file) {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(file));
			
			//write header
			out.print(trimmedNames[0]);
			for (int i=1; i< trimmedNames.length; i++){
				out.print("\t");
				out.print(trimmedNames[i]);
			}
			out.println();
			
			//write maxtrix lines
			//for each row
			for (int j=0; j< matrix[0].length; j++){
				out.print(matrix[0][j]);
				for (int i=1; i< trimmedNames.length; i++){
					out.print("\t");
					out.print(matrix[i][j]);
				}
				out.println();
			}
			out.close();
		} catch (IOException e) {
			System.err.println("\nError: problem saving matrix file.");
			e.printStackTrace();
		}
		
	}

	private void flagDatasets() {
		int index = coverageThreshold -1;
		for (int i=0; i< names.length; i++){
			if (matrix[i][index] < minimumFraction) failedNames.add(names[i]);
		}
	}

	private void loadMatrix() {
		matrix = new float[names.length][100];
		for (int i=0; i< dataFiles.length; i++){
			matrix[i] = loadReadCoverage(dataFiles[i]);
		}
	}
	
	private void makeExcelDoc() {
		try {
			/* Create a Workbook object that will hold the final chart */
			XSSFWorkbook workbook = new XSSFWorkbook();
			/* Create a worksheet object for the line chart. This worksheet will contain the chart */
			XSSFSheet worksheet = workbook.createSheet("Read Coverage");

			//add description
			CellStyle style = workbook.createCellStyle();
			Font font = workbook.createFont();
			font.setBoldweight(Font.BOLDWEIGHT_BOLD);
			style.setFont(font);
			
			XSSFRow desRow = worksheet.createRow((short)0);
			XSSFCell desCell = desRow.createCell((short)0);
			desCell.setCellValue("Fraction of interrogated bps with X or greater fold coverage");
			desCell.setCellStyle(style);
			
			//any failing samples
			XSSFRow sampleRow = worksheet.createRow((short)1);
			XSSFCell sampleCell = sampleRow.createCell((short)1);
			String des = "Samples failing thresholds (min frac "+Num.formatNumber(minimumFraction, 2)+", min fold "+coverageThreshold+") "+failedNames.size()+": ";
			sampleCell.setCellValue(des);
			if (failedNames.size()!=0){
				XSSFCell badSampleNames = sampleRow.createCell((short)6);
				badSampleNames.setCellValue(Misc.stringArrayListToString(failedNames, ", "));
				badSampleNames.setCellStyle(style);
			}
			
			//add trimmed names row
			XSSFRow nameRow = worksheet.createRow((short)3);
			XSSFCell nameCell = nameRow.createCell((short)0);
			nameCell.setCellValue("FoldCoverage");
			for (int i=0; i< trimmedNames.length; i++){
				XSSFCell cell = nameRow.createCell((short)(i+1));
				cell.setCellValue(trimmedNames[i]);
			}
			//create table
			for (int rowIndex = 0; rowIndex < matrix[0].length; rowIndex++){
				/* Add a row that contains the chart data */
				XSSFRow my_row = worksheet.createRow((short)rowIndex+4);
				//add coverage 1x,2x,3x ...
				XSSFCell cov = my_row.createCell((short)0);
				cov.setCellValue((short)(rowIndex +1));
				for (int colIndex = 0; colIndex < matrix.length; colIndex++){
					/* Define column values for the row that is created */
					XSSFCell cell = my_row.createCell((short)(colIndex+1));
					cell.setCellValue(matrix[colIndex][rowIndex]);
				}
			}

			/* At the end of this step, we have a worksheet with test data, that we want to write into a chart */
			/* Create a drawing canvas on the worksheet */
			XSSFDrawing xlsx_drawing = worksheet.createDrawingPatriarch();
			/* Define anchor points in the worksheet to position the chart */
			XSSFClientAnchor anchor = xlsx_drawing.createAnchor(0, 0, 0, 0, 2, 5, 20, 40);
			/* Create the chart object based on the anchor point */
			XSSFChart chart = xlsx_drawing.createChart(anchor);
			/* Define legends for the line chart and set the position of the legend */
			if (trimmedNames.length<11){
				XSSFChartLegend legend = chart.getOrCreateLegend();
				legend.setPosition(LegendPosition.RIGHT); 
			}
			/* Create data for the chart */
			LineChartData data = chart.getChartDataFactory().createLineChartData();     
			/* Define chart AXIS */
			ChartAxis bottomAxis = chart.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis = chart.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			//leftAxis.setCrosses(AxisCrosses.AUTO_ZERO); 
			leftAxis.setMaximum(1.0);
			leftAxis.setMinimum(0.4);
			
			/* Define Data sources for the chart */
			/* Set the right cell range that contain values for the chart */
			/* Pass the worksheet and cell range address as inputs */
			/* Cell Range Address is defined as First row, last row, first column, last column */
			//make source for the x axis
			ChartDataSource<Number> xs = DataSources.fromNumericCellRange(worksheet, new CellRangeAddress(4, 53, 0, 0));
			
			//for each dataset
			for (int i=0; i< trimmedNames.length; i++){
				ChartDataSource<Number> sample = DataSources.fromNumericCellRange(worksheet, new CellRangeAddress(4, 53, 1+i, 1+i));
				data.addSerie(xs, sample).setTitle(trimmedNames[i]);
			}

			/* Plot the chart with the inputs from data and chart axis */
			chart.plot(data, new ChartAxis[] { bottomAxis, leftAxis });
			/* Finally define FileOutputStream and write chart information */               
			FileOutputStream fileOut = new FileOutputStream("/Users/u0028003/Desktop/delme3.xlsx");
			workbook.write(fileOut);
			fileOut.close();
		} catch (Exception e){
			System.err.println("\nError writing excel spreadsheet.\n");
			e.printStackTrace();
		}
	}

	/*Parses the following from each file
BaseCoverage	ObservedBasesWithGivenCoverage	FractionObserved	FractionObservedWithGivenOrMoreCoverage
0	9785.0	0.000	1.000
1	8252.0	0.000	1.000
2	11706.0	0.000	0.999
3	15316.0	0.001	0.999
4	19162.0	0.001	0.998
5	23075.0	0.001	0.997
6	27094.0	0.001	0.996
	 */
	private float[] loadReadCoverage(File file) {
		float[] fracCov = new float[100];
		try {
			BufferedReader in = IO.fetchBufferedReader(file);
			String line;
			String[] fields;
			//find header line
			boolean failed = true;
			while ((line = in.readLine()) !=null){
				if (line.contains(headerWord)) {
					failed = false;
					break;
				}
			}
			if (failed) Misc.printErrAndExit("\nError: failed to find header in "+file.getName()+". Is this a read coverage output file?\n");
			
			//skip first
			line = in.readLine();
			
			//extract next 100
			for (int i=0; i< fracCov.length; i++) {
				line = in.readLine();
				fields = Misc.TAB.split(line);
				fracCov[i] = Float.parseFloat(fields[3]);
			}

			in.close();
		} catch (IOException e) {
			System.err.println("\nError: problem parsing "+file);
			e.printStackTrace();
		}
		return fracCov;
	}

	private void loadNames() {
		// TODO Auto-generated method stub
		names = new String[dataFiles.length];
		for (int i=0; i< names.length; i++) names[i] = Misc.removeExtension(dataFiles[i].getName());
		trimmedNames = Misc.trimCommon(names);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ReadCoverageParser(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		File dir = null;
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': dir = new File(args[++i]); break;
					case 's': saveFile = new File(args[++i]); break;
					case 'c': coverageThreshold = Integer.parseInt(args[++i]); break;
					case 'm': minimumFraction = Double.parseDouble(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (dir == null || dir.isDirectory() == false || dir.canRead() == false) Misc.printExit("\nError: cannot find your data file directory?\n");
		
		//check save file
		if (saveFile == null) saveFile = new File (dir,"readCoverageReport.xlxs");
		if (saveFile.exists()) saveFile.delete();
		
		dataFiles = IO.extractOnlyFiles(dir); 
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find any data files to parse?\n");


		//print info
		System.out.println(coverageThreshold+ "\tRead coverage flag");
		System.out.println(minimumFraction+ "\tMinimum fraction interrogated bases at the read coverage to avoid flagging");
		System.out.println(saveFile.getName() +"\tSave file");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Read Coverage Parser: Nov 2014                          **\n" +
				"**************************************************************************************\n" +
				"Parses and plots read coverage output from Sam2USeq. Also flags datasets that fail the\n"+
				"coverage thresholds. Generates a native xlsx spreadsheet with interactive graph.\n"+
				
				"\nOptions:\n"+
				"-d Directory containing Sam2USeq log files with read coverage information.\n" +

				"\nDefault Options:\n"+
				"-s Save file, defaults writing log to -d.\n"+
				"-c Minimum coverage threshold, defaults to 10.\n"+
				"-m Minimum fraction interrogated bases at the coverage threshold, defaults to 0.95\n"+
				"\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/ReadCoverageParser -d /QC/Sam2USeqLogs/\n" +
				"     -c 15 -m 0.9  \n\n" +

				"**************************************************************************************\n");

	}	

}
