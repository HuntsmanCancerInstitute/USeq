package edu.utah.seq.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
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

public class CollectBamStats {

	//fields
	private File[] dataFiles;
	private File saveFile;
	private String[] names;
	private String[] trimmedNames;
	private float[][] readCoverageMatrix;
	private float[][] alignmentStatMatrix;
	private LinkedHashMap <String, ArrayList<String>> failedNames = new LinkedHashMap <String, ArrayList<String>>();
	private String readCoverageHeaderWord = "FractionObservedWithGivenOrMoreCoverage";
	private String mergePairedHeaderWord = "be suspicious of zero read catagories";

	//alignmentStats
	private int numberAlignmentStats;
	private HashMap<String, Integer> nameIndexAlignmentStats = null;
	private String[] alignmentStatColumnNames = null;
	private int unmappedIndex;
	private int duplicateIndex;
	private int passAllFiltersIndex;
	private int overlappingBasesIndex;
	private float maximumUnmapped;
	private float maximumDuplicate;
	private float minimumPassAll;
	private float maximumOverlappingBases;

	//read coverage
	private int coverageThreshold;
	private float minimumCoverage;

	//constructors
	public CollectBamStats(String[] args){

		long startTime = System.currentTimeMillis();
		processArgs(args);
		loadNames();
		loadAlignmentStatArrays();
		

		loadAlignmentStatMatrix();
		flagAlignmentStatDatasets();
		printAlignmentStats();

		loadReadCoverageMatrix();
		flagReadCoverageDatasets();
		printReadCoverageMatrix();

		printThresholds();
		
		System.out.println("\nFailed samples: ");
		if (failedNames.size()!=0) Misc.printHashMap(failedNames, ", ");
		else System.out.println("\tAll passed");

		//saveMatrix(saveFile);
		//makeExcelDoc();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" Sec\n");
	}


	private void printThresholds() {
		System.out.println("\nThresholds for failing samples:");
		System.out.println(maximumUnmapped+"\tMaximum fraction unmapped reads");
		System.out.println(maximumDuplicate+"\tMaximum fraction duplicate reads");
		System.out.println(minimumPassAll+"\tMinimum fraction passing alignments");
		System.out.println(maximumOverlappingBases+"\tMaximum fraction overlapping bps");
		System.out.println(coverageThreshold+"\tCoverage threshold");
		System.out.println(minimumCoverage+"\tMinimum fraction exonic bases that meet the Coverage threshold");
	}


	private void printAlignmentStats() {
		System.out.println("Alignment statistics:");
		System.out.println("Names\t"+Misc.stringArrayToString(alignmentStatColumnNames, "\t"));
		for (int i=0; i< names.length; i++){
			System.out.println(names[i]+"\t"+Num.floatArrayToString(alignmentStatMatrix[i], "\t"));
		}

	}


	private void flagAlignmentStatDatasets() {
		//for each sample
		for (int i=0; i<alignmentStatMatrix.length; i++ ){
			ArrayList<String> al = new ArrayList<String>();
			float[] stats = alignmentStatMatrix[i];
			//check unmapped
			if (stats[unmappedIndex] > maximumUnmapped) al.add("Unmapped");
			//check duplicate
			if (stats[duplicateIndex] > maximumDuplicate) al.add("Duplicates");
			//check passAllFilters
			if (stats[passAllFiltersIndex] < minimumPassAll) al.add("Passing Alignments");
			//check overlappingBases
			if (stats[overlappingBasesIndex] > maximumOverlappingBases) al.add("Overlapping Bases");
			//save it?
			if (al.size() !=0){
				if (failedNames.containsKey(names[i]) == false){
					failedNames.put(names[i], al);
				}
				else failedNames.get(names[i]).addAll(al);
			}
		}

	}


	private void loadAlignmentStatArrays(){
		numberAlignmentStats = 10;

		String[] columnNames = {
				"Total # alignments",
				"Unmapped reads",
				"Alignments failing alignment score",
				"Alignments failing mapping quality score",
				"Adapter alignments",
				"PhiX alignments",
				"Duplicate alignments",
				"Alignments passing all filters",
				"Fraction overlapping bases in paired alignments",
				"Median Coverage"
		};
		nameIndexAlignmentStats = new HashMap<String, Integer>();
		for (int i=0; i<columnNames.length; i++) nameIndexAlignmentStats.put(columnNames[i], i);

		//short version
		alignmentStatColumnNames = new String[] {
				"Total",
				"Unmapped",
				"Failing alignment score",
				"Failing mapping quality score",
				"Adapter",
				"PhiX",
				"Duplicate",
				"Passing all filters",
				"Overlapping bases",
				"Exonic Coverage"
		};

		//indexes for the scores;
		unmappedIndex = 1;
		duplicateIndex = 6;
		passAllFiltersIndex = 7;
		overlappingBasesIndex = 8;

		//thresholds for alignment stats
		maximumUnmapped = 0.01f;
		maximumDuplicate = 0.15f;
		minimumPassAll = 0.8f;
		maximumOverlappingBases = 0.1f;

		//read coverage thresholds
		coverageThreshold = 10;
		minimumCoverage = 0.95f;

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
			for (int j=0; j< readCoverageMatrix[0].length; j++){
				out.print(readCoverageMatrix[0][j]);
				for (int i=1; i< trimmedNames.length; i++){
					out.print("\t");
					out.print(readCoverageMatrix[i][j]);
				}
				out.println();
			}
			out.close();
		} catch (IOException e) {
			System.err.println("\nError: problem saving matrix file.");
			e.printStackTrace();
		}

	}

	private void printReadCoverageMatrix() {
		//write header
		System.out.println("\nRead coverage statistics:");
		System.out.print("XCoverage");
		for (int i=0; i< trimmedNames.length; i++){
			System.out.print("\t");
			System.out.print(trimmedNames[i]);
		}
		System.out.println();

		//write maxtrix lines
		//for each row
		for (int j=0; j< readCoverageMatrix[0].length; j++){
			System.out.print(j+1);
			float total = 0;
			for (int i=0; i< trimmedNames.length; i++){
				System.out.print("\t");
				System.out.print(Math.abs(readCoverageMatrix[i][j]));
				total += readCoverageMatrix[i][j];
			}
			System.out.println();
			if (total == 0f) break;
			
		}
	}

	private void flagReadCoverageDatasets() {
		int index = coverageThreshold -1;
		for (int i=0; i< names.length; i++){
			if (readCoverageMatrix[i][index] < minimumCoverage) {
				ArrayList<String> al = failedNames.get(names[i]);
				if (al == null) {
					al = new ArrayList<String>();
					failedNames.put(names[i], al);
				}
				al.add("Read Coverage");
			}
		}
	}

	private void loadReadCoverageMatrix() {
		readCoverageMatrix = new float[names.length][100];
		for (int i=0; i< dataFiles.length; i++){
			readCoverageMatrix[i] = loadReadCoverage(dataFiles[i]);
		}
	}

	private void loadAlignmentStatMatrix() {
		alignmentStatMatrix = new float[names.length][numberAlignmentStats];
		for (int i=0; i< dataFiles.length; i++){
			alignmentStatMatrix[i] = loadAlignmentStats(dataFiles[i]);
		}
	}

	private void makeExcelDoc() {
		try {
			/* Create a Workbook object that will hold the final chart */
			XSSFWorkbook workbook = new XSSFWorkbook();
			/* Create a worksheet object for the line chart. This worksheet will contain the chart */
			XSSFSheet worksheet = workbook.createSheet("Read Coverage");

			//add description
			CellStyle styleBold = workbook.createCellStyle();
			Font font = workbook.createFont();
			font.setBoldweight(Font.BOLDWEIGHT_BOLD);
			styleBold.setFont(font);

			XSSFRow desRow = worksheet.createRow((short)0);
			XSSFCell desCell = desRow.createCell((short)0);
			desCell.setCellValue("Fraction of interrogated bps with X or greater fold coverage");
			desCell.setCellStyle(styleBold);

			//any failing samples, need to fix
			XSSFRow sampleRow = worksheet.createRow((short)1);
			XSSFCell sampleCell = sampleRow.createCell((short)1);
			String des = "Samples failing thresholds (min frac "+Num.formatNumber(minimumCoverage, 2)+", min fold "+coverageThreshold+") "+failedNames.size()+": ";
			sampleCell.setCellValue(des);
			if (failedNames.size()!=0){
				XSSFCell badSampleNames = sampleRow.createCell((short)6);
				//badSampleNames.setCellValue(Misc.stringArrayListToString(failedNames, ", "));
				//badSampleNames.setCellStyle(style);
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
			for (int rowIndex = 0; rowIndex < readCoverageMatrix[0].length; rowIndex++){
				/* Add a row that contains the chart data */
				XSSFRow my_row = worksheet.createRow((short)rowIndex+4);
				//add coverage 1x,2x,3x ...
				XSSFCell cov = my_row.createCell((short)0);
				cov.setCellValue((short)(rowIndex +1));
				for (int colIndex = 0; colIndex < readCoverageMatrix.length; colIndex++){
					/* Define column values for the row that is created */
					XSSFCell cell = my_row.createCell((short)(colIndex+1));
					cell.setCellValue(readCoverageMatrix[colIndex][rowIndex]);
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
				if (line.contains(readCoverageHeaderWord)) {
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

	/*Parses the following from the file

	 */
	private float[] loadAlignmentStats(File file) {
		float[] stats = new float[numberAlignmentStats];
		try {
			BufferedReader in = IO.fetchBufferedReader(file);
			String line;
			String[] fields;
			//find header line
			boolean failed = true;
			while ((line = in.readLine()) !=null){
				if (line.contains(mergePairedHeaderWord)) {
					failed = false;
					break;
				}
			}
			if (failed) Misc.printErrAndExit("\nError: failed to find header in "+file.getName()+". Does this contain a MergePairedAlignment log output?\n");

			//parse out remainder
			boolean lastNotFound = true;
			while ((line = in.readLine()) !=null){
				line = line.trim();
				if (line.length()==0) continue;
				fields = Misc.TAB.split(line);
				if (fields.length < 2 || fields.length > 3) continue;
				//look for a match
				Integer index = nameIndexAlignmentStats.get(fields[0]);
				if (index != null) {
					float last = Float.parseFloat(fields[fields.length-1]);
					stats[index.intValue()] = last;
				}
				//look for last
				if (fields[0].equals("Median Coverage")){
					lastNotFound = false;
					break;
				}
			}

			//last one found?
			if (lastNotFound){
				Misc.printErrAndExit("\nError, failed to find the 'Median Coverage' stat in the log file? Does this not contain the output of the Sam2USeq app? "+file.getName());
			}


			in.close();
		} catch (IOException e) {
			System.err.println("\nError: problem parsing "+file);
			e.printStackTrace();
		}
		return stats;
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
		new CollectBamStats(args);
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
					case 'l': dir = new File(args[++i]); break;
					//case 's': saveFile = new File(args[++i]); break;
					case 'x': coverageThreshold = Integer.parseInt(args[++i]); break;
					case 'c': minimumCoverage = Float.parseFloat(args[++i]); break;
					case 'u': maximumUnmapped = Float.parseFloat(args[++i]); break;
					case 'd': maximumDuplicate = Float.parseFloat(args[++i]); break;
					case 'p': minimumPassAll = Float.parseFloat(args[++i]); break;
					case 'o': maximumOverlappingBases = Float.parseFloat(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (dir == null || dir.isDirectory() == false || dir.canRead() == false) Misc.printExit("\nError: cannot find your log directory?\n");

		//check save file
		if (saveFile == null) saveFile = new File (dir,"readCoverageReport.xlsx");
		if (saveFile.exists()) saveFile.delete();

		dataFiles = IO.extractOnlyFiles(dir); 
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find any data files to parse?\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Collect Bam Stats: Dec 2014                         **\n" +
				"**************************************************************************************\n" +
				"Parses and plots bam alignment quality statistics from a log file containing the\n"+
				"output of the MergePairedAlignments and Sam2USeq apps. Will flag datasets that fail\n"+
				"the set thresholds.\n"+// Generates a native xlsx spreadsheet with interactive graph.\n"+

				"\nOptions:\n"+
				"-l Directory containing a combine log file of MergePairedAlignments and Sam2USeq,\n"+
				"      one per sample.\n" +

				"\nDefault Options:\n"+
				//"-s Save file, defaults writing log to -d.\n"+
				"-x Minimum alignment coverage threshold, defaults to 10.\n"+
				"-c Minimum fraction interrogated bases at the coverage threshold, defaults to 0.95\n"+
				"-u Maximum fraction unmapped reads, defaults to 0.01\n"+
				"-d Maximum fraction duplicate reads, defaults to 0.15\n"+
				"-p Minimum fraction passing alignments, defaults to 0.8\n"+
				"-o Maximum fraction overlapping bps in paired alignments, defaults to 0.1\n"+
				"\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/CollectBamStats -l /QC/Sam2USeqLogs/\n" +
				"     -x 15 -c 0.9  \n\n" +

				"**************************************************************************************\n");

	}	
}
