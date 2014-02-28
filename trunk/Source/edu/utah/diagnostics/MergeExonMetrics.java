package edu.utah.diagnostics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

import org.apache.poi.hssf.usermodel.HSSFClientAnchor;
import org.apache.poi.hssf.usermodel.HSSFPicture;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.ss.util.CellRangeAddress;
import org.apache.poi.ss.util.RegionUtil;
import org.apache.poi.util.IOUtils;
import org.apache.poi.xssf.usermodel.XSSFShape;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;


public class MergeExonMetrics {
	private File metricsDir = null;
	private ArrayList<File> metricsList = null;
	private ArrayList<HashMap<String,File>> fileList = null;
	private String prefix = null;
	private ArrayList<HashMap<String,String>> lookupTables = new ArrayList<HashMap<String,String>>();
	
	private CellStyle integerStyle = null;
	private CellStyle generalStyle = null;
	private CellStyle floatStyle = null;
	private CellStyle headerStyle = null;
	private CellStyle decimalStyle = null;
	
	private Workbook wb = null;
	private Sheet sheet1 = null;
	
	

	public MergeExonMetrics(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		
		//Process command line arguments
		this.processArgs(args);
		
		//Read in dictionary Files	
		this.readDictionaryFiles();
		
		//Create the excel file
		//mergeFiles();
		createMetricsExcel();
		
		for (File m: metricsList) {
			m.delete();
		}
		
		System.out.println("Done!");
		
		
	}
	
	private void readDictionaryFiles() {
		for (File mf: this.metricsList) {
			HashMap<String,String> fileData = new HashMap<String,String>();
			
			try {
				BufferedReader br = new BufferedReader(new FileReader(mf));
				String line = null;
				
				while ((line = br.readLine()) != null) {
					String[] items = line.split("\t");
					if (fileData.containsKey(items[0])) {
						System.out.println("Already encountered this key!!!. " + items[0]);
						System.exit(1);
					} else {
						fileData.put(items[0], items[1]);
					}
				}
				
				br.close();
				lookupTables.add(fileData);
				
			} catch (IOException ioex) {
				System.out.println("Error reading dictionary file, exiting");
				ioex.printStackTrace();
				System.exit(1);
			}
		}
	}
	
	private void createMetricsExcel() {
		wb = new XSSFWorkbook();
		sheet1 = wb.createSheet("Alignment Metrics");
		
		for (int i=0;i<metricsList.size()*5 + 1;i++) {
			sheet1.setColumnWidth(i, 255*12);
		}
		
		//RowIndex
		int rIdx = 0;
	
		/*******************************************
		 * Create cell Styles
		 ******************************************/
		
		//Create styles
		DataFormat format = wb.createDataFormat();
		
		integerStyle = wb.createCellStyle();
		integerStyle.setDataFormat(format.getFormat("#,##0"));
		
		floatStyle = wb.createCellStyle();
		floatStyle.setDataFormat((short)10);
		
		generalStyle = wb.createCellStyle();
		generalStyle.setDataFormat((short)1);
		
		decimalStyle = wb.createCellStyle();
		decimalStyle.setDataFormat(format.getFormat("0.0000"));
		
		headerStyle = wb.createCellStyle();
		headerStyle.setDataFormat((short)1);
		Font hsFont = wb.createFont();
		hsFont.setBoldweight(Font.BOLDWEIGHT_BOLD);
		headerStyle.setFont(hsFont);
		
		
		//Slap borders on those styles!
		CellStyle[] csa = {integerStyle,floatStyle,generalStyle,headerStyle,decimalStyle};
		for (CellStyle cs: csa) {
			cs.setBorderBottom(CellStyle.BORDER_THIN);
			cs.setBorderTop(CellStyle.BORDER_THIN);
			cs.setBorderRight(CellStyle.BORDER_THIN);
			cs.setBorderLeft(CellStyle.BORDER_THIN);
			cs.setTopBorderColor(IndexedColors.BLACK.getIndex());
			cs.setBottomBorderColor(IndexedColors.BLACK.getIndex());
			cs.setLeftBorderColor(IndexedColors.BLACK.getIndex());
			cs.setRightBorderColor(IndexedColors.BLACK.getIndex());
		}
		
		/*******************************************
		 * Create Sheet header
		 ******************************************/
		
		//Create the header
		CellStyle styleMain = wb.createCellStyle();
		Font fontMain = wb.createFont();
		fontMain.setBoldweight(Font.BOLDWEIGHT_BOLD);
		fontMain.setFontHeightInPoints((short)24);
		styleMain.setFont(fontMain);
		
		Row headerRow = sheet1.createRow(rIdx++);
		Cell headerCell = headerRow.createCell(0);
		headerCell.setCellValue("Sample Metrics: " + this.prefix);
		headerCell.setCellStyle(styleMain);
		sheet1.addMergedRegion(new CellRangeAddress(0,0,0,2));
		
		
		//Space between header and section 1
		rIdx++;
		
		/*******************************************
		 * Create refcounts Section
		 ******************************************/
		
		//Section header
		CellStyle styleHeader = wb.createCellStyle();
		Font fontHeader = wb.createFont();
		fontHeader.setBoldweight(Font.BOLDWEIGHT_BOLD);
		fontHeader.setFontHeightInPoints((short)18);
		styleHeader.setFont(fontHeader);
		
		Row refCountTitleR = sheet1.createRow(rIdx++);
		Cell refCountTitleC = refCountTitleR.createCell(0);
		refCountTitleC.setCellValue("Reference Counts");
		refCountTitleC.setCellStyle(styleHeader);
		sheet1.addMergedRegion(new CellRangeAddress(rIdx-1,rIdx-1,0,2));
		
		//Reference counts table 1a
		
		String[] metricNames = {"Metric","Total Reads","Aligned Reads","Standard Reads","Non-Standard Reads","PhiX Reads","Adapter Reads"};
		String[] keyNames = {"sampleName","totalReads","alignReads","standardReads","nonStandardReads","phiXReads","adapterReads"};
		String[] typeList = {"header","integer","integer","integer","integer","integer","integer"};
		
		rIdx = createTable(rIdx,metricNames,keyNames,typeList,"Table 1a: Read Counts for Each Reference Sequence",lookupTables);
		
		//Reference counts table 1b
		
		String[] mn2 = {"Metric","Aligned Reads","Standard Reads","Non-Standard Reads","PhiX Reads","Adapter Reads"};
		String[] kn2 = {"sampleName","alignReadsP","standardReadsP","nonStandardReadsP","phiXReadsP","adapterReadsP"};
		String[] tl2 = {"header","float","float","float","float","float"};
		
		rIdx = createTable(rIdx,mn2,kn2,tl2,"Table 1b: Read Percentages for Each Reference Sequence",lookupTables);
		
		/*******************************************
		 * Create Duplicates section
		 ******************************************/
		
		//Section header
		Row dupTitleRow = sheet1.createRow(rIdx++);
		Cell dupTitleCell = dupTitleRow.createCell(0);
		dupTitleCell.setCellValue("Duplicates");
		dupTitleCell.setCellStyle(styleHeader);
		sheet1.addMergedRegion(new CellRangeAddress(rIdx-1,rIdx-1,0,2));
		
		//Duplicates table 2a
		
		String[] mn3 = {"Metric","Duplicate Unpaired Reads","Duplicate Paired Reads"};
		String[] kn3 = {"sampleName","dupUnpaired","dupPaired"};
		String[] tl3 = {"header","integer","integer"};
		
		rIdx = createTable(rIdx,mn3,kn3,tl3,"Table 2a: Duplicate Read Counts",lookupTables);
		
		//Duplicates table 2b
		
		String[] mn4 = {"Metric","Duplicate Unpaired Reads","Duplicate Paired Reads","Percent Duplications"};
		String[] kn4 = {"sampleName","dupUnpairedP","dupPairedP","percentDup"};
		String[] tl4 = {"header","float","float","float"};
		
		rIdx = createTable(rIdx,mn4,kn4,tl4,"Table 2b: Duplicate Reads Percentages",lookupTables);
		
		/*******************************************
		 * Create 'error rate' section
		 ******************************************/
		Row errorTitleRow = sheet1.createRow(rIdx++);
		Cell errorTitleCell = errorTitleRow.createCell(0);
		errorTitleCell.setCellValue("Error Rate");
		errorTitleCell.setCellStyle(styleHeader);
		sheet1.addMergedRegion(new CellRangeAddress(rIdx-1,rIdx-1,0,2));
		
		//Space!!
		rIdx++;
		
		//Insert images
		
		rIdx =  imageDriver("error",rIdx);
		
		/*******************************************
		 * Create 'Alignment' section
		 ******************************************/
		
		Row alignmentTitleRow = sheet1.createRow(rIdx++);
		Cell alignmentTitleCell = alignmentTitleRow.createCell(0);
		alignmentTitleCell.setCellValue("Alignment");
		alignmentTitleCell.setCellStyle(styleHeader);
		sheet1.addMergedRegion(new CellRangeAddress(rIdx-1,rIdx-1,0,2));
		
		//Create new lookup tables with read1 and read2 split out as separate samples
		String[] kn5a = {"sampleName","Read1","aTotalReads1","aAlignedReads1","aAlignedReadsP1","aPairedReads1","aPairedReadsP1","aStrand1","aError1"};
		String[] kn5b = {"sampleName","Read2","aTotalReads2","aAlignedReads1","aAlignedReadsP2","aPairedReads2","aPairedReadsP2","aStrand2","aError2"};
		String[] kn5c = {"sampleName","Read","totalReads","alignedReads","alignedPercent","pairedReads","pairedPercent","strand","error"};
		ArrayList<HashMap<String,String>> lk2 = expandDictionaries(lookupTables,kn5a,kn5b,kn5c);
		
		//Create read1/read2 alignment table
		String[] mn5 = {"","Metric","Total Reads","Aligned Reads","%Aligned","Paired","%Paired","Strand Balance","Error Rate"};
		String[] tl5 = {"header","header","integer","integer","float","integer","float","float","decimal"};
		
		rIdx = createTable(rIdx,mn5,kn5c,tl5,"Table 3a: Individual Read Alignment Statistics for Standard hg19 Chromosomes After Recalibration, Realignment and De-duplication",lk2);
		
		String[] kn6 = {"sampleName","aTotalReadsC","aAlignedReadsC","aAlignedReadsPC","aPairedReadsC",
				"aPairedReadsPC","aErrorC","aNonProper","aNonProperP","aSingleton","aSingleTonP"};
		String[] mn6 = {"Metric","Total Reads","Aligned Reads","%Aligned","Paired","%Paired",
				"Error Rate","Non-Proper Pairs","%Non-Proper Pairs","Singletons","%Singletons"};
		String[] tl6 = {"header","integer","integer","float","integer","float","decimal","integer","float","integer","float"};
		
		rIdx = createTable(rIdx,mn6,kn6,tl6,"Table 3b: Combined Read Alignment Statistics for Standard hg19 Chromosomes After Recalibration, Realignment and De-duplication",lookupTables);
		
		
		/*******************************************
		 * Create Insert Size section
		 ******************************************/
		Row insertTitleRow = sheet1.createRow(rIdx++);
		Cell insertTitleCell = insertTitleRow.createCell(0);
		insertTitleCell.setCellValue("Insert Size");
		insertTitleCell.setCellStyle(styleHeader);
		sheet1.addMergedRegion(new CellRangeAddress(rIdx-1,rIdx-1,0,2));
		
		//Space!!
		rIdx++;
		
		//Insert images
		
		rIdx =  imageDriver("insert",rIdx);
		
		/*******************************************
		 * Create Depth of Coverage section
		 ******************************************/
		Row docTitleRow = sheet1.createRow(rIdx++);
		Cell docTitleCell = docTitleRow.createCell(0);
		docTitleCell.setCellValue("Depth of Coverage");
		docTitleCell.setCellStyle(styleHeader);
		sheet1.addMergedRegion(new CellRangeAddress(rIdx-1,rIdx-1,0,2));
		
		//Space!!
		rIdx++;
		
		//Insert images
		
		rIdx =  imageDriver("coverage2",rIdx);
		
		

		
		try {
			FileOutputStream fileOut = new FileOutputStream(this.prefix + ".xlsx");
			wb.write(fileOut);
			fileOut.close();
			
		} catch (FileNotFoundException fnfe) {
			System.out.println("Could not find the output file, which makes no sense");
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Could not write to the output file");
			System.exit(1);
		}
		
	}
	
	private ArrayList<HashMap<String,String>> expandDictionaries(ArrayList<HashMap<String,String>> luTables, String[] set1, String[] set2, String[] newKeys) {
		ArrayList<HashMap<String,String>> nluTable = new ArrayList<HashMap<String,String>>();
		for (int j=0;j<luTables.size();j++) {
			HashMap<String,String> nHash1 = new HashMap<String,String>();
			HashMap<String,String> nHash2 = new HashMap<String,String>();
			
			for (int i=0; i<set1.length;i++) {
				
				if (luTables.get(j).containsKey(set1[i])) {
					nHash1.put(newKeys[i], luTables.get(j).get(set1[i]));
				} else {
					nHash1.put(newKeys[i], set1[i]);
				}
				
				if (luTables.get(j).containsKey(set2[i])) {
					nHash2.put(newKeys[i], luTables.get(j).get(set2[i]));
				} else {
					nHash2.put(newKeys[i], set2[i]);
				}
		
			}
			nluTable.add(nHash1);
			nluTable.add(nHash2);
		}
			
		return nluTable;
		
	}
	
	private int imageDriver(String imageName,int rIdx) {
		int cIdx = 0;
		int maxImagesRow = 3;
		int currImagesRow = 0;
		
		for (int i=0; i<fileList.size();i++) {
			if (currImagesRow == maxImagesRow ) {
				cIdx = 0;
				rIdx += 21;
				currImagesRow = 0;
			}
			
			createImages(cIdx, rIdx, cIdx+5, rIdx+20,fileList.get(i).get(imageName)); 
			
			cIdx += 6;
			currImagesRow++;

		}
		
		rIdx+=21;
		return rIdx;
		
	}
	
	private void createImages(int startCol, int startRow,int endCol,int endRow, File file) {
		try {
			InputStream is = new FileInputStream(file);
			byte[] bytes = IOUtils.toByteArray(is);
			int pictureIdx = wb.addPicture(bytes, Workbook.PICTURE_TYPE_PNG);
			
			
			Drawing drawing = sheet1.createDrawingPatriarch();
			ClientAnchor anchor = drawing.createAnchor(0, 0, 0, 0, startCol, startRow, endCol, endRow); 
			
			Picture pict = drawing.createPicture(anchor, pictureIdx);
			//pict.resize();
			
			
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException ioex) {
			ioex.printStackTrace();
		}
		
		
	}
	
	private int createTable(int rIdx,String[] metricNames, String[] keyNames, String[] typeList, String header, ArrayList<HashMap<String,String>> luTables) {
		
		//Space between section header and data
		rIdx++;
		
		//Create table caption row, fill rest in later
		int refCountCaptionIdx = rIdx;
		Row refCountCaptionR = sheet1.createRow(rIdx++);
		rIdx++;
		
		//Metric names
		HashMap<String,Row> table1a = new HashMap<String,Row>();
		HashMap<String,String> typeDict = new HashMap<String,String>();
		for (int i=0; i<metricNames.length; i++) {
			String mn = metricNames[i];
			String key = keyNames[i];
			String type = typeList[i];
			Row tempRow = sheet1.createRow(rIdx++);
			table1a.put(key, tempRow);
			typeDict.put(key, type);
			Cell tempCell = tempRow.createCell(0);
			tempCell.setCellValue(mn);
			tempCell.setCellStyle(headerStyle);
		}
		
		int cIdx = 1;
		for (HashMap<String,String> lt: luTables) {
			for (String key: keyNames) {
				Cell tc = table1a.get(key).createCell(cIdx);
				
				if (typeDict.get(key).equals("general")) {
					tc.setCellValue(lt.get(key));
					tc.setCellStyle(generalStyle);
				} else if (typeDict.get(key).equals("integer")) {
					tc.setCellValue(Integer.parseInt(lt.get(key)));
					tc.setCellStyle(integerStyle);
				} else if (typeDict.get(key).equals("header")) {
					tc.setCellValue(lt.get(key));
					tc.setCellStyle(headerStyle);
				} else if (typeDict.get(key).equals("float")) {
					tc.setCellValue(Float.parseFloat(lt.get(key))/100);
					tc.setCellStyle(floatStyle);
				} else if (typeDict.get(key).equals("decimal")) {
					tc.setCellValue(Float.parseFloat(lt.get(key)));
					tc.setCellStyle(decimalStyle);
				}
				
			}
			
			cIdx++;
		}
		
		//Create table caption
		Cell captionCell = refCountCaptionR.createCell(0);
		captionCell.setCellValue(header);
		sheet1.addMergedRegion(new CellRangeAddress(refCountCaptionIdx,refCountCaptionIdx,0,cIdx-1));
		
		//Space at the end of the table
		rIdx++;
		rIdx++;
		
		return rIdx;
	}
	
	private void mergeFiles() {
		try {
					
			ArrayList<String> sectionNames = new ArrayList<String>();
			ArrayList<ArrayList<ArrayList<String>>> sectionData = new ArrayList<ArrayList<ArrayList<String>>>();
			
			Pattern sectionStart = Pattern.compile("<section id='(.+?)'>");
			Pattern sectionEnd = Pattern.compile("</section>");
			Pattern headingPatt = Pattern.compile("<h2>.+?</h2>");
			
			for (int i=0; i<metricsList.size();i++) {
				BufferedReader br = new BufferedReader(new FileReader(metricsList.get(i)));
				String line = null;
				boolean slurp = false;
			
				ArrayList<ArrayList<String>> sample = new ArrayList<ArrayList<String>>();
				ArrayList<String> section = null;
				
				while ((line = br.readLine()) != null) {
					Matcher matcherStart = sectionStart.matcher(line);
					Matcher matcherEnd = sectionEnd.matcher(line);
					Matcher headingMatch = headingPatt.matcher(line);
					if (matcherStart.matches()) {
						slurp = true;
						if (i==0) {
							
							
							sectionNames.add(matcherStart.group(1));
						}
						section = new ArrayList<String>();
					} else if (matcherEnd.matches()) {
						slurp = false;
						sample.add(section);
						
					} else if  (headingMatch.matches()) {
						continue;
					} else if (slurp) {
						section.add(line);
					} 
				}
				
				sectionData.add(sample);
				br.close();
			}
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.metricsDir,prefix + ".metrics.html")));
			
			bw.write("<!DOCTYPE html>\n");
			bw.write("<html>\n");
			bw.write("<head>\n");
			bw.write("<title>Sample Metrics: " + prefix + "</title>\n");
			bw.write("<meta name='author' content='USeq ParseExonMetrics'>\n");
			bw.write("</head>\n");
			bw.write("<body>\n");
			bw.write("<h1>Sample Metrics: " + this.prefix + "</h1>\n");
			
			for (int i=0; i<sectionNames.size(); i++) {
				bw.write("<section id='" + sectionNames.get(i) + "'>\n");
				bw.write("<h2>" + sectionNames.get(i) + "</h2>\n");
				for (ArrayList<ArrayList<String>> sample: sectionData) {
					for (String line: sample.get(i)) {
						bw.write(line);
					}
					bw.write("<br />\n");
				}
				bw.write("</section>\n");
			}
			
			
			bw.write("</body>\n");
			bw.write("</html>\n");
			bw.close();
		} catch (IOException ioex) {
			System.out.println("Problem reading/writing metrics files");
			ioex.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static void main(String[] args) {
		new MergeExonMetrics(args);
	}
	
	private void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          MergeExonMetrics : June 2013                              **\n" +
				"**************************************************************************************\n" +
				"This app simply merges the output from several metrics html files.\n\n\n" +

				"Required:\n"+
				"-f Directory containing metrics dictionary files and a image directory\n"+
				"-o Name of the combined metrics file\n\n" +

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/MergeExonMetrics -f metrics -o 9908_metrics \n" +
		        "**************************************************************************************\n");

	}
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
	
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': this.metricsDir = new File(args[++i]); break;
					case 'o': this.prefix = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (this.prefix == null) {
			System.out.println("Must specify an output name");
			System.exit(1);
		}
		
		if (this.metricsDir == null) {
			System.out.println("Must specify and metrics directory");
			System.exit(1);
		}
		
		if (!metricsDir.exists()) {
			System.out.println("The specified metrics directory does not exist");
			System.exit(1);
		}
		
		File imageDir = new File(metricsDir,"images");
		
		if (!imageDir.exists()) {
			System.out.println("An image directory does not exist in the metrics directory folder");
			System.exit(1);
		}
		
		//Grab file names
		File[] mFiles = metricsDir.listFiles();
		
		//Convert to arrayList and sort
		ArrayList<File> mFileList = new ArrayList<File>(Arrays.asList(mFiles));
		Collections.sort(mFileList);
		
		ArrayList<String> names = new ArrayList<String>();
		Pattern p = Pattern.compile("(.+?)\\.dict\\.txt");
		metricsList = new ArrayList<File>();
		
		
		//Grab all metrics files in the directory and store their names
		for (File file: mFileList) {
			Matcher m = p.matcher(file.getName());
			if (m.matches()) {
				this.metricsList.add(file);
				names.add(m.group(1));
			}
		}
		
		//Make sure html files are found
		if (names.size() == 0) {
			System.out.println("No files matching *dict.txt were found");
			System.exit(1);
		}
		
		
		//Create patterns for each prefix
		ArrayList<Pattern> namePatterns = new ArrayList<Pattern>();
		fileList = new ArrayList<HashMap<String,File>>();
		
		
		for (int i=0;i<names.size();i++) {
			namePatterns.add(Pattern.compile(names.get(i) + "_(.+?)\\.png"));
			fileList.add(new HashMap<String,File>());
		}
		
		
		//Check the images directory to see if there are matching image files
		//Check that the number of images files is equal for all samples. 
		//Don't want to specify the exact image name in case metrics change in the future.
		File[] iFiles = imageDir.listFiles();
		
		
		for (File file: iFiles) {
			for (int i=0;i<namePatterns.size();i++) {
				Pattern np = namePatterns.get(i);
				Matcher m = np.matcher(file.getName());
				if (m.matches()) {
					fileList.get(i).put(m.group(1), file);
				}
			}
		}
		
		Integer size = null;
		
		for (int i=0;i<fileList.size();i++) {
			if (fileList.get(i).size() == 0) {
				System.out.println("No image files matching prefix: " + names.get(i));
				System.exit(1);
			}
			
			if (size == null) {
				size = fileList.get(i).size();
			} else if (size != fileList.get(i).size()) {
				System.out.println("Number of image files for each prefix doesn't match");
				System.exit(1);
			}
			
		}
		
		//Sort the filelist
		

	}
		

}
