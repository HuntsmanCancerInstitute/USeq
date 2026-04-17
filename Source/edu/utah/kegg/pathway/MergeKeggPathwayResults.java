package edu.utah.kegg.pathway;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.CreationHelper;
import org.apache.poi.ss.usermodel.FillPatternType;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Hyperlink;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFColor;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import util.gen.*;

public class MergeKeggPathwayResults {

	//user fields
	private File[] pathwayFilesToMerge = null;
	private File spreadsheetFile = null;
	private double maxPathwayAdjPVal = 0.1;
	
	//internal
	//pathwayName, hash of fileName and KeggPathwayToMerge obj with that datasets info
	private TreeMap<String, HashMap<String, KeggPathwayToMerge>> pathwayInfo = new TreeMap<String, HashMap<String,KeggPathwayToMerge>>();
	private HashMap<String, String> pathwayNakedUrl = new HashMap<String, String>();
	private ArrayList<String> fileNames = new ArrayList<String>();

	private CreationHelper createHelper = null;
	private CellStyle boldCenterStyle = null;
	private XSSFCellStyle lightRedStyle = null;
	private CellStyle linkStyle = null;
	private ArrayList<String> header = null;
	
	//constructor for cmd line
	public MergeKeggPathwayResults(String[] args) {
		try {

			processArgs(args);
			
			parsePathwayFiles();
			
			printResultsToXlsxFile();
			
		} catch (Exception e) {
			e.printStackTrace();
			IO.el("ERROR running the MergeKeggPathwayResults");
			System.exit(1);
		}
	}
	
	private void printResultsToXlsxFile() throws Exception {
		IO.pl("\nPrinting xlsx spreadsheet results...");
		
		Workbook workbook = new XSSFWorkbook();
		createHelper = workbook.getCreationHelper();
		
		// Create cell styles 
		// Bold center for header
        Font boldFont = workbook.createFont();
        boldFont.setBoldweight(boldFont.BOLDWEIGHT_BOLD);
        boldCenterStyle = workbook.createCellStyle();
        boldCenterStyle.setFont(boldFont);
        boldCenterStyle.setAlignment(CellStyle.ALIGN_CENTER);
        boldCenterStyle.setWrapText(true);
        
        // Light red for sig pathway FDRs
        byte[] rgb = new byte[] {(byte)255, (byte)230, (byte)230}; 
        XSSFColor lightRed = new XSSFColor(rgb);
        lightRedStyle = (XSSFCellStyle) workbook.createCellStyle();
        lightRedStyle.setFillForegroundColor(lightRed);
        lightRedStyle.setFillPattern(FillPatternType.SOLID_FOREGROUND);
        
        // Link style
        linkStyle = workbook.createCellStyle();
        Font linkFont = workbook.createFont();
        linkFont.setUnderline(Font.U_SINGLE);
        linkFont.setColor(IndexedColors.BLUE.getIndex());
        linkStyle.setFont(linkFont);
        
        // Make header
        header = new ArrayList<String>();
        header.add("Num Sig Datasets");
        header.add("Pathway Desc HyperLink");
        header.addAll(fileNames);
        
        addPValueSheet(workbook);
        
        addLoadedPathwaySheet(workbook);
        
        addNetworksSheet(workbook);
        
        addGenesSheet(workbook);

		//save it
		FileOutputStream fileOut = new FileOutputStream(spreadsheetFile);
        workbook.write(fileOut);
        fileOut.close();
	}

	private void addPValueSheet(Workbook workbook) throws Exception {
		Sheet sheet = workbook.createSheet("AdjPValues");

		//add header
		int counter = 0;
		Row row = sheet.createRow(counter++);
		for (int i=0; i< header.size(); i++) {
			Cell cell = row.createCell(i);
			cell.setCellValue(header.get(i));
			cell.setCellStyle(boldCenterStyle);
		}
		// and freeze the header row
		sheet.createFreezePane(0, 1);

		//for each pathway
		for (String pd: pathwayInfo.keySet()) {
			Row nRow = sheet.createRow(counter++);
			
			//add placeholder cell for sig count
			Cell numSigCell = nRow.createCell(0);

			//add the pathway hyperlink
			Cell netHyper = nRow.createCell(1);
			Hyperlink link = createHelper.createHyperlink(Hyperlink.LINK_URL);
			link.setAddress(pathwayNakedUrl.get(pd));
			netHyper.setHyperlink(link);
			netHyper.setCellValue(pd);
			netHyper.setCellStyle(linkStyle);

			//for each fileName
			int numSig = 0;
			HashMap<String, KeggPathwayToMerge> nameData = pathwayInfo.get(pd);
			
//if (pd.contains("Calcium signaling pathway")) IO.pl("Here "+nameData.size());
			
			int cellCounter = 2;
			for (String fileName: fileNames) {
				Cell info = nRow.createCell(cellCounter++);
				//any pathway info
				KeggPathwayToMerge p = nameData.get(fileName);
				if (p!=null) {
					info.setCellValue(p.adjPVal);
					if (p.adjPVal<= maxPathwayAdjPVal) {
						info.setCellStyle(lightRedStyle);
						numSig++;
					}
				}
			}
			
			//set numSig
			numSigCell.setCellValue(numSig);
		}
	}
	
	private void addGenesSheet(Workbook workbook) throws Exception {
		Sheet sheet = workbook.createSheet("Genes");

		//add header
		int counter = 0;
		Row row = sheet.createRow(counter++);
		for (int i=0; i< header.size(); i++) {
			Cell cell = row.createCell(i);
			cell.setCellValue(header.get(i));
			cell.setCellStyle(boldCenterStyle);
		}
		// and freeze the header row
		sheet.createFreezePane(0, 1);

		//for each pathway
		for (String pd: pathwayInfo.keySet()) {
			Row nRow = sheet.createRow(counter++);

			//add placeholder cell for sig count
			Cell numSigCell = nRow.createCell(0);

			//add the pathway hyperlink
			Cell netHyper = nRow.createCell(1);
			Hyperlink link = createHelper.createHyperlink(Hyperlink.LINK_URL);
			link.setAddress(pathwayNakedUrl.get(pd));
			netHyper.setHyperlink(link);
			netHyper.setCellValue(pd);
			netHyper.setCellStyle(linkStyle);

			//for each fileName
			int numSig = 0;
			HashMap<String, KeggPathwayToMerge> nameData = pathwayInfo.get(pd);
			int cellCounter = 2;
			for (String fileName: fileNames) {
				Cell info = nRow.createCell(cellCounter++);
				//any pathway info
				KeggPathwayToMerge p = nameData.get(fileName);
				if (p!=null) {
					if (p.adjPVal<= maxPathwayAdjPVal) {
						info.setCellValue(p.allGenes);
						numSig++;
					}
					else info.setCellValue(".");
				}
				else info.setCellValue(".");
			}
			
			//set numSig
			numSigCell.setCellValue(numSig);
		}
	}

	private void addNetworksSheet(Workbook workbook) throws Exception {
		Sheet sheet = workbook.createSheet("Networks");

		//add header
		int counter = 0;
		Row row = sheet.createRow(counter++);
		for (int i=0; i< header.size(); i++) {
			Cell cell = row.createCell(i);
			cell.setCellValue(header.get(i));
			cell.setCellStyle(boldCenterStyle);
		}
		// and freeze the header row
		sheet.createFreezePane(0, 1);

		//for each pathway
		for (String pd: pathwayInfo.keySet()) {
			Row nRow = sheet.createRow(counter++);

			//add placeholder cell for sig count
			Cell numSigCell = nRow.createCell(0);

			//add the pathway hyperlink
			Cell netHyper = nRow.createCell(1);
			Hyperlink link = createHelper.createHyperlink(Hyperlink.LINK_URL);
			link.setAddress(pathwayNakedUrl.get(pd));
			netHyper.setHyperlink(link);
			netHyper.setCellValue(pd);
			netHyper.setCellStyle(linkStyle);

			//for each fileName
			int numSig = 0;
			HashMap<String, KeggPathwayToMerge> nameData = pathwayInfo.get(pd);
			int cellCounter = 2;
			for (String fileName: fileNames) {
				Cell info = nRow.createCell(cellCounter++);
				//any pathway info
				KeggPathwayToMerge p = nameData.get(fileName);
				if (p!=null) {
					if (p.adjPVal<= maxPathwayAdjPVal) {
						info.setCellValue(Misc.stringArrayListToString(p.networkNames, ", "));
						numSig++;
					}
					else info.setCellValue(".");
				}
				else info.setCellValue(".");
			}
			//set numSig
			numSigCell.setCellValue(numSig);
		}
	}
	
	private void addLoadedPathwaySheet(Workbook workbook) throws Exception {
		Sheet sheet = workbook.createSheet("AnnotatedPathways");

		//add header
		int counter = 0;
		Row row = sheet.createRow(counter++);
		for (int i=0; i< header.size(); i++) {
			Cell cell = row.createCell(i);
			cell.setCellValue(header.get(i));
			cell.setCellStyle(boldCenterStyle);
		}
		// and freeze the header row
		sheet.createFreezePane(0, 1);

		//for each pathway
		for (String pd: pathwayInfo.keySet()) {
			Row nRow = sheet.createRow(counter++);

			//add placeholder cell for sig count
			Cell numSigCell = nRow.createCell(0);

			//add the pathway hyperlink
			Cell netHyper = nRow.createCell(1);
			netHyper.setCellValue(pd);

			//for each fileName
			int numSig = 0;
			HashMap<String, KeggPathwayToMerge> nameData = pathwayInfo.get(pd);
			int cellCounter = 2;
			for (String fileName: fileNames) {
				Cell info = nRow.createCell(cellCounter++);
				//any pathway info
				KeggPathwayToMerge p = nameData.get(fileName);
				if (p!=null) {
					if (p.adjPVal<= maxPathwayAdjPVal) {
						Hyperlink url = createHelper.createHyperlink(Hyperlink.LINK_URL);
						url.setAddress(p.pathwayURL);
						info.setHyperlink(url);
						info.setCellValue(p.pathwayName.replace(" pathway", ""));
						info.setCellStyle(linkStyle);
						numSig++;
					}
					else info.setCellValue(".");
				}
				else info.setCellValue(".");
			}
			//set numSig
			numSigCell.setCellValue(numSig);
		}
	}

	public void parsePathwayFiles() throws Exception {
		IO.pl("Parsing pathway xls files...");
		for (int i=0; i< pathwayFilesToMerge.length; i++) {
			String fileName = Misc.removeExtension(pathwayFilesToMerge[i].getName());
			fileNames.add(fileName);
			IO.pl("\t"+pathwayFilesToMerge[i].getName());

			String[] lines = IO.loadFile(pathwayFilesToMerge[i]);
			String[] fields = null;
			for (int x=0; x< lines.length; x++) {
				//advance to next PATHWAY
				if (lines[x].startsWith("PATHWAY")) {
					fields = Misc.TAB.split(lines[x]);
					//pull or create new entry for the pathway
					String pathwayName = fields[2];
					HashMap<String, KeggPathwayToMerge> al = pathwayInfo.get(pathwayName);
					if (al == null) {
						al = new HashMap<String, KeggPathwayToMerge>();
						// =HYPERLINK("https://www.kegg.jp/entry/hsa04261","https://www.kegg.jp/entry/hsa04261")
						String[] split = MergeKeggNetworkResults.quoteCommaQuote.split(fields[3]);
						if (split.length != 2) throw new Exception("ERROR: failed to split url in 2 on comma quotes "+fields[3]);
						split[0] = split[0].substring(12);
						split[1] = split[1].substring(0, split[1].length()-2);
						pathwayInfo.put(pathwayName, al);
						pathwayNakedUrl.put(pathwayName, split[0]);
					}
					//add in new KeggPathwayToMerge obj representing the file pathway dataset
					KeggPathwayToMerge kp = new KeggPathwayToMerge(fields[2]);
					al.put(fileName, kp);
					
					
					//load in network info, KeggLink, and PathwayFDR etc to the KeggPathwayToMerge obj
					for (int k = x+1; k<lines.length; k++) {
						//skip Networks line
						if (lines[k].startsWith("Networks")) continue;
						fields = Misc.TAB.split(lines[k]);
//if (lines[k].contains("Diff Exp")) IO.pl("Here:"+fields[0]+"---"+fields[2]); 
						//last line? PathwayFDR?
						if (lines[k].startsWith("PathwayFDR")) {
							double fdr = Double.parseDouble(fields[1]);
							kp.adjPVal = fdr;
							x=k;
							break;
						}
						//AllGenes?
						else if (lines[k].startsWith("AllGenes")) kp.allGenes = fields[1];
						//KeggLink?
						else if (lines[k].startsWith("KeggLink")) kp.pathwayURL = fields[1];
						//Network from combine
						else if (fields[0].startsWith("Diff")) {
							kp.networkNames.add(fields[2].replace(" pathway", ""));
						}
						//Network from gene or variant
						else if (fields[0].startsWith("N0")) {
							kp.networkNames.add(fields[1].replace(" pathway", ""));
						}
					}
				}
			}
		}
	}

	private class KeggPathwayToMerge {
		String pathwayName;
		double adjPVal;
		String allGenes;
		String pathwayURL;
		ArrayList<String> networkNames = new ArrayList<String>();
		
		public KeggPathwayToMerge(String pathwayName) {
			this.pathwayName = pathwayName;
		}
	}

	public static void main(String[] args) throws IOException {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergeKeggPathwayResults(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': pathwayFilesToMerge = IO.extractFiles(new File(args[++i]), ".xls"); break;
					case 's': spreadsheetFile = new File(args[++i]); break;
					case 'm': maxPathwayAdjPVal = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (pathwayFilesToMerge == null || pathwayFilesToMerge.length == 0) {
			Misc.printErrAndExit("\nERROR: failed to find any xxx.xls Kegg network files to merge in your -n directory.\n");
		}
		if (spreadsheetFile == null) {
			Misc.printErrAndExit("\nERROR: failed to find your xxx.xlsx output file.\n");
		}
		if (spreadsheetFile.getName().endsWith(".xlsx") == false) {
			Misc.printErrAndExit("\nERROR: your output spreadsheed file doesn't end in  '.xlsx'\n");
		}
		
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Merge Kegg Pathway Results : Feb 2026                      **\n" +
				"**************************************************************************************\n" +
				"MKPR merges USeq KEGG 'pathway' spreadsheet xxx.xls (not xlsx) files from the\n"+
				"USeq Kegg analysis applications into a four tab xlsx spreadsheet. Useful for \n"+
				"comparing between multiple pathway analysis.\n"+
				
				"\nNOTE: If have multiple KeggGeneAndVariantPathwayAnalyzer datasets, ONLY run this\n"+
				"tool on the combineGeneVariantPathwaysMinGenXXXMaxFdrXXX.xls files. Don't put in the\n"+
				"intermediate genePathwayXXX.xls or variantPathwayXXX.xls files.\n"+
				
				"\nApp Parameters:\n\n"+
				
				"-p Directory containing the combine or gene and/or variant pathway xxx.xls files from\n"+
				"      the USeq Kegg analysis tools to merge. See the NOTE above.\n"+
				"-s Path to an xxx.xlsx file for saving the merged results.\n"+
				"-m Maximum pathway FDR / adjusted p-value to consider significant, defaults to 0.1\n"+
				
				"\nExample:\n\n"+ 
				"java -Xmx1G -jar pathTo/USeq/Apps/MergeKeggPathwayResults -m 0.15 -s mergedPaths.xlsx\n"+
				"   -p USeqKeggPathwayFiles/\n"+
				
				"\n**************************************************************************************\n");

	}

}
