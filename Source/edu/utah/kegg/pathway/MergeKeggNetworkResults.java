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

public class MergeKeggNetworkResults {

	//user fields
	private File[] networkFilesToMerge = null;
	private File spreadsheetFile = null;
	private double maxAdjPVal = 0.1;
	
	//internal
	private TreeMap<String, HashMap<String, KeggNetworkToMerge>> netInfo = new TreeMap<String, HashMap<String,KeggNetworkToMerge>>();
	private ArrayList<String> fileNames = new ArrayList<String>();
	private CreationHelper createHelper = null;
	private CellStyle boldCenterStyle = null;
	private XSSFCellStyle lightRedStyle = null;
	private CellStyle linkStyle = null;
	private ArrayList<String> header = null;
	private Pattern quoteCommaQuote = Pattern.compile("\",\"");
	private Pattern spaceColon = Pattern.compile(" : ");
	
	//constructor for cmd line
	public MergeKeggNetworkResults(String[] args) {
		try {

			processArgs(args);
			
			parseNetworkFiles();
			
			//printMergedNetworksAll();
			//printMergedNetworksAdjPVals();
			//printMergedNetworksGenes();
			//printMergedNetworksPathways();
			printResultsToXlsxFile();
			
		} catch (Exception e) {
			e.printStackTrace();
			IO.el("ERROR running the MergeKeggNetworkResults");
			System.exit(1);
		}
	}
	
	private void printResultsToXlsxFile() throws Exception {
		IO.pl("\nPrinting xlsx spreadsheet results...");
		
		Workbook workbook = new XSSFWorkbook();
		createHelper = workbook.getCreationHelper();
		
		// Create cell styles 
		// Bold center
        Font boldFont = workbook.createFont();
        boldFont.setBoldweight(boldFont.BOLDWEIGHT_BOLD);
        boldCenterStyle = workbook.createCellStyle();
        boldCenterStyle.setFont(boldFont);
        boldCenterStyle.setAlignment(CellStyle.ALIGN_CENTER);
        boldCenterStyle.setWrapText(true);
        
        // Light red 
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
        header.add("Network Desc HyperLink");
        header.addAll(fileNames);
        
        addPValueSheet(workbook);
        
        addGenesSheet(workbook);
        
        addPathwaysSheet(workbook);
        
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

		//for each network
		for (String nd: netInfo.keySet()) {
			Row nRow = sheet.createRow(counter++);

			//add placeholder cell for sig count
			Cell numSigCell = nRow.createCell(0);

			//add the network hyperlink
			String[] linkUrlName = parseLinkNameUrl(nd);
			Cell netHyper = nRow.createCell(1);
			Hyperlink link = createHelper.createHyperlink(Hyperlink.LINK_URL);
			link.setAddress(linkUrlName[0]);
			netHyper.setHyperlink(link);
			netHyper.setCellValue(linkUrlName[1]);
			netHyper.setCellStyle(linkStyle);

			//for each fileName
			int numSig = 0;
			HashMap<String, KeggNetworkToMerge> nameData = netInfo.get(nd);
			int cellCounter = 2;
			for (String fileName: fileNames) {
				Cell info = nRow.createCell(cellCounter++);
				//any network info
				KeggNetworkToMerge n = nameData.get(fileName);
				if (n!=null) {
					info.setCellValue(n.adjPVal);
					if (n.adjPVal<= maxAdjPVal) {
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

		//for each network
		for (String nd: netInfo.keySet()) {
			Row nRow = sheet.createRow(counter++);

			//add placeholder cell for sig count
			Cell numSigCell = nRow.createCell(0);

			//add the network hyperlink
			String[] linkUrlName = parseLinkNameUrl(nd);
			Cell netHyper = nRow.createCell(1);
			Hyperlink link = createHelper.createHyperlink(Hyperlink.LINK_URL);
			link.setAddress(linkUrlName[0]);
			netHyper.setHyperlink(link);
			netHyper.setCellValue(linkUrlName[1]);
			netHyper.setCellStyle(linkStyle);

			//for each fileName
			int numSig = 0;
			HashMap<String, KeggNetworkToMerge> nameData = netInfo.get(nd);
			int cellCounter = 2;
			for (String fileName: fileNames) {
				Cell info = nRow.createCell(cellCounter++);
				//any network info
				KeggNetworkToMerge n = nameData.get(fileName);
				if (n!=null) {
					if (n.adjPVal<= maxAdjPVal) {
						info.setCellValue(n.intersectingGenes);
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

	private void addPathwaysSheet(Workbook workbook) throws Exception {
		Sheet sheet = workbook.createSheet("Pathways");

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

		//for each network
		for (String nd: netInfo.keySet()) {
			Row nRow = sheet.createRow(counter++);

			//add placeholder cell for sig count
			Cell numSigCell = nRow.createCell(0);

			//add the network hyperlink
			String[] linkUrlName = parseLinkNameUrl(nd);
			Cell netHyper = nRow.createCell(1);
			Hyperlink link = createHelper.createHyperlink(Hyperlink.LINK_URL);
			link.setAddress(linkUrlName[0]);
			netHyper.setHyperlink(link);
			netHyper.setCellValue(linkUrlName[1]);
			netHyper.setCellStyle(linkStyle);

			//for each fileName
			int numSig = 0;
			HashMap<String, KeggNetworkToMerge> nameData = netInfo.get(nd);
			int cellCounter = 2;
			for (String fileName: fileNames) {
				Cell info = nRow.createCell(cellCounter++);
				//any network info?
				KeggNetworkToMerge n = nameData.get(fileName);
				if (n!=null) {
					if (n.adjPVal<= maxAdjPVal) {
						//info.setCellValue(n.intersectingGenes);
						ArrayList<String> pathways = n.pathwayLinks;
						
						//just one?
						if (pathways.size() == 1) {
							String[] urlName = parseLinkNameUrl(pathways.get(0));
							Hyperlink lnk = createHelper.createHyperlink(Hyperlink.LINK_URL);
							lnk.setAddress(urlName[0]);
							info.setHyperlink(lnk);
							info.setCellValue(urlName[1]);
							info.setCellStyle(linkStyle);
						}
						//greater than one
						else if (pathways.size()>1) {
							StringBuilder sb = new StringBuilder();
							boolean addComma = false;
							for (String p: pathways) {
								if (addComma) sb.append(", ");
								String[] urlName = parseLinkNameUrl(p);
								sb.append(urlName[1]);
								addComma = true;
							}
							info.setCellValue(sb.toString());
						}
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

	private String[] parseLinkNameUrl(String hyperLink) throws Exception {
//What about the other type		
		// =HYPERLINK("https://www.kegg.jp/entry/N00009","N00009: TRK fusion kinase to RAS-ERK signaling pathway")
		//Apoptosis : https://www.kegg.jp/kegg-bin/show_pathway?map=hsa04210&multi_query=317%20%23FFFFAC%0A572%20%23FFFFAC%0A578%20%23FFFFAC%0A581%20%23FFFFAC%0A27113%20%23D6B4FC%0A596%20%23FFD580%0A10018%20%23FFFFAC%0A836%20%23FFFFAC%0A840%20%23FFFFAC%0A54205%20%23FFFFAC&network=N00098
		String[] split = null;
		if (hyperLink.startsWith("=")) {
			split = quoteCommaQuote.split(hyperLink);
			if (split.length != 2) throw new Exception("ERROR: failed to split url in 2 on comma quotes "+hyperLink);
			split[0] = split[0].substring(12);
			split[1] = split[1].substring(0, split[1].length()-2);
		}
		else {
			split = spaceColon.split(hyperLink);
			String name = split[0];
			String url = split[1];
			split[0] = url;
			split[1] = name;
		}
		split[1] = split[1].replace("pathway", "");
		return split;
	}

	private void printMergedNetworksAdjPVals() {
		//for each networkDesc
		IO.pl("NumSigDatasets\tNetworkDescHyperLink\t"+Misc.stringArrayListToString(fileNames, "\t"));
		for (String nd: netInfo.keySet()) {
			int numSig = 0;
			StringBuilder sb = new StringBuilder();
			sb.append(nd); 
			
			//for each fileName
			HashMap<String, KeggNetworkToMerge> nameData = netInfo.get(nd);
			for (String fileName: fileNames) {
				sb.append("\t");
				//any network info
				KeggNetworkToMerge n = nameData.get(fileName);
				if (n!=null) {
					sb.append(n.adjPVal);
					if (n.adjPVal<= maxAdjPVal) numSig++;
				}
				else sb.append("NA");
			}
			IO.pl(numSig+"\t"+sb);
		}
	}
	
	private void printMergedNetworksGenes() {
		//for each networkDesc
		IO.pl("NumSigDatasets\tNetworkDescHyperLink\t"+Misc.stringArrayListToString(fileNames, "\t"));
		for (String nd: netInfo.keySet()) {
			int numSig = 0;
			StringBuilder sb = new StringBuilder();
			sb.append(nd); 
			
			//for each fileName
			HashMap<String, KeggNetworkToMerge> nameData = netInfo.get(nd);
			for (String fileName: fileNames) {
				sb.append("\t");
				//any network info
				KeggNetworkToMerge n = nameData.get(fileName);
				if (n!=null) {
					//only print if sig
					if (n.adjPVal<= maxAdjPVal) {
						sb.append(n.intersectingGenes);
						numSig++;
					}
					else sb.append(".");
				}
				else sb.append(".");
			}
			IO.pl(numSig+"\t"+sb);
		}
	}
	
	private void printMergedNetworksPathways() throws IOException {
		//for each networkDesc
		IO.pl("NumSigDatasets\tNetworkDescHyperLink\t"+Misc.stringArrayListToString(fileNames, "\t"));
		for (String nd: netInfo.keySet()) {
			int numSig = 0;
			StringBuilder sb = new StringBuilder();
			sb.append(nd); 
			
			//for each fileName
			HashMap<String, KeggNetworkToMerge> nameData = netInfo.get(nd);
			for (String fileName: fileNames) {
				sb.append("\t");
				//any network info
				KeggNetworkToMerge n = nameData.get(fileName);
				if (n!=null) {
					//only print if sig
					if (n.adjPVal<= maxAdjPVal) {
						sb.append(n.getPathwayNames());
						numSig++;
					}
					else sb.append(".");
				}
				else sb.append(".");
			}
			IO.pl(numSig+"\t"+sb);
		}
	}

	private void printMergedNetworksAll() {
		//for each networkDesc
		for (String nd: netInfo.keySet()) {
			IO.pl("-------------------------------------");
			IO.pl(nd);
			
			//for each fileName
			HashMap<String, KeggNetworkToMerge> nameData = netInfo.get(nd);
			for (String fileName: fileNames) {
				IO.pl(fileName);
				//any network info
				KeggNetworkToMerge n = nameData.get(fileName);
				if (n!=null) {
					IO.pl("\tPVal\t"+n.adjPVal);
					IO.pl("\tGenes\t"+n.intersectingGenes);
					IO.pl("\tPath\t"+n.pathwayLinks);
				}
			}
		}
	}


	public void parseNetworkFiles() {
		IO.pl("Parsing xls files...");
		for (int i=0; i< networkFilesToMerge.length; i++) {
			String fileName = Misc.removeExtension(networkFilesToMerge[i].getName());
			fileNames.add(fileName);
			IO.pl("\t"+networkFilesToMerge[i].getName());
			Boolean isGene = null;
			String[] lines = IO.loadFile(networkFilesToMerge[i]);
			for (String line : lines) {
				//header?
				if (line.contains("Network Desc Link")) {
					//what file type is this?
					if (line.contains("ANoHits")) isGene = false;
					else if (line.contains("Int Gene Symbols")) isGene = true;
					continue;
				}
				if (line.trim().length() == 0) continue;
				
				String[] f = Misc.TAB.split(line);
				// Gene
				// NetworkName NetworkDescLink Pval AdjPval AllNetworkGenes FoundNetworkGenes SelectGenes Intersect Int/Net IntGeneSymbols IntGenesLgRtos PathwayMapLinks...
				//      0              1         2     3           4                5              6           7       8           9         10              11 12 13 ...          
				if (f.length < 10) continue;

				// Variant
				// #NetworkName(s) NetworkDescLink Pval AdjPval AHits ANoHits FracAHits AGeneHits BHits BNoHits FracBHits BGeneHits Log2(fracA/fracB) AllGeneHits PathwayMapLinksWithTopMapDescription...													
				//       0                1          2     3      4      5        6         7       8      9       10         11            12             13               14 15 16 ...
				
				HashMap<String, KeggNetworkToMerge> al = netInfo.get(f[1]);
				if (al == null) {
					al = new HashMap<String, KeggNetworkToMerge>();
					netInfo.put(f[1], al);
				}
				String genesInvolved = null;
				int pathwayIndex = -1;
				if (isGene) {
					genesInvolved = f[9];
					pathwayIndex = 11;
				}
				else {
					genesInvolved = f[13];
					pathwayIndex = 14;
				}
				
				KeggNetworkToMerge m = new KeggNetworkToMerge(fileName, Double.parseDouble(f[3]), genesInvolved);
				
				for (int x=pathwayIndex; x< f.length; x++) m.pathwayLinks.add(f[x]);
				al.put(fileName, m);
			}
		}

	}

	private class KeggNetworkToMerge {
		String fileName;
		double adjPVal;
		String intersectingGenes;
		ArrayList<String> pathwayLinks = new ArrayList<String>();
		
		public KeggNetworkToMerge(String fileName, double adjPVal, String intersectingGenes) {
			this.fileName = fileName;
			this.adjPVal = adjPVal;
			this.intersectingGenes = intersectingGenes;
		}

		
		public String getPathwayNames() throws IOException {
			StringBuilder sb = new StringBuilder();
			//=HYPERLINK("https://www.kegg.jp/kegg-bin/show_pathway?map=hsa04062&multi_query=4773%20%2387CEEB%0A7852%20%23FFB6C1&network=N00401","Chemokine signaling pathway")
			//Apoptosis : https://www.kegg.jp/kegg-bin/show_pathway?map=hsa04210&multi_query=317%20%23FFFFAC%0A572%20%23FFFFAC%0A578%20%23FFFFAC%0A581%20%23FFFFAC%0A27113%20%23D6B4FC%0A596%20%23FFD580%0A10018%20%23FFFFAC%0A836%20%23FFFFAC%0A840%20%23FFFFAC%0A54205%20%23FFFFAC&network=N00098
			boolean addComma = false;
			String p = null;
			for (String hyperLink : pathwayLinks) {
				if (addComma) sb.append(", ");
				if (hyperLink.startsWith("=")) {
					String[] split = quoteCommaQuote.split(hyperLink);
					if (split.length!=2) throw new IOException("ERROR: failed to split the hyperlink in two: "+hyperLink);
					p = split[1].substring(0, split[1].length()-2);
				}
				else {
					String[] split = spaceColon.split(hyperLink);
					if (split.length!=2) throw new IOException("ERROR: failed to split the hyperlink in two: "+hyperLink);
					p = split[0];
				}
				p = p.replace("pathway", "");
				sb.append(p.trim());
				addComma = true;
			}
			return sb.toString();
		}
	}

	public static void main(String[] args) throws IOException {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergeKeggNetworkResults(args);
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
					case 'n': networkFilesToMerge = IO.extractFiles(new File(args[++i]), ".xls"); break;
					case 's': spreadsheetFile = new File(args[++i]); break;
					case 'm': maxAdjPVal = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (networkFilesToMerge == null || networkFilesToMerge.length == 0) {
			Misc.printErrAndExit("\nERROR: failed to find any xxx.xls Kegg network files to merge in your -n directory.\n");
		}
		if (spreadsheetFile == null) {
			Misc.printErrAndExit("\nERROR: failed to find your xxx.xlsx output file.\n");
		}
		
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Merge Kegg Network Results : Feb 2026                      **\n" +
				"**************************************************************************************\n" +
				"MKNR merges gene and variant network spreadsheet xxx.xls (not xlsx) files from the\n"+
				"USeq Kegg analysis applications into a three tab xlsx spreadsheet. Useful for \n"+
				"comparing between multiple pathway analysis.\n"+
				
				"\nApp Parameters:\n\n"+
				
				"-n Directory containing gene and or variant network xxx.xls files from the USeq Kegg\n"+
				"      analysis tools to merge.\n"+
				"-s Path to an xxx.xlsx file for saving the merged results.\n"+
				"-m Maximum adjusted p-value to consider significant, defaults to 0.1\n"+
				
				"\nExample:\n\n"+ 
				"java -Xmx1G -jar pathTo/USeq/Apps/MergeKeggNetworkResults -m 0.15 -s mergedNets.xlsx\n"+
				"   -n USeqKeggNetworkFiles/\n"+
				
				"\n**************************************************************************************\n");

	}

}
