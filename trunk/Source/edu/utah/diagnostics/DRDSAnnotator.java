package edu.utah.diagnostics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.poi.openxml4j.exceptions.InvalidFormatException;
import org.apache.poi.openxml4j.opc.OPCPackage;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import util.gen.Misc;

public class DRDSAnnotator {
	private File inputFile = null;
	private File outputFile = null;
	private File annotationFile = null;
	private HashMap<String,String[]> annotations = null;
	private String[] header = null;
	private int entries = 0;
	private Workbook wb = null;
	

	public DRDSAnnotator(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		parseArgs(args);
		
		System.out.println("Parsing biomart annotations");
		this.parseBiomart();
		
		System.out.println("Adding annotations to file");
		this.addAnnotations();
		
		System.out.println("Writing workbook to file");
		this.writeWorkBook();
		
		System.out.println("Finished");
	}
	
	public static void main(String[] args) {
		new DRDSAnnotator(args);

	}
	
	private void addAnnotations() {
		try {
			
			//Open up xlsx file
			wb = WorkbookFactory.create(new File(this.inputFile.getAbsolutePath()));	
			
			//Grab appropriate sheet and bail if it doesn't exist
			Sheet oldSheet = wb.getSheet("Analyzed Genes");
			if (oldSheet == null) {
				System.out.println("Could not find sheet 'Analyzed Genes', are you sure you are using the correct input file?  Exiting: "
						+ this.inputFile.toString());
				System.exit(1);
			}
			
			//Create a new sheet.  It looks like POI doesn't support inserting columns, so we need to recreate the sheet :(
			Sheet newSheet = wb.createSheet("Analyzed Genes Annotated");
			
			
			//CellStyle container
			CellStyle[] csList = new CellStyle[oldSheet.getRow(0).getPhysicalNumberOfCells()];
			
			
			//Iterate through rows
			for (int r = 0; r<oldSheet.getPhysicalNumberOfRows(); r++) {
				newSheet.createRow(r);
				int spacer = 0;
				
				if (r % 1000 == 0 && r != 0) {
					System.out.println(String.format("Annotated %d lines.",r));
				}
				
				for (int c=0;c < oldSheet.getRow(c).getPhysicalNumberOfCells(); c++) {
					
					//Get old cell
					Cell oldCell = oldSheet.getRow(r).getCell(c);
					
					//Simply duplicate old cell without fanfare
					Cell newCell = newSheet.getRow(r).createCell(c + spacer);
					try {
						newCell.setCellValue(oldCell.getStringCellValue());
					} catch (IllegalStateException ise) {
						newCell.setCellValue(oldCell.getNumericCellValue());
					}
					
					//Grab existing style, greatly speeds up the annotation
					CellStyle newCellStyle;
					if (r == 1 || r == 0) {
						newCellStyle = wb.createCellStyle();
						newCellStyle.cloneStyleFrom(oldCell.getCellStyle());
							csList[c] = newCellStyle;
					} else {
						newCellStyle = csList[c];
					}
					
					newCell.setCellStyle(newCellStyle);
					
					//If first column, run matching
					if (c == 0) {
						if (r == 0) {
							//First row is header, so don't match, just add headers
							for (String h: this.header) {
								spacer++; //
								Cell ns = newSheet.getRow(r).createCell(c + spacer);
								ns.setCellValue(h); //Add header information
								ns.setCellStyle(newCellStyle); //just use existing cell style
							}	
						} else {
							//Subsequent rows need to be matched
							
							//Grab annotations if there is match
							String[] annList = null;
							if (this.annotations.containsKey(oldCell.getStringCellValue())) {
						    	annList = this.annotations.get(oldCell.getStringCellValue());
						    }
							
							for (int i=0;i<entries;i++) {
								spacer++;
								
								//set new value if there is match
								String newValue = "NA";
								if (annList != null) {
									newValue = annList[i];
								}
							    
								//Create and add new cell
								Cell ns = newSheet.getRow(r).createCell(c + spacer);
								ns.setCellValue(newValue);
							} //end for
						} //end if else
					} //end if c
				} //end for c
			} //end for r
			
			//Duplicate freeze panes
			newSheet.createFreezePane(0,1);
			
			//Remove old sheet, rename new sheet and move into position
			wb.removeSheetAt(wb.getSheetIndex("Analyzed Genes"));
			wb.setSheetName(wb.getSheetIndex("Analyzed Genes Annotated"), "Analyzed Genes");
			wb.setSheetOrder("Analyzed Genes", 0);
			
			
		} catch (InvalidFormatException e) {
			System.out.println("DRDS xlsx file is not in the correct format, exiting");
			e.printStackTrace();
			System.exit(1);
		} catch (FileNotFoundException e) {
			System.out.println("DRDS xlsx file not found, exiting: " + this.inputFile.toString());
			e.printStackTrace();
			System.exit(1);
		} catch (IOException e) {
			System.out.println("Error reading DRDS xlsx file, exiting");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	private void writeWorkBook() {
		//Write workbook
		FileOutputStream fileOut;
		try {
			fileOut = new FileOutputStream(this.outputFile);
			wb.write(fileOut);
			fileOut.close();
		} catch (FileNotFoundException e) {
			System.out.println("Could not find output file, exiting: " + this.outputFile.getAbsolutePath());
			e.printStackTrace();
			System.exit(1);
		} catch (IOException e) {
			System.out.println("Error writing to file, exiting: " + this.outputFile.getAbsolutePath());
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	private void parseBiomart() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.annotationFile));
			
			String[] header = br.readLine().split("\t");
			if (header.length <= 2) {
				System.out.println(String.format("The annotation file only has %s columns, which isn't enough for annotation.  "
						+ "Please check to make sure your file is in tab-delimited format. Exiting.", header.length));
				System.exit(1);
			}
			if (!header[0].equals("Ensembl Gene ID")) {
				System.out.println("WARNING! The first column of the annotation file is not 'Ensembl Gene ID', if Ensembl Gene ID "
						+ "information isn't in the first column, annotation might fail.\n");
			} 
			if (!header[1].equals("Ensembl Transcript ID")) {
				System.out.println("WARNING! This second column of the annotaiton file is not 'Ensembl Transcript ID', the second"
						+ " column of the annotation file is skipped, so make sure you don't have important information in this "
						+ " column");
			}
			
			//Set header
			this.header = Arrays.copyOfRange(header, 2, header.length);
			
			//Set entries
			this.entries = header.length - 2;
			
			//Initialize annotation container
			this.annotations = new HashMap<String,String[]>();
			
			//Read through file and add annotiation information
			String temp = null;
			while ((temp = br.readLine()) != null)  {
				String[] items = temp.split("\t",-1);
				if (!this.annotations.containsKey(items[0])) {
					this.annotations.put(items[0], Arrays.copyOfRange(items, 2, items.length));
				}
				
			}
			
			
			br.close();
		} catch (FileNotFoundException fnfe) {
			System.out.println("Could not find annotation file, exiting: " + this.annotationFile.toString());
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Error reading annotation file, exiting: " + this.annotationFile.toString());
			ioex.printStackTrace();
			System.exit(1);
		}
	}
	
	private void parseArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'o': outputFile = new File(args[++i]); break;
					case 'a': annotationFile = new File(args[++i]); break;
					case 'i': inputFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		
		//Make sure all arguments were specified
		if (outputFile.getAbsolutePath().equals(inputFile.getAbsolutePath())) {
			exitMessage("Output file path and input file path are the same, exiting");
		}
		
		if (outputFile == null) {
			exitMessage("Annotated output file not specified (-o), exiting.");
		}
		
		if (inputFile == null) {
			exitMessage("DRDS xlsx file not specified (-i), exiting.");
		}
		
		if (annotationFile == null) {
			exitMessage("Biomart annotation file not specified (-a), exiting");
		}	
		
		//Make sure files exist
		if (!inputFile.exists()) {
			exitMessage("DRDS xlsx file specifed does not exist, exiting: " + this.inputFile.toString());
		}
		
		if (!annotationFile.exists()) {
			exitMessage("Biomart annotation file specified does not exist, exiting: " + this.annotationFile.toString());
		}
		
	    //Makes sure the suffixes are correct	
		if (!inputFile.toString().endsWith(".xlsx")) {
			exitMessage("The DRDS xlsx file specied does not end with xlsx, exiting: " + this.inputFile.toString());
		}
		
		if (!outputFile.toString().endsWith(".xlsx")) {
			exitMessage("The annotated output file specified does not end with xlsx, exiting" + this.outputFile.toString());
		}
	}
	
	private void exitMessage(String message) {
		printDocs();
		System.out.println(message);
		System.exit(0);
	}
	
	private static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                DRDS Annotator: January 2014                      **\n" +
				"**************************************************************************************\n" +
				"This application annotates DefinedRegionDifferentialSeq xlsx files using Ensembl \n" + 
				"biomart tab-delimited annotation files. By default, ensembl biomart output files will \n" +
				"list the Ensembl gene id in the first column and Ensembl transcript id in the second \n" +
				"column.  This application assumes these defaults.  It will match the gene id in the \n" +
			    "first column of the biomart file to the name listed in the 'IGB HyperLink' column \n" +
				"found in the 'Analyzed Genes' tab of the DRDS xlxs output. All biomart columns after \n" +
			    "the transcript id column are added to the output file.  The data is inserted between \n" +
				"the 'Alt Name' and locus columns in the 'Analyzed Genes' tab.\n\n" +
				"The biomart output files can have multiple annotation lines for each gene id.  \n" +
				"Currently, this app uses the first annotation line encountered.\n\n" +
				
				"\nRequired Arguments:\n\n"+
				"-i Input file. Path to DRDS xlsx output file you wish to annotate \n" +
				"-a Annotation file. Path to biomart annotation file. \n" +
				"-o Annotated output file. Path to the annotated output file\n" +
				
				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/DRDSAnnotator -i geneStats.xlsx \n" +
				"               -a mm10.biomart.txt -o geneStats.ann.xlsx\n\n" +

				"**************************************************************************************\n");

	}

}
