package trans.main;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import trans.tpmap.MapSplitter;
import util.gen.*;

/**
 * Prints out the oligo intensity values, split by chromosome, as '.sgr' files for import into Affy's IGB.
 */
public class OligoIntensityPrinter {
	
	//fields
	private File tpmapFile; 
	private File infoFile;
	private File treatmentDirFile; 
	private File controlDirFile;
	private File tempFile;
	
	public OligoIntensityPrinter (String[] args){
		processArgs(args); 
		if (treatmentDirFile.isDirectory()) tempFile = new File(treatmentDirFile,"OIPDelMe");
		else tempFile = new File (treatmentDirFile.getParent(), "OIPDelMe");
		
		
		//fetch float[][] of processed intensity values for treatment and controls
		float[][] treatments = ScanChip.fetchFloatArrays(IO.getFullPathName(treatmentDirFile),false);
		float[] treatment = Num.averageFloatArraysFlippedToFloats(treatments);
		treatments = null;
		
		float[][] controls = ScanChip.fetchFloatArrays(IO.getFullPathName(controlDirFile), false);
		float[] control = Num.averageFloatArraysFlippedToFloats(controls);
		controls = null;
		
		//calculate stats
		System.out.println("\nCalculating statistics...");
		double meanTreatment = Num.mean(treatment);
		double meanControl = Num.mean(control);
		double stndTreatment = Num.standardDeviation(treatment,meanTreatment);
		double stndControl = Num.standardDeviation(control, meanControl);
		System.out.println("\t"+(float)meanTreatment+"\tMean Treatment");
		System.out.println("\t"+(float)stndTreatment+"\tStandard Deviation Treatment");
		System.out.println("\t"+(float)meanControl+"\tMean Control");
		System.out.println("\t"+(float)stndControl+"\tStandard Deviation Control\n");
		
		//get .tpmap file info
		ArrayList info = (ArrayList)IO.fetchObject(infoFile);
		
		//break both apart by chromosome using the tpmapInfo ArrayList, save to disk, 
		//	might be able to get away with keeping em in memory?  Only need to run once.

		String uniqueIdT = MapSplitter.breakSaveIntensityValues(info, treatment, IO.getFullPathName(tempFile));
		String uniqueIdC = MapSplitter.breakSaveIntensityValues(info, control, IO.getFullPathName(tempFile));
		
		//for each chromosome, safe to parallel process on cluster
		int num = info.size();
		String chromosome;
		int[] positions;
		
		//make print writers
		PrintWriter outT = null;
		PrintWriter outC = null; 
		PrintWriter outR = null; 
		try {
			outT = new PrintWriter(new FileWriter( new File (tempFile.getParentFile(), "aveTreatment.sgr")));
			outC = new PrintWriter(new FileWriter( new File (tempFile.getParentFile(), "aveControl.sgr")));
			outR = new PrintWriter(new FileWriter( new File (tempFile.getParentFile(), "aveRatio.sgr")));
		}catch(IOException e){
			e.printStackTrace();
			System.exit(1);
		}
		
		//for each chromosome
		for (int i=1; i<num; i+=4){
			//fetch windows to scan, the int[][] made in the WindowMaker
			chromosome = (String)info.get(i);
			System.out.println("\tTesting chromosome: "+chromosome);
			
			//fetch treatment and control normalized, transformed, and chromosome divided int[]'s
			treatment = (float[])IO.fetchObject(new File(tempFile+chromosome+uniqueIdT));
			control = (float[])IO.fetchObject(new File(tempFile+chromosome+uniqueIdC));

			//fetch bp positions to use in converting window indexes to real base pairs
			positions = (int[])IO.fetchObject(new File(tpmapFile+chromosome));

			//run oligos
			int numberOligos = treatment.length;
			float ratio;
			String pos;
			for (int j=0; j<numberOligos; j++){
				//calculate intensity difference and ratio 
				if (control[j] !=0) ratio = treatment[j]/control[j];
				else ratio = 0;
				pos = chromosome+"\t"+ positions[j]+"\t";
				//print treatment
				outT.println(pos+treatment[j]);
				//print control
				outC.println(pos+control[j]);
				//print ratio
				outR.println(pos+ratio);
			}
		}
		//close PrintWriters
		outT.close();
		outC.close();
		outR.close();
		
		//delete .delMeTmp files
		IO.deleteFiles(tempFile.getParent(), uniqueIdT);
		IO.deleteFiles(tempFile.getParent(), uniqueIdC);
	}
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File directory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': directory = new File(args[i+1]); i++; break;
					case 't': treatmentDirFile = new File(args[i+1]); i++; break;
					case 'c': controlDirFile = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check to see if they entered required params
		if (directory==null || directory.isDirectory() == false){
			System.out.println("\nCannot find your tpmap directory!\n");
			System.exit(0);
		}
		else {
			tpmapFile = new File(directory.toString(),"tpmap.fa"); 
			infoFile = new File(directory.toString(),"tpmap.faInfo");
		}
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Oligo Intensity Printer: Dec 2005                         **\n" +
				"**************************************************************************************\n" +
				"OIP prints to file the average treatment, control, and ratio oligo intensity scores\n" +
				"      (no windowing) in a .sgr text format for direct import into Affy's IGB.\n" +
				"\n" +
				"Use the following options when running OIP:\n\n" +
				"-b Full path directory text for the 'xxxTpmapFiles' generated by the TPMapProcessor\n" +
				"-t The full path file text for a directory or file containing serialized float[] \n" +
				"      'xxx.celp' treatment file(s).\n" +
				"-c The full path file text for a directory or file containing serialized float[]\n" +
				"      'xxx.celp' control file(s).\n" +
				"\n" +
				"Example: java -Xmx256M -jar pathTo/T2/Apps/OligoIntensityPrinter -b\n" +
				"      /affy/675bp3MinFsTPMapFiles/ -t /affy/t.celp -c /affy/c.celp \n" +
				"\n" +
				
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length <3){
			printDocs();
			System.exit(0);
		}
		new OligoIntensityPrinter(args);
	}
	
}
