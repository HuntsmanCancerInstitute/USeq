package edu.expr;
import java.io.*;
import util.gen.*;
import java.util.*;

public class ParseManipulateColumnDataForGeneSifter {

	public static void main (String[] args){
		File dataFile1 = new File (args[0]);
		File dataFile2 = new File (args[1]);
		File posIdsFile = new File (args[2]);
		File geneNamesFile = new File (args[3]);
		boolean convertLog10ToLog2 = true;
		boolean extraColumns = true;

		//load patient ids that are six1 positive others are assume to be six1 negative
		String[] patientIDs = IO.loadFileIntoStringArray(posIdsFile);

		//load gene names
		String[] geneNames = IO.loadFileIntoStringArray(geneNamesFile);

		//load each data file, should have all doubles except for top row which contains patientIDs
		HashMap idValues = Num.loadDoubles(dataFile1, convertLog10ToLog2);
		idValues.putAll(Num.loadDoubles(dataFile2, convertLog10ToLog2));
		System.out.println("\nNumber of data columns " +idValues.size());
		System.out.println("Number of positives "+ patientIDs.length);
		System.out.println("Number of negatives "+ (idValues.size() - patientIDs.length));

		//Sort into different classes by patient type
		double[][] patientData = new double[idValues.size()][];
		String[] names = new String[idValues.size()];
		//load patients
		int index = 0;
		for (; index < patientIDs.length; index++){
			patientData[index] = (double[]) idValues.remove(patientIDs[index]);
			if (extraColumns) names[index] = patientIDs[index] +"\tMSP+";
			else names[index] = patientIDs[index] +"_MSP+";
		}

		//System.out.println("Starting " +index);
		//load remainder
		Iterator it = idValues.keySet().iterator();
		while (it.hasNext()){
			//System.out.println("\t"+ index);
			String name = ((String)it.next());
			if (extraColumns) names[index] = name  +"\tMSP-";
			else names[index] = name  +"_MSP-";
			patientData[index] = (double[]) idValues.get(name);
			index++;
		}

		//print to file
		File sortedData = new File (posIdsFile.getParentFile(), "merged.txt");
		try {
			PrintWriter out = new PrintWriter( new FileWriter(sortedData));
			//print header
			StringBuffer sb = new StringBuffer();
			//leave two blank columns for gene names
			if (extraColumns) sb.append("\t\t");
			else sb.append("\t");
			sb.append(names[0]);
			for (int i=1; i< names.length; i++){
				sb.append("\t");
				sb.append(names[i]);
			}
			out.println(sb);
			
			//print data with extra blank columns
			//for each row
			for (int j=0; j< patientData[0].length; j++){
				sb = new StringBuffer();
				//append text line
				sb.append (geneNames[j]);
				sb.append ("\t");
				//append first
				sb.append(patientData[0][j]);
				//for each column
				for (int i=1; i< patientData.length; i++){
					if (extraColumns) sb.append ("\t\t");
					else sb.append ("\t");
					sb.append(patientData[i][j]);
				}
				out.println(sb);
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}


}
