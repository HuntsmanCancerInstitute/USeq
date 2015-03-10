package edu.utah.seq.cnv;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class PoReDataChunk extends Thread{

	//fields
	private GeneExonSample[] exonSample;
	private PoReCNV poisRegCNV;
	private float significanceThreshold;
	private String name;
	
	//0 waiting to start, 1 running, 2 complete w/o errors, 3 errors
	private int status = 0;
		
	//temp files
	private File countTableFile;
	private File scriptFile;
	private File rOut;
	private File residualsFile;
	private File obsExpLg2RtosFile;
	private File sigLevelFile;
	
	//constructor
	public PoReDataChunk (String name, ArrayList<GeneExonSample>[] s, PoReCNV poisRegCNV){
		//get total counts
		int num = 0;
		for (int i=0; i< s.length; i++) num+= s[i].size();
		//convert to standard array
		exonSample = new GeneExonSample[num];
		int index = 0;
		for (int i=0; i< s.length; i++){
			for (GeneExonSample e: s[i]) exonSample[index++] = e;
		}
		this.poisRegCNV = poisRegCNV;
		this.name = name;
	}

	//methods
	public void run(){
		status = 1;
		writeDataTable();
		if (status !=3) executeCaseAnalysisScript();
		if (status !=3) loadCaseAnalysisResults();
		if (status !=3) {
			if (poisRegCNV.isDeleteTempFiles()){
				countTableFile.delete();
				countTableFile.delete();
				rOut.delete();
				scriptFile.delete();
				residualsFile.delete();
				obsExpLg2RtosFile.delete();
				sigLevelFile.delete();
			}
			status = 2;
		}
	}
	
	private void writeDataTable() {
		countTableFile = new File(poisRegCNV.getSaveDirectory(), name+"CountTable.txt");
		try {
			PrintWriter out = new PrintWriter( new FileWriter(countTableFile));
			for (int i=0; i< exonSample.length; i++) out.println(exonSample[i].fetchRDataLine());
			out.close();
		}
		catch (Exception e){
			System.err.println("\nProblem writing out count table for "+name);
			e.printStackTrace();
			status = 3;
		}
	}

	private void executeCaseAnalysisScript() {
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("source('"+ poisRegCNV.getAlunRScript()+ "')\n");
			sb.append("x = readdata('"+countTableFile+"',3,3)\n");
			sb.append("x = fitmodel(x)\n");
			sb.append("exportCNVData(x, '"+poisRegCNV.getSaveDirectory() 
					+"', '"+name+"', testsig="+poisRegCNV.getMinimumAdjPVal()
					+", ntests="+poisRegCNV.getNumExonsProcessed()+" )\n");

			//write script to file
			scriptFile = new File (poisRegCNV.getSaveDirectory(), name+"Cnv_RScript.R");
			rOut = new File(poisRegCNV.getSaveDirectory(), name+"Cnv_RScript.Rout");
			IO.writeString(sb.toString(), scriptFile);

			//make command
			String[] command = new String[] {
					poisRegCNV.getFullPathToR().getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};

			//execute command
			IO.executeCommandLine(command);

			//look for results files
			residualsFile = new File (poisRegCNV.getSaveDirectory(), name+"Residuals.txt");
			obsExpLg2RtosFile = new File (poisRegCNV.getSaveDirectory(), name+"ObsExpLg2Rtos.txt");
			sigLevelFile = new File (poisRegCNV.getSaveDirectory(), name+"SigLevel.txt");
			if (residualsFile.exists() == false || obsExpLg2RtosFile.exists() == false || sigLevelFile.exists() == false){
				throw new IOException("\nError: cannot find the R results files? Check the "+rOut.getName()+" log file for errors.\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Error: failed to execute case only cnv analysis in R for "+name+"\n");
			status = 3;
			System.exit(1);
		}
	}

	private void loadCaseAnalysisResults() {
		try {
			//load data
			float[] residuals = Num.loadFloats(residualsFile);
			float[] obsExpRtos = Num.loadFloats(obsExpLg2RtosFile);
			float[] sigLevelArray = Num.loadFloats(sigLevelFile);
			//check
			if (residuals == null || residuals.length != exonSample.length || obsExpRtos == null || obsExpRtos.length != exonSample.length || sigLevelArray == null || sigLevelArray.length !=1){
				throw new Exception("\nError: cannot load appropriate data from R results, check logs.\n");
			}
			//load em
			for (int i=0; i< exonSample.length; i++){
				exonSample[i].setResidual(residuals[i]);
				exonSample[i].setObsExpLgRto(obsExpRtos[i]);
			}
			significanceThreshold = sigLevelArray[0];

		} catch (Exception e){
			e.printStackTrace();
			System.err.println("\nError: problem parsing results from "+name);
			status = 3;
		}
	}

	//getters and setters
	public float getSignificanceThreshold() {
		return significanceThreshold;
	}
	public int getStatus() {
		return status;
	}

	public String getRealName() {
		return name;
	}

}
