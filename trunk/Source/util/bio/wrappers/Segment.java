package util.bio.wrappers;
import java.io.*;
import util.gen.*;

/**This wraps Huber's segment() function in his 'tilingArray' R package.
 * Use it to divide each column of float[]s, in a tab delimited file, into 3 segments, 2 indexes.
 * The first row in the file should contain the number of subsequent floats in each column to use during segmentation.
 * Must have installed the BioConductor tilingArray package and it's 1/2 dozen or so dependancies.*/
public class Segment {
	//fields
	private File rApp;
	private File dataFile;
	private File tempDirectory;
	private File tempRScriptFile;
	private File tempResultsFile;
	private File tempROutputFile;
	/**The segment slice occurs immediately before a given index, eg 5 means a cut between indexes 4 and 5.
	 * int[columns][2 indexes]*/
	private int [][] segmentIndexes;
		
	//constructor
	public Segment (File rApp, File dataFile, File tempDirectory){
		this.rApp = rApp;
		this.dataFile = dataFile;
		this.tempDirectory = tempDirectory;
	}
	
	//methods
	public void segment() throws Exception{
		//check files
		checkFiles();
		
		//make temp files
		String randomWord = Passwords.createRandowWord(8);
		tempRScriptFile = new File (tempDirectory, "tempFile_RScript_"+randomWord+".txt");
		tempROutputFile = new File (tempDirectory, "tempFile_RScript_"+randomWord+".txt.Rout"); 
		tempResultsFile = new File (tempDirectory, "tempFile_Results_"+randomWord+".txt");
		
		//write R script file
		writeRScriptFile();
		
		//execute script
		executeRScript();
		
		//parse results
		parseResults();
	}
	
	public void deleteTempFiles(){
		tempRScriptFile.delete();
		tempResultsFile.delete();
		tempROutputFile.delete();
	}
	
	public void parseResults() throws Exception{
		segmentIndexes = IO.loadTableOfInts(tempResultsFile);
		if (segmentIndexes == null || segmentIndexes[0].length !=2) throw new Exception ("Problem with R results");
		deleteTempFiles();
		//subtract one to get into java zero based array coordinates, R is 1 based
		for (int i=0; i< segmentIndexes.length; i++){
			for (int j=0; j< segmentIndexes[i].length; j++){
				segmentIndexes[i][j]--;
			}
		}
	}
	
	public void executeRScript() throws Exception{
		//make command
		String[] command = new String[] {
				rApp.getCanonicalPath(),
				"CMD",
				"BATCH",
				tempRScriptFile.getCanonicalPath(),
				tempROutputFile.getCanonicalPath()};			
		//execute
		String[] results = IO.executeCommandLine(command);
		if (results == null || results.length > 1 && results[0].trim().length()!=0) throw new Exception ("Problem executing batch R command");
	}
	
	public void checkFiles() throws Exception{
		if (rApp == null || rApp.canRead() == false) throw new Exception ("Cannot read R app "+rApp);
		if (dataFile == null || dataFile.canRead() == false) throw new Exception ("Cannot read data file "+dataFile);
		if (tempDirectory == null || tempDirectory.isDirectory() == false || tempDirectory.canWrite() == false) throw new Exception ("Cannot find or write to temp results directory "+tempDirectory);
	}
	
	public void writeRScriptFile() throws Exception{
		String script = 
			"#to launch on the command line '"+rApp+" CMD BATCH "+dataFile+"'\n"+
			"\n"+
			"#set data and results files\n"+
			"dataFile = \""+dataFile+"\"\n"+
			"resultsFile = \""+tempResultsFile+"\"\n"+
			"\n"+
			"#load library\n"+
			"library(\"tilingArray\")\n"+
			"\n"+
			"#read in data\n"+
			"data <- read.table(dataFile, sep=\"\\t\")\n"+
			"\n"+
			"#how many columns\n"+
			"numColumns = ncol(data)\n"+
			"\n"+
			"#make matrix to hold results\n"+
			"segMatrix = matrix(nrow = 2, ncol = numColumns)\n"+
			"\n"+
			"#set max segments\n"+
			"maxseg = 3\n"+
			"\n"+
			"#for each column\n"+
			"for (x in 1:numColumns){\n"+
			"	#read number of values and assign to maxk, data[rows,columns]\n"+
			"	maxk = data[1,x]\n"+
			"	finalRow = maxk + 1\n"+
			"	#launch app\n"+
			"	seg = segment (data[2:finalRow,x], maxseg, maxk)\n"+
			"	#assign the results, they are the index immediately before the breakpoint\n"+
			"	#thus a breakpoint of 9 would occur between 8 and 9\n"+
			"	segMatrix[1,x] = seg@breakpoints[3][[1]][1]\n"+
			"	segMatrix[2,x] = seg@breakpoints[3][[1]][2]\n"+
			"}\n"+
			"\n"+
			"#save results to file\n"+
			"write.table(segMatrix, row.names=FALSE, col.names=FALSE, sep=\"\\t\", file = resultsFile)\n";
		//System.out.println(script);

		if (IO.writeString(script, tempRScriptFile) == false) {
			tempRScriptFile.delete();
			throw new Exception ("Cannot write temp R script file "+tempRScriptFile);
		}
	}
	
	/*
	public static void main (String[] args){
		File rApp = new File ("/usr/bin/R");
		File dataFile = new File ("/Users/nix/Desktop/delme.txt");
		File tempDir = new File ("/Users/nix/Desktop/");
		Segment seg = new Segment (rApp, dataFile, tempDir);
		try{
			seg.segment();
			Misc.printArray(seg.segmentIndexes);
		} catch (Exception e){
			e.printStackTrace();
		}
	} */


	public int[][] getSegmentIndexes() {
		return segmentIndexes;
	}
}
