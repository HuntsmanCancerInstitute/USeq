package edu.utah.ames.bioinfo;

import java.io.*;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class test {
	
	private static String x = "/Users/darren/Desktop/novoalignerTestDir/reports/autoalign_2013-01-02.txt";
	private static String freshDataReport = "/tomato/job/autoaligner/reports/autoalign_*.txt";
	private static String dest = "/Users/darren/Desktop/novoalignerTestDir/processedReports/";
	
	public static void parseSomething() throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(x));
		br.readLine();
		String line;
		
		while ((line = br.readLine()) != null) {

			String dataValue[] = line.split("\t");
			if (dataValue.length < 15) {
			//	System.out.println("you have issues.");
				continue;
			}
			else {
				Sample s = new Sample(dataValue);
			//	System.out.println("you are awesome.");
			}	
		}
		br.close();
	//	System.out.println("you are horribly lost.");
		
	//	System.out.println(freshDataReport.substring(32));
	}
	
	public static String[] fileNames(String directoryPath) throws IOException {
		
		String sampleID = "9786X2";
		String destin = "/Users/darren/Desktop/novoalignerTestDir/linkedFiles/";
		//TODO build path from sample data
		File dir = new File("/Users/darren/Desktop/novoalignerTestDir/testStuff/");
		//TODO build regex from sample data
		String fileMatch = sampleID + "[_0-9A-Z]+{6,}.txt.gz";
		Pattern p = Pattern.compile(fileMatch);
		Collection<String> files = new ArrayList<String>();
		
		if (dir.isDirectory()) {
			File[] listFiles = dir.listFiles();
			
			for (File file : listFiles) {
				if (file.isFile()) {
					files.add(file.getName());
				}
			}
			
			for (int i = 0; i < listFiles.length; i++) {
				//System.out.println(listFiles[i]);
				
				Matcher m = p.matcher(listFiles[i].toString());
				if (m.find()) {
					//System.out.println(m.group());
					String dirs = dir.toString() + "/";
					System.out.println(dirs);
					String targetFiles = dirs + m.group();
					//TODO use found filename as input name for soft linking
					Process process = Runtime.getRuntime().exec(new String[] {"ln", "-s", targetFiles, destin});
				}
			else {
					//System.out.println("it didn't work.");
				}
			}
		}
		return files.toArray(new String[]{});
	}
	
	public static void main(String[] args) throws Exception {
		fileNames("/Users/darren/Desktop/novoalignerTestDir/testStuff/");
		
		//parseSomething();
		GregorianCalendar gc = new GregorianCalendar();
		int year = gc.get(Calendar.YEAR);
	//	System.out.println(year);
		
		File file = new File(x);
	//	System.out.println(file.getName());
		
		boolean success = file.renameTo(new File(dest + file.getName()));
		if (!success) {
	//		System.out.println("\nit's not working\n");
		}
		
		String targetFiles = "/Users/darren/Desktop/novoalignerTestDir/testStuff/";
		String destin = "/Users/darren/Desktop/novoalignerTestDir/linkedFiles/";
		
		
		//soft link files
		//Process process = Runtime.getRuntime().exec(new String[] {"ln", "-s", targetFiles, destin});
		
		
		//get new analysis report number and create dir
		String line;
		InputStream stderr = null;
		InputStream stdout = null;
		String analysisNum = null;
		String analysisPath = null;
		
		//String cmd[] = {"./httpclient_create_analysis.sh", "-properties", "gnomex_httpclient.properties",
		//		"-serverURL", "https://bioserver.hci.utah.edu", "-name", "foo", "-folderName",
		//		"test", "-lab", "Brad Cairns", "-organism", "Human", "-genomeBuild", "hg19",
		//		"-analysisType", "Alignment", "-seqLane", "9769F1_1"};
		
		//Process p = Runtime.getRuntime().exec(cmd, null, new File("/Users/darren/Desktop/novoalignerTestDir/create_analysis/"));
		//stderr = p.getErrorStream();
		//stdout = p.getInputStream();
		
		//BufferedReader brCleanup = new BufferedReader(new InputStreamReader(stdout));
		
		//while (((line = brCleanup.readLine()) != null)) {
		//	System.out.println("[stdout] " + line);
		//	if (line.indexOf("idAnalysis=") > 0) {
				//grab the new analysis number in output
		//		String matchANum = "([A][0-9]+)";
		//		Pattern p1 = Pattern.compile(matchANum);
		//		Matcher m = p1.matcher(line);
		//		if (m.find()) {
					//assign analysis number to var
		//			analysisNum = m.group();
		//			System.out.println(analysisNum);
		//		}
		///	}
		//	String matchPath = "[/\\w]+[A][0-9]+";
		//	Pattern p2 = Pattern.compile(matchPath);
		//	Matcher m = p2.matcher(line);
		//	if (m.find()) {
				//set the path in sample object
		//		analysisPath = m.group() + "/";
		//		System.out.println(analysisPath);
		//	}
	//	}
	}
}
