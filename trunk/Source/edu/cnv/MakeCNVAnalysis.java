package edu.cnv;

import java.io.*;
import java.util.*;
import util.gen.*;

/**Utility app to pull classes into folders and extract menus to html.*/
public class MakeCNVAnalysis {
	//fields
	
	private File jarAppDir = null;
	private File srcDir = null;
	
	//For the USeq package
	/**Relative directories to add to a library code jar*/
	private String[] dirForLibJar = {
			"edu",
			"trans",
			"expr",
			"util", 
			"cnv"
	};
	
	/**Classes to turn into jar files*/
	private String[] classesToJar = {
			"trans/misc/Bar2Gr",
			"cnv/CNVScanner",
			"util/bio/parsers/FetchGenomicSequences",
			"trans/anno/FindNeighboringGenes",
			"expr/FileCrossFilter",
			"expr/FileMatchJoiner",
			"util/apps/FileJoiner",
			"util/apps/FileSplitter",
			"trans/misc/Gr2Bar",
			"cnv/IntersectCNVs",
			"util/apps/IntersectLists",
			"trans/anno/IntersectRegions",
			"util/bio/wrappers/Primer3Wrapper",
			"util/apps/PrintSelectColumns",
			"trans/roc/ScoreParsedBars",
			"trans/misc/Sgr2Bar",
			"trans/misc/Wig2Bar",
	};

	//constructor
	public MakeCNVAnalysis (){
		//make bioToolsCodeLibrary.jar		
		makeLibrary();
		//make source code
		makeSource();
		//make command menu html doc
		makeMenuHTMLDoc();
		//make jar files for each class		
		makeMockJars();
	}
	
	//methods
	/**Makes mock Jar files.*/
	public void makeMockJars(){
		File manFile = null;
		//for each class
		for (int i=0; i< classesToJar.length; i++){
			//check if class file exists
			File classFile = new File(classesToJar[i]+".class");
			if (classFile.exists() == false) Misc.printExit("\nError: class file does not exist "+classFile);
			//make and write manifest
			String manifest = "Manifest-Version: 1.0\nMain-Class: " + classesToJar[i]+ "\nClass-Path: ../LibraryJars/bioToolsCodeLibrary.jar\n";
			manFile = new File ("MANIFEST.MF");
			IO.writeString(manifest, manFile);
			//make dummy jar
			String[] tokens = classesToJar[i].split("/");
			String[] jarCmd = {"jar", "cmf0", manFile.getName(), tokens[tokens.length-1]};
			String[] results = IO.executeCommandLineReturnAll(jarCmd);
			Misc.printArray(results);
			if (results.length != 0){
				Misc.printArray(results);
				Misc.printExit("\nError: problem with making dummy jar for "+classesToJar[i]);
			}
			//move into Apps
			File jarApp = new File (tokens[tokens.length-1]);
			File jarAppInDir = new File (jarAppDir, jarApp.getName());
			jarApp.renameTo(jarAppInDir); 
		}
		manFile.delete();
	}
	 
	/**Generates the cmdLnMenus.html doc.*/
	public void makeMenuHTMLDoc(){
		StringBuffer names = new StringBuffer();
		StringBuffer menus = new StringBuffer();

		for (int i=0; i< classesToJar.length; i++){
			String[] cmd = {"java", classesToJar[i]};
			String[] menu = IO.executeCommandLineReturnAll(cmd);
			if (menu[0].startsWith("Exception")) {
				Misc.printArray(menu);
				Misc.printExit("\nError: problem lauching "+classesToJar[i]);
			}
			String[] tokens = classesToJar[i].split("/");
			String jarName = tokens[tokens.length-1];
			names.append("<a href=\"#"+  jarName +"\">"+ jarName +"</a><br>\n");
			menus.append("<a text=\""+ jarName +"\"><pre>\n");
			for (int j=0; j<menu.length; j++){
				menus.append(menu[j]);
				menus.append("\n");
			}
			menus.append("</pre><p>\n");
		}
		
		StringBuffer sb = new StringBuffer(menuHeader);
		sb.append(names);
		sb.append("<p>");
		sb.append(menus);
		sb.append("</body></html>");
		File menuFile = new File ("CNVAnDocumentation/cmdLnMenus.html");
		IO.writeString(sb.toString(), menuFile); 
	}
	
	/**Makes the code library jar from the directoriesForLibJar array.
	 * Saves as bioToolsCodeLibrary.jar*/
	public boolean makeLibrary(){
		
		jarAppDir = new File ("Apps");
		if (jarAppDir.exists()) IO.deleteDirectory(jarAppDir);
		jarAppDir.mkdir();
		for (int i=0; i< dirForLibJar.length; i++){
			File source = new File (dirForLibJar[i]);
			if (source.exists()== false) Misc.printExit("\nError: cannot find source directory "+source);
			if (IO.copyDirectoryRecursive(source, new File (jarAppDir,source.getName()),"class$|properties$") == false) Misc.printExit("\nError: directory copy failed for "+source);
		}
		//jar it		
		StringBuffer sb = new StringBuffer();
		sb.append("\n**** Type on the command line****\n\ncd ");
		sb.append(IO.getFullPathName(jarAppDir));
		sb.append("\njar cf0 bioToolsCodeLibrary.jar ");
		sb.append(Misc.stringArrayToString(dirForLibJar," "));
		sb.append("\ncd ..\nrm -rf ");
		for (int i=0; i<dirForLibJar.length; i++){
			sb.append(" Apps/");
			sb.append(dirForLibJar[i]);
		}
		sb.append(
				"\n"+
				"mkdir CNVAn\n"+
				"mkdir CNVAn/Documentation CNVAn/LibraryJars \n"+ // USeq/OthersCode \n"+
				"zip -r -q SourceCode SourceCode\n"+
				"mv SourceCode.zip CNVAn/\n"+
				"rm -rf SourceCode\n"+
				"cp -R CNVAnDocumentation/ CNVAn/Documentation/\n"+
				"mv Apps/bioToolsCodeLibrary.jar CNVAn/LibraryJars/\n"+
				"mv Apps/ CNVAn/\n" +
				"zip -r -q CNVAn.zip CNVAn/\n\n");
		
		System.out.println(sb);
		  
		return true;
	}
	
	/**Makes the code library jar from the directoriesForLibJar array.
	 * Saves as bioToolsCodeLibrary.jar*/
	public void makeSource(){
		srcDir = new File ("SourceCode");
		if (srcDir.exists()) IO.deleteDirectory(srcDir);
		srcDir.mkdir();
		for (int i=0; i< dirForLibJar.length; i++){
			File source = new File (dirForLibJar[i]);
			if (source.exists()== false) Misc.printExit("\nError: cannot find source directory "+source);
			if (IO.copyDirectoryRecursive(source, new File (srcDir,source.getName()),".java") == false) Misc.printExit("\nError: directory copy failed for "+source);
		}
	}
	
	public static void main (String[] args){
		new MakeCNVAnalysis();
	}
	
	public static final String menuHeader = "<html><head><title>CNVAn Command Line Menus</title>" +
			"<style type=\"text/css\">#rt{text-align:right; color: #000000; font-weight: bold}" +
			"#grBk {background-color: #CC9966;}TD {font-family: Verdana, Arial, Helvetica, sans-serif; " +
			"font-size:12;}H1 {color: #996633; font:arial; font-size:16;}H2 {color: #996633; " +
			"font:arial; font-size:12;}BODY {color:black; background-color:white; font-family: " +
			"Verdana, Arial, Helvetica, sans-serif; font-size:12;}A:link    {text-decoration: none; " +
			"color: #000000; font-weight: bold}  A:visited {text-decoration: none; color: #000000; " +
			"font-weight: bold}   A:hover   {text-decoration: none; color: #FFCC66; font-weight: bold} " +
			"A:active  {text-decoration: none; color: #000000; font-weight: bold}   </style></head><body>" +
			"<H1>Command Line Menus</H1> <p>";
	
}
