package util.apps;
import java.io.*;
import java.util.ArrayList;

import util.gen.*;

/**Utility app to pull classes into folders and extract menus to html.*/
public class MakeBioTools {
	//fields
	
	private File jarAppDir = null;
	private File srcDir = null;
	
	//For the T2 package
	/**Relative directories to add to a library code jar for T2*/
	private String[] t2DirForLibJar = {
			"trans",
			"util", 
			"igb",
			"expr",
			"edu",
			"gata",
			"meme",
			"selex",
			"com",
			"javax",
			"org"	
	};
	/**Classes to turn into jar files*/
	private String[] t2ClassesToJar = {
			"edu/utah/seq/analysis/AggregatePlotter",
			"trans/anno/AnnotateRegions",
			"trans/misc/Bar2Gr",
			"trans/cel/CelFileConverter",
			"trans/qc/CelFileQualityControl",
			"trans/graphics/CelMasker",
			"trans/cel/CelProcessor",
			"trans/qc/CoordinateExtractor1lq",
			"trans/misc/ConvertAgilentData",
			"trans/misc/ConvertGeoData",
			"trans/misc/ConvertNimblegenNDF2TPMap",
			"trans/misc/ConvertNimblegenPAIR2Cela",
			"expr/CorrelationMaps",
			"util/apps/Correlate",
			"trans/misc/BestWindowScoreExtractor",
			"util/bio/annotation/ExportIntergenicRegions",
			"util/bio/annotation/ExportIntronicRegions",
			"trans/main/FDRWindowConverter",
			"util/bio/parsers/FetchGenomicSequences",
			"trans/anno/FindNeighboringGenes",
			"expr/FileCrossFilter",
			"util/apps/FileJoiner",
			"util/apps/FileSplitter",
			"trans/tpmap/FilterTPMapByRegions",
			"trans/main/FindSubBindingRegions",
			"trans/misc/Gr2Bar",
			"util/bio/cluster/HierarchicalClustering",
			"util/bio/seq/IndexFastas",
			"trans/main/IntensityPrinter",
			"util/apps/IntersectLists",
			"trans/anno/IntersectRegions",
			"trans/main/IntervalFilter",
			"trans/main/IntervalGFFPrinter",
			"trans/main/IntervalGraphPrinter",
			"trans/main/IntervalMaker",
			"trans/graphics/IntervalPlotter",
			"trans/main/IntervalReportPrinter",
			"util/apps/JQSub",
			"trans/main/LoadChipSetIntervalOligoInfo",
			"trans/main/LoadIntervalOligoInfo",
			"trans/cel/MakeChromosomeSets",
			"trans/main/MergeWindowArrays",
			"util/apps/MergeRegions",
			"trans/main/MultiWindowIntervalMaker",
			"trans/tpmap/MummerMapper",
			"trans/main/OligoIntensityPrinter",
			"util/bio/seq/OligoTiler",
			"trans/main/OverlapCounter",
			"trans/roc/ScoreParsedBars",
			"util/bio/wrappers/Primer3Wrapper",
			"util/apps/PrintSelectColumns",
			"util/apps/RandomizeTextFile",
			"trans/graphics/RankedSetAnalysis",
			"trans/main/ScanChip",		
			"trans/main/ScanChromosomes",
			"trans/main/ScanGenes",
			"util/apps/ScatterPlot",
			"util/apps/ScoreChromosomes",
			"trans/main/ScoreIntervals",
			"util/apps/ScoreSequences",
			"trans/main/SetNumberIntervalMaker",
			"trans/misc/Sgr2Bar",
			"util/apps/SubtractRegions",
			"util/bio/parsers/SynonymMatching",
			"trans/main/T2",
			//"edu.utah.seq/Tag2Bar",
			"trans/tpmap/TPMapOligoBlastFilter",
			"trans/tpmap/TPMapProcessor",
			"trans/tpmap/TPMapSort",
			"trans/graphics/VirtualCel",
			"igb/util/Windows2HeatMapSgr",
	};
	



	//constructor
	public MakeBioTools (){
		//make bioToolsCodeLibrary.jar
System.out.println("MakingLibrary");		
		makeLibrary();
		//make source code
System.out.println("MakingSource");
		makeSource();
		//make command menu html doc
System.out.println("MakingMenuHTML");
		makeMenuHTMLDoc();
		//make jar files for each class
System.out.println("MakingMockJars");		
		makeMockJars();
		//clean up and pack
		/*
		File t2 = new File ("T2");
		if (t2.exists()) t2.delete();
		t2.mkdir();
		//move SourceCode Dir to T2
		File move = new File (t2,srcDir.getName());
		srcDir.renameTo(move);
		System.out.println("cp *html T2; cp *txt T2; cp -R doc/ T2/; cp -R OthersCode/RBSymPTest T2/\n");
		*/
		
		
	}
	
	//methods
	/**Makes mock Jar files.*/
	public void makeMockJars(){
		File manFile = null;
		//for each class
		for (int i=0; i< t2ClassesToJar.length; i++){
			//check if class file exists
			File classFile = new File(t2ClassesToJar[i]+".class");
			if (classFile.exists() == false) Misc.printExit("\nError: class file does not exist "+classFile);
			//make and write manifest
			String manifest = "Manifest-Version: 1.0\nMain-Class: " + t2ClassesToJar[i]+ "\nClass-Path: ../LibraryJars/bioToolsCodeLibrary.jar\n";
			manFile = new File ("MANIFEST.MF");
			IO.writeString(manifest, manFile);
			//make dummy jar
			String[] tokens = t2ClassesToJar[i].split("/");
			String[] jarCmd = {"jar", "cmf0", manFile.getName(), tokens[tokens.length-1]};
			String[] results = IO.executeCommandLineReturnAll(jarCmd);
			Misc.printArray(results);
			if (results.length != 0){
				Misc.printArray(results);
				Misc.printExit("\nError: problem with making dummy jar for "+t2ClassesToJar[i]);
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

		for (int i=0; i< t2ClassesToJar.length; i++){
			String[] cmd = {"java", t2ClassesToJar[i]};
			String[] menu = IO.executeCommandLineReturnAll(cmd);
			if (menu[0].startsWith("Exception")) {
				Misc.printArray(menu);
				Misc.printExit("\nError: problem lauching "+t2ClassesToJar[i]);
			}
			String[] tokens = t2ClassesToJar[i].split("/");
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
		File menuFile = new File ("T2Documentation/cmdLnMenus.html");
		IO.writeString(sb.toString(), menuFile); 
	}
	
	/**Makes the code library jar from the directoriesForLibJar array.
	 * Saves as bioToolsCodeLibrary.jar*/
	public boolean makeLibrary(){
		
		jarAppDir = new File ("Apps");
		if (jarAppDir.exists()) IO.deleteDirectory(jarAppDir);
		jarAppDir.mkdir();
		for (int i=0; i< t2DirForLibJar.length; i++){
			File source = new File (t2DirForLibJar[i]);
			if (source.exists()== false) Misc.printExit("\nError: cannot find source directory "+source);
			if (IO.copyDirectoryRecursive(source, new File (jarAppDir,source.getName()),"class$|properties$") == false) Misc.printExit("\nError: directory copy failed for "+source);
		}
		//jar it		
		StringBuffer sb = new StringBuffer();
		sb.append("\n**** Type on the command line****\n\ncd ");
		sb.append(IO.getFullPathName(jarAppDir));
		sb.append("\njar cf0 bioToolsCodeLibrary.jar ");
		sb.append(Misc.stringArrayToString(t2DirForLibJar," "));
		sb.append("\ncd ..\nrm -rf ");
		for (int i=0; i<t2DirForLibJar.length; i++){
			sb.append(" Apps/");
			sb.append(t2DirForLibJar[i]);
		}
		sb.append(
				"\n"+
				"mkdir T2\n"+
				"mkdir T2/Documentation T2/LibraryJars T2/OthersCode \n"+
				"zip -r -q SourceCode SourceCode\n"+
				"mv SourceCode.zip T2/\n"+
				"rm -rf SourceCode\n"+
				"cp -R T2Documentation/ T2/Documentation/\n"+
				"cp -R OthersCode/ T2/OthersCode/\n"+
				"mv Apps/bioToolsCodeLibrary.jar T2/LibraryJars/\n"+
				"mv Apps/ T2/\n");
		
		System.out.println(sb);
		  
		return true;
	}
	
	/**Makes the code library jar from the directoriesForLibJar array.
	 * Saves as bioToolsCodeLibrary.jar*/
	public void makeSource(){
		srcDir = new File ("SourceCode");
		if (srcDir.exists()) IO.deleteDirectory(srcDir);
		srcDir.mkdir();
		for (int i=0; i< t2DirForLibJar.length; i++){
			File source = new File (t2DirForLibJar[i]);
			if (source.exists()== false) Misc.printExit("\nError: cannot find source directory "+source);
			if (IO.copyDirectoryRecursive(source, new File (srcDir,source.getName()),".java") == false) Misc.printExit("\nError: directory copy failed for "+source);
		}
	}
	
	public static void main (String[] args){
		new MakeBioTools();
	}
	
	public static final String menuHeader = "<html><head><title>TiMAT2 Command Line Menus</title>" +
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
