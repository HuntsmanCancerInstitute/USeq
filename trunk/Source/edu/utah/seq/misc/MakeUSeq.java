package edu.utah.seq.misc;

import java.io.*;
import java.util.*;
import util.gen.*;

/**Utility app to build the USeq package. This has been hacked to death.  Anyone want to make an ant build script?*/
public class MakeUSeq {
	//fields

	private File jarAppDir = null;
	private String version;
	private boolean makeWebStart = false;

	//For the USeq package
	/**Relative directories to add to a library code jar for USeq*/
	private String[] uSeqDirForLibJar = {
			"trans",
			"util", 
			"edu",
			"meme",
			"com",
			"javax",
			"org",
			"info", 
			"net" 
	};

	/**Apps to turn into jar files*/
	private String[] uSeqClassesToJar = {
			"util/apps/ABITraceTCPeakCalculator",
			"edu/utah/seq/analysis/AggregatePlotter",
			"edu/expr/Alleler",
			"edu/utah/seq/analysis/AllelicMethylationDetector",
			"edu/utah/seq/base/BamIntensityJoiner",
			"edu/utah/seq/base/BamNMerIntensityParser",
			"trans/misc/Bar2Gr",
			"edu/utah/seq/useq/apps/Bar2USeq",
			"edu/utah/seq/base/BaseClassifier",
			"trans/misc/Bed2Bar",
			"edu/utah/seq/parsers/BedStats",
			"edu/utah/seq/analysis/BisSeq",
			"edu/utah/seq/analysis/BisSeqAggregatePlotter",
			"edu/utah/seq/data/BisSeqErrorAdder",
			"edu/utah/seq/analysis/BisStat",
			"edu/utah/seq/analysis/BisStatRegionMaker",
			"edu/utah/seq/analysis/ChIPSeq",
			"util/bio/wrappers/CHPCAligner",
			"util/apps/CompareIntersectingRegions",
			"edu/utah/seq/parsers/CompareParsedAlignments",
			"util/bio/seq/ConcatinateFastas",
			"edu/utah/seq/data/CorrelatePointData",
			"util/bio/seq/BisulfiteConvertFastas",
			"edu/expr/CorrelationMaps",
			"util/bio/seq/ConvertFastaA2G",
			"util/bio/seq/ConvertFastqA2G",
			"util/bio/seq/ConvertFasta2GCBoolean",
			"util/bio/seq/ConvertFasta2GCBarGraph",
			"edu/utah/seq/analysis/DefinedRegionBisSeq",
			"edu/utah/seq/analysis/DefinedRegionDifferentialSeq",
			"edu/utah/seq/analysis/DefinedRegionScanSeqs",
			"edu/utah/seq/analysis/EnrichedRegionMaker",
			"edu/utah/seq/parsers/ElandMultiParser",
			"edu/utah/seq/parsers/ElandParser",
			"edu/utah/seq/parsers/ElandSequenceParser",
			"util/bio/annotation/ExportExons",
			"util/bio/annotation/ExportIntergenicRegions",
			"util/bio/annotation/ExportIntronicRegions",
			"util/bio/annotation/ExportTrimmedGenes",
			"util/bio/parsers/FetchGenomicSequences",
			"trans/anno/FindNeighboringGenes",
			"util/apps/FindOverlappingGenes",
			"util/apps/FindSharedRegions",
			"edu/expr/FileCrossFilter",
			"edu/expr/FileMatchJoiner",
			"util/apps/FileJoiner",
			"util/apps/FileSplitter",
			"edu/utah/seq/parsers/FilterDuplicateAlignments",
			"edu/utah/seq/data/Graph2Bed",
			"util/apps/FilterIntersectingRegions",
			"edu/utah/seq/data/FilterPointData",
			"trans/misc/Gr2Bar",
			"util/apps/gwrap/GWrap_GUI_ClickMe",
			"util/apps/IntersectLists",
			"trans/anno/IntersectKeyWithRegions",
			"trans/anno/IntersectRegions",
			"edu/expr/KeggPathwayEnrichment",
			"util/bio/parsers/MaqSnps2Bed",
			"util/bio/seq/MakeSpliceJunctionFasta",
			"util/bio/seq/MakeTranscriptome",
			"util/bio/seq/MaskExonsInFastaFiles",
			"util/bio/seq/MaskRegionsInFastaFiles",
			"edu/utah/seq/parsers/MergePairedSamAlignments",
			"edu/utah/seq/data/MergePointData",
			"util/apps/MergeRegions",
			"util/bio/parsers/MergeUCSCGeneTable",
			"edu/utah/seq/analysis/multi/MultipleConditionRNASeq",
			"edu/utah/seq/analysis/MultipleReplicaDefinedRegionScanSeqs",
			"edu/utah/seq/analysis/MultipleReplicaScanSeqs",
			"edu/utah/seq/parsers/NovoalignBisulfiteParser",
			"edu/utah/seq/parsers/NovoalignIndelParser",
			"edu/utah/seq/parsers/NovoalignParser",
			"edu/utah/seq/parsers/NovoalignPairedParser",
			"util/bio/seq/OligoTiler",
			"edu/utah/seq/analysis/OverdispersedRegionScanSeqs",
			"edu/utah/seq/parsers/ParseIntersectingAlignments",
			"edu/utah/seq/data/ParsePointDataContexts",
			"edu/utah/seq/analysis/PeakShiftFinder",
			"edu/utah/seq/data/PointDataManipulator",
			"util/bio/wrappers/Primer3Wrapper",
			"util/apps/PrintSelectColumns",
			"edu/utah/seq/analysis/QCSeqs",
			"edu/utah/seq/data/Qseq2Fastq",
			"util/apps/RandomizeTextFile",
			"trans/graphics/RankedSetAnalysis",
			"edu/utah/seq/data/ReadCoverage",
			"edu/utah/seq/parsers/RNAEditingPileUpParser",
			"edu/utah/seq/analysis/RNAEditingScanSeqs",
			"edu/utah/seq/analysis/RNASeq",
			"edu/utah/seq/data/RNASeqSimulator",
			"edu/utah/seq/parsers/Sam2Fastq",
			"edu/utah/seq/data/Sam2USeq",
			"edu/utah/seq/base/SamAlignmentExtractor",
			"edu/utah/seq/parsers/SamParser",
			"edu/utah/seq/parsers/SamTranscriptomeParser",
			"edu/utah/seq/parsers/SamFixer",
			"edu/utah/seq/analysis/ScanSeqs",
			"edu/utah/seq/analysis/ScoreMethylatedRegions",
			"util/bio/parsers/ShiftAnnotationPositions",
			"edu/utah/seq/parsers/SoapV1Parser",
			"util/apps/SubtractRegions",
			"util/apps/ScoreChromosomes",
			"trans/roc/ScoreParsedBars",
			"util/apps/ScoreSequences",
			"trans/misc/Sgr2Bar",
			"edu/utah/seq/data/Simulator",
			"edu/utah/seq/base/SNPComparator",
			"edu/utah/seq/analysis/StrandedBisSeq",
			"util/bio/parsers/SRAProcessor",
			"edu/utah/seq/data/SubSamplePointData",
			"edu/utah/seq/parsers/Tag2Point",
			"edu/utah/seq/useq/apps/Text2USeq",
			"edu/utah/seq/useq/apps/UCSCBig2USeq",
			"edu/utah/seq/useq/apps/USeq2UCSCBig",
			"edu/utah/seq/useq/apps/USeq2Text",
			"trans/misc/Wig2Bar",
			"edu/utah/seq/useq/apps/Wig2USeq",
	};

	//constructor
	public MakeUSeq (String[] args){
		version = args[0];
		if (args.length > 1) makeWebStart = true;

		//make bioToolsCodeLibrary.jar		
		makeLibrary();

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
		for (int i=0; i< uSeqClassesToJar.length; i++){
			//check if class file exists
			File classFile = new File(uSeqClassesToJar[i]+".class");
			if (classFile.exists() == false) Misc.printExit("\nError: class file does not exist "+classFile);
			//make and write manifest
			String path;
			if (uSeqClassesToJar[i].contains("GWrap_GUI_ClickMe")) path = "\nClass-Path: LibraryJars/bioToolsCodeLibrary.jar\n";
			else path = "\nClass-Path: ../LibraryJars/bioToolsCodeLibrary.jar\n";
			String manifest = "Manifest-Version: 1.0\nMain-Class: " + uSeqClassesToJar[i]+ path+"\nImplementation-Version: USeq_"+version;
			manFile = new File ("MANIFEST.MF");
			manFile.deleteOnExit();
			IO.writeString(manifest, manFile);
			//make dummy jar
			String[] tokens = uSeqClassesToJar[i].split("/");
			String[] jarCmd = {"jar", "cmf0", manFile.getName(), tokens[tokens.length-1]};
			String[] results = IO.executeCommandLineReturnAll(jarCmd);

			if (results.length != 0){
				Misc.printArray(results);
				Misc.printExit("\nError: problem with making dummy jar for "+uSeqClassesToJar[i]);
			}
			//move into Apps
			File jarApp = new File (tokens[tokens.length-1]);
			File jarAppInDir = new File (jarAppDir, jarApp.getName());
			jarApp.renameTo(jarAppInDir); 
		}
	}

	/**Generates the cmdLnMenus.html doc.*/
	public void makeMenuHTMLDoc(){
		StringBuffer names = new StringBuffer();
		StringBuffer menus = new StringBuffer();

		for (int i=0; i< uSeqClassesToJar.length; i++){
			//skip GWrap_GUI_ClickMe
			if (uSeqClassesToJar[i].contains("GWrap_GUI_ClickMe")) continue;
			String[] cmd = {"java", uSeqClassesToJar[i]};
			String[] menu = IO.executeCommandLineReturnAll(cmd);
			if (menu[0].startsWith("Exception")) {
				Misc.printArray(menu);
				Misc.printExit("\nError: problem lauching "+uSeqClassesToJar[i]);
			}
			String[] tokens = uSeqClassesToJar[i].split("/");
			String jarName = tokens[tokens.length-1];
			names.append("<a href=\"#"+  jarName +"\">"+ jarName +"</a><br>\n");
			menus.append("<a name=\""+ jarName +"\"><pre>\n");
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
		File menuFile = new File ("../Documentation/USeqDocumentation/cmdLnMenus.html");
		IO.writeString(sb.toString(), menuFile); 
	}

	/**Makes the code library jar from the directoriesForLibJar array.
	 * Saves as bioToolsCodeLibrary.jar*/
	public boolean makeLibrary(){
		
		//copy over and uncompress com, org, info, etc
		String cmd = "rm -rf com info javax net org; cp ../Misc/JavaxComOrgInfoNetClasses.zip . ; unzip -qo JavaxComOrgInfoNetClasses.zip ; mv -f JavaxComOrgInfoNetClasses/* . ; rm -rf JavaxComOrgInfoNetClasses* __MACOSX";
		String[] errors = IO.executeShellScript(cmd, new File("."));
		for (String error: errors) if (error.length()!=0) Misc.printErrAndExit("\nProblem copying and uncompressing JavaxComOrgInfoNetClasses.zip. Aborting! -> \n"+Misc.stringArrayToString(errors, "\n"));

		//copy over class files
		jarAppDir = new File ("Apps");
		if (jarAppDir.exists()) IO.deleteDirectory(jarAppDir);
		jarAppDir.mkdir();
		for (int i=0; i< uSeqDirForLibJar.length; i++){
			File source = new File (uSeqDirForLibJar[i]);
			if (source.exists()== false) Misc.printExit("\nError: cannot find source directory "+source);
			if (IO.copyDirectoryRecursive(source, new File (jarAppDir,source.getName()),"class$|properties$") == false) Misc.printExit("\nError: directory copy failed for "+source);
		}
		String dirName = "USeq_"+version;
		
		//make manifiest file
		String manifest = "Manifest-Version: 1.0\nImplementation-Vendor: University of Utah Bioinformatics Shared Resource (http://bioserver.hci.utah.edu)\nImplementation-Version: USeq_"+version+"\n";
		File manFile = new File (jarAppDir, "MANIFEST.MF");
		IO.writeString(manifest, manFile);

		//jar it		
		StringBuffer sb = new StringBuffer();
		sb.append("\n**** Type on the command line****\n\n");
		
		//standard sourceforge
			sb.append("cd "+IO.getFullPathName(jarAppDir)+"\n");
			sb.append("jar cmf0 MANIFEST.MF bioToolsCodeLibrary.jar ");
			sb.append(Misc.stringArrayToString(uSeqDirForLibJar," "));
			sb.append("\nrm -rf MANIFEST.MF \n");
			sb.append("cd ..\nrm -rf ");
			for (int i=0; i<uSeqDirForLibJar.length; i++){
				sb.append(" Apps/");
				sb.append(uSeqDirForLibJar[i]);
			}

			sb.append(
					"\nrm -rf USeq_*\n"+
					"mkdir "+dirName+"\n"+
					"mkdir "+dirName+"/Documentation "+dirName+"/LibraryJars \n"+ // USeq/OthersCode \n"+
					"touch "+dirName+ "/Documentation/version_"+version+"\n"+
					"zip -r -q "+dirName+"/SourceCode.zip ../Source\n"+
					"cp -R ../Documentation/USeqDocumentation/ "+dirName+"/Documentation/\n"+
					"mv -f Apps/bioToolsCodeLibrary.jar "+dirName+"/LibraryJars/\n"+
					"mv -f Apps/GWrap_GUI_ClickMe "+dirName+"/GWrap_GUI_ClickMe.jar\n"+
					"mv -f Apps/ "+dirName+"/\n" +
					"zip -r -q "+dirName+".zip "+dirName+"\n"+
					"rm -rf ../Releases/"+dirName+"* \n"+
					"mv -f "+dirName+"* ../Releases/\n\n"+
					
					"#Nix specific configurations\n"+
					"rm -rf ~/AppsUSeq\n"+
					"ln -s /Users/davidnix/Code/USeq/Releases/"+dirName+"/Apps/ ~/AppsUSeq\n"+
					"rsync -r /Users/davidnix/Code/USeq/Releases/"+dirName+"/* $moab:/home/BioApps/USeq/\n\n");
		

		System.out.println(sb);
		return true;
	}

	public static void main (String[] args){
		if (args.length==0) Misc.printErrAndExit("\nPlease enter a version number, ie '3.7'\n");
		new MakeUSeq(args);
	}

	public final String menuHeader = "<html><head><title>USeq Command Line Menus "+version+"</title>" +
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
