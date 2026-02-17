package edu.utah.seq.maf;

import java.io.*;
import java.util.ArrayList;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;


/**
 * Parses variant maf files collecting stats, fast and loose, get er done! 
 * @author davidnix*/
public class CBioFiltForWontak {

	//user defined fields
	private File[] mafDataFiles;
	private File outputDir;
	
	//internal
	PrintWriter out = null; 

	public CBioFiltForWontak(String[] args) throws IOException{
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		String[] pDots = {"p.G12R", "p.G12C", "p.G12D", "p.G12A", "p.G12V", "p.G12S", "p.G13D", "p.G13R", "p.G13V", "p.Q61R", "p.Q61K", "p.Q61L", "p.Q61H", "p.Q61E", "p.Q61P"};
		TreeSet<String> rasCanPDots = new TreeSet<String>();
		for (String pd: pDots) rasCanPDots.add(pd);
		IO.pl("RAS canonical pDots to find "+rasCanPDots);
		
		pDots = new String[]{"p.K117N", "p.A146T", "p.A146V", "p.L19F","p.T158A"}; 
		TreeSet<String> rasNonCanPDots = new TreeSet<String>();
		for (String pd: pDots) rasNonCanPDots.add(pd);
		IO.pl("RAS non canonical pDots to find "+rasNonCanPDots+"\n");
		
		out = new PrintWriter (new FileWriter(new File(outputDir, "cBioStats.xls")));
		boolean printHeader = true;
		
		File tempDir = new File(outputDir,"TmpDelme");
		tempDir.mkdirs();
		File r = new File("/usr/local/bin/R");
		IntersectListsHypergeometric ih = new IntersectListsHypergeometric(tempDir, r);
		
		for (int i=0; i< mafDataFiles.length; i++) {
			String fileName = Misc.removeExtension(mafDataFiles[i].getName());
			IO.pl(fileName);
			
			//make container for results
			CBioProjStatRes res = new CBioProjStatRes(fileName);
			
			//make parser
			MafParser2 mafParser = new MafParser2(mafDataFiles[i]);
			
			//number samples
			int numSamples = mafParser.fetchUniqueColumnValues("Tumor_Sample_Barcode").size();
			IO.pl((int)numSamples+"\tNum samples");
			res.setNumberSamples(numSamples);
			
			//fetch lines with MAP2K1 or MAP2K2
			TreeSet<String> map2k = new TreeSet<String>();
			map2k.add("MAP2K1");
			map2k.add("MAP2K2");
			ArrayList<String[]> map2kLines = mafParser.fetchGeneLines(map2k);
			res.setMafMAP2KLines(map2kLines);
			TreeSet<String> sampleIdsWithMap2k = mafParser.fetchUniqueColumnValues(map2kLines, "Tumor_Sample_Barcode");
			IO.pl(sampleIdsWithMap2k.size()+"\tNum samples with MAP2K1/2 non synonymous mutations");
			res.setNumberSamplesMAP2K(sampleIdsWithMap2k.size());
			TreeSet<String> mapPDots = mafParser.fetchUniqueColumnValues(map2kLines, "HGVSp_Short");
			res.setObservedMAP2KpDots(mapPDots);
			IO.pl("\tObserved MAP2K1/2 pDots\t"+mapPDots);
			
			RasStats kRS = intersectWithRas("KRAS", sampleIdsWithMap2k, mafParser, rasCanPDots, rasNonCanPDots);
			RasStats nRS = intersectWithRas("NRAS", sampleIdsWithMap2k, mafParser, rasCanPDots, rasNonCanPDots);
			RasStats hRS = intersectWithRas("HRAS", sampleIdsWithMap2k, mafParser, rasCanPDots, rasNonCanPDots);
			res.getRasStats().add(kRS);
			res.getRasStats().add(nRS);
			res.getRasStats().add(hRS);
			
			res.calculateIntersectionPValues(ih);
			
			if (printHeader) {
				out.println(res.getHeader());
				printHeader = false;
			}
			
			out.println(res.toString());
		}
		out.close();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		IO.pl("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	/*
	tree sets for pDots
	spreadsheet output
	oncoKB for muts
	pDots for Mapk
	*/
	
	private RasStats intersectWithRas(String ras, TreeSet<String> sampleIdsWithMap2k, MafParser2 mafParser, TreeSet<String> rasCanPDots, TreeSet<String> rasNonCanPDots) throws IOException {
		IO.pl("\n------- "+ras+" -------");
		RasStats rs = new RasStats(ras);
		
		//pull lines for the ras
		TreeSet<String> rasName = new TreeSet<String>();
		rasName.add(ras);
		ArrayList<String[]> rasLines = mafParser.fetchGeneLines(rasName);
		rs.setMafRasLines(rasLines);
		TreeSet<String> rasMutSamples = mafParser.fetchUniqueColumnValues(rasLines, "Tumor_Sample_Barcode");
		IO.pl(rasMutSamples.size()+"\tNum "+ras+" mutant samples, any");
		rs.setNumberSamplesRasAny(rasMutSamples.size());
		TreeSet<String> rasPDots = mafParser.fetchUniqueColumnValues(rasLines, "HGVSp_Short");
		IO.pl("\tObserved "+ras+" pDots\t"+rasPDots);
		rs.setObservedRaspDots(rasPDots);
		
		//filter for canonical
		ArrayList<String[]> rasLinesCanonical = mafParser.filterLines(rasLines, rasCanPDots, "HGVSp_Short");
		TreeSet<String> sampleIdsWithRasCan = mafParser.fetchUniqueColumnValues(rasLinesCanonical, "Tumor_Sample_Barcode");
		IO.pl(sampleIdsWithRasCan.size()+"\tNum "+ras+" mutant samples, canonical\n");
		rs.setNumberSamplesRasCanonical(sampleIdsWithRasCan.size());
		
		ArrayList<String> subCan = findIntersectingValues (sampleIdsWithRasCan, sampleIdsWithMap2k);
		IO.pl(subCan.size()+"\tNum MAP2K1/2 and "+ras+" canonical mut samples: "+subCan+"\n");
		rs.setNumberSamplesMAP2K_RasCanonical(subCan.size());
		
		//filter for non canonical
		ArrayList<String[]> rasLinesNonCanonical = mafParser.filterLines(rasLines, rasNonCanPDots, "HGVSp_Short");
		TreeSet<String> sampleIdsWithNonRasCan = mafParser.fetchUniqueColumnValues(rasLinesNonCanonical, "Tumor_Sample_Barcode");
		IO.pl(sampleIdsWithNonRasCan.size()+"\tNum "+ras+" mutant samples, non canonical\n");
		rs.setNumberSamplesRasNonCanonical(sampleIdsWithNonRasCan.size());
		ArrayList<String> subNon = findIntersectingValues (sampleIdsWithNonRasCan, sampleIdsWithMap2k);
		IO.pl(subNon.size()+"\tNum MAP2K1/2 and "+ras+" non canonical mut samples: "+subNon);
		rs.setNumberSamplesMAP2K_RasNonCanonical(subNon.size());
		
		return rs;
		
	}

	public ArrayList<String> findIntersectingValues(TreeSet<String> a, TreeSet<String> b){
		ArrayList<String> common = new ArrayList<String>();
		int sizeA = a.size();
		int sizeB = b.size();
		if (sizeA < sizeB) {
			for (String s: a) if (b.contains(s)) common.add(s);
		}
		else {
			for (String s: b) if (a.contains(s)) common.add(s);
		}
		return common;
	}
	
	public static void main(String[] args) throws IOException {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CBioFiltForWontak(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': forExtraction = new File(args[++i]); break;
					case 'o': outputDir = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a maf file or directory containing such to parse.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".maf");
		tot[1] = IO.extractFiles(forExtraction,".maf.gz");
		tot[2] = IO.extractFiles(forExtraction,".maf.zip");
		mafDataFiles = IO.collapseFileArray(tot);
		if (mafDataFiles == null || mafDataFiles.length ==0 || mafDataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find any xxx.maf(.txt/.gz OK) file(s)!\n");

		//outputDir
		if (outputDir == null) Misc.printExit("\nError: please provide an output directory for writing the fixed files.\n");
		outputDir.mkdirs();
		
		
	}	
	


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                               cBioFiltForWontak: Feb 2026                        **\n" +
				"**************************************************************************************\n" +
				"\n\n"+


		        "**************************************************************************************\n");
	}

}
