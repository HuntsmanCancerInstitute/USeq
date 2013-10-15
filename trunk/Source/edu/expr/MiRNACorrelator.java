package edu.expr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class MiRNACorrelator {
	
	//user defined fields
	private File miRNAData;
	private File allMiRNANames;
	private File targetData;
	private File geneExpressionData;
	private File results;
	
	//internal fields
	private static final Pattern TAB = Pattern.compile("\t");
	private HashSet<String> namesAllMiRNAs;
	/*nameMiRNA : valuesMiRNA*/
	private HashMap<String, MiRNA> miRNAs = new HashMap<String, MiRNA>();
	private String[] queriedMiRNANames;
	/*nameGene : geneExpressionValuesLg2FDR*/
	private HashMap<String, double[]> genes = new HashMap<String, double[]>();
	/*nameGene : setOfMiRNANames*/
	private HashMap<String, HashSet<String>> targetGeneNames2MiRNA = new HashMap<String, HashSet<String>>();
	private HashMap<String, HashSet<String>>  miRNANames2TargetNames = new HashMap<String, HashSet<String>>();
	
	
	//constructors
	public MiRNACorrelator(String[] args){
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		System.out.println("Loading gene expression data... "+geneExpressionData);
		parseGeneExpression();
		//System.out.println(genes);
		
		System.out.println("\nLoading all miRNA names... "+allMiRNANames);
		loadMiRNANames();
		//System.out.println(namesAllMiRNAs);
		
		System.out.println("\nLoading gene target 2 miRNA data... "+targetData);
		parseTargets();
		//System.out.println(targetGeneNames2MiRNA);
		//System.out.println(miRNANames2TargetNames);
		
		System.out.println("\nLoading select miRNA data... "+miRNAData);
		parseMiRNAs();
		//System.out.println(miRNAs);
		
		System.out.println("\nPrinting out merged spreadsheet... "+results);
		outputTable();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}
	
	private void loadMiRNANames() {
		namesAllMiRNAs = IO.loadFileIntoHashSet(allMiRNANames);
	}

	/*GeneName MiRNAName1 MiRNAName2 ...*/
	private void outputTable() {
		try {
			PrintWriter out = new PrintWriter( new FileWriter(results));
			
			//write out header
			out.print("GeneName");
			for (String mN: queriedMiRNANames){
				out.print("\t");
				out.print(mN);
			}
			out.println();
			
			//add miRNA log2Rtos
			out.print("Log2(2^(ctA-ctB))");
			for (String mN: queriedMiRNANames){
				out.print("\t");
				MiRNA m = miRNAs.get(mN);
				out.print(m.log2Rto);
			}
			out.println();
			
			//print gene lines
			//for each gene name
			for (String geneName: genes.keySet()){
				//is it ever a miRNA target?
				HashSet<String> miRNAsHittingThisGene = targetGeneNames2MiRNA.get(geneName);
				if (miRNAsHittingThisGene == null) continue;
				
				//were any of the miRNAs that hit this gene interrogated in this comparison?
				boolean skip = true;
				for (String mN: queriedMiRNANames){
					if (miRNAsHittingThisGene.contains(mN)){
						skip = false;
						break;
					}
				}
				if (skip) continue;
				
				double[] lg2FDR = genes.get(geneName);
				String val = "\t"+lg2FDR[0]+" "+lg2FDR[1];
				
				//print name
				out.print(geneName);
				//for each interrogated miRNA
				for (String mN: queriedMiRNANames){
					//does this miRNA target this gene?
					if (miRNAsHittingThisGene.contains(mN)){
						//print expression info
						out.print(val); 
					}
					else out.print("\t");
				}
				out.println();
			}
			
			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	/*Name NormCtA NormCtB Log2Rto*/
	private void parseMiRNAs() {
		String[] lines = IO.loadFileIntoStringArray(miRNAData);
		ArrayList<String> queriedMiRNANamesAL = new ArrayList<String>();
		//skip first header line
		for (int i=1; i< lines.length; i++) {
			MiRNA m = new MiRNA(lines[i]);
			//any targets?
			HashSet<String> targets = miRNANames2TargetNames.get(m.name);
			if (targets == null){
				System.out.println("\tSkipping "+m.name+", no gene targets");
				continue;
			}
			//are targets in gene expression data?
			boolean found = false;
			for (String geneName: targets){
				if (genes.containsKey(geneName)){
					found = true;
					break;
				}
			}
			
			if (found){
				miRNAs.put(m.name, m);
				queriedMiRNANamesAL.add(m.name);
			}
			else {
				System.out.println("\tSkipping "+m.name+", no gene targets in expression data");
			}
		}
		queriedMiRNANames = Misc.stringArrayListToStringArray(queriedMiRNANamesAL);
	}
	
	private class MiRNA{
		String name;
		double normCtA;
		double normCtB;
		double log2Rto;
		
		public MiRNA(String line){
			String[] t = TAB.split(line);
			name = t[0];
			if (namesAllMiRNAs.contains(name) == false){
				Misc.printErrAndExit("\nThis miRNA's name was not found in the all miRNAs file, aborting. -> "+name+"\n");
			}
			normCtA = Double.parseDouble(t[1]);
			normCtB = Double.parseDouble(t[2]);
			log2Rto = Double.parseDouble(t[3]);
		}
		
		public String toString(){
			return name+"\t"+normCtA+"\t"+normCtB+"\t"+log2Rto;
		}
	}
	
	/*GeneName Log2 FDR*/
	private void parseGeneExpression(){
		String[] lines = IO.loadFileIntoStringArray(geneExpressionData);
		//skip header line
		for (int i=1; i< lines.length; i++){
			String[] t = TAB.split(lines[i]);
			String geneName = t[0];
			double[] lgFdr = new double[]{Double.parseDouble(t[1]),Double.parseDouble(t[2])};
			if (genes.containsKey(geneName)) Misc.printErrAndExit("\nDuplicate gene name found! "+lines[i]+"\n");
			genes.put(geneName, lgFdr);
		}
	}

	/*TargetGeneName MiRNAName*/
	private void parseTargets() {
		String[] lines = IO.loadFileIntoStringArray(targetData);
		//skip first header line
		for (int i=1; i< lines.length; i++) {
			String[] t = TAB.split(lines[i]);
			String geneName = t[0];
			String miRNAName = t[1];
			//is the miRNA name in the all list?
			if (namesAllMiRNAs.contains(miRNAName) == false){
				Misc.printErrAndExit("\nThis miRNA's name was not found in the all miRNAs file, aborting. -> "+t[1]+"\n");
			}
			//load hash
			HashSet<String> m = targetGeneNames2MiRNA.get(geneName);
			if (m == null){
				m = new HashSet<String>();
				targetGeneNames2MiRNA.put(geneName, m);
			}
			m.add(miRNAName);
			
			//load hash
			HashSet<String> g = miRNANames2TargetNames.get(miRNAName);
			if (g == null){
				g = new HashSet<String>();
				miRNANames2TargetNames.put(miRNAName, g);
			}
			g.add(geneName);
		} 
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MiRNACorrelator(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': miRNAData = new File(args[++i]); break;
					case 'a': allMiRNANames = new File(args[++i]); break;
					case 't': targetData = new File(args[++i]); break;
					case 'e': geneExpressionData = new File(args[++i]); break;
					case 'r': results = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		if (miRNAData.exists() == false) Misc.printErrAndExit("\nMissing MiRNA data! "+miRNAData);
		if (allMiRNANames.exists() == false) Misc.printErrAndExit("\nMissing all miRNA name file! "+allMiRNANames);
		if (targetData.exists() == false) Misc.printErrAndExit("\nMissing gene target to miRNAs file! "+targetData);
		if (geneExpressionData.exists() == false) Misc.printErrAndExit("\nMissing gene expression data! "+geneExpressionData);
		if (results == null) Misc.printErrAndExit("\nMissing results file! "+results);
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             MiRNA Correlator: Oct 2013                           **\n" +
				"**************************************************************************************\n" +
				"Generates a spreadsheet to use in comparing changing miRNA levels to changes in gene\n"+
				"expression.\n"+

				"\nOptions:\n"+
				"-r Results file.\n"+
				"-a All miRNA name file (single column of miRNA names).\n"+
				"-m MiRNA data (two columns: miRNA name, miRNA log2Rto).\n"+
				"-t Gene target to miRNA data (two columns: gene target name, miRNA name).\n"+
				"-e Gene expression data (three columns: gene name, log2Rto, FDR).\n"+
				
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MiRNACorrelator -m miRNA_CLvsMOR.txt -a \n"+
				"allMiRNANamesNoPs.txt -t targetGene2MiRNA.txt -e geneExp_CLvsMOR.txt -r results.xls\n\n" +

		"**************************************************************************************\n");

	}


}
