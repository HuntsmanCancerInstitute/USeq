package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**
 * @author Nix
 * */
public class ROCPlotsMulti {

	//user defined fields
	private File[] dirsToParse;
	private File[] toParse;
	private HashMap<String, FdrSet> fdrSets = new HashMap<String, FdrSet>();
	private File pdfFile = null;
	private PrintWriter out = null;
	private int numRows = 2;
	private int numColumns = 3;
	private File fullPathToR = new File ("/usr/local/bin/R");

	//constructors
	/**Stand alone.*/
	public ROCPlotsMulti(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//for each dir, say 0.1, 0.05, 0.001
		for (File d: dirsToParse) {
			IO.pl("\nProcessing "+d);
			toParse = IO.extractFiles(d, ".txt");
			if (toParse.length > (numRows * numColumns)) Misc.printErrAndExit("Adjust numRows or numColumns");
			
			//parse the FdrSets for each of the comps fdrTprSummary.txt
			fdrSets.clear();
			for (File f: toParse) {
				IO.pl("\tParsing "+f);
				parse(f);
			}

			//print script file to generate 5 plots, one for each caller app, Strelka, ssc, lofreq, ...
			try {
				File scriptFile = new File(d, "rScript.txt");
				out = new PrintWriter(new FileWriter(scriptFile));
				pdfFile = new File(d, d.getParentFile().getName()+"_"+d.getName()+".pdf");
				//start printer and set grid
				out.println("pdf(file='"+pdfFile+"', width=12, height=7)");
				out.println("old.par = par(mfrow=c("+numRows+", "+numColumns+"))");

				//walk hash map
				for (String setName: fdrSets.keySet()) {
					//IO.pl(setName);
					FdrSet fs = fdrSets.get(setName);
					fs.printPlot(setName);
				}
				//close the grid and printer
				out.println("par(old.par)");
				out.println("dev.off()");
				out.close();

				//make command
				File rOut = new File(d, "rOutput.txt");
				String[] command = new String[] {
						fullPathToR.getCanonicalPath(),
						"CMD",
						"BATCH",
						"--no-save",
						"--no-restore",
						scriptFile.getCanonicalPath(),
						rOut.getCanonicalPath()};			
				//execute
				//Misc.printArray(command);
				if (IO.executeCommandLine(command) == null) Misc.printErrAndExit("\nError executing R script.\n"+scriptFile);
				String[] output = IO.loadFile(rOut);
				for (String s: output) if (s.toLowerCase().contains("error")) Misc.printErrAndExit("\nError executing R script.\n"+scriptFile);
				
				
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}




	private void parse(File f) {
		BufferedReader in;
		try {
			in = IO.fetchBufferedReader(f);

			//pull first line
			//dFDR	12_Snvs_0.1.lofreq	dFDR	12_Snvs_0.1.mutect	dFDR	12_Snvs_0.1.ssc	dFDR	12_Snvs_0.1.strelka	dFDR	12_Snvs_0.1.tnscope
			//dFDR	15352X1_Indels_0.1_Lofreq_low	dFDR	15352X1_Indels_0.1_Mutect_low	dFDR	15352X1_Indels_0.1_SSC_low	dFDR	15352X1_Indels_0.1_Strelka_low	dFDR	15352X1_Indels_0.1_TNscope_low


			//fix graph names, remove 1534X? _low _moderate?
			String hl = in.readLine();
			hl = hl.replaceAll("15352X", "");
			//hl = hl.replaceAll("_mod", "");

			String[] header = Misc.TAB.split(hl);

			//make ArrayLists
			ArrayList<Float>[] fAL = new ArrayList[header.length];
			for (int i=0; i< header.length; i++) fAL[i] = new ArrayList<Float>();

			//for each line parse values
			String line;
			String[] values;
			while ((line = in.readLine())!= null) {
				values = Misc.TAB.split(line, -1);
				if (values.length != header.length) Misc.printErrAndExit("Mismatche in data column length \n\t"+line);
				for (int i=0; i< header.length; i++) {
					if (values[i].length()!=0) {
						fAL[i].add(Float.parseFloat(values[i]));
					}
				}
			}

			//load hash
			for (int i=0; i< header.length; i++) {
				ArrayList<Float> fdr = fAL[i];
				i++;
				ArrayList<Float> tpr = fAL[i];

				//parse set name, 12_Snvs_0.1.lofreq to Snvs_0.1.lofreq
				String name = header[i];
				String setName = name.substring(name.indexOf("_")+1);

				FdrSet fdrSet = fdrSets.get(setName);
				if (fdrSet == null) {
					fdrSet = new FdrSet();
					fdrSets.put(setName, fdrSet);
				}
				fdrSet.add(name, fdr, tpr);
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}



	}

	private class FdrSet {
		ArrayList<String>  individualNames = new ArrayList<String>();
		ArrayList<float[]> fdr = new ArrayList<float[]>();
		ArrayList<float[]> tpr = new ArrayList<float[]>();

		public void add(String name, ArrayList<Float> fdrAL, ArrayList<Float> tprAL ) {
			individualNames.add(Misc.UNDERSCORE.split(name)[0]);
			fdr.add(Num.arrayListOfFloatToArray(fdrAL));
			tpr.add(Num.arrayListOfFloatToArray(tprAL));
		}

		public void print() {
			for (int i=0; i< individualNames.size(); i++) {
				IO.pl("\t"+individualNames.get(i));
				IO.pl("\tFDR "+Misc.floatArrayToString(fdr.get(i), ","));
				IO.pl("\tTPR "+Misc.floatArrayToString(tpr.get(i), ","));
			}
		}


		String[] colors = {"#FE2712", "#FB9902",  "#FEFE33", "#B2D732", "#66B032", "#347C98", "#0247FE", "#FCCC1A", "#4424D6", "#8601AF", "#C21460", "#FC600A"};
		int[] symbols = {15,16,17,18,19,20,21,22,23,24,25,1,2,3,4,5,6};

		public void printPlot(String title) throws IOException {
			TreeMap<Integer, String[]> ordered = fetchOrderedData();
			Set<Integer> nameNumber = ordered.keySet();
			ArrayList<Integer> al = new ArrayList<Integer>();
			al.addAll(nameNumber);

			//save for legend
			ArrayList<String> lNames = new ArrayList<String>();
			ArrayList<String> lColors = new ArrayList<String>();
			ArrayList<String> lSymbols = new ArrayList<String>();

			//find first good one
			String[] first = null;
			int i=0;
			for (; i< al.size(); i++) {
				String[] data = ordered.get(al.get(i));
				lColors.add(colors[i]);
				lSymbols.add(symbols[i]+"");
				if (data != null) {
					first = data;
					lNames.add(al.get(i)+"");
					break;
				}
				else lNames.add(al.get(i)+" NA");
			}

			if (first == null) IO.pl("\t\t# "+title+ " All Failed!");
			else {
				//print out data
				for (String[] d: ordered.values()) if(d != null) out.println(d[0]);

				//print plot
				String parsedTitle = title.replaceAll("_", " ");
				out.println("plot("+first[1]+",pch="+symbols[i]+",col='"+colors[i]+"',main='"+parsedTitle+"',type='b',bty='l',xlab='dFDR',ylab='TPR',lwd=3,xlim=c(0,1),ylim=c(0,1))");



				//print lines
				i++;
				for (; i< al.size(); i++) {
					String[] data = ordered.get(al.get(i));
					lColors.add(colors[i]);
					lSymbols.add(symbols[i]+"");
					if (data != null) {
						out.println("lines("+data[1]+",pch="+symbols[i]+",col='"+colors[i]+ "',lwd=3,type='b')" );
						lNames.add(al.get(i)+"");
					}
					else lNames.add(al.get(i)+" NA");
				}

				//print grid
				out.println("abline(h=seq(0,1,0.1),v=seq(0,1,0.1),col='grey',lwd=0.8)");

				//print the legend
				out.println("legend('bottomright',bty ='n',pt.cex=1,cex=1,text.col='black',horiz=F,inset=c(0.05, 0.05),");
				out.println("legend=c('"+ Misc.stringArrayListToString(lNames, "','")+ "'),");
				out.println("col=c('"+ Misc.stringArrayListToString(lColors, "','")+ "'),");
				out.println("pch=c("+ Misc.stringArrayListToString(lSymbols, ",")+ ") )");
			}
		}

		public TreeMap<Integer, String[]> fetchOrderedData() {
			TreeMap<Integer, String[]> ordered = new TreeMap<Integer, String[]>();
			for (int i=0; i< individualNames.size(); i++) {
				int numberName = Integer.parseInt(individualNames.get(i));
				//is it a failed array
				float[] iFdr = fdr.get(i);
				if (iFdr.length==3 && Num.sumArray(iFdr)==3.0f) ordered.put(numberName, null);
				else {
					String[] s = {"d"+numberName+ "Fdr=c("+ Misc.floatArrayToString(fdr.get(i), ",")+ ")\n"+
							"d"+numberName+ "Tpr=c("+ Misc.floatArrayToString(tpr.get(i), ",")+ ")", "d"+numberName+"Tpr~d"+numberName+"Fdr"};
					ordered.put(numberName, s);
				}
			}
			return ordered;
		}

		public double numEmptySets() {
			double num = 0;
			for (int i=0; i< individualNames.size(); i++) {
				float[] iFdr = fdr.get(i);
				if (iFdr.length==3 && Num.sumArray(iFdr)==3.0f) num++;
			}

			return num/(double)individualNames.size();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ROCPlotsMulti(args);
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
					case 'd': dirsToParse = IO.fetchDirsRecursively(new File(args[++i])); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group()+", Try -h");
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		Arrays.sort(dirsToParse);
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Sam 2 Fastq: May 2018                              **\n" +
				"**************************************************************************************\n" +

				"Example: java -Xmx2G -jar pathTo/USeq/Apps/Sam2Fastq -a myQNSorted.bam -s S2F/\n\n"+

				"**************************************************************************************\n");

	}


}
