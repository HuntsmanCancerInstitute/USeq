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
 * Generates a single combine ROC plot for each directory
 * @author Nix
 * */
public class ROCPlotsSingle {

	//user defined fields
	private File[] dirsToParse;
	private File[] toParse;
	private File nameSwitcher;
	private HashMap<String, String> renamer = null;
	private FdrSet fdrSet = null;
	private File pdfFile = null;
	private PrintWriter out = null;
	private File fullPathToR = new File ("/usr/local/bin/R");
	private double xMax = 0.3; //or 1
	private int numRows = 2;
	private int numColumns = 3;

	//constructors
	/**Stand alone.*/
	public ROCPlotsSingle(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//check num rows and columns
		if (dirsToParse.length > (numRows * numColumns)) Misc.printErrAndExit("Adjust numRows or numColumns");

		try {

			//start printer and set grid for composite pdf
			File scriptFile = new File(dirsToParse[0].getParentFile(), "rScript.txt");
			out = new PrintWriter(new FileWriter(scriptFile));
			pdfFile = new File(dirsToParse[0].getParentFile(), dirsToParse[0].getParentFile().getName()+"_pROC.pdf");
			out.println("pdf(file='"+pdfFile+"', width=12, height=7)");
			out.println("old.par = par(mfrow=c("+numRows+", "+numColumns+"))");


			//for each dir, say 0.1, 0.05, 0.001
			for (int i=dirsToParse.length-1; i>=0; i--) {
				File d = dirsToParse[i];
				IO.pl("\nProcessing "+d);
				toParse = IO.extractFiles(d, ".txt");

				//parse the FdrSets for each of the comps fdrTprSummary.txt, combining all
				fdrSet = new FdrSet();
				for (File f: toParse) {
					IO.pl("\tParsing "+f);
					parse(f);
				}

				//print script file to generate 1

				/*File scriptFile = new File(d, "rScript.single");
				//scriptFile.deleteOnExit();
				//out = new PrintWriter(new FileWriter(scriptFile));
				//pdfFile = new File(d, d.getParentFile().getName()+"_"+d.getName()+".single.pdf");
				//start printer
				//out.println("pdf(file='"+pdfFile+"', width=12, height=7)");
				 */

				//define title and print it
				String title = d.getParentFile().getName()+"_"+d.getName();
				fdrSet.printPlot(title);

				//close the printer
				/*out.println("dev.off()");
				//out.close();

				//make command
				File rOut = new File(d, "rOutput.single");
				rOut.deleteOnExit();
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
				 */



			}
			
			//close the grid and printer
			out.println("par(old.par)");
			out.println("dev.off()");
			out.close();
			
			//make command
			File rOut = new File(dirsToParse[0].getParentFile(), "rOutput.single");
			rOut.deleteOnExit();
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


			//fix all graph names, remove 1534X? 15352X1_Indels_0.0075.1st to 1_Indels_0.0075.1st
			String hl = in.readLine();
			//hl = hl.replaceAll("15352X", "");
			//hl = hl.replaceAll("Indels_", "");
			//hl = hl.replaceAll("Snvs_", "");

			String[] header = Misc.TAB.split(hl);
			//change the names?
			if (renamer !=null) {
				for (int i=1; i< header.length; i++) {
					String newName = renamer.get(header[i]);
					if (newName == null) throw new IOException("Couldn't rename "+header[i]);
					header[i] = newName;
					i++;
				}
			}

			//make ArrayLists
			ArrayList<Float>[] fAL = new ArrayList[header.length];
			for (int i=0; i< header.length; i++) fAL[i] = new ArrayList<Float>();

			//for each line parse values
			String line;
			String[] values;
			while ((line = in.readLine())!= null) {
				values = Misc.TAB.split(line, -1);
				if (values.length != header.length) Misc.printErrAndExit("Mismatch in data column length \n\t"+line);
				for (int i=0; i< header.length; i++) {
					if (values[i].length()!=0) {
						fAL[i].add(Float.parseFloat(values[i]));
					}
				}
			}

			//load fdrSet
			for (int i=0; i< header.length; i++) {
				ArrayList<Float> fdr = fAL[i];
				i++;
				ArrayList<Float> tpr = fAL[i];
				fdrSet.add(header[i], fdr, tpr);
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
			individualNames.add(name);
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


		String[] colors = {"#FE2712", "#FB9902",  "#FEFE33", "#B2D732", "#66B032", "#347C98", "#0247FE", "#FCCC1A", "#4424D6", "#8601AF", "#C21460", "#FC600A", "#000000", "#98f5ff"};
		int[] symbols = {15,16,17,18,19,20,21,22,23,24,25,1,2,3,4,5,6};

		public void printPlot(String title) throws IOException {
			TreeMap<String, String[]> ordered = fetchOrderedData();			
			Set<String> nameNumber = ordered.keySet();
			ArrayList<String> al = new ArrayList<String>();
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
				out.println("plot("+first[1]+",pch="+symbols[i]+",col='"+colors[i]+"',main='"+parsedTitle+
						"',type='b',bty='l',xlab='dFDR',ylab='TPR',lwd=3,xlim=c(0,"+xMax+"),ylim=c(0,1))");



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

		public TreeMap<String, String[]> fetchOrderedData() {
			TreeMap<String, String[]> ordered = new TreeMap<String, String[]>();
			for (int i=0; i< individualNames.size(); i++) {
				//is it a failed array
				float[] iFdr = fdr.get(i);
				if (iFdr.length==3 && Num.sumArray(iFdr)==3.0f) ordered.put(individualNames.get(i), null);
				else {
					String[] s = {"d"+i+ "Fdr=c("+ Misc.floatArrayToString(fdr.get(i), ",")+ ")\n"+
							"d"+i+ "Tpr=c("+ Misc.floatArrayToString(tpr.get(i), ",")+ ")", "d"+i+"Tpr~d"+i+"Fdr"};
					ordered.put(individualNames.get(i), s);
				}
			}
			return ordered;
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ROCPlotsSingle(args);
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
					case 'n': nameSwitcher = new File(args[++i]); break;
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
		if (nameSwitcher != null) renamer = IO.loadFileIntoHashMap(nameSwitcher);
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
