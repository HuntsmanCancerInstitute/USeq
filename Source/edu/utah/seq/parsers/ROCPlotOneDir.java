package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**
 * Generates a single combine ROC plot for the directory
 * @author Nix
 * */
public class ROCPlotOneDir {

	//user defined fields
	private File dirToParse = new File ("/Users/u0028003/HCI/Labs/Underhill/Hg38cfDNASims/NewTestTrainAnalysis/Replicas/VCFComp/ForRoc/Snv/Snv_0.05");
	private File[] toParse;
	private FdrSets fdrSets = null;
	private File pdfFile = null;
	private PrintWriter out = null;
	private File fullPathToR = new File ("/usr/local/bin/R");
	private double xMax = 0.15; //or 1
	private double yMin = 0.4; //or 0
	private double yMax = 1; // or 1
	//private double xMax = 1; //or 1
	//private double yMin = 0; //or 0
	//private double yMax = 0.25; // or 1
	private int numRows = 1;
	private int numColumns = 1;
	private String[] colors = {"#FE2712", "#FB9902",  "#FEFE33", "#B2D732", "#66B032", "#347C98", "#0247FE", "#FCCC1A", "#4424D6", "#8601AF", "#C21460", "#FC600A", "#000000", "#98f5ff",      "#FE2712", "#FB9902",  "#FEFE33", "#B2D732", "#66B032", "#347C98", "#0247FE", "#FCCC1A", "#4424D6", "#8601AF", "#C21460", "#FC600A", "#000000", "#98f5ff"};
	private int[] symbols = {15,16,17,18,19,20,21,22,23,24,25,1,2,3,4,5,6,      15,16,17,18,19,20,21,22,23,24,25,1,2,3,4,5,6};


	//constructors
	/**Stand alone.*/
	public ROCPlotOneDir(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		//processArgs(args);

		try {

			//pull files to parse
			IO.pl("\nProcessing "+dirToParse);
			toParse = IO.extractFiles(dirToParse, ".txt");

			//start printer and set grid for composite pdf
			File scriptFile = new File(dirToParse, "rScript.txt");
			scriptFile.deleteOnExit();
			out = new PrintWriter(new FileWriter(scriptFile));
			pdfFile = new File(dirToParse, dirToParse.getName()+"_pROC.pdf");
			out.println("pdf(file='"+pdfFile+"', width=12, height=7)");
			out.println("old.par = par(mfrow=c("+numRows+", "+numColumns+"))");

			//parse the FdrSets for each of the comps fdrTprSummary.txt, combining all
			fdrSets = new FdrSets();
			for (File f: toParse) {
				IO.pl("\tParsing "+f);
				parse(f);
			}
			
//HERE: define the title
			String title = dirToParse.getName();
			fdrSets.printPlot(title.replaceAll("_", " "));

			//close the grid and printer
			out.println("par(old.par)");
			out.println("dev.off()");
			out.close();

			//make command
			File rOut = new File(dirToParse, "rOutput.single");
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


//HERE: fix all graph names, remove 1534X? 15352X1_Indels_0.0075.1st to 1_Indels_0.0075.1st
			String hl = in.readLine();
			IO.pl("Ori "+hl);
			//Pattern pat = Pattern.compile("_Indels_0\\.\\d+");
			Pattern pat = Pattern.compile("_Snvs_0\\.\\d+");
			hl = pat.matcher(hl).replaceAll("");
			hl = hl.replaceAll("merged", "Merged");
			hl = hl.replaceAll("replica", "Replica");
			hl = hl.replaceAll("single", "Single");
			IO.pl("Pst "+hl+"\n");
			//hl = hl.replaceAll("15352X", "");
			//hl = hl.replaceAll("Indels_", "");
			//hl = hl.replaceAll("Snvs_", "");

			String[] header = Misc.TAB.split(hl);

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
				String legendName = header[i];
//HERE: modify legend name
String[] split = header[i].split("\\.");
legendName = split[1]+ " "+ split[0];

				fdrSets.addNewGraph(legendName, fdr, tpr);
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}



	}
	
	/*private class FdrSet {
		String legendName;
		float[] fdr;
		float[] tpr;
		String color;
		int symbol;
		
		public FdrSet(String name, float[] fdr, float[] tpr) {
			legendName = name;
			this.fdr = fdr;
			this.tpr = tpr;
		}
	}*/

	private class FdrSets {

		ArrayList<String>  individualNames = new ArrayList<String>();
		ArrayList<float[]> fdr = new ArrayList<float[]>();
		ArrayList<float[]> tpr = new ArrayList<float[]>();
		
		//ArrayList<FdrSet> fdrSetAL = new ArrayList<FdrSet>();

		public void addNewGraph(String name, ArrayList<Float> fdrAL, ArrayList<Float> tprAL) {
			individualNames.add(name);
			fdr.add(Num.arrayListOfFloatToArray(fdrAL));
			tpr.add(Num.arrayListOfFloatToArray(tprAL));
			//fdrSetAL.add(new FdrSet(name, Num.arrayListOfFloatToArray(fdrAL), Num.arrayListOfFloatToArray(tprAL)));
		}

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
			String color = null;
			for (; i< al.size(); i++) {
				String[] data = ordered.get(al.get(i));
				String legendName = al.get(i);
//HERE: modify graph colors, do the same thing below			
				if (legendName.startsWith("Merged")) color = "#FF0000";
				else if (legendName.startsWith("Replica")) color = "#008000";
				else if (legendName.startsWith("Single")) color = "#0000FF";
				else color = colors[i];
				lColors.add(color);
				
				lSymbols.add(symbols[i]+"");
				if (data != null) {
					first = data;
					lNames.add(legendName+"");
					break;
				}
				else lNames.add(legendName+" NA");
			}

			if (first == null) IO.pl("\t\t# "+title+ " All Failed!");
			else {
				//print out data
				for (String[] d: ordered.values()) if(d != null) out.println(d[0]);

				//print plot
				String parsedTitle = title.replaceAll("_", " ");
				out.println("plot("+first[1]+",pch="+symbols[i]+",col='"+color+"',main='"+parsedTitle+
						"',type='b',bty='l',xlab='dFDR',ylab='TPR',lwd=3,xlim=c(0,"+xMax+"),ylim=c("+yMin+","+yMax+"))");
				out.println("grid()");

				//print lines
				i++;
				for (; i< al.size(); i++) {
					String[] data = ordered.get(al.get(i));
					String legendName = al.get(i);
//HERE: modify graph colors				
					if (legendName.startsWith("Merged")) color = "#FF0000";
					else if (legendName.startsWith("Replica")) color = "#008000";
					else if (legendName.startsWith("Single")) color = "#0000FF";
					else color = colors[i];
					lColors.add(color);
					
					lSymbols.add(symbols[i]+"");
					if (data != null) {
						out.println("lines("+data[1]+",pch="+symbols[i]+",col='"+color+ "',lwd=3,type='b')" );
						lNames.add(legendName+"");
					}
					else lNames.add(legendName+" NA");
				}

				//print grid
				out.println("abline(h=seq(0,1,0.1),v=seq(0,1,0.1),col='grey',lwd=0.8)");

//HERE: where to print the legend bottomright topleft
				out.println("legend('bottomright',bty ='n',pt.cex=1,cex=1,text.col='black',horiz=F,inset=c(0.05, 0.05),");
				out.println("legend=c('"+ Misc.stringArrayListToString(lNames, "','")+ "'),");
				out.println("col=c('"+ Misc.stringArrayListToString(lColors, "','")+ "'),");
				out.println("pch=c("+ Misc.stringArrayListToString(lSymbols, ",")+ ") )");
			}
		}

		public TreeMap<String, String[]> fetchOrderedData() {
			TreeMap<String, String[]> ordered = new TreeMap<String, String[]>();
			//for each set
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
		new ROCPlotOneDir(args);
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
					case 'd': dirToParse = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group()+", Try -h");
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
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
