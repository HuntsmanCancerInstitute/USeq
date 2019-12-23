package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**
  * @author Nix
 * */
public class ParseFdrColumns {

	//user defined fields
	private File[] toParse;
	private HashMap<String, FdrSet> fdrSets = new HashMap<String, FdrSet>();
	private File pdfDirectory = null;
	private Gzipper out = null;

	//constructors
	/**Stand alone.*/
	public ParseFdrColumns(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		for (File f: toParse) {
			IO.pl("Parsing "+f);
			parse(f);
		}

		try {
			out = new Gzipper(new File(pdfDirectory, "rScript.txt"));
			//walk hash map
			for (String setName: fdrSets.keySet()) {
				//IO.pl(setName);
				FdrSet fs = fdrSets.get(setName);
				//fs.print();
				//IO.pl("\t"+fs.numEmptySets());
				fs.printPlot(setName, pdfDirectory);
			}
			out.close();
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
			String[] header = Misc.TAB.split(in.readLine());
			
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
		
		public void printPlot(String title, File saveDir) throws IOException {
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
			
			if (first == null) out.println("# "+title+ " All Failed!\n");
			else {
				out.println("pdf('"+saveDir+"/"+title+".pdf')");
				//print out data
				for (String[] d: ordered.values()) if(d != null) out.println(d[0]);
				
				//print plot
				out.println("plot("+first[1]+",pch="+symbols[i]+",col='"+colors[i]+"',main='"+title+"',type='b',bty='l',xlab='dFDR',ylab='TPR',lwd=3,xlim=c(0,1),ylim=c(0,1))");
				
				
				
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
				out.println("legend('bottomright',bty ='n',pt.cex=2,cex =1.2,text.col='black',horiz=F,inset=c(0.1, 0.1),");
				out.println("legend=c('"+ Misc.stringArrayListToString(lNames, "','")+ "'),");
				out.println("col=c('"+ Misc.stringArrayListToString(lColors, "','")+ "'),");
				out.println("pch=c("+ Misc.stringArrayListToString(lSymbols, ",")+ ") )");
				out.println("dev.off()\n");
				
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
		new ParseFdrColumns(args);
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
					case 'd': toParse = IO.fetchFilesRecursively(new File(args[++i]), ".txt"); break;
					case 'p': pdfDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group()+", Try -h");
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (pdfDirectory.exists() == false) pdfDirectory.mkdirs();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Sam 2 Fastq: May 2018                              **\n" +
				"**************************************************************************************\n" +
				"Given a query name sorted alignment file, S2F writes out fastq data for paired and \n"+
				"unpaired alignemnts. Any non-primary, secondary, or supplemental \n"+
				"alignments are written to a failed sam file.  This app doesn't have the memory leak\n"+
				"found in Picard, writes gzipped fastq, and error checks the reads. Provide an\n"+
				"unfiltered fastq-bam to use in retrieving missing mates. S2F will remove pre and\n"+
				"post naming info from the BamBlaster restoring the original fragment names.\n"+

				"\nOptions:\n"+
				"-s Path to a directory for saving parsed data.\n"+
				"-a Path to a query name sorted bam/sam alignment file. \n"+
				"-u (Optional) Path to a query name sorted unfiltered bam/sam alignment file for use\n"+
				"      in fetching missing mates of the first bam. Convert fastq to bam, then qn sort.\n"+

				"\n"+

				"Example: java -Xmx2G -jar pathTo/USeq/Apps/Sam2Fastq -a myQNSorted.bam -s S2F/\n\n"+

				"**************************************************************************************\n");

	}


}
