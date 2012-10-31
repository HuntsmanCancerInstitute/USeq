package trans.main;
import util.gen.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;

import util.gen.Misc;

/**Associates an empirical FDR estimate with each window provided a mock IP was performed.*/
public class FDRWindowConverter {
	//fields
	private File mockFile;
	private File[] ipFiles;
	private int scoreIndex = 0;
	private int numWindowsMock;
	private int numWindowsIP;
	private boolean printSgrs = false;
	private int sizeOfOligoMinusOne= 24;
	
	public FDRWindowConverter(String[] args){
		processArgs(args);
		WindowComparator comparator = new WindowComparator();
		
		//load mock Window[], set score, sort
		System.out.println("\nLoading mock Window[]...");
		Window[] mock = (Window[])IO.fetchObject(mockFile);
		numWindowsMock = mock.length;
		setScoreAndSort(mock);
		
		//for each ipFile
		for (int x=0; x<ipFiles.length; x++){
			//load windows set sortby and sort, big to small
			System.out.println("Loading "+ipFiles[x].getName()+" Window[]...");
			Window[] ip = (Window[])IO.fetchObject(ipFiles[x]);
			numWindowsIP = ip.length;
			
			setScoreAndSort(ip);
			//for each different window score calc FDR
			double score = -10000;
			int mockIndex = 0;
			double fdr = 0;
			double log10 = Math.log(10);
			double maxFdr = -1;
			for (int y=0; y< numWindowsIP; y++){
				//different score therefore recalculate fdr
				if (ip[y].getSortBy() != score){					
					//new score calculate FDR
					double numIPWins = y+1;					
					//find index of mock windows with a score < ip score
					mockIndex = thresholdWindows(mock, ip[y].getSortBy(), mockIndex);
					//watch out for mockIndex = numWindowsMock or zero
					if (mockIndex ==0 || mockIndex == numWindowsMock) fdr = 0;
					else {
						fdr = ((double)mockIndex) / numIPWins;
						if (fdr >= 1) fdr = 0.000001;
						else fdr = (-10.0) * (Math.log(fdr)/log10);
					}
					score = ip[y].getSortBy();
					if (fdr > maxFdr) maxFdr = fdr;
				}
				//set fdr
				double[] added = Num.appendDouble(ip[y].getScores(),fdr);
				ip[y].setScores(added);
			}

			//set max fdrs in non measurable windows as max
			int finalIndex = ip[0].getScores().length -1;
			for (int y=0; y< numWindowsIP; y++){
				double[] modScores = ip[y].getScores();
				if (modScores[finalIndex] ==0) modScores[finalIndex] = maxFdr;
			}
			
			//sort by chrom, pos, length
			Arrays.sort(ip, comparator);
			//save 
			IO.saveObject(new File(IO.getFullPathName(ipFiles[x])+"EmpFdr"), ip);
			
			//print sgr files?
			if (printSgrs){
				try{
				File pointSgrFile  = new File(ipFiles[x].getCanonicalPath()+"EmpFdr.sgr");
				PrintWriter out = new PrintWriter( new FileWriter (pointSgrFile));
				//for each window print sgr line
				for (int i=0; i< numWindowsIP; i++){
					double diff = sizeOfOligoMinusOne + ip[i].getStartLastOligo() - ip[i].getStart1stOligo();
					String pos = ip[i].getChromosome() +"\t"+ ( (int)Math.round(diff/2.0) + ip[i].getStart1stOligo() )+"\t";
					float trunkScore = new Double(ip[i].getScores()[finalIndex]).floatValue();
					out.println(pos+trunkScore);
				}
				out.close();
				IO.zipAndDelete(pointSgrFile);
				
				//print heat map
				WindowBlockMaker bm = new WindowBlockMaker(sizeOfOligoMinusOne +1);
				bm.setScoreIndex(finalIndex);
				bm.makeHeatMapSgrFile(ip, new File(ipFiles[x].getCanonicalPath()+"EmpFdrHM.sgr.zip"), true);
				
				} catch (IOException e){
					e.printStackTrace();
				}
			}
		}
		System.out.println("\nDone!\n");
	}
	
	public int thresholdWindows(Window[] w, double score, int startingIndex){
		for (int x=startingIndex; x< numWindowsMock; x++){			
			if (w[x].getSortBy() <= score) {				
				return x;
			}
		}		
		return numWindowsMock;
	}
	
	public void setScoreAndSort(Window[] w){
		int num = w.length;
		for (int y=0; y< num; y++){
			double[] scores = w[y].getScores();
			w[y].setSortBy((float)scores[scoreIndex]);
		}
		Arrays.sort(w);
	}
	
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'r': ipFiles = IO.extractFiles(new File(args[i+1])); i++; break;
					case 's': scoreIndex = Integer.parseInt(args[i+1]); i++; break;
					case 'z': sizeOfOligoMinusOne = Integer.parseInt(args[i+1]) - 1;i++; break;
					case 'm': mockFile = new File(args[i+1]); i++; break;
					case 'p': printSgrs = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group()+"\n");
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (mockFile == null || ipFiles[0] == null || mockFile.canRead()==false || ipFiles[0].canRead() == false){
			Misc.printExit("\nCannot find or read one of your Mock or IP Window[] files!\n");
		}
	}	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        FDR Window Converter: May 2006                            **\n" +
				"**************************************************************************************\n" +
				"Given a Window[] from a mock IP analysis (ie IgG) and Window[]s from real IPs, FDRWC\n" +
				"will associate an empirical FDR (# mock/ # real windows at each threshold) for each\n" +
				"real IP Window. The FDR will be appended on to each Window's score array as a\n" +
				"-10Log10(fdr) transformed value. \n\n" +
				
				"-m Full path file text for the mock IP serialized Window[] array.\n"+
				"-r Full path file text for a directory or file containing real IP serialized Window[]s.\n" +
				"-s Score index to use converting to FDRs, see ScanChip or ScanChromosome.\n" +
				"-z Size of oligo, defaults to 25.\n"+
				"-p Print sgr files for the empirical FDR scores.\n\n"+
				
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/FDRWindowConverter -m /affy/mockWins -r \n" +
				"      /affy/realIPs/ -s 1\n\n" +
				
		"**************************************************************************************\n");		
	}
	
	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);	
		}
		new FDRWindowConverter(args);
	}
	
}
