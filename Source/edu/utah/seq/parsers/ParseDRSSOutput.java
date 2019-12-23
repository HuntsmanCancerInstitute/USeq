
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import util.gen.*;
import java.util.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class ParseDRSSOutput{

	//user defined fields
	private File drssFile;
	private double treatmentCounts = 0;
	private double controlCounts = 0;
	
	

	//constructors
	public ParseDRSSOutput(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

						
			//load and sort
			String[] lines = IO.loadFile(drssFile);
			TreeMap<Integer, String[]> sortedLines = new TreeMap<Integer, String[]>();
			for (String l: lines){
				if (l.startsWith("#")) {
					//String[] t = Misc.TAB.split(l);
					//for (int i=0; i< t.length; i++){
					//	System.out.println(t[i]+"\t"+i);
					//}
					continue;
				}
				String[] cells = Misc.TAB.split(l);
				//rip number
				//"=HYPERLINK("http://localhost:7085/UnibrowControl?version=B37&seqid=17&start=56393245&end=56428242","5016")"
				String orderNumer = Misc.COMMA.split(cells[0])[1];
				orderNumer = orderNumer.substring(1, orderNumer.length()-2);
				sortedLines.put(Integer.parseInt(orderNumer), cells);	
			}
			
			//incoming #DisplayName Name Chr Strand Start Stop pVal qValFDR eFDR pValSkew pValDiffDist TotalRegionBPs Log2((sumT+1)/(sumC+1)) tSumPlus tSumMinus tRPKM cSumPlus cSumMinus cRPKM GenomeVersion=B37, TotalTreatObs=20150215, TotalCtrlObs=18248714
			//               0        1   2    3      4    5     6     7      8      9         10             11                 12              13        14      15      16       17      18        19
			String name = Misc.removeExtension(drssFile.getName());
			File outFile = new File(drssFile.getParentFile(), "parsed_"+name+".txt.gz");
			Gzipper out = new Gzipper(outFile);
			out.println("Chr\tStart\tStop\tqValFDR\tTotalRegionBPs\ttTotal\tcTotal\t"+name+"_Std_tRPKM\t"+name+"_Std_cRPKM\t"+name+"_Std_log2Rto\t"+name+"_Dmel_tRPKM\t"+name+"_Dmel_cRPKM\t"+name+"_Dmel_log2Rto");
			//outgoing 
			//Chr Start Stop qValFDR TotalRegionBPs tTotal cTotal tRPKM cRPKM log2Rto tRPKM cRPKM log2Rto
			for (String[] cells: sortedLines.values()){
				if (cells.length != 19) {
					IO.pl("Incorrect number of fields, skipping "+Misc.stringArrayToString(cells, "\t"));
					continue;
				}
				
				ArrayList<String> p = new ArrayList<String>();
				//Chr
				p.add(cells[2]);
				//Start
				p.add(cells[4]);
				//Stop
				p.add(cells[5]);
				//qValFDR
				p.add(cells[7]);
				//TotalRegionBPs
				p.add(cells[11]);
				//tTotal
				double tT = Double.parseDouble(cells[13]) + Double.parseDouble(cells[14]);
				p.add((int)tT+"");
				//cTotal
				double cT = Double.parseDouble(cells[16]) + Double.parseDouble(cells[17]);
				p.add((int)cT+"");
				//tRPKM - std
				p.add(cells[15]);
				//cRPKM - std
				p.add(cells[18]);
				//log2Rto - std
				double rto = Double.parseDouble(cells[15]) / Double.parseDouble(cells[18]);
				p.add( Num.formatNumber( Num.log2(rto), 4));
				//tRPKM - dmel
				double size = Double.parseDouble(cells[11]);
				double tR = (double)tT / size / treatmentCounts * 1000000000;
				p.add(Num.formatNumber(tR, 4));
				//cRPKM - dmel
				double cR = (double)cT / size / controlCounts * 1000000000;
				p.add(Num.formatNumber(cR, 4));
				//log2Rto - dmel
				double rtoDmel = tR/cR;
				double lgDmel = Num.log2(rtoDmel);
				p.add(Num.formatNumber(lgDmel, 4));
				out.println(Misc.stringArrayListToString(p, "\t"));
				
			}
			out.close();


			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			System.err.println("Prob!");
			e.printStackTrace();
		}
	}

	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ParseDRSSOutput(args);
	}		


	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		String dmelCountString = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': drssFile = new File(args[++i]); break;
					case 'd': dmelCountString = args[++i]; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//check
		if (drssFile == null || drssFile.exists() == false) Misc.printErrAndExit("\nFailed to find your drss file!");
		if (dmelCountString == null) Misc.printErrAndExit("\nFailed to find your treatment control count string!");
		int[] dmel = Num.parseInts(dmelCountString, Misc.COMMA);
		treatmentCounts = dmel[0];
		controlCounts = dmel[1];



	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                                **\n" +
				"**************************************************************************************\n" +
				

		"**************************************************************************************\n");

	}

}
