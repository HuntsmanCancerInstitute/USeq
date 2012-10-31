package util.apps;
import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.*;
import util.gen.*;
import util.bio.parsers.*;

public class FindSharedRegions {

	private File bedFile1 = null;
	private File bedFile2 = null;
	private File saveFile;
	private int minimumLength = 0;
	HashMap<String,RegionScoreText[]> regionsHash1;
	HashMap<String,RegionScoreText[]> regionsHash2;


	public  FindSharedRegions(String[] args) {

		processArgs(args);

		//for each chromosome
		try {
			PrintWriter bedOut = new PrintWriter( new FileWriter (saveFile));
			for (String chrom : regionsHash1.keySet()){

				//get regions 
				RegionScoreText[] sub1 = regionsHash1.get(chrom);
				RegionScoreText[] sub2 = regionsHash2.get(chrom);
				if (sub2 == null) continue;
				System.out.println("\t"+chrom);

				//make boolean of covered bps
				int lastBase1 = sub1[sub1.length -1].getStop();
				boolean[] coveredBps1 = new boolean[lastBase1 +1];
				Arrays.fill(coveredBps1, false);
				
				//for each region flip to false
				for (int i=0; i< sub1.length; i++){
					int start = sub1[i].getStart();
					if (start < 0) start = 0;
					int stop = sub1[i].getStop();
					for (int k=start; k< stop; k++) coveredBps1[k] = true;
				}
				
				boolean[] commonBPs = new boolean[ coveredBps1.length];
				Arrays.fill(commonBPs, true);
				
				//for each region flip to false
				for (int i=0; i< sub2.length; i++){
					int start = sub2[i].getStart();
					if (start < 0) start = 0;
					int stop = sub2[i].getStop();
					for (int k=start; k< stop; k++) {
						if (k > lastBase1) break;
						if (coveredBps1[k]) commonBPs[k] = false;
					}
				}
				
				//write out blocks
				int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(commonBPs, 0, minimumLength);
				for (int j=0; j< blocks.length; j++){
					bedOut.println(chrom+"\t"+ blocks[j][0]+ "\t"+ (blocks[j][1]+1));
				}
			}
			bedOut.close();
			System.out.println("\nDone!\n");
		} catch (Exception e){
			System.out.println();
			e.printStackTrace();
			System.exit(1);
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FindSharedRegions(args);
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
					case 'f': bedFile1 = new File (args[++i]); break;
					case 's': bedFile2 = new File (args[++i]); break;
					case 'r': saveFile = new File (args[++i]); break;
					case 'm': minimumLength = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check params
		if (bedFile1 == null || bedFile1.canRead() == false) Misc.printExit("\nCannot find or read your first bed file?!\n");
		if (bedFile2 == null || bedFile2.canRead() == false) Misc.printExit("\nCannot find or read your second bed file?!\n");
		if (saveFile == null) Misc.printExit("\nPlease enter a file name for saving your results.\n");

		//load regions, no strand
		regionsHash1 = Bed.parseBedFile(bedFile1, true);
		regionsHash2 = Bed.parseBedFile(bedFile2, true);

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Find Shared Regions: Dec 2011                         **\n" +
				"**************************************************************************************\n" +
				"Writes out a bed file of shared regions. Interbase coordinates.\n\n"+

				"Options:\n"+
				"-f First bed file (tab delimited: chr start stop ...).\n"+
				"-s Second bed file.\n"+
				"-r Results file.\n"+
				"-m Minimum length, defaults to 0.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/FindSharedRegions -f \n" +
				"      /Res/firstBedFile.bed -s /Res/secondBedFile.bed -r /Res/common.bed -m 100\n\n"+ 
				

		"************************************************************************************\n");

	}

}
