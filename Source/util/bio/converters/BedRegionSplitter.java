package util.bio.converters;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.annotation.Bed;
import util.gen.*;

/**
 * For each region, splits it into xxx bp chunks if necessary. >
 */
public class BedRegionSplitter {

	private int chunkSize = 2000;
	private File[] toParse;
	private int bpPadding = 0;

	public BedRegionSplitter(String[] args) {
		processArgs(args);
		run();
		System.out.println("\nDone!\n");
	}

	public void run(){
		//for each file
		System.out.println("File\t#Start\t#End");
		for (int i=0; i< toParse.length; i++){
			if (toParse[i].isDirectory()) continue;
			parseIt(toParse[i]);
		}
	}

	private void parseIt(File file) {
		Gzipper out = null;
		File outputBed = null;
		try {
			//make a writer for the output
			String name = "_split.bed.gz";
			if (bpPadding !=0) name = "_splitPadded"+bpPadding+".bed.gz";
			outputBed = new File(file.getParentFile(), Misc.removeExtension(file.getName())+name);
			out = new Gzipper(outputBed);

			//load the bed regions
			Bed[] allRegions = Bed.parseFile(file, 0, 0);
			Arrays.sort(allRegions);
			HashMap<String, Bed[]> chromBed = Bed.splitBedByChrom(allRegions);

			int numStartingRegions = 0;
			int numFinalRegions = 0;

			//for each chromosome
			IO.p("\t");
			for (String chr: chromBed.keySet()) {
				IO.p(chr+" ");
				ArrayList<Bed> splitRegions = new ArrayList<Bed>();

				//for each bed region
				for (Bed b: chromBed.get(chr)) {
					numStartingRegions++;
					int[][] chunks = Num.chunkRegionExact(chunkSize, b.getStart(), b.getStop());
					numFinalRegions+= chunks.length;
					for (int i=0; i< chunks.length; i++){
						if (b.getName()==null) splitRegions.add(new Bed(b.getChromosome(), chunks[i][0], chunks[i][1]));
						else splitRegions.add(new Bed(b.getChromosome(), chunks[i][0], chunks[i][1], b.getName(), b.getScore(), b.getStrand()));
					}
				}

				if (bpPadding !=0) padIt(splitRegions);

				//print it
				for (Bed b: splitRegions) {
					out.println(b.toString());
				}
			}
			IO.pl("\n"+file.getName()+"\t"+numStartingRegions+" -> "+numFinalRegions);
			out.close();
		} catch (Exception e) {
			if (out!=null) out.closeNoException();
			if (outputBed!= null) outputBed.delete();
			e.printStackTrace();
			Misc.printErrAndExit("ERROR processing "+file+", ABORTING!");
		}

	}

	private void padIt(ArrayList<Bed> splitRegions) {
		//find max bp
		int maxBp = 0;
		for (Bed b: splitRegions) if (b.getStop()> maxBp) maxBp = b.getStop();
		
//IO.pl("\nMax bp "+maxBp);

		//make array
		int[] bps = new int[maxBp+bpPadding+10];
		for (int i=0; i< bps.length; i++) bps[i] = -1;
		for (int i=0; i< splitRegions.size(); i++) {
			Bed b = splitRegions.get(i);
			int start = b.getStart();
			int stop = b.getStop();
			//load the array
			for (int x = start; x< stop; x++) bps[x] = i;
		}
		
		//for (int i=0; i< bps.length; i++) if (bps[i]!=12) IO.pl("Array "+i+" : "+bps[i]);

		//for every bp to expand, advance the bps one base
		for (int x=0; x< bpPadding; x++) {
			for (int i=0; i< bps.length; i++) {
				//is it not -1, thus a bp from a split region
				if (bps[i] != -1) {
					//IO.pl("Not -1 "+bps[i]);
					//is there a -1 to the left? expand
					if (i !=0 && bps[i-1] == -1) {
						bps[i-1] = bps[i];
						//IO.pl("Adjusting left");
					}
					//is there a -1 to the right?
					if ((i+1) < bps.length && bps[i+1] == -1) {
						bps[i+1] = bps[i];
						//IO.pl("Adjusting right");
						i++;
					}
				}
			}
		}
		
		//for (int i=0; i< bps.length; i++) if (bps[i]!=-1) IO.pl("Padded "+i+" : "+bps[i]);
		
		boolean inRegion = false;
		int start = 0;
		int bedIndex = 0;
		for (int i=0; i< bps.length; i++) {
			
			//are we in a region
			if (inRegion) {
				//is this new base a -1
				if (bps[i]==-1) {
					//save old
					Bed oldBed = splitRegions.get(bedIndex);
					oldBed.setStart(start);
					oldBed.setStop(i);
//IO.pl("New -1 "+oldBed.toString());
					//reset
					inRegion = false;
				}
				//has the bedIndex changed? e.g. they are adjacent
				else if (bps[i] != bedIndex) {
					//save old
					Bed oldBed = splitRegions.get(bedIndex);
					oldBed.setStart(start);
					oldBed.setStop(i);
//IO.pl("New bp "+oldBed.toString());
					//reset
					bedIndex = bps[i];
					start = i;
				}
			}
			//OK we are not in a region
			else {
				//is this new base a region base?
				if (bps[i]!=-1) {
					//start a new region
					inRegion = true;
					bedIndex = bps[i];
					start = i;
				}
			}
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BedRegionSplitter(args);
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
					case 'd': toParse = IO.extractFiles(new File(args[++i])); break;
					case 'c': chunkSize = Integer.parseInt(args[++i]); break;
					case 'p': bpPadding = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (toParse == null || toParse.length ==0) Misc.printExit("\nError: please provide a directory containing bed region files to split.\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Bed Region Splitter : Nov 2022                       **\n" +
				"**************************************************************************************\n" +
				"Regions exceeding the chunk size are split into multiple parts, each the minimum or\n"+
				"larger. The split regions can also be padded up to a provided length without \n"+
				"overlapping adjacent regions.\n"+

				"\nRequired:\n"+
				"-d Path to a file or directory containing just bed files to split and or pad.\n"+

				"\nOptional:\n" +
				"-c Minimum BP chunk size, defaults to 2000.\n" +
				"-p Pad the split regions up to the pad size without intersecting adjacent split\n"+
				"     regions. Note the input bed files must not have overlapping regions for this to\n"+
				"     work. If in doubt run the MergeRegions app first.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/BedRegionSplitter -d ToSplit/ -c 150 -p 200 \n\n"+

				"**************************************************************************************\n");

	}


}
