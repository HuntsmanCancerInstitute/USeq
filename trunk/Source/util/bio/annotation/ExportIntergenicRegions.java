package util.bio.annotation;
import java.io.*;
import java.util.regex.*;

import util.bio.parsers.gff.Gff3Feature;
import util.bio.parsers.gff.Gff3Parser;
import util.gen.*;
import java.util.*;

public class ExportIntergenicRegions {

	private File[] gffFiles;
	private int minSize = 60;
	private boolean subtractOne = false;
	private int bpToTrim = 0;	

	public ExportIntergenicRegions (String[] args){
		processArgs(args);
		System.out.println("\nLaunching...\n\tMinimum Size "+minSize+"\n\tBP to trim "+
				bpToTrim+"\n\tSubtract one? "+subtractOne+"\n");
		exportIntergenicRegions();
	}
	
	public void exportIntergenicRegions(){
		//for each file
		for (int i=0; i< gffFiles.length; i++){
			//parse file
			Gff3Parser parser = new Gff3Parser();
			parser.setRegExTypes(".+");
			parser.setRelax(true);
			if (parser.parseIt(gffFiles[i]) == false) Misc.printExit("Problem parsing gff! "+gffFiles[i]);

			//subtract 1 from start and stop
			if (subtractOne) parser.subtractOneFromFeatures();

			//fetch by chromosome
			Gff3Feature[][] gffLines = parser.getChromSplitFeatures();

			//for each chromosome
			for (int x=0; x< gffLines.length; x++){
				//mask boolean array, those that are false are not part of annotation
				boolean[] masked = maskBoolean(gffLines[x]);
				//fetch blocks of false
				int[][] startStop = fetchFalseBlocks(masked, bpToTrim, minSize);
				//print blocks
				printBlocks(startStop, gffLines[x][0].getSeqId(), gffFiles[i]);
			}

		}

	}

	/**Prints a file containing the block regions. Strips off the extension from the file and appends the chromosome text.*/
	public boolean printBlocks(int[][] startStop, String chromosome, File file){
		try{
			String trunk = Misc.removeExtension(file.getCanonicalPath());
			PrintWriter out = new PrintWriter( new FileWriter (new File (trunk+"_"+chromosome+".bed")));
			for (int y=0; y< startStop.length; y++){
				out.println(chromosome+"\t"+startStop[y][0]+ "\t" + startStop[y][1]);
			}
			out.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}

	}

	public static int[][] fetchFalseBlocks(boolean[] b, int bpToTrim, int minSize){
		ArrayList<int[]> al = new ArrayList<int[]>();
		int startFalse = 0;
		boolean inTrue;
		//check first base and set params
		if (b[0] == true) inTrue = true;
		else {
			inTrue = false;
			startFalse = 0;
		}

		//find blocks
		for (int i=1; i< b.length; i++){
			//true found?
			if (b[i]) {
				if (inTrue == false){					
					int size = i - startFalse;
					int left = startFalse + bpToTrim;
					int right = (i-1) - bpToTrim;
					size = right - left + 1;
					if (right >= left && size >= minSize){
						//make block, new true found
						int[] block = {left, right};
						al.add(block);
					}
					inTrue = true;
				}
				//otherwise continue
			}
			//false found
			//if inTrue is true then set the startFalse
			else if (inTrue) {
				startFalse = i;
				inTrue = false;
			}
		}
		//last one?
		if (inTrue == false){
			int size = b.length - startFalse;
			int left = startFalse + bpToTrim;
			int right = (b.length-1) - bpToTrim;
			size = right - left + 1;
			if (right >= left && size >= minSize){
				//make block, new true found
				int[] block = {left, right};
				al.add(block);
			}
		}
		//convert to int[][]
		int[][] blocks = new int[al.size()][2];
		for (int i=0; i< blocks.length; i++){
			blocks[i] = (int[]) al.get(i);
		}
		return blocks;
	}

	/**Uses a boolean[] to score whether a base is part of an annotation (true) or not (false).*/
	public boolean[] maskBoolean(Gff3Feature[] features){
		//find largest base
		int largest = 0;
		for (int i=0; i< features.length; i++){
			int end = features[i].getEnd();
			if (end > largest) largest = end;
		}
		//make boolean to hold
		boolean[] bps = new boolean[largest+1];
		//for each region set booleans to true
		for (int i=0; i< features.length; i++){
			int end = features[i].getEnd() + 1;
			int start = features[i].getStart();
			for (int j= start; j< end; j++){
				bps[j] = true;
			}
		}
		return bps;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ExportIntergenicRegions(args);
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
					case 'g': gffFiles = IO.extractFiles(new File (args[i+1])); i++; break;
					case 'm': minSize = Integer.parseInt(args[i+1]); i++; break;
					case 't': bpToTrim = Integer.parseInt(args[i+1]); i++; break;
					case 's': subtractOne = true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (gffFiles == null || gffFiles.length ==0) Misc.printExit("\nError: cannot find your gff file(s)!\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Export Intergenic Regions    May 2007                     **\n" +
				"**************************************************************************************\n" +
				"EIR takes a gff file and uses it to mask a boolean array.  Parts of the boolean array\n" +
				"that are not masked are returned and represent integenic sequences. Be sure to put in\n" +
				"a gff line at the stop of each chromosome noting the last base so you caputure the last\n" +
				"intergenic region. (eg chr1 GeneDB lastBase 3600000 3600001 . + . lastBase). Base\n" +
				"coordinates are assumed to be stop inclusive, not interbase.\n\n"+
				
				"Parameters:\n"+
				"-g Full path file text for a gff file or directory containing such.\n"+
				"-t Base pairs to trim from the ends of each intergenic region, defaults to 0.\n"+
				"-m Minimum acceptable intergenic size, those smaller will be tossed, defaults to 60bp\n"+
				"-s Subtract one from the start and stop coordinates.\n\n"+
				
				"Example: java -Xmx1000M -jar pathTo/T2/Apps/ExportIntergenicRegions -s -m 100 -g\n"+
				"                 /user/Jib/GffFiles/Pombe/sanger.gff\n\n"+
		        "**************************************************************************************\n");		
	}
}
