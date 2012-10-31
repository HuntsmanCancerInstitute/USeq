package trans.qc;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;

/**
 * Extracts control an pm lines from a text 1lq file.
 * A serialized int[][][] for the xy coordinates is saved. 
 */
public class CoordinateExtractor1lq {
	private File[] onelqFiles = null;
	private int numRowsOrColumns;
	private int numRowsMinusOne;
	private boolean saveVirtualCelFiles = false;
	private File resultsDirectory = null;
	private boolean make7G = true;
	
	//constructors
	public CoordinateExtractor1lq(String[] args){
		processArgs(args);
		if (make7G) System.out.println("\nExtracting coordinates for 7G cel files.");
		else System.out.println("\nExtracting coordinates for pre 7G, unrotated, cel files.");
		extract1lqFiles();
		System.out.println("\nDone!");
	}
	
	//main methods
	/**Fetches coordinates for DESTYPEs of 0, 3, -3, 6, or -6 (Control) and PM -111, 111.
	 * Lumps sense and antisense together -6/6, -3/-3, -111/111
	 * int[probe index][x,y] order noSynth, dim, bright, and pm.
	 * Returns null if a problem is encountered.*/
	public int[][][] extract1lqFile(File onelqFile, File outputFile){	
		try{
			ArrayList dimAl = new ArrayList(10000);
			ArrayList brightAl = new ArrayList(10000);
			ArrayList noSynthAl = new ArrayList(10000);
			ArrayList pmAl = new ArrayList(100000);
			BufferedReader in = new BufferedReader( new FileReader(onelqFile));
			PrintWriter out = new PrintWriter( new FileWriter(outputFile));

			String line;
			//parse header and pull number of rows
			while (true){
				line = in.readLine();
				if (line == null) return null;
				if (line.startsWith("COLS")){
					String[] tokens = line.split("\\s+");
					numRowsOrColumns = Integer.parseInt(tokens[1]);
					numRowsMinusOne = numRowsOrColumns -1;
					System.out.println("\tNumber rows/ columns = "+numRowsOrColumns);
				}
				else if (line.startsWith("X")) break;
			}
			//did number of rows get parsed?
			if (numRowsOrColumns == 0) Misc.printExit("\nFailed to parse the number of rows/ columns. There should " +
					"be a line in the header like 'COLS/ROWS=2560  2560    Dummy'. Fix and rerun.");
			//check DESTYPE 
			Pattern controlPattern = Pattern.compile("-?[3|6|0]");
			Pattern pmPattern = Pattern.compile("-?111"); 
			Matcher m;
			String[] tokens;
			while ((line = in.readLine()) !=null){
				tokens = line.split("\\s+");
				int[] xy = new int[2]; 
				//rotate it?
				if (make7G){
					xy[0] = Integer.parseInt(tokens[1]);
					xy[1] = numRowsMinusOne - Integer.parseInt(tokens[0]);
				}
				else {
					xy[0] = Integer.parseInt(tokens[0]);
					xy[1] = Integer.parseInt(tokens[1]);
				}
				//is it a PM oligo
				m = pmPattern.matcher(tokens[3]);
				if (m.matches()){
					pmAl.add(xy);
				}
				else {
					m = controlPattern.matcher(tokens[3]);
					if (m.matches()) {
						out.println(line);
						//which DESTYPE			
						if (tokens[3].equals("0")) {
							noSynthAl.add(xy);
						}
						else if (tokens[3].charAt(0)== '-') {
							
							dimAl.add(xy);
						}
						else {
							brightAl.add(xy);
						}
					}
				}
			}
			in.close();
			out.close();
			//convert ArrayLists to int[][]
			int[][] dim = new int[dimAl.size()][2];
			int[][] bright = new int[brightAl.size()][2];
			int[][] noSynth = new int[noSynthAl.size()][2];
			int[][] pm = new int[pmAl.size()][2];
			
			dimAl.toArray(dim);
			brightAl.toArray(bright);
			noSynthAl.toArray(noSynth);
			pmAl.toArray(pm);
			int[][][] ndb = new int[4][][];
			ndb[0] = noSynth;
			ndb[1] = dim;
			ndb[2] = bright;
			ndb[3] = pm;
			return ndb;
		} catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}
	
	/**Extracts control oligo information from 1lq files.*/
	public void extract1lqFiles(){
			for (int i=0; i< onelqFiles.length; i++){
				System.out.println("\t"+onelqFiles[i].getName());
				File output = new File (resultsDirectory, onelqFiles[i].getName()+".controls");
				int[][][] ndb = extract1lqFile(onelqFiles[i], output);
				int[] numCoor = countNumberIntensities(ndb);
				System.out.println("\t\t# No Synth (0)\t"+numCoor[0]);
				System.out.println("\t\t# Dim (-3,-6)\t"+numCoor[1]);
				System.out.println("\t\t# Bright (3,6)\t"+numCoor[2]);
				System.out.println("\t\t# PM (-111,111)\t"+numCoor[3]);
				//save an int[][][] of xy coordinates for noSynth, dim, bright
				File arrays = new File (resultsDirectory,onelqFiles[i].getName()+".SerCont");
				IO.saveObject(arrays, ndb);
				
				//save virtual celas for each class
				if (saveVirtualCelFiles) {
					float[][] vc = makeVirtualCel(ndb[0],numRowsOrColumns);
					File vcFile = new File (resultsDirectory,onelqFiles[i].getName()+"NoSynth.cela");
					IO.saveObject(vcFile, vc);
					vc = makeVirtualCel(ndb[1],numRowsOrColumns);
					vcFile = new File (resultsDirectory,onelqFiles[i].getName()+"Dim.cela");
					IO.saveObject(vcFile, vc);
					vc = makeVirtualCel(ndb[2],numRowsOrColumns);
					vcFile = new File (resultsDirectory,onelqFiles[i].getName()+"Bright.cela");
					IO.saveObject(vcFile, vc);
					vc = makeVirtualCel(ndb[3],numRowsOrColumns);
					vcFile = new File (resultsDirectory,onelqFiles[i].getName()+"PM.cela");
					IO.saveObject(vcFile, vc);
				}
			}

	}
	
	/**Make a virtual cel file loading xyCoors with an intensity of 100.*/
	public static float[][] makeVirtualCel(int[][] xyCoor, int numRows){
		float[][] vc = Num.loadFloatArray(numRows, numRows,1); //new float[numRows][numRows];
		for (int i=0; i< xyCoor.length; i++){
			int[] xy = xyCoor[i];
			vc[xy[0]][xy[1]] = 100;
		}
		return vc;
	}
	
	/**Calculate number of noSynth, dim, bright controls, and pm probes.*/
	public static int[] countNumberIntensities(int[][][] coordinates){
		return new int[]{coordinates[0].length, coordinates[1].length, coordinates[2].length, coordinates[3].length};
	}
	
	//main
	public static void main(String[] args) {
		if (args.length!=0) {
			new CoordinateExtractor1lq(args);
		}
		else {
			printDocs();
		}
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File dirFile = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': dirFile = new File(args[i+1]); i++; break;
					case 'r': resultsDirectory = new File(args[i+1]); i++; break;
					case 's': saveVirtualCelFiles = true; break;
					case 'p': make7G = false; break;
					case 'h': printDocs(); System.exit(0);
					default: {
						System.out.println("\nError: unknown option! -> " + mat.group()+"\n");
						System.exit(0);
					}
					}
				}
				catch (Exception e){
					System.out.print("\nError: something doesn't look right with this parameter request: -"+test+"\n");
					System.exit(0);
				}
			}
		}
		//check for tempDirectory
		if (dirFile == null || dirFile.exists()== false){
			System.out.println("\nError: cannot find your 1lq file or directory?! -> "+dirFile+"\n");
			System.exit(0);
		}
		//fetch 1lq files
		onelqFiles = IO.extractFiles(dirFile, "1lq");
		if (onelqFiles.length == 0){
			System.out.println("\nError: cannot find any 'xxx.1lq' files in -> "+dirFile+"\n");
			System.exit(0);
		} 
		Arrays.sort(onelqFiles);
		//check for an alternative results directory
		if (resultsDirectory !=null){ 
			if (resultsDirectory.isDirectory()== false) Misc.printExit("\nError: cannot find or read your alternative results directory -> "+resultsDirectory+"\n");
		}
		else resultsDirectory = onelqFiles[0].getParentFile();
		
	}
	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Coordinate Extractor 1lq: Sept  2007                        **\n" +
				"**************************************************************************************\n" +
				"Extracts control and pm coordinates from an unrotated text 1lq file.  An int[][][]\n" +
				"is saved to disk for each 1lq file. Also creates virtual cel files for each class for\n" +
				"visualization using the VirtualCel app.\n\n"+
				
				"-f Full path to a 'xxx.1lq' file or directory containing such for extraction.\n" +
				"-r Alternative results directory, defaults to the xxx.1lq parent directory.\n"+
				"-s Save virtual 'xxx.cela' files for each set of coordinates, each is assigned an\n" +
				"      intensity of 100, default is no.\n"+
				"-p Extract coordinates for pre 7G, unrotated, cel files, default is to rotate 1lq.\n\n"+
				
				"Example: java -Xmx512M -jar pathTo/T2/Apps/CoordinateExtractor1lq -f\n" +
				"      /data/hmn1lqFiles/ \n\n"+
				
		"**************************************************************************************\n");
	}
	
}
