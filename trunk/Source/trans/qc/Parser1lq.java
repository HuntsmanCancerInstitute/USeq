package trans.qc;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;

/**
 * Extracts control an pm lines from a text 1lq file.
 * A serialized int[][][] for the xy coordinates is saved. 
 */
public class Parser1lq {
	private File[] onelqFiles = null;
	private int numRowsOrColumns = 2560;
	private boolean saveVirtualCelFiles = false;
	
	//constructors
	public Parser1lq(String[] args){
		processArgs(args);
		System.out.println("\nProcessing...");
		extract1lqFiles();
		System.out.println("\nDone!");
	}
	
	//main methods
	/**Fetches coordinates for DESTYPEs of 0, 3, -3, 6, or -6 (Control) and PM -111, 111.
	 * Lumps sense and antisense together -6/6, -3/-3, -111/111
	 * int[probe index][x,y] order noSynth, dim, bright, and pm.
	 * Returns null if a problem is encountered.*/
	public static int[][][] extract1lqFile(File onelqFile, File outputFile){	
		try{
			ArrayList zero = new ArrayList(10000);
			ArrayList three = new ArrayList(10000);
			ArrayList negThree = new ArrayList(10000);
			ArrayList six = new ArrayList(10000);
			ArrayList negSix = new ArrayList(10000);
			ArrayList pmAl = new ArrayList(100000);
			BufferedReader in = new BufferedReader( new FileReader(onelqFile));
			boolean inHeader = true;
			String line;
			//skip header
			while (inHeader){
				line = in.readLine();
				if (line == null) return null;
				if (line.startsWith("X")) inHeader=false;
			}
			//check DESTYPE 
			Pattern controlPattern = Pattern.compile("-?[3|6|0]");
			Pattern pmPattern = Pattern.compile("-?111"); 
			Matcher m;
			String[] tokens;
			while ((line = in.readLine()) !=null){
				tokens = line.split("\\s+");
				//is it a PM oligo
				m = pmPattern.matcher(tokens[3]);
				if (m.matches()){
					int[] xy = {Integer.parseInt(tokens[0]), Integer.parseInt(tokens[1])};
					pmAl.add(xy);
				}
				else {
					m = controlPattern.matcher(tokens[3]);
					if (m.matches()) {
						//which DESTYPE
						int[] xy = {Integer.parseInt(tokens[0]), Integer.parseInt(tokens[1])};
						int number = Integer.parseInt(tokens[3]);
						if (number == 0) {
							zero.add(xy);
						}
						else if (number == 3) three.add(xy);
						else if (number == -3) negThree.add(xy);
						else if (number == 6) six.add(xy);
						else if (number == -6) negSix.add(xy);
						else {
							Misc.printExit("\nWarning no match of destype "+number+" orig token "+tokens[3]+" line "+line);
							
						}
					}
				}
			}
			in.close();
			//convert ArrayLists to int[][]
			int[][] zeroA = new int[zero.size()][2];
			int[][] threeA = new int[three.size()][2];
			int[][] negThreeA = new int[negThree.size()][2];
			int[][] sixA = new int[six.size()][2];
			int[][] negSixA = new int[negSix.size()][2];
			int[][] pm = new int[pmAl.size()][2];
			
			zero.toArray(zeroA);
			three.toArray(threeA);
			negThree.toArray(negThreeA);
			six.toArray(sixA);
			negSix.toArray(negSixA);
			pmAl.toArray(pm);
			
			int[][][] ndb = new int[6][][];
			ndb[0] = zeroA;
			ndb[1] = threeA;
			ndb[2] = negThreeA;
			ndb[3] = sixA;
			ndb[4] = negSixA;
			ndb[5] = pm;
			return ndb;
		} catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}
	
	/**Extracts control oligo information from 1lq files.*/
	public void extract1lqFiles(){
		try {
			for (int i=0; i< onelqFiles.length; i++){
				System.out.println("\t"+onelqFiles[i].getName());
				File output = new File (onelqFiles[i].getCanonicalPath()+".controls");
				int[][][] ndb = extract1lqFile(onelqFiles[i], output);
				int[] numCoor = countNumberIntensities(ndb);
				System.out.println("\t\t# 0\t"+numCoor[0]);
				System.out.println("\t\t# 3\t"+numCoor[1]);
				System.out.println("\t\t# -3\t"+numCoor[2]);
				System.out.println("\t\t# 6\t"+numCoor[3]);
				System.out.println("\t\t# -6\t"+numCoor[4]);
				System.out.println("\t\t# PM (-111,111)\t"+numCoor[5]);
				//save an int[][][] of xy coordinates for noSynth, dim, bright
				File arrays = new File (onelqFiles[i].getCanonicalPath()+".SerCont");
				IO.saveObject(arrays, ndb);
				
				//save virtual celas for each class
				if (saveVirtualCelFiles) {
					float[][] vc = makeVirtualCel(ndb[0],numRowsOrColumns);
					File vcFile = new File (onelqFiles[i].getCanonicalPath()+"NoSynth.cela");
					IO.saveObject(vcFile, vc);
					vc = makeVirtualCel(ndb[1],numRowsOrColumns);
					vcFile = new File (onelqFiles[i].getCanonicalPath()+"Dim.cela");
					IO.saveObject(vcFile, vc);
					vc = makeVirtualCel(ndb[2],numRowsOrColumns);
					vcFile = new File (onelqFiles[i].getCanonicalPath()+"Bright.cela");
					IO.saveObject(vcFile, vc);
					vc = makeVirtualCel(ndb[3],numRowsOrColumns);
					vcFile = new File (onelqFiles[i].getCanonicalPath()+"PM.cela");
					IO.saveObject(vcFile, vc);
				}
			}
		}catch(IOException e){
			e.printStackTrace();
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
		return new int[]{coordinates[0].length, coordinates[1].length, coordinates[2].length, coordinates[3].length, coordinates[4].length, coordinates[5].length};
	}
	
	//main
	public static void main(String[] args) {
		if (args.length!=0) {
			new Parser1lq(args);
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
					case 'n': numRowsOrColumns = Integer.parseInt(args[i+1]); i++; break;
					case 's': saveVirtualCelFiles = true; break;
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
		
	}
	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Parser 1lq: Jan  2006                           **\n" +
				"**************************************************************************************\n" +
				"Extracts control and pm coordinates from a text 1lq file.  A serialized int[][][]\n" +
				"is saved to disk for each 1lq file. Also creates virtual cel files for each class for\n" +
				"visualization using the VirtualCel app.\n\n"+
				
				"-f Full path to a 'xxx.1lq' file or directory containing such for extraction.\n" +
				//"-s Save virtual 'xxx.cela' files for each set of coordinates, default is no.\n"+
				"-n Number of rows/ columns on the array, defaults to 2560.\n\n"+
				
				"Example: java -Xmx512M CoordinateExtractor1lq -f /data/hmn1lqFiles/ -n 1280\n\n"+
				
				"You may need to increase the heap size of the java virtual machine using the -Xmx\n" +
				"flag. Adjust accordingly in response to out of memory errors and the available\n" +
				"resources. \n\n" +
				
				"Questions, comments, suggestions? Contact Gingeras Group or David_Nix@Affymetrix.com\n" +
		"**************************************************************************************\n");
	}
	
}
