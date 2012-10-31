package trans.cel;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import trans.misc.*;
import java.util.*;

import trans.tpmap.*;

public class MakeChromosomeSets {
	
	//fields
	private File[] setDirectories = null;	//directories of directories containing split chromosome IntensityFeature[]s, one for each set of chips
	private boolean makeOligoPositions = true;
	private File oligoPositionsDir = null;
	private String[] chromosomes;// = Util.getNumberChromosome();
	private boolean verbose = true;
	private File saveDirectory;
	
	//constructor
	public MakeChromosomeSets (String[] args){
		System.out.println("\nLaunching...");
		processArgs(args);	
		//oligo positions?
		if (makeOligoPositions) {
			oligoPositionsDir = new File (saveDirectory, "OligoPositions");
			if (oligoPositionsDir.exists()) IO.deleteDirectory(oligoPositionsDir);
			oligoPositionsDir.mkdir();
		}
		System.out.println("\tProcessing "+saveDirectory.getName());
		
		//find all chromosomes
		HashSet chroms = new HashSet();
		for (int x=0; x< setDirectories.length; x++){
			File[] files = IO.extractFiles(setDirectories[x]);
			for (int y=0; y< files.length; y++) chroms.add(files[y].getName());
		}
		chromosomes = Misc.hashSetToStringArray(chroms);
		System.out.println("\tChromosomes found: "+Misc.stringArrayToString(chromosomes, ", "));
		
		//for each chromosome attempt to collect IntensityFeature[]s from each directory, combine, sort, and save
		long totalOligos = makeSets(setDirectories, saveDirectory);
		System.out.println("\tTotal\t"+totalOligos + " oligo values\n");
		
		//save total for comparison 
		if (totalOligos !=0 ){
			File total = new File (saveDirectory, "totalNumberOligos.txt");
			IO.writeString(totalOligos+"", total);
		}
	}
	
	/**Merges chromosome info together and saves as whole chromosomes, will also save the oligoPositions if oligoPositionsDirectory isn't null.
	 * Returns the total number of oligos found across all chromosomes.*/
	public long makeSets(File[] chipSetDirectories, File saveDirectory){

		long totalOligos = 0;
		//for each chromosome attempt to collect IntensityFeature[]s, sort save
		for (int j=0; j< chromosomes.length; j++){
			File[] chromInts = fetchChromosomeFiles(chipSetDirectories, chromosomes[j]);			
			if (chromInts.length != 0 ) {
				IntensityFeature[] ints = fetchIntensityFeatures(chromInts);
				if (verbose) System.out.println("\t\t"+chromosomes[j]+"\t"+ints.length + " oligo values");
				totalOligos += ints.length;
				Arrays.sort(ints);
				File chromFile = new File (saveDirectory, chromosomes[j]);
				splitAndSaveIntensities(ints, chromFile);
				if (oligoPositionsDir !=null){
					File pos = new File(oligoPositionsDir, chromosomes[j]);
					splitAndSavePositions(ints, pos);
				}
			}
		}
		return totalOligos;
	}
	//methods
	/**Saves the intensities to file as a int[]*/
	public static void splitAndSaveIntensities(IntensityFeature[] ints, File f){
		float[] values = new float[ints.length];
		for (int i=0; i< ints.length; i++){
			values[i] = ints[i].intensity;
		}
		IO.saveObject(f, values);
	}
	/**Saves the positions to file as a int[]*/
	public static void splitAndSavePositions(IntensityFeature[] ints, File f){
		int[] pos = new int[ints.length];
		for (int i=0; i< ints.length; i++){
			pos[i] = ints[i].start;
		}
		IO.saveObject(f, pos);
	}
	
	/**Loads IntensityFeature[]s into a single master.*/
	public static IntensityFeature[] fetchIntensityFeatures(File[] f){
		IntensityFeature[][] ints = new IntensityFeature[f.length][];
		for (int i=0; i<f.length; i++){
			ints[i] = (IntensityFeature[])IO.fetchObject(f[i]);
		}
		IntensityFeature[] ifs = collapseIntensityFeatureArray(ints);
		return ifs;
	}
	
	/**Collapses an uneven array of IntensityFeature[][] to IntensityFeature[].*/
	public static IntensityFeature[] collapseIntensityFeatureArray (IntensityFeature[][] fs){
		int numI = fs.length;
		int total = 0;
		for (int i=0; i< numI; i++){
			for (int j=0; j< fs[i].length; j++){
				total++;
			}
		}
		IntensityFeature[] combine = new IntensityFeature[total];
		int k =0;
		for (int i=0; i< numI; i++){
			for (int j=0; j< fs[i].length; j++){
				combine[k++] = fs[i][j];
			}
		}
		return combine;
	}
	
	/**Extracts all the Files with a given chromosome text from a list of directories*/
	public static File[] fetchChromosomeFiles(File[] chipDir, String chromosomeName){
		ArrayList files = new ArrayList();
		for (int i=0; i< chipDir.length; i++){
			File test = new File(chipDir[i], chromosomeName);
			if (test.exists()) files.add(test);
		}
		File[] x = new File[files.size()];
		files.toArray(x);
		return x;
	}
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		try{
			String dir = null;
			Pattern pat = Pattern.compile("-[a-z]");
			for (int i = 0; i<args.length; i++){
				String lcArg = args[i].toLowerCase();
				Matcher mat = pat.matcher(lcArg);
				if (mat.matches()){
					char test = args[i].charAt(1);
					switch (test){
					case 'd': dir = args[i+1]; i++; break;
					case 'n': saveDirectory = new File (args[i+1]); i++; break;
					case 's': makeOligoPositions = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
					
				}
			}
			if (dir == null) Misc.printExit("\nPlease enter a comma delimited list, no spaces, of directories containing split chromosome intensity files.\n");
			setDirectories = IO.extractFiles(dir);
			if (setDirectories == null) Misc.printExit("\nProblem parsing your directories -> "+dir+"\n");
			if (saveDirectory == null) Misc.printExit("\nPlease provide a full path directory text to use in saving your joined chromosomes.\n");
			if (saveDirectory.exists() == false) saveDirectory.mkdir();
		}
		catch (Exception e){
			System.out.println("\nSorry, something doesn't look right with your command line options!\n");
			e.printStackTrace();
			System.exit(1);
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Make Chromosome Sets: Feb  2007                         **\n" +
				"**************************************************************************************\n" +
				"MCS takes directories containing split chromosome intensity files from the\n" +
				"CelProcessor app and combines them into a master set of chromosome split intensity\n" +
				"files for use by ScanChromosomesCNV.\n" +
				"\n" +
				"Use the following options when running CP:\n\n" +
				"-d Comma delimited list of directories containing split chromosome files from\n" +
				"      CelProcessor, no spaces. These should represent one complete replica.\n" +
				"-n Full path directory text to save the combine chomosomes.\n"+
				"-s Skip making chromosome oligo positions (you already have them for this chip set).\n" +
				"\n" +
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/MakeChromosomeSets -d\n" +
				"     /ProcCelFiles/Chip1RepA,/ProcCelFiles/Chip2RepA,/ProcCelFiles/Chip3RepA -n\n" +
				"     /ProcCelFiles/RepA\n" +
				
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MakeChromosomeSets(args);
	}
	
}
