package trans.misc;
import java.util.*;
import java.io.*;
import java.util.regex.*;
import edu.utah.seq.parsers.*;
import util.gen.*;
import trans.tpmap.*;

/**Reads in multiple directories of bar files and makes xxx.cela and xxx.tpmap files. */
public class ConvertBars2TPMapCel {

	//fields
	private File[] barDirectories;
	private File resultsDirectory;
	private boolean delog = false;
	private float defaultValue = 1;
	
	//internal fields
	private HashMap<String,HashSet<Integer>> chromPos = new HashMap<String,HashSet<Integer>>();
	String[] chromosomes;
	int[][] positions;
	int[] baseIndexes;
	int totalPositions = 0;
	
	//constructor
	public ConvertBars2TPMapCel( String[] args){

		//process arguments
		processArgs(args);

		//scan all bar files to get master list of chromosomes and positions
		scanBarDirectories();
		
		//build the tpmap file from the master
		buildTPMAP();
		
		//build cela files (float[1][] arrays)
		buildCelas();
		makeBlankCela();

		System.out.println("\nDone!");
	}
	
	public void makeBlankCela(){
		float[][] cela = new float[1][totalPositions];
		Arrays.fill(cela[0], defaultValue);
		IO.saveObject(new File(resultsDirectory, "blank_1.cela"), cela);
	}
	
	public void buildCelas(){
		//for each barDirectory
		for (int i=0; i<barDirectories.length; i++){
			System.out.println("\tMaking cela file for "+barDirectories[i].getName());
			
			//make cela
			float[][] cela = new float[1][totalPositions];
			Arrays.fill(cela[0], defaultValue);
			//for each chromosome attempt to load cela
			BarParser bp = new BarParser();
			for (int j=0; j<chromosomes.length; j++){
				File barFile = new File(barDirectories[i], chromosomes[j]+".bar");
				if (barFile.exists() == false) barFile = new File(barDirectories[i], chromosomes[j]+".bar.zip");
				if (barFile.exists() == false) {
					System.out.println("\tNo chrom found for "+chromosomes[j]);
					continue;
				}
				//fetch positions and values
				bp.readBarFile(barFile, true);
				int[] pos = bp.getBasePositions();
				float[] val = bp.getValues();
				if (delog) val = Num.antiLog(val, 2);
				//check vals
				for (int x=0; x<val.length; x++) if (Float.isNaN(val[x])) val[x] = defaultValue;
				//for each base position find index in master and add value
				for (int k=0; k< pos.length; k++){
					int index = Arrays.binarySearch(positions[j], pos[k]);
					int extendedIndex = index + baseIndexes[j];
					cela[0][extendedIndex] = val[k];
					//System.out.println("\t\t"+index+"\t"+extendedIndex+"\t"+cela[0][extendedIndex]);
				}
				//save it				
				IO.saveObject(new File(resultsDirectory, barDirectories[i].getName()+".cela"), cela);
			}
		}
	}
	
	public void buildTPMAP(){
		//make chromosomes and positions
		ArrayList<String> chromsAL = new ArrayList<String>();
		chromsAL.addAll(chromPos.keySet());
		chromosomes = Misc.stringArrayListToStringArray(chromsAL);
		Arrays.sort(chromosomes);
		positions = new int[chromosomes.length][];
		baseIndexes = new int[chromosomes.length];
		//for each chrom
		try{
			PrintWriter out = new PrintWriter ( new FileWriter ( new File(resultsDirectory,"bar.tpmap")));
			out.println("#seq_group_name\n#version");
			for (int i=0; i< chromosomes.length; i++){
				//set base index
				baseIndexes[i] = totalPositions;
				//get positions
				HashSet<Integer> posHS = chromPos.get(chromosomes[i]);
				int[] pos = Num.hashSetToInt(posHS);
				Arrays.sort(pos);
				positions[i] = pos;
				//write mock tpmap line
				String pre = "x\tf\t"+chromosomes[i]+"\t";
				for (int j=0; j< pos.length; j++){
					out.println(pre + positions[i][j]+"\t0\t"+totalPositions);
					totalPositions++;
				}
			}
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}
			
		
	}
	
	public void scanBarDirectories(){
		//for each barDirectory
		for (int i=0; i<barDirectories.length; i++){
			System.out.println(barDirectories[i].getName());
			//get bar files
			File[] bars = IO.extractFiles(barDirectories[i], ".bar");
			if (bars == null || bars.length == 0) bars = IO.extractFiles(barDirectories[i], ".bar.zip");
			if (bars == null || bars.length == 0) {
				System.out.println("\tSkipping, no xxx.bar or xxx.bar.zip files.");
				continue;
			}
			//for each bar file load and add positions to master chromPos
			BarParser bp = new BarParser();
			for (int j=0; j< bars.length; j++){
				bp.readBarFile(bars[j], true);
				String chrom = bp.getChromosome();
				//get hashSet or make it
				HashSet<Integer> hash;
				if (chromPos.containsKey(chrom)) hash = chromPos.get(chrom);
				else {
					hash = new HashSet<Integer>();
					chromPos.put(chrom, hash);
				}
				//add positions
				int[] positions = bp.getBasePositions();
				for (int k=0; k< positions.length; k++) hash.add(new Integer(positions[k]));
			}
		}
	}


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		//load
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': barDirectories = IO.extractOnlyDirectories(new File(args[i+1]));  i++; break;
					case 'd': delog = true; break;
					case 'r': resultsDirectory = new File(args[i+1]);  i++; break;
					case 's': defaultValue = Float.parseFloat(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check
		if (resultsDirectory == null) Misc.printExit("\nError making/ finding results directory, aborting.\n");
		if (barDirectories == null || barDirectories.length ==0) Misc.printExit("\nCannot find your bar directories, aborting.\n");
		resultsDirectory.mkdir();

	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Convert Bars2TPMapCel: April 2007                        **\n" +
				"**************************************************************************************\n" +
				"Reads in multiple directories of bar files and makes xxx.cela and xxx.tpmap files.\n" +
				"\n"+

				"Options:\n"+
				"-b Full path directory containing directories of xxx.bar files.\n" +
				"-r Full path directory text for saving the results.\n" +
				"-d Delog2 values.\n" +
				"-s Score to set for missing values, defaults to 1.\n\n" +

				"Example: java -Xmx1500M -jar pathTo/T2/Apps/ConvertBars2TPMapCell -f /data/ -d -r\n" +
				"              /data/converted\n\n" +

		"**************************************************************************************\n");		
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new ConvertBars2TPMapCel(args);
	}

}
