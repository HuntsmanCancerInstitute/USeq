package trans.misc;
import java.util.*;
import java.io.*;
import java.util.regex.*;

import util.gen.*;
import trans.tpmap.*;

/**Reads and parses in a GEO data file into a tpmap and xxx.cela files for TiMAT2 processing. ie
#ID_REF = 
#Rval = Cy5 (red) signal, Input total DNA from WI38 cells
#Gval = Cy3 (green) signal, MeDIP DNA from WI38 cells
#VALUE = log2(Gval/Rval)
ID_REF	Rval	Gval	VALUE
CHR1P000009757_HSAP0406S00000003	2934.44	6303.33	1.1030
CHR1P000009637_HSAP0406S00000003	4742.78	13346.00	1.4926
CHR1P000009627_HSAP0406S00000003	3311.56	8825.89	1.4142
CHR1P000009617_HSAP0406S00000003	3295.78	8769.11	1.4118
CHR1P000009607_HSAP0406S00000003	3812.56	10092.00	1.4044
 *
 *Saves xxx.cela files, and a xxx.tpmap
 */
public class ConvertGeoData {

	//fields
	private File[] dataFiles;
	private int indexChrom = 0;
	private int indexRed = 0;
	private int indexGreen = 0;
	private int lengthOligo = 60;
	private String mockOligo;
	private boolean switchRatioToRedGreen = false;
	private int minNumberColumns = 3;
	
	//internal fields
	private String trimmedFileName;
	private TPMapLine[] tpmapLines;
	private float[][] red;
	private float[][] green;
	private HashSet chromNames = new HashSet();
	
	
	//constructor
	public ConvertGeoData( String[] args){
		//process arguments
		processArgs(args);
		
		//make mock oligo
		StringBuffer sb = new StringBuffer();
		for (int i=0; i< lengthOligo; i++) sb.append("x");
		mockOligo = sb.toString();

		//for each file parse 
		for (int i =0; i< dataFiles.length; i++){
			System.out.println("\nParsing "+dataFiles[i]+"...");
			//parse file saving sgr
			if (parseIt(dataFiles[i]) ==false) Misc.printExit("\nError: problem parsing file, aborting.\n");
			//save cel files
			File cela = new File (dataFiles[i].getParentFile(), trimmedFileName+"_R.cela");
			IO.saveObject(cela, red);
			red = null;
			cela = new File (dataFiles[i].getParentFile(), trimmedFileName+"_G.cela");
			IO.saveObject(cela, green);
			green = null;
			//save tpmap file
			Arrays.sort(tpmapLines);
			File tpmapFile = new File (dataFiles[i].getParentFile(), trimmedFileName + ".tpmap");
			if (TPMapLine.saveTPMap(tpmapLines, tpmapFile) == false) Misc.printExit("\nError: problem writing tpmap file, aborting.\n");

		}
		System.out.println("\nChromosomes found: "+chromNames);

		System.out.println("\nDone!");
	}
	


	public boolean parseIt(File file){
		try {
			trimmedFileName = Misc.removeExtension(file.getName());
			//count data lines
			int numLines = countDataLines(file);
			if (numLines == 0){
				System.out.println("No data lines found? Skipping.");
				return false;
			}
			System.out.println(numLines+" # data lines found.");
			//instantiate objects
			red = new float[1][numLines];
			green = new float[1][numLines];
			tpmapLines = new TPMapLine[numLines];
			File sgrFile = new File (file.getParentFile(), trimmedFileName+".sgr");
			PrintWriter outSgr = new PrintWriter (new FileWriter (sgrFile));
			
			//skip header
			BufferedReader in = new BufferedReader( new FileReader(file));
			String line;
			while ((line = in.readLine()) != null){
				if (line.startsWith("#")) continue;
				else break;
			}
			
			//parse data
			Pattern tab = Pattern.compile("\\t+");
			Pattern chromPosPat = Pattern.compile("CHR(.+)P(\\d+)_.+");
			int counter = 0;
			while ((line = in.readLine()) != null){
				
				//parse columns
				String[] cols = tab.split(line, 0);
				if (cols.length < minNumberColumns) {
					System.out.println("\nError: data line found with < "+minNumberColumns+" columns! "+line);
					return false;
				}
				red[0][counter] = Float.parseFloat(cols[indexRed]);
				green[0][counter] = Float.parseFloat(cols[indexGreen]);
				
				//calc ratio
				float ratio;
				if (switchRatioToRedGreen) ratio = Num.log2(red[0][counter]/green[0][counter]);
				else ratio = Num.log2(green[0][counter]/red[0][counter]);
				
				//parse chromosome and base position
				Matcher mat = chromPosPat.matcher(cols[indexChrom]);
				if (mat.matches()== false){
					System.out.println("\nError: cannot parse chromosome and position information from "+line);
					return false;
				}
				//record chromosome
				chromNames.add(mat.group(1));
				String chromosome = "chr"+mat.group(1);
				
				//switch those ending in R with _random
				if (chromosome.endsWith("R")){
					chromosome = chromosome.substring(0, chromosome.length()-1) + "_random";
				}
				int position = Integer.parseInt(mat.group(2));
				
				//write sgr line
				outSgr.println(chromosome +"\t"+position+"\t"+ratio);
				
				//instantiate TPMapLine
				tpmapLines[counter] = new TPMapLine(chromosome, position, 0, counter, mockOligo);
				
				//increment counter
				counter++;
			}
			
			//close print writers
			in.close();
			outSgr.close();
			
			//zip sgr file
			IO.zipAndDelete(sgrFile);
			
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}

	public int countDataLines(File file){
		try{
			//skip # and first non #, damn this is bad form, what's with GEO!
			BufferedReader in = new BufferedReader( new FileReader(file));
			String line;
			while ((line = in.readLine()) != null){
				if (line.startsWith("#")) continue;
				else break;
			}
			//count lines
			int numDataLines = 0;
			while (in.readLine() != null) numDataLines ++;

			in.close();
			return numDataLines;
		}
		catch (Exception e){
			e.printStackTrace();
		}
		return -1;
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
					case 'f': dataFiles = IO.extractFiles(new File(args[i+1]));  i++; break;
					case 's': switchRatioToRedGreen = true; break;
					case 'c': indexChrom = Integer.parseInt(args[i+1]); i++; break;
					case 'o': lengthOligo = Integer.parseInt(args[i+1]); i++; break;
					case 'r': indexRed = Integer.parseInt(args[i+1]); i++; break;
					case 'g': indexGreen = Integer.parseInt(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check agilentFiles
		if (dataFiles == null || dataFiles.length ==0) Misc.printExit("\nCannot find your data file(s)?\n");
		System.out.println("\nIndexes: "+indexChrom+"-ChromPos "+indexRed+"-Red "+indexGreen+"-Green");
		System.out.println("Switching ratio to red/green? "+switchRatioToRedGreen);

	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Convert Geo Data: April 2007                           **\n" +
				"**************************************************************************************\n" +
				"Parses a Gene Expression Omnibus individual data file for a particular sample, not a\n" +
				"Download family file.  The SOFT, MINiMl and TXT formats are useless since they don't\n" +
				"give access to the raw data only the processed ratios. This is very free form. Check\n" +
				"your results carefully! Assumes two color data files, makes a tpmap and cela files\n" +
				"for use in TiMAT2. Also creates log2(Green/Red) sgr files for IGB. Column indexes\n" +
				"start with zero! Your GEO data file should look something like the following:\n\n" +

				"   #ID_REF = \n" +
				"   #Rval = Cy5 (red) signal, Input total DNA from WI38 cells\n" +
				"   #Gval = Cy3 (green) signal, MeDIP DNA from WI38 cells\n" +
				"   #VALUE = log2(Gval/Rval)\n" +
				"              ID_REF                 Rval    Gval     VALUE\n" +
				"   CHR1P000009757_HSAP0406S00000003 2934.45 6303.33   1.1030\n" +
				"   CHR1P000009637_HSAP0406S00000003 4742.78 13346.00  1.4926\n" +
				"   CHR1P000009627_HSAP0406S00000003 3311.56 8825.89   1.4142\n" +
				"   ...\n\n" +

				"Options:\n"+
				"-f Full path file text or directory containing text files.\n" +
				"-o Length of oligo.\n"+
				"-c Column index of chromosome position.\n" +
				"-r Column index of green intensity values.\n" +
				"-g Column index of red intensity values.\n" +
				"-s Switch ratio to Red/Green for making sgr files.\n\n" +

				"Example: java -Xmx1500M -jar pathTo/T2/Apps/ConvertGeoData -f /data/ -c 0 -r 2 -g 1\n\n" +

		"**************************************************************************************\n");		
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new ConvertGeoData(args);
	}

}
