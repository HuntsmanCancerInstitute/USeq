package trans.misc;
import java.util.*;
import java.io.*;
import java.util.regex.*;
import util.gen.*;
import trans.tpmap.*;

/**Reads in an Agilent tiling data file and outputs a tpmap and four cela files.
 * Cela files only contain intensities for chromosome mapped oligos, control intensities are set to zero.
 * Agilent uses zero based inter digit numbering, just like UCSC. Not necessarily! I've found stop base included.
 * Assumes each oligo is a 60mer. Some arrays are isothermal so the actual oligo length varies.*/
public class ConvertAgilentData {

	//fields
	private File[] agilentFiles;
	private String fileName;
	private int fileIndex;
	private int numberRows;
	private int numberColumns;
	private int numberSkippedControls;
	private int numToSubtractFromOligoStart = 0;
	private float[][] rawGreen;
	private float[][] rawRed;
	private float[][] processedGreen;
	private float[][] processedRed;
	private TPMapLine[] tpmapLines;
	private String[] controlTypes;
	private boolean switchRatio = false;
	private boolean redPresent = true;
	private boolean excludeControlProbes = false;
	private boolean sequencePresent = false;
	private Pattern dataLineMatchFilter = Pattern.compile(".*");
	private Pattern chrOligoStartRegEx = null;

	//constructor
	public ConvertAgilentData(String[] args){
		//process arguments
		processArgs(args);

		//for each file parse into a cela file and a tpmap
		for (fileIndex =0; fileIndex< agilentFiles.length; fileIndex++){
			System.out.println("\nParsing "+agilentFiles[fileIndex]+"...");

			//parse it
			if (parseAgilentFile() == false){
				Misc.printExit("\nError parsing "+fileName+"\n");
			}
			System.out.println("\n\t"+numberRows+"\tRows\t"+numberColumns+"\tColumns");
			System.out.println("\tParsed "+tpmapLines.length+" chromosome data lines.");
			if (tpmapLines.length == 0) Misc.printExit("\nError! No data lines parsed?\n");
			System.out.println("\tSkipped "+numberSkippedControls+" control data lines.");

			//save 1lq? must do before sorting tpmap array!
			if (sequencePresent){
				System.out.println("\tWriting mock 1lq file");
				save1lq();
			}

			//save tpmap and sgr files?
			if (chrOligoStartRegEx != null){
				System.out.println("\tWriting tpmap file...");
				Arrays.sort(tpmapLines);
				saveTPMap();
				//make sgr files
				if (redPresent){
					if (switchRatio) System.out.print("\tWriting log2(R/G)");
					else System.out.print("\tWriting log2(G/R)");
					System.out.println(" sgr files, score assigned to center of 60mer ...");
					writeSgrFiles();
				}
			}

			//write out cela files
			System.out.println("\tWriting cela files...");
			writeCelaFiles();
		}

		System.out.println("Done!");
	}

	public boolean parseAgilentFile (){
		fileName = Misc.removeExtension(agilentFiles[fileIndex].getName());
		String line = "";
		ArrayList tpmapLinesAL = new ArrayList();
		ArrayList controlTypesAL = new ArrayList();

		try {
			BufferedReader in = IO.fetchBufferedReader(agilentFiles[fileIndex]);
			//parse number of rows and columns
			//skip first line
			in.readLine();
			//parse 2nd to find index of Grid_NumRows and Grid_NumCols
			line = in.readLine();
			String[] tokens = line.split("\\t");
			int indexRows = -1;
			for (int i=0; i< tokens.length; i++) {
				if (tokens[i].equals("Grid_NumRows")) {
					indexRows = i;
					break;
				}
			}
			int indexColumns = -1;
			for (int i=0; i< tokens.length; i++) {
				if (tokens[i].equals("Grid_NumCols")) {
					indexColumns = i;
					break;
				}
			}
			//find values
			line = in.readLine();
			if (line.startsWith("DATA") == false || indexRows == -1 || indexColumns == -1) {
				System.err.println("\nError: Problem trying to parse rows and column numbers from "+agilentFiles[fileIndex]+"\n\tLine "+line);
				return false;
			}
			tokens = line.split("\\t");
			numberRows = Integer.parseInt(tokens[indexRows]);
			numberColumns = Integer.parseInt(tokens[indexColumns]);

			//make array to hold values
			rawGreen = new float[numberRows][numberColumns];
			rawRed = new float[numberRows][numberColumns];
			processedGreen = new float[numberRows][numberColumns];
			processedRed = new float[numberRows][numberColumns];

			//advance header to FEATURES line and find indexes, these change with the type of array
			int rowIndex =-1;
			int colIndex =-1;
			int controlIndex =-1;
			int systematicNameIndex =-1;
			int gProcessedSignalIndex =-1;
			int rProcessedSignalIndex =-1;
			int gMedianSignalIndex =-1;
			int rMedianSignalIndex =-1;
			int sequenceIndex =-1;

			boolean featuresFound = false;
			while ((line = in.readLine()) != null){
				if (line.startsWith("FEATURES")){
					HashMap hash = new HashMap();
					String[] features = line.split("\\t");
					for (int i=0; i< features.length; i++){
						hash.put(features[i], new Integer(i));
					}
					//assign indexes
					boolean missedIndex = false;
					if (hash.containsKey("Row")) rowIndex = ((Integer)hash.get("Row")).intValue();
					else missedIndex = true;
					if (hash.containsKey("Col")) colIndex = ((Integer)hash.get("Col")).intValue();
					else missedIndex = true;
					if (hash.containsKey("ControlType")) controlIndex = ((Integer)hash.get("ControlType")).intValue();
					else missedIndex = true;
					if (hash.containsKey("SystematicName")) systematicNameIndex = ((Integer)hash.get("SystematicName")).intValue();
					else missedIndex = true;
					if (hash.containsKey("gProcessedSignal")) gProcessedSignalIndex = ((Integer)hash.get("gProcessedSignal")).intValue();
					else missedIndex = true;
					if (hash.containsKey("rProcessedSignal")) rProcessedSignalIndex = ((Integer)hash.get("rProcessedSignal")).intValue();
					else missedIndex = true;
					if (hash.containsKey("gMedianSignal")) gMedianSignalIndex = ((Integer)hash.get("gMedianSignal")).intValue();
					else missedIndex = true;
					if (hash.containsKey("rMedianSignal")) rMedianSignalIndex = ((Integer)hash.get("rMedianSignal")).intValue();
					else missedIndex = true;
					if (missedIndex){
						if (rProcessedSignalIndex == -1 && rMedianSignalIndex == -1) {
							System.out.println("\nNo red channel...");
							redPresent = false;
							processedRed = null;
							rawRed = null;
						}
						else {
							System.err.println("\nError: Failed to assign all column indexes from FEATURES line ?! "+agilentFiles[fileIndex]);
							return false;
						}
					}
					//attempt to find sequence index, not in all array types
					if (hash.containsKey("Sequence")) {
						sequenceIndex = ((Integer)hash.get("Sequence")).intValue();
						sequencePresent = true;
					}
					
					featuresFound = true;
					break;
				}
			}
			if (featuresFound == false){
				System.err.println("\nError: FEATURES line not found in file?! "+agilentFiles[fileIndex]);
				return false;
			}

			//run through remaining
			numberSkippedControls = 0;
			while ((line = in.readLine()) != null){
				if (line.startsWith("DATA") == false) continue;
				String[] items = line.split("\\t");

				//skip based on ControlType?
				if (excludeControlProbes && items[controlIndex].equals("0") == false){
					numberSkippedControls++;
					continue;
				}

				//skip base on regular expression
				Matcher cont = dataLineMatchFilter.matcher(items[systematicNameIndex]);
				if (cont.matches() == false){
					numberSkippedControls++;
					continue;
				}

				//parse row and column, subtracting one to put in zero coordiantes
				int row = Integer.parseInt(items[rowIndex]) -1;
				int column = Integer.parseInt(items[colIndex]) -1;

				//parse sequence?
				String sequence = "xxxxxxxxxxxxxxxxxxxxxxxxxxxx";
				if (sequenceIndex != -1) sequence = items[sequenceIndex];

				//parse tpmap line?
				TPMapLine tpmapLine;
				if (chrOligoStartRegEx != null){
					String chrNum = items[systematicNameIndex];
					Matcher mat = chrOligoStartRegEx.matcher(chrNum);
					if (mat.matches() == false) Misc.printExit("\nERROR: cannot parse your chromosome and start postion from -> "+chrNum+"\n\tLine -> "+line+"\n\tPattern "+chrOligoStartRegEx.pattern());
					String chromosome = "chr"+mat.group(1);
					//put bp and indexes in zero coordinates
					int start = Integer.parseInt(mat.group(2)) - numToSubtractFromOligoStart;
					int end = Integer.parseInt(mat.group(3));
					//check length and adjust so oligo is centered in a mock 60mer
					double length = end - start;
					if (length < 59){
						int halfDiff = (int)Math.round((60.0 - length)/2.0);
						start = start - halfDiff;
					}
					//String chr, int pos, int row, int column, String sequence
					tpmapLine = new TPMapLine(chromosome, start, row, column, sequence);

					if (tpmapLine.getStart() == -1){
						System.err.println("\nError: Problem parsing tpmap line from Agilent data file "+agilentFiles[fileIndex]+"\n\tLine "+line);
						System.out.println("\tSN\t"+chrNum);
						System.out.println("\tchr\t"+chromosome);
						System.out.println("\tstart\t"+start);
						System.out.println("\trow\t"+row);
						System.out.println("\tcolumn\t"+column);
						System.out.println("\tseq\t"+sequence);
						return false;
					}
				}
				//make mock
				else tpmapLine = new TPMapLine("chr", 0, row, column, sequence);
				
				//set tpmapLine
				tpmapLinesAL.add(tpmapLine);
				
				//save controlType for 1lq
				controlTypesAL.add(items[controlIndex]);

				//parse for green
				float rawMedianG = Float.parseFloat(items[gMedianSignalIndex]);
				float processedG = Float.parseFloat(items[gProcessedSignalIndex]);
				rawGreen[tpmapLine.getPmX()][tpmapLine.getPmY()] = rawMedianG;
				processedGreen[tpmapLine.getPmX()][tpmapLine.getPmY()] = processedG;

				//parse for red?
				if (redPresent){
					float rawMedianR = Float.parseFloat(items[rMedianSignalIndex]);
					float processedR = Float.parseFloat(items[rProcessedSignalIndex]);
					//load float[][]
					rawRed[tpmapLine.getPmX()][tpmapLine.getPmY()] = rawMedianR;
					processedRed[tpmapLine.getPmX()][tpmapLine.getPmY()] = processedR;
				}
			}


			in.close();
			//convert AL to TPMap[]
			tpmapLines = new TPMapLine[tpmapLinesAL.size()];
			tpmapLinesAL.toArray(tpmapLines);

			//convert controlTypesAL to String[]
			controlTypes = new String[controlTypesAL.size()];
			controlTypesAL.toArray(controlTypes);

			return true;
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("\nError: Problem trying to parse Agilent data file "+agilentFiles[fileIndex]+"\n\tLine "+line);
			return false;
		}
	}

	public void writeCelaFiles(){
		File rawGFile = new File (agilentFiles[fileIndex].getParentFile(), fileName+"_RawG.cela");
		File procGFile = new File (agilentFiles[fileIndex].getParentFile(), fileName+"_ProcG.cela");
		IO.saveObject(rawGFile, rawGreen);
		IO.saveObject(procGFile, processedGreen);
		if (redPresent){
			File procRFile = new File (agilentFiles[fileIndex].getParentFile(), fileName+"_ProcR.cela");
			File rawRFile = new File (agilentFiles[fileIndex].getParentFile(), fileName+"_RawR.cela");
			IO.saveObject(rawRFile, rawRed);
			IO.saveObject(procRFile, processedRed);
		}
	}

	public boolean writeSgrFiles(){
		try {
			File rawSgrFile = new File (agilentFiles[fileIndex].getParentFile(), fileName+"_Raw.sgr");
			File procSgrFile = new File (agilentFiles[fileIndex].getParentFile(), fileName+"_Proc.sgr");
			PrintWriter outRaw = new PrintWriter(new FileWriter(rawSgrFile));
			PrintWriter outProc = new PrintWriter(new FileWriter(procSgrFile));
			for (int i=0; i<tpmapLines.length; i++){
				int pos = tpmapLines[i].getStart()+30;
				double rawRatio = rawGreen[tpmapLines[i].getPmX()][tpmapLines[i].getPmY()]/ rawRed[tpmapLines[i].getPmX()][tpmapLines[i].getPmY()];
				double procRatio = processedGreen[tpmapLines[i].getPmX()][tpmapLines[i].getPmY()]/ processedRed[tpmapLines[i].getPmX()][tpmapLines[i].getPmY()];
				rawRatio = Math.log(rawRatio)/ Num.log2;
				procRatio = Math.log(procRatio)/ Num.log2;
				if (switchRatio){
					rawRatio = -1 * rawRatio;
					procRatio = -1 * procRatio;
				}
				outRaw.println(tpmapLines[i].getChromosome()+"\t"+pos+"\t"+ rawRatio);
				outProc.println(tpmapLines[i].getChromosome()+"\t"+pos+"\t"+ procRatio);
			}
			outRaw.close();
			outProc.close();
			IO.zipAndDelete(rawSgrFile);
			IO.zipAndDelete(procSgrFile);
			return true;
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("\nError: Problem writing sgr files for "+agilentFiles[fileIndex]);
			return false;
		}
	}

	public boolean saveTPMap(){
		try {
			File tpmapFile = new File (agilentFiles[fileIndex].getParentFile(), fileName+".tpmap");
			PrintWriter outRes = new PrintWriter(new FileWriter(tpmapFile));
			for (int i=0; i<tpmapLines.length; i++){
				outRes.println(tpmapLines[i].getLine());
			}
			outRes.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("\nError: Problem trying to save a tpmap for "+agilentFiles[fileIndex]);
			return false;
		}
	}

	public boolean save1lq(){
		try {
			File file = new File (agilentFiles[fileIndex].getParentFile(), fileName+".1lq");
			PrintWriter outRes = new PrintWriter(new FileWriter(file));
			//write header
			outRes.println("X\tY\tSeq\tDestype");
			for (int i=0; i<tpmapLines.length; i++){
				//X	Y	Seq	Destype
				StringBuffer seq = new StringBuffer(tpmapLines[i].getSequence());
				outRes.println(tpmapLines[i].getPmX()+ "\t" +tpmapLines[i].getPmY()+ "\t" +seq.reverse()+ "\t" +controlTypes[i]);
			}
			outRes.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("\nError: Problem trying to save a 1lq for "+agilentFiles[fileIndex]);
			return false;
		}
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
					case 'f': agilentFiles = IO.extractFiles(new File(args[i+1]));  i++; break;
					case 's': switchRatio = true; break;
					case 'e': excludeControlProbes = true; break;
					case 'n': numToSubtractFromOligoStart = Integer.parseInt(args[i+1]); i++; break;
					case 'c': chrOligoStartRegEx = Pattern.compile(args[i+1]); i++; break;
					case 'd': dataLineMatchFilter = Pattern.compile(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					e.printStackTrace();
				}
			}
		}
		//check agilentFiles
		if (agilentFiles == null || agilentFiles.length ==0) Misc.printExit("\nCannot find your Agilent data file(s)?\n");

	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Convert Agilent Data: Nov 2007                          **\n" +
				"**************************************************************************************\n" +
				"Parses Agilent two color data files to map and cela files for use in TiMAT2.\n" +
				"Assumes each oligo is a 60mer. Also creates log2(Green/Red) sgr files for IGB.\n" +
				"WARNING, Agilent file formats are so relaxed as to make this parser near useless.\n" +
				"Each raw txt data file must be tested extensively! \n\n" +

				"-f Full path file text or directory containing text Agilent two color tiling data.\n" +
				"-e Exclude any 'ControlType' data lines that are not equal to zero.\n"+
				"-d To filter data lines based on their 'SystematicName,' provide a Java RE that will\n" +
				"      match. Non matches will be dropped. For example, '.*Spombe.+' .\n"+
				"-c To generate a tpmap and sgr files, provide a Java RE to extract the chromosome #\n" +
				"      and the oligo start and stop positions, ie for 'Spombe|complement_chr2:10-70' one would\n" +
				"      provide '.+chr(.+):(\\d+)-(\\d+)' from the 'SystematicName' column.\n"+
				"-n Number to subtract from oligo start position, defaults to zero.  Use to change\n"+
				"      coordinates to interbase/ zero based if needed.\n"+
				"-s Switch ratio to Red(cy5)/Green(cy3) for making sgr files.\n\n" +

				"Example: java -Xmx1500M -jar pathTo/T2/Apps/ConvertAgilentData -f /badData/ -c\n" +
				"      '.+chr(.+):(\\d+)-(\\d+)' -s -d '.*Spombe.+' -e\n\n" +

		"**************************************************************************************\n");		
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new ConvertAgilentData(args);
	}

}
