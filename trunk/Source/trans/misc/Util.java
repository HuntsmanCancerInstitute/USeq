package trans.misc;
import java.io.*;

import trans.main.*;
import trans.roc.*;
import util.gen.*;

import java.util.*;

/**
 * Static helper methods for all things Affy
 */
public class Util {
	
	
	public static void main(String[] args)  {
		//write a bar file
		int[] positions = {1,5,10,15,20,25,30};
		float[] values = {10,50,100,150,200,250,300};
		String chromosome = "chrMT";
		String groupVersion = "hg17";
		String strand = "+";
		File fil = new File("/Users/nix/Desktop/test.bar");
		writeSimpleBarFile(chromosome, groupVersion, strand, positions, values, fil);
		
		readSimpleGrBarFile(fil);
	}
	
	/**Sorts an array of Interval by the median ratio of the best SubWindow.  If a sub is not
	 * found then the score is set to best window score index 1.
	 * Modifies the original array.
	 * Sets the sortBy field to the median ratio score.*/
	public static void sortIntervalsBySubWindowMedianRatio(Interval[] intervals){
		int numIntervals = intervals.length;
		for (int i=0; i< numIntervals; i++){
			double score = 0;
			if (intervals[i].getBestSubWindow()!=null){
				score = intervals[i].getBestSubWindow().getMedianRatio();
			}
			else score = intervals[i].getBestWindow().getScores()[1];
			intervals[i].setSortBy(score);
		}
		Arrays.sort(intervals);
	}
	
	/**Splits an ArrayList containing Interval by chromosome into a HashMap.*/
	public static HashMap splitIntervalArrayListByChromosome (ArrayList intervals){
		HashMap map = new HashMap();
		int num = intervals.size();
		for (int i=0; i< num; i++){
			Interval x = (Interval) intervals.get(i);
			if (map.containsKey(x.getChromosome())){
				ArrayList al = (ArrayList) map.get(x.getChromosome());
				al.add(x);
			}
			else {
				ArrayList al = new ArrayList();
				al.add(x);
				map.put(x.getChromosome(), al);
			}
		}
		return map;
	}

	/**Finds min max score*/
	public static double[] minMaxWindowScores(Window[] win, int scoreIndex){
		double min=0;
		double max=0;
		for (int i=0; i< win.length; i++){
			double score = win[i].getScores()[scoreIndex];
			if (score < min) min = score;
			else if (score > max) max = score;
		}
		return new double[] {min, max};
	}

	/**This is a complete hack of Gregg Helt's IGB code to parse bar files.  It only works on single graph xxx.gr type files, 
	 * ie two columns, ints and floats.
	 * Is compatable with TiMAT's Gr2Bar writer.*/
	public static GrGraph readSimpleGrBarFile(File simpleBarFile){
		try {
			FileInputStream fis = new FileInputStream(simpleBarFile);
			DataInputStream dis = new DataInputStream(new BufferedInputStream(fis));
			//read header
			byte[] headbytes = new byte[8];
			dis.readFully(headbytes);
			float barVersion = dis.readFloat();       
			dis.readInt();
			int vals_per_point = dis.readInt();
			
			int[] val_types = new int[vals_per_point];
			for (int i=0; i<vals_per_point; i++) {
				val_types[i] = dis.readInt();
			}
			int tvcount = dis.readInt(); 
			HashMap tagValues = readTagValPairs(dis, tvcount);
			//read seq
			int namelength = dis.readInt();
			byte[] barray = new byte[namelength];
			dis.readFully(barray);
			String chromosome = new String(barray);
			boolean bar2 = (barVersion >= 2);
			if (bar2) {
				int grouplength = dis.readInt();
				barray = new byte[grouplength];
				dis.readFully(barray);
			}
			int verslength = dis.readInt();			
			barray = new byte[verslength];
			dis.readFully(barray);
			String seqversion = new String(barray);
			// hack to extract seq version and seq text from seqname field for bar files that were made
			//   with the version and text concatenated (with ";" separator) into the seqname field
			int sc_pos = chromosome.lastIndexOf(";");
			if (sc_pos >= 0) {
				seqversion = chromosome.substring(0, sc_pos);
				chromosome = chromosome.substring(sc_pos+1);
			}
			
			if (bar2) {
				int seq_tagval_count = dis.readInt();
				tagValues.putAll(readTagValPairs(dis, seq_tagval_count));
			}
			
			int total_points = dis.readInt();
			//build GrGraph
			int[] basePositions = new int[total_points];
			float[] values = new float[total_points];
			
			for (int i= 0; i<total_points; i++) {
				basePositions[i] = dis.readInt();
				values[i] = dis.readFloat();
				//System.out.println(chromosome+"\t"+basePositions[i]+"\t"+values[i]);
			}
			dis.close();
			//System.out.println(tagValues);
			return new GrGraph(chromosome, seqversion, basePositions, values);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	/**Methodfor parsing tag values from a binary stream.*/
	public static HashMap<String,String> readTagValPairs(DataInput dis, int pair_count) throws IOException  {
		HashMap<String,String> tagValues = new HashMap();
		for (int i=0; i<pair_count; i++) {
			int taglength = dis.readInt();
			byte[] barray = new byte[taglength];
			dis.readFully(barray);
			String tag = new String(barray);
			int vallength = dis.readInt();
			barray = new byte[vallength];
			dis.readFully(barray);
			String value = new String(barray);
			tagValues.put(tag, value);
		}
		return tagValues;
	}
	
	/**
	 * Writes a simple bar format graph for a xxx.gr type graph.  Modified from IGB code.
	 */
	public static boolean writeSimpleBarFile(GrGraph gr, File barFile){
		return writeSimpleBarFile(gr.getChromosome(), gr.getGenomeVersion(), ".", gr.getBasePositions(), gr.getValues(), barFile);
	}
	
	/**
	 * Writes a simple bar format graph for a xxx.gr type graph.  Modified from IGB code.
	 * @param chromosome - ie 'chr21'
	 * @param versionedGenome - ie 'hg18' or 'H_sapiens_Mar_2006' where possible follow UCSC
	 * @param strand - '+, -, or .' no exceptions. This is written out as a sequence tag value pair.
	 */
	public static boolean writeSimpleBarFile(String chromosome, String versionedGenome, String strand, int[] basePositions, float[] intensityValues, File barFile) {
		
		int BYTE4_FLOAT = 1;
		int BYTE4_SIGNED_INT = 2;
		boolean success = false;
		DataOutputStream dos;
		
		//groupId and version set to same - ie hg17
		String groupid = versionedGenome;	//chipset?
		String version = groupid;		//genome version?
		//chromosome - ie chr2
		String seqid = chromosome;
		try {
			dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(barFile)));
			dos.writeBytes("barr\r\n\032\n");  // char  "barr\r\n\032\n"
			dos.writeFloat(2.0f);       // version of bar format = 2.0
			dos.writeInt(1);  // number of seq data sections in file -- if single graph, then 1
			dos.writeInt(2);  // number of columns (dimensions) per data point
			dos.writeInt(BYTE4_SIGNED_INT);  // int  first column/dimension type ==> 4-byte signed int
			dos.writeInt(BYTE4_FLOAT);  // int  second column/dimension type ==> 4-byte float
			
			// should write out all properties from group and/or graphs as tag/vals?  For now just saying no tag/vals
			dos.writeInt(0);
			
			// assuming one graph for now, so only one seq section
			dos.writeInt(seqid.length());
			dos.writeBytes(seqid);
			dos.writeInt(groupid.length());
			dos.writeBytes(groupid);
			dos.writeInt(version.length());
			dos.writeBytes(version);
			
			// write out strand as a tag value, the tag is "strand" the value "+, -, or ."
			dos.writeInt(1);	//Integer	4	Number of tag/value pairs.
			dos.writeInt("strand".length());	//Integer	4	The number of characters in the text of the tag. Referred to as TAGNAMELEN.
			dos.writeBytes("strand");	//Char	TAGNAMELEN	The text of the tag.  Used only in version 2.0 or greater.
			dos.writeInt(strand.length());	//Integer	4	The number of characters in the value part of the tag/value pair. Referred to as TAGVALLEN.
			dos.writeBytes(strand);	//Char	TAGVALLEN	The value of the tag/value pair.
			
			
			int total_points = basePositions.length;
			
			dos.writeInt(total_points);
			for (int i=0; i<total_points; i++) {
				dos.writeInt(basePositions[i]);
				dos.writeFloat(intensityValues[i]);
				//System.out.println(chromosome+"\t"+basePositions[i]+"\t"+values[i]);				
			}
			dos.close(); 
			success = true;
		}
		catch (Exception ex) {
			ex.printStackTrace();
			success = false;
		}
		return success;
	}
	
	/**Attempts to parse a chromosome text. Returns null if it fails.
	 * Looks for chr1-22, chrX, chrY, chrM, and chrMT.*/
	public static String parseChromosomeName(String name){
		String[] chroms = getNumberChromosome();
		for (int x=chroms.length -1; x>=1; x--){
			if (name.indexOf(chroms[x]) !=-1) return chroms[x];
		}
		return null;
	}
	
	/**Creates a HashMap of "chr1"=Byte(1) through "chr22"=Byte(22) plus 
	 * chrX=23, chrY=24, and chrMT=25. */
	public static HashMap getChromosomeNumber(){
		HashMap chromNumber = new HashMap(26);
		for (byte i=1; i<23; i++){
			chromNumber.put("chr"+i, new Byte(i));
		}
		chromNumber.put("chrX",new Byte((byte)23));
		chromNumber.put("chrY",new Byte((byte)24));
		chromNumber.put("chrM",new Byte((byte)25));
		chromNumber.put("chrMT",new Byte((byte)26));
		return chromNumber;
	}
	/**Creates a String[] containing [1]"chr1", [2]"chr2" ... [23]"chrX", 
	 * [24]"chrY", [25]"chrM", [26]"chrM", [27]"chrMT", and [28]"chrCtrls".
	 * [0] is empty. */
	public static String[] getNumberChromosome(){
		String[] numChrom = new String[28];
		for (byte i=1; i<23; i++){
			numChrom[i] = "chr"+i;
		}
		numChrom[23] ="chrX";
		numChrom[24] ="chrY";
		numChrom[25] ="chrM";
		numChrom[26] ="chrMT";
		numChrom[27] ="chrCtrls";
		return numChrom;
	}
	
	/**Returns the min and max values of a given score index in a Window[].
	 * @return double[]{min, max}*/
	public static double[] minMaxWindowArray(Window[] win, int scoreIndex){
		double min = 10000;
		double max = -10000;
		for (int i=0; i< win.length; i++){
			double score = win[i].getScores()[scoreIndex];
			if (score < min) min = score;
			if (score > max) max = score;
		}
		return new double[]{min, max};
	}
	
	/**Removes windows that contain < minNumberOligos.*/
	public static Window[] removeLowOligoWindows(Window[] win, int minNumberOligos){
		ArrayList al = new ArrayList(win.length);
		for (int i=0; i< win.length; i++){
			if (win[i].getNumberOligos() >= minNumberOligos) al.add(win[i]);
		}
		Window[] win2 = new Window[al.size()];
		al.toArray(win2);
		al = null;
		return win2;
	}
	
	/**Reads in an sgr file.*/
	public static Sgr[] loadSgrFile(File sgrFile){
		ArrayList sgrs = new ArrayList(10000);
		try{
			String line;
			BufferedReader in = new BufferedReader(new FileReader(sgrFile));
			while ((line=in.readLine())!=null){
				if (line.trim().length()==0) continue;
				//make Sgr to hold line info
				sgrs.add(new Sgr(line));
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		Sgr[] x = new Sgr[sgrs.size()];
		sgrs.toArray(x);
		return x;
	}
	
	/**Reads in an sgr file.*/
	public static Sgr[] loadGrFile(File grFile, String chromName){
		ArrayList sgrs = new ArrayList(100000);
		try{
			String line;
			BufferedReader in = new BufferedReader(new FileReader(grFile));
			while ((line=in.readLine())!=null){
				if (line.trim().length()==0) continue;
				//make Sgr to hold line info
				sgrs.add(new Sgr(line, chromName));
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		Sgr[] x = new Sgr[sgrs.size()];
		sgrs.toArray(x);
		return x;
	}
	
	
	/**Fetches intensity values from a text version cel file.*/
	public static float[] fetchTextCelFileIntensities (File celFile){
		ArrayList al = new ArrayList(6553600); //2560x2560
		try{
			BufferedReader in = new BufferedReader(new FileReader(celFile));
			//skip header
			boolean inHeader = true;
			String line;
			while (inHeader) {
				line = in.readLine();
				if (line.indexOf("CellHeader") != -1) inHeader = false;
			}
			//parse xy coords and intensity and add to float array
			String[] tokens;
			while ((line = in.readLine()) !=null){
				tokens = line.trim().split("\\s+");
				//check that a data line was found and not stop of file or junk
				if (tokens.length!=5){
					break;
				}
				//add to intensity array
				Float f = new Float(tokens[2]);
				al.add(f);
			}
			in.close();
		} catch (IOException e){
			e.printStackTrace();
		}
		return Num.arrayListOfFloatToArray(al);
	}
	
	
	/**Loads a text version cel file into a zeroed virtual chip.*/
	public static float[][] createVirtualCel (File celFile){
		int numRows = 0;
		float[][] intensities;
		//load virtual slide with .cel file intensities
		try{
			BufferedReader in = new BufferedReader(new FileReader(celFile));
			//skip header
			boolean inHeader = true;
			String line;
			int counter = 0;
			while (inHeader) {
				line = in.readLine();
				//set number of rows
				int index = line.indexOf("Rows=");
				if (index !=-1 ) {
					numRows = Integer.parseInt(line.substring(index+"Rows=".length()));
				}
				if (line.indexOf("CellHeader") != -1) inHeader = false;
			}
			
			//rows set?
			if (numRows == 0) {
				Misc.printExit("\nProblem parsing number of rows from cel file header!\n");
			}
			
			//create a float[][] loaded with 0's
			intensities = Num.zeroedFloatArray(numRows, numRows);

			//parse xy coords and intensity and add to float array
			String[] tokens;
			int x;
			int y;
			float inten;
			while ((line = in.readLine()) !=null){
				counter++;
				tokens = line.trim().split("\\s+");
				//check that a data line was found and not stop of file or junk
				if (tokens.length!=5){
					//System.out.println ("\tProcessed "+counter+" cel file data lines.");
					break;
				}
				//add to intensity array
				x = Integer.parseInt(tokens[0]);
				y = Integer.parseInt(tokens[1]);
				inten = Float.parseFloat(tokens[2]);
				intensities[x][y]=inten;
			}
			in.close();
		} catch (IOException e){
			e.printStackTrace();
			return null;
		}
		return intensities;
	}
	
	
	/**Loads a text version cel file into a zeroed virtual chip.*/
	public static float[][] createVirtualCel (String[] lines){
		int numRows = 0;
		
		//load virtual slide with .cel file intensities
		int numberLines = lines.length;
		int index = 0;
		//skip header
		boolean inHeader = true;
		String line;
		int counter = 0;
		while (inHeader) {
			line = lines[index++];
			//set number of rows
			int indexSub = line.indexOf("Rows=");
			if (indexSub !=-1 ) {
				numRows = Integer.parseInt(line.substring(indexSub+"Rows=".length()));
			}
			if (line.indexOf("CellHeader") != -1) inHeader = false;
		}
		
		//rows set?
		if (numRows == 0) {
			Misc.printExit("\nProblem parsing number of rows from cel file header!\n");
		}
		
		//create a float[][] loaded with 0's
		float[][] intensities = Num.zeroedFloatArray(numRows, numRows);
		
		//parse xy coords and intensity and add to float array
		String[] tokens;
		int x;
		int y;
		float inten;
		//while ((line = in.readLine()) !=null){
		for (; index < numberLines; index++){
			counter++;
			tokens = lines[index].trim().split("\\s+");
			//check that a data line was found and not stop of file or junk
			if (tokens.length!=5){
				//System.out.println ("\tProcessed "+counter+" cel file data lines.");
				break;
			}
			//add to intensity array
			x = Integer.parseInt(tokens[0]);
			y = Integer.parseInt(tokens[1]);
			inten = Float.parseFloat(tokens[2]);
			intensities[x][y]=inten;
		}
		
		return intensities;
	}
	
	/**Given an array of int[][start, stop] returns the total number of integers, stop included.*/
	public static int countNumberOfOligoPositions(int[][] startStops){
		int numberOligos = 0;
		for (int x=0; x< startStops.length; x++){
			int windowSize = 1+startStops[x][1]-startStops[x][0];
			numberOligos += windowSize;
		}
		return numberOligos;
	}
	
	public static Window[][] splitWindowsByChromosome(Window[] win){
			String currentChrom = win[0].getChromosome();
			ArrayList windowArrays = new ArrayList();
			ArrayList chromWinAL = new ArrayList();
			for (int j=0; j< win.length; j++){
				if (win[j].getChromosome().equals(currentChrom)) chromWinAL.add(win[j]);
				else {
					//convert to Window[] and save in ArrayList
					Window[] chromWin = new Window[chromWinAL.size()];
					chromWinAL.toArray(chromWin);
					windowArrays.add(chromWin);

					//begin anew
					currentChrom = win[j].getChromosome();
					chromWinAL.clear();
					chromWinAL.add(win[j]);
				}
			}
			//save last
			Window[] chromWin = new Window[chromWinAL.size()];
			chromWinAL.toArray(chromWin);
			windowArrays.add(chromWin);
			
			//convert to Array
			int num = windowArrays.size();
			Window[][] split = new Window[num][];
			for (int i=0; i< num; i++) split[i] = (Window[]) windowArrays.get(i);
			return split;
	}
	
}
