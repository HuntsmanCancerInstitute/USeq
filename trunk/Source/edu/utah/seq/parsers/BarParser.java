package edu.utah.seq.parsers;
import java.io.*;
import java.util.*;
import java.util.zip.*;

import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;

import util.gen.*;
import trans.misc.*;


/**Class for reading and writing xxx.bar files.
 * @author Nix
 *
 <pre>
   Bar format definition:

   MAIN HEADER FOR ENTIRE FILE

   1	Char	8	The file type identifier. This is always set to "barr\r\n\032\n".
   2	Float	4	The file version number.  Valid versions are 1.0 and 2.0
   3	Integer	4	The number of sequences stored in the file. Referred to as NSEQ.
   4	Integer	4	The number of columns per data point. Referred to as NCOL.
   5	Integer	4*NCOL	The field types, one per column of data. The possible values are:
             0 - Double
             1 - Float
             2 - 4 byte signed integer
             3 - 2 byte signed integer
			 4 - 1 byte signed integer
			 5 - 4 byte unsigned integer
			 6 - 2 byte unsigned integer
			 7 - 1 byte unsigned integer
   6	Integer	4	Numbern of tag/value pairs.
   7	Integer	4	The number of characters in the text of the tag. Referred to as TAGNAMELEN.
   8	Char	TAGNAMELEN	The text of the tag.
   9	Integer	4	The number of characters in the value part of the tag/value pair. Referred to as TAGVALLEN.
   10	Char	TAGVALLEN	The value of the tag/value pair.


   SECTION HEADER FOR EACH BAR SEQ/DATA 

   11	Integer	4	The number of characters in the text of the sequence. Referred to as SEQNAMELEN.
   12	Char	SEQNAMELEN	The sequence text (ie chr21).
   13	Integer	4	The number of characters in the text of the sequence group.  Referred to as SEQGROUPNAMELEN.  Used only in version 2.0 or greater.
   14	Char	SEQGROUPNAMELEN	The text of the group of which the sequence is a member (for example, often specifies organism, ie Homo sapiens).  Referred to as SEQGROUPNAME.  Used only in version 2.0 or greater.
   15	Integer	4	The number of characters in the sequence version text. Referred to as SEQVERLEN.
   16	Char	SEQVERLEN	The sequence version. (ie hg18)
   17	Integer	4	Number of tag/value pairs.  Used only in version 2.0 or greater.
   18	Integer	4	The number of characters in the text of the tag. Referred to as TAGNAMELEN.  Used only in version 2.0 or greater.
   19	Char	TAGNAMELEN	The text of the tag.  Used only in version 2.0 or greater.
   20	Integer	4	The number of characters in the value part of the tag/value pair. Referred to as TAGVALLEN.  Used only in version 2.0 or greater.
   21	Char	TAGVALLEN	The value of the tag/value pair.  Used only in version 2.0 or greater.
   22	Integer	4	The number of data points defined in the sequence. Each data point will contain NCOL column values.
   23			The next set of values in the file is the data points for the sequence. Each data point contains NCOL column values. The type, thus the size, of each column is defined above in the field types section.
 </pre> 

 *Note, several reserve xxx_TAG names are used to put info into the tagValue pair HashMap, see below, these include:
 *Strand information can be recorded '+, -, or .'
 *Read length, for signature sequencers such as SOLiD and Solexa.
 *
 *Loads xxx.bar.zip compressed files too!
 * */
public class BarParser {

	//fields
	private File barFile;
	private HashMap <String,String> tagValues = new HashMap <String,String> ();
	private String chromosome;
	private String versionedGenome;
	private String strand = ".";
	private int[] basePositions;
	private float[] values;
	private int numberPositionValues;
	private double scoreTotal = 0;
	private boolean zipCompress = false;
	private boolean loadPositionValues = false;
	private boolean mergeIdenticalPositions = false;
	private boolean closeDataInputStream = true;
	private boolean replaceWithHitCount = false;
	private DataInputStream dis;
	
	/**Reserved tag names for writing to bar file tagValue pairs.
	 * Don't use these keys for other purposes when putting extra info into the tagValue HashMap.
	 * These have specific meaning to the IGB browser*/
	//what is the strand? +, -, or .
	public static final String STRAND_TAG = "strand";
	//what is the size of the individual reads in the file?
	public static final String READ_LENGTH_TAG = "readLength";
	//where is the data from? Parsed from file x?
	public static final String SOURCE_TAG = "source";
	//what kind of graph type should be used to display
	public static final String GRAPH_TYPE_TAG = "initialGraphStyle";
	public static final String GRAPH_TYPE_BAR = "Bar";
	public static final String GRAPH_TYPE_DOT = "Dot";
	public static final String GRAPH_TYPE_LINE = "Line";
	public static final String GRAPH_TYPE_MINMAXAVE = "Min/Max/Ave";
	public static final String GRAPH_TYPE_STAIRSTEP = "Stairstep";
	public static final String GRAPH_TYPE_HEATMAP = "Heat Map";
	//what color, hex color values #0000FF
	public static final String GRAPH_TYPE_COLOR_TAG = "initialColor";

	//bound Y
	public static final String GRAPH_MIN_Y_TAG = "initialMinY";
	public static final String GRAPH_MAX_Y_TAG = "initialMaxY";
	//description?
	public static final String DESCRIPTION_TAG = "description";
	//what units are the bar float values?
	public static final String UNIT_TAG = "unit";
	//bp shift off original alignment position
	public static final String BP_3_PRIME_SHIFT = "bp3PrimeShift";
	//window size?
	public static final String WINDOW_SIZE = "windowSize";
	//score total
	public static final String SCORE_TOTAL = "scoreTotal";
	
	
	//methods
	/**Writes a xxx.bar file to disk. 
	 * @param barFile - should stop with '.bar'. If zip compressing, the barFile text will be appended with the '.zip' extension.
	 * @param strand - either +, -, or .
	 * @param tagValues - optional input, can be null.*/
	public boolean writeBarFile(File barFile, String chromosome, String versionedGenome, char strand, 
			int[] basePositions, float[] values, HashMap <String, String> tagValues){
		this.barFile = barFile;
		this.chromosome = chromosome;
		this.versionedGenome = versionedGenome;
		this.strand = new String (new char[]{strand});
		this.basePositions = basePositions;
		this.values = values;
		if (tagValues != null) {
			this.tagValues = tagValues;
		}
		try {
			writeSimpleBarFile();
			return true;
		} catch (IOException e){
			e.printStackTrace();
			return false;
		}
	}

	/**Reads a bar file and sets the appropriate fields in the BarParser.
	 * Set loadPositionValues = true to load the position values or false to
	 * just load the header info.*/
	public boolean readBarFile(File barFile, boolean loadPositionValues){
		try {
			this.barFile = barFile.getCanonicalFile();
			this.loadPositionValues = loadPositionValues;
			loadSimpleBarFile();
			return true;
		} catch (IOException e){
			System.err.println("\nProblem reading bar file -> "+barFile);
			System.err.println("\tLoad data boolean? "+loadPositionValues);
			System.err.println("\tExist? "+barFile.exists());
			System.err.println("\tReadable? "+barFile.canRead());
			e.printStackTrace();
			System.exit(1);
		}
		return false;
	}
	
	public Point fetchNextPoint() {
		try {
			int pos = dis.readInt();
			float val = dis.readFloat();
			return new Point(pos,val);
		} catch (Exception e){
			return null;
		}
	}
	

	/**This is a complete hack of Gregg Helt's IGB code to parse bar files.  
	 * It only works on single graph type files, ie two columns, int and float,
	 * one sequence. Reads xxx.bar or zip compressed xxx.bar.zip files.
	 * @param loadPositionValues - set to true if you want to load the actual 
	 * positions and values, otherwise just the header is loaded.*/
	public void loadSimpleBarFile() throws IOException {
		//reset HashMap
		tagValues.clear();
		//get streams
		if (barFile.getName().endsWith(".zip")){
			ZipFile zf = new ZipFile(barFile);
			ZipEntry ze = (ZipEntry) zf.entries().nextElement();
			dis = new DataInputStream (new BufferedInputStream (zf.getInputStream(ze)));
		}
		else if (barFile.getName().endsWith(".gz")){
			dis = new DataInputStream(new GZIPInputStream(new FileInputStream(barFile)));
		}
		else {
			FileInputStream fis = new FileInputStream(barFile);
			dis = new DataInputStream(new BufferedInputStream(fis));
		}

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
		readTagValPairs(dis, tvcount);
		//read seq
		int namelength = dis.readInt();
		byte[] barray = new byte[namelength];
		dis.readFully(barray);
		chromosome = new String(barray);
		boolean bar2 = (barVersion >= 2);
		if (bar2) {
			int grouplength = dis.readInt();
			barray = new byte[grouplength];
			dis.readFully(barray);
		}
		int verslength = dis.readInt();
		barray = new byte[verslength];
		dis.readFully(barray);
		versionedGenome = new String(barray);
		// hack to extract seq version and seq text from seqname field for bar files that were made
		//   with the version and text concatenated (with ";" separator) into the seqname field
		int sc_pos = chromosome.lastIndexOf(";");
		if (sc_pos >= 0) {
			versionedGenome = chromosome.substring(0, sc_pos);
			chromosome = chromosome.substring(sc_pos+1);
		}
		if (bar2) {
			int seq_tagval_count = dis.readInt();
			readTagValPairs(dis, seq_tagval_count);
		}
		numberPositionValues = dis.readInt();
		if (loadPositionValues) {
			if (mergeIdenticalPositions) loadAndMerge(dis);
			else loadPositionValues(dis);
		}
		if (closeDataInputStream) dis.close();
		//set strand
		if (tagValues.containsKey(STRAND_TAG)) strand = tagValues.get(STRAND_TAG);
		else strand = ".";
		//set totalScore
		if (tagValues.containsKey(SCORE_TOTAL)) scoreTotal = Double.parseDouble(tagValues.get(SCORE_TOTAL));
	}
	
	private void loadPositionValues(DataInputStream dis) throws IOException{
		//load basePositions and values
		basePositions = new int[numberPositionValues];
		values = new float[numberPositionValues];
		for (int i= 0; i<numberPositionValues; i++) {
			basePositions[i] = dis.readInt();
			values[i] = dis.readFloat();
		}
	}
	
	private void loadAndMerge(DataInputStream dis) throws IOException{
		//load basePositions and values
		ArrayList<Point> pts = new ArrayList<Point>();
		int pos = dis.readInt();
		float val = dis.readFloat();
		if (replaceWithHitCount) val = 1;
		for (int i= 1; i<numberPositionValues; i++) {
			//read next
			int testPos = dis.readInt();
			float testVal = dis.readFloat();
			//add?
			if (testPos == pos) {
				if (replaceWithHitCount) val += 1;
				else val += testVal;
			}
			else {
				pts.add(new Point(pos, val));
				pos = testPos;
				if (replaceWithHitCount) val = 1;
				else val = testVal;
			}
		}
		//add last
		pts.add(new Point(pos, val));
		Point[] p = new Point[pts.size()];
		pts.toArray(p);
		PointData pd = Point.extractPositionScores(p);
		basePositions = pd.getPositions();
		values = pd.getScores();
		numberPositionValues = values.length;
	}

	/**Method for parsing tag values from a binary stream.*/
	public void readTagValPairs(DataInput dis, int pair_count) throws IOException  {
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
	}

	/** Writes a simple bar format graph using the BarParser field values.
	 *  Based on Gregg Helt's IGB code.*/
	public void writeSimpleBarFile() throws IOException {
		//defined in spec above
		int BYTE4_FLOAT = 1;
		int BYTE4_SIGNED_INT = 2;

		//write header
		DataOutputStream dos = new DataOutputStream(new BufferedOutputStream( new FileOutputStream(barFile)));
		dos.writeBytes("barr\r\n\032\n");  // char  "barr\r\n\032\n"
		dos.writeFloat(2.0f);       // version of bar format = 2.0
		dos.writeInt(1);  // number of seq data sections in file -- if single graph, then 1
		dos.writeInt(2);  // number of columns (dimensions) per data point
		dos.writeInt(BYTE4_SIGNED_INT);  // int  first column/dimension type ==> 4-byte signed int
		dos.writeInt(BYTE4_FLOAT);  // int  second column/dimension type ==> 4-byte float

		// write out tag value pairs that apply to whole document. Not implemented
		writeTagValuePairs (dos, null);
		//for each sequence (only one at present)
		//write chromosome
		dos.writeInt(chromosome.length());
		dos.writeBytes(chromosome);
		//write species
		dos.writeInt(versionedGenome.length());
		dos.writeBytes(versionedGenome);
		//write genome version		
		dos.writeInt(versionedGenome.length());
		dos.writeBytes(versionedGenome);
		//set strand in tagValues
		tagValues.put(STRAND_TAG, strand);
		//set scoreTotal
		tagValues.put(SCORE_TOTAL, Num.sumArrayReturnDouble(values)+"");
		//write out sequence specific tag values?
		writeTagValuePairs (dos, tagValues);
		//write points
		int total_points = basePositions.length;		
		dos.writeInt(total_points);
		for (int i=0; i<total_points; i++) {
			dos.writeInt(basePositions[i]);
			dos.writeFloat(values[i]);				
		}
		
		//close stream
		dos.close();
		
		//zip compress?
		if (zipCompress){
			if (IO.zipAndDelete(barFile)==false) throw new IOException("Failed to zip compress "+barFile);
			//rename barFile to barFile.zip
			barFile = new File(barFile.getCanonicalPath()+".zip");
		}
	}
	
	public static void writeTagValuePairs(DataOutputStream dos,HashMap<String,String> tagValuePairs) throws IOException{
		//any tag value pairs?
		if (tagValuePairs == null || tagValuePairs.size() == 0) dos.writeInt(0);
		else {
			//write number of pairs
			dos.writeInt(tagValuePairs.size());
			//write tag values preceded by their respective lengths
			Iterator<String> it = tagValuePairs.keySet().iterator();
			while (it.hasNext()){
				String tag = it.next();
				dos.writeInt(tag.length());
				dos.writeBytes(tag);
				String value = tagValuePairs.get(tag);
				dos.writeInt(value.length());
				dos.writeBytes(value);
			}
		}
	}

	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("Bar File:\t"+barFile);
		sb.append("\nVersioned Genome:\t"+versionedGenome);
		sb.append("\nChromosome:\t"+chromosome);
		sb.append("\nStrand:\t"+strand);
		if (tagValues.size()>1) sb.append("\nTag Values:\t"+tagValues);
		int length = basePositions.length;
		sb.append("\n# Base Values:\t"+length);
		sb.append("\nBasePostions\tValues");
		if (length > 25) length = 25;
		for (int i=0; i< length; i++){
			sb.append("\n");
			sb.append(basePositions[i]);
			sb.append("\t");
			sb.append(values[i]);
		}
		sb.append("\n...");
		return sb.toString();
	}
	
	/*
	public static void main (String[] args){
		BarParser bp = new BarParser();
		bp.readBarFile(new File ("/Users/nix/Desktop/chr13.bar.zip"));
		bp.setStrand("-");
		bp.getTagValues().put("readLength", "35");
		bp.setZipCompress(true);
		File f = new File ("/Users/nix/Desktop/chr13Mod.bar");
		bp.setBarFile(f);
		try {bp.writeSimpleBarFile();} catch(Exception e){}
		System.out.println(bp.getBarFile());
		bp.readBarFile(bp.getBarFile());
		
		System.out.println(bp.getTagValues()+"\n"+bp);
	}*/
	

	public File getBarFile() {
		return barFile;
	}

	public void setBarFile(File barFile) {
		this.barFile = barFile;
	}

	public int[] getBasePositions() {
		return basePositions;
	}

	public void setBasePositions(int[] basePositions) {
		this.basePositions = basePositions;
	}

	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public String getStrand() {
		return strand;
	}

	public void setStrand(String strand) {
		this.strand = strand;
	}

	public HashMap<String, String> getTagValues() {
		return tagValues;
	}

	public void setTagValues(HashMap<String, String> tagValue) {
		this.tagValues = tagValue;
	}

	public float[] getValues() {
		return values;
	}

	public void setValues(float[] values) {
		this.values = values;
	}

	public String getVersionedGenome() {
		return versionedGenome;
	}

	public void setVersionedGenome(String versionedGenome) {
		this.versionedGenome = versionedGenome;
	}

	public boolean isZipCompress() {
		return zipCompress;
	}

	public void setZipCompress(boolean zipCompress) {
		this.zipCompress = zipCompress;
	}

	public int getNumberPositionValues() {
		return numberPositionValues;
	}

	public boolean isLoadPositionValues() {
		return loadPositionValues;
	}

	public void setLoadPositionValues(boolean loadPositionValues) {
		this.loadPositionValues = loadPositionValues;
	}

	public boolean isMergeIdenticalPositions() {
		return mergeIdenticalPositions;
	}

	public void setMergeIdenticalPositions(
			boolean mergeIdenticalPositions) {
		this.mergeIdenticalPositions = mergeIdenticalPositions;
	}

	public DataInputStream getDis() {
		return dis;
	}

	public boolean isCloseDataInputStream() {
		return closeDataInputStream;
	}

	public void setCloseDataInputStream(boolean closeDataInputStream) {
		this.closeDataInputStream = closeDataInputStream;
	}

	public boolean isReplaceWithHitCount() {
		return replaceWithHitCount;
	}

	public void setReplaceWithHitCount(boolean replaceWithHitCount) {
		this.replaceWithHitCount = replaceWithHitCount;
	}
}
