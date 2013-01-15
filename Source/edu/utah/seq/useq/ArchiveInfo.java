package edu.utah.seq.useq;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.*;

/**This contains information related to all of the data slices in a USeq archive. It is always called archiveReadMe.xxx and is the first ZipEntry in the archive. 
 * The format of the archiveReadMe.txt version is simply comment lines beginning with '#' that are not parsed and key=values delimited by a return, thus one per line.
 * The first '=' sign in each key=value is used to split the tokens. Key's must not contain '=' signs or white space. White space before and after the '=' are permitted.
 * Where possible do use the reserve key names and add new ones as needed. At some point a archiveReadMe.xml version with a DTD should be created, volunteers?
 * 
 * @author david.nix@hci.utah.edu*/
public class ArchiveInfo {

	//required fields
	/**Version of the USeq Archive, current version is 1.0*/
	public static final String ARCHIVE_VERSION_KEY = "useqArchiveVersion";
	public static final String ARCHIVE_VERSION_VALUE_ONE = "1.0";
	public static final String ARCHIVE_README_NAME = "archiveReadMe.txt";

	/**Preferably of the DAS/2 form (e.g. H_sapiens_Mar_2006, C_elegans_May_2008)*/
	public static final String VERSIONED_GENOME_KEY = "versionedGenome";
	public static final Pattern DAS2_VERSIONED_GENOME_FORM = Pattern.compile("^\\w_\\w+_\\w+_\\d+$");

	/**Currently just two types exist graph and region.  This is used as to indicate how the data should be displayed (e.g. as a graph, or as 
	 * a regions/blocks). Others to add?*/
	public static final String DATA_TYPE_KEY = "dataType";
	public static final String DATA_TYPE_VALUE_GRAPH = "graph";
	public static final String DATA_TYPE_VALUE_REGION = "region";

	//reserved optional fields
	/**Where is the data from? Parsed from file x?*/
	public static final String ORIGINATING_DATA_SOURCE_KEY = "originatingDataSource";

	/**hat kind of graph style should be used to display, using IGB formats*/
	public static final String GRAPH_STYLE_KEY = "initialGraphStyle";
	public static final String GRAPH_STYLE_VALUE_BAR = "Bar";
	public static final String GRAPH_STYLE_VALUE_DOT = "Dot";
	public static final String GRAPH_STYLE_VALUE_LINE = "Line";
	public static final String GRAPH_STYLE_VALUE_MINMAXAVE = "Min_Max_Ave";
	public static final String GRAPH_STYLE_VALUE_STAIRSTEP = "Stairstep";
	public static final String GRAPH_STYLE_VALUE_HEATMAP = "HeatMap";

	/**What color, hex color values #0000FF*/
	public static final String COLOR_KEY = "initialColor";
	public static final String BACKGROUND_COLOR_KEY = "initialBackground";
	public static final Pattern COLOR_HEX_FORM = Pattern.compile("#\\w{6}");

	/**Initial minimum Y and maximum Y values to set for the data.*/
	public static final String MIN_Y_KEY = "initialMinY";
	public static final String MAX_Y_KEY = "initialMaxY";

	/**Free text description*/
	public static final String DESCRIPTION_KEY = "description";

	/**What units are the float values?*/
	public static final String UNIT_KEY = "units";

	/**USeq archive creation date. This will be added automatically.*/
	public static final String ARCHIVE_CREATION_DATE = "archiveCreationDate";

	//standard fields
	/**Comment lines, each beginning with '#' */
	private String[] commentLines = null;
	/**Container for all of the keyValues.*/
	private LinkedHashMap<String,String> keyValues = null;
	public static final Pattern KEY_VALUE_SPLITTER = Pattern.compile("\\s*([^=\\s]+)\\s*=\\s*(.+)\\s*");

	//constructors
	public ArchiveInfo(String versionedGenome, String dataType, boolean printGenomeVersionWarning){
		//instantiate keyValues Hash and add required fields
		keyValues = new LinkedHashMap<String,String>();
		keyValues.put(ARCHIVE_VERSION_KEY, ARCHIVE_VERSION_VALUE_ONE);
		//set date
		keyValues.put(ARCHIVE_CREATION_DATE, new Date().toString());
		//set dataType
		keyValues.put(DATA_TYPE_KEY, dataType);
		//set versioned genome
		keyValues.put(VERSIONED_GENOME_KEY, versionedGenome);
		//check to see if verisonedGenome follows form 
		if (printGenomeVersionWarning == true && DAS2_VERSIONED_GENOME_FORM.matcher(versionedGenome).matches() == false) System.err.println("\nWARNING: Versioned genome does not follow recommended form (e.g. H_sapiens_Mar_2006) correct -> "+versionedGenome);
	}
	public ArchiveInfo(File readMeTxtFile) throws IOException{
		loadTextArchiveReadMeFile(readMeTxtFile);
		//look for required fields
		if (keyValues.containsKey(ARCHIVE_VERSION_KEY) == false || keyValues.containsKey(VERSIONED_GENOME_KEY) == false || keyValues.containsKey(DATA_TYPE_KEY) == false){
			throw new IOException ("Error: text archiveReadMe.txt file does not contain required keys.  Add '"+ARCHIVE_VERSION_KEY+"' and or '"+VERSIONED_GENOME_KEY+"' and or '"+DATA_TYPE_KEY+"' to "+readMeTxtFile);
		}
		//check archive version
		if (keyValues.get(ARCHIVE_VERSION_KEY).equals(ARCHIVE_VERSION_VALUE_ONE) == false ){
			throw new IOException ("Error: this ArchiveInfo parser only supports "+ARCHIVE_VERSION_KEY+" = "+ARCHIVE_VERSION_VALUE_ONE);
		}
	}
	/**One way to get this is by zipFile.getInputStream(zipEntry)*/
	public ArchiveInfo(InputStream is, boolean closeStreams) {
		InputStreamReader isr = null;
		BufferedReader br = null;
		try {
			isr = new InputStreamReader(is);
			br = new BufferedReader(isr);
			loadTextArchiveReadMeFile(br);
		} catch (Exception e) {
			e.printStackTrace();
			USeqUtilities.safeClose(br);
			USeqUtilities.safeClose(isr);
		} finally {
			if (closeStreams) {
				USeqUtilities.safeClose(br);
				USeqUtilities.safeClose(isr);
			}
		}
	}

	//methods
	@SuppressWarnings("unchecked")
	public static ArchiveInfo fetchArchiveInfo(File useqArchive, boolean closeStream) {
		InputStream is = null;
		ArchiveInfo ai = null;
		try {
			if (USeqUtilities.USEQ_ARCHIVE.matcher(useqArchive.getName()).matches() == false) return null;
			ZipFile zf = new ZipFile(useqArchive);
			Enumeration e = zf.entries();
			ZipEntry ze = (ZipEntry) e.nextElement();
			is =  zf.getInputStream(ze);
			ai = new ArchiveInfo(is, closeStream);
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (closeStream) USeqUtilities.safeClose(is);
		}
		return ai;
	}


	/**Writes the comment lines and 'key = value' to file, one per line. Will overwrite.
	 * Returns the archiveReadMe.txt File object*/
	public File writeReadMeFile (File saveDirectory){
		//set date
		keyValues.put(ARCHIVE_CREATION_DATE, new Date().toString());
		PrintWriter out = null;
		try{
			File readme = new File (saveDirectory, ARCHIVE_README_NAME);
			out = new PrintWriter (new FileWriter (readme));
			//any comment lines?
			if (commentLines!= null){
				for (int i=0; i< commentLines.length; i++) out.println(commentLines[i]);
				out.println();
			}
			//print key values, spaces flanking = are permitted
			Iterator<String> it = keyValues.keySet().iterator();
			while (it.hasNext()){
				String key = it.next();
				String value = keyValues.get(key);
				out.println(key +" = "+value);
			}
			return readme;
		} catch (IOException e){
			e.printStackTrace();
			return null;
		} finally {
			USeqUtilities.safeClose(out);
		}
	}
	/**For appending Archive into onto a text file.*/
	public void appendCommentedKeyValues (PrintWriter out){
		//any comment lines?
		if (commentLines!= null){
			for (int i=0; i< commentLines.length; i++) out.println(commentLines[i]);
			out.println();
		}
		//print key values, spaces flanking = are permitted
		Iterator<String> it = keyValues.keySet().iterator();
		while (it.hasNext()){
			String key = it.next();
			String value = keyValues.get(key);
			out.println("# "+key +" = "+value);
		}
	}

	/**This does not close the BufferedReader.*/
	public void loadTextArchiveReadMeFile (BufferedReader in) {
		try {
			keyValues = new LinkedHashMap<String,String>();
			String line;
			ArrayList<String> comments = new ArrayList<String>();
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				if (line.startsWith("#")) comments.add(line);
				else {
					//split line
					Matcher mat = KEY_VALUE_SPLITTER.matcher(line);
					if (mat.matches() == false) throw new IOException("Error in parsing archiveReadMe.txt file. Found a non comment and non key = value line. Bad line -> '"+line);
					keyValues.put(mat.group(1), mat.group(2));
				}
			}
		} catch (IOException e) {
			e.printStackTrace();			
			USeqUtilities.safeClose(in);
		}
	}

	/**This does close the streams.*/
	public void loadTextArchiveReadMeFile (File readMeTxt) {
		FileReader fr = null;
		BufferedReader br = null;
		try {
			fr = new FileReader (readMeTxt);
			br = new BufferedReader (fr);
			loadTextArchiveReadMeFile (br);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} finally {
			USeqUtilities.safeClose(fr);
			USeqUtilities.safeClose(br);
		}
	}


	//getters and setters
	public void setArchiveVersion(String archiveVersion) {
		keyValues.put(ARCHIVE_VERSION_KEY, archiveVersion);
	}
	public void setOriginatingDataSource(String originatingDataSource) {
		keyValues.put(ORIGINATING_DATA_SOURCE_KEY, originatingDataSource);
	}
	public void setDataType(String dataType) {
		keyValues.put(DATA_TYPE_KEY, dataType);
	}
	public void setInitialGraphStyle(String initialGraphStyle) {
		keyValues.put(GRAPH_STYLE_KEY, initialGraphStyle);
	}
	public void setInitialColor(String initialColor) throws IOException{
		if (COLOR_HEX_FORM.matcher(initialColor).matches()== false) throw new IOException ("Error: initial color does not follow hex form (e.g. #B2B300)! "+initialColor);
		keyValues.put(COLOR_KEY, initialColor);
	}
	public void setInitialBackgroundColor(String initialBackgroundColor) throws IOException{
		if (COLOR_HEX_FORM.matcher(initialBackgroundColor).matches()== false) throw new IOException ("Error: initial background color does not follow hex form (e.g. #B2B300)! "+initialBackgroundColor);
		keyValues.put(BACKGROUND_COLOR_KEY, initialBackgroundColor);
	}
	public void setInitialMinY(String initialMinY) {
		keyValues.put(MIN_Y_KEY, initialMinY);
	}
	public void setInitialMaxY(String initialMaxY) {
		keyValues.put(MAX_Y_KEY, initialMaxY);
	}
	public void setDescription(String description) {
		keyValues.put(DESCRIPTION_KEY, description);
	}
	public void setUnits(String units) {
		keyValues.put(UNIT_KEY, units);
	}
	public String[] getCommentLines() {
		return commentLines;
	}
	public void setCommentLines(String[] commentLines) {
		this.commentLines = commentLines;
	}
	public String getValue(String key) {
		return keyValues.get(key);
	}
	public void setKeyValue(String key, String value){
		keyValues.put(key, value);
	}
	public LinkedHashMap<String,String> getKeyValues(){
		return keyValues;
	}
	/**Must contain ARCHIVE_VERSION_KEY, VERSIONED_GENOME_KEY, and DATA_TYPE_KEY .*/
	public void setKeyValues (LinkedHashMap<String,String> keyValues) throws IOException {
		if (keyValues.containsKey(ARCHIVE_VERSION_KEY) == false || keyValues.containsKey(VERSIONED_GENOME_KEY) == false || keyValues.containsKey(DATA_TYPE_KEY) == false){
			throw new IOException ("Error: keyValues do not contain required keys.  Add '"+ARCHIVE_VERSION_KEY+"' and or '"+VERSIONED_GENOME_KEY+"' and or '"+DATA_TYPE_KEY+"'");
		}
		this.keyValues = keyValues;
	}
	public String getVersionedGenome(){
		return keyValues.get(VERSIONED_GENOME_KEY);
	}
	public String getArchiveVersion(){
		return keyValues.get(ARCHIVE_VERSION_KEY);
	}
	public String getDataType(){
		return keyValues.get(DATA_TYPE_KEY);
	}
	public boolean isGraphData(){
		if (keyValues.get(DATA_TYPE_KEY).equals(DATA_TYPE_VALUE_GRAPH)) return true;
		return false;
	}
	public boolean isRegionData(){
		if (keyValues.get(DATA_TYPE_KEY).equals(DATA_TYPE_VALUE_REGION)) return true;
		return false;
	}
}
