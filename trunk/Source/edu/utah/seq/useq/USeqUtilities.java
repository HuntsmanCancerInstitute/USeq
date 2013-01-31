package edu.utah.seq.useq;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;
import java.util.*;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Passwords;

/** Utility class to hold constants for working with the USeq binary data type
 * @author david.nix@hci.utah.edu*/
public class USeqUtilities {
	//for building file types
	public static final String BYTE = "b";
	public static final String SHORT = "s";
	public static final String INT = "i";
	public static final String FLOAT = "f";
	public static final String DOUBLE = "d";
	public static final String TEXT = "t";
	public static final String BOOLEAN = "o";
	public static final String USEQ_EXTENSION_NO_PERIOD = "useq";
	public static final String USEQ_EXTENSION_WITH_PERIOD = ".useq";

	//possible binary file types
	//Position
	public static final Pattern POSITION = Pattern.compile("["+INT+SHORT+"]");
	public static final Pattern POSITION_INT = Pattern.compile(INT);
	public static final Pattern POSITION_SHORT = Pattern.compile(SHORT);
	//PositionScore
	public static final Pattern POSITION_SCORE = Pattern.compile("["+INT+SHORT+"]"+FLOAT);
	public static final Pattern POSITION_SCORE_INT_FLOAT = Pattern.compile(INT+FLOAT);
	public static final Pattern POSITION_SCORE_SHORT_FLOAT = Pattern.compile(SHORT+FLOAT);
	//PositionText
	public static final Pattern POSITION_TEXT = Pattern.compile("["+INT+SHORT+"]"+TEXT);
	public static final Pattern POSITION_TEXT_INT_TEXT = Pattern.compile(INT+TEXT);
	public static final Pattern POSITION_TEXT_SHORT_TEXT = Pattern.compile(SHORT+TEXT);
	//PositionScoreText
	public static final Pattern POSITION_SCORE_TEXT = Pattern.compile("["+INT+SHORT+"]"+FLOAT+TEXT);
	public static final Pattern POSITION_SCORE_TEXT_INT_FLOAT_TEXT = Pattern.compile(INT+FLOAT+TEXT);
	public static final Pattern POSITION_SCORE_TEXT_SHORT_FLOAT_TEXT = Pattern.compile(SHORT+FLOAT+TEXT);
	//Region
	public static final Pattern REGION = Pattern.compile("["+INT+SHORT+"]{2}");
	public static final Pattern REGION_INT_INT = Pattern.compile(INT+INT);
	public static final Pattern REGION_INT_SHORT = Pattern.compile(INT+SHORT);
	public static final Pattern REGION_SHORT_INT = Pattern.compile(SHORT+INT);
	public static final Pattern REGION_SHORT_SHORT = Pattern.compile(SHORT+SHORT);
	//RegionScore
	public static final Pattern REGION_SCORE = Pattern.compile("["+INT+SHORT+"]{2}"+FLOAT);
	public static final Pattern REGION_SCORE_INT_INT_FLOAT = Pattern.compile(INT+INT+FLOAT);
	public static final Pattern REGION_SCORE_INT_SHORT_FLOAT = Pattern.compile(INT+SHORT+FLOAT);
	public static final Pattern REGION_SCORE_SHORT_INT_FLOAT = Pattern.compile(SHORT+INT+FLOAT);
	public static final Pattern REGION_SCORE_SHORT_SHORT_FLOAT = Pattern.compile(SHORT+SHORT+FLOAT);
	//RegionText
	public static final Pattern REGION_TEXT = Pattern.compile("["+INT+SHORT+"]{2}"+TEXT);
	public static final Pattern REGION_TEXT_INT_INT_TEXT = Pattern.compile(INT+INT+TEXT);
	public static final Pattern REGION_TEXT_INT_SHORT_TEXT = Pattern.compile(INT+SHORT+TEXT);
	public static final Pattern REGION_TEXT_SHORT_INT_TEXT = Pattern.compile(SHORT+INT+TEXT);
	public static final Pattern REGION_TEXT_SHORT_SHORT_TEXT = Pattern.compile(SHORT+SHORT+TEXT);
	//RegionScoreText
	public static final Pattern REGION_SCORE_TEXT = Pattern.compile("["+INT+SHORT+"]{2}"+FLOAT+TEXT);
	public static final Pattern REGION_SCORE_TEXT_INT_INT_FLOAT_TEXT = Pattern.compile(INT+INT+FLOAT+TEXT);
	public static final Pattern REGION_SCORE_TEXT_INT_SHORT_FLOAT_TEXT = Pattern.compile(INT+SHORT+FLOAT+TEXT);
	public static final Pattern REGION_SCORE_TEXT_SHORT_INT_FLOAT_TEXT = Pattern.compile(SHORT+INT+FLOAT+TEXT);
	public static final Pattern REGION_SCORE_TEXT_SHORT_SHORT_FLOAT_TEXT = Pattern.compile(SHORT+SHORT+FLOAT+TEXT);

	//misc
	public static final Pattern USEQ_ARCHIVE = Pattern.compile("(.+)\\.(useq)$");

	//for GenoViz DAS/2
	public static final List<String> USEQ_FORMATS = new ArrayList<String>();
	static {
		USEQ_FORMATS.add(USEQ_EXTENSION_NO_PERIOD);
	}

	//static helper methods from the USeq project util.gen.Misc, util.gen.IO
	/**Rounds to an int and will set to 0 if neg or 1000 if > 1000 according to bed format requirements.*/
	public static int fixBedScore(float value){
		int score = (int)Math.round(value);
		if (score < 0) score = 0;
		else if (score > 1000) score = 1000;
		return score;
	}
	/**Prints message to screen, then exits.*/
	public static void printErrAndExit (String message){
		System.err.println (message);
		System.exit(0);
	}
	/**Removes an extension if found xxx.txt -> xxx.
	 * If none found returns the original.
	 */
	public static String removeExtension(String txt) {
		int index = txt.lastIndexOf(".");
		if (index != -1)  return txt.substring(0,index);
		return txt;
	}
	/**Returns a String separated by commas for each bin.*/
	public static String stringArrayToString(String[] s, String separator){
		if (s==null) return "";
		int len = s.length;
		if (len==1) return s[0];
		if (len==0) return "";
		StringBuffer sb = new StringBuffer(s[0]);
		for (int i=1; i<len; i++){
			sb.append(separator);
			sb.append(s[i]);
		}
		return sb.toString();
	}
	/**Returns a gz zip or straight file reader on the file based on it's extension compression.*/
	public static BufferedReader fetchBufferedReader( File txtFile) {
		BufferedReader in = null;
		try {
			String name = txtFile.getName().toLowerCase();
			if (name.endsWith(".zip")) {
				ZipFile zf = new ZipFile(txtFile);
				ZipEntry ze = (ZipEntry) zf.entries().nextElement();
				in = new BufferedReader(new InputStreamReader(zf.getInputStream(ze)));
			}
			else if (name.endsWith(".gz")) {
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(txtFile))));
			}
			else in = new BufferedReader (new FileReader (txtFile));
		} catch (Exception e) {
			e.printStackTrace();
			safeClose(in);
		} 
		return in;
	}
	/**Prints message to screen, then exits.*/
	public static void printExit (String message){
		System.out.println (message);
		System.exit(0);
	}
	/**Given a String of ints delimited by something, will parse or return null.*/
	public static int[] stringArrayToInts(String s, String delimiter) throws NumberFormatException{
		String[] tokens = s.split(delimiter);
		int[] num = new int[tokens.length];
		for (int i=0; i< tokens.length; i++){
			num[i] = Integer.parseInt(tokens[i]);
		}
		return num;
	}
	/**Extracts the full path file names of all the files and directories in a given directory. If a file is given it is
	 * returned as the File[0].
	 * Skips files starting with a '.'*/
	public static File[] extractFiles(File directory){
		File[] files = null;	
		String[] fileNames;
		if (directory.isDirectory()){
			fileNames = directory.list();
			int num = fileNames.length;
			ArrayList<File> al = new ArrayList<File>();
			try{
				String path = directory.getCanonicalPath();
				Pattern pat = Pattern.compile("^\\w+.*");
				Matcher mat; 
				for (int i=0; i< num; i++)  {
					mat = pat.matcher(fileNames[i]);
					if (mat.matches()) al.add(new File(path, fileNames[i]));
				}
				//convert arraylist to file[]
				if (al.size() != 0){
					files = new File[al.size()];
					al.toArray(files);
				}
			}catch(IOException e){
				System.out.println("Problem extractFiles() "+directory);
				e.printStackTrace();
				return null;
			}
		}
		if (files == null){
			files = new File[1];
			files[0] = directory;
		}
		Arrays.sort(files);
		return files;
	}
	
	/**Fetches all files with a given extension in a directory recursing through sub directories.
	 * Will return a file if a file is given with the appropriate extension, or null.*/
	public static File[] fetchFilesRecursively (File directory, String extension){
		if (directory.isDirectory() == false){
			return extractFiles(directory, extension);
		}
		ArrayList<File> al = fetchAllFilesRecursively (directory, extension);
		File[] files = new File[al.size()];
		al.toArray(files);
		return files;
	}
	
	/**Fetches all files with a given extension in a directory recursing through sub directories.*/
	public static ArrayList<File> fetchAllFilesRecursively (File directory, String extension){
		ArrayList<File> files = new ArrayList<File>(); 
		File[] list = directory.listFiles();
		for (int i=0; i< list.length; i++){
			if (list[i].isDirectory()) {
				ArrayList<File> al = fetchAllFilesRecursively (list[i], extension);
				files.addAll(al);
			}
			else{
				if (list[i].getName().endsWith(extension)) files.add(list[i]);
			}
		}
		return files;
	}
	
	/**Extracts the full path file names of all the files in a given directory with a given extension (ie txt or .txt).
	 * If the dirFile is a file and ends with the extension then it returns a File[] with File[0] the
	 * given directory. Returns null if nothing found. Case insensitive.*/
	public static File[] extractFiles(File dirOrFile, String extension){
		if (dirOrFile == null) return null;
		File[] files = null;
		Pattern p = Pattern.compile(".*"+extension+"$", Pattern.CASE_INSENSITIVE);
		Matcher m;
		if (dirOrFile.isDirectory()){
			files = dirOrFile.listFiles();
			int num = files.length;
			ArrayList<File> chromFiles = new ArrayList<File>();
			for (int i=0; i< num; i++)  {
				m= p.matcher(files[i].getName());
				if (m.matches()) chromFiles.add(files[i]);
			}
			files = new File[chromFiles.size()];
			chromFiles.toArray(files);
		}
		else{
			m= p.matcher(dirOrFile.getName());
			if (m.matches()) {
				files=new File[1];
				files[0]= dirOrFile;
			}
		}
		if (files != null) Arrays.sort(files);
		return files;
	}
	/**Takes a file, capitalizes the text, strips off .gz or .zip and any other extension, then makes a directory of the file and returns it.*/
	public static File makeDirectory (File file, String extension){
		String name = capitalizeRemoveExtension(file, extension);
		File dir = new File (file.getParentFile(), name+ extension);
		dir.mkdir();
		return dir;
	}
	public static String capitalizeRemoveExtension(File file, String extension){
		String name = file.getName();
		name = capitalizeFirstLetter(name);
		name = name.replace(".gz", "");
		name = name.replace(".zip", "");
		name = removeExtension(name);
		return name;
	}
	/**Capitalizes the first letter in a String.*/
	public static String capitalizeFirstLetter(String s){
		char[] first = s.toCharArray();		
		if (Character.isLetter(first[0])){
			first[0] = Character.toUpperCase(first[0]);
			return new String(first);
		}
		return s;
	}
	/**Deletes a directories contents, and then the directory.
	 * Returns false if at any time a file cannot be deleted or the directory is null.*/
	public static boolean deleteDirectory(File directory){
		if (directory == null) return false;
		File[] files = directory.listFiles();
		int num = files.length;
		for (int i=0; i<num; i++){
			if (files[i].isDirectory()) {
				if (deleteDirectory(files[i]) == false) return false;
			}
			else if (files[i].delete() == false) return false;
		}
		if (directory.delete() == false)return false;
		return true;
	}
	/**Zip compresses an array of Files, be sure to text your zipFile with a .zip extension!*/
	public static boolean zip(File[] filesToZip, File zipFile ){
		byte[] buf = new byte[2048];
		ZipOutputStream out = null;
		FileInputStream in = null;
		try {
			out = new ZipOutputStream(new FileOutputStream(zipFile));		
			// Compress the files
			for (int i=0; i<filesToZip.length; i++) {
				in = new FileInputStream(filesToZip[i]);
				out.putNextEntry(new ZipEntry(filesToZip[i].getName()));
				int len;
				while ((len = in.read(buf)) != -1) {
					out.write(buf, 0, len);
				}
				out.closeEntry();
				in.close();
			}
			out.close();
		} catch (IOException e) {	
			System.err.println("Can't zip()");
			e.printStackTrace();
			safeClose(out);
			safeClose(in);
			return false;
		} 
		return true;
	}
	/**Merges all files in File[][] to a File[].*/
	public static File[] collapseFileArray(File[][] f){
		ArrayList<File> al = new ArrayList<File>();
		for (int i=0; i< f.length; i++){
			if (f[i] != null){
				for (int j=0; j< f[i].length; j++){
					al.add(f[i][j]);
				}
			}
		}
		File[] files = new File[al.size()];
		al.toArray(files);
		return files;
	}
	/**Using interbase coordinates so length = stop - start.*/
	public static int calculateMiddleIntergenicCoordinates(int start, int end){
		if (start == end) return start;
		double length = end - start;
		double halfLength = length/2.0;
		return (int)Math.round(halfLength) + start;
	}
	
	/**Converts #3333CC to RGB values 51 51 204*/
	public static String convertHexadecimal2RGB(String hexColor, String divider){
		String hex = hexColor.replaceFirst("#", "");
		if (hex.length() != 6) return null;
		StringBuilder sb = new StringBuilder();
		String toCon = hex.substring(0,2);
		int color= Integer.parseInt(toCon,16);
		sb.append(color);
		for (int i=2; i<6; i+=2){
			toCon = hex.substring(i,i+2);
			color= Integer.parseInt(toCon,16);
			sb.append(divider);
			sb.append(color);
		}
		return sb.toString();
	}
	
	/**Executes tokenized params on command line, use full paths.
	 * Put each param in its own String.  
	 * Returns null if a problem is encountered.
	 * Returns both error and data from execution.
	 */
	public static String[] executeCommandLineReturnAll(String[] command){
		ArrayList<String> al = new ArrayList<String>();		
		try {
			Runtime rt = Runtime.getRuntime();
			rt.traceInstructions(true); //for debugging
			rt.traceMethodCalls(true); //for debugging
			Process p = rt.exec(command);
			//Process p = rt.exec(Misc.stringArrayToString(command, " "));
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			String line;
			while ((line = data.readLine()) != null){
				al.add(line);
			}
			while ((line = error.readLine()) != null){
				al.add(line);
			}
			data.close();
			error.close();

		} catch (Exception e) {
			System.out.println("Problem executing -> "+Misc.stringArrayToString(command," ")+" "+e.getLocalizedMessage());
			e.printStackTrace();
			return null;
		}
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}
	
	/**Executes a String of shell script commands via a temp file.  Only good for Unix.
	 * @throws IOException */
	public static String[] executeShellScript (String shellScript, File tempDirectory) throws IOException{
		//make shell file
		File shellFile = new File (tempDirectory, "tempFile_"+Passwords.createRandowWord(10) +".sh");
		String fullPath = IO.getFullPathName(shellFile);
		//write to shell file
		IO.writeString(shellScript, shellFile);
		//set permissions for execution
		String[] cmd = {"chmod", "777", fullPath};
		String[] res = IO.executeCommandLineReturnAll(cmd);
		if (res == null || res.length !=0 ) {
			shellFile.delete();
			throw new IOException ("Error found making shellscript executable : "+shellScript);
		}
		//execute
		cmd = new String[]{fullPath};
		res = IO.executeCommandLineReturnAll(cmd);
		shellFile.delete();
		return res; 
	}
	
	/**
	 * Safely close a Closeable object.  If it doesn't exist, return.
	 */
	public static <S extends Closeable> void safeClose(S s) {
		if (s == null) {
			return;
		}
		try {
			s.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}
