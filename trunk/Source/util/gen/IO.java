package util.gen;
import java.io.*;

import javax.servlet.http.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.*;
import java.util.jar.JarFile;
import java.util.jar.Manifest;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.*;
import org.apache.tools.bzip2.*;
import java.net.*;


/**
 * Static methods for Input Output related tasks.
 */
public class IO {

	public static final Pattern COMMA = Pattern.compile(",");

	/**Writes a String, max length 256.*/
	public static void writeStringShort(String x, DataOutputStream out) throws IOException{
		int length = x.length();
		if (length > 256) throw new IOException("String length exceeds 256 characters. See "+x);
		length = length -128;
		//write length - 128
		out.writeByte(length);
		//write string
		out.writeBytes(x);
	}

	/**Reads a String, max length 256.*/
	public static String readStringShort(DataInputStream dis) throws IOException{
		//read length
		byte[] barray = new byte[dis.readByte() + 128];
		dis.readFully(barray);
		return new String(barray);
	}

	/**Attempts to fetch the 'Implementation-Version' from the calling jar file's manifest.*/
	public static String fetchUSeqVersion(){
		String version = "";
		try {
			String jarPath = IO.class.getProtectionDomain().getCodeSource().getLocation().getPath();
			jarPath = URLDecoder.decode(jarPath, "UTF-8");
			if (jarPath.endsWith(".jar") == false) return version;
			JarFile jar = new JarFile(new File (jarPath));
			Manifest manifest = jar.getManifest();
			version = manifest.getMainAttributes().getValue("Implementation-Version");
			if (version == null) version = "";
		} catch (Exception x) {
			System.err.println("\nProblem fetching version information from the jar file.\n");
			x.printStackTrace();  
		}
		return version;
	}

	/**Returns the names of the files or directories.*/
	public static String[] fetchFileNames(File[] files){
		String[] names = new String[files.length];
		for (int i=0; i< files.length; i++){
			names[i] = files[i].getName();
		}
		return names;
	}

	/**Returns the names of the files or directories sans their extension*/
	public static String[] fetchFileNamesNoExtensions(File[] files){
		String[] names = new String[files.length];
		for (int i=0; i< files.length; i++){
			names[i] = Misc.removeExtension(files[i].getName());
		}
		return names;
	}

	/**Returns the amount of memory being used.*/
	public static String memoryUsed(){
		System.gc();
		Runtime rt = Runtime.getRuntime();
		System.gc();
		long usedMB = (rt.totalMemory() - rt.freeMemory()) / 1024 / 1024;
		return usedMB+"MB @ "+System.currentTimeMillis()/1000+"Sec";
	}

	/**Checks to see if the current java is 1.6, 1.7, or 1.8 by executing and parsing 'java -version'*/
	public static boolean checkJava(){
		String[] cmd = {"java","-version"};
		String[] results = IO.executeCommandLineReturnAll(cmd);
		if (results == null || results.length ==0) return false;
		for (int i=0; i< results.length; i++){
			if (results[i].startsWith("java version")){
				if (results[i].contains("1.6.") || results[i].contains("1.7.")|| results[i].contains("1.8.")) return true;
				return false;
			}
		}
		return false;
	}

	/**Returns null is something bad happened, otherwise returns the r output lines containing 'error' case insensitive or returns '' if no error found.*/
	public static String runRCommandLookForError(String rCmdLine, File rApp, File tempDir){
		//look for R
		if (rApp.canExecute() == false) return "Cannot execute R check -> "+rApp;
		String rndWrd = Passwords.createRandowWord(7);
		File inFile = new File (tempDir, "RCmdLine_"+rndWrd);
		File outFile = new File (tempDir, "ROutput_"+rndWrd);
		String[] command = {
				rApp.toString(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore", 
				inFile.toString(),
				outFile.toString()
		};
		if (IO.writeString(rCmdLine, inFile) == false) return null;
		String[] messages = IO.executeCommandLineReturnAll(command);
		inFile.delete();
		if (messages == null || messages.length !=0){
			Misc.printArray(messages);
			outFile.delete();
			return null;
		}
		//load results
		String[] rLines = IO.loadFile(outFile);
		outFile.delete();
		//look for error
		boolean errorFound = false;
		StringBuilder errors = new StringBuilder();
		Pattern errorPat = Pattern.compile("error", Pattern.CASE_INSENSITIVE);
		for (int i=0; i< rLines.length; i++){
			if (errorFound) {
				errors.append(" ");
				errors.append(rLines[i]);
			}
			else if (errorPat.matcher(rLines[i]).lookingAt()){
				errors.append(rLines[i]);
				errorFound = true;
			}
		}
		if (errorFound) return errors.toString();
		return "";
	}

	/**Takes a file, capitalizes the text, strips off .gz or .zip and any other extension, then makes a directory of the file and returns it.
	 * returns null if the new directory exists.*/
	public static File makeDirectory (File file, String extension){
		String name = file.getName();
		name = Misc.capitalizeFirstLetter(name);
		name = name.replace(".gz", "");
		name = name.replace(".zip", "");
		name = Misc.removeExtension(name);
		File dir = new File (file.getParentFile(), name+ extension);
		if (dir.exists()) return null;
		dir.mkdir();
		return dir;
	}

	/**Takes a file strips off .gz or .zip and returns it's extension without the leading . or "".*/
	public static String fetchFileExtension (File file){
		String name = file.getName();
		name = name.replace(".gz", "");
		name = name.replace(".zip", "");
		int index = name.lastIndexOf(".");
		if (index != -1)  return name.substring(index+1);
		return "";
	}

	/**Uses Ants implementation of the bzip2 compression algorithm.
	 * Be sure to text your file xxx.bz2*/
	public static boolean bzip2(File toZip, File zipped) {
		try {
			//output
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(zipped));
			bos.write('B');
			bos.write('Z');
			CBZip2OutputStream zOut = new CBZip2OutputStream(bos);
			//input
			FileInputStream in = new FileInputStream(toZip);
			//zip it
			byte[] buffer = new byte[8 * 1024];
			int count = 0;
			do {
				zOut.write(buffer, 0, count);
				count = in.read(buffer, 0, buffer.length);
			} while (count != -1);
			//close streams
			zOut.close();
			in.close();
			return true;
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return false;
		}
	}

	/**Uncompressed bzip2 files using the Ant implementation.
	 * Be sure to text your file xxx.bz2*/
	public static boolean bunzip2(File zippedFile, File decompressedFile) {
		try {
			FileOutputStream out = new FileOutputStream(decompressedFile);
			FileInputStream fis = new FileInputStream(zippedFile);
			BufferedInputStream bis = new BufferedInputStream(fis);
			int b = bis.read();
			int z = bis.read();
			if (b != 'B' && z != 'Z') {
				System.out.println("Invalid bz2 file! "+ zippedFile);
				return false;
			}
			CBZip2InputStream zIn = new CBZip2InputStream(bis);
			byte[] buffer = new byte[8 * 1024];
			int count = 0;
			do {
				out.write(buffer, 0, count);
				count = zIn.read(buffer, 0, buffer.length);
			} while (count != -1);
			out.close();
			fis.close();
			return true;
		} catch (IOException ioe) {
			System.out.println("Problem bunzip2ing your file! ");
			ioe.printStackTrace();
		}
		return false;
	}

	/**Converts a String of "grp1=/my/Dir1,grp2=/my/Dir2, etc that contains 
	 * directories into a LinkedHashMap.
	 * Returns null and prints an error message if an error was thrown, ie no directory, not a directory
	 */
	public static LinkedHashMap buildDirectoryFileGroups (String map){
		//convert key=value into a hashmap
		String[] lines = map.split(",");
		LinkedHashMap hash = new LinkedHashMap();
		for (int i=0; i<lines.length; i++){
			String[] keyValue= lines[i].split("=");
			//check if keyValue exists
			if (keyValue.length != 2) {
				System.out.println("\nError: key=value pair not found? -> "+lines[i]+"\n");
				return null;
			}
			//check if file exists and is readable
			File file = new File (keyValue[1]);
			if (file.canRead() == false || file.isDirectory() == false){
				System.out.println("\nError: could not read this directory -> "+file+"\n");
				return null;
			}
			//check to see if key already exists
			if (hash.containsKey(keyValue[0])){
				System.out.println("\nError: duplicate key found -> "+lines[i]+"\n");
				return null;
			}
			//add
			hash.put(keyValue[0], file);

		}
		return hash;
	}

	/**Converts a String of "grp1=/my/file1,grp1=/my/old/file2,grp2=/my/file2,
	 * grp2=/my/new/file3,grp3=/my/new/dir/,/my/default etc"
	 * to a File[][] where the files have been broken by grouping into different File[].
	 * If no key is provided, the files are assigned to a common 'default' grouping.
	 * If a directory is provided all files within are grouped together.
	 * Mixing of keyed and unkeyed is permitted.
	 * Returns null if an error was thrown.
	 * @param directoryFileExtensionFilter - if directories are given in the map, you can 
	 * filter which files to extract (ie 'cel'), case insensitive, set to "." for everything.
	 * @return - ArrayList containing String[] groupNames and a File[][] of grouped files, and
	 * a LinkedHashMap of the two.*/
	public static ArrayList buildFileGroups (String map, String directoryFileExtensionFilter){
		//convert key=value into a hashmap
		if (map == null){
			System.out.println("\nError: please enter files to check.\n");
			return null;
		}
		String[] lines = map.split(",");
		LinkedHashMap hash = new LinkedHashMap();
		for (int i=0; i<lines.length; i++){
			String[] keyValue= lines[i].split("=");
			//check if keyValue exists
			if (keyValue.length == 1) {
				keyValue = new String[]{"Default", lines[i]};
			}
			else if (keyValue.length > 2){
				System.out.println("\nError: splitting of your grouping=filename failed?!  Problem line -> "+lines[i]+"\n");
				return null;
			}
			//check if file exists and is readable
			File file = new File (keyValue[1]);
			if (file.canRead() == false){
				System.out.println("\nError: could not read this file/ directory -> "+file+"\n");
				return null;
			}
			//is it a directory?
			if (file.isDirectory()){
				//get cel files
				File[] all = IO.extractFiles(file,directoryFileExtensionFilter);
				Arrays.sort(all);
				//load hash
				//does the grouping already exist?
				if (hash.containsKey(keyValue[0])){
					ArrayList files = (ArrayList)hash.get(keyValue[0]);
					for (int x=0; x<all.length; x++) files.add(all[x]);
				}
				else {
					ArrayList files = new ArrayList();
					for (int x=0; x<all.length; x++) files.add(all[x]);
					hash.put(keyValue[0],files);
				}
			}
			//no it is not a directory
			else {
				//does the grouping already exist?
				if (hash.containsKey(keyValue[0])){
					ArrayList files = (ArrayList)hash.get(keyValue[0]);
					files.add(file);
				}
				else {
					ArrayList files = new ArrayList();
					files.add(file);
					hash.put(keyValue[0],files);
				}
			}
		}

		//convert hash to File[][]
		File[][] files = new File[hash.size()][];
		String[] groupNames = new String[hash.size()];
		Iterator it = hash.keySet().iterator();
		int counter = 0;
		while (it.hasNext()){
			String key = (String)it.next();
			groupNames[counter] = key;
			ArrayList al = (ArrayList)hash.get(key);
			File[] fileGrp = new File[al.size()];
			al.toArray(fileGrp);	
			//sort
			Arrays.sort(fileGrp);
			files[counter++] = fileGrp;
			//set in hash replacing original
			hash.put(key, fileGrp);
		}
		ArrayList results = new ArrayList(3);
		results.add(groupNames);		
		results.add(files);
		results.add(hash);
		return results;
	}

	/**
	 * Use to read in xxx.bz2 files directly to a String[].
	 * Returns null if a problem encountered.
	 * Must have bunzip2 installed.
	 */
	public static String[] bunZip2(File bunZip2App, File bz2FileToUncompress) {
		String[] command = {IO.getFullPathName(bunZip2App), "--stdout", "--decompress", IO.getFullPathName(bz2FileToUncompress)};
		return IO.executeCommandLine(command);
	}

	/**Check to see that a file exists and is not a directory.*/
	public static void checkFile(String par){
		File f = new File(par);
		if (f.exists()) return;
		if (f.isDirectory()) return;
		System.out.println("\nSorry, I can't find your file or directory! ->"+par);
		System.exit(0);
	}

	/**Merges all files in File[][] to a File[].*/
	public static File[] collapseFileArray(File[][] f){
		ArrayList al = new ArrayList();
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

	/**Concatinates the full path file names for the given File array.*/
	public static String concatinateFileFullPathNames (File[] f, String seperator){
		try{
			String[] names = new String[f.length];
			for (int i=0; i< f.length; i++) names[i] = f[i].getCanonicalPath();
			return Misc.stringArrayToString(names, seperator);
		}
		catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}
	/** Fast & simple file copy. From GForman http://www.experts-exchange.com/M_500026.html
	 * Hit an odd bug with a "Size exceeds Integer.MAX_VALUE" error when copying a vcf file. -Nix.*/
	public static boolean copyViaFileChannel(File source, File dest){
		FileChannel in = null, out = null;
		try {
			in = new FileInputStream(source).getChannel();
			out = new FileOutputStream(dest).getChannel();
			long size = in.size();
			MappedByteBuffer buf = in.map(FileChannel.MapMode.READ_ONLY, 0, size);
			out.write(buf);
			if (in != null) in.close();
			if (out != null) out.close();
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	/** Copy via read line. Slower than channel copy but doesn't have Integer.MAX_Value problem.
	 * Source file can be txt, txt.gz, or txt.zip. So decompresses if needed.
	 * @author Nix*/
	public static boolean copyViaReadLine(File source, File dest){
		BufferedReader in = null;
		PrintWriter out = null;
		try {
			in = IO.fetchBufferedReader(source);
			out = new PrintWriter (new FileWriter (dest));
			String line;
			while ((line = in.readLine()) != null) out.println(line);
			return true;
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		} finally {
				try {
					if (in != null) in.close();
					if (out != null) out.close();
				} catch (IOException e) {}
			
		}
	}
	
	
	
	/**Attempts to uncompress a xxx.gz or xxx.zip file and write it to the same location without the extension.  
	 * Returns null if any issues are encountered. If the uncompressed file already exists, it is returned.
	 * If uncompressed file is null then it is created from the gzipOrZipCompressed's name.*/
	public static File uncompress (File gzipOrZipCompressed, File uncompressed) {

		File f = null;
		PrintWriter out = null;
		BufferedReader in = null;

		try {
			String name = gzipOrZipCompressed.getName();
			if (name.endsWith(".zip")) name = name.substring(0, name.length()-4);
			else if (name.endsWith(".gz")) name = name.substring(0, name.length()-3);
			else return null;
			
			if (uncompressed != null) f = uncompressed;
			else {
				f = new File (gzipOrZipCompressed.getParentFile(), name);
				if (f.exists() && f.length() > 40) return f;
			}
			out = new PrintWriter( new FileWriter(f));
			in = IO.fetchBufferedReader(gzipOrZipCompressed);
			String line;
			while ((line = in.readLine()) != null) out.println(line);
			
			in.close();
			out.close();
			return f;
		} catch (IOException e){
			e.printStackTrace();
			return null;
		} finally {
			try {
				if (out != null) out.close();
				if (in != null) in.close();
			} catch (IOException e){
				e.printStackTrace();
			}
		}
	}

	/**Copies a given directory and it's contents to the destination directory.
	 * Use a extension (e.g. "class$|properties$") to limit the files copied over or set to null for all.*/
	public static boolean copyDirectoryRecursive (File sourceDir, File destDir, String extension){
		Pattern pat = Pattern.compile(extension);
		if (destDir.exists() == false) destDir.mkdir();
		//for each file in source copy to destDir		
		File[] files = IO.extractFiles(sourceDir);
		for (int i=0; i< files.length; i++){
			if (files[i].isDirectory()) {
				copyDirectoryRecursive(files[i], new File (destDir, files[i].getName()), extension);
			}
			else {
				Matcher mat = pat.matcher(files[i].getName());
				if (extension == null || mat.find()){
					File copied = new File (destDir, files[i].getName());					
					if (copyViaFileChannel(files[i], copied) == false ) return false;
				}
			}
		}
		return true;
	}

	/**Counts the number of lines in a file skipping blanks.*/
	public static long countNonBlankLines(File file){
		long num =0;
		try {
			BufferedReader in = new BufferedReader(new FileReader(file));
			String line;
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() !=0)num++;
			}
		}
		catch (IOException e){
			System.out.println("\nProblem counting the number of lines in the file: "+file);
			e.printStackTrace();
		}
		return num;
	}

	/**Counts the number of lines in a file*/
	public static long countNumberOfLines(File file){
		long num =0;
		try {
			BufferedReader in = new BufferedReader(new FileReader(file));
			while ((in.readLine()) !=null) {
				num++;
			}
		}
		catch (IOException e){
			System.out.println("\nProblem counting the number of lines in the file: "+file);
			e.printStackTrace();
		}
		return num;
	}

	/**Counts the number of File in the array that actually exist.*/
	public static int countNumberThatExist(File[] f){
		int num =0;
		for (int i=0; i< f.length; i++) if (f[i].exists()) num++;
		return num;
	}

	/**Attempts to delete a directory and it's contents.
	 * Returns false if all the file cannot be deleted or the directory is null.
	 * Files contained within scheduled for deletion upon close will cause the return to be false.*/
	public static void deleteDirectory(File dir){
		if (dir == null || dir.exists() == false) return;
		if (dir.isDirectory()) {
			File[] children = dir.listFiles();
			for (int i=0; i<children.length; i++) {
				deleteDirectory(children[i]);
			}
			dir.delete();
		}
		dir.delete();
	}

	public static void deleteDirectoryViaCmdLine(File dir){
		try {
			IO.executeCommandLine(new String[]{"rm","-rf",dir.getCanonicalPath()});
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**Attempts to delete a file.*/
	public static boolean deleteFile(File file) {
		String name = ""; 
		try {
			name = file.getCanonicalPath();
			if (!file.delete()) {
				System.out.println("Warning: could not delete the file " + name);
				return false;
			}
			return true;
		}

		catch (Exception ex){
			System.out.println("Problem with not deleteFile() " + name);
			ex.printStackTrace();
			return false;
		}
	}

	/**Attempts to delete a file.*/
	public static boolean deleteFile (String fullPathFileName){
		return deleteFile(new File(fullPathFileName));
	}


	/**Deletes files in a given directory with a given extension.*/
	public static boolean deleteFiles(File directory, String extension){		
		boolean deleted = true;
		if (directory.isDirectory()){
			String[] files = directory.list();
			int num = files.length;
			try{
				String path = directory.getCanonicalPath();
				for (int i=0; i< num; i++)  {
					if (files[i].endsWith(extension)) {
						deleted = new File(path,files[i]).delete();
						if (deleted ==false) return false;
					} 
				}
			}catch(IOException e){
				System.out.println("Prob deleteFiles, dir ->"+directory+" "+extension);
				e.printStackTrace();
			}
		}
		else deleted = false;
		return deleted;
	}

	/**Deletes files in a given directory with a given prefix, that are older cutoff.*/
	public static void deleteFiles(File directory, String prefix, int minutes){		
		if (directory.isDirectory()){
			long cutoff = minutes*1000*60;
			long current = System.currentTimeMillis();
			cutoff = current - cutoff;
			File[] files = directory.listFiles();
			int num = files.length;
			for (int i=0; i< num; i++)  {
				if (files[i].getName().startsWith(prefix) && files[i].lastModified()< cutoff) {
					files[i].delete();
				} 
			}
		}
	}

	/**Attempts to delete the array of File.*/
	public static boolean deleteFiles (File[] files){
		if (files == null) return false;
		for (int i=0; i< files.length; i++){
			boolean deleted = files[i].delete();
			if (deleted == false) return false;
		}
		return true;
	}

	/**Looks in the given directory for files ending with given String and deletes them.*/
	public static void deleteFiles (String directory, String endsWith){
		File dir = new File (directory);
		String[] files = dir.list();
		int numFiles = files.length;
		for (int i=0; i< numFiles; i++){
			if (files[i].endsWith(endsWith)){
				new File (directory, files[i]).delete();
			}
		}
	}

	/**Deletes files in a given directory that don't stop in a given extension.*/
	public static void deleteFilesNotEndingInExtension(File directory, String extension){		
		File[] files = directory.listFiles();
		for (int i=0; i < files.length; i++){
			if (files[i].getName().endsWith(extension) == false) files[i].delete();
		}
	}

	/**Executes tokenized params on command line, use full paths.
	 * Put each param in its own String.  
	 * Returns null if a problem is encountered.
	 */
	public static String[] executeCommandLine(String[] command){
		ArrayList al = new ArrayList();
		try {
			Runtime rt = Runtime.getRuntime();
			//rt.traceInstructions(true); //for debugging
			//rt.traceMethodCalls(true); //for debugging
			Process p = rt.exec(command);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			//BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			String line;
			while ((line = data.readLine()) != null){
				al.add(line);
			}
			//while ((line = error.readLine()) != null){
			//	System.out.println(line);
			//}
			data.close();
			//error.close();

		} catch (Exception e) {
			System.out.println("Problem executingCommandLine(), command -> "+Misc.stringArrayToString(command," "));
			e.printStackTrace();
			return null;
		}
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}

	/**Executes tokenized params on command line, use full paths.
	 * Put each param in its own String.  
	 * Returns null if a problem is encountered. Does not print any stack trace. Returns the output to standard out and standard error.
	 */
	public static String[] executeCommandLineNoError(String[] command){
		ArrayList<String> al = new ArrayList<String>();
		try {
			Runtime rt = Runtime.getRuntime();
			Process p = rt.exec(command);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
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
			return null;
		}
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}

	/**Executes tokenized params on command line, use full paths.
	 * Put each param in its own String.  
	 * Writes lines to a file.
	 */
	public static String[] executeCommandLine(String[] command, File txtOutPutFile){
		ArrayList al = new ArrayList();
		try {
			Runtime rt = Runtime.getRuntime();
			//rt.traceInstructions(true); //for debugging
			//rt.traceMethodCalls(true); //for debugging
			Process p = rt.exec(command);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			//BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			PrintWriter out = new PrintWriter( new FileWriter(txtOutPutFile));
			String line;
			while ((line = data.readLine()) != null){
				out.println(line);
			}
			//while ((line = error.readLine()) != null){
			//	System.out.println(line);
			//}
			out.close();
			data.close();
			//error.close();

		} catch (Exception e) {
			System.out.println("Problem executingCommandLine(), command -> "+Misc.stringArrayToString(command," "));
			e.printStackTrace();
		}
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}

	/**Executes tokenized params on command line, use full paths.
	 * Put each param in its own String.  
	 * Writes lines to a file.
	 */
	public static String[] executeCommandLine(String[] command, File txtOutPutFile, String[] env){
		ArrayList al = new ArrayList();
		try {
			Runtime rt = Runtime.getRuntime();
			rt.traceInstructions(true); //for debugging
			rt.traceMethodCalls(true); //for debugging
			Process p = rt.exec(command, env);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			PrintWriter out = new PrintWriter( new FileWriter(txtOutPutFile));
			String line;
			int counter = 0;
			while ((line = error.readLine()) != null){
				System.out.println(line);
			}

			while ((line = data.readLine()) != null){
				out.println(line);
				if (counter++ == 10) break;
			}

			out.close();
			data.close();
			error.close();

		} catch (Exception e) {
			System.out.println("Problem executingCommandLine(), command -> "+Misc.stringArrayToString(command," "));
			e.printStackTrace();
		}
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
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

	/**Executes a String of shell script commands via a temp file.  Only good for Unix.*/
	public static String[] executeShellScript (String shellScript, File tempDirectory){
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
			return null;
		}
		//execute
		cmd = new String[]{fullPath};
		res = IO.executeCommandLineReturnAll(cmd);
		shellFile.delete();
		return res; 
	}

	/**Takes a String, possibly comma delimited, no spaces.  Looks for URLs defined by the http prefix, otherwise assumes it's a file,
	 * fetches all files recursively if it is a directory.  Set the extension to null if you want everything.  "." hidden files are ignored.
	 * Returns null if something bad happened or no files were returned.*/
	public static String[] fetchFileURLStrings(String s, String extension){
		//attempt to split on ,
		String[] names = COMMA.split(s);
		//add or attempt to pull file names
		ArrayList<String> toReturn = new ArrayList<String>();
		for (int i=0; i< names.length; i++){
			//is it a url?
			if (names[i].startsWith("http")) {
				if (extension !=null && names[i].endsWith(extension)) toReturn.add(names[i]);
				else toReturn.add(names[i]);
			}
			else {
				File[] f = null;
				if (extension != null) f = fetchFilesRecursively(new File(names[i]), extension);
				else f = fetchFilesRecursively(new File(names[i]));
				if (f != null){
					for (int x=0; x< f.length; x++) toReturn.add(f[x].toString());
				}
			}
		}
		int num = toReturn.size();
		if (num == 0) return null;
		names = new String[num];
		toReturn.toArray(names);
		return names;
	}

	/**Extracts the full path file names of all the files and directories in a given directory. If a file is given it is
	 * returned as the File[0].
	 * Skips files starting with a '.'*/
	public static File[] extractFiles(File directory){
		try{
			directory = directory.getCanonicalFile();
			File[] files = null;	
			String[] fileNames;
			if (directory.isDirectory()){
				fileNames = directory.list();
				int num = fileNames.length;
				ArrayList<File> al = new ArrayList<File>();

				String path = directory.getCanonicalPath();
				for (int i=0; i< num; i++)  {
					if (fileNames[i].startsWith(".") == false) al.add(new File(path, fileNames[i])); 
				}
				//convert arraylist to file[]
				if (al.size() != 0){
					files = new File[al.size()];
					al.toArray(files);
				}
			}
			if (files == null){
				files = new File[1];
				files[0] = directory;
			}
			Arrays.sort(files);
			return files;

		}catch(IOException e){
			System.out.println("Problem extractFiles() "+directory);
			e.printStackTrace();
			return null;
		}
	}

	/**Returns directories or null if none found. Not recursive.
	 * Skips those beginning with a period.*/
	public static File[] extractOnlyDirectories(File directory){
		if (directory.isDirectory() == false) return null;
		File[] fileNames = directory.listFiles();
		ArrayList al = new ArrayList();
		Pattern pat = Pattern.compile("^\\w+.*");
		Matcher mat; 
		for (int i=0; i< fileNames.length; i++)  {
			if (fileNames[i].isDirectory() == false) continue;
			mat = pat.matcher(fileNames[i].getName());
			if (mat.matches()) al.add(fileNames[i]);
		}
		//convert arraylist to file[]
		if (al.size() != 0){
			File[] files = new File[al.size()];
			al.toArray(files);
			Arrays.sort(files);
			return files;
		}
		else return new File[]{directory};
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
			ArrayList chromFiles = new ArrayList();
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


	/**Extracts all files from a comma delimited text of files and directories.
	 * No spaces. Returns null if file.canRead() is false. Not recursive.*/
	public static File[] extractFiles(String commaSeparList){
		if (commaSeparList == null) return null;
		String[] items = commaSeparList.split(",");
		File[] files = new File[items.length];
		for (int i=0; i<items.length; i++){
			files[i] = new File(items[i]);
			if (files[i].canRead()==false) return null;
		}
		Arrays.sort(files);
		return files;
	}

	/**Extracts all files with a given extension from a comma delimited text of files and directories.
	 * No spaces.*/
	public static File[] extractFiles(String commaSeparList, String extension){
		ArrayList filesAL = new ArrayList();
		String[] items = commaSeparList.split(",");
		for (int i=0; i<items.length; i++){
			File test = new File(items[i]);
			if (test.canRead()==false) return null;
			File[] files = extractFiles(test, extension);
			if (files == null) return null;
			for (int j=0; j<files.length; j++){
				filesAL.add(files[j]);
			}
		}
		File[] collection = new File[filesAL.size()];
		filesAL.toArray(collection);
		return collection;
	}

	/**Given full path file or directory names, extracts all files given an extension.*/
	public static File[] extractFiles(String[] dirFiles, String extension){
		int numDirectories = dirFiles.length;
		ArrayList celFilesAL = new ArrayList();
		String[] fileNames;
		try{
			for (int i=0; i<numDirectories; i++){
				File dir = new File(dirFiles[i]);
				if (dir.isDirectory()){
					String path = dir.getCanonicalPath();
					fileNames = dir.list();
					int num = fileNames.length;
					for (int j=0; j< num; i++)  {
						if (fileNames[j].endsWith(extension)) celFilesAL.add(new File(path,fileNames[j]));
					}
				}
				else if (dir.isFile() && dir.getName().endsWith(extension)) celFilesAL.add(dir);
			}

		}catch(IOException e){
			System.out.println("Problem extractFiles() ");
			e.printStackTrace();
		}
		File[] celFiles = new File[celFilesAL.size()];
		celFilesAL.toArray(celFiles);
		return celFiles;
	}

	/**Extracts the full path file names of all the files in a given directory with a given extension.*/
	public static File[] extractFilesReturnFiles(File directory, String extension){
		String[] files = extractFilesStringNames(directory, extension);
		int numFiles = files.length;
		File[] finalFiles = new File[numFiles];
		for (int i=0; i<numFiles; i++){
			finalFiles[i]= new File(files[i]);
		}
		return finalFiles;
	}

	/**Extracts the full path file names of all the files in a given directory with a given extension.*/
	public static String[] extractFilesStringNames(File directory, String extension){
		String[] files = null;			
		if (directory.isDirectory()){
			files = directory.list();
			int num = files.length;
			ArrayList chromFiles = new ArrayList();
			try{
				String path = directory.getCanonicalPath() + File.separator;
				for (int i=0; i< num; i++)  {
					if (files[i].endsWith(extension)) chromFiles.add(path+files[i]);
				}
				files = new String[chromFiles.size()];
				chromFiles.toArray(files);
			}catch(IOException e){
				System.out.println("Problem extractFilesStringNames() "+directory+" "+extension);
				e.printStackTrace();}
		}
		return files;
	}

	/**Extracts the full path file names of all the files not directories in a given directory. 
	 * Skips files starting with a non word character (e.g. '.', '_', etc).
	 * Returns null if no files found.*/
	public static File[] extractOnlyFiles(File directory){
		File[] fileNames = directory.listFiles();
		if (fileNames == null) return null;
		ArrayList<File> al = new ArrayList<File>();
		Pattern pat = Pattern.compile("^\\w+.*");
		Matcher mat; 
		for (int i=0; i< fileNames.length; i++)  {
			if (fileNames[i].isDirectory()) continue;
			mat = pat.matcher(fileNames[i].getName());
			if (mat.matches()) al.add(fileNames[i]);
		}
		//convert arraylist to file[]
		if (al.size() != 0){
			File[] files = new File[al.size()];
			al.toArray(files);
			Arrays.sort(files);
			return files;
		}
		else return null;
	}


	/**Fetches files that don't start with a '.' from a directory recursing through sub directories.*/
	public static ArrayList<File> fetchAllFilesRecursively (File directory){
		ArrayList<File> files = new ArrayList<File>(); 
		File[] list = directory.listFiles();
		if (list != null){
			for (int i=0; i< list.length; i++){
				if (list[i].isDirectory()) {
					ArrayList<File> al = fetchAllFilesRecursively (list[i]);
					int size = al.size();
					for (int x=0; x< size; x++){
						File test = al.get(x);
						if (test.getName().startsWith(".") == false) files.add(test);
					}
				}
				else if (list[i].getName().startsWith(".") == false) files.add(list[i]);				
			}
		}
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

	/**Fetches an ArrayList stored as a serialize object file.*/
	public static ArrayList fetchArrayList(File file) {
		return (ArrayList)fetchObject(file);
	}

	/**For each base text (ie /home/data/treat1 ), makes a File object using the chromExtenstion
	 * (ie /home/data/treat1.chromExtension ), if it actually exists, it is saved and a File[]
	 * returned.*/
	public static File[] fetchFiles(String[] baseNames, String chromExtension){
		ArrayList al = new ArrayList();
		//for each base text look for an associated file with the chromExtension
		for (int x=0; x<baseNames.length; x++){
			File f = new File(baseNames[x]+"."+chromExtension);
			if (f.exists())al.add(f);
		}
		File[] files = new File[al.size()];
		al.toArray(files);
		return files;
	}
	/**Fetches all files in a directory recursing through sub directories.*/
	public static File[] fetchFilesRecursively (File directory){
		ArrayList<File> al = fetchAllFilesRecursively (directory);
		File[] files = new File[al.size()];
		al.toArray(files);
		return files;
	}



	/**Loads float[][]s from disk.*/
	public static float[][][] fetchFloatFloatArrays(File[] files){
		int i=0;
		try {
			int num = files.length;
			float[][][] f = new float[num][][];
			for (;i<num; i++){
				f[i] = (float[][])IO.fetchObject(files[i]);
			}
			return f;
		} catch (Exception e){
			System.out.println("\nError: One of your 'xxx.cela' files does not appear" +
					" to be a serialized java float[][] array?! -> "+files[i]+"\n");
			System.exit(0);
		}
		return null;
	}

	/**Given a directory, returns a HashMap<String, File> of the containing directories' names and a File obj for the directory. */
	public static HashMap<String, File> fetchNamesAndDirectories(File directory){
		HashMap<String, File> nameFile = new HashMap<String, File>();
		File[] files = IO.extractFiles(directory);	
		for (int i=0; i< files.length; i++){
			if (files[i].isDirectory()) nameFile.put(files[i].getName(), files[i]);
		}
		return nameFile;
	}
	/**Given a directory, returns a HashMap<String, File> of the containing files names and a File obj for the directory. */
	public static HashMap<String, File> fetchNamesAndFiles(File directory){
		HashMap<String, File> nameFile = new HashMap<String, File>();
		File[] files = IO.extractFiles(directory);	
		for (int i=0; i< files.length; i++){
			if (files[i].isDirectory() == false) nameFile.put(files[i].getName(), files[i]);
		}
		return nameFile;
	}
	/**Fetches an Object stored as a serialized file.
	 * Can be zip/gz compressed too.*/
	public static Object fetchObject(File file) {
		Object a = null;
		try {
			ObjectInputStream in;
			if (file.getName().endsWith(".zip")){
				ZipFile zf = new ZipFile(file);
				ZipEntry ze = (ZipEntry) zf.entries().nextElement();
				in = new ObjectInputStream( zf.getInputStream(ze));
			}
			else if (file.getName().endsWith(".gz")) {
				in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(file)));
			}
			else in = new ObjectInputStream(new FileInputStream(file));
			a = in.readObject();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Problem fetchObject() "+file);
		}
		return a;
	}

	/**Fetches a BufferedReader from a url or file, zip/gz OK.*/
	public static BufferedReader fetchBufferedReader(String s){
		try {
			if (s.startsWith("http")) return fetchBufferedReader (new URL(s));
			return fetchBufferedReader (new File (s));
		} catch (Exception e) {
			System.out.println("Problem fetching buffered reader fro -> "+s);
			e.printStackTrace();
		} 
		return null;

	}

	/**Fetches a BufferedReader from a url, zip/gz OK.*/
	public static BufferedReader fetchBufferedReader(URL url) throws IOException{
		BufferedReader in = null;
		InputStream is = url.openStream();
		String name = url.toString();
		if (name.endsWith(".gz")) {
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(is)));
		}
		else if (name.endsWith(".zip")){
			ZipInputStream zis = new ZipInputStream(is);
			zis.getNextEntry();
			in = new BufferedReader(new InputStreamReader(zis));
		}
		else in = new BufferedReader(new InputStreamReader(is));
		return in;
	}

	/**Returns a gz zip or straight file reader on the file based on it's extension.
	 * @author davidnix*/
	public static BufferedReader fetchBufferedReader( File txtFile) throws IOException{
		BufferedReader in;
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
		return in;
	}

	/**Returns a gz zip or straight input strean on the file based on it's extension compression.*/
	public static InputStream fetchInputStream( File txtFile) throws IOException{
		InputStream in;
		String name = txtFile.getName().toLowerCase();
		if (name.endsWith(".zip")) {
			ZipFile zf = new ZipFile(txtFile);
			ZipEntry ze = (ZipEntry) zf.entries().nextElement();
			in = zf.getInputStream(ze);
		}
		else if (name.endsWith(".gz")) {
			in = new GZIPInputStream(new FileInputStream(txtFile));
		}
		else in = new FileInputStream(txtFile);
		return in;
	}

	/**Returns a BufferedReader from which you can call readLine() directly from a single entry zipped file without decompressing.
	 * Returns null if there is a problem.*/
	public static BufferedReader fetchReaderOnZippedFile (File zippedFile) {
		try {
			ZipFile zf = new ZipFile(zippedFile);
			ZipEntry ze = (ZipEntry) zf.entries().nextElement();
			return new BufferedReader(new InputStreamReader(zf.getInputStream(ze)));
		} catch(Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**Returns a BufferedReader from which you can call readLine() directly from a gzipped (.gz) file without decompressing.
	 * Returns null if there is a problem.*/
	public static BufferedReader fetchReaderOnGZippedFile (File gzFile) {
		try {
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gzFile))));
		} catch(Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**Fires qsub with the commands.
	 * @return String beginning with 'Error: ' if the error stream or java bombed. 
	 * otherwise it returns a tab delimited list of the output and the qsub parameters.
	 * Will leave the shell file in the temp directory, renamed as the jobNumber.sh
	 * @param qsubParams = '-q queue@loc -l nodes=1' */
	public static String fireQSub(String qsubParams, String commands, File tempDir){
		//make random word
		String name = "N"+Passwords.createRandowWord(9);
		//Create Results directory and files
		File resultsDirectory = new File (tempDir, "Results");
		if (resultsDirectory.exists() == false) resultsDirectory.mkdir();
		File launched = new File (resultsDirectory, "L_"+name);
		File results = new File (resultsDirectory, "E_"+ name);
		File returned = new File (resultsDirectory, "R_"+name);
		//create shell script wrapper file
		StringBuffer shell = new StringBuffer();
		shell.append ("echo running >> "); shell.append(launched); shell.append("\n");
		shell.append ("date > "); shell.append(results); shell.append("\n");
		shell.append (commands); shell.append(" >> "); shell.append(results);shell.append(" 2>> "); shell.append(results); shell.append("\n");
		shell.append ("date >> "); shell.append(results); shell.append("\n");
		shell.append ("rm "); shell.append(launched); shell.append("\n");
		shell.append ("mv "); shell.append(results); shell.append(" "); shell.append(returned); shell.append("\n");

		//QSub messages
		File qsubDir = new File (tempDir, "QSub");
		if (qsubDir.exists() == false) qsubDir.mkdir();
		String qsub = "qsub "+"-N "+name+" -j oe -o "+qsubDir+" "+qsubParams;
		//Save in shell directory
		File shellDirectory = new File(tempDir,"ShellScripts");
		if (shellDirectory.exists() == false) shellDirectory.mkdir();
		File shellFile = new File(shellDirectory,name +".sh");
		writeString(shell.toString(), shellFile);
		//write command line to launched
		writeString(qsub+" "+shellFile+"\n", launched); 
		try {
			//fire shell script
			Runtime rt = Runtime.getRuntime();
			//System.out.println("\nFiring: "+qsub+shellFile+ "\n");
			Process p = rt.exec(qsub +shellFile);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
			String line;

			//collect potential error
			ArrayList errorAL = new ArrayList();
			while ((line = error.readLine()) != null){
				errorAL.add(line);
			}
			if (errorAL.size()!=0){
				data.close();
				error.close();
				Misc.printExit("Error: "+Misc.stringArrayListToString(errorAL, "\n")  +"\n"+ commands +"\n"+ qsub+" "+shellFile);
			}
			//collect output of firing qsub
			ArrayList dataAL = new ArrayList();
			while ((line = data.readLine()) != null){
				dataAL.add(line);
			}
			//close handles
			data.close();
			error.close();
			//get job text
			String job = ((String)dataAL.get(0)).replaceAll("\\.cluster.+","");

			return name+"\t"+job;
		} catch (Exception e){
			return "Error: "+ e;
		} 
	}

	/**Fires qsub with the commands, will wait and return results and delete x.sh and the output files if you like.
	 * Be sure to put an stop of line '\n' after each command.*/
	public static String[] fireQSub(String queue, String commands, String fullPathToTempDir, boolean waitAndClean){
		//create shell script wrapper file
		String name = "tmp"+Passwords.createRandowWord(6);
		String shell = "#!/bin/sh \n" +
		"#$ -N "+name+" \n" +
		"#$ -j oe \n" +
		"#$ -o "+fullPathToTempDir+" \n"+
		"#$ -q "+queue+" \n"+
		"#$ -l nodes=1 \n"+
		commands+" \n" +
		"echo \"Finished!\"\n";
		System.out.print(shell);
		//System.exit(0);
		File shellFile = new File(fullPathToTempDir,name +".sh");
		writeString(shell, shellFile);

		//fire shell script
		try {
			Runtime rt = Runtime.getRuntime();
			Process p = rt.exec("qsub "+shellFile);
			System.out.println("Fired qsub");

			//check if user wants to wait for results and clean up files
			if (waitAndClean == false) return null;

			//Wait 5min then check for file and for Finished!
			boolean wait = true;
			File dir = new File(fullPathToTempDir);
			File outputFile = null;
			long counter =0;
			long milSec = 30000;  //checking every 30 sec
			ArrayList dataArrayList = new ArrayList(1000);

			System.out.print("Waiting");
			while (wait){
				Thread.sleep(milSec);
				//find file
				if (outputFile==null){				
					//fetch contents of directory
					String[] fileNames = dir.list();
					//run through each file
					for (int i= fileNames.length-1; i>=0; i--){
						if (fileNames[i].startsWith(name+".o")) { //the qsub output file MemeR1081900495429.o50151
							outputFile = new File(fullPathToTempDir, fileNames[i]);
							System.out.println("\nFound results file: "+outputFile.getName());
							break;
						}
					}
				}
				if (outputFile!=null){ //check if null again
					//check if the last line is "Finnished!"
					String lastLine = "";
					String line;
					dataArrayList.clear();
					BufferedReader in = new BufferedReader(new FileReader(outputFile));
					while ((line = in.readLine()) !=null) {
						lastLine= line;
						dataArrayList.add(line);
					}
					in.close();
					//System.out.println("\nLast line is : "+lastLine);			
					if (lastLine.startsWith("Finished!") || lastLine.startsWith("cd: Too many arguments")) wait=false;  //stop found exit wait loop

				}

				//put break in so thread will stop after a very long time.
				counter++;	
				System.out.print(".");

				if (counter>57600000) { //16hrs 57600000
					System.out.println("\n    Error: shell script failed to return from qsub after "+counter/60000+
							" minutes.\nFind and kill the job: (text Nix, process number is the last digits after the .o in -> "+outputFile.getName());
					System.exit(1);	
				}

			}


			System.out.println("\nResults are ready!");

			//delete files
			String[] fileNames = dir.list();
			//run through each file
			for (int i= fileNames.length-1; i>=0; i--){
				if (fileNames[i].startsWith(name)) { //the qsub output file MemeR1081900495429.o50151
					new File(fullPathToTempDir, fileNames[i]).delete();
					System.out.println("Deleting-> "+fileNames[i]);
				}
			}
			System.out.println("Time for run: "+(counter*(milSec/1000))+" seconds\n");

			//return program output
			dataArrayList.trimToSize();
			String[] x = new String[dataArrayList.size()];
			dataArrayList.toArray(x);
			return x;
		} catch (Exception e){
			System.out.println("Problem with fireQSub()");
			e.printStackTrace();
		}
		return null;
	}
	/**Gets full path text using getCanonicalPath() on a File*/
	public static String getFullPathName(File fileDirectory){
		String full = null;
		try {
			full = fileDirectory.getCanonicalPath();
		}catch (IOException e){
			System.out.println("Problem with getFullPathtName(), "+fileDirectory);
			e.printStackTrace();
		}
		return full;
	}

	/**Fetches and truncates common chars, from file names, front and back, uses getName() so no path info is involved.*/
	public static String[] getTruncatedNames(File[] files){
		String[] names  = new String[files.length];
		for (int i=0; i< names.length; i++){
			names[i] = files[i].getName();
		}
		names = Misc.trimCommon(names);
		return names;
	}

	/**Fetches and truncates common chars, from file names, front and back, uses getName() so no path info is involved.*/
	public static String[][] getTruncatedNames(File[][] files){
		//collapse to single array
		ArrayList names = new ArrayList();
		for (int i=0; i< files.length; i++){
			for (int j=0; j< files[i].length; j++){
				names.add(files[i][j].getName());
			}
		}
		//trim edges and break into original structure
		String[] trimmed = new String[names.size()];
		names.toArray(trimmed);
		trimmed = Misc.trimCommon(trimmed);
		String[][] finalTrimmed = new String[files.length][];
		int index = 0;
		for (int i=0; i< files.length; i++){
			finalTrimmed[i] = new String[files[i].length];
			for (int j=0; j< files[i].length; j++){
				finalTrimmed[i][j] = trimmed[index++];
			}
		}
		return finalTrimmed;
	}

	/**Loads a file's lines into a String[], will save blank lines. gz/zip OK*/
	public static String[] loadFile(File file){
		ArrayList<String> a = new ArrayList<String>();
		try{
			BufferedReader in = IO.fetchBufferedReader(file);
			String line;
			while ((line = in.readLine())!=null){
				line = line.trim();
				a.add(line);
			}
		}catch(Exception e){
			System.out.println("Prob loadFileInto String[]");
			e.printStackTrace();
		}
		String[] strings = new String[a.size()];
		a.toArray(strings);
		return strings;
	}

	/**Loads a file's lines into a hash first column is the key, second the value.
	 * */
	public static HashMap<String,String> loadFileIntoHashMap(File file){
		HashMap<String,String> names = new HashMap<String,String>(1000);
		try{
			BufferedReader in = IO.fetchBufferedReader(file);
			String line;
			String[] keyValue;
			while ((line = in.readLine())!=null){
				keyValue = line.split("\\s+");
				if (keyValue.length !=2 || keyValue[0].startsWith("#")) continue;
				names.put(keyValue[0].trim(), keyValue[1].trim());
			}
		}catch(Exception e){
			System.out.println("Prob loadFileInttoHash()");
			e.printStackTrace();
		}
		return names;
	}

	/**Loads a file's lines into a hash first column is the key, second the value.
	 * */
	public static HashMap<String,Integer> loadFileIntoHashMapStringInteger(File file){
		HashMap<String,Integer> names = new HashMap<String,Integer>(1000);
		try{
			BufferedReader in = IO.fetchBufferedReader(file);
			String line;
			String[] keyValue;
			while ((line = in.readLine())!=null){
				keyValue = line.split("\\s+");
				if (keyValue.length !=2 || keyValue[0].startsWith("#")) continue;
				names.put(keyValue[0].trim(), new Integer(keyValue[1].trim()));
			}
		}catch(Exception e){
			System.out.println("Prob loadFileInttoHash()");
			e.printStackTrace();
		}
		return names;
	}

	/**Loads a file's lines into a hash set, keys only.*/
	public static HashSet<String> loadFileIntoHashSet(File file){
		HashSet<String> names = new HashSet(10000);
		try{
			BufferedReader in = fetchBufferedReader(file);
			String line;
			while ((line = in.readLine())!=null){
				line = line.trim();
				if (line.length() ==0) continue;
				//if (names.contains(line)) {
				//	System.out.println("\tDuplicate line found while loading hash -> "+line);
				//}
				names.add(line);
			}
		}catch(Exception e){
			System.out.println("Prob loadFileInttoHash()");
			e.printStackTrace();
		}
		return names;
	}



	/**Loads a file's lines into a hash.*/
	public static HashSet loadFileIntoHashSplitLines(File file){
		HashSet names = new HashSet(10000);
		try{
			BufferedReader in = fetchBufferedReader(file);
			String line;
			String[] tokens;
			while ((line = in.readLine())!=null){
				line = line.trim();
				tokens = line.split("\\s+");
				if (tokens.length ==0) continue;
				for (int i=0; i< tokens.length; i++) names.add(tokens[i]);
			}
		}catch(Exception e){
			System.out.println("Prob loadFileInttoHashSplitLines()");
			e.printStackTrace();
		}
		return names;
	}

	/**Loads a file's lines into a linked hash set, keys only.*/
	public static LinkedHashSet loadFileIntoLinkedHashSet(File file){
		LinkedHashSet names = new LinkedHashSet(10000);
		try{
			BufferedReader in = fetchBufferedReader(file);
			String line;
			while ((line = in.readLine())!=null){
				line = line.trim();
				if (line.length() ==0) continue;
				names.add(line);
			}
		}catch(Exception e){
			System.out.println("Prob loadFileIntoLinkedHashSet()");
			e.printStackTrace();
		}
		return names;
	}



	/**Loads a file's lines into a String[] trimming blank lines and beginning and ending whitespace.*/
	public static String[] loadFileIntoStringArray(File file){
		ArrayList a = new ArrayList();
		try{
			BufferedReader in = fetchBufferedReader(file);
			String line;
			while ((line = in.readLine())!=null){
				line = line.trim();
				if (line.length() ==0) continue;
				a.add(line);
			}
		}catch(Exception e){
			System.out.println("Prob loadFileInto String[]");
			e.printStackTrace();
		}
		String[] strings = new String[a.size()];
		a.toArray(strings);
		return strings;
	}
	/**Loads a tab delimited file into a LinkedHashMap where each line is broken into cells and loaded into
	 * a LinkedHashMap where the cells are the keys and the value is the line.*/
	public static LinkedHashMap loadHash (File tabDelimited){
		LinkedHashMap hash = new LinkedHashMap();
		try{
			BufferedReader in = fetchBufferedReader(tabDelimited);
			String line;
			String[] cells;
			Pattern pat = Pattern.compile("\\s");
			//run through each line
			while ((line = in.readLine()) != null){
				String trimmedLine = new String (line.trim());
				if (trimmedLine.length() == 0 || trimmedLine.startsWith("#")) continue;
				cells = pat.split(trimmedLine);
				//only add cells that are unique to big hash
				HashSet unique = new HashSet();
				for (int i=0; i< cells.length; i++){
					//is it blank?
					cells[i] = cells[i].trim();
					if (cells[i].length() == 0) continue;
					//is it unique?
					if (unique.contains(cells[i]) == false){
						//add to hash
						unique.add(cells[i]);
						//add to big hash only if it isn't present, throw warning if present
						if (hash.containsKey(cells[i])) {
							System.err.println("Warning: duplicate key found!\n\tLine -> "+line +"\n\tKey -> "+cells[i]);
						}
						else hash.put(cells[i], trimmedLine);
					}
				}
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return hash;
	}

	/**Loads a key = value containing text file into a LinkedHashMap.
	 * Strips all white space, ignores lines beginning with / or #.
	 * Ignores blank lines.*/
	public static LinkedHashMap loadKeyValueFile(File file){
		LinkedHashMap map = new LinkedHashMap();
		try{
			BufferedReader in = fetchBufferedReader(file);
			String line;
			Pattern comment = Pattern.compile("^[/|#].*");
			Pattern whiteSpace = Pattern.compile("\\s+");
			Pattern equal = Pattern.compile("=");
			while ((line = in.readLine())!=null){
				line = line.trim();
				//ignore blank lines and comment lines '// or #'
				if (line.length() == 0 || comment.matcher(line).matches()) continue;
				//strip white space
				line = whiteSpace.matcher(line).replaceAll("");
				//split and put
				String[] keyValue = equal.split(line);
				if (keyValue.length !=2) continue;
				map.put(keyValue[0], keyValue[1]);
			}
		}catch(Exception e){
			System.out.println("Prob with loadKeyValueFile()");
			e.printStackTrace();
		}
		return map;
	}

	/**Loads a key value containing text file into a LinkedHashMap.
	 * Ignores lines beginning with / or #.
	 * Ignores blank lines.
	 * Splits on spaces/ tabs thus no spaces within key or value.
	 * Returns null if split keyValue number is != 2.*/
	public static LinkedHashMap loadKeyValueSpaceTabDelimited(File file, boolean printErrorMessages){
		LinkedHashMap map = new LinkedHashMap();
		try{
			BufferedReader in = fetchBufferedReader(file);
			String line;
			Pattern comment = Pattern.compile("^[/|#].*");
			Pattern whiteSpace = Pattern.compile("\\s+");
			while ((line = in.readLine())!=null){
				line = line.trim();
				//ignore blank lines and comment lines '// or #'
				if (line.length() == 0 || comment.matcher(line).matches()) continue;
				//split and put
				String[] keyValue = whiteSpace.split(line);
				if (keyValue.length !=2) {
					if (printErrorMessages) System.err.println("Cannot parse line -> "+line);
					return null;
				}
				map.put(keyValue[0], keyValue[1]);
			}
		}catch(Exception e){
			if (printErrorMessages) e.printStackTrace();
		}
		return map;
	}

	/**Loads a tab delimited file containing ints into an int[columns][rows], all lines included, assumes no missing values.*/
	public static int[][] loadTableOfInts (File file){
		int[][] columnsRows = null;
		try{
			BufferedReader in = fetchBufferedReader(file);
			//read first line 
			String line=in.readLine();
			if (line == null) return null;
			String[] tokens = line.split("\\t");
			//find number columns
			int numColumns = tokens.length;
			//find number of rows
			int numRows = (int)IO.countNumberOfLines(file);
			columnsRows = new int[numColumns][numRows];
			//assign first row to table
			int y=0;
			for (int x =0; x< numColumns; x++) columnsRows[x][y] = Integer.parseInt(tokens[x]);
			y++;
			//load other rows to table
			for (; y< numRows; y++){
				line=in.readLine();
				tokens = line.split("\\t");
				for (int x =0; x< numColumns; x++) columnsRows[x][y] = Integer.parseInt(tokens[x]);
			}
			in.close();
		} catch (Exception e){
			System.out.println("Problem loading table");
			e.printStackTrace();
		}
		return columnsRows;
	}
	
	/**Loads a tab delimited file containing double into an double[columns][rows], all lines included, assumes no missing values. No headers please.*/
	public static double[][] loadTableOfDoubles (File file){
		double[][] columnsRows = null;
		try{
			BufferedReader in = fetchBufferedReader(file);
			//read first line 
			String line=in.readLine();
			if (line == null) return null;
			String[] tokens = line.split("\\t");
			//find number columns
			int numColumns = tokens.length;
			//find number of rows
			int numRows = (int)IO.countNumberOfLines(file);
			columnsRows = new double[numColumns][numRows];
			//assign first row to table
			int y=0;
			for (int x =0; x< numColumns; x++) columnsRows[x][y] = Double.parseDouble(tokens[x]);
			y++;
			//load other rows to table
			for (; y< numRows; y++){
				line=in.readLine();
				tokens = line.split("\\t");
				for (int x =0; x< numColumns; x++) columnsRows[x][y] = Double.parseDouble(tokens[x]);
			}
			in.close();
		} catch (Exception e){
			System.out.println("Problem loading table");
			e.printStackTrace();
		}
		return columnsRows;
	}

	/**Attempts to make a file and return it's canonical path text (full path text).*/
	public static String makeFullPathName(String resultsDirectory, String fileName){
		File dump = new File(resultsDirectory, fileName);
		String fullPathName = "";
		try{
			fullPathName = dump.getCanonicalPath();
		}
		catch (IOException e){
			System.out.println("Problem building text to write files!");
			e.printStackTrace();
		}
		return fullPathName;
	}

	/**Returns a HashMap of the file text sans extension and the file.
	 * Returns null if duplicate text is found.*/
	public static HashMap makeNameFileHash(File[] files){
		HashMap hash = new HashMap();
		for (int i=0; i< files.length; i++){
			String name = Misc.removeExtension(files[i].getName());
			if (hash.containsKey(name)) return null;
			hash.put(name, files[i]);
		}
		return hash;
	}

	/**Returns the number of files that stop with the given extension.*/
	public static int numberFilesExist(File directory, String extension){
		String[] fileNames = directory.list();
		int num = 0;
		for (int i=0; i< fileNames.length; i++){
			if (fileNames[i].endsWith(extension)) num++;
		}
		return num;
	}

	/**Returns the number of files that exist in the array.*/
	public static int numberFilesExist(File[] f){
		int num =0;
		for (int i=0; i< f.length; i++) if (f[i].exists()) num++;
		return num;
	}

	/**Parses the indexed column of a tab delimited file, all lines included.
	 * Will return an empty String if index column doesn't exist*/
	public static String[] parseColumn (File file, int index){
		ArrayList al = new ArrayList();
		try{
			BufferedReader in = new BufferedReader (new FileReader(file));
			String line;
			String[] tokens;
			while ((line=in.readLine()) != null){
				tokens = line.split("\\t");
				if (tokens.length <= index) al.add("");
				else al.add(tokens[index].trim());
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printExit("\nError: problem parsing matcher file, aborting.\n");
		}
		String[] col = new String[al.size()];
		al.toArray(col);
		return col;
	}

	/**Parses a tab delimited file, the indexed column is used as the key, 
	 * the entire line as the value, blank lines skipped, returns null if 
	 * a duplicate key is found.*/
	public static HashMap parseFile (File file, int index, boolean ignoreDuplicateKeys){
		HashMap al = new HashMap();
		String line = null;
		try{
			BufferedReader in = fetchBufferedReader(file);

			String[] tokens;
			while ((line=in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				tokens = line.split("\\t");
				String item = tokens[index].trim();
				if (ignoreDuplicateKeys) al.put(item, line);
				else if (al.containsKey(item)) throw new Exception("Duplicate found for -> "+tokens[index]);
				else al.put(item, line);
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("\nError: problem parsing file to hash. Line -> "+line);
			return null;
		}
		return al;
	}

	/**Parses a tab delimited file, the indexed column is used as the key, 
	 * the entire line as the value, blank lines skipped, removes duplicate keys.*/
	public static HashMap parseUniqueKeyFile (File file, int index){
		HashMap al = new HashMap();
		String line = null;
		try{
			BufferedReader in = fetchBufferedReader(file);
			HashSet<String> badKeys = new HashSet<String>();
			String[] tokens;
			while ((line=in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0) continue;
				tokens = line.split("\\t");
				String item = tokens[index].trim();
				if(al.containsKey(item)) badKeys.add(item);
				else al.put(item, line);
			}
			//remove bad keys
			for (String baddie: badKeys) al.remove(baddie);
			if (badKeys.size()!=0) System.err.println("Dropped duplicate keys -> "+badKeys);
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("\nError: problem parsing file to hash.");
			return null;
		}
		return al;
	}


	/**Removes the extension from each File in the directory.*/
	public static void removeExtension(File directory, String extension){		
		File[] files = directory.listFiles();
		for (int i=0; i < files.length; i++){
			String name = files[i].getName();
			if (name.endsWith(extension)) {
				int index = name.lastIndexOf(extension);
				name = name.substring(0,index);
				File renamed = new File (directory, name);
				files[i].renameTo(renamed);
			}
		}
	} 

	/**Renames files containing the toReplace with replaceWith using String.replaceAll(String, String)*/
	public static void renameFiles(File directory, String toReplace, String replaceWith){
		File[] files = directory.listFiles();
		for (int i=0; i< files.length; i++){
			String name = files[i].getName();
			String changed = name.replaceAll(toReplace,replaceWith);
			if (name.equals(changed) == false) files[i].renameTo(new File (files[i].getParent()+File.separator+changed));
		}
	}

	public static boolean saveObject(File file, Object ob) {
		try {
			ObjectOutputStream out =
				new ObjectOutputStream(new FileOutputStream(file));
			out.writeObject(ob);
			out.close();
			return true;
		} catch (Exception e) {
			System.out.println("There appears to be a problem with saving this file: "+ file);
			e.printStackTrace();
		}
		return false;
	}

	/**Writes each element of an Arraylist on a seperate line to a file.*/
	public static boolean writeArrayList(ArrayList al, File file){
		try{
			PrintWriter out  = new PrintWriter(new FileWriter(file));
			int num = al.size();
			for (int i=0; i<num; i++){
				out.println(al.get(i));
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;	
	}

	/**Writes each key tab value on a separate line to the file.*/
	public static boolean writeHashMap(HashMap hashMap, File file){
		try{
			PrintWriter out  = new PrintWriter(new FileWriter(file));
			Iterator it = hashMap.keySet().iterator();
			while (it.hasNext()){
				Object obj = it.next();
				out.println(obj+"\t"+hashMap.get(obj));
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;	
	}
	
	/**Writes each key tab value on a separate line to the file.*/
	public static boolean writeHashSet(HashSet set, File file){
		try{
			PrintWriter out  = new PrintWriter(new FileWriter(file));
			Iterator it = set.iterator();
			while (it.hasNext()){
				out.println(it.next());
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;	
	}

	/**Checks to see if the browser sending the request is from a windows machine.
	 * If so it converts the file text dividers to unix.
	 * In either case, replaces any spaces with underscores.*/
	public static String winFileConvert(String fileName, HttpServletRequest request){
		String name = fileName;
		if (request.getHeader("user-agent").indexOf("Win")!=-1){
			int index = fileName.lastIndexOf("\\");
			name = fileName.substring(index+1);
		}
		return name.replaceAll(" ","_");
	}

	/**Writes a String to disk. */
	public static boolean writeString(String data, File file) {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(file));
			out.print(data);
			out.close();
			return true;
		} catch (IOException e) {
			System.out.println("Problem writing String to disk!");
			e.printStackTrace();
			return false;
		}
	}

	/**Writes a String to disk.*/
	public static boolean writeString(String data, String fullPathFileName) {
		return writeString(data, new File(fullPathFileName));
	}	
	/**Writes a String onto the stop of a file.
	 * No need to add a \n onto the final stop. */
	public static boolean writeStringAppend(String data, File file) {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(file, true));
			out.print(data);
			out.close();
			return true;
		} catch (IOException e) {
			System.out.println("Problem writing String to disk!");
			e.printStackTrace();
			return false;
		}
	}

	/**Zip compress a single file.*/
	public static boolean zip(File fileToZip, File zipFile){
		File[] file = new File[1];
		file[0]= fileToZip;
		return zip(file, zipFile);
	}

	/**Zip compresses an array of Files, be sure to text your zipFile with a .zip extension!*/
	public static boolean zip(File[] filesToZip, File zipFile ){
		byte[] buf = new byte[2048];	
		try {
			ZipOutputStream out = new ZipOutputStream(new FileOutputStream(zipFile));		
			// Compress the files
			for (int i=0; i<filesToZip.length; i++) {
				FileInputStream in = new FileInputStream(filesToZip[i]);
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
			return false;
		}
		return true;
	}

	/**Zip compresses the file and if sucessful, then deletes the original file.*/
	public static boolean zipAndDelete(File f){
		File zip = null;
		try {
			//make new file with .zip extension
			zip = new File(f.getCanonicalPath()+".zip");
			//zip it and delete original
			if (IO.zip(f, zip)) {
				f.delete();
				return true;
			}
			else {
				zip.delete();
				return false;
			}
		} catch (Exception e){
			e.printStackTrace();
			zip.delete();
			return false;
		}
	}

	/**Returns the contents of a zipped file. Returns null if there is a problem.*/
	public static String[] zippedFile2StringArray(File zippedFile){
		try {
			BufferedReader in = fetchReaderOnZippedFile(zippedFile);
			ArrayList al = new ArrayList();
			String line = null;
			while ((line= in.readLine()) !=null) al.add(line);
			return Misc.stringArrayListToStringArray(al);
		} catch(Exception e) {
			e.printStackTrace();
		}
		return null;
	}

}
