package edu.utah.ames.bioinfo;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * This class is used to perform operations with directories and files present in the file system.
 */
public class Directory {

	private String filePath;
	private File file = null;

	public String getFilePath() {
		return filePath;
	}

	public void setFilePath(String filePath) {
		this.filePath = filePath;
	}
	
	//constructor
	public Directory(String path) {
		this.filePath = path;
		file = new File(filePath);
	}

	/**
	 * Returns a File Array containing list of sub directories in the given directory.
	 * 
	 * @return Sub-Directories List
	 * @throws Exception
	 */
	public File[] getSubDirectories() throws Exception {
		List<File> fileLst = new ArrayList<File>();

		if (file.isDirectory()) {
			for (File f : file.listFiles()) {
				if (f.isDirectory()) {
					fileLst.add(f);
					continue;
				}
			}
			File[] array = (File[]) fileLst.toArray(new File[fileLst.size()]);
			return array;
		}
		else {
			throw new Exception("The given path is not a Directory");
		}
	}

	/**
	 * Returns a File Array containing list of sub directories in the given directory. The returned array will be in
	 * sorted order based on the boolean arguments.
	 * 
	 * @param sort
	 * @param sortDescending
	 * @return Sub-Directories List
	 * @throws Exception
	 */
	public File[] getSubDirectories(boolean sort, boolean sortDescending) throws Exception {

		File[] flArr = getSubDirectories();
		if (sort == true) {
			if (sortDescending == true)
			{
				Arrays.sort(flArr, new FileComp());
			}
			else {
				Arrays.sort(flArr);
			}
		}
		return flArr;
	}

	/**
	 * Returns a File Array containing list of files in the directory
	 * 
	 * @return
	 * @throws Exception
	 */
	public File[] getFiles() throws Exception {
		List<File> fileLst = new ArrayList<File>();

		if (file.isDirectory()) {
			for (File f : file.listFiles()) {
				if (!f.isDirectory()) {
					fileLst.add(f);
					continue;
				}
			}
			File[] array = (File[]) fileLst.toArray(new File[fileLst.size()]);
			return array;
		}
		else {
			throw new Exception("The given path is not a Directory");
		}
	}

	/**
	 * Returns a File Array in sorted order based on the boolean argument.
	 * 
	 * @param sort
	 * @param sortDescending
	 * @return File List
	 * @throws Exception
	 */
	public File[] getFiles(boolean sort, boolean sortDescending) throws Exception {
		File[] fileArr = getFiles();

		if (sort == true) {
			if (sortDescending == true) {
				Arrays.sort(fileArr, new FileComp());
			}
			else {
				Arrays.sort(fileArr);
			}
		}
		return fileArr;
	}

	/**
	 * Returns a File object for the given directory path.
	 * 
	 * @return File
	 * @throws Exception
	 */
	public File getFile() throws Exception {

		if (file.isDirectory()) {
			throw new Exception("The given path is not a Directory");
		}
		else {
			return file;
		}
	}

	/**
	 * Returns a list of file extensions in the given directory
	 * 
	 * @return List of File extensions
	 * @throws Exception
	 */
	public List<String> getFileExtensionList() throws Exception {
		File[] fileArr = getFiles();
		HashSet<String> s = new HashSet<String>();
		for (File f : fileArr) {
			String fileName = f.getName();
			s.add(fileName.substring(fileName.lastIndexOf('.') + 1, fileName.length()));
		}
		return new ArrayList<String>(s);
	}

	/**
	 * Returns a boolean value depending on the availability of the given extension as argument in the directory.
	 * 
	 * @param extn
	 * @return
	 * @throws Exception
	 */
	public boolean isExtensionAvailable(String extn) throws Exception {
		List<String> extnLst = getFileExtensionList();
		return extnLst.contains(extn);
	}

	/**
	 * Returns a list of files satisfying the given regular expressions
	 * 
	 * @param regExp
	 * @return
	 * @throws Exception
	 */
	public void getFiles(Sample sample, String regExp) throws Exception {
		File[] fileArr = getFiles();

		for (File f : fileArr) {
			String fileName = f.getName();
			Matcher m = Pattern.compile(regExp).matcher(fileName);
			if (m.find()) {
				sample.setNovoindex(fileName);
			}
		}
	}

	/**
	 * Displays all files and directories in the given directory path.
	 * 
	 * @return
	 * @throws Exception
	 */
	public void displayFilesAndDircectories() throws Exception {

		for (File f : file.listFiles()) {
			if (f.isDirectory()) {
				System.out.println("   Directory -> " + f.getName());
				continue;
			}
			System.out.println("   File -> " + f.getName());
		}
	}

	public static void main(String[] args) throws Exception {
		// set the directory path
		Directory f = new Directory("C:\\Project\\java");
		System.out.println("Display all Files and Directories in the Given Path -> ");
		f.displayFilesAndDircectories();
		System.out.println("Get Files in the Given Directory -> ");
		for (File f1 : f.getFiles()) {
			System.out.println("   " + f1.getName());
		}
		System.out.println("Get Files in the Given Directory in Descending Order-> ");
		for (File f1 : f.getFiles(true, true)) {
			System.out.println("   " + f1.getName());
		}
		System.out.println("Get SubDirectories in the Given Directory in Order-> ");
		for (File f1 : f.getSubDirectories()) {
			System.out.println("   " + f1.getName());
		}
		System.out.println("Get SubDirectories in the Given Directory in Descending Order-> ");
		for (File f1 : f.getSubDirectories(true, true)) {
			System.out.println("   " + f1.getName());
		}
		System.out.println("Get files in the Given Directory satisfying the regular expression '^h.' ");
		//for (File f1 : Arrays.asList(f.getFiles("^h.", filePath))) {
			//System.out.println("   " + f1.getName());
		//}
		//System.out.println("List of the file extensions in the given directory ");
		for (String f1 : f.getFileExtensionList()) {
			System.out.println("   " + f1);
		}
		System.out.println("Does the given directory contains jar files -> " + f.isExtensionAvailable("jar"));
	}
}

class FileComp implements Comparator<File> {

	@Override
	public int compare(File f1, File f2) {
		String a1 = f1.getName().toLowerCase();
		String b1 = f2.getName().toLowerCase();
		return b1.compareTo(a1);
	}
}