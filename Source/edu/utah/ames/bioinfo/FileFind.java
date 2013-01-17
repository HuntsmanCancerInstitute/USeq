package edu.utah.ames.bioinfo;
import java.io.File;

/**
 * Recursively finds cmd.txt file in user directories
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class FileFind {

	private File fileObject;
	private User user;
 
	public static void main(String[] args) {
		String path = "/Users/darren/Desktop/testDir/";
		User user = new User();
		FileFind ff = new FileFind(new File(path), user);
		ff.find();
	}
	
	//constructor
	public FileFind(File fileObject, User user) {
		this.fileObject = fileObject;
		this.user = user;
	}
 
	public void find() {
		recursiveFind(fileObject);
	}
	
	//recursively search dir
	public void recursiveFind(File fileObject) {
		//System.out.println(fileObject + " RT");
		if (fileObject.isDirectory()) {
			//System.out.println(indent + fileObject.getName());
			File allFiles[] = fileObject.listFiles();
			for (File aFile : allFiles) {
				recursiveFind(aFile);
			}
		} 
		
		if (fileObject.getName().equalsIgnoreCase("cmd.txt")) {
			user.getCmdFiles().add(fileObject);
		}
	}
}