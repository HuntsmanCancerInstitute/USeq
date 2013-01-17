package edu.utah.ames.bioinfo;
import java.io.File;

public class FolderTraversarTest {

	public static void main(String[] args) {
		String folderPath = "/Users/darren/Desktop/testDir/";
		FileFind traversal = new FileFind(new File(folderPath), null);
		traversal.find();
	}
	
}
