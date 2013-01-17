package edu.utah.ames.bioinfo;
import java.io.File;
import java.io.IOException;

//Finds cmd.txt file in specified directory, using FindFile

public class FindFileTest {
	
    public static void main(String[] args) throws IOException {
        new FindFileTest();
    }
    
    //define file search location
    public FindFileTest() throws IOException {
        File dir = new File("/Users/darren/Desktop/testDir/");
        
        //list files found using FileFilter
        File[] files = dir.listFiles(new FindFile());
        
        //print resulting found files
        for (File f : files) {
            System.out.println("file: " + f.getCanonicalPath());
        }
    }
    
}
