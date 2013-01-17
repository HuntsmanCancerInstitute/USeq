package edu.utah.ames.bioinfo;
import java.io.File;
import java.io.IOException;

//Finds cmd.txt file in specified directory, using FindFile

/**
 *
 * @author darren.ames@hci.utah.edu
 */
public class FindFileTest2 {
	
	
	FileInfo fileInfo = null;
	File dir = null;
    
    //define file search location
    public FindFileTest2() throws IOException {
        
    	this.dir = dir;
    	
        //create new FileInfo object to hold found file name
        fileInfo = new FileInfo();
       
        //list files found using FileFilter in FindFile
        File[] files = dir.listFiles(new FindFile());
         
        //print resulting found files
        for (File f : files) {
            //System.out.println("file: " + f.getCanonicalPath());
            if (f.getName().equals("cmd.txt")) fileInfo.setCmdFile(f);
            //System.out.println(fileInfo.getCmdFile());
        }
    }


	public FileInfo getFileInfo() {
		return fileInfo;
	}


	public void setFileInfo(FileInfo fileInfo) {
		this.fileInfo = fileInfo;
	}
}
