package edu.utah.ames.bioinfo;
import java.io.File;

/**
 * simple class with methods that report several bits of info about a file
 */
/**
 *
 * @author darren.ames@hci.utah.edu
 */
public class FileInfo {

    //fields
    private boolean warn;
    private boolean delete;
    private int ageDays;
    private int fileSize;
    private File cmdFile;
    
    //no-arg constructor
    public FileInfo() {
    	cmdFile = null;
        warn = false;
        delete = false;
        ageDays = 0;
        fileSize = 0;
    }
    
    //regular constructor
    public FileInfo(File cmdFile, boolean w, boolean d, int a, int s) {
        this.cmdFile = cmdFile;
        warn = w;
        delete = d;
        ageDays = a;
        fileSize = s;
    }
    
    //qc method to test if file exists
    public boolean fileExists() {
        String path = null;
        File file = new File(path);
        boolean exists = false;
        if (file.exists() == true)
            exists = true;
        return exists;
    }
    
    //assigns its arg to the warn field
    public void setWarn(boolean w) {
        warn = w;
    }
    
    //assigns its arg to the delete field
    public void setDelete(boolean d) {
        delete = d;
    }
    
    //assigns its arg to the ageDays field
    public void setAgeDays(int a) {
        ageDays = a;
    }
    
    //assigns its arg to the fileSize field
    public void setFileSize(int s) {
        fileSize = s;
    }
    
    
    //returns the boolean in the warn field
    public boolean getWarn() {
        return warn;
    }
    
    //returns the boolean in the delete field
    public boolean getDelete() {
        return delete;
    }
    
    //returns the int value in the ageDays field
    public int getAgeDays() {
        return ageDays;
    }
    
    //returns the int value in the fileSize field
    public int getFileSize() {
        return fileSize;
    }
    
    public File getCmdFile() {
    	return cmdFile;
    }
    
    public void setCmdFile(File cmdFile) {
    	this.cmdFile = cmdFile;
    }
}
