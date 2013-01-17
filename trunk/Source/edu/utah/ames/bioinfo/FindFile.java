package edu.utah.ames.bioinfo;
import java.io.File;
import java.io.FileFilter;

/**
 *
 * @author darren.ames@hci.utah.edu
 */

/*
 * A class implementing the Java FileFilter interface
 */

public class FindFile implements FileFilter {
    
    private final String[] fileName = new String[] {"cmd.txt"};
    
    @Override
    public boolean accept(File file) {
        
        for (String name : fileName) {
            if (file.getName().toLowerCase().equals(name)) {
                return true;
            }
        }
        return false;
    }
    
}
