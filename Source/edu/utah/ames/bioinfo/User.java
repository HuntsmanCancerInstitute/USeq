package edu.utah.ames.bioinfo;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

/*
 * This class holds information about users

/**
 *
 * @author darren.ames@hci.utah.edu
 */
public class User {
    
    //fields
    private String name;
    private File directory;
    private ArrayList<File> cmdFiles = new ArrayList<File>();
    private HashSet<String> emails = new HashSet<String>();
    //private ArrayList<String> bigOldFiles = new ArrayList<String>();
    private ArrayList<File> warnFiles = new ArrayList<File>();
    private ArrayList<File> deleteFiles = new ArrayList<File>();
    
    //no-arg constructor
    public User() {
    	name = null;
    	directory = null;
    	cmdFiles = null;
    	emails = null;
    	//bigOldFiles = null;
    	warnFiles = null;
    	deleteFiles = null;
    }
    
    public User(File dir) {
        directory = dir;
        name = dir.getName();
    }
    
    //constructor accepts a String arg assigned to the name field
    public User(String n) {
        name = n;
    }
    
    //constructor accepts a String arg assigned to the name field
    //and a String arg assigned to the emailAddress field
    public User(String n, HashSet<String> e) {
        name = n;
        emails = e;
    }
    
    public boolean isGood() {
        boolean rollWithIt = false;
        
        //instantiate EmailAddress
        String e = null;
        EmailAddress email = new EmailAddress(e);
        FileInfo f = new FileInfo();
        if (email.isValidEmailAddress(e) && f.fileExists() == true)
            rollWithIt = true;
            return rollWithIt;
    }
    
    //assigns its arg to the name field
    public void setName(String n) {
        name = n;
    }
    
    //returns the String in the name field
    public String getName() {
        return name;
    }
    
    public void setEmail(HashSet<String> em) {
    	emails = em;
    }
    
    public HashSet<String> getEmail() {
    	return emails;
    }

	public ArrayList<File> getCmdFiles() {
		return cmdFiles;
	}

	public void setCmdFiles(ArrayList<File> cmdFiles) {
		this.cmdFiles = cmdFiles;
	}

	public File getDirectory() {
		return directory;
	}

	public void setDirectory(File directory) {
		this.directory = directory;
	}

	public ArrayList<File> getWarnFiles() {
		return warnFiles;
	}

	public void setWarnFiles(ArrayList<File> warnFiles) {
		this.warnFiles = warnFiles;
	}

	public ArrayList<File> getDeleteFiles() {
		return deleteFiles;
	}

	public void setDeleteFiles(ArrayList<File> deleteFiles) {
		this.deleteFiles = deleteFiles;
	}
	
}
