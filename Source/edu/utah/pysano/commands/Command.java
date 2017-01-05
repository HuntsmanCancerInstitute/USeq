package edu.utah.pysano.commands;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import edu.utah.pysano.utils.Logger;
import edu.utah.pysano.daemon.CThread;


public class Command {
	//File paths
	private File jobDirectory;
	private File finalDirectory;
	private File projectDirectory;
	private ArrayList<File> finalDirectoryList = new ArrayList<File>();

	private ArrayList<File> inputFiles = new ArrayList<File>();
	private ArrayList<File> externalFiles = new ArrayList<File>();
	private ArrayList<File> outputFiles = new ArrayList<File>();
	private ArrayList<File> protectList = new ArrayList<File>();
	private ArrayList<File> deleteList = new ArrayList<File>();
	private ArrayList<File> skipList = new ArrayList<File>(); //List of files that no longer need to be checked for output
	private HashMap<File,ArrayList<File>> conditionalSkipList = new HashMap<File,ArrayList<File>>(); //if key file is found, no longer check for the files contained in the arrayList
	
	private ArrayList<String> templateStringList = new ArrayList<String>();
	
	private HashMap<String,String> globalProperties = new HashMap<String,String>();
	private ArrayList<HashMap<String,String>> localProperties = new ArrayList<HashMap<String,String>>();
	
	//Command options
	private String cluster;
	private String email;
	private Integer wallTime;
	private boolean suppress;

	//SampleInfo
	private String sampleName = "NA";
	private String sampleID = "NA";
	
	
	//Logger
	Logger logFile = null;
	
	public Command(File projectDirectory, String cluster, String email, Integer wallTime, boolean suppress, Logger logFile) {
		this.cluster = cluster;
		this.email = email;
		this.wallTime = wallTime;
		this.suppress = suppress;
		this.logFile = logFile;
		this.projectDirectory = projectDirectory;
	}

	public void postProcess() {
		
		//Create final directories
		for (File d: finalDirectoryList) {
			if (!d.exists()) {
				d.mkdirs();
			}
		}
		
		//Move output files, skipping things that should no longer exist
		for (File outputFile: outputFiles) {
			if (skipList.contains(outputFile)) {
				continue;
			}
			String directory = outputFile.getParent();
			String basename = outputFile.getName();
			File workingFile = new File(jobDirectory,basename);
			File finalFile = new File(directory,basename);
			if (workingFile.isDirectory()) {
				finalFile.mkdirs();
				logFile.writeInfoMessage("[Command] copying contents of " + workingFile.getAbsolutePath() + " to " + finalFile.getAbsolutePath());
				for (File f: workingFile.listFiles()) {
					cpFile(f,finalFile);
				}
			} else {
				logFile.writeInfoMessage("[Command] moving " + workingFile.getAbsolutePath() + " to " + finalFile.getAbsolutePath());
				moveFile(workingFile,finalFile, logFile);	
			}	
		}
		
		//Create Job log directories.
		File finalJobDirectory = new File(projectDirectory,"Jobs");
		if (!finalJobDirectory.exists()) {
			finalJobDirectory.mkdir();
		}
		logFile.writeInfoMessage("[Command cleaning up directory " + jobDirectory.getAbsolutePath());
		
		//Delete files
		for (File df: deleteList) {
			this.deleteFile(df);
		}
		
		//Tomato sometimes relaunches old files, rename cmd.txt to avoid this
		File cmdFileOrig = new File(jobDirectory,"cmd.txt");
		if (cmdFileOrig.exists()) {
			File cmdFileRename = new File(jobDirectory,"cmd_done.txt");
			moveFile(cmdFileOrig, cmdFileRename, logFile);
		}
		
		logFile.writeInfoMessage("[Command] moving job directory " + jobDirectory.getAbsolutePath() + " to " + finalJobDirectory.getAbsolutePath());
		moveFile(jobDirectory,finalJobDirectory,logFile);
		
	}
	
	public CThread prepareJobDirectory(int threadNumber, int heartbeat, int failMax) {
		if (!jobDirectory.exists()) {
			jobDirectory.mkdir();
		}
		
		for (File inputFile: inputFiles) {
			File destFile = new File(jobDirectory,inputFile.getName());
			createLink(inputFile,destFile);
		}
		
		for (File externalFile: externalFiles) {
			File destFile = new File(jobDirectory,externalFile.getName());
			createLink(externalFile,destFile);
		}	
		
		createCmd();
		
		CThread thread = new CThread(jobDirectory, failMax, threadNumber, heartbeat, getProtectList(), logFile, email);
		return thread;
	}
	
	protected void merge(Command c2) {
		addOutputFiles(c2.getOutputFiles());
		addToLocalProperties(c2.getLocalProperties());
		addToTemplateString(c2.getTemplateString());
		addProtectList(c2.getProtectList());
		addDeleteList(c2.getDeleteList());
		addExternalFileList(c2.getExternalFiles());
		addSkipList(c2.getSkipFiles());
		addConditionalSkipList(c2.getConditionalSkipList());
		this.finalDirectoryList.add(c2.getFinalDirectory());
	}
	
    public boolean doesFinalExist() {
    	boolean doesExist = true;
    	
    	//Check for conditional output files and update skip list
    	for (File keyFile: conditionalSkipList.keySet()) {
    		if (keyFile.exists()) {
    			for (File newSkip: conditionalSkipList.get(keyFile)) {
    				skipList.add(newSkip);
    			}
    		}
    	}
    	
    	//Check for the existance of all output files, ignoring files that are marked as skipped.
    	for (File outputFile: outputFiles) {
    		if (skipList.contains(outputFile)) {
    			continue;
    		}
    		if (!outputFile.exists()) {
    			doesExist = false;
    		}
    	}
    	return doesExist;
    }
    
    public boolean doesWorkingExist() {
    	boolean doesExist = true;
    	for (File outputFile: outputFiles) {
    		String basename = outputFile.getName();
    		File workingFile = new File(jobDirectory,basename);
    		if (!workingFile.exists()) {
    			doesExist = false;
    		}
    	}
    	return doesExist;
    }
    
    
    
	//Getters and Setters
    public File getFinalDirectory() {
    	return this.finalDirectory;
    }
    
    public File getProjectDirectory() {
    	return this.projectDirectory;
    }
    
    public void setFinalDirectory(File finalDirectory) {
    	this.finalDirectory = finalDirectory;
    	this.finalDirectoryList.add(finalDirectory);
    }
    
    public void setJobDirectory(File f) {
    	this.jobDirectory = f;
    }
    
    public String getSampleName() {
    	return this.sampleName;
    }
    
  
    public String getSampleId() {
    	return this.sampleID;
    }
    
    public void setSampleName(String sampleName) {
    	this.sampleName = sampleName;
    }
    
    public void setSampleID(String sampleID) {
    	this.sampleID = sampleID;
    }
    
 
	protected ArrayList<File> getOutputFiles() {
		return outputFiles;
	}
	
	protected void addOutputFiles(ArrayList<File> files) {
		outputFiles.addAll(files);
	}
	
	protected HashMap<String,String> getLocalProperties() {
		return localProperties.get(0);
	}
	
	private void addToLocalProperties(HashMap<String,String> properties) {
		localProperties.add(properties);
	}
	
	protected void addLocalProperty(String key, String value) {
		if (this.localProperties.isEmpty()) {
			HashMap<String,String> newProperties = new HashMap<String,String>();
			this.localProperties.add(newProperties);
		}
		this.localProperties.get(0).put(key, value);
	}
	
	protected void setGlobalProperties(HashMap<String,String> properties) {
		this.globalProperties = properties;
	}
	
	protected String getTemplateString() {
		return this.templateStringList.get(0);
	}
	
	protected void setTemplateFile(File template) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(template));
			String temp = "";
			String templateString = "";
			while ((temp = br.readLine()) != null) {
				templateString += temp + "\n";
			}
			
			templateStringList.add(templateString);
			br.close();
		} catch (IOException ex) {
			logFile.writeErrorMessage("Error parsing template file: " + ex.getMessage(), false);
		}
	}
	
	protected void addToTemplateString(String template) {
		this.templateStringList.add(template);
	}
	
	protected ArrayList<File> getProtectList() {
		return protectList;
	}
	
	protected ArrayList<File> getDeleteList() {
		return deleteList;
	}
	
	protected void addProtect(File protect) {
		String name = protect.getName();
		File destFile = new File(jobDirectory,name);
		if (!this.protectList.contains(destFile)) {
			this.protectList.add(destFile);
		}
		
	}
	
	protected void addProtectList(ArrayList<File> protectList) {
		for (File f: protectList) {
			addProtect(f);
		}
	}
	
	protected void addConditionalSkip(File keyFile, File skipFile) {
		if (!conditionalSkipList.containsKey(keyFile)) {
			conditionalSkipList.put(keyFile, new ArrayList<File>());
		}
		conditionalSkipList.get(keyFile).add(skipFile);
	}
	
	protected HashMap<File,ArrayList<File>> getConditionalSkipList() {
		return this.conditionalSkipList;
	}
	
	protected void addSkipFile(File skipFile) {
		this.skipList.add(skipFile);
	}
	
	protected void addSkipList(ArrayList<File> skipFiles) {
		this.skipList.addAll(skipFiles);
	}
	
	protected ArrayList<File> getSkipFiles() {
		return this.skipList;
	}
	
	protected void addExternalFile(File external) {
		if (!this.externalFiles.contains(external)) {
			this.externalFiles.add(external);
		}
	}
	
	protected ArrayList<File> getExternalFiles() {
		return externalFiles;
	}
	
	protected void addExternalFileList(ArrayList<File> externalFiles) {
		for (File ef: externalFiles) {
			this.externalFiles.add(ef);
		}
		
	}
	
	protected void addConditionalSkipList(HashMap<File,ArrayList<File>> skipList) {
		for (File keyFile: skipList.keySet()) {
			if (!this.conditionalSkipList.containsKey(keyFile)) {
				this.conditionalSkipList.put(keyFile, new ArrayList<File>());
			}
			this.conditionalSkipList.get(keyFile).addAll(skipList.get(keyFile));
		}
	}
	
	protected void addDelete(File delete) {
		String name = delete.getName();
		File destFile = new File(jobDirectory,name);
		if (!this.deleteList.contains(destFile)) {
			this.deleteList.add(destFile);
		}
		
	}
	
	protected void addDeleteList(ArrayList<File> deleteList) {
		for (File f: deleteList) {
			addDelete(f);
		}
	}
	
	protected void addInputFiles(ArrayList<File> inputFiles) {
		this.inputFiles.addAll(inputFiles);
		
	}
	
	protected void addOutputFile(File outputFile) {
		this.outputFiles.add(outputFile);
	}

    
	protected static HashMap<String,File> matchFiles(Pattern pattern, File[] contents) {
		HashMap<String, File> matchingFiles = new HashMap<String,File>();
		for (File f: contents) {
			Matcher m = pattern.matcher(f.getName());
			if (m.matches()) {
				matchingFiles.put(m.group(1), f);
			}
		}
		return matchingFiles;
	}

	    
    //Methods shared by all commands
    protected void sortVCF(File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[sortVCF] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			ProcessBuilder pb = new ProcessBuilder(globalProperties.get("VCF_PATH_LOCAL") + "/vcf-sort",source.getAbsolutePath());
			
			Process p = pb.start();
			
			BufferedInputStream bis = new BufferedInputStream(p.getInputStream());
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(dest));
			
			
			byte[] buffer = new byte[1024*1024*10];
			int n = -1;
			
			while((n = bis.read(buffer))!=-1) {
			  bos.write(buffer,0,n);
			}
		

			int val = p.waitFor();
			bos.close();
			bis.close();

			
			if (val != 0) {
				logFile.writeErrorMessage("[sortVcf] Error while sorting the VCF file: " + dest.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[sortVcf] IO Exception while trying to move your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[sortVcf] Process was interrupted while trying to move your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
    
    protected void deleteFile(File source) {
    	if (source.exists()) {
    		try {
    		
    			ProcessBuilder pb = new ProcessBuilder("rm",source.getAbsolutePath());
    			Process p = pb.start();
    			
    			int val = p.waitFor();
    			
    			if (val != 0) {
    				logFile.writeErrorMessage("[deleteFile] System could not delete your file " + source.getAbsolutePath(),true);
    				System.exit(1);
    			}
    			
    		} catch (IOException ioex) {
    			logFile.writeErrorMessage("[deleteFile] IO Exception while trying to delete your file: " + source.getAbsolutePath(),true);
    			System.exit(1);
    		} catch (InterruptedException ieex) {
    			logFile.writeErrorMessage("[deleteFile] Process was interrupted while trying to delete your file: " + source.getAbsolutePath(),true);
    			System.exit(1);
    		}
    	} else {
			
    	}
    	
    }
    
    protected static void moveFile(File source, File dest, Logger logFile) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[moveFile] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			ProcessBuilder pb = new ProcessBuilder("mv",source.getAbsolutePath(),dest.getAbsolutePath());
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[moveFile] System could not move your file " + source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[moveFile] IO Exception while trying to move your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[moveFile] Process was interrupted while trying to move your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
    
    protected void cpFile(File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[copyFile] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			ProcessBuilder pb = new ProcessBuilder("cp","-r",source.getAbsolutePath(),dest.getAbsolutePath());
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[copyFile] System could not copy your file " + source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[copyFile] IO Exception while trying copy your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[copyFile] Process was interrupted while trying to copy your file: " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
	
	protected void createLink(File filename, File linkname) {
		try {
			if (linkname.exists()) {
				this.deleteFile(linkname);
			}
				
			ProcessBuilder pb = new ProcessBuilder("ln","-s",filename.getAbsolutePath(),linkname.getAbsolutePath());
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[createLink]  System could not create a link to your file " + filename.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
				
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[createLink]  IO Exception while trying to create a link to your file: " + filename.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[createLink]  Process was interrupted while trying to create a link to your file: " + filename.getAbsolutePath(),true);
			System.exit(1);
		}
		
	}
	
	
	
	protected void mergeVcf(ArrayList<File> vcfList, File dest) {
		try {
			ArrayList<String> command  = new ArrayList<String>();
			command.add(globalProperties.get("VCF_PATH_LOCAL") + "/vcf-concat");

			for (File bam: vcfList) {
				if (!bam.exists()) {
					logFile.writeErrorMessage("[mergeVcf] Expected file does not exist: " + bam.getAbsolutePath(),true);
					System.exit(1);
				}
				command.add(bam.getAbsolutePath());
			}
			
			
			ProcessBuilder pb = new ProcessBuilder(command);
			Process p = pb.start();
			
			BufferedInputStream bis = new BufferedInputStream(p.getInputStream());
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(dest));
			
			
			byte[] buffer = new byte[1024*1024*10];
			int n = -1;
			
			while((n = bis.read(buffer))!=-1) {
			  bos.write(buffer,0,n);
			}
		

			int val = p.waitFor();
			bos.close();
			bis.close();

			
			if (val != 0) {
				logFile.writeErrorMessage("[mergeVcf] Error while merging the VCF files: " + dest.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
			
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[mergeVcf] IO Exception while trying to merge the VCF files: " + dest.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[mergeVcf] Process was interrupted while trying to merge the VCF files: " + dest.getAbsolutePath(),true);
			System.exit(1);
		}
	}
	
	protected void bgzip(File source) {
		try {
			
			if (!source.exists()) {
				logFile.writeErrorMessage("[bgzip] Expected file does not exist: " + source.getAbsolutePath(),true);
			}
			
			ProcessBuilder pb = new ProcessBuilder(globalProperties.get("TABIX_PATH_LOCAL") + "/bgzip",source.getAbsolutePath());
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[bgzip] Error while bgzipping the vcf: " + source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[bgzip] IO Exception while trying to bgzip the VCF files: " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[bgzip] Process was interrupted while trying to bgzip the VCF files: " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
	
	protected void tabix(File source) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[tabix] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			ProcessBuilder pb = new ProcessBuilder(globalProperties.get("TABIX_PATH_LOCAL") + "tabix","-p","vcf",source.getAbsolutePath());
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[tabix] Error while indexing the VCF files: " + source.getAbsolutePath(),true);
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[tabix] IO Exception while trying to index the VCF files: " + source.getAbsolutePath(),true);
			System.exit(1);
		} catch (InterruptedException ieex) {
			logFile.writeErrorMessage("[tabix] Process was interrupted while trying to indexthe VCF files: " + source.getAbsolutePath(),true);
			System.exit(1);
		}
	}
	
	
	protected void deleteFolder(File folder) {
	    File[] files = folder.listFiles();
	    if(files!=null) { //some JVMs return null for empty dirs
	        for(File f: files) {
	            if(f.isDirectory()) {
	                deleteFolder(f);
	            } else {
	                f.delete();
	            }
	        }
	    }
	    folder.delete();
	}
   
	
    protected void createCmd() {
		//Write command file
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(jobDirectory,"cmd.txt")));
			
			
			if (this.suppress) {
				bw.write(String.format("#e %s -abe\n", this.email));
				//bw.write(String.format("#e %s -n\n#c ember\n", this.email));
			} else {
				bw.write(String.format("#e %s \n", this.email));
				//bw.write(String.format("#e %s \n#c ember\n", this.email));
			}
			
			if (this.cluster != null) {
				bw.write(String.format("#c " + this.cluster));
			}
			
			if (wallTime != null) {
				bw.write(String.format("#w %s\n", this.wallTime));
			}
			
			bw.write("\n\n");
			
			//Determine node metrics
			bw.write("MEMTOTAL=`free | grep Mem | awk '{ print $2 }'`\n");
			bw.write("MEMGB=`expr $MEMTOTAL / 1048576`\n");
			bw.write("SMGB=`expr $MEMGB - 2`\n");
			bw.write("SMGB2=`expr $MEMGB / 2`\n");
			bw.write("NCPU=`nproc`\n");
			bw.write("GCT=$NCPU\n");
			bw.write("DISKAVAIL=`df -hP . | tail -n 1 | awk '{print $4}'`\n");
			
			bw.write("if [ $NCPU -gt 8 ]\n");
			bw.write("then\n");
			bw.write("   GCT=`expr $NCPU \\* 5 / 8 + 3`\n");
			bw.write("fi\n");
			bw.write("HOST=`hostname`\n");
			bw.write("echo \"Hostname: \" $HOST\n");
			bw.write("echo \"Total cpu: \" $NCPU\n");
			bw.write("echo \"GC threads: \" $GCT\n");
			bw.write("echo \"Total memory: \" $MEMGB\n");
			bw.write("echo \"Java memory: \" $SMGB\n");
			bw.write("echo \"Disk available: \" $DISKAVAIL\n");
			bw.write("\n\n");
			
		

			for (int i=0; i<templateStringList.size();i++) {
				String edited = templateStringList.get(i);
				for (String key: localProperties.get(i).keySet()) {
					edited = edited.replaceAll(key, localProperties.get(i).get(key));
				}
				for (String key: globalProperties.keySet()) {
					edited = edited.replaceAll(key, globalProperties.get(key));
				}
				bw.write(edited + "\n");
				
			}
		
	
			bw.close();
			
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[TFCommand] Error writing command file, exiting: " + ioex.getMessage(),true);
			System.exit(1);
		}	
		
	}
    
  
    protected static void validateFastqSet(File file1, File file2, Logger logFile) {
		logFile.writeInfoMessage("[TFCommand] Validating paired-end file 1: " + file1.getAbsolutePath());
		int lines1 = validateFastq(file1, logFile);
		if (lines1 == 0) {
			logFile.writeErrorMessage("Paired-end fastq file 1 is empty",false);
			System.exit(1);
		}
		logFile.writeInfoMessage("[TFCommand] Validated " + lines1 + " records");
		logFile.writeInfoMessage("[TFCommand] Validating paired-end file 2");
		int lines2 = validateFastq(file2, logFile);
		if (lines2 == 0) {
			logFile.writeErrorMessage("[TFCommand] Paired-end fastq file 2 is empty",false);
			System.exit(1);
		}
		logFile.writeInfoMessage("[TFCommand] Validated " + lines2 + " records");
		if (lines2 != lines1) {
			logFile.writeErrorMessage("[TFCommand] Paired-end fastq files have unequal numbers of lines: " + lines1 + " " + lines2,false);
			System.exit(1);
		}
	}
	
	protected static int validateFastq(File file, Logger logFile) {
		int retval = 0;
		try {
			GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(file));
			BufferedReader br = new BufferedReader(new InputStreamReader(gzip));
			int counter = 0;
			int readLength = 0;
			Pattern p = Pattern.compile("^[ACTGN]+$");
			String line = null;
			int min = 100000000;
			int max = 0;
			int rl = 0;
			
		
			while((line = br.readLine()) != null) {
				if (counter % 4 == 0) {
					if (!line.startsWith("@")) {
						logFile.writeErrorMessage("[TFCommand] Fastq line 1 does not start with @,  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
				} else if (counter % 4 == 1) {
					Matcher m = p.matcher(line);
					readLength = line.length();
					if (readLength > rl) {
						rl = readLength;
					}
					if (!m.matches()) {
						logFile.writeErrorMessage("[TFCommand] Fastq line 2 has characters other than A,C,G,T or N.  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
				} else if (counter % 4 == 2) {
					if (!line.startsWith("+")) {
						logFile.writeErrorMessage("[TFCommand] Fastq line 3 does not start with +.  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
				} else if (counter % 4 == 3) {
					if (line.length() != readLength) {
						logFile.writeErrorMessage("[TFCommand] Fastq line 4 length does not match fastq line 2,  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
					//Quality check
					for (int i = 0; i < line.length(); i++){
					    char c = line.charAt(i);        
					  
					    int val = (int)c;
					    if (val > max) {
					    	max = val;
					    } else if (val < min) {
					    	min = val;
					    }
					    
					}
				} else {
					logFile.writeErrorMessage("[TFCommand] Should not reach this line",true);
					System.exit(1);
				}
				counter++;
			}
			
			if (min-33 >= 31) {
				logFile.writeInfoMessage("[TFCommand] Fastq detected as ASCII-64.  Min: " + (min-64) + ". Max: " + (max-64));
				
			} else {
				logFile.writeInfoMessage("[TFCommand] Fastq detected as ASCII-33.  Min: " + (min-33) + ". Max: " + (max-33));
				
			}
			
			
			if (counter % 4 != 0) {
				logFile.writeErrorMessage("[TFCommand] Fastq file number not divible by 4, exiting",false);
				System.exit(1);
			}
			
			br.close();
			retval = counter / 4;
			
			
		} catch (Exception ioex) {
			logFile.writeErrorMessage("[TFCommand] Error reading fastq file, exiting: " + ioex.getMessage(),true) ;
			ioex.printStackTrace();
			System.exit(1);
		}
		return retval;
	}
}