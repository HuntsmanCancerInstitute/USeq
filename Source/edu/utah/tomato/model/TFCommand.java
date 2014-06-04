package edu.utah.tomato.model;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import edu.utah.tomato.daemon.TFThreadDaemon;
import edu.utah.tomato.model.command.TFCommandExomeAlignUgp;
import edu.utah.tomato.model.command.TFCommandExomeAlignHci;
import edu.utah.tomato.model.command.TFCommandExomeMetrics;
import edu.utah.tomato.model.command.TFCommandExomeVariant;
import edu.utah.tomato.util.TFConstants;
import edu.utah.tomato.util.TFLogger;


public abstract class TFCommand {
	//Type of command
	protected String commandString = null;
	protected String commandType = null;
	
	//Directory locations
	protected ArrayList<File> templateFiles = null;
	protected File rootDirectory = null;
	
	//Log file
	protected TFLogger logFile = null;

	//task list
	protected List<Callable<Boolean>> taskList = new ArrayList<Callable<Boolean>>();
	
	//cmd file options
	protected Integer wallTime = 0;
	protected String email = null;
	protected boolean suppress = false;
	protected boolean isFull = false;
	protected Integer jobs= null;
	
	//Tomato thread options
	protected Integer heartbeat = null;
	protected Integer failmax = null;
	
	//Thread daemon
	protected TFThreadDaemon daemon = null;

	//Variables
	protected HashMap<String,String> properties = null;
	
	//Job directories
	protected File finalDirectory = null;
	protected File jobDirectory = null;
	
	public static TFCommand getCommandObject(File rootDirectory, String commandLine, TFLogger logFile, 
			String email, Integer wallTime, Integer heartbeat, Integer failmax, Integer jobs, boolean suppress,
			boolean deleteMetricsBams, boolean deleteReducedBams, boolean isFull, boolean use1KGenomes, boolean validateFastq, String study, String splitType,  
			File targetFile, HashMap<String,String> properties) {
		
		//Validate template file
		File templateDir = new File(properties.get("TEMPLATE_PATH"));
		String[] commandParts = commandLine.split(":");
		if (commandParts.length < 2) {
			logFile.writeErrorMessage("[TFCommand] Must be at least two parts for every command", true);
			System.exit(1);
		}
		ArrayList<File> templateFiles = new ArrayList<File>();
		File templateFile = new File(templateDir,"templates/" + commandParts[0] + "/" + commandParts[0] + "." + commandParts[1] + ".txt");
		if (!templateFile.exists()) {
			logFile.writeErrorMessage("[TFCommand] Configuration file invalid, specified template file doesn not exist: " + commandLine,false);
			System.exit(1);
		}
		templateFiles.add(templateFile);
		
		//Check for dependencies
		File depFile = new File(templateDir,"templates/" + commandParts[0] + "/deps." + commandParts[1] + ".txt");
		if (depFile.exists()) {
			logFile.writeInfoMessage("[TFCommand] Found dependencies: ");
			try {
				BufferedReader br = new BufferedReader(new FileReader(depFile));
				String line = "";
				while ((line = br.readLine()) != null) {
					File extraFile = new File(templateDir,"templates/" + commandParts[0] + "/" + line);
					if (!extraFile.exists()) {
						logFile.writeErrorMessage("[TFCommand] Configuration file invalid, specified template file doesn not exist: " + extraFile.toString(),false);
						System.exit(1);
					}
					logFile.writeInfoMessage("[TFCommand] \t" + extraFile.getName());
					templateFiles.add(extraFile);
				}
				br.close();
			} catch (FileNotFoundException fnfe) {
				logFile.writeErrorMessage("[TFCommand] Count not find deps file", true);
			} catch (IOException ioex) {
				logFile.writeErrorMessage("[TFCommand]  Error reading deps file",true);
			}
		}
		
		TFCommand returnCommand = null;
		if (commandParts[0].equals(TFConstants.ANALYSIS_EXOME_ALIGN_NOVO) || commandParts[0].equals(TFConstants.ANALYSIS_EXOME_ALIGN_BWA) ) {
			if (failmax == null) {
				failmax = 3;
			}
			returnCommand = new TFCommandExomeAlignUgp(templateFiles,rootDirectory, commandLine, commandParts[0], logFile,email, wallTime, 
					heartbeat, failmax, jobs, suppress, isFull, validateFastq, properties);
		} else if (commandParts[0].equals(TFConstants.ANALYSIS_EXOME_ALIGN_BEST)) {
			if (failmax == null) {
				failmax = 3;
			}
			returnCommand = new TFCommandExomeAlignHci(templateFiles,rootDirectory, commandLine, commandParts[0], logFile,email, wallTime, 
					heartbeat, failmax, jobs, suppress, isFull, validateFastq, properties);
		}
		
		else if (commandParts[0].equals(TFConstants.ANALYSIS_EXOME_METRICS)) {
			if (failmax == null) {
				failmax = 3;
			}
			returnCommand = new TFCommandExomeMetrics(templateFiles,rootDirectory, commandLine, commandParts[0], logFile,email, wallTime, 
					heartbeat, failmax, jobs, suppress, deleteMetricsBams, isFull, study, targetFile, properties);
		} else if (commandParts[0].equals(TFConstants.ANALYSIS_EXOME_VARIANT_RAW) || commandParts[0].equals(TFConstants.ANALYSIS_EXOME_VARIANT_VQSR) || commandParts[0].equals(TFConstants.ANALYSIS_EXOME_VARIANT_BEST)) {
			if (failmax == null) {
				failmax = 5;
			}
			returnCommand = new TFCommandExomeVariant(templateFiles,rootDirectory, commandLine, commandParts[0], logFile,email, wallTime, 
					heartbeat, failmax, jobs, suppress, deleteReducedBams, isFull, use1KGenomes, study, splitType, targetFile, properties);
		}
		
		else {
			logFile.writeErrorMessage("[TFCommand] Command type not recognized: " + commandParts[0], true);
			System.exit(1);
		}
		
		logFile.writeInfoMessage("[TFCommand] " + commandLine + " is allowed " + failmax + " failures");
		
		return returnCommand;
	}
	

	public TFCommand(ArrayList<File> templateFile, File rootDirectory, String commandString, String commandType, TFLogger logFile, 
			String email, Integer wallTime, Integer heartbeat, Integer failmax, Integer jobs, boolean suppress, boolean isFull, HashMap<String,String> properties) {
		this.templateFiles = templateFile;
		this.rootDirectory  = rootDirectory;
		this.commandString = commandString;
		this.commandType = commandType;
		this.logFile = logFile;
		this.wallTime = wallTime;
		this.email = email;
		this.heartbeat = heartbeat;
		this.failmax = failmax;
		this.jobs = jobs;
		this.suppress  = suppress;
		this.properties = properties;
		this.isFull = isFull;
	}
	
	
	/** Create file links needed based on the si object */
	public abstract ArrayList<TFSampleInfo> run(ArrayList<TFSampleInfo> sampleList);
	
	protected abstract ArrayList<TFSampleInfo> validateSampleSet(ArrayList<TFSampleInfo> sampleList);
	protected abstract ArrayList<TFSampleInfo> findPrereqsExisting(ArrayList<TFSampleInfo> sampleList);
	protected abstract ArrayList<TFSampleInfo> findPrereqsNew(ArrayList<TFSampleInfo> sampleList);
	
    public String getCommandType() {
    	return this.commandType;
    }
    
    public String getCommandString() {
    	return this.commandString;
    }
    
    public void shutdown() {
    	if (this.daemon != null && this.daemon.isAlive()) {
    		logFile.writeWarningMessage("[TFCommand] Recieved termination signal, interrupting daemon");
    		this.daemon.interrupt();
    		while(this.daemon.isAlive()){}
    	}
    }
    
    protected ArrayList<TFSampleInfo> findPatternsExisting(ArrayList<TFSampleInfo> sampleList, File directory, ArrayList<TFMatchObject> dependantPatterns) {
    	//Iterate through file in the directory
    	File[] contents = directory.listFiles();
    	
    	for (File file: contents) {
			if (file.isDirectory()) {
				continue;
			}
			
			for (TFMatchObject pattern: dependantPatterns) {
				if (pattern.testMatch(file.getName())) {
					sampleList = this.matchDependantPatterns(file,sampleList,pattern,directory);
				}
			}
		}
		
		sampleList = validateSampleSet(sampleList);
		
		return sampleList;
    }
    
    /**
     * This method searches the directory for matching files.  MasterPatterns are used to initialize SampleInfo objects.
     * DependantPatterns are then added to already-initialized SampleInfo objects.  The validateSampleSet method
     * is called at the end to make sure the necessary files were found.
     * 
     * @param directory
     * @param masterPatterns
     * @param dependantPatterns
     * @return validatedSampleSets
     */
    protected ArrayList<TFSampleInfo> findPatternsNew (File directory, ArrayList<TFMatchObject> masterPatterns, ArrayList<TFMatchObject> dependantPatterns) {
		//Create containers for valid/ignored files
		ArrayList<TFSampleInfo> sampleList = new ArrayList<TFSampleInfo>();
				
		//Iterate through file in the directory
		File[] contents = directory.listFiles();
		
		for (File file: contents) {
			//Skip directories
			if (file.isDirectory()) {
				continue;
			}
		
			for (TFMatchObject pattern: masterPatterns) {
				if (pattern.testMatch(file.getName())) {
					sampleList = this.matchMasterPatterns(file, sampleList, pattern, directory);
				}
			}
		}
		
		for (File file: contents) {
			if (file.isDirectory()) {
				continue;
			}
			
			for (TFMatchObject pattern: dependantPatterns) {
				if (pattern.testMatch(file.getName())) {
					sampleList = this.matchDependantPatterns(file,sampleList,pattern,directory);
				}
			}
		}
		
		sampleList = validateSampleSet(sampleList);
		
		return sampleList;
	}
    
    /** This method finds files that match the specified patterns.  If the pattern wasn't previously observed
     * a new SampleInfo object is created. 
     * 
     * @param f
     * @param sampleList
     * @param p
     * @param directory
     */
    private ArrayList<TFSampleInfo> matchMasterPatterns(File f, ArrayList<TFSampleInfo> sampleList, TFMatchObject p, File directory) {
 
		TFSampleInfo used = null;
		
		for (TFSampleInfo sample: sampleList) {
			if (sample.comparePrefix(p.getMatchPrefix())) {
				used = sample;
			}
		}
		if (used == null) {
			String sampleID = "unknown";
			String puID = "unknown";
			String sampleName = "unknown";
			if (p.getPrefixType().equals(TFConstants.PREFIX_SAMPLEID)) {
				sampleID = p.getMatchPrefix();
				String[] idParts = sampleID.split("_",2);
				
				if (idParts.length != 2) {
			    	logFile.writeErrorMessage("[TFCommand] File name improperly formed. Should have an underscore separating the sample name and the flowcell name: SAMPLE_FLOWCELL_1.fastq and SAMPLE_FLOWCELL_2.fastq",false);
			    	System.exit(1);
			    }
				
				sampleName = idParts[0];
				puID = idParts[1];
			} else if (p.getPrefixType().equals(TFConstants.PREFIX_SAMPLENAME)) {
				sampleName = p.getMatchPrefix();
			} else if (p.getPrefixType().equals(TFConstants.PREFIX_PUID)) {
				puID = p.getMatchPrefix();
			}
			
			
			//Create a new SampleInfoObject
			used = new TFSampleInfo(sampleID,sampleName,puID,logFile);
			sampleList.add(used);
		} 
		
		TFFileObject tempFO = new TFFileObject(p.getMatchFull(),directory,directory);
		used.setFileObject(p.getMatchName(), tempFO);
		
		return sampleList;
		
	}
    
    /** This method adds the files that match specified patterns to an existing sampleList.
     * 
     * @param f
     * @param sampleList
     * @param p
     * @param directory
     */
    private ArrayList<TFSampleInfo> matchDependantPatterns(File f, ArrayList<TFSampleInfo> sampleList, TFMatchObject p, File directory) {
    	TFSampleInfo used = null;
    	
    	for (TFSampleInfo sample: sampleList) {
    		if (p.getPrefixType().equals(TFConstants.PREFIX_SAMPLEID)) {
    			if (sample.getSampleID().equals(p.getMatchPrefix())) {
    				used = sample;
    			}
    		} else if (p.getPrefixType().equals(TFConstants.PREFIX_SAMPLENAME)) {
    			if (sample.getSampleName().equals(p.getMatchPrefix())) {
    				used = sample;
    			}
    			
    		} else if (p.getPrefixType().equals(TFConstants.PREFIX_PUID)) {
    			if (sample.getPuID().equals(p.getMatchPrefix())) {
    				used = sample;
    			}
    		} else {
    			logFile.writeErrorMessage(String.format("[TFCommand] The specified StringPattern matchtype isn't recognized: %s", p.getPrefixType()), true);
    			System.exit(1);
    		}
    	}
    	
    	if (used != null) {
    		TFFileObject tempFO = new TFFileObject(p.getMatchFull(),directory,directory);
			used.setFileObject(p.getMatchName(), tempFO);
    	}
    	
    	return sampleList;
    }
	
	
	    
    //Methods shared by all commands
    protected void sortVCF(File source, File dest) {
		try {
			if (!source.exists()) {
				logFile.writeErrorMessage("[sortVCF] Expected file does not exist: " + source.getAbsolutePath(),true);
				System.exit(1);
			}
			
			ProcessBuilder pb = new ProcessBuilder(properties.get("VCF_PATH_LOCAL") + "/vcf-sort",source.getAbsolutePath());
			
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
			//logFile.writeErrorMessage("[deleteFile] Expected file does not exist: " + source.getAbsolutePath(),true);
			//System.exit(1);
    	}
    	
    }
    
    protected void moveFile(File source, File dest) {
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
//				if (linkname.getCanonicalFile().equals(filename.getAbsoluteFile())) {
//					logFile.writeWarningMessage("[createLink] An existing link matches what you requested, using existing link");
//				} else {
//					logFile.writeErrorMessage("[createLink] An existing file exists at link location and it doesn't match request, refusing to overwrite existing file"
//							+ "\nExisting: " + linkname.getCanonicalPath()
//							+ "\nNew: " + filename.getAbsolutePath(),false);
//					System.exit(1);
//				}
			}
				
			ProcessBuilder pb = new ProcessBuilder("ln","-s",filename.getAbsolutePath(),linkname.getAbsolutePath());
			Process p = pb.start();
			
			int val = p.waitFor();
			
			if (val != 0) {
				logFile.writeErrorMessage("[createLink]  System could not create a link to you file " + filename.getAbsolutePath(),true);
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
			command.add(properties.get("VCF_PATH_LOCAL") + "/vcf-concat");

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
			
			ProcessBuilder pb = new ProcessBuilder(properties.get("TABIX_PATH_LOCAL") + "/bgzip",source.getAbsolutePath());
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
			
			ProcessBuilder pb = new ProcessBuilder(properties.get("TABIX_PATH_LOCAL") + "tabix","-p","vcf",source.getAbsolutePath());
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
   
	
    protected void createCmd(HashMap<String,String> replacements, File cmdFile, int index) {
		//Write command file
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.templateFiles.get(index)));
			BufferedWriter bw = new BufferedWriter(new FileWriter(cmdFile));
			
			
			if (this.suppress) {
				bw.write(String.format("#e %s -n\n", this.email));
				//bw.write(String.format("#e %s -n\n#c ember\n", this.email));
			} else {
				bw.write(String.format("#e %s \n", this.email));
				//bw.write(String.format("#e %s \n#c ember\n", this.email));
			}
			
			if (this.wallTime != null) {
				bw.write(String.format("#w %s\n", this.wallTime));
			}
			
			bw.write("\n\n");
		
			
			//Write command file
			String line = null;

			
			while((line = br.readLine()) != null) {
				String edited = line;
				for (String key: replacements.keySet()) {
				
					edited = edited.replaceAll(key, replacements.get(key));
				}
				bw.write(edited + "\n");
			}
			
			
			br.close();
			bw.close();
			
			
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[TFCommand] Error writing command file, exiting: " + ioex.getMessage(),true);
			System.exit(1);
		}	
		
	}
    
    protected void cleanup(HashSet<File> runDirectoryList, HashSet<File> deleteList) {
		//clean up unwanted files
		for (File df: deleteList) {
			this.deleteFile(df);
		}
		
		//move job directories
		for (File rd: runDirectoryList) {
			//Tomato sometimes relaunches old files, rename cmd.txt to avoid this
			File cmdFileOrig = new File(rd,"cmd.txt");
			if (cmdFileOrig.exists()) {
				File cmdFileRename = new File(rd,"cmd_done.txt");
				this.moveFile(cmdFileOrig, cmdFileRename);
			}
			
			
			File existDir = new File(this.jobDirectory,rd.getName());
			if (existDir.exists()) {
				deleteFolder(existDir);
			}
			this.moveFile(rd, this.jobDirectory);
			
		}
	}
    
    protected void validateFileSet(TFSampleInfo fi) {
		logFile.writeInfoMessage("[TFCommand] Validating fastq files for sample: " + fi.getSampleName());
		logFile.writeInfoMessage("[TFCommand] Validating paired-end file 1");
		int lines1 = validateFastq(fi,TFConstants.FILE_FASTQ1);
		if (lines1 == 0) {
			logFile.writeErrorMessage("Paired-end fastq file 1 is empty",false);
			System.exit(1);
		}
		logFile.writeInfoMessage("[TFCommand] Validated " + lines1 + " records");
		logFile.writeInfoMessage("[TFCommand] Validating paired-end file 2");
		int lines2 = validateFastq(fi,TFConstants.FILE_FASTQ2);
		if (lines2 == 0) {
			logFile.writeErrorMessage("[TFCommand] Paired-end fastq file 2 is empty",false);
			System.exit(1);
		}
		logFile.writeInfoMessage("[TFCommand] Validated " + lines2 + " records");
		if (lines2 != lines1) {
			logFile.writeErrorMessage("[TFCommand] Paired-end fastq files have unequal numbers of lines: " + lines1 + " " + lines2,false);
			System.exit(1);
		}
		fi.setRecordCount(lines1);
	}
	
	protected int validateFastq(TFSampleInfo fi, String name) {
		int retval = 0;
		try {
			GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(fi.getFileObject(name).getFinalPath()));
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
				fi.setQual64();
			} else {
				logFile.writeInfoMessage("[TFCommand] Fastq detected as ASCII-33.  Min: " + (min-33) + ". Max: " + (max-33));
				
			}
			
			
			if (counter % 4 != 0) {
				logFile.writeErrorMessage("[TFCommand] Fastq file number not divible by 4, exiting",false);
				System.exit(1);
			}
			
			br.close();
			retval = counter / 4;
			
			fi.setReadLength(rl);
		} catch (Exception ioex) {
			logFile.writeErrorMessage("[TFCommand] Error reading fastq file, exiting: " + ioex.getMessage(),true) ;
			ioex.printStackTrace();
			System.exit(1);
		}
		return retval;
	}
}
