package util.apps;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

/*** Zip archives particular folders, deletes particular files. Use to clean up analysis results directories prior to cloud upload.*/
public class JobCleaner {

	//user defined fields
	private File rootDirectory = null;
	private String[] fileExtensionsToDelete = null;
	private String[] directoryNamesToZip = null;
	private boolean dryRun = false;
	private boolean mergeInDirectory = false;
	
	//internal
	private boolean lookForFiles;
	private boolean lookForDirectories;
	private ArrayList<File> toDelete = new ArrayList<File>();
	private ArrayList<File> toZip = new ArrayList<File>();

	public JobCleaner (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			IO.pl("Walking root directory...");
			walkDir(rootDirectory);
			
			deleteFiles();
			
			zipDirectories();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
			
		} catch (Exception e) {
			IO.el("\nERROR running the JobCleaner, aborting");
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void zipDirectories() throws IOException {
		if (lookForDirectories == false) return;
		if (dryRun) {
			IO.pl("\nThe following files would be zip archived and then deleted if this wasn't a dry run...");
			for (File f: toZip) IO.pl("\t"+f.toString());
		}
		else {

			IO.pl("\nZip archiving and then deleting the following files and directories...");
			if (mergeInDirectory) {
				//split by parent directory
				HashMap<String, ArrayList<File>> parentDirs = new HashMap<String, ArrayList<File>>();
				for (File d: toZip) {
					IO.pl("\t"+d.toString());
					if (d.exists()== false) IO.pl("\t\tMissing, already zipped?! Skipping "+d.toString());
					else {
						String parent = d.getParentFile().getCanonicalPath();
						ArrayList<File> al = parentDirs.get(parent);
						if (al == null) {
							al = new ArrayList<File>();
							parentDirs.put(parent, al);
						}
						al.add(d);
					}
				}
				//for each parent dir, zip the children dirs into a combine archive
				for (String par: parentDirs.keySet()) {
					ArrayList<File> children = parentDirs.get(par);
					File[] toCombine = new File[children.size()];
					children.toArray(toCombine);
					Arrays.sort(toCombine);
					StringBuilder comboName = new StringBuilder(par);
					comboName.append("/");
					for (File f: toCombine) comboName.append(f.getName());
					comboName.append(".zip");
					File zipFile = new File (comboName.toString());
					boolean zipped = IO.zipDirectoriesInSameParentDirectory(toCombine, zipFile);
					if (zipped == false) throw new IOException("Failed to zip: "+comboName);

					//delete the individual dirs
					for (File f: toCombine) {
						boolean deleted = IO.deleteDirectorySimple(f);
						if (deleted == false) throw new IOException("Failed to delete directory after zipping: "+f);
					}
					
				}
			}
			else {
				for (File d: toZip) {
					IO.pl("\t"+d.toString());
					if (d.exists()== false) IO.pl("\t\tMissing, already zipped?! "+d.toString());
					else {
						boolean zipped = IO.zipDirectory(d);
						if (zipped == false) throw new IOException("Failed to zip: "+d);

						boolean deleted = IO.deleteDirectorySimple(d);
						if (deleted == false) throw new IOException("Failed to delete directory after zipping: "+d);
					}
				}
			}
		}
	}

	private void deleteFiles() throws IOException {
		if (lookForFiles == false) return;
		if (dryRun) {
			IO.pl("\nThe following files would be deleted if this wasn't a dry run...");
			for (File f: toDelete) IO.pl("\t"+f.toString());
		}
		else {
			
			IO.pl("\nDeleting the following files...");
			for (File f: toDelete) {
				String name = f.toString();
				IO.pl("\t"+name);
				f.delete();
				//check if it exists, on Macs, the .delete() is returning false!
				if (f.exists())  throw new IOException("ERROR: failed to delete "+name);
			}
		}
	}

	public void walkDir (File directory){ 
		File[] list = directory.listFiles();
		if (list != null){
			for (int i=0; i< list.length; i++){
				//is it a directory?
				if (list[i].isDirectory()) {
					if (lookForDirectories) checkDir(list[i]);
					walkDir (list[i]);
				}
				//must be a file, do they want to delete anything?
				else if (lookForFiles) checkFile(list[i]);				
			}
		}
	}

	private void checkFile(File file) {
		String fileName = file.getName();
		for (String e: fileExtensionsToDelete) {
			if (fileName.endsWith(e)) {
				toDelete.add(file);
				return;
			}
		}
	}

	private void checkDir(File dir) {
		String dirName = dir.getName();
		for (String n: directoryNamesToZip) {
			if (dirName.equals(n)) {
				toZip.add(dir);
				return;
			}
		}
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new JobCleaner(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
			IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
			Pattern pat = Pattern.compile("-[a-zA-Z]");
			for (int i = 0; i<args.length; i++){
				Matcher mat = pat.matcher(args[i]);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'r': rootDirectory = new File(args[++i]); break;
						case 'e': fileExtensionsToDelete = Misc.COMMA.split(args[++i]); break;
						case 'n': directoryNamesToZip = Misc.COMMA.split(args[++i]); break;
						case 'd': dryRun = true; break;
						case 'm': mergeInDirectory = true; break;
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
						}
					}
					catch (Exception e){
						e.printStackTrace();
						Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					}
				}
			}
			
			lookForDirectories = directoryNamesToZip!=null;
			lookForFiles = fileExtensionsToDelete!=null;

			
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                Job Cleaner : Jan 2024                            **\n" +
				"**************************************************************************************\n" +
				"Zip archives particular folders, deletes particular files. Use to clean up analysis \n"+
				"result directories prior to cloud upload.\n"+
				
				"\nOptions:\n"+
				"-r Root directory to recursively look for files and folders.\n"+
				"-e File extensions and file names to delete, comma delimited, no spaces.\n" +
				"-n Directory names to zip archive and then delete, comma delimited, no spaces.\n"+
				"-m Create a merged zip archive for directories defined in -n that exist in the same\n"+
				"      parent directory, e.g. LogsRunScripts.zip instead of Logs.zip and RunScripts.zip\n"+
				"-d Dry run, just list the files and directories.\n"+


				"\nExample: java -jar pathToUSeq/Apps/JobCleaner -d -n 'Logs,RunScripts' -r CJobs/ -e \n"+
				"    '.tbi,.crai,.bai,COMPLETE' -m \n"+

				"\n**************************************************************************************\n");
	}

}
