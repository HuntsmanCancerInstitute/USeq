package util.bio.parsers;
import org.apache.commons.net.ftp.*;

import util.gen.*;
import edu.utah.seq.useq.USeqUtilities;
import java.io.*;
import java.util.*;
import java.util.regex.*;

/**Downloads and converts xxx.sra files from the SRA to fastq.gz, then launches tomato.  So SRAToDamnSAM!
 * @author davidnix*/
public class SRAProcessor {
	
	/* ftp://ftp-trace.ncbi.nlm.nih.gov

	## by run SRR016669
	/sra/sra-instant/reads/ByRun/sra/SRR/SRR016/SRR016669/SRR016669.sra

	## by project SRP000401
	/sra/sra-instant/reads/ByStudy/sra/SRP
	/sra/sra-instant/reads/ByStudy/sra/SRP/SRP000/SRP000401 many others possible
	/sra/sra-instant/reads/ByStudy/sra/SRP/SRP000/SRP000401/SRR016669/SRR016669.sra
	 */
	
	//fields
	private String ipAddress = "ftp-trace.ncbi.nlm.nih.gov";
	private String srrPath = "/sra/sra-instant/reads/ByRun/sra/SRR/";
	private String srpPath = "/sra/sra-instant/reads/ByStudy/sra/SRP/";
	private String email = "useq-users@lists.sourceforge.net";
	private FTPClient ftp;
	private File saveDirectory;
	private File fastqDump;
	private File cmdFile = null;
	private ArrayList<String> srpToProcess = new ArrayList<String>();
	private ArrayList<String> srrToProcess = new ArrayList<String>();
	private ArrayList<String> srrFailedToProcess = new ArrayList<String>();
	private ArrayList<String> srrProcessed = new ArrayList<String>();
	private static final Pattern JUST_SRR = Pattern.compile("(SRR\\d+).*");

	//constructors
	public SRAProcessor(String[] args) {
		long startTime = System.currentTimeMillis();
		processArgs(args);

		//connect
		connect();
		
		//any projects?
		if (srpToProcess.size() != 0){
			System.out.println("Fetching SRR names from project(s)...");
			for (String project: srpToProcess) fetchProjectSRRNames(project);
		}
		
		//look for any already converted
		System.out.println("\nLooking for existing converted SRR's, the following will be skipped, delete to reprocess...");
		removeAlreadyPresent();
		
		//process the runs
		String[] srrNamesToProcess = Misc.stringArrayListToStringArray(srrToProcess);
		System.out.println("\nProcessing "+srrNamesToProcess.length+" SRRs...");
		for (int i=0; i< srrNamesToProcess.length; i++){
			fetchConvertDeleteSRR(srrNamesToProcess[i]);
		}
		
		//disconnect
		disconnect();
		
		
		
		//stats
		printProcessed();
		
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(60000.0 * 60.0);
		System.out.println("\nRun complete "+Num.formatNumber(diffTime, 3)+"hrs");
	}
	
	public void removeAlreadyPresent(){
		//fetch all the files in the save directory
		File[] files = saveDirectory.listFiles();
		HashMap<String, ArrayList<File>> fileNames = new HashMap<String, ArrayList<File>>();
		for (File f: files) {
			String name = f.getName();
			Matcher mat = JUST_SRR.matcher(name);
			if (mat.matches()) {
				name = mat.group(1);
				ArrayList<File> al ;
				if (fileNames.containsKey(name)) al = fileNames.get(name);
				else {
					al = new ArrayList<File>();
					fileNames.put(name, al);
				}
				al.add(f);
			}
		}
		
		//for each srr look to see if it is already present
		System.out.print("\t");
		for (int i=0; i< srrToProcess.size(); i++){
			String srr = srrToProcess.get(i);
			//System.out.println("Checking "+srr);
			if (fileNames.containsKey(srr)){
				//System.out.println("Found it, any gz's?");
				ArrayList<File> al = fileNames.get(srr);
				//look for gz file
				for (File f: al){
					File[] gz = IO.extractFiles(f, ".gz");
					if (gz !=null && gz.length !=0 && gz[0].length() > 1000){
						System.out.print(srr+" ");
						srrToProcess.remove(i);
						i--;
						break;
					}
				}
			}
		}
		System.out.println();
	}

	//methods
	public boolean convertSRA(File sra){
		try {
			String[] cmd = {
					fastqDump.getCanonicalPath(),
					"--split-3",
					"--gzip",
					"--origfmt",
					"--outdir",
					sra.getParentFile().getCanonicalPath(),
					sra.getCanonicalPath(),
			};

			System.out.println("\tconverting "+sra);

			String[] output = IO.executeCommandLineReturnAll(cmd);
			//look for error
			for (String line: output){
				if (line.contains("error")) {
					System.err.println("\nProblem executing fastq conversion for "+sra);
					for (String l: output) System.err.println(l);
					return false;
				}
			}

			return true;
		} catch (IOException e) {
			System.err.println("\nProblem converting to fastq for -> "+sra);
			e.printStackTrace();
			return false;
		}



	}

	/**Fetches file, return null if something bad happened.*/
	public File fetchRun(String srr){
		String path = buildSRRPath(srr);
		File file = new File(saveDirectory, srr+".sra");

		FileOutputStream out = null;
		try {
			out = new FileOutputStream(file);
			System.out.println("\tfetching   "+path);
			boolean fetched = ftp.retrieveFile(path,out);
			if (fetched == false) return null;
			return file;
		} catch (FTPConnectionClosedException e){
			//attempt to reconnect connect, this will exit if it fails
			System.out.println("\twarning    connection dropped attempting reconnect");
			connect();
			//if it gets here then it is connected
			return fetchRun(srr);
		} catch (Exception e){
			System.err.println("Error fetching "+srr);
			e.printStackTrace();
			return null;
		} finally {
			if (out != null)
				try {
					out.close();
				} catch (IOException e) {}
		}
	}

	public boolean fetchConvertDeleteSRR(String srr){
		srrToProcess.remove(srr);
		//fetch it
		File srrFile = fetchRun(srr);
		if (srrFile == null) {
			srrFailedToProcess.add(srr);
			return false;
		}

		//make own folder?
		File dir = null;
		if (cmdFile != null){
			dir = new File (saveDirectory, srr);
			dir.mkdir();
			File cmdCopy = new File (dir, cmdFile.getName());
			IO.copyViaFileChannel(cmdFile, cmdCopy);
			File moved = new File (dir, srrFile.getName());
			srrFile.renameTo(moved);
			srrFile = moved;
		}

		//convert it to fastq.gz
		if (convertSRA(srrFile) == false) {
			srrFailedToProcess.add(srr);
			return false;
		}

		//deletes only if both passed
		srrFile.delete();
		srrProcessed.add(srr);

		//launch tomato
		if (cmdFile != null){
			File b = new File(dir,"b");
			try {
				b.createNewFile();
				//System.out.println("\tlaunching\t"+dir);
			} catch (IOException e) {
				System.err.println("\nFailed to make a b file in "+dir);
				e.printStackTrace();
			}
		}

		return true;
	}

	public String buildSRRPath(String srr){
		//build path
		StringBuilder path = new StringBuilder(srrPath);
		//add first six
		path.append(srr.substring(0, 6));
		path.append("/");
		//append dir
		path.append(srr);
		path.append("/");
		//append file
		path.append(srr);
		path.append(".sra");
		return path.toString();
	}

	public String[] fetchProjectSRRNames( String srpName){
		String projectPath = buildSRPPath(srpName);
		System.out.print("\t"+projectPath);
		ArrayList<String> names = fetchSRRNames(projectPath);
		System.out.println("\t"+Misc.stringArrayListToString(names, " "));
		srrToProcess.addAll(names);
		return Misc.stringArrayListToStringArray(names);
	}

	public String buildSRPPath(String srpName){
		//build path
		StringBuilder path = new StringBuilder(srpPath);
		//add first six
		path.append(srpName.substring(0, 6));
		path.append("/");
		//append dir
		path.append(srpName);
		path.append("/");
		return path.toString();
	}

	public void connect(){
		try {
			//any old connections?
			if (ftp != null) ftp.disconnect();
			//make new
			ftp = new FTPClient();
			ftp.connect(ipAddress);
			ftp.enterLocalPassiveMode();
			//connect
			if(FTPReply.isPositiveCompletion(ftp.getReplyCode()) == false) throw new Exception ("Failed to connect, aborting.");
			//login
			if (ftp.login("anonymous",email) == false) throw new Exception("Failed to login, aborting");
			//set binary
			if (ftp.setFileType(FTPClient.BINARY_FILE_TYPE) == false) throw new Exception("Failed to set binary file type, aborting");
			ftp.setControlKeepAliveTimeout(15);
			
		} catch (Exception ex) {
			System.out.println("\nFailed to connect! Aborting.");
			printProcessed();
			System.out.println();
			IO.deleteFiles(saveDirectory, ".sra");
			ex.printStackTrace();
			System.exit(1);
		}
	}

	public void printProcessed(){
		System.out.println("\nProcessing results:");
		System.out.println("\tProcessed: "+Misc.stringArrayListToString(srrProcessed, ","));
		System.out.println("\tFailed processing: "+Misc.stringArrayListToString(srrFailedToProcess, ","));
		System.out.println("\tWaiting processing: "+Misc.stringArrayListToString(srrToProcess, ","));
	}


	/**Fetches recursively all xxx.sra file names*/
	public ArrayList<String> fetchSRRNames(String path){
		//System.out.println("FetchSRRNames for \t"+path);
		ArrayList<String> names = new ArrayList<String>();
		try {

			//fetch files
			FTPFile[] files = ftp.listFiles(path);
			
			if (files != null) {
				for (FTPFile f: files){
					String name = f.getName();
					//System.out.println("Testing\t"+name);
					if (name.startsWith("SR") == false || name.endsWith(".md5")) {
						//System.out.println("\tSkipping "+name);
						continue;
					}
					if (name.endsWith(".sra")) {
						//System.out.println("xxxxxxx adding "+name);
						names.add(name.substring(0, name.length()-4));
					}
					else {
						String newPath = path + name + "/";
						//System.out.println ("\t\t\tNewPath\t"+ newPath);
						ArrayList<String> newNames = fetchSRRNames(newPath);
						names.addAll(newNames);
						//System.out.println("yyyyyyy adding "+newNames);
						
					}
					//System.out.println("CurrListOfNames "+names);
				}
			}
			//System.out.println("returning "+names);
			return names;
		} catch (IOException e) {
			System.err.println("\nProblem fetching srr file names from path "+path);
			e.printStackTrace();
		}
		return names;
	}


	public void disconnect(){
		if (ftp != null)
			try {
				ftp.disconnect();
			} catch (IOException e) {}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SRAProcessor(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+USeqUtilities.stringArrayToString(args, " ")+"\n");
		String[] names = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'n': names =  args[++i].split(","); break;
					case 'f': fastqDump = new File (args[++i]); break;
					case 's': saveDirectory = new File (args[++i]); break;
					case 'c': cmdFile = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: USeqUtilities.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					USeqUtilities.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check params
		if (names == null || names.length == 0) USeqUtilities.printErrAndExit("\nError: please enter SRP or SRR names to fetch.\n");
		for (String n: names){
			if (n.contains(".")) Misc.printErrAndExit("\nError: SRP and SRR names cannot contain a . fix -> "+n);
			if (n.startsWith("SRP")) srpToProcess.add(n);
			else if (n.startsWith("SRR")) srrToProcess.add(n);
			else Misc.printErrAndExit("\nError: SRP and SRR names must start as such, fix -> "+n);
		}
		if (fastqDump.exists() == false || fastqDump.canExecute() == false) USeqUtilities.printErrAndExit("\nError: Cannot find or execute your fastq-dump application.\n");
		if (saveDirectory.exists() == false) {
			if (saveDirectory.mkdirs() == false) Misc.printErrAndExit("\nError: Failed to make your save directory, aborting.\n");
		}



	}

	public static void printDocs(){


		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               SRA Processor: August 2012                         **\n" +
				"**************************************************************************************\n" +
				"Fetchs SRA files from the Sequence Read Archive and converts them to gzipped fastq.\n" +
				"Use in conjunction with Tomato to align these on the ember cluster. Be sure the SRA\n" +
				"archives you want are really in fastq format.\n\n" +

				"Required Parameters:\n"+
				"-n Names of SRRs (runs) or SRPs (projects) to fetch, comma delimited, no spaces.\n" +
				"       (e.g. SRR016669 or SRP000401).\n"+
				"-f Fastq-dump executable, full path, from the SRA Toolkit, download from\n" +
				"       http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software\n"+
				"-s Save directory, full path.\n"+

				"\nOptional Parameters:\n"+
				"-c Full path to a cmd.txt file to copy into converted SRA folders. If the save\n" +
				"       directory is scanned by tomato, a tomato job is then launched,\n" +
				"       see http://bioserver.hci.utah.edu/BioInfo/index.php/Software:Tomato\n" +

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/SRAProcessor -n SRP000401 /\n" +
				"      -s /tomato/job/Nix/SRP000401/ -f ~/sratoolkit.2.1.8-centos_linux64/fastq-dump\n" +
				"      -c /tomato/job/Nix/SRP000401/cmd.txt \n" +

				"\n"+

		"**************************************************************************************\n");

	}


}
