package util.bio.seq;
import java.io.*;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.*;

/**
 * Converts fasta and sequence A's to G's for ADAR edited alignment.
 */
public class ConvertFastqA2G {
	//fields
	private File[] fastqFiles;
	private File saveDirectory;
	private Pattern a = Pattern.compile("a", Pattern.CASE_INSENSITIVE);

	
	public ConvertFastqA2G(String[] args){
		
		processArgs(args);
		System.out.print("\nConverting");
		
		
		for (File f: fastqFiles){
			System.out.print(" "+f.getName());
			try{
				File convertedFile = new File (saveDirectory, "converted_"+f.getName()); 
				Gzipper out = new Gzipper(convertedFile);
				BufferedReader in = IO.fetchBufferedReader(f);
				
				String line;
				while ((line = in.readLine())!=null){
					if (line.startsWith("@")) {
						//print header line
						out.println(line);
						//load and modify seq line
						line = in.readLine();
						line = a.matcher(line).replaceAll("G");
						out.println(line);
						//print quality lines
						out.println(in.readLine());
						out.println(in.readLine());
					}
				}
				in.close();
				out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		System.out.println("\n\nDone!\n");
		
	}

	/**This method will process each argument and assign any new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fastqFiles = IO.extractFiles(new File(args[++i])); break;
					case 's': saveDirectory = new File(args[++i]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		//check to see if they entered required params
		if (fastqFiles==null || fastqFiles.length ==0){
			Misc.printErrAndExit("\nError: please provide some fastq files to parse.\n");
		}
		if (saveDirectory == null) {
			saveDirectory = fastqFiles[0].getParentFile();
			if (saveDirectory == null) saveDirectory = new File(".");
		}
		saveDirectory.mkdirs();
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Convert Fastq A 2 G: Mar 2012                             **\n" +
				"**************************************************************************************\n" +
				"Converts all the sequence A's to G's, case insensitive.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path for the fastq file or directory containing such. xxx.gz/.zip OK.\n" +
				"-s Optional, full path directory to save the converted files.\n" +
	
				"\n" +
				"Example: java -Xmx2G -jar pathTo/Apps/ConvertFastqA2G -f /IllData/Fastq/ \n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new ConvertFastqA2G(args);
	}

}

