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
public class ConvertFastaA2G {
	//fields
	private File directory;
	private File saveDirectory;
	private Pattern lowerCase = Pattern.compile("a");
	private Pattern upperCase = Pattern.compile("A");
	
	public ConvertFastaA2G(String[] args){
		
		processArgs(args);
		System.out.print("\nConverting");
		
		HashMap<String, File> chromFiles = Seq.fetchChromosomeFastaFileHashMap(directory);
		if (chromFiles == null || chromFiles.keySet().size()==0) Misc.printErrAndExit("\nERROR: Cannot find any fasta files in "+directory);
		
		for (String chrom: chromFiles.keySet()){
			System.out.print(" "+chrom);
			try{
				File convertedFile = new File (saveDirectory, chrom+".fasta.gz"); 
				Gzipper out = new Gzipper(convertedFile);
				BufferedReader in = IO.fetchBufferedReader(chromFiles.get(chrom));
				
				String line;
				while ((line = in.readLine())!=null){
					if (line.length() == 0 || line.startsWith(">")) out.println(line);
					else {
						line = lowerCase.matcher(line).replaceAll("g");
						line = upperCase.matcher(line).replaceAll("G");
						out.println(line);
					}
				}
				in.close();
				out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		System.out.println("\n\nDone!");
		
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
					case 'f': directory = new File(args[i+1]); i++; break;
					case 's': saveDirectory = new File(args[i+1]); i++; break;
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
		if (directory==null){
			System.out.println("\nCannot find your fasta directory!\n");
			System.exit(0);
		}
		if (saveDirectory == null) Misc.printErrAndExit("\nERROR: please enter a directory in which to save the converted files.");
		if (saveDirectory.equals(directory)) Misc.printErrAndExit("\nERROR: the save directory cannot be the same as the fasta directory.");
		saveDirectory.mkdirs();
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Convert Fasta A 2 G: Mar 2012                             **\n" +
				"**************************************************************************************\n" +
				"Converts all the a/A's to g/G's in fasta file(s) maintaining case.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path for the fasta file (.fa/.fasta/.gz/.zip OK) or directory containing such.\n" +
				"-s Full path directory to save the converted files.\n" +
	
				"\n" +
				"Example: java -Xmx2G -jar pathTo/Apps/ConvertFastaA2G -f /mm9/Fastas/ -s\n" +
				"      /mm9/AGConvertedFastas/\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new ConvertFastaA2G(args);
	}

}

