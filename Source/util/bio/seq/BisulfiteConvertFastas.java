package util.bio.seq;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.*;

/**
 * Converts fasta file(s) into serialized boolean[]s where every base g or c is true all others false.
 */
public class BisulfiteConvertFastas {
	//fields
	private File directory;
	private Pattern lowerCaseC = Pattern.compile("c");
	private Pattern upperCaseC = Pattern.compile("C");
	
	public BisulfiteConvertFastas(String[] args){
		
		processArgs(args);
		System.out.print("Converting");
		
		File[] files = IO.extractFiles(directory, "fasta.gz");
		
		if (files == null || files.length == 0) {
			files = IO.extractFiles(directory, "fa.gz");
		}
		
		for (int i=0; i< files.length; i++){
			System.out.print(".");
			try{
				String name = Misc.removeExtension(files[i].getName())+".bs.fasta";
				File outFile = new File (files[i].getParentFile(),name);
				BufferedReader in = IO.fetchBufferedReader(files[i]);
				PrintWriter out = new PrintWriter( new FileWriter(outFile));
				String line;
				Matcher mat;
				while ((line = in.readLine())!=null){
					if (line.length() == 0 || line.startsWith(">")) out.println(line);
					else {
						line = lowerCaseC.matcher(line).replaceAll("t");
						line = upperCaseC.matcher(line).replaceAll("T");
						out.println(line);
					}
				}
				
				
				in.close();
				out.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		System.out.println("\nDone!");
		
	}

	/**This method will process each argument and assign any new varibles*/
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
			System.out.println("\nCannot find your directory!\n");
			System.exit(0);
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Bisulfite Convert Fastas: Dec 2008                        **\n" +
				"**************************************************************************************\n" +
				"Converts all the c/C's to t/T's in fasta file(s) maintaining case.\n\n"+
				
				"Required Parameters:\n"+
				"-f Full path text for the xxx.fasta file or directory containing such.\n" +
	
				"\n" +
				"Example: java -Xmx2000M -jar pathTo/Apps/BisulfiteConvertFastas -f /affy/Fastas/\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new BisulfiteConvertFastas(args);
	}

}

