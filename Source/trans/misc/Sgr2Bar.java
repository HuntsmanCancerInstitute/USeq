package trans.misc;
import util.gen.*;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.parsers.*;

public class Sgr2Bar {

	private File[] sgrFiles;

	private String genomeVersion = null;
	private boolean stairStep = false;
	private String strand = ".";
	private HashMap <String, String> tagValues = new HashMap <String, String>();
	private BarParser barParser = new BarParser();


	//constructor
	public Sgr2Bar(String[] args) {
		try {
			//check for args 
			processArgs(args);

			System.out.println("Genome version -> "+genomeVersion);
			System.out.println("Strand -> "+strand);
			System.out.println("Stair Step? -> "+stairStep);
			System.out.println();

			//load tagValues
			if (stairStep) tagValues.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
			else tagValues.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);

			for (int x=0; x< sgrFiles.length; x++){
				System.out.println("\tLoading -> "+sgrFiles[x]);
				GrGraph[] grs = GrGraph.loadSgrFile(sgrFiles[x]);

				//make save directory
				String dirName;

				dirName = Misc.removeExtension(sgrFiles[x].getCanonicalPath());

				File dir = new File (dirName);
				dir.mkdir();
				//print bar files
				System.out.println("\tSaving...");
				for (int i=0; i< grs.length; i++){
					File barFile = new File (dir, grs[i].getChromosome()+".bar");
					barParser.writeBarFile(barFile, grs[i].getChromosome(), genomeVersion, strand.charAt(0), grs[i].getBasePositions(), grs[i].getValues(), tagValues);
				}
			}
			System.out.println("\nDone!\n");
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nFailed to parse file!");
		}
	}



	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}	
		new Sgr2Bar(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File dir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': dir = new File(args[i+1]); i++; break;
					case 'v': genomeVersion = args[i+1]; i++; break;
					case 's': strand = args[++i]; break;
					case 't': stairStep = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (dir == null || dir.canRead() == false) Misc.printExit("\nError: cannot find or read your sgr file/ directory.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(dir,".sgr");
		tot[1] = IO.extractFiles(dir,".sgr.zip");
		tot[2] = IO.extractFiles(dir,".sgr.gz");
		sgrFiles = IO.collapseFileArray(tot);

		if (sgrFiles == null || sgrFiles.length ==0) Misc.printExit("\nError: cannot find your xxx.sgr.zip file(s)");
		if (genomeVersion == null) Misc.printExit("\nError: you must supply a genome version. Goto http://genome.ucsc.edu/cgi-" +
		"bin/hgGateway load your organism to find the associated genome version.\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Sgr2Bar: Jan 2012                                  **\n" +
				"**************************************************************************************\n" +
				"Converts xxx.sgr(.zip) files to chromosome specific bar files.\n\n" +

				"-f The full path directory/file text for your xxx.sgr(.zip or .gz) file(s).\n" +
				"-v Genome version (ie H_sapiens_Mar_2006, M_musculus_Jul_2007), get from UCSC Browser.\n" +
				"-s Strand, defaults to '.', use '+', or '-'\n"+
				"-t Graphs should be viewed as a stair-step, defaults to bar\n\n"+


				"Example: java -Xmx1500M -jar pathTo/Apps/Sgr2Bar -f /affy/sgrFiles/ -s + -t\n" +
				"      -v D_rerio_Jul_2006\n\n" +

		"**************************************************************************************\n");		
	}

}
