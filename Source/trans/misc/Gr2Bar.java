package trans.misc;
import util.gen.*;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import trans.misc.*;
import trans.roc.Sgr;

public class Gr2Bar {
	
	private File[] files;
	private String genomeVersion;
	
	public Gr2Bar(String[] args) {
		//check for args 
		processArgs(args);
		
		for (int x=0; x< files.length; x++){
			//attempt to extract a chrom text
			String chrom = Util.parseChromosomeName(files[x].getName().replaceAll("Chr", "chr"));
			if (chrom == null) System.out.println("\tCould not parse a chromosome text! Skipping -> "+files[x].getName());
			else {
				//does it exist?
				File bar = new File (files[x].getParentFile(), Misc.replaceEnd(files[x].getName(), ".gr.zip", ".bar"));
				if (bar.exists()) {
					System.out.println("\t"+bar+" already exists! Skipping...");
					continue;
				}
				System.out.println("\tProcessing -> "+files[x]+"  Parsed chrom -> "+chrom +"\n\tLoading...");
				
				//make arrays to hold gr values
				ArrayList positions = new ArrayList(10000);
				ArrayList values = new ArrayList(10000);
				
				//load file
				try{
					String line;
					BufferedReader in = IO.fetchReaderOnZippedFile(files[x]);
					while ((line=in.readLine())!=null){
						if (line.trim().length()==0) continue;
						String[] tokens = line.split("\\s");
						Integer position = new Integer(tokens[0]);
						Float value = new Float(tokens[1]);
						positions.add(position);
						values.add(value);
					}
					in.close();
				} catch (Exception e){
					e.printStackTrace();
				}
				
				//convert ArrayList
				int[] intPositions = Num.arrayListOfIntegerToInts(positions);
				positions = null;
				float[] floatValues = Num.arrayListOfFloatToArray(values);
				values = null;
				
				//write xxx.bar file?
				if (intPositions != null && intPositions.length != 0){
					System.out.println("\tWriting...");
					Util.writeSimpleBarFile(chrom, genomeVersion, ".", intPositions, floatValues, bar);
				}
				else System.out.println("\tProblem making bar file for "+files[x]);
			}
		}
		
	}
	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}	
		new Gr2Bar(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String filesString = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': filesString = args[i+1]; i++; break;
					case 'v': genomeVersion = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//parse files and genome version
		if (filesString == null ) Misc.printExit("\nError: cannot find your xxx.sgr.zip file(s)?");
		files = IO.fetchFilesRecursively(new File(filesString), ".gr.zip");
		if (genomeVersion == null) Misc.printExit("\nError: you must supply a genome version. Goto http://genome.ucsc.edu/cgi-" +
				"bin/hgGateway load your organism to find the associated genome version.\n");
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 Gr2Bar: Nov 2006                                 **\n" +
				"**************************************************************************************\n" +
				"Converts xxx.gr.zip files to chromosome specific bar files.\n\n" +

				"-f The full path directory/file text for your xxx.gr.zip file(s).\n" +
				"-v Genome version (ie H_sapiens_Mar_2006), get from UCSC Browser,\n" +
				"      http://genome.ucsc.edu/FAQ/FAQreleases\n" +
				
				"\nExample: java -Xmx1500M -jar pathTo/T2/Apps/Gr2Bar -f /affy/GrFiles/ -v hg17 \n\n" +
				
		"**************************************************************************************\n");		
	}
	
}
