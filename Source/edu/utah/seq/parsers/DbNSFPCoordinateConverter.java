package edu.utah.seq.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class DbNSFPCoordinateConverter {

	//fields
	private File[] chrFiles;
	private File saveDirectory;
	private HashMap<String, Gzipper> chrWriters = new HashMap<String, Gzipper>();
	private long linesParsed = 0;
	private String header = null;
	private int columnNumber = 0;
	
	public DbNSFPCoordinateConverter(String[] args) {
		processArgs(args);

		try {

			System.out.println("Parsing:");
			for (File toParse: chrFiles) {
				System.out.print("\t" +toParse.getName()+"\t");
				long num = parse(toParse);
				System.out.println(num);
				linesParsed+= num;
			}
			System.out.println("\tTotal\t"+linesParsed);

			//close the zippers
			for (Gzipper g: chrWriters.values()) g.close();
			
			//load each, sort and write out composite
			System.out.println("\nSorting:");
			sort();
			
			System.out.println("\nDone!");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void sort() throws Exception{
		PrintWriter out = new PrintWriter(new FileWriter (new File(saveDirectory, "corrSwappedSortedDbNSFP.txt")));
		out.println(header);
		
		//for each chrom file
		for (Gzipper g: chrWriters.values()) {
			BufferedReader in = IO.fetchBufferedReader(g.getGzipFile());
			System.out.println("\t" +g.getGzipFile().getName()+"\t");
			String line;
			String[] columns;
			ArrayList<StringSort> al = new ArrayList<StringSort>();
			
			//for each line, load into memory
			while ((line = in.readLine())!= null){
				columns = Misc.TAB.split(line);
				int position = Integer.parseInt(columns[1]);
				al.add( new StringSort(line, position));
			}
			in.close();
			
			//sort it
			StringSort[] ss = new StringSort[al.size()];
			al.toArray(ss);
			Arrays.sort(ss);
			
			//write it out
			for (StringSort s: ss) out.println(s.line);
		}
		//close final writer
		out.close();
	}
	
	private class StringSort implements Comparable<StringSort> {
		String line;
		int position;
		public StringSort(String line, int position){
			this.line = line;
			this.position = position;
		}
		public int compareTo(StringSort o) {
			if (position > o.position) return 1;
			if (position < o.position) return -1;
			return 0;
		}
	}

	private long parse(File toParse) throws IOException {
		long numParsed = 0;
		String[] columns;
		String line;
		BufferedReader in;
		String currChrom = "";
		Gzipper out = null;

		in = IO.fetchBufferedReader(toParse);
		while ((line=in.readLine())!= null){

			//header line?
			if (line.startsWith("#")){
				if (header == null) header = line;
			}
			else{
				columns = Misc.TAB.split(line);
				//check number columns
				if (columnNumber ==0) columnNumber = columns.length;
				else if (columnNumber != columns.length) Misc.printErrAndExit("\nError, incorrect number columns in \n"+ line);
				//replace columns 0 and 1 with the hg19 coordinates
				//#chr	pos(1-based)	ref	alt	aaref	aaalt	rs_dbSNP150	hg19_chr	hg19_pos(1-based)	hg18_chr	hg18_pos(1-based
				
				//watch out for no converted coordinates
				if (columns[7].equals(".")) continue;
				
				//swap
				columns[0] = columns[7];
				columns[1] = columns[8];
				
				//fetch writer
				if (columns[0].equals(currChrom) == false){
					currChrom = columns[0];
					out = chrWriters.get(columns[0]);
					if (out == null) {
						File f = new File(saveDirectory, "temp_"+columns[0]+".gz");
						f.deleteOnExit();
						out = new Gzipper(f);
						chrWriters.put(columns[0], out);
					}
				}
				out.println(columns, "\t");
				numParsed++;
			}
		}
		in.close();


		return numParsed;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new DbNSFPCoordinateConverter(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();

		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': chrFiles = IO.extractOnlyFiles(new File(args[++i])); break;
					case 's': saveDirectory = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (chrFiles == null || chrFiles.length ==0) Misc.printErrAndExit("\nCan't find your files to parse?\n");
		if (saveDirectory == null) Misc.printErrAndExit("\nCan't find your save directory? "+saveDirectory);
		saveDirectory.mkdirs();
	}
		
		public static void printDocs(){
			System.out.println("\n" +
					"**************************************************************************************\n" +
					"**                         DbNSFP Coordinate Converter: Dec 2017                    **\n" +
					"**************************************************************************************\n" +
					"Walks a directory of dbNSFP files swapping the B38 coordinates with the B37, splits\n" +
					"by chromosome, sorts, and writes out the final composite. "+

					"\nOptions:\n"+
					"-d Path to a directory of dbNSFP files to parse.\n" +
					"-s Path to a directory for saving the results.\n\n"+

				
					"\nExample: java -Xmx20G -jar pathToUSeq/Apps/DbNSFPChrSplitter -d DbNSFP3.5a \n" +
					"     -s B37_DbNSFP3.5a \n\n" +

					"**************************************************************************************\n");

		}

}
