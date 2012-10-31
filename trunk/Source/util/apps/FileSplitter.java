package util.apps;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import util.gen.IO;

/**
 * Application for splitting a large text file into smaller files given a maximum number of lines.
 * 
 * @author nix
 *
 */
public class FileSplitter {
	//fields
	private File[] files;
	private int numLinesPerFile = 0;
	private int fileNumber = 1;
	private boolean gzip = false;

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
					case 'f': files = IO.extractFiles(new File(args[++i])); break;
					case 'n': numLinesPerFile=Integer.parseInt(args[i+1]); i++; break;
					case 'g': gzip = true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check to see if they entered required params
		if (files==null || files[0].canRead() == false){
			System.out.println("\nCannot find your file!\n");
			System.exit(0);
		}
		if (numLinesPerFile == 0){
			System.out.println("\nEnter a maximum number of lines to place in each split file!\n");
			System.exit(0);
		}

	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          File Splitter: July 2010                                **\n" +
				"**************************************************************************************\n" +
				"Splits a big text file into smaller files given a maximum number of lines.\n\n"+

				"Required Parameters:\n"+
				"-f Full path file text or directory for the text file(s) (.zip/.gz OK).\n" +
				"-n Maximum number of lines to place in each.\n" +
				"-g GZip split files.\n"+

				"\n" +
				"Example: java -Xmx256M -jar pathTo/T2/FileSplitter -f /affy/bpmap.txt -n 50000\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		System.out.println("Launching FileSplitter...");
		new FileSplitter(args);
	}

	public FileSplitter(String[] args){
		processArgs(args);

		String line;
		try{
			for (int i=0; i< files.length; i++){
				fileNumber = 1;
				BufferedReader in = IO.fetchBufferedReader(files[i]);

				if (gzip){
					byte[] cr = "\n".getBytes(); 
					GZIPOutputStream out = getGZipOutputStream (files[i]);
					int counter = 0;
					byte[] b = null;
					while ((line = in.readLine()) !=null) { 
						if (counter < numLinesPerFile) {
							counter++;
						}
						else {
							counter =1;
							out.finish();
							out.close();
							out = getGZipOutputStream (files[i]);
						}
						b=line.getBytes();
						out.write(b);
						out.write(cr);
					}
					//close final writer
					out.finish();
					out.close();
				}
				else {
					PrintWriter out = getPrintWriter(files[i]);
					int counter = 0;
					while ((line = in.readLine()) !=null) { 
						if (counter < numLinesPerFile) {
							counter++;
						}
						else {
							counter =1;
							out.close();
							out = getPrintWriter(files[i]);
						}
						out.println(line);
					}
					//close final PrintWriter
					out.close();
				}
			}
		}
		catch (IOException e){
			e.printStackTrace();
		}
	}

	public PrintWriter getPrintWriter(File file){
		PrintWriter out =null;
		try {
			String name = file.getName();
			name = name.replaceAll(".gz", "");
			name = name.replaceAll(".zip", "");
			out = new PrintWriter(new FileWriter(new File(file.getParentFile(),fileNumber+"_"+name)));
			fileNumber++;
		}
		catch (IOException e){
			e.printStackTrace();
		}
		return out;
	}

	public GZIPOutputStream getGZipOutputStream(File file){
		GZIPOutputStream out =null;
		try {
			String name = file.getName();
			name = name.replaceAll(".gz", "");
			name = name.replaceAll(".zip", "");
			name = name+".gz";
			out = new GZIPOutputStream(new FileOutputStream(new File(file.getParentFile(),fileNumber+"_"+name)));
			fileNumber++;
		}
		catch (IOException e){
			e.printStackTrace();
		}
		return out;
	}
}

