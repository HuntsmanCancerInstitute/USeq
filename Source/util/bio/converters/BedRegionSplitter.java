package util.bio.converters;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**
 * For each region, splits it into xxx bp chunks if necessary. >
 */
public class BedRegionSplitter {
	
	private int chunkSize = 2000;
	private File[] toParse;
	

	public BedRegionSplitter(String[] args) {
		processArgs(args);
		run();
		System.out.println("\nDone!\n");
	}
	
	public void run(){
		//for each file
		System.out.println("File\t#Start\t#End");
		for (int i=0; i< toParse.length; i++){
			if (toParse[i].isDirectory()) continue;
			System.out.print(toParse[i].getName());
			parse(toParse[i]);
		}
	}

	private void parse(File file) {
		try {
			BufferedReader in = IO.fetchBufferedReader(file);
			File p = new File(file.getParentFile(), Misc.removeExtension(file.getName())+"_split.bed");
			Gzipper out = new Gzipper(p);
			int[] startEnd = parseIt(in, out);
			out.close();
			in.close();
			System.out.println("\t"+startEnd[0]+ "\t"+ startEnd[1]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private int[] parseIt(BufferedReader in, Gzipper out) throws IOException {
		String[] t;
		String line;
		int numStartingRegions = 0;
		int numFinalRegions = 0;
		while ((line = in.readLine())!=null){
			t = Misc.TAB.split(line);
			if (line.trim().length() == 0 || line.startsWith("#") || t.length < 3) out.println(line);
			else{
				int start = Integer.parseInt(t[1]);
				int stop = Integer.parseInt(t[2]);
				numStartingRegions++;
				//small enough?
				int size = stop - start;
				if (size <= chunkSize) {
					out.println(line);
					numFinalRegions++;
				}
				else {
					int[][] chunks = Num.chunkRegion(chunkSize, start, stop);
					numFinalRegions+= chunks.length;
					String extra =null;
					if (t.length > 3) {
						StringBuilder x = new StringBuilder();
						for (int i=3; i< t.length; i++){
							x.append("\t");
							x.append(t[i]);
						}
						extra = x.toString();
					}
					String chr = t[0]+"\t";
					for (int i=0; i< chunks.length; i++){
						out.print(chr);
						out.print(chunks[i][0]);
						out.print("\t");
						out.print(chunks[i][1]);
						if (extra != null) out.println(extra);
						else out.println();
					}
				}
			}
		}
		return new int[]{numStartingRegions , numFinalRegions};
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BedRegionSplitter(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': toParse = IO.extractFiles(new File(args[++i])); break;
					case 'c': chunkSize = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (toParse == null || toParse.length ==0) Misc.printExit("\nError: please provide a directory containing bed region files to split.\n");
				
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Bed Region Splitter : June 2017                      **\n" +
				"**************************************************************************************\n" +
				"Regions exceeding the chunk size are split into multiple parts.\n"+

				"\nRequired:\n"+
				"-d Path to a file or directory containing such to chunk.\n"+
								
				"\nOptional:\n" +
				"-c BP chunk size, defaults to 2000.\n" +

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/BedRegionSplitter -d ToSplit/ -c 5000 \n\n"+

		"**************************************************************************************\n");

	}


}
