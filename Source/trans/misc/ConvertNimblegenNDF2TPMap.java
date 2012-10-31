package trans.misc;
import java.io.*;
import java.util.regex.*;
import trans.tpmap.TPMapLine;
import util.gen.*;
import java.util.*;

/** 
 * Converts a Nimblegen NDF txt file to a tpmap.!  
 */
public class ConvertNimblegenNDF2TPMap {
	//fields
	private File[] ndfFiles;
	private TPMapLine[] tpmapLines;
	private int fileIndex;
	private int maxX = 0;
	private int maxY = 0;
	private int numberSkippedControls = 0;
	private boolean parseAll = false;
	
	//constructor
	public ConvertNimblegenNDF2TPMap (String[] args){
		processArgs(args);
		
		//for each file
		for (fileIndex =0; fileIndex< ndfFiles.length; fileIndex++){
			System.out.println("\nParsing "+ndfFiles[fileIndex]+"...");
			
			//parse it
			if (parseNDF() == false)	Misc.printExit("\nError parsing ndf file!\n");
			System.out.println("\n\t"+maxY+" Rows\t"+maxX+" Columns");
			System.out.println("\tParsed "+tpmapLines.length+" chromosome data lines.");
			System.out.println("\tSkipped "+numberSkippedControls+" control data lines.");
			
			//save tpmap
			System.out.println("\tWriting tpmap file...");
			saveTPMap();
		}
		System.out.println("\nDone!\n");
	}
	
	public boolean parseNDF(){
		ArrayList tpmapLinesAL = new ArrayList();
		String line = null;
		try{
			BufferedReader in = new BufferedReader (new FileReader( ndfFiles[fileIndex]));
			while ((line = in.readLine()) != null){
				if (line.length()== 0 || line.startsWith("PROBE")) continue;
				String[] tokens = line.split("\\t");
				//parse x, y and possibly set max values
				int pmX = Integer.parseInt(tokens[10]);
				int pmY = Integer.parseInt(tokens[9]);
				if (maxX< pmX) maxX = pmX;
				if (maxY< pmY) maxY = pmY;
				//is it a chromosome line?
				if (parseAll || tokens[4].startsWith("CHR") || tokens[4].startsWith("chr")  ) {
					TPMapLine tpmapLine = new TPMapLine(tokens, true);
					if (parseAll == false && tpmapLine.getStart() == -1){
						throw (new Exception ("Could not parse chr start from ndf line!"));
					}
					tpmapLinesAL.add(tpmapLine);
				}
				else {
					numberSkippedControls++;
					continue;
				}
				
			}
			in.close();
			//convert AL to TPMap[]
			tpmapLines = new TPMapLine[tpmapLinesAL.size()];
			tpmapLinesAL.toArray(tpmapLines);
			Arrays.sort(tpmapLines);
			return true;
		} catch (Exception e){
			System.err.println("\nError: problem parsing NDF File -> "+ndfFiles[fileIndex]+"\n\tLine -> "+line);
			e.printStackTrace();
			return false;
		}
	}
	
	public boolean saveTPMap(){
		try {
			String fileName = Misc.removeExtension(ndfFiles[fileIndex].getName()) +".tpmap";
			File tpmapFile = new File (ndfFiles[fileIndex].getParentFile(), fileName);
			PrintWriter outRes = new PrintWriter(new FileWriter(tpmapFile));
			for (int i=0; i<tpmapLines.length; i++){
				outRes.println(tpmapLines[i].getLine());
			}
			outRes.close();
			return true;
		} catch (Exception e){
			System.err.println("\nError: Problem trying to save a tpmap for "+ndfFiles[fileIndex]);
			e.printStackTrace();
			return false;
		}
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': ndfFiles = IO.extractFiles(new File(args[i+1]), ".ndf");  i++; break;
					case 'a': parseAll = true; break; 
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check agilentFiles
		if (ndfFiles == null || ndfFiles.length ==0) Misc.printExit("\nCannot find your xxx.ndf Nimblegen mapping file(s)?\n");
		
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Convert Nimblegen NDF 2 TPMap: June 2008                 **\n" +
				"**************************************************************************************\n" +
				"Converts a Nimblegen NDF text file to a tpmap.\n\n" +
				
				"-f Full path file text or directory containing text xxx.ndf file(s).\n" +
				"-a Parse all lines, default is to look for linese beginning with 'chr' or 'CHR'\n\n" +
				
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/ConvertNimblegenNDF2TPMap -f /badData/\n\n" +
				
		"**************************************************************************************\n");		
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new ConvertNimblegenNDF2TPMap(args);
	}

}
