package trans.misc;
import java.io.*;
import java.util.regex.*;
import util.gen.*;


/** 
 * Converts Nimblegen PAIR txt files to a cela files.!  
 */
public class ConvertNimblegenPAIR2Cela {
	//fields
	private File[] pairFiles;
	private float[][] intensities;
	private int fileIndex;
	private String fileName;
	private int maxX =0;
	private int maxY =0;
	
	//constructor
	public ConvertNimblegenPAIR2Cela (String[] args){
		processArgs(args);
		
		//for each file
		for (fileIndex =0; fileIndex< pairFiles.length; fileIndex++){
			System.out.println("\nParsing "+pairFiles[fileIndex]+"...");
			
			//set maxX and Y?
			if (fileIndex == 0) {
				System.out.print("\tFinding max indexes...");
				if (setMaxValues() == false) System.exit(0);
				System.out.println(" "+maxY+" Rows, "+maxX+" Columns");
			}
			
			//parse it
			if (parsePAIR() == false)	Misc.printExit("\nError parsing pair file!\n");
			
			//save cela
			System.out.println("\tWriting cela file...");
			File celaFile = new File (pairFiles[fileIndex].getParentFile(), fileName+".cela");
			IO.saveObject(celaFile, intensities);
		}
		System.out.println("\nDone!\n");
	}
	
	public boolean setMaxValues(){
		String line = null;
		try{
			BufferedReader in = new BufferedReader (new FileReader( pairFiles[0]));
			//skip header
			while ((line = in.readLine()) != null){
				if (line.startsWith("IMAGE_ID")) break;
			}
			while ((line = in.readLine()) != null){
				String[] tokens = line.split("\\t");
				//parse x, y and possibly set max values
				int pmX = Integer.parseInt(tokens[5]);
				int pmY = Integer.parseInt(tokens[6]);
				if (maxX< pmX) maxX = pmX;
				if (maxY< pmY) maxY = pmY;
			}
			in.close();
			return true;
		} catch (Exception e){
			System.err.println("\nError: problem parsing PAIR File -> "+pairFiles[fileIndex]+"\n\tLine -> "+line);
			e.printStackTrace();
			return false;
		}
	}

	
	public boolean parsePAIR(){
		fileName = Misc.removeExtension(pairFiles[fileIndex].getName());
		intensities = new float[maxX][maxY];
		String line = null;
		try{
			BufferedReader in = new BufferedReader (new FileReader( pairFiles[fileIndex]));
			//skip header
			while ((line = in.readLine()) != null){
				if (line.startsWith("IMAGE_ID")) break;
			}
			//load intensity array
			while ((line = in.readLine()) != null){
				if (line.length()== 0 ) continue;
				String[] tokens = line.split("\\t");
				//parse x, y, subtract 1 to get to zero index based
				int pmX = Integer.parseInt(tokens[5]) -1;
				int pmY = Integer.parseInt(tokens[6]) -1;
				//parse PM value
				float pm = Float.parseFloat(tokens[9]);
				//set
				intensities[pmX][pmY] = pm;
				//System.out.println (intensities[pmX][pmY] +"\t"+line);
			}
			in.close();
			return true;
		} catch (Exception e){
			System.err.println("\nError: problem parsing PAIR File -> "+pairFiles[fileIndex]+"\n\tLine -> "+line);
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
					case 'f': pairFiles = IO.extractFiles(new File(args[i+1]), ".PAIR");  i++; break;
					//case 's': switchXYColumns = true; break; 
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
		if (pairFiles == null || pairFiles.length ==0) Misc.printExit("\nCannot find your xxx.PAIR Nimblegen mapping file(s)?\n");
		
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Convert Nimblegen PAIR 2 Cela: Jan 2007                     **\n" +
				"**************************************************************************************\n" +
				"Converts Nimblegen PAIR text data file(s) to cela.\n\n" +
				
				"-f Full path file text or directory containing text xxx.pair file(s).\n\n" +
				
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/ConvertNimblegenPAIR2Cela -f /badData/\n\n" +
				
		"**************************************************************************************\n");		
	}
	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		else new ConvertNimblegenPAIR2Cela(args);
	}

}
