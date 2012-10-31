package trans.cel;

import util.gen.*;

import java.io.*;
import java.util.regex.*;

/**Saves text cel files as serialized float[][]s for downstream use by TiMAT applications.*/
public class CelFileConverter {

	//fields
	private File[] celFiles;
	private File bunZip2App = null;
	private File saveDirectory = null;
	private boolean rotateMatrix = false;
	private boolean saveConcatenatedArray = false;

	public CelFileConverter (String[] args){
		processArgs(args);
		convert();
		System.out.println("\tDone!");
	}

	public void convert(){
		//write float[][] to disk
		for (int i=0; i<celFiles.length; i++){
			System.out.println("\tLoading "+celFiles[i].getName());
			float[][] f;
			if (bunZip2App == null) f= trans.misc.Util.createVirtualCel(celFiles[i]);
			else {
				String[] lines = IO.bunZip2(bunZip2App, celFiles[i]);
				if (lines == null || lines.length == 0) Misc.printExit("\nError, problem with loading bz2 file?!\n");
				f= trans.misc.Util.createVirtualCel(lines);
			}
			//rotate cel file?
			if (rotateMatrix){
				System.out.println("\t\tRotating...");
				f = Num.rotateClockwise(f);
			}
			//save cela file
			String name;
			if (bunZip2App == null) name = Misc.replaceEnd(celFiles[i].getName(),"cel","cela");
			else name = Misc.replaceEnd(celFiles[i].getName(),"cel.bz2","cela");
			File serCel = new File(saveDirectory, name);
			System.out.println("\tSaving "+serCel);
			IO.saveObject(serCel,f);
			//save celp file?
			if (saveConcatenatedArray){
				float[] exp = Num.concatinate(f);
				File expFile = new File(saveDirectory, name+"Con");
				System.out.println("\tSaving "+expFile);
				IO.saveObject(expFile,exp);}
		}
	}

	//main
	public static void main(String[] args) {
		if (args.length!=0) {
			new CelFileConverter(args);
		}
		else {
			printDocs();
		}
	}

	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File directory = null;
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
					case 'b': bunZip2App = new File(args[i+1]); i++; break;
					case 'r': rotateMatrix = true; break;
					case 'c': saveConcatenatedArray = true; break;
					case 'h': printDocs(); System.exit(0);
					default: {
						System.out.println("\nProblem, unknown option! " + mat.group());
						System.exit(0);
					}
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test+"\n");
					System.exit(0);
				}
			}
		}
		//check for tempDirectory
		if (directory == null || directory.exists()== false){
			System.out.println("\nCannot find your 'xxx.cel' file or directory?! -> "+directory);
			System.exit(0);
		}
		String extension = "cel";
		if (bunZip2App !=null) {
			extension = "cel.bz2";
			if (bunZip2App.exists()==false )Misc.printExit("\nError, cannot find your bunzip2 application -> "+bunZip2App+"\n");
		}
		celFiles = IO.extractFiles(directory,extension);
		//any files
		if (celFiles.length == 0){
			System.out.println("\nCannot find any 'xxx."+extension+"' files?! -> "+directory);
			System.exit(0);
		}
		//check for an alternative save directory
		if (saveDirectory !=null){ 
			if (saveDirectory.isDirectory()== false) Misc.printExit("\nError: cannot find or read your alternative results directory -> "+saveDirectory+"\n");
		}
		else saveDirectory = celFiles[0].getParentFile();

	}



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Cel File Converter: Feb  2008                          **\n" +
				"**************************************************************************************\n" +
				"Converts text version cel files into serialized java float[][]s for use by TiMAT2\n" +
				"applications.\n\n" +

				"Parameters:\n"+
				"-f Full path to a text version 'xxx.cel' file or directory containing the same.\n" + 
				"-b Optional, full path to bunzip2. Convert compressed cel files 'xxx.cel.bz2'\n"+
				"-s Full path to alternative save directory, defaults to cel file parent directory.\n"+
				"-c Save a float[] by concatinating the lines of the float[][]\n"+
				"-r Rotate cel file 90 degrees clockwise.\n\n"+

				"Example: java -Xmx512M -jar pathTo/T2/Apps/CelFileConverter -f /data/cels/ \n\n"+

		"**************************************************************************************\n");
	}	

}
