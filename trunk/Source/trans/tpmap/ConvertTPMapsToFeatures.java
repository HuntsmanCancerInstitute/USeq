package trans.tpmap;

import java.util.regex.*;
import java.io.*;
import util.gen.*;
import java.util.*;

/**
 * Converts TPMaps to MapFeature[] and saves them to disk.
 */
public class ConvertTPMapsToFeatures {
	//fields
	private File[] tpmapFiles;
	
	public ConvertTPMapsToFeatures(String[] args){
		processArgs(args);
		System.out.println("\nLaunching...");
		
		//save array
		for (int x=0; x<tpmapFiles.length; x++){
			MapFeature[] f = convert(tpmapFiles[x]);
			System.out.println("\tSaving MapFeature[]...");
			String name = IO.getFullPathName(tpmapFiles[x]);
			name = name.replaceFirst(".tpmap", ".fa");
			File newFile = new File(name);
			IO.saveObject(newFile, f);
		}
		System.out.println("\nDone!\n");
	}
	
	/**Converts a text bpmap to a MapFeature[]*/
	public static MapFeature[] convert(File tpmapFile){
		MapFeature[] f = null;
		ArrayList al = null;
		try {
			//count lines
			System.out.print("\tCounting # tpmap lines in "+tpmapFile.getName());
			int numLines = (int)IO.countNonBlankLines(tpmapFile);
			System.out.println(" -> "+(numLines-2));
			al = new ArrayList(numLines);
			//read in and make array
			System.out.println("\tMaking MapFeature[]...");
			String line;
			BufferedReader in = new BufferedReader(new FileReader(tpmapFile));
			
			while ((line = in.readLine()) !=null) {
				if (line.startsWith("#") == false) al.add(new MapFeature(line));
			}
			in.close();
			
		} catch (Exception e) {e.printStackTrace();}
		f = new MapFeature[al.size()];
		al.toArray(f);
		return f;
	}
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File dir = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': dir = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nError: problem, unknown option! " + mat.group()+"\n"); System.exit(0);
					}
				}
				catch (Exception e){
					System.out.print("\nError: something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check to see if they entered required params
		if (dir==null || dir.canRead() == false){
			System.out.println("\nError: cannot find your tpmap file or directory! ->"+dir);
			System.exit(0);
		}
		tpmapFiles = IO.extractFiles(dir, ".tpmap");
		if (tpmapFiles.length ==0){
			System.out.println("\nError: didn't find any xxx.tpmap files !\n");
			System.exit(0);
		}
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                   Convert TPMaps To Features: Oct 2005                     **\n" +
				"**************************************************************************************\n" +
				"Converts TPMap files to MapFeature[] arrays and saves them to disk. Tpmap lines:\n" +
				"seq ori chrom start pmX pmY mmX mmY. The mm coords are optional.\n\n"+
				
				"Parameters:\n"+
				"-f Full path text for the 'xxx.tpmap' file or directory containing the same.\n" +
				
				"\n" +
				"Example: java -Xmx500M ConvertTPMapsToFeatures -f /affy/dmel.tpmap\n" +
				"\n" +
				
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
		}
		else new ConvertTPMapsToFeatures(args);
	}
	
	
	
	
	
	
}
