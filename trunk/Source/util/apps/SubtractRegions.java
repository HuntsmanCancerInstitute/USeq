package util.apps;
import java.io.*;

import util.gen.*;
import util.bio.annotation.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**Intersects two lists of regions subtracting the first from the second.*/
public class SubtractRegions {

	private HashMap maskerCoor;
	private File masker;
	private File directory;
	private File[] toMask;

	public SubtractRegions(String[] args){
		//process args
		processArgs(args);
		
		//load regions to use in subtraction
		Coordinate[] masks = Coordinate.parseFile( masker, 0, 0);
		Arrays.sort(masks);
		maskerCoor = Coordinate.splitByChromosome(masks);

		//for each file to be reduced/ masked
		for (int i=0; i< toMask.length; i++){
			Coordinate[] coor = Coordinate.parseFile( toMask[i], 0, 0);
			if (directory.isDirectory() == false) directory = directory.getParentFile();
			File subtracted = new File (directory, Misc.removeExtension(toMask[i].getName()) + "_Sub.bed");
			mask (coor, subtracted);
		}
	}

	public void mask (Coordinate[] coor, File toSave){		
		try {
			PrintWriter out = new PrintWriter (new FileWriter (toSave));
			Arrays.sort(coor);
			HashMap chrCoor = Coordinate.splitByChromosome(coor);
			//for each chromosome
			Iterator it = chrCoor.keySet().iterator();
			Coordinate[] mask;
			Coordinate[] forReduction;
			while (it.hasNext()){
				String chromosome = (String) it.next();
				mask = (Coordinate[])maskerCoor.get(chromosome);
				forReduction = (Coordinate[])chrCoor.get(chromosome);
				//convert to arraylist
				ArrayList forReductionAL = new ArrayList(forReduction.length);
				for (int i=0; i< forReduction.length; i++) forReductionAL.add(forReduction[i]);
				
				//for each region to be reduced
				while (forReductionAL.size() !=0){
					Coordinate region = (Coordinate) forReductionAL.remove(0);
					if (mask !=null){
						//for each masked region
						for (int j=0; j< mask.length; j++){
							Coordinate[] pair = Coordinate.subtract (region, mask[j]);
							//completely covered?
							if (pair == null) {
								region = null;
								break;
							}
							//trimmed
							else if (pair.length == 1){
								region = pair[0];
							}
							//center punched
							else {
								//add second to arraylist
								forReductionAL.add(pair[1]);
								region = pair[0];
							}
						}
					}
					if (region != null) {
						out.println(region);
						//System.out.println(region);
					}
				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SubtractRegions(args);
	}		

	/**This method will process each argument and assign new varibles*/
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
					case 'm': masker = new File (args[++i]); break;
					case 'd': directory = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (masker == null || masker.canRead()==false){
			Misc.printErrAndExit("\nPlease enter a file to use in masking.\n");
		}
		if (directory == null || directory.isDirectory() == false){
			Misc.printErrAndExit("\nPlease enter a directory containing files to mask.\n");
		}
		toMask = IO.extractFiles(directory);
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Subtract Regions: May 2009                           **\n" +
				"**************************************************************************************\n" +
				"Removes regions and parts there of that intersect the masking region file.  Provide\n" +
				"tab delimited bed files (chr start stop ...). Assumes interbase coordinates.\n" +

				"\nOptions:\n"+
				"-m Bed file to use in subtracting/ masking.\n"+
				"-d Directory containing bed files to mask.\n"+

				"\nExample: java -Xmx4000M -jar pathTo/Apps/SubtractRegions -d /Anno/TilingDesign/\n" +
				"       -m /Anno/repeatMaskerHg18.bed\n\n" +

		"************************************************************************************\n");
	}
}
