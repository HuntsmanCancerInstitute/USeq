package trans.tpmap;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Misc;
import trans.anno.*;

/**
 * Filters a TPMap file by regions (tab delimited: chr start stop).  Any oligo not contained within one of the regions is pitched.
 */
public class FilterTPMapByRegions {
	//fields
	private File tpmap;
	private File regions;
	
	public FilterTPMapByRegions(String[] args){
		System.out.println("\nLaunching TPMap filter...");
		processArgs(args);
		
		//load regions
		GenomicRegion[] rgns = GenomicRegion.parseRegions(regions);
		int numRegions = rgns.length;
		
		//read in tpmap, make tpmap line, check against regions, if found print out to new file, if not pitch
		try{
			PrintWriter out = new PrintWriter (new FileWriter(new File (tpmap.getCanonicalPath()+"Filtered")));
			BufferedReader in = new BufferedReader (new FileReader(tpmap));
			String line;
			int droppedLines = 0;
			int savedLines = 0;
			while ((line=in.readLine()) != null){
				boolean noIntersection = true;
				//skip comment lines
				if (line.startsWith("#")) out.println(line);
				else {
					TPMapLine tLine = new TPMapLine(line);
					for (int i=0; i<numRegions; i++){
						if (rgns[i].intersect(tLine)){
							out.println(line);
							savedLines++;
							noIntersection = false;
							break;
						}
					}
					if (noIntersection) droppedLines++;
				}
			}
			
			//print stats
			System.out.println("\n"+savedLines+" TPMap lines intersected the regions.");
			System.out.println(droppedLines+ " were dropped.\n");
			
			
			in.close();
			out.close();
			
		}catch (Exception e){
			e.printStackTrace();
		} 
	}
	
	
	
	
	public static void main (String[] args){
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new FilterTPMapByRegions(args);
	}
	
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
					case 't': tpmap = new File(args[i+1]); i++; break;
					case 'r': regions = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		//check to see if they entered required params
		if (tpmap == null || regions == null || tpmap.canRead()== false || regions.canRead() == false){
			Misc.printExit("\nCannot find or read one of your files?! Check and restart.\n");
		}	
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Filter TPMap By Regions: Jan 2006                         **\n" +
				"**************************************************************************************\n" +
				"FTBR strips any oligos not contained within one or more of the user defined regions.\n\n" +

				"    -t Full path file text for the tpmap file.\n" +
				"    -r Full path file text for the tab delimited regions file (chr start stop \n" +
				"           (inclusive, zero based)).\n\n"+
				
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/FilterTPMapByRegions -t /affy/enc.tpmap\n" +
				"           -r /affy/encRgns.txt\n\n" +
				
		"**************************************************************************************\n");
	}
}
