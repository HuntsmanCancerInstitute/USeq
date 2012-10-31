package trans.roc;
import java.io.*;
import java.util.*;

import util.gen.*;

/**
 * Returns oligo values that overlap particular regions.
 */
public class ParseSgrsForParticularRegions {
	public ParseSgrsForParticularRegions(String[] args){
		if (args.length==0) {
			System.out.println("\nTo extract feature values that overlap particular regions, enter two file names: SgrFile/Dir and a RegionsFile (chrom, start, stop)\n");
			System.exit(1);
		}
		//load RegionsFile
		Positive[] regions = parseRegionFile(new File(args[1]));
		//run thru .sgr files
		File[] sgrFiles = IO.extractFiles(new File(args[0]), "sgr");
		for (int x=0; x<sgrFiles.length; x++){
			String line;
			BufferedReader in = null;
			Sgr sgr = null;
			try{
				in = new BufferedReader(new FileReader(sgrFiles[x]));
				while ((line=in.readLine())!=null){
					if (line.trim().length()==0) continue;
					//make Sgr to hold line info
					sgr = new Sgr(line);
					//load it
					loadRegions(sgr, regions);
					//find whether it is in a region
					/*int regionIndex = inRegion(sgr, regions);
					if (regionIndex != -1){
						//out.println(sgr);
						out.println(regions[regionIndex].getStart()+"\t"+sgr);
					}
					*/
				}
				//print results
				for (int a=0; a<regions.length; a++){
					System.out.println(regions[a].toStringSimple());
					ArrayList al = regions[a].getScores();
					if (al.size() !=0){
						System.out.println(al);
					}
				}
				in.close();

			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}
	
	/**Checks to see if an sgr object is in one of the positives, if not returns -1, 
	 * otherwise returns index number.*/
	public static int inRegion(Sgr sgr, Positive[] pos){
		int numPos = pos.length;
		for (int i=0; i< numPos; i++){
			if (pos[i].matches(sgr)) return i;
		}
		return -1;
	}
	
	/**Checks to see if an sgr object is in  the positives, if so, adds the sgr object to the positive.*/
	public static void loadRegions(Sgr sgr, Positive[] pos){
		int numPos = pos.length;
		for (int i=0; i< numPos; i++){
			if (pos[i].matches(sgr)) {
				ArrayList al = pos[i].getScores();
				al.add(sgr);
			}
		}
	}
	
	public static Positive[] parseRegionFile(File picksFile){
		Positive[] regions =null;
		try{
			BufferedReader in = IO.fetchBufferedReader(picksFile);
			String line;
			String[] tokens;
			ArrayList regionsAL = new ArrayList();
			String chromosome;
			int start;
			int stop;
			int index = 0;
			//chrom, start, stop
			//  0      1      2         
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.equals("") || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length<3) {
					System.err.println("Problem parsing this line, skipping -> "+line);
					continue;
				}
				chromosome = tokens[0];
				start = Integer.parseInt(tokens[1]);
				stop = Integer.parseInt(tokens[2]);
				Positive pos = new Positive(index++, chromosome, start, stop);
				if (tokens.length > 3 && tokens[3].equals("-")) pos.setSenseStrand(false);
				regionsAL.add(pos);
			}
			regions = new Positive[regionsAL.size()];
			regionsAL.toArray(regions);
			
		}catch (IOException e){
			e.printStackTrace();
		}
		return regions;
	}
	
	public static void main(String[] args){
		new ParseSgrsForParticularRegions(args);
	}
	
}
