package trans.misc;
import util.gen.*;
import java.io.*;

/**Converts a file format to something usable by Excel from the joinRegionsWithScores.sh script.*/
public class ScoredRegionsFormatFixer {

	public static void main (String[] args){
		if (args.length ==0) Misc.printExit("\nEnter the full path file text to convert format for excel sorting, from the joinRegionsWithScores.sh script.\n");
		
		//fetch file
		File txtFile = new File(args[0]);
		if (txtFile.canRead()==false) Misc.printExit("\nCannot read or find your 'xxx.txt' file?\n");
		
		try{
			//read in txt file and write out chromosome specific gr files
			BufferedReader in = new BufferedReader (new FileReader(txtFile));
			PrintWriter out = new PrintWriter(new FileWriter(txtFile+".etxt"));
			String line;
			String chrom;
			//find text
			while ((line = in.readLine()) !=null){
				if (line.startsWith("chr")){
					chrom = line.trim();
					while ((line = in.readLine()) !=null){
						//is it empty?
						if (line.trim().length() ==0) break;
						//nope it contains a region and score
						String[] tokens = line.split("\\s+");
						out.println(chrom+"\t"+tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]);
					}
				}
			}
			in.close();
			out.close();
	
		}catch (Exception e){
			e.printStackTrace();
		}
		System.out.println("Done!");
		
	}
	
}
