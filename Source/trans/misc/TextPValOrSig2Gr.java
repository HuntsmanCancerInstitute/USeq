package trans.misc;
import util.gen.*;
import java.io.*;
import java.util.*;

/**Converts a bed file to multiple gr files.*/
public class TextPValOrSig2Gr {

	public static void main (String[] args){
		if (args.length ==0) Misc.printExit("\nEnter the full path 'xxx_pvalue.txt' or 'xxx_signal.txt' " +
				"file text to convert to multiple chromosome specific 'xxx.gr' files.\n");
		
		//fetch file
		File txtFile = new File(args[0]);
		if (txtFile.canRead()==false || args[0].endsWith(".txt") == false) Misc.printExit("\nCannot read or find your 'xxx.txt' file?\n");
		
		try{
			System.out.println("Loading txt file...");
			//read in txt file and write out chromosome specific gr files
			BufferedReader in = new BufferedReader (new FileReader(txtFile));
			String line;
			//find text
			System.out.print("Writing ");
			while ((line = in.readLine()) !=null){
				if (line.startsWith("# Name")){
					String[] tokens = line.split("\\s+");
					String chrom = tokens[2];
					System.out.print(chrom+" ");
					//make printwriter
					String fileName = Misc.replaceEnd(args[0], ".txt", chrom+".gr");
					PrintWriter out = new PrintWriter(new FileWriter(fileName));
					//skip next # comment line 
					in.readLine();
					//print each line till hit next comment line
					while ((line = in.readLine()) !=null){
						//is it a comment?
						if (line.startsWith("#")) {
							out.close();
							break;
						}
						//print if not blank
						else if (line.trim().length()!=0) {
							out.println(line);
						}
					}
					//close last PrintWriter
					out.close();
				}
			}
			in.close();
			System.out.println(" gr files.");
	
		}catch (Exception e){
			e.printStackTrace();
		}
		System.out.println("Done!");
		
	}
	
}
