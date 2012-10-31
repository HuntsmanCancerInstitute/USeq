package trans.cel;
import util.gen.*;
import java.util.*;
import java.io.*;

public class RotateCelFile {
	
	public static void main (String[] args){
		if (args.length == 0) Misc.printExit("\nEnter the full path file text for a text xxx.cel file and the number of rows on the chip.\n");
		int num = Integer.parseInt(args[1]);
		File file = new File (args[0]);
		rotateCelFile(file, num);
	}
	
	/**Rotates a text version cel file 90 degrees clockwise.*/
	public static void rotateCelFile (File celFile, int numRows){
		try{ 
			//first read in cel file converting coordinates and making an ArrayList of CelLines
			int widthMinusOne = numRows -1;
			BufferedReader in = new BufferedReader(new FileReader(celFile));
			String rotatedFile = "rotated_"+celFile.getName();
			PrintWriter out = new PrintWriter (new FileWriter( new File(celFile.getParentFile(), rotatedFile) ));
			
			//print header
			boolean inHeader = true;
			String line;
			int counter = 0;
			while (inHeader) {
				line = in.readLine();
				if (line.indexOf("CellHeader") != -1) inHeader = false;
				out.println(line);
			}
			
			//parse xy coords and make array for sorting
			ArrayList al = new ArrayList (numRows * numRows / 2);
			String[] tokens;
			while ((line = in.readLine()) !=null){
				counter++;
				tokens = line.trim().split("\\s+");
				//check that a data line was found and not stop of file or junk
				if (tokens.length!=5){
					//System.out.println ("\tProcessed "+counter+" cel file data lines.");
					break;
				}
				//add to intensity array
				int x = Integer.parseInt(tokens[0]);
				int y = Integer.parseInt(tokens[1]);
				int transX = widthMinusOne - y;
				int transY = x;
				tokens[0] = transX+"";
				tokens[1] = transY+"";
				//make line
				String celLine = Misc.stringArrayToString(tokens,"\t");
				//save in array
				al.add(new CelLine(transX, transY, celLine));
			}
			in.close();
			
			//convert arraylist to array
			CelLine[] lines = new CelLine[al.size()];
			al.toArray(lines);
			Arrays.sort(lines);
			
			//print lines
			for (int i=0; i< lines.length; i++) out.println(lines[i].getLine());
			
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
}
