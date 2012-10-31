package trans.roc;
import java.io.*;
import java.util.*;

/**Assignes values from an sgr file to a collection of {@link RocWindow}s.*/
public class ExtractScoresForRoc {
	public ExtractScoresForRoc(String[] args){
		
		//fetch roc windows
		System.out.println("\n\tReading pos and neg windows...");
		RocWindow[] pos = ParsePatternedWindows.loadTxtRocFile(new File (args[0]));
		RocWindow[] neg = ParsePatternedWindows.loadTxtRocFile(new File (args[1]));
		
		//make HashMap
		HashMap win = new HashMap(pos.length+neg.length);
		for (int i=0; i<pos.length; i++){
			pos[i].setType(1); //for positives
			pos[i].setScore(0);
			win.put(pos[i].getChromosome()+pos[i].getMiddle(), pos[i]);
		}
		for (int i=0; i<neg.length; i++){
			neg[i].setType(-1); //for negatives
			neg[i].setScore(0);
			win.put(neg[i].getChromosome()+neg[i].getMiddle(), neg[i]);
		}
		
		//run thru sgr file and assign score to RocWindow
		System.out.println("\tReading Sgr file and assigning scores...");
		File sgrFile = new File(args[2]);
		try{
			BufferedReader in = new BufferedReader(new FileReader(sgrFile));
			String line;
			Sgr sgr;
			RocWindow roc;
			Object obj;
			while ((line = in.readLine()) != null){
				sgr = new Sgr(line);
				obj = null;
				obj = win.get(sgr.getChromosome()+sgr.getPosition());
				//was it present?
				if (obj != null) {
					roc = (RocWindow)obj;
					//check for score
					if (roc.getScore() != 0){
						System.out.println("\nError: assigning a second score!");
						System.exit(1);
					}
					roc.setScore(sgr.getScore());
				}
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		
		//check if all were scored
		System.out.println("\tWriting Results...");
		try{
			PrintWriter out = new PrintWriter(new FileWriter( new File(sgrFile.getCanonicalPath()+".rScrd")));
			for (int i=0; i<pos.length; i++){
				out.println(pos[i].getScore()+"\t+");
				if (pos[i].getScore() == 0) {
					System.out.println("\nError: score not set!");
					System.exit(0);
				}
			}
			for (int i=0; i<neg.length; i++){
				out.println(neg[i].getScore()+"\t-");
				if (neg[i].getScore() == 0) {
					System.out.println("\nError: score not set!");
					System.exit(0);
				}
			}	
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		System.out.println("\tDone!\n");
		
		
	}
	
	
	
	public static void main(String[] args) {
		if (args.length == 0) System.out.println("\n Pos, Neg, Sgr files\n");
		else new ExtractScoresForRoc(args);
		
	}
	
}
