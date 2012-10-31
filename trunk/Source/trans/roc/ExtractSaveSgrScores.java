package trans.roc;
import java.io.*;
import java.util.*;

import util.gen.*;
/**
 * Application for extracting and saving the scores from a .sgr file as a serialized float[].
 */
public class ExtractSaveSgrScores {
	//fields
	private File[] sgrFiles;
	
	public ExtractSaveSgrScores(String[] args){
		if (args.length==0) {
			System.out.println("Enter a full path for an Sgr file.");
			System.exit(1);
		}
		//parse params
		File sgr1 = new File(args[0]);
		saveScores(sgr1);
		

	}

	public static void main(String[] args){
		new ExtractSaveSgrScores(args);
	}
	
	public static void saveScores(File sgrFile){
		String line;
		BufferedReader in = null;
		Sgr sgr = null;
		ArrayList al = new ArrayList();
		try{
			in = new BufferedReader(new FileReader(sgrFile));
			while ((line=in.readLine())!=null){
				if (line.trim().length()==0) continue;
				//make Sgr to hold line info
				sgr = new Sgr(line);
				al.add(new Float(new Double(sgr.getScore()).floatValue()));
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		float[] f = Num.arrayListOfFloatToArray(al);
		IO.saveObject(new File(sgrFile+".scrs"), f);
	}
}
