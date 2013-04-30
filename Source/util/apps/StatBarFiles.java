package util.apps;
import util.gen.*;
import java.io.*;
import edu.utah.seq.parsers.*;

public class StatBarFiles {
	public static void main (String[] args){
		if (args.length ==0 ) Misc.printExit("\nEnter a file or directory containing log2 xxx.bar files to calculate stats on.\n");
	
		BarParser parser = new BarParser();
		
		File[] barFiles = IO.extractFiles(new File(args[0]), "bar");
		
		for (int i=0; i< barFiles.length; i++){
			parser.readBarFile(barFiles[i], true);
			float[] delogged = Num.antiLog(parser.getValues(), 2);
			//float[] delogged = parser.getValues();
			System.out.println("\n"+barFiles[i].getName());
			Num.statFloatArray(delogged, false);
		}
	
	
	}
}
