package util.apps;
import util.gen.*;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import trans.anno.GenomicRegion;
import trans.cel.MultiSetQuantile;
import trans.main.Window;
import trans.roc.Sgr;

public class TASSignalChromFileSplitter {
	
	public static void main (String[] args){
		try{
			File file = new File(args[0]);
			PrintWriter out = null;
			String[] lines = IO.loadFileIntoStringArray(file);
			String chr = null;
			for (int i=0; i<lines.length; i++){
				//hit new line
				if (lines[i].indexOf("# Name")!= -1){
					System.out.println(lines[i]);
					if (out !=null) out.close();
					int index = lines[i].indexOf("c");
					chr = lines[i].substring(index);
					String name = file.getCanonicalPath()+chr+".txt";
					out = new PrintWriter (new FileWriter (name));
				}
				else if (lines[i].startsWith("#") == false && Misc.isNotEmpty(lines[i])) out.println(chr+"\t"+lines[i]);
			}
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
	
	
	
	
}
