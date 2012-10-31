package trans.main;

import java.io.*;
import util.gen.*;
import java.util.*;

public class ReplaceWindowScores {
	
	//fields
	private int[] cyberTColumnsToParse= {18,26};
	private boolean[] convertScores = {false,false};	//-10*log10(value) transform?
	private CyberTLine[] cyberTLines;
	private Window[] windows;
	
	//constructor
	public ReplaceWindowScores (String[] args){
		//fetch cyber t lines, ordered one per window
		File cyberTFile = new File (args[0]);
		parseCyberT(cyberTFile);
		//fetch Window[]
		File windowsFile = new File (args[1]);
		windows = (Window[])IO.fetchObject(windowsFile);
		//check same number
		if (cyberTLines.length != windows.length) Misc.printExit("\nError: number of cyber t lines and windows differ!");
		
		//replace window scores with cyber t scores
		for (int i=0; i< windows.length; i++){
			if (i<5){
				System.out.println(cyberTLines[i].name+"\t"+ Misc.doubleArrayToString(cyberTLines[i].scores, "\t")+
						"\tWin "+ windows[i].getChromosome()+"\t"+
						Misc.doubleArrayToString(windows[i].getScores(),"\t"));
			}
			windows[i].setScores(cyberTLines[i].scores);
		}
		
		//save window array
		File modWindows = new File(IO.getFullPathName(windowsFile)+".rpl");
		IO.saveObject(modWindows, windows);
		
	}
	
	public static void main (String[] args){
		if (args.length==0) Misc.printExit("\nEnter full path file names for the cyberTFile and Window[] file.\n");
		new ReplaceWindowScores(args);
	}
	
	/**Parses a file for a tab delimited list of at minimum, chrom, start, stop, (optional) notes into a
	 * GenomicRegion[].*/
	public CyberTLine[] parseCyberT(File file){
		String line = null;
		ArrayList ctLinesAL = new ArrayList();
		try{
			BufferedReader in = new BufferedReader(new FileReader(file));    
			while ((line = in.readLine()) !=null) {
				line = line.trim().replaceAll("\"","");
				if (line.length() == 0) continue;
				ctLinesAL.add(new CyberTLine(line));
			}
			cyberTLines = new CyberTLine[ctLinesAL.size()];
			ctLinesAL.toArray(cyberTLines);
			Arrays.sort(cyberTLines);
			return cyberTLines;
			
		}catch (IOException e){
			System.out.println ("Problem parsing this line-> "+line);
			e.printStackTrace();
		}
		return null;
	}
	
	private class CyberTLine implements Comparable{
		int index;
		double[] scores;
		String name;
		
		public CyberTLine(String line){
			String[] tokens = line.split("\t");
			index = Integer.parseInt(tokens[0]);
			name = tokens[1];
			//parse scores
			double log10 = Math.log(10);
			scores = new double[cyberTColumnsToParse.length];
			for (int i=0; i<cyberTColumnsToParse.length; i++){
				scores[i] = Double.parseDouble(tokens[cyberTColumnsToParse[i]]);
				//transform
				if (convertScores[i]) scores[i] = -10*(Math.log(scores[i])/log10);
			}
		}
		//methods
		public int compareTo(Object obj){
			CyberTLine other = (CyberTLine)obj;
			if (other.index>index) return -1;
			if (other.index<index) return 1;
			return 0;
		}
	}
}