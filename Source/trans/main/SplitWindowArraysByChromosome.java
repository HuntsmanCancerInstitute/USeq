package trans.main;
import java.util.regex.*;
import java.util.*;
import util.gen.*;
import java.io.*;

public class SplitWindowArraysByChromosome {
	
	private File[] windowFiles;
	
	public SplitWindowArraysByChromosome(String[] args){
		System.out.println("\nLaunched:");
		processArgs(args);
		WindowComparator wc = new WindowComparator();
		
		for (int i=0; i< windowFiles.length; i++){
			System.out.println("\tLoading "+windowFiles[i].getName());
			Window[] win = (Window[])IO.fetchObject(windowFiles[i]);
			Arrays.sort(win, wc);
			String currentChrom = win[0].getChromosome();
			ArrayList winAL = new ArrayList();
			System.out.println("\t\t"+win.length+"\tTotal Windows");
			for (int j=0; j< win.length; j++){
				if (win[j].getChromosome().equals(currentChrom)) winAL.add(win[j]);
				else {
					//convert to Window[]
					Window[] chromWin = new Window[winAL.size()];
					winAL.toArray(chromWin);
					File f = new File (windowFiles[i].getParentFile(), currentChrom+"_"+windowFiles[i].getName());
					System.out.println("\t\t"+chromWin.length+"\tSaved for "+f.getName());
					IO.saveObject(f, chromWin);
					chromWin = null;
					//begin anew
					currentChrom = win[j].getChromosome();
					winAL.clear();
					winAL.add(win[j]);
				}
			}
			//save last
			Window[] chromWin = new Window[winAL.size()];
			winAL.toArray(chromWin);
			File f = new File (windowFiles[i].getParentFile(), currentChrom+"_"+windowFiles[i].getName());
			System.out.println("\t\t"+chromWin.length+"\tSaved for "+f.getName());
			IO.saveObject(f, chromWin);
			chromWin = null;
			winAL = null;
		}
		System.out.println("\nDone!\n");
	}
	

	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new SplitWindowArraysByChromosome(args);
	}
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		File file = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': file = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (file ==null || file.isDirectory() == false){
			Misc.printExit("\nPlease enter a directory containing window arrays.\n");
		}
		windowFiles = IO.extractFiles(file);
	}	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Split Window Arrays by Chromosome: Nov 2006                  **\n" +
				"**************************************************************************************\n" +
				"Splits and saves serialized Window[]s by chromosome.\n\n" +
				
				"-f Full path file or directory text for the serialized Window[]s.\n\n" +
				
				"Example: java -Xmx5000M -jar pathTo/T2/Apps/SplitWindowArraysByChromosome -f /af/win/ \n\n" +
				
		"**************************************************************************************\n");		
	}

}
