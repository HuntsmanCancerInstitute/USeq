package trans.tpmap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;
import util.gen.*;

/**
 * This reversed the oligo sequence and orientation in a tpmap file.
 */
public class OligoFlipper {
	//fields
	private File[] tpmapFiles;
	private File resultsDirectory;
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': tpmapFiles = IO.extractFiles(new File(args[i+1]), ".tpmap"); i++; break;
					case 'r': resultsDirectory = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            OligoFlipper: Dec 2005                                **\n" +
				"**************************************************************************************\n\n" +
				
				"-f Full path file text for the text bpmap file.\n" +
				"-r ResultsDirectory.\n"+
				
				"\n" +
				
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length ==0 ){
			printDocs();
			System.exit(0);
		}
		System.out.println("Launching OligoFlipper...");
		new OligoFlipper(args);
	}
	
	
	public OligoFlipper(String[] args){
		processArgs(args);
		
		try {
			long start = System.currentTimeMillis();
			String line;
			for (int i=0; i< tpmapFiles.length; i++){
			BufferedReader in = new BufferedReader(new FileReader(tpmapFiles[i]));
			PrintWriter out = new PrintWriter ( new FileWriter (new File (resultsDirectory, tpmapFiles[i].getName())  ) );
			String[] tokens;
			while ((line = in.readLine()) !=null) { 
				if (line.startsWith("#")) out.println(line);
				else{
					tokens = line.split("\\s+");
					tokens[0] = new StringBuffer(tokens[0]).reverse().toString();
					if (tokens[1].equals("f")) tokens[1] = "t";
					else tokens[1] = "f";
					out.println(Misc.stringArrayToString(tokens, "\t"));
				}
			}
			in.close();
			out.close();
			}
			
			int elapse = (int)(System.currentTimeMillis()-start)/1000;
			System.out.println("Finished "+elapse+" seconds");
		} catch (Exception e) {e.printStackTrace();}
		
	}
	
	
	
}
