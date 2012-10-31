package util.bio.seq;
import util.gen.*;
import java.io.*;
import java.util.regex.*;

/**Converts a directory of files that contain lowercase masked sequences to Xs.*/
public class LowerCase2Xs {

	public static void main(String[] args) {
		//check args
		if (args.length ==0) Misc.printExit("\nEnter the full path directory or file with lowercase " +
				"sequence to convert and optionally, the character to replace them with, defaults to 'X'.\n");
		
		//fetch mfa files
		File[] mfaFiles = IO.extractFiles(new File(args[0]));
		
		//character to use in replacement
		String replacementString = "X";
		if (args.length == 2) replacementString = args[1];
		System.out.println("Replacement character -> "+replacementString);
		
		//for each file
		Pattern pat = Pattern.compile("[gatc]");
		Matcher mat; 
		for (int i=0; i< mfaFiles.length; i++){
			System.out.println("Processing "+mfaFiles[i].getName());
			try {
				BufferedReader in = IO.fetchBufferedReader(mfaFiles[i]);
				String name = mfaFiles[i].getCanonicalPath();
				PrintWriter out = new PrintWriter(new FileWriter(new File(name+".Xed")));
				String line;
				//for each line
				while ((line = in.readLine()) !=null){
					//check it aint a header line
					if (line.startsWith(">") == false){
						mat = pat.matcher(line);
						line = mat.replaceAll(replacementString);
					}
					out.println(line);
				}
				out.close();
				in.close();
			} catch (IOException e){
				e.printStackTrace();
			}
		}
		
		System.out.println("\nDone!\n");
	}

}
